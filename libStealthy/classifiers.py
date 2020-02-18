#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from abc import ABCMeta, abstractmethod
from scipy.stats import ttest_ind, mannwhitneyu, entropy, chisquare

from libStealthy.CCE import libCCE
from libStealthy.util import window, equalAreaBins
from libStealthy.srl import calcProbLow, statisticNu0IPD, statisticSignedNu0IPD

import numpy as np
import collections
import sys


class baseClassifier(object, metaclass=ABCMeta):
    '''
    This abstract class defines the basic structure of a classifier,
    `__init__` should be use to setup any classifier specfic
    arguments.
    '''
    @abstractmethod
    def setup(self, refIPDs):
        '''
        Used to "calibrate" the classifer with the refference IPDs.
        This method is called once per ref file.

        Args:
            refIPDs (iterable):
                An iterable containing the refference IPDs to be
                used.
        '''
        pass

    @abstractmethod
    def classify(self, obsIPDs):
        '''
        Calculate and return the test stat for the classifier.


        Args:
            obsIPDs (iterable):
                An iterable containing the observed IPDs to be used.
        
        Returns:
            (tuple of floats):A tuple containing one or more floats
            that describe the output of the test stat
        '''
        pass

    @abstractmethod
    def getMetaData(self):
        '''
        Get metadata related to this specfic intance of the
        clssifier.

        Returns:
            (dict):A dict containing classifier specfic metadata 
        '''
        pass

class welcht(baseClassifier):
    detectorStr = 'welcht'
    def __init__(self):
        pass
    def setup(self, refIPDs):
        self.refIPDs = refIPDs

    def classify(self, obsIPDs):
        t, p = ttest_ind(obsIPDs, self.refIPDs, equal_var=False)
        return (t.item(), p)
    def getMetaData(self):
        return {
            'detectorStr':self.detectorStr
        }

class mwu(baseClassifier):
    detectorStr = 'mwu'
    def __init__(self):
        pass

    def setup(self, refIPDs):
        self.refIPDs = refIPDs

    def classify(self, obsIPDs):
        U, p = mannwhitneyu(obsIPDs, self.refIPDs, alternative='two-sided')
        return (U.item(), p)
    def getMetaData(self):
        return {
            'detectorStr':self.detectorStr
        }


class chisq(baseClassifier):
    detectorStr = 'chisq'
    def __init__(self, numOfBins):
        self.numBins = numOfBins

    def setup(self, refIPDs):
        self.refIPDs = refIPDs
        self.maxEnt = max(self.refIPDs)
        self.minEnt = min(self.refIPDs)
        # Use bins of uniform width
        self.bin_width = (self.maxEnt - self.minEnt) / self.numBins

        tmp_ranges = []
        for iter in range(0, self.numBins):
            # Floats are rounded to 6 decimal places
            tmp_ranges.append(float("{0:.6f}".format(self.minEnt + iter * self.bin_width)))
        tmp_ranges.append(float("{0:.6f}".format(self.maxEnt)))
        tmp_ranges = np.array(tmp_ranges)

        ref_hist, ref_bin_edges = np.histogram(self.refIPDs, bins=tmp_ranges)

        # If histogram has any 0 entries, merge those bins into the bins
        # into the bins to their left.
        # (Leftmost bin includes min. element, so this is OK.)
        if 0 in ref_hist:        
            to_delete = []
            for bin_num in range(0, len(ref_hist)):
                if ref_hist[bin_num] == 0:
                    to_delete.append(bin_num)
            ref_hist = np.delete(ref_hist, to_delete)
            tmp_ranges = np.delete(tmp_ranges, to_delete)
        
        self.ref_hist = ref_hist
        self.bin_edges = tmp_ranges


    def classify(self, obsIPDs):
        '''
        Return the chi-squared statistic and p value.
        
        Compare the observed sequence obsIPDs with the reference sequence
        used in setting up the classifier using the chi-squared statistic.
        
        Args:
            obsIPDs (iterable):
                The observed IPD sequence for which the statistic
                should be calculated.

        Returns:
            (chisq, p):
                A pair containing the chi-squared statistic chisq and the
                associated p-value p.
        '''

        # Expand first and last bins, if needed, to include all obsIPDs
        obs_edges = self.bin_edges
        if min(obsIPDs) < obs_edges[0]:
            obs_edges[0] = min(obsIPDs)
        if max(obsIPDs) > obs_edges[-1]:
            obs_edges[-1] = max(obsIPDs)
        obs_hist, obs_bin_edges = np.histogram(obsIPDs, bins=obs_edges)

        # Rescale the reference histogram to get the values we would expect
        # if there were as many refIPDs as there are obsIPDs
        scale = len(obsIPDs) / len(self.refIPDs)
        scaled_ref_hist = [scale * hist_val for hist_val in self.ref_hist]
        
        C, p = chisquare(obs_hist, scaled_ref_hist)
        return (C.item(), p)

    def getMetaData(self):
        return {
            'detectorStr':self.detectorStr
        }



class filler(baseClassifier):
    detectorStr = 'filler'
    def __init__(self, fillerP = None):
        if fillerP is not None:
            self.probLow = fillerP
        else:
            self.fillerP = fillerP

    def setup(self, refIPDs):
        if self.fillerP is None:
            self.probLow = calcProbLow(tuple(refIPDs))

    def classify(self, obsIPDs):
        nu = statisticNu0IPD(obsIPDs, self.probLow)
        return (nu,)
    def getMetaData(self):
        return {
            'probLow':self.probLow
            ,'detectorStr':self.detectorStr
        }

class sgfiller(baseClassifier):
    detectorStr = 'sgfiller'
    def __init__(self, fillerP = None):
        if fillerP is not None:
            self.probLow = fillerP
        else:
            self.fillerP = fillerP

    def setup(self, refIPDs):
        if self.fillerP is None:
            self.probLow = calcProbLow(tuple(refIPDs))

    def classify(self, obsIPDs):
        sgnu = statisticSignedNu0IPD(obsIPDs, self.probLow)
        return (sgnu,)

    def getMetaData(self):
        return {
            'probLow':self.probLow
            ,'detectorStr':self.detectorStr
        }

class cce(baseClassifier):
    '''
    A classifier which calcs Corrected Conditional Entropy to decide
    if a given sequence is embedded or plain.

    Args:
        seqLen (int):
            The length of the sub-sequences you would like to
            break the observed input sequences into.  Making this
            long will increase compute time sigficantly.
        bins (int or (list or tuple)):
            If bins is an int it is used to set how many
            equal area bins the input will be binned into.

            If bins is a list or tuple then we use it as the
            boundries for the binning process
    '''
    def __init__(self, seqLen, bins):
        self.detectorStr = 'cce'
        if type(bins) is int:
            self.binCount = bins
        elif len(bins) == 1:
            self.binCount = int(bins[0])
        else:
            self.bins = [float(bin) for bin in bins]
            self.binCount = None
        self.seqLen = seqLen

    def setup(self, refIPDs):
        '''
        Use the refIPDs to calculate the bins for the obsIPDs unless
        we already have bins set from init

        Args:
            refIPDs (tuple or list):
                The refference IPDs you want to use when setting up
                bins.
        '''
        # After init and this, self.bins has the left endpoints of our bins
        if self.binCount is not None: # I.e., single value given for binCount in init
            # create self.binCount equal area bins (each contains 
            # same number of refIPDs)
            self.bins = equalAreaBins(refIPDs, self.binCount)
            self.bins[0] = 0.0

    def classify(self, obsIPDs):
        '''
        Compute the CCE of the observed IPDs and return the minimum.

        Args:
            obsIPDs (iterable):
                The un-binned observed IPD sequence you want to
                calculate CCE for.

        Returns:
            (tuple with a float):
                A Tuple containing the minimum CCE value of the
                observed IPDs
        '''
        # -1 ensures that all values are between 0 and the number of bins -1
        obsIPDs = np.digitize(obsIPDs, self.bins) - 1
        # That should take sequence of bin values and compute CCE
        obsCCE = libCCE.CCE(self.binCount, obsIPDs, self.seqLen)
        oCCE = obsCCE.calculateCCE()
        return (min(oCCE),)

    def getMetaData(self):
        '''
        Get metadata related to this specfic intance of the
        clssifier.
        This detector has the metadata of:

        * binCount
        * bins
        * sequenceLen
        * detectorStr

        Returns:
            (dict):A dict containing detector specfic metadata 
        '''
        return {
            'binCount': len(self.bins)
            ,'bins': [str(bin_) for bin_ in self.bins]
            ,'sequenceLen': self.seqLen
            ,'detectorStr':self.detectorStr
        }

class localMinCCE(cce):
    detectorStr = 'localMinCCE'
    
    def __init__(self, seqLen, bins):
        self.detectorStr = 'localMinCCE'
        if type(bins) is int:
            self.binCount = bins
        elif len(bins) == 1:
            self.binCount = int(bins[0])
        else:
            self.bins = [float(bin) for bin in bins]
            self.binCount = None
        self.seqLen = seqLen
    
    
    def classify(self, obsIPDs):
        obsIPDs = np.digitize(obsIPDs, self.bins) - 1
        obsCCE = libCCE.localMinCCE(tuple(obsIPDs), self.seqLen)
        return (obsCCE,)

#
# This was an example for writing a new classifier; not used anywhere
#
class Entropy(baseClassifier):
    '''
    Class that was added as an example.  (This is not used anywhere.)
    '''
    detectorStr = 'entropy'
    def __init__(self, bins):
        if type(bins) is int:
            self.binCount = bins
        elif len(bins) == 1:
            self.binCount = int(bins[0])
        else:
            self.bins = [float(bin) for bin in bins]
            self.binCount = None

    def setup(self, refIPDs):
        if self.binCount is not None:
            # create self.binCount equal area bins (each contains 
            # same number of refIPDs)
            self.bins = equalAreaBins(refIPDs, self.binCount)
            self.bins[0] = 0.0
    
    def classify(self, obsIPDs):
        bindedIPDs = np.digitize(obsIPDs, self.bins)
        ent = entropy(bindedIPDs)
        return (ent,)
    
    def getMetaData(self):
            return {'detectorStr':self.detectorStr}

# Here be introspective dragons... 
# We use some of Python's inrospection features to build
# references to classes from this module's internal dict.
#
def getClassifier(classifier):
    '''
    Get a reference to a classifier.

    Args:
        rngName (str):
            A string containing the name of the classifier you want.
    Returns:
        A reference to the constructor of the requested class
    
    Raises:
        NotImplementedError:
            Raised if the provided string doesn't correspond to
            anything in this module
    '''
    # get a ref to the current module (classifiers)
    currentModule = sys.modules[__name__]
    # check to see if the current module has anything matchng the
    # supplied name.
    if hasattr(currentModule, classifier):
        # return a ref to the object if so
        return getattr(currentModule, classifier)
    else:
        # raise NotImplmentedError if we don't have that classifier
        raise NotImplementedError('Classifier %d is not implmented.' % classifier)

def getClassifiers():
    '''
    Get a tuple of the names of the classfiers implmented by this
    module.

    Returns:
        (tuple):A tuple of the names of the classifiers implmented by
        this module.
    '''
    # get a ref to the current module (classifiers)
    currentModule = sys.modules[__name__]
    classifiers = []
    # create a generator to get memebers from the module by name.
    moduleMembers = ((getattr(currentModule, name), name) for name in dir(currentModule))
    for member, name in moduleMembers:
        # remove things that aren't classes
        if isinstance(member, type):
            # remove things that aren't subclassed from
            # baseClassifier
            if issubclass(member, baseClassifier):
                classifiers.append(name)
    return tuple(classifiers)
