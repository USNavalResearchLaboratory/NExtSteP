#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from math import sqrt, log10, log
from Crypto.Random import random as HQrandom
from base64 import b64decode, b64encode
from itertools import islice
from pathlib import Path
from array import array
import configparser
import numpy as np
import datetime
import json
import sys
import re
import os


# These should all be able to take a float; rounding will be done later.
rateF = {}
rateF['const'] = lambda N : 1.0
rateF['log'] = lambda N : log10(N)
rateF['cubrt'] = lambda N :pow(N, 1./3.)
rateF['sqrt'] = lambda N : pow(N, 1./2.)
rateF['2rt3'] = lambda N : pow(N, 2./3.)
rateF['linear'] = lambda N : N


def calc_rate_constants(n, k):
    '''
    Determine rate constants so that k out of n IPDs are changed.
    '''
    rateC = {}
    for fname in rateF:
        rateC[fname] = float(k) / rateF[fname](n)
    return rateC

def calc_mssage_length(N, rateType, rateConstant):
    '''
    Determine how many bits to change (out of N total) for a given embedding rateType.
    '''
    return int(round(rateConstant * rateF[rateType](N)))


def calc_mssage_prob(N, rateType, rateConstant):
    '''
    Determine probability of changing a  given embedding rateType.
    '''
    # Do computation analogous to calc_mssage_length without rounding
    # Divide result by total cover length N
    return rateConstant * rateF[rateType](N) / N


def hqrand():
    '''
    Generate a high-quality random number in [0,1]

    Returns:
        (float):A random floating point number between [0, 1]
        
    '''
    # FUTURE: The number of random bits used here should be changeable

    return float(HQrandom.getrandbits(31))/(2**31 - 1)


def aur(classifier_outputs, unrev=False):
    '''
    Compute the area under the ROC given the outputs of a classifier.

    classifier_outputs is a 0-1 list whose entries indicate whether the
    tested condition is true; these are sorted according to the outputs
    of a classifier.  E.g., [1,0,0,1,0,1] captures six samples with three
    negative and three positive such that the classifier ranked one of the
    positive samples as most likely to be positive, then two of the negative
    samples, etc.

    This returns a value between 0.5 and 1.0 corresponding to the area
    under the ROC curve specified by these classifier outputs.  (If the value
    would be less than 0.5, it is first subtracted from 1.0.)
    '''

    ones_so_far = 0
    zeros = 0
    temp_AUR = 0

    for x in classifier_outputs:
        if x == 0:
            temp_AUR += ones_so_far
            zeros += 1
        elif x == 1:
            ones_so_far += 1
        else:
            raise ValueError('Value other than 0 or 1 encountered!')

    rect = ones_so_far * zeros

    if 2 * temp_AUR >= rect or unrev:
        return temp_AUR/float(rect)
    else:
        return (rect - temp_AUR)/float(rect)


def pe(classifier_outputs):
    '''
    Compute the 1-P_E value (following Ker et al.) of a binary classifier.

    classifier_outputs is a 0-1 list whose entries indicate whether the
    tested condition is true; these are sorted according to the outputs
    of a classifier.  E.g., [1,0,0,1,0,1] captures six samples with three
    negative and three positive such that the classifier ranked one of the
    positive samples as most likely to be positive, then two of the negative
    samples, etc.

    This returns 1 - P_E, where P_E is one half of the minimum error fraction
    (following Ker et al.).  I.e., we minimize the sum of the number of false
    positives plus the number of false negatives and divide this by 2 * the
    total number of sequences (plain plus stego); the factor of 2 comes from
    Ker et al.
    
    Here, we compute the error for both the standard detector and its
    complement, i.e., in which the interpretation of above/below the
    threshold is reversed.  We return the smaller of these errors; this is
    analogous to complementing the detector if the AUR is less than 0.5.
    '''
    # Go through list
    # At each point, count:
    # #zeros before plus #ones after
    # #ones before plus #zeros after
    # Find the minimum of each sum
    # Modify (divide by two, subtract from one) to get measure of interest
    # Return measure of interest

    num_zeros = classifier_outputs.count(0)
    num_ones = classifier_outputs.count(1)

    zeros_so_far = 0
    ones_so_far = 0
    zeros_left = num_zeros
    ones_left = num_ones
    
    # 01_sum is 0s before, 1s after
    # 10_sum is 1s before, 0s after
    min_01_sum = zeros_so_far + ones_left
    cur_01_sum = min_01_sum

    min_10_sum = zeros_so_far + ones_left
    cur_10_sum = min_10_sum

    for next_true in classifier_outputs:
        if next_true == 0:
            zeros_so_far += 1
            zeros_left -= 1
        elif next_true == 1:
            ones_so_far += 1
            ones_left -= 1
        else:
            raise ValueError('Classifier output that is neither 0 nor 1!')
        cur_01_sum = zeros_so_far + ones_left
        cur_10_sum = ones_so_far + zeros_left
        
        if cur_01_sum < min_01_sum:
            min_01_sum = cur_01_sum
        if cur_10_sum < min_10_sum:
            min_10_sum = cur_10_sum

    if min_01_sum <= min_10_sum:
        min_sum = min_01_sum
    else:
        min_sum = min_10_sum

    PsubE = float(min_sum) / (2. * len(classifier_outputs))
    return 1 - PsubE


def pe_pair(classifier_outputs):
    '''
    Compute the 1-P_E values of a binary classifier and its complement.

    classifier_outputs is a 0-1 list whose entries indicate whether the
    tested condition is true; these are sorted according to the outputs
    of a classifier.  E.g., [1,0,0,1,0,1] captures six samples with three
    negative and three positive such that the classifier ranked one of the
    positive samples as most likely to be positive, then two of the negative
    samples, etc.

    This returns a pair of 1 - P_E values, where P_E is one half of the
    minimum error fraction.  I.e., we minimize the sum of the number of false
    positives plus the number of false negatives and divide this by 2 * the
    total number of sequences (plain plus stego); the factor of 2 comes from
    Ker et al.
    
    Here, we compute the error for both the standard detector and its
    complement, i.e., in which the interpretation of above/below the
    threshold is reversed.  We return a pair in which the first 1-P_E value
    corresponds to 0 before 1 as an error and the second value corresponds to
    1 before 0 as an error.
    '''
    # Go through list
    # At each point, count:
    # #zeros before plus #ones after
    # #ones before plus #zeros after
    # Find the minimum of each sum
    # Modify (divide by two, subtract from one) to get measure of interest
    # Return measure of interest

    num_zeros = classifier_outputs.count(0)
    num_ones = classifier_outputs.count(1)

    zeros_so_far = 0
    ones_so_far = 0
    zeros_left = num_zeros
    ones_left = num_ones
    
    # 01_sum is 0s before, 1s after
    # 10_sum is 1s before, 0s after
    min_01_sum = zeros_so_far + ones_left
    cur_01_sum = min_01_sum

    min_10_sum = zeros_so_far + ones_left
    cur_10_sum = min_10_sum

    for next_true in classifier_outputs:
        if next_true == 0:
            zeros_so_far += 1
            zeros_left -= 1
        elif next_true == 1:
            ones_so_far += 1
            ones_left -= 1
        else:
            raise ValueError('Classifier output that is neither 0 nor 1!')
        cur_01_sum = zeros_so_far + ones_left
        cur_10_sum = ones_so_far + zeros_left
        
        if cur_01_sum < min_01_sum:
            min_01_sum = cur_01_sum
        if cur_10_sum < min_10_sum:
            min_10_sum = cur_10_sum


    PsubE_01 = float(min_01_sum) / (2. * len(classifier_outputs))
    PsubE_10 = float(min_10_sum) / (2. * len(classifier_outputs))
    return (1 - PsubE_01, 1 - PsubE_10)


def window(seq, n=3):
    '''
    Returns a sliding window (of width n) over data from the iterable"
    ::

       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...

    Args:
        seq (iterable):
            The sequence you would like to get a sliding window over.
        n (int):
            The number of elements in each window.
    Yields:
        (tuple):A tuple of length n containing the elements of the list.
    
    '''
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it: # it is iterator
        result = result[1:] + (elem,)
        yield result


class DataFile(object):
    '''
    Abstracts creating and reading datafiles.

    Args:
        ifile (str):
            The path to a DataFile on disk
        fast (bool):
            Sould fast encoding be the default.  By default this is
            `False` because fast encoding should only be used with
            lists of numbers.


    The DataFile uses a JSON schema like the following:
    ::

        {
            //Metadata describing the file's creation.
            "metadata":{
                // The datetime the file was created.
                "created": "2018-07-30 10:35:09.257966"
                // Does the file use fast output encoding?
                ,"fast": true
                // Other assorted metadata like the detector string or
                //embedding rate.
            },
            # A list or a string depending on if fast encoding was used
            "data":[[1.3,1], [2.4, 0]]
        }

    Raises:
        TypeError:
            Raised if the file is not a DataFile.
    '''
    fileNameEndngRegex = r'\d{3,}.json'
    '''
    Default regex for the end of the file, looks for at least three
    digits followed by the literal string ".json".
    '''
    __slots__ = ['dataFile', 'fast']
    def __init__(self, ifile = None, fast = False):
        if ifile is not None:
            self.load(ifile)
        else:
            self.dataFile = {
                    'metadata':{
                        'created':str(datetime.datetime.now())
                        ,'fast':fast
                        ,'type':'DataFile'
                    },
                    'data':[]
                }
        self.fast = fast
    
    def load(self, fname):
        '''
        Load a datafile

        Args:
            fname (str):
                The path to the datafile
        Raises:
            TypeError:
                Raised if the file is not a DataFile.
        '''
        with open(fname) as fp:
            loaded = json.load(fp)
        # check to see if this is an older datafile and handle it.
        if isinstance(loaded, dict):
            # check to see if this is actually a DataFile
            if 'metadata' in loaded and 'type' in loaded:
                raise TypeError('File\'s %s structure does not match any known filetype' % fname)
            if loaded['metadata']['type'] != 'DataFile':
                raise TypeError('File %s is a %s not a DataFile' % (fname, str(loaded['metadata']['type'])))
            self.dataFile = loaded
        elif isinstance(loaded, (list, tuple)):
            self.dataFile = {
                'metadata':{
                    'created':str(datetime.datetime.now())
                    ,'fast':False
                    ,'type':'DataFile'
                },
                'data':loaded
            }
        else:
            raise TypeError('File\'s %s structure does not match any known filetype' % fname)
        # handle a fast encoded datafile
        if self.getMetaData('fast'):
            arr = array('d')
            encoded = self.dataFile['data']
            decoded = b64decode(encoded)
            # convert to bytes object
            if isinstance(decoded, str):
                decoded = decoded.encode('ascii')
            arr.frombytes(decoded)
            self.setData(list(arr))
            self.removeMetaData('fast')

    def write(self, fname, fast=None):
        '''
        Serilize this DataFile to a file.
        
        Args:
            fname (str): 
                The path to the datafile
            fast (bool):
                Sould fast encoding be used, by default uses whatever
                encoding the DataFile was created with but can be
                overridden here.
        '''
        if fast is None:
            fast = self.fast
        
        # update timestamp
        self.addMetaData('created',str(datetime.datetime.now()))

        if fast:
            data = self.getData()
            self.addMetaData('fast', True)
            arr = array('d', data)
            encoded = b64encode(arr.tobytes())
            # convert from bytes to str or json encoding.
            encoded = encoded.decode('ascii')
            self.setData(encoded)
            with open(fname, 'w') as fp:
                json.dump(self.dataFile, fp)
            self.setData(data)
        else:
            with open(fname, 'w') as fp:
                json.dump(self.dataFile, fp)

    def getData(self):
        '''
        Get the data stored in this DataFile

        Returns:
            list:
                A list containing the data stored in the DataFie.
        '''
        return self.dataFile['data']

    def setData(self, data):
        '''
        Set the contents of this DataFile

        Args:
            data:
                Some data for the DataFile to store, if `self.fast`
                is `False` any JSON serilizable data can be set, if
                `self.fast` is `True` then only arrays of numbers can
                be stored.
        '''
        self.dataFile['data'] = data

    def getMetaData(self, attribute = None):
        '''
        Get the metadata from this instance, if attrubute is provided
        only that attrbute of the metadata will be returned.

        Args:
            attribute (str):
                The name of the attrubute you want.
        '''
        if attribute is None:
            return self.dataFile['metadata']
        else:
            return self.dataFile['metadata'][attribute]

    def addMetaData(self, attribute, value):
        '''
        Add some data to the metadata of the current DataFile

        Args:
            attribute (str):
                The name of the attrubute you want to set
            value (int, float, str):
                The value to be stored, may only be types that can be
                directly mapped to json
        '''
        self.dataFile['metadata'][attribute] = value

    def removeMetaData(self, attribute):
        '''
        Remove a specfic attrubute of the metadata.

        Args:
            attribute (str):
                The name of the attrubute you would like to remove.
        
        Raises:
            KeyError:
                When the provided attribute is not in the metadata
        '''
        del self.dataFile['metadata'][attribute]
        

    def mergeMetaData(self, metaDataToMerge):
        '''
        Add the metadata from an external dict to this DataFile. The
        external dict's values will replace any tags that already 
        happen to exist in this DataFile.

        Args:
            metaDataToMerge (dict):
                A dict conaining the metadata you want to merge into
                this DataFile.
        '''
        self.dataFile['metadata'].update(metaDataToMerge)

    @staticmethod
    def getValidDataFiles(inputDirectory, inputPrefix, matchingExpression = None):
        '''
        Yields paths to valid input files for the given input directory
        and prefix.

        Args:
            inputDirectory (str):
                The path to the data directory to search
            inputPrefix (str):
                The name of the data file before the index and .json file
                type.
            
            matchingExpression (str or None):
                A regular expression to match the remainder of the
                filename after prefix, by default the matching
                expression consists of the input prefix followed by
                upto three digits followd by the literal ".json"

        Yields:
            (str):Path to a matching data file'''
        if matchingExpression is None:
            matchingExpression = DataFile.fileNameEndngRegex

        inputPrefix = re.escape(inputPrefix)
        regexString = inputPrefix + matchingExpression
        regex = re.compile(regexString)


        for inputFile in sorted(filter(regex.search, os.listdir(inputDirectory))):
            yield os.path.join(inputDirectory, inputFile)


class SeedBank(object):
    '''
    Handles storage and retrieval of seeds from a "SeedBank" file.
    

    Args:
        bankPath (str):
            The path to the seedbank file
        runID (str or None):
            The ID of a previous run you want to re-use the seeds
            from.  If None, new seeds are generated. If the runID
            doesn't exist in this seedbank generate new seeds and
            store them under that runID.
    
    The SeedBank follows the following schema::

        {
            "data": {
                "runID1": {
                    "owner1": "seed1"
                    ,"owner2": "seed2"
                }
                ,"runID2": {
                    "subRunID1": {
                        "owner1": "seed3"
                        ,"owner2": "seed4"
                    }
                }
            }
            ,"metadata": {
                "created": "2018-08-06 10:58:21.473026"
                ,"type": "SeedBank"
            }
        }


    Where runID is a name given to that particular generation of
    seeds, by default this is the date and time the new seeds were
    requested.  subRunID is some sub-specifier for the run, for
    instance the cover lengths.
    
    An example:

    ::

        with SeedBank("SeedBank.json", 'myTest') as bank:
            files = ["a.txt", "b.txt", "c.txt"]
            for file, seed in bank.getNewSeeds(files):
                print(file, seed)

    '''
    def __init__(self, bankPath, runID = None):
        self.bankPath = bankPath
        if runID is None:
            self.runID = str(datetime.datetime.now()).replace(' ', '_')
        else:
            self.runID = runID

    def __enter__(self):
       self.open()
       return self

    def __exit__(self, type, value, traceback):
        self.write()
    
    def open(self, path = None):
        '''
        Opens the BankFile specfied on creation.  If no such BankFile
        exists we create one.

        Args:
            path (str):
                The path to read the bankfile from.  If not provided
                the bankpath provided in __init__ will be used. If a
                path is provided self.bankPath will be updated to
                match.

        Raises:
            TypeError:
                Raised if the metadata type is not "SeedBank".
        '''
        if path is not None:
            self.bankPath = path
        try:
            with open(self.bankPath, 'r') as fp:
                bank = json.load(fp)
            if 'metadata' in bank and 'type' in bank:
                raise TypeError('File %s is not any known filetype' % self.bankPath)
            if bank['metadata']['type'] != 'SeedBank':
                raise TypeError('[!!!] %s is a %s not a SeedBank')
            self.bank = bank
        except FileNotFoundError:
            # populate the internal dict
            self.bank = {
                'metadata':{
                    'created':str(datetime.datetime.now())
                    ,'type':'SeedBank'
                }
                ,'data':{}
            }
            with open(self.bankPath, 'w') as fp:
                json.dump(self.bank, fp)
        
        if self.runID not in self.bank:
            self.bank['data'][self.runID] = {}

    def write(self, path = None):
        '''
        Write all seeds to the bankfile and close the file handle to
        it.
        
        Args:
            path (str):
                The path to write the bankfile to, uses self.bankpath
                by default.
        '''
        if path is None:
            path = self.bankPath
        with open(path, 'w') as fp:
            json.dump(self.bank, fp)
        

    def getSeedsByID(self, subRunID = None):
        '''
        Get existing seeds from SeedBank.

        Args:
            subRunID (str):
                The subrun ID of the specfic group of seeds you want.
        
        Yields:
            (tuple):A tuple containing the name of the file the seed
            was used for and the seed.
        '''
        if subRunID is not None:
            seeds = self.bank['data'][self.runID][str(subRunID)]
        else:
            seeds = self.bank['data'][self.runID]
        for owner, seed in seeds.items():
            yield (owner, seed)

    def getNewSeeds(self, owners, subRunID = None, seedBits = 256):
        '''
        Generate new seeds and add them to the bank.

        Args:
            owners (list or tuple):
                A list or tuple of values representing the owner
                of the seed to be generated. In our use case the
                values of the list or tuple are usually filenames.
            subRunID (str, int or tuple):
                An ID to store a part of a run under.  In our case
                this is usaully a cover length.
            seedBits (int):
                The number of bits to get from the randomness source.

        Yields:
            (tuple):A tuple containing a seed and its owner in the
            format (owner, seed).
        '''
        seeds = {}
        for owner in owners:
            seed = HQrandom.getrandbits(seedBits)
            seeds[owner] = seed
            yield (owner, seed)
        if subRunID is None:
            self.bank['data'][self.runID] = seeds
        else:
            self.bank['data'][self.runID][subRunID] = seeds
    
    def getMetaData(self, attribute=None):
        '''
        Get the metadata from this instance, if attrubute is provided
        only that attrbute of the metadata will be returned.

        Args:
            attribute (str):
                The name of the attrubute you want.
        '''
        if attribute is None:
            return self.bank['metadata']
        return self.bank['metadata'][attribute]
        

def printProgress(completed, count, prefix='', size=40):
    '''
    Print a progress bar and reset the cursor to he start of the
    line.

    Args:
        completed (int):
            The number of items completed
        count (int):
            The total number of items to process.
        prefix (str):
            A message to print before the progress bar
        size (int):
            The number of character to make the progress bar long.
    
    Returns:
        (int):The number of characters printed.

    '''
    x = int(size*completed/count)
    outStr = "%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), completed, count)
    sys.stdout.write(outStr)
    sys.stdout.flush()
    return len(outStr)


def progressBarIterable(lst, prefix="", size=40, count=None):
    '''
    Create a generator that updates progress bar whenever it yeilds
    an item.  Once the generator is exausted it will print whitespace
    and reset the cursor to the start of the line to clean up after
    itself.

    Args:
        lst (iterable):
            An iterable containing the values to be looped over.
        count (int):
            The total number of items to process. This is only
            required if the iterable provided doesn't support len.
        prefix (str):
            A message to print before the progress bar
        size (int):
            The number of character to make the progress bar long.
    
    Raises:
        TypeError:
            When lst is a generator and a count is not provided.

    Yeilds:
        (object):Elements of the list you provided

    '''
    if count is None:
        if not hasattr(lst, '__len__'):
            raise TypeError('If lst is a generator count must be supplied.')
        count = len(lst)
    progressLen = printProgress(0, count, prefix, size)
    for i, item in enumerate(lst):
        yield item
        progressLen = printProgress(i+1, count, prefix, size)
        sys.stdout.flush()
    sys.stdout.write((' ' * progressLen) + '\r')


def equalAreaBins(data, binCount = 5):
    '''
    Return binCount bins with equal area under the curve.
    
    Args:
        data (iterable):
            The data you would like binned
        binCount (int):
            The number of bins you would like to use

    Returns:
        (tuple):The boundries of the bins.  (Left endpoints)
    '''
    percentiles=[100/binCount*i for i in range(binCount)]
    bins = np.percentile(data, percentiles)
    return list(bins)


def chunks(l, n):
    """
    Yield successive `n`-sized chunks from `l`.

    Args:
        l (iterable):
            A list or tuple you want to make chunks out of.
        n (int):
            The size of the chunks made from `l`

    Yields:
        (iterable):Slices of the orignal list of size `n`
    """
    for i in range(0, len(l), n):
        yield islice(l, i, i + n)


def setArgsFromConfig(cfgPath, args, section='general'):
    '''
    Set config options that are in the config file but aren't set
    by the current commad line arguments.

    Note:
        Don't use config files from un-trusted sources as this
        function uses `eval` as a hacky cast.

    Args:
        cfgPath (str): 
            The path to the input config file.
        args (argparse.Namespace):
            The argprase namespace you want to update with values
            from the config file.

    Return:
        (argparse.Namespace):The argparse namespace containing the data
        from the config file and the command line.

    '''
    cfg = configparser.ConfigParser()
    # Override the optionxform method of RawConfigParser so it does
    # not conert our keys to lowercase.
    cfg.optionxform = str
    cfg.read(cfgPath)
    for cfgKey, cfgVal in cfg[section].items():
        if not hasattr(args, cfgKey) or getattr(args, cfgKey) is None:
            setattr(args, cfgKey, eval(cfgVal))
    return args