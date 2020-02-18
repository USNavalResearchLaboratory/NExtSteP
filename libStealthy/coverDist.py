#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
'''
Classes (BaseCoverDist and its subclasses) for producing a
sequences of IPDs for cover traffic. These all have the
generateIPDseq method, which returns a list of IPDs.  In some
cases, there is also an IPDdist method, which returns the
distribution function captured by the object.
'''

from scipy.stats import expon, halfnorm, foldnorm, norm, ttest_ind, lomax
from markovify.chain import accumulate, BEGIN
from Crypto.Random import random as HQrandom
from abc import ABCMeta, abstractmethod
from markovify import Chain

from libStealthy.exceptions import NoDistributionError, EncodingError, DecodingError
from libStealthy.util import window, equalAreaBins, chunks
from libStealthy.pRandom import PRNG

import numpy as np
import argparse
import bisect
import sys

class BaseCoverDist(object, metaclass=ABCMeta):
    '''
    The abstract base class describing the interface for all cover
    distribution classes.
    
    '''
    distStr = 'ABC'
    @abstractmethod
    def generateIPDseq(self, numt):
        '''
        Generate a sequence according to the distribution captured by this
        object.
        
        Args:
            numt (int):
                The number of IPDs you would like to generate.
        
        Returns:
            (list):A list of numbers representing the IPDs.
        '''
        pass

    @abstractmethod
    def IPDdist(self):
        '''
        Distribution captured by this object. Calling the returned
        method will return a value distributed according to the
        distribution.

        Returns:
            (function):A function that when called returns a value
            from the distrubuiton.
        '''
        pass

    @abstractmethod
    def getMetaData(self):
        '''
        Get distribution-specific metadata.

        Returns:
            (dict):A dict containing distribution-specific metadata.
        '''
        pass

    @staticmethod
    @abstractmethod
    def generateSubparser(subparsers, alias = 'base'):
        '''
        Create a subparser with distribution-specific arguments.

        Args:
            subparsers (argparse._SubParsersAction):
                An argparse SubParsersAction object to which we will
                attach our subparser.
            alias (str):
                The string to be used invoke the subparser from the
                command line.
        Returns:
            (argparse.ArgumentParser):The handle to the subparser
            created.

        '''
        pass

class Exponential(BaseCoverDist):
    '''
    Cover distribution in which IPDs follow exponential distribution
    with parameter lam.
    
    Args:
        lam (int or float):
            lam... whatever that is...
    '''
    distStr = 'Exponential'
    def __init__(self, lam):
        self.lam = lam
        self.dist = expon(scale=1./self.lam).rvs

    def generateIPDseq(self, numt):
        IPDseq = [self.dist() for i in range(numt)]
        return IPDseq

    def IPDdist(self):
        return self.dist
    
    def getMetaData(self):
        return {
            'lam':self.lam
            ,'coverType':self.distStr
        }

    @staticmethod
    def generateSubparser(subparsers, alais = 'exponential'):
        sp = subparsers.add_parser(alais, help='Cover distribution in which IPDs follow exponential distribution.')
        sp.add_argument('--lam', type=float)
        sp.set_defaults(parserID='exponential')
        return sp

# # Preliminary implementation; not currently used.
#
# class FoldNormInt(BaseCoverDist):
#     '''
#     Cover distribution IPDs follow a folded-normal distribution
#     '''
#     distStr = 'FoldNormInt'
#     def __init__(self, c, loc=0, scale=1):
#         self.c = c
#         self.loc = loc
#         self.scale = scale
#         self.dist = foldnorm(self.c, loc=self.loc, scale=self.scale).rvs
# 
#     def generateIPDseq(self, numt):
#         IPDseq = [int(self.dist()) for i in range(numt)]
#         return IPDseq
# 
#     def IPDdist(self):
#         return lambda : int(self.dist())
#     
#     def getMetaData(self):
#         return {
#             'c':self.c
#             ,'loc':self.loc
#             ,'scale':self.scale
#             ,'coverType':self.distStr
#         }

class BasePRNGCover(object):
    '''
    A base class to handle common patterns related to the use of PRNG

    Args:
        PRNGSource (str):
            A string with the name of a PRNG source.
        seed (bytes or int):
            A seed for the PRNG.  By default this uses a random one.
    '''
    def __init__(self, PRNGSource = 'chacha20', seed = None):
        self.seed = seed
        self.PRNGSource = PRNGSource
        self.PRNGClass = PRNG.getRng(PRNGSource)

    def prepRNG(self, seed = None):
        '''
        Prepare an instance of the PRNG class specfied during creation.

        Args:
            seed (int or bytes):
                A seed value for the new PRNG to use.  If not supplied
                the one supplied at creation will be used, if one was not
                supplied at creation one will be generated.

        Returns:
            An instance of the PRNG class specified by PRNGSource.
        '''
        if seed is None:
            if self.seed is None:
                self.seed = HQrandom.getrandbits(self.PRNGClass.seedSize)
            seed = self.seed
        return self.PRNGClass(seed)


class Normal(BasePRNGCover, BaseCoverDist):
    '''
    Cover distribution that follows a (truncated) normal distribution

    Args:
        mean (int):
            The mean value for the distribution
        sd (float):
            The standard deviation for the distrubuiton
        PRNGsorce (str):
            The name of a PRNG source to use.
        seed (bytes or int):
            A seed for the PRNG.  By default this uses a random one.
        minIPD (int):
            The maximum delay between pakets allowed.
        maxIPD (int):
            The maximum delay between pakets allowed.

    '''
    distStr = 'Normal'
    def __init__(self, mean, sd, PRNGSource = 'chacha20', seed = None, minIPD =1000, maxIPD=5000):
        # minIPD and maxIPD bound the values that are returned
        # These incorporate the effects of resFactor in generateIPDseq
        super().__init__(PRNGSource, seed)
        self.mean = mean
        self.sd = sd
        self.minIPD = minIPD
        self.maxIPD = maxIPD


    def generateIPDseq(self, numt, seed = None):
        ra = self.prepRNG(seed)
        resFactor = 10 # 1 gives ms, 10 gives 10^{-4} s

        rv = norm(loc=self.mean, scale=self.sd)

        F = rv.cdf
        Finv = rv.ppf

        normals = []
        while len(normals) < numt:
            r = ra.random()
            ipd = int(resFactor*Finv(r))
            if self.minIPD <= ipd <= self.maxIPD:
                normals.append(ipd)

        return normals
    
    def IPDdist(self, seed = None):
        return lambda: self.generateIPDseq(1, seed)[0]

    def getMetaData(self):
        return {
            'seed':self.seed
            ,'mean':self.mean
            ,'sd':self.sd
            ,'PRNG':self.PRNGSource
            ,'coverType':self.distStr
        }

    @staticmethod
    def generateSubparser(subparsers, alias='normal'):
        sp = subparsers.add_parser(alias, help='Use a user-specfied normal distribution.  The mean and standard deviation are set using --mean and --sd respectively.', description='Use a user specfied Pareto distribution.')
        sp.add_argument('--mean', metavar='MU', type=int, help='Mean value for normal distribution')
        sp.add_argument('--sd', metavar='SIGMA', type=float, help='Standard deviation for normal distribution')
        sp.set_defaults(parserID = 'normal')
        return sp

class Pareto(BasePRNGCover, BaseCoverDist):
    '''
    Cover distribution that follows a Pareto distribution

    Args:
        alpha (float):
            The alpha value for the distribution (lower bound of support)
        beta (float):
            The beta value for the distrubuiton
        PRNGsorce (str):
            The name of a PRNG source to use.
        seed (bytes or int):
            A seed for the PRNG.  By default this uses a random one.
        maxIPD (int):
            The maximum delay between pakets allowed.

    '''
    distStr = 'Pareto'
    def __init__(self, alpha, beta, PRNGSource = 'chacha20', seed = None, maxIPD=500):
        super().__init__(PRNGSource, seed)
        self.alpha = alpha
        self.beta = beta
        self.maxIPD = maxIPD

    def generateIPDseq(self, numt, seed = None):
        ra = self.prepRNG(seed)
        resFactor = 10 # 1 gives ms, 10 gives 10^{-4} s

        F = lambda x : 1 - (float(self.alpha)/float(x))**self.beta
        maxRand = F(self.maxIPD + 1.0/resFactor)

        Finv = lambda x : self.alpha * ((1 / (1.0 - x))**(1.0/self.beta))

        rands = [ra.random() for i in range(numt) if i <= maxRand]

        while len(rands) < numt:
            r = ra.random()
            if r <= maxRand:
                rands.append(r)
        paretos = [int(resFactor*Finv(x)) for x in rands]
        return paretos
    
    def IPDdist(self, seed = None):
        return lambda: self.generateIPDseq(1, seed)[0]

    def getMetaData(self):
        return {
            'seed':self.seed
            ,'alpha':self.alpha
            ,'beta':self.beta
            ,'PRNG':self.PRNGSource
            ,'coverType':self.distStr
        }

    @staticmethod
    def generateSubparser(subparsers, alais='pareto'):
        sp = subparsers.add_parser(alais, help='Use a user specfied Pareto distribution.  alpha and beta are set using --paretoa, and --paretob respectively.', description='Use a user specfied Pareto distribution.')
        sp.add_argument('--paretoa', metavar='ALPHA', type=int, help='Alpha value (minimum value possible) for Pareto distribution')
        # Process beta as a string to facilitate standardization of floating-point argument
        sp.add_argument('--paretob', metavar='BETA', type=float, help='Beta value (shape) for Pareto distribution')
        sp.set_defaults(parserID = 'pareto')
        return sp

class ParetoNarrow(Pareto):
    '''
    Cover distribution that follows a Pareto distribution 
    where alpha=100 beta=10.

    Args:
        PRNGsorce (str):
            The name of a PRNG source to use.
        seed (bytes or int):
            A seed for the PRNG.  By default this uses a random one.
        maxIPD (int):
            The maximum delay between pakets allowed.

    '''
    distStr = 'ParetoNarrow'
    def __init__(self, PRNGSource = 'chacha20', seed = None, maxIPD=500):
        super().__init__(100, 10, PRNGSource, seed, maxIPD)

    @staticmethod
    def generateSubparser(subparsers, alais = 'narrow'):
        sp = subparsers.add_parser(alais, help='Use a pre-defined Pareto distribution of alpha=100 beta=10', description='Use a pre-defined Pareto distribution of alpha=100 beta=10')
        sp.set_defaults(parserID = 'narrow')
        return sp

class ParetoSellke(Pareto):
    '''
    Cover distribution that follows a Pareto distribution 
    where alpha=100 beta=0.95.

    Args:
        PRNGsorce (str):
            The name of a PRNG source to use.
        seed (bytes or int):
            A seed for the PRNG.  By default this uses a random one.
        maxIPD (int):
            The maximum delay between pakets allowed.

    '''
    distStr = 'ParetoSellke'
    def __init__(self, PRNGSource = 'chacha20', seed = None, maxIPD=500):
        super().__init__(100, 0.95, PRNGSource, seed, maxIPD)

    @staticmethod
    def generateSubparser(subparsers, alias = 'sellke'):
        sp = subparsers.add_parser(alias, help='Use a pre-defined Pareto distribution of alpha=100 beta=0.95', description='Use a pre-defined Pareto distribution of alpha=100 beta=0.95')
        sp.set_defaults(parserID='sellke')
        return sp


class hqChain(Chain):
    def __init__(self, corpus, state_size, model=None, PRNGSource = None, seed = None):
        """
        `corpus`: A list of lists, where each outer list is a "run"
        of the process (e.g., a single sentence), and each inner list
        contains the steps (e.g., words) in the run. If you want to simulate
        an infinite process, you can come very close by passing just one, very
        long run.
        `state_size`: An integer indicating the number of items the model
        uses to represent its state. For text generation, 2 or 3 are typical.
        """
        
        # setup PRNG
        if PRNGSource is None:
            self.prng = PRNG.ChaCha20(seed)
        else:
            self.prng = PRNGSource(seed)
        # Create markov chain
        super().__init__(corpus, state_size, model)

    def move(self, state):
        """
        Given a state, choose the next item at random.
        """
        if state == tuple([ BEGIN ] * self.state_size):
            choices = self.begin_choices
            cumdist = self.begin_cumdist
        else:
            choices, weights = zip(*self.model[state].items())
            cumdist = list(accumulate(weights))
        r = self.prng.random() * cumdist[-1]
        selection = choices[bisect.bisect(cumdist, r)]
        return selection

class MarkovIngest(BasePRNGCover, BaseCoverDist):
    '''
    Create a markov based traffic IPD generator from a corpus of input data. 
    
    Args:
        corpus (list):
            A list of values to bulid the model from.
        state_size (int):
            The number of output types you want.
        lookbackSize (int):
            The lengths of the runs to use when building the
            chain.
        bins (list):
            A list containing the edges of the bins used to
            represent diffrent possible output ranges. By default
            bins of equal area are calculated automatically.
        disjoint (bool):
            Use dsijoint runs when building the model.
        PRNGsorce (str):
            A string with the name of a PRNG source.
        seed (bytes or int):
            A seed for the PRNG.  By default this uses a random one.

    '''
    distStr = 'MarkovIngest'
    
    def __init__(self, corpus, state_size, lookbackSize = 2, bins = None, disjoint = False, PRNGSource = 'chacha20', seed = None):
        super().__init__(PRNGSource, seed)
        
        # setup bins if we aren't provided any.
        if bins is None:
            bins = equalAreaBins(corpus, state_size)

        # Convert to python ints from np.int64
        bindedCorpus = [int(digit) for digit in np.digitize(corpus, bins) - 1 ]

        # Setup IPDs to draw from when when retruning values
        self.knownIPDs = {bin_:[] for bin_ in range(len(bins))}
        for digit, ipd in zip(bindedCorpus, corpus):
            self.knownIPDs[digit].append(ipd)

        # break corpus into runs.
        if disjoint:
            runs = chunks(bindedCorpus, lookbackSize)
        else:
            runs = window(bindedCorpus, lookbackSize)
        
        # filter out short runs and convert tuples to lists
        bindedCorpus = [list(run) for run in runs if len(run) == lookbackSize]
        
        # create the chain
        self.chain = hqChain(bindedCorpus, state_size, PRNGSource=self.PRNGClass, seed=seed)

    def generateIPDseq(self, numt, seed = None):
        cover = []
        ra = self.prepRNG(seed)
        self.chain.prng = ra
        while len(cover) < numt:
            # Walk the chain and convert the resulting bins to known IPDs
            ipds = [ra.choice(self.knownIPDs[bin_]) for bin_ in self.chain.walk()]
            cover.extend(ipds)
        cover = cover[:numt]
        return cover

    def IPDdist(self, seed = None):
        def IPDGen():
            while True:
                for ipd in self.generateIPDseq(1000, seed):
                    yield ipd
        ipds = IPDGen()
        return lambda: next(ipds)

    def getMetaData(self):
        return {
            'seed':self.seed
            ,'model':self.chain.to_json()
            ,'PRNG':self.PRNGSource
            ,'coverType':self.distStr
        }
    
    @staticmethod
    def generateSubparser(subparsers, alias='markov-ingest'):
        sp = subparsers.add_parser(alias, help='A cover generator based on a Markov Chain')
        sp.add_argument('--corpusPath', type=str, nargs='+', help='The path to the corpus of of the traffic whose distribution you want to approximate.')
        sp.add_argument('--lookbackSize', type=int, default=2, help='The number of IPDs to lookback to generate runs for the model.')
        sp.add_argument('--stateSize', type=int, default=2, help='The size of the interal state the Markov Chain will use, default is 2')
        sp.add_argument('--bins', type=float, nargs='*', default=None, help='Specify the specific boundaries of the bins.')
        sp.add_argument('--disjointLookback', action='store_true', help='Use disjoint (no overlap) segments from the orignial cover.')
        sp.set_defaults(parserID='markov-ingest')
        return sp

# Here be introspective dragons...
# We use some of Python's inrospection features to build
# references to classes from this module's internal dict.
def getCoverDist(coverDist):
    '''
    Get a reference to a CoverDist.

    Args:
        rngName (str):
            A string containing the name of the CoverDist you want.
    Returns:
        A reference to the constructor of the requested class
    
    Raises:
        NotImplementedError:
            Raised if the provided string doesn't correspond to
            anything in this module
    '''
    # get a ref to the current module (CoverDists)
    currentModule = sys.modules[__name__]
    # check to see if the current module has anything matchng the
    # supplied name.
    if hasattr(currentModule, coverDist):
        # return a ref to the object if so
        return getattr(currentModule, coverDist)
    else:
        # raise NotImplmentedError if we don't have that CoverDist
        raise NotImplementedError('CoverDist %d not implmented.' % coverDist)


def getCoverDists():
    '''
    Get a tuple of the names of the CoverDists implmented by this
    module.

    Returns:
        (tuple):A tuple of the names of the CoverDists implmented by
        this module.
    '''
    # get a ref to the current module (coverDist)
    currentModule = sys.modules[__name__]
    coverDists = []
    # create a generator to get memebers from the module by name.
    moduleMembers = ((getattr(currentModule, name), name) for name in dir(currentModule))
    for member, name in moduleMembers:
        # remove things that aren't classes
        if isinstance(member, type):
            # remove things that aren't subclassed from BaseCoverDist
            # like BasePRNGCover, and the markov subclasses.
            if issubclass(member, BaseCoverDist):
                coverDists.append(name)
    return tuple(coverDists)
