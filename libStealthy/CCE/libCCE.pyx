#cython: embedsignature=True
# this import should be refactored and split up into seprate files
from libCCE cimport llRoot, createLinkedList, insertSequence, freeLL, tNode, createTree, printLayer, printList, calcCCEs, applyToList, llNode, isnan, resetArena
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.stdlib cimport free, malloc
from libc.stdio cimport printf
from collections import Counter
from scipy.stats import entropy
from itertools import islice
import numpy as np
import cython


@cython.embedsignature(True)
cdef class CCE(object):
    """
    Calculate CCE of a given sequence.  The inital body of samples
    can be provided to `__init__` and they will be populated into the
    tree using a faster pure C method.  Any new data can be provided
    to `insertSeqence` which will insert it into the tree.

    This class is a wrapper for the cCCE implementation which is based
    on the work done by S. Gianvecchio and H. Wang's in "Detecting
    Covert Timing Channels: An Entropy-Based Approach" section 3.4.

    As a general rule this implementation will try and fail via 
    assertion before returing compromised information to the user.

    Args:
        branchingFactor (int): 
            The branching factor of the tree, this is equivalent to
            the number of bins you have.
        initalSequences (iterable):
            A python iterable which contains IPDs for insertion.
        subSeqLen (int): 
            The length you would like the subsequences devided into.
    Raises:
        ValueError:
            if subSeqLen is less than the length of the original 
            sequence or if the branching factor is less than or equal
            to zero.
    """
    # instance level c attrubutes 
    cdef tRoot* root
    cdef unsigned int branchingFactor
    
    def __cinit__(self, unsigned int branchingFactor=5, initalSequences=None, int subSeqLen=50):
        cdef int seqLen 
        cdef int* cSeq
        cdef int i = 0
        cdef int offset

        if branchingFactor <= 0:
            raise ValueError('branching factor must be > 0.')

        self.root = createTree(branchingFactor)
        self.branchingFactor = branchingFactor

        if initalSequences is not None:
            # allocate memory to pass the values from the python
            # object to the CCE code
            seqLen = len(initalSequences)
            if subSeqLen > seqLen:
                raise ValueError("subSeqLen must be less than the length of the initial sequance")
            cSeq = <int*>PyMem_Malloc(seqLen * sizeof(int))
            
            # Copy the python object to the memory
            for seq in initalSequences:
                cSeq[i] = seq
                # try our darndest to prevent foot shooting via heap 
                # based buffer overflows in the branching array.
                if cSeq[i] >= branchingFactor:
                    raise ValueError("Bin value %d at index %d is not in range [0, %d)" % (seq, i, self.branchingFactor))
                i += 1
            
            # insert sequences into the tree using pointer math to
            # set the section we are intrested in.
            for i in range((seqLen - subSeqLen)):
                insertSequence(self.root, cSeq + i, subSeqLen)
            
            PyMem_Free(cSeq)


    def __init__(self, branchingFactor=5, initalSequence=None, subSeqLen=50):
        # All real work is done in __cinit__
        pass

    def __dealloc__(self):
        """ 
        Used to cleanup memory on destruction of object. Automaticlly
        called should NOT be manually called it WILL break everything.
        """
        freeTree(self.root)

    def insertSequence(self, pySeq):
        """
        Wrapper for insertSequence, converts a list or tuple into an
        intiger array of the same length.

        Args:
            pySeq (list or tuple):
                A python iterable that supports indexing and len()
        """
        cdef int* cSeq
        cdef int i
        cdef int seqLen

        seqLen = len(pySeq)
        cSeq = <int*>PyMem_Malloc(seqLen * sizeof(int))
        # populate the C buffer with values from the python object.
        for i in range(seqLen):
            cSeq[i] = pySeq[i]
            # try our darndest to prevent foot shooting via heap 
            # based buffer overflows
            if cSeq[i] >= self.branchingFactor:
                raise ValueError("Bin value %d at index %d is not in range [0, %d)" % (cSeq[i], i, self.branchingFactor))
        insertSequence(self.root, cSeq, seqLen)
        PyMem_Free(cSeq)

    def calculateCCE(self):
        """
        Calculate the CCE for all inserted sequences and return the
        CCE for all sub-sequence lengths upto the maximum sub-
        sequence length.

        Returns:
            (list): A list of floats containing the CCE for sub-
            sequence length.
        """
        cdef double* CCEArr
        cdef llNode* n
        cdef int i = 0

        CCEs = []

        CCEArr = calcCCEs(self.root)
        
        # The array returned from calcCCEs is NaN terminated so we
        # loop till we find the NaN.
        while not isnan(CCEArr[i]):
            CCEs.append(CCEArr[i])
            i += 1

        # free memory alocated in calcCCEs to hold the return values.
        free(CCEArr)

        return CCEs

    def resetTree(self):
        '''
        This funciton can be used to reset the tree and reuse it 
        of freeing all the memory assocated with the tree.

        This can be super handy when you want to calculate new data
        that is a of a similar size to the previously analized data.
        '''
        resetArena(self.root.arenaManagement)


def _window(seq, int n=2):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

cdef double _precentUnique(frequencys):
    f = [x for x in frequencys if x == 1] 
    return <double>(len(f) / sum(frequencys))

@cython.embedsignature(True)
def localMinCCE(seq, int maxLen):
    '''
    Calculate the CCE of a sequence using the local minimum as an
    early stoping metric.

    Args:
        seq (iterable):
            An iterable containing the pre-bined values of your IPD
            sequence.
        maxLen (int):
            The maximum subsequence length to be tried.  This is
            more of an upper bound as the local minimum is usually
            found before this is reached.
    
    Returns:
        (float):The local minimum CCE of the provided sequence.
    '''
    cdef int seqLen 
    cdef double firstOrderEntropy
    cdef double correctionFactor
    cdef double previousEntropy
    cdef double uniques
    cdef double minCCE
    cdef double CCE
    cdef double ent

    previousEntropy = 0.0
    firstOrderEntropy = 0.0
    minCCE = float('inf')
    # print('CCE, Conditional Entropy, Correction Factor') #debug info
    
    # Prime the first order entropy.
    c = Counter(_window(seq[:len(seq)-maxLen], 1))
    frequency = list(c.values())
    ent = entropy(frequency)
    firstOrderEntropy = ent

    for seqLen in range(2, maxLen+1):
        conditionalEntroy = ent - previousEntropy
        previousEntropy = ent
        uniques = _precentUnique(frequency)
        correctionFactor = firstOrderEntropy * uniques
        CCE = conditionalEntroy + correctionFactor
        # print("%f, %f, %f" % (CCE, conditionalEntroy, correctionFactor)) # debug info
        if CCE > minCCE:
            break
        else:
            minCCE = CCE
        if uniques == 1.0: 
            break
        # This funky slicing ensures that the output of the tree
        # based version and this one match.  You lose out a few
        # samples per length, but they would have been from an
        # imcompleate segment anyway.
        c = Counter(_window(seq[:len(seq)-(maxLen-(seqLen-1))], seqLen))
        frequency = list(c.values())
        # print(frequency)
        ent = entropy(frequency)

    return minCCE