#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from math import modf, sqrt
from scipy.stats import ttest_ind, mannwhitneyu
from os.path import isfile, isdir
from os import makedirs
from functools import lru_cache

from libStealthy.stego.fromIPDs import SCTimeReplayFromIPDsUE
from libStealthy.exceptions import EncodingError

import numpy as np
import warnings
import random
import json


@lru_cache()
def calcProbLow(refIPDs):
    ''' 
    Go through refIPDs and determine probLow
    That will be number of below-median values
    divided by the adjusted (dropping median values)
    length of refIPDs
    '''
    count = 0
    adjN = 0
    med = np.median(refIPDs)
    for ipd in refIPDs:
        if ipd < med:
            count += 1
            adjN += 1
        elif ipd > med:
            adjN += 1
    probLow = count/float(adjN)
    return probLow


def srlFillerEmbedding(beta, p, coverMsg):
    '''
    Flip each bit of coverMsg independently with probability beta * p.

    Take cover source coverMsg, a 0-1 list, and flip each bit independently with probability beta * p.  Return a new 0-1 list.
    '''

    # Validate inputs
    if p < 0.0 or p > 1.0:
        raise RuntimeError('p must be in [0,1]!')
    if beta < 0.0 or beta > 1.0:
        raise RuntimeError('beta must be in [0,1]!')


    ###

    outList = []
    flipProb = beta * p

    for bit in coverMsg:
        if (bit != 0) and (bit != 1):
            raise RuntimeError('coverMsg is not a 0-1 iterable!')

        flipRand = np.random.random()
        if flipRand < flipProb:
            outList.append(1 - bit)
        else:
            outList.append(bit)

    return outList


def skewedTREmbedding(beta, p, q, coverLen):
    '''
    Generate random 0-1 list with certain parameters.

    This produces a 0-1 list of length coverLen.  The probability that an element of this list is 0 equals q + beta * p - 2 * q * beta * p, and the probability that an element of this list is 1 equals 1 - q - beta * p + 2 * q * beta * p.  beta, p, and q must be in [0,1].
    '''

    ### Validate inputs

    ## Do we want to check abs(q - 0.5) < eps instead?
    if q == 0.5:
        raise RuntimeError('q must not be 0.5!')
    elif abs(q - 0.5) < 0.001:
        warnings.warn('q is very close to (within 0.001 of) 0.5!')
        # Produce a warning that we're close to 0.5

    if q < 0.0 or q > 1.0:
        raise RuntimeError('q must be in [0,1]!')
    if p < 0.0 or p > 1.0:
        raise RuntimeError('p must be in [0,1]!')
    if beta < 0.0 or beta > 1.0:
        raise RuntimeError('beta must be in [0,1]!')

    i = 0
    outList = []
    zeroProb = q + (beta * p) - (2 * q * beta * p)
    while i < coverLen:
        r = np.random.random()

        if r < zeroProb:
            outList.append(0)
            i += 1
        else:
            outList.append(1)
            i += 1

    return outList


def computeNu0Binary(seq, prob0):
    '''
    Compute nu for the input binary sequence given probability of a 0 in the input.

    Note that nu is the same whether 0 or 1 is used (for two-letter case).
    '''

    n = len(seq)

    count = 0
    for bit in seq:
        if bit == 0:
            count += 1
        elif bit != 1:
            raise RuntimeError('seq is not a 0-1 iterable!')

    nu = float(np.sqrt(n) * np.absolute((count / float(n)) - prob0))

    return nu


def detectorNu0Binary(seq, PFAstar, prob0):
    '''
    Detector following notes (given threshold and probability of 0)

    For the given 0-1 sequence seq, evaluate the test statistic nu (given the probability prob0 of a 0 in the sequence-generation).  If this is greater than the threshold T (determined as a function of PFAstar), return True (seq is a stego. sequence).  Otherwise, return False (seq is not a stego. sequence).
    '''

    # Value of C from notes
    C = 41

    # Threshold based on C and P_{FA}^*
    T = np.sqrt(C/float(PFAstar))

    nu = computeNu0Binary(seq, prob0)

    if nu > T:
        return True
    else:
        return False


def statisticNu0Binary(seq, prob0):
    '''
    Return statistic nu following notes (given threshold and probability of 0)

    For the given 0-1 sequence seq, compute the test statistic nu (given the probability prob0 of a 0 in the sequence-generation).
    '''

    nu = computeNu0Binary(seq, prob0)

    return nu


def computeNu0IPD(seq, probLow, signed = False):
    '''
    Compute nu for the input IPD sequence given probability of a below-median value in the input.
    '''

    med = float(np.median(seq))

    count = 0
    n = 0

    for ipd in seq:
        if ipd < med:
            count += 1
            n += 1
        elif ipd > med:
            n += 1


    nu = float(np.sqrt(n) * np.absolute((count / float(n)) - probLow))

    return nu


def computeSignedNu0IPD(seq, probLow):
    '''
    Compute signed version of nu for the input IPD sequence given probability of a below-median value in the input.

    This version checks whether IPDs are *strictly* below the median.
    '''

    med = np.median(seq)

    count = 0
    n = 0

    for ipd in seq:
        if ipd < med:
            count += 1
            n += 1
        elif ipd > med:
            n += 1


    sgnu = float(np.sqrt(n) * (count / float(n) - probLow))

    return sgnu


def detectorNu0IPD(seq, PFAstar, probLow):
    '''
    Detector following notes (given threshold and probability of 0)

    For the given IPD sequence seq, evaluate the test statistic nu (given the probability probLow of a below-median value in the sequence-generation distribution).  If this is greater than the threshold T (determined as a function of PFAstar), return True (seq is a stego. sequence).  Otherwise, return False (seq is not a stego. sequence).
    '''

    # Value of C from notes
    C = 41

    # Threshold based on C and P_{FA}^*
    T = np.sqrt(C/float(PFAstar))

    nu = computeNu0IPD(seq, probLow)

    if nu > T:
        return True
    else:
        return False



def statisticNu0IPD(seq, probLow, signed=False):
    '''
    Return statistic nu following notes (given threshold and probability of below-median IPD)

    For the given IPD sequence seq, compute the test statistic nu (given the probability probLow of a below-median value in the sequence-generation).
    '''

    nu = computeNu0IPD(seq, probLow, signed)

    return nu


def statisticSignedNu0IPD(seq, probLow):
    '''
    Return signed version of statistic nu following notes (given threshold and probability of below-median IPD)

    For the given IPD sequence seq, compute signed version the test statistic nu (given the probability probLow of a below-median value in the sequence-generation).

    This version tests for *strictly* below-median values.
    '''

    sgnu = computeSignedNu0IPD(seq, probLow)

    return sgnu

