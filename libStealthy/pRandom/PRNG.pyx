#cython: embedsignature=True
#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from SFMT cimport sfmt_genrand_uint64, sfmt_genrand_uint32, sfmt_t, sfmt_init_gen_rand
from PCG cimport pcg32_random_t, pcg32_srandom_r, pcg32_random_r, pcg32_boundedrand_r
from xoshiro cimport next as next_
from libc.stdint cimport uint32_t, uint64_t, uint8_t, uint16_t
from libc.stdio cimport printf
from random import Random
from Crypto.Random import random as HQrandom
from Crypto.Cipher import ChaCha20 as chacha20
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
import sys


cdef class cpcg:
    '''
    Python cdef class wrapper for C functions relted to the PCG PRNG.

    This class should not be instantiated directly from python, only from
    the related pure python wrapper class. 
    '''
    cdef pcg32_random_t state
    cdef uint64_t seedVal
    def __cinit__(self, uint64_t seedVal, uint64_t rounds = 5):
        pcg32_srandom_r(&self.state, seedVal, rounds)
    def __init__(self, seedVal, rounds=5):
        pass
    def seed(self, uint64_t seedVal, uint64_t rounds = 5):
        pcg32_srandom_r(&self.state, seedVal, rounds)
    def randint8(self):
        return pcg32_random_r(&self.state) >> 24
    def randint16(self):
        return pcg32_random_r(&self.state) >> 16
    def randint32(self):
        return pcg32_random_r(&self.state)

cdef class csfmt:
    '''
    Python cdef class wrapper for C functions relted to the SFMT PRNG.

    This class should not be instantiated directly from python, only from
    the related pure python wrapper class. 
    '''
    cdef sfmt_t state
    cdef uint32_t seedVal
    def __cinit__(self, uint32_t seedVal):
        sfmt_init_gen_rand(&self.state, seedVal)
        self.seedVal = seedVal
    def __init__(self, seedVal):
        pass
    def seed(self, uint32_t seedVal):
        sfmt_init_gen_rand(&self.state, seedVal)
    def randint8(self):
        return sfmt_genrand_uint32(&self.state) >> 24
    def randint16(self):
        return sfmt_genrand_uint32(&self.state) >> 16
    def randint32(self):
        return sfmt_genrand_uint32(&self.state)
    def randint64(self):
        return sfmt_genrand_uint64(&self.state)

cdef class cxoshiro:
    '''
    Python cdef class wrapper for C functions related to the xoshiro PRNG.

    This class should not be instantiated directly from python, only from
    the related pure python wrapper class. 
    '''
    cdef uint8_t* state
    cdef uint8_t* seedVal
    def __cinit__(self, seedVal):
        cdef int i
        self.state = <uint8_t*>PyMem_Malloc(32 * sizeof(uint8_t))
        self.seedVal = <uint8_t*>PyMem_Malloc(32 * sizeof(uint8_t))
        for i in range(32):
            self.state[i] = 0
            self.seedVal[i] = 0
    def ___init__(self, seedVal):
        self.seed(seedVal)

    def __dealloc__(self):
        PyMem_Free(self.state)
        PyMem_Free(self.seedVal)

    def seed(self, seedVal):
        cdef int i
        if len(seedVal) != 32:
            raise ValueError('Seed must be exactly 32 bytes (256 bits) long.')
        for i in range(32):
            self.state[i] = seedVal[i]
            self.seedVal[i] = seedVal[i]
        
    def randint64(self):
        return next_(<uint64_t*>self.state)
    def randint32(self):
        return self.randint64() >> 32
    def randint16(self):
        return self.randint64() >> 48
    def randint8(self):
        return self.randint64() >> 56

class cRandom(Random):
    '''
    Base class for pure Python wrapper classes, __init__ is meant to be
    implmented by the classes which subclass it.
    '''
    def __init__(self, seed):
        super().__init__(seed)

    # version is unused and only exists to maintain compatiblity
    # with random API
    def seed(self, seed, version = 2):
        '''
        Args:
            seed (int or bytes):
                A 32 bit value used to seed the state of the PRNG.
            version (int):
                This is unused nd only exists to maintain
                compatablity with the python Random API
        '''
        if isinstance(seed, bytes):
            seed = int.from_bytes(seed, 'little')
        # truncate seed to 64 bits
        seed = 0xFFFFFFFF & seed
        self.seedVal = seed
        self.cRand.seed(seed)

    def random(self):
        '''
        Get a random float, the float has 64 bits of randomness

        Returns:
            (float):A number between 0 and 1
        '''
        return self.cRand.randint64() / ((2 ** 64) - 1)

    def getrandbits(self, int bits):
        '''
        Returns an int composed of a user specfied number of random
        bits. This is implmented so that all features of the Random
        parent class will funciton using this PRNG.

        Args:
            bits(int):
                The number of bits you would like form the PRNG
        Returns:
            (int):An int composed of random bits.
        '''
        cdef int i
        cdef int numBytes
        cdef uint64_t randVal
        cdef uint8_t* bitBuffer

        if bits <= 0:
            raise ValueError('number of bits must be greater than 0')
        # if we need <= 32 bits, get them and shift away what isn't needed
        if bits <= 32:
            randVal = self.cRand.randint32()
            return randVal >> (32 - bits)

        numBytes = (bits + 7) // 8
        # setup a buffer to hold the output
        bitBuffer = <uint8_t*>PyMem_Malloc(numBytes)
        # fill the buffer with random bits
        for i in range(1, numBytes-1):
            bitBuffer[i] = self.cRand.randint8()

        # trim what isn't needed off the top
        bitBuffer[0] = self.cRand.randint8() >> ((numBytes * 8) - bits)
        return int.from_bytes(bitBuffer, 'little')

class pcg(cRandom):
    '''
    Create an instance of the PCG PRNG

    Underlying C implmentation by Melissa E. O'Neill:
    http://www.pcg-random.org

    Args:
        SeedVal (int):
            Whatever 32 bit value you would like to seed with. If no
            value is supplied one will be generated for you.
    '''
    seedSize = 32
    def __init__(self, seed = None):
        self.seedVal = HQrandom.getrandbits(32)
        self.cRand = cpcg(self.seedVal)
        if seed is not None:
            self.seed(seed)
        super().__init__(self.seedVal)

    def random(self):
        '''
        Get a random float, float has 32bits or randomness

        Returns:
            (float):A number between 0 and 1
        '''
        return self.cRand.randint32() / ((2 ** 32) - 1)

class sfmt(cRandom):
    '''
    Create an instance of the SFMT PRNG.

    Underlying C implmentation by Makoto Matsumoto and Mutsuo Saito:
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/

    Args:
        SeedVal (int or bytes):
            Whatever 32 bit value you would like to seed with. If no
            value is supplied one will be generated for you.
    '''
    seedSize = 32
    def __init__(self, seed = None):
        self.seedVal = HQrandom.getrandbits(32)
        self.cRand = csfmt(self.seedVal)
        if seed is not None:
            self.seed(seed)
        super().__init__(self.seedVal)

class xoshiro(cRandom):
    '''
    Create an instance of the xoshiro PRNG

    Underlying C implmentation by David Blackman and Sebastiano Vigna:
    http://xoshiro.di.unimi.it

    Args:
        SeedVal:
            Whatever value you would like to seed with. If no
            value is supplied one will be generated for you.
    '''
    seedSize = 256
    def __init__(self, seed = None):
        if seed is not None:
            self.seedVal = seed
        else:
            self.seedVal = HQrandom.getrandbits(256)
        self.cRand = cxoshiro(self.seedVal)
        self.seed(self.seedVal)
        super().__init__(self.seedVal)

    def seed(self, seedVal, version=2):
        '''
        Seed the PRNG using a 256 bit value.  Any value not long
        enough will be deterministically padded little endian to 
        avoid zero state in the PRNG.  Any value too long will be
        truncated to its uppper 256 bits.

        Args:
            seedVal:
                the int or bytes object containing the 
        '''
        cdef int i
        if type(seedVal) != int and type(seedVal) != bytes:
            raise TypeError('The seed must either be a bytes like object or an int.')
        if type(seedVal) == int:
            # use bitwise and to truncate the int to 256 bits
            seedVal = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF & seedVal
            seedVal = seedVal.to_bytes(32, 'little')
        # if we are less than 32 pad with determinstic non-zero padding.
        elif len(seedVal) < 32 and type(seedVal) == bytes:
            padLen = 32 - len(seedVal)
            for i in range(padLen):
                seedVal += bytes([i])
        self.cRand.seed(seedVal)

class ChaCha20(Random):
    '''
    Create an instance of the ChaCha20 PRNG

    Args:
        SeedVal (bytes or int):
            Whatever value you would like to seed with. If
            no value is supplied one will be generated for you.
    '''
    seedSize = 256
    def __init__(self, seedVal = None):
        if seedVal is None:
            seedVal = HQrandom.getrandbits(256)
        self.seed(seedVal)
        super().__init__(seedVal)
        
    def seed(self, seedVal, version = 2):
        '''
        Seed the PRNG using a 256 bit vaue.  Any value not long
        enough will be deterministically padded little endian to 
        avoid zero state in the PRNG.  Any value too long will be
        truncated to its uppper 256 bits.

        Args:
            seedVal (int or bytes):
                The initialization vector for the PRNG
            version (int):
                This is unused and only exists to maintain
                compatabliy with the Python random API.
        '''
        cdef int i
        # if seedVal is None:
        #     return
        if type(seedVal) != int and type(seedVal) != bytes:
            raise TypeError('The seed must either be a bytes like object or an int not %s.' % str(type(seedVal)))
        if type(seedVal) == int:
            # use bitwise and to truncate the int to 256 bits
            seedVal = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF & seedVal
            seedVal = seedVal.to_bytes(32, 'little')
        # if we are less than 32 pad with determinstic non-zero padding.
        elif len(seedVal) < 32 and type(seedVal) == bytes:
            padLen = 32 - len(seedVal)
            for i in range(padLen):
                seedVal += bytes([i])
        self.seedVal = int.from_bytes(seedVal, 'little')
        # setup PRNG
        self.cipher = chacha20.new(key=seedVal, nonce=seedVal[:8])
        self.feedback = self.cipher.encrypt(seedVal)

    def randint256(self):
        '''
        Get a random 256bit int

        Returns:
            (int):a 256 bit random int.
        '''
        self.feedback = self.cipher.encrypt(self.feedback)
        return int.from_bytes(self.feedback, 'little')
        
    def random(self):
        '''
        Get a random float, the float draws from 256 bits of
        randomness.

        Returns:
            (float):A number between 0 and 1
        '''
        return self.randint256() / ((2**256) -1)

    def getrandbits(self, bits):
        '''
        Returns an int composed of a user specified number of random
        bits. This is implemented so that all features of the Random
        parent class will funciton using this PRNG

        Args:
            bits (int):
                The number of bits you would like from the PRNG
        Returns:
            (int):An int composed of the random bits.
        '''
        if bits <= 0:
            raise ValueError('number of bits must be greater than 0')
        # if we need <= 256 bits we can just get them and bit shift
        # away the ones we don't want.
        if bits <= 256:
            retVal = self.randint256()
            return retVal >> (256 - bits)
        bitBuffer = b''
        # construct a bytestring 256bits at a time untill we have at
        # least enough bits
        while len(bitBuffer)*8 <= bits:
            self.feedback = self.cipher.encrypt(self.feedback)
            bitBuffer += self.feedback
        # Convert the bytes object to an int and bitshift away the
        # extra bits.
        return int.from_bytes(bitBuffer) >> (len(bitBuffer)*8 - bits)

# Here be introspective dragons... seriously its kind of a mess and
# I'm sorry.  We use some of Python's inrospection features to build
# refferences to classes from this module's internal dict.
def getRng(rngName):
    '''
    Get a refference to an PRNG class from string

    Args:
        rngName (str):
            A string containing the name of the PRNG class you want.
    Returns:
        A reference to the requested class
    
    Raises:
        NotImplementedError:
            Raised if the provided string doesn't correspond to anything in this module
    '''
    # get a ref to the current module (PRNG)
    currentModule = sys.modules[__name__]
    # check to see if the current module has anything matchng the
    # supplied name.
    if hasattr(currentModule, rngName):
        # return a ref to the object if so
        return getattr(currentModule, rngName)
    else:
        # raise NotImplmentedError if we don't have that PRNG
        raise NotImplementedError('RNG source %d not implmented.' % rngName)


def getListOfRNGSources():
    '''
    Get a tuple of the names of the RNG sources implmented by this
    module which python is alowed to call.

    Returns:
        (tuple):A tuple of the names of the RNG sources implmented by
        this module.
    '''
    # get a ref to the current module (PRNG)
    currentModule = sys.modules[__name__]
    RNGClasses = []
    # create a generator to get memebers from the module by name.
    moduleMembers = ((getattr(currentModule, name), name) for name in dir(currentModule))
    for member, name in moduleMembers:
        # remove things that aren't classes
        if isinstance(member, type):
            # remove things that aren't subclassed from Random, like
            # cdef RNG sources (csfmt, cxoshiro, cpcg)
            if issubclass(member, Random):
                RNGClasses.append(name)
    return tuple(RNGClasses)