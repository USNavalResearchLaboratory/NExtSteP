#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from libStealthy.exceptions import EncodingError
from abc import ABCMeta, abstractmethod
from bisect import bisect_left

from libStealthy.util import hqrand

import numpy as np

class BaseStegoFromIPDs(object, metaclass=ABCMeta):
    '''
    Construct a fresh IPD sequence that encodes msg, given as an
    iterable from {'0','1'}*.  This uses an IPD sequence that is
    given when the BaseStegoFromIPDs (or sub-class) object is
    constructed.

    A typical approach is to first construct a distribution from the
    input IPD sequence and then use this distributiion to construct
    the new IPD sequence.

    Subclasses could also use the input IPD sequence (given at
    construction time)
    '''

    def __init__(self, IPDseq, RNGSource, seed = None):
        self.initIPDs = [x for x in IPDseq]

    @abstractmethod
    def generateIPDseq(self, msg, numt):
        pass


class SCTimeReplayFromIPDs(BaseStegoFromIPDs):
    '''
    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use it
    to construct an IPD sequence containing a stego message inserted
    using time replays.

    Pad the length out to the requested length numt if specified
    (using the constructed distribution to independently generate
    each additional IPD).
    '''

    def __init__(self, IPDseq, RNGSource, seed = None):
        '''
        Construct the distribution upon creating the object
        '''
        self.HQrand = RNGSource(seed)
        self.initIPDs = [x for x in IPDseq]
        self.dist = lambda : self.HQrand.choice(self.initIPDs)
        self.medInCover = np.median(self.initIPDs)
        self.vecForZero = [x for x in self.initIPDs if x < self.medInCover]
        self.vecForOne = [x for x in self.initIPDs if x > self.medInCover]

        if len(self.vecForZero) == 0  or len(self.vecForOne) == 0:
            raise ValueError('Need more distinct IPD cover values!')

        self.distForZero = lambda : self.HQrand.choice(self.vecForZero)
        self.distForOne = lambda : self.HQrand.choice(self.vecForOne)

    def generateIPDseq(self, msg, numt=-1):
        '''
        Construct the stego traffic
        '''
        if numt < 0:
            numt = len(msg)

        if numt < len(msg):
            raise EncodingError('Message too long for requested number of IPDs!')

        stegoIPDs = []
        for i in range(len(msg)):
            if msg[i] == '0':
                stegoIPDs.append(self.distForZero())
            else:
                stegoIPDs.append(self.distForOne())

        # Pad out to numt IPDs
        if len(msg) < numt:
            for i in range(numt - len(msg)):
                stegoIPDs.append(self.dist())

        return stegoIPDs


class SCTimeReplayFromIPDsUE(BaseStegoFromIPDs):
    '''
    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use
    it to construct an IPD sequence containing a stego message
    inserted using time replays.

    If the message is shorter than the requested number of IPDs,
    embed it uniformly (but with constant offset) into a cover
    sequence.  If the requested number of IPDs is at most the number
    given at instantiation, use an initial prefix of that IPD
    sequence.  If it is more than the number of IPDs given at
    initialization, then generate the cover sequence using the
    distribution defined by the initialization IPDs.
    '''

    def __init__(self, IPDseq, RNGSource, seed = None):
        '''
        Construct the distribution upon creating the object
        '''
        self.HQrand = RNGSource(seed)
        # self.HQrand.random = hqrand
        self.initIPDs = [x for x in IPDseq]
        self.dist = lambda : self.HQrand.choice(self.initIPDs)
        self.medInCover = np.median(self.initIPDs)
        self.vecForZero = [x for x in self.initIPDs if x < self.medInCover]
        self.vecForOne = [x for x in self.initIPDs if x > self.medInCover]

        if len(self.vecForZero) == 0  or len(self.vecForOne) == 0:
            raise ValueError('Need more distinct IPD cover values!')

        self.distForZero = lambda : self.HQrand.choice(self.vecForZero)
        self.distForOne = lambda : self.HQrand.choice(self.vecForOne)

    def generateIPDseq(self, msg, numt=-1):
        '''
        Construct the stego traffic
        '''
        if numt < 0:
            numt = len(msg)

        if numt < len(msg):
            raise EncodingError('Message too long for requested number of IPDs!')

        stegoIPDs = []

        if len(msg) < numt:
            if numt <= len(self.initIPDs):
                stegoIPDs = self.initIPDs[0:numt]
            else:
                stegoIPDs = [self.HQrand.choice(self.initIPDs) for i in range(numt)]

            di = int(numt/len(msg))
            if di < 1:
                raise ValueError('Trying to embed with di = 0!')

            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs[i * di] = self.distForZero()
                else:
                    stegoIPDs[i * di] = self.distForOne()

                ### Modify stegoIPDs
        else:
            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs.append(self.distForZero())
                else:
                    stegoIPDs.append(self.distForOne())

            # Pad out to numt IPDs
            if len(msg) < numt:
                for i in range(numt - len(msg)):
                    stegoIPDs.append(self.dist())


        return stegoIPDs


class SCTimeReplayPercentileFromIPDsUE(BaseStegoFromIPDs):
    '''
    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use it
    to construct an IPD sequence containing a stego message inserted
    using skewed time replays.  The parameter q (in [0,100]) used to
    construct the object gives the percentile in IPD values used to
    represent a 0.  The IPD sequence used to construct the object is
    sorted, and the value at percentile q (0 = min., 100 = max.) is
    identified.  To encode a 0, an IPD from those below the
    identified value is selected at random from those observed (with
    duplicates increasing the probability); to encode a 1, an IPD
    from those above the identified is selected at random.

    If the message is shorter than the requested number of IPDs,
    embed it uniformly (but with constant offset) into a cover
    sequence.  If the requested number of IPDs is at most the number
    given at instantiation, use an initial prefix of that IPD
    sequence.  If it is more than the number of IPDs given at
    initialization, then generate the cover sequence using the
    distribution defined by the initialization IPDs.
    '''

    def __init__(self, IPDseq, RNGSource, q, seed = None):
        '''
        Construct the distribution upon creating the object
        '''
        self.HQrand = RNGSource(seed)
        # self.HQrand.random = hqrand  
        self.initIPDs = [x for x in IPDseq]
        self.dist = lambda : self.HQrand.choice(self.initIPDs)
        self.percInCover = np.percentile(self.initIPDs, q)
        self.vecForZero = [x for x in self.initIPDs if x < self.percInCover]
        self.vecForOne = [x for x in self.initIPDs if x > self.percInCover]

        if len(self.vecForZero) == 0  or len(self.vecForOne) == 0:
            raise ValueError('Need more distinct IPD cover values!')

        self.distForZero = lambda : self.HQrand.choice(self.vecForZero)
        self.distForOne = lambda : self.HQrand.choice(self.vecForOne)

    def generateIPDseq(self, msg, numt=-1):
        '''
        Construct the stego traffic
        '''
        if numt < 0:
            numt = len(msg)

        if numt < len(msg):
            raise EncodingError('Message too long for requested number of IPDs!')

        stegoIPDs = []

        if len(msg) < numt:
            if numt <= len(self.initIPDs):
                stegoIPDs = self.initIPDs[0:numt]
            else:
                stegoIPDs = [self.HQrand.choice(self.initIPDs) for i in range(numt)]

            di = int(numt/len(msg))
            if di < 1:
                raise ValueError('Trying to embed with di = 0!')

            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs[i * di] = self.distForZero()
                else:
                    stegoIPDs[i * di] = self.distForOne()

                ### Modify stegoIPDs
        else:
            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs.append(self.distForZero())
                else:
                    stegoIPDs.append(self.distForOne())

            # Pad out to numt IPDs
            if len(msg) < numt:
                for i in range(numt - len(msg)):
                    stegoIPDs.append(self.dist())


        return stegoIPDs


class SCTimeReplayPercentileFromIPDsUEHQ(BaseStegoFromIPDs):
    '''
    High-quality version; also adds <= in computing vecForZero

    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use it
    to construct an IPD sequence containing a stego message inserted
    using skewed time replays.  The parameter q (in [0,100]) used to
    construct the object gives the percentile in IPD values used to
    represent a 0.  The IPD sequence used to construct the object is
    sorted, and the value at percentile q (0 = min., 100 = max.) is
    identified.  To encode a 0, an IPD from those below the
    identified value is selected at random from those observed (with
    duplicates increasing the probability); to encode a 1, an IPD
    from those above the identified is selected at random.

    If the message is shorter than the requested number of IPDs,
    embed it uniformly (but with constant offset) into a cover sequence.
    If the requested number of IPDs is at most the number given at
    instantiation, use an initial prefix of that IPD sequence.  If
    it is more than the number of IPDs given at initialization, then
    generate the cover sequence using the distribution defined by the
    initialization IPDs.
    '''

    def __init__(self, IPDseq, RNGSource, q, seed = None):
        '''
        Construct the distribution upon creating the object
        '''

        self.randHQ = RNGSource(seed)
        # self.randHQ.random = hqrand  
        self.initIPDs = [x for x in IPDseq]
        self.dist = lambda : self.randHQ.choice(self.initIPDs)
        self.percInCover = np.percentile(self.initIPDs, q)
        self.vecForZero = [x for x in self.initIPDs if x <= self.percInCover]
        self.vecForOne = [x for x in self.initIPDs if x > self.percInCover]

        if len(self.vecForZero) == 0  or len(self.vecForOne) == 0:
            raise ValueError('Need more distinct IPD cover values!')

        self.distForZero = lambda : self.randHQ.choice(self.vecForZero)
        self.distForOne = lambda : self.randHQ.choice(self.vecForOne)

    def generateIPDseq(self, msg, numt=-1):
        '''
        Construct the stego traffic
        '''
        if numt < 0:
            numt = len(msg)

        if numt < len(msg):
            excStr = 'Message too long for requested number of IPDs! numt = '+str(numt)+', len(msg) = '+str(len(msg))
            print(excStr)
            raise EncodingError('Message too long for requested number of IPDs!')

        stegoIPDs = []

        if len(msg) < numt:
            if numt <= len(self.initIPDs):
                stegoIPDs = self.initIPDs[0:numt]
            else:
                stegoIPDs = [self.randHQ.choice(self.initIPDs) for i in range(numt)]

            di = int(numt/len(msg))
            if di < 1:
                raise ValueError('Trying to embed with di = 0!')

            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs[i * di] = self.distForZero()
                else:
                    stegoIPDs[i * di] = self.distForOne()

                ### Modify stegoIPDs
        else:
            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs.append(self.distForZero())
                else:
                    stegoIPDs.append(self.distForOne())

            # Pad out to numt IPDs
            if len(msg) < numt:
                for i in range(numt - len(msg)):
                    stegoIPDs.append(self.dist())


        return stegoIPDs


class SCTimeReplayPercentileFromOtherIPDsUEHQ(BaseStegoFromIPDs):
    '''
    High-quality version; also adds <= in computing vecForZero

    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use it
    to embed a stego message into a (potentially different) IPD
    sequence using skewed time replays.  The parameter q (in [0,100])
    used to construct the object gives the percentile in IPD values
    used to represent a 0.  The IPD sequence used to construct the
    object is sorted, and the value at percentile q
    (0 = min., 100 = max.) is identified.  To encode a 0, an IPD from
    those below the identified value is selected at random from those
    observed (with duplicates increasing the probability); to encode
    a 1, an IPD from those above the identified is selected at random.

    If the message is shorter than the sequence into which it is being
    embedded, embed it uniformly (but with constant offset) into the
    cover sequence.
    '''

    def __init__(self, IPDseq, RNGSource, q, seed = None):
        '''
        Construct the distribution upon creating the object

        The IPD sequence given at construction is used to generate
        the IPDs that encode the message.  These are embedded into a
        (potentially different) IPD sequence.
        '''

        self.randHQ = RNGSource(seed)
        self.initIPDs = [x for x in IPDseq]
        self.dist = lambda : self.randHQ.choice(self.initIPDs)
        self.percInCover = np.percentile(self.initIPDs, q)
        self.vecForZero = [x for x in self.initIPDs if x <= self.percInCover]
        self.vecForOne = [x for x in self.initIPDs if x > self.percInCover]

        if len(self.vecForZero) == 0  or len(self.vecForOne) == 0:
            raise ValueError('Need more distinct IPD cover values!')

        self.distForZero = lambda : self.randHQ.choice(self.vecForZero)
        self.distForOne = lambda : self.randHQ.choice(self.vecForOne)

    def generateIPDseq(self, msg, hostSeq):
        '''
        Embed the stego traffic into the host sequence; return the result
        '''


        if len(hostSeq) < len(msg):
            excStr = 'Message too long for requested number of IPDs! len(hostSeq) = '+str(len(hostSeq))+', len(msg) = '+str(len(msg))
            print(excStr)
            raise EncodingError('Message too long for requested number of IPDs!')

        elif len(msg) == 0:
            stegoIPDs = hostSeq[:]

        elif len(msg) < len(hostSeq):

            stegoIPDs = hostSeq[:]

            di = int(len(hostSeq)/len(msg))
            if di < 1:
                raise ValueError('Trying to embed with di = 0!')

            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs[i * di] = self.distForZero()
                else:
                    stegoIPDs[i * di] = self.distForOne()


        else: # Message is exactly as long as host sequence
            stegoIPDs = []
            for i in range(len(msg)):
                if msg[i] == '0':
                    stegoIPDs.append(self.distForZero())
                else:
                    stegoIPDs.append(self.distForOne())

        return stegoIPDs


class SCTimeReplayProbFromIPDsUEHQ(BaseStegoFromIPDs):
    '''
    High-quality version; also adds <= in computing vecForZero

    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use it
    to construct an IPD sequence containing a stego message inserted
    using probabalistic skewed time replays.  The parameter q
    (in [0,100]) used to construct the object gives the percentile in
    IPD values used to represent a 0.  The IPD sequence used to
    construct the object is sorted, and the value at percentile q
    (0 = min., 100 = max.) is identified.  To encode a 0, an IPD from
    those below the identified value is selected at random from those
    observed (with duplicates increasing the probability); to encode
    a 1, an IPD from those above the identified is selected at random

    If the message is shorter than the requested number of IPDs,
    embed it uniformly (but with constant offset) into a cover
    sequence.  If the requested number of IPDs is at most the number
    given at instantiation, use an initial prefix of that IPD
    sequence as cover sequence.  If it is more than the number of
    IPDs given at initialization, then generate the cover sequence
    using the distribution defined by the initialization IPDs.
    '''

    def __init__(self, IPDseq, RNGSource, q, seed = None):
        '''
            Construct the distribution upon creating the object

            '''
        self.randHQ = RNGSource(seed)
        self.initIPDs = [x for x in IPDseq]
        self.dist = lambda: self.randHQ.choice(self.initIPDs)
        self.percInCover = np.percentile(self.initIPDs, q)
        self.vecForZero = [x for x in self.initIPDs if x <= self.percInCover]
        self.vecForOne = [x for x in self.initIPDs if x > self.percInCover]

        if len(self.vecForZero) == 0  or len(self.vecForOne) == 0:
            raise ValueError('Need more distinct IPD cover values!')

        self.distForZero = lambda : self.randHQ.choice(self.vecForZero)
        self.distForOne = lambda : self.randHQ.choice(self.vecForOne)


    def generateIPDseq(self, msg, prob,numt=-1):
        '''
        Construct the stego traffic
        '''

        if prob <0 or prob>1.0:
            excStr='Probability must be between 0 and 1'
            print(excStr)
            raise EncodingError('Invalid probability used')

        if numt < 0:
            numt = len(msg)

        if numt < len(msg):
            excStr = 'Message too long for requested number of IPDs! numt = '+str(numt)+', len(msg) = '+str(len(msg))
            print(excStr)
            raise EncodingError('Message too long for requested number of IPDs!')

        stegoIPDs = []

        if numt <= len(self.initIPDs):
            stegoIPDs = self.initIPDs[0:numt]
            # use IPD from construction of object to initialize stego
        else:
            stegoIPDs = [self.randHQ.choice(self.initIPDs) for i in range(numt)] 
            #sample from distribution to intialize stego.

        msgbitnum=0
        # embedd_word is the word exchanged by Alice and Bob before
        # embedding to decide which cover bits message will be
        # embedded into
        embedd_word=[]
        for i in range(numt):
            randnum=self.randHQ.random()
            if(randnum <prob):
                embedd_word.append('1')
                if (msg[msgbitnum]=='0' and (stegoIPDs[i]>self.percInCover)):
                    stegoIPDs[i]=self.distForZero()
                if (msg[msgbitnum] == '1' and (stegoIPDs[i] <= self.percInCover)):
                    stegoIPDs[i] = self.distForOne()

                msgbitnum+=1
            else:
                embedd_word.append('0')
        return stegoIPDs, msgbitnum


class SCTimeReplayProbFromIPDsUEHQ_alphabet(BaseStegoFromIPDs):
    '''
    High-quality version; also adds <= in computing vecForZero

    Given an IPD sequence and a message (given as an iterable from
    {'0','1'}*), generate a distribution from the sequence and use it
    to construct an IPD sequence containing a stego message inserted
    using probabalistic skewed time replays.  The parameter q
    (in [0,100]) used to construct the object gives the percentile in
    IPD values used to represent a 0.  The IPD sequence used to
    construct the object is sorted, and the value at percentile q
    (0 = min., 100 = max.) is identified.  To encode a 0, an IPD from
    those below the identified value is selected at random from those
    observed (with duplicates increasing the probability); to encode
    a 1, an IPD from those above the identified is selected at random

    If the message is shorter than the requested number of IPDs,
    embed it uniformly (but with constant offset) into a cover
    sequence.  If the requested number of IPDs is at most the number
    given at instantiation, use an initial prefix of that IPD
    sequence as cover sequence.  If it is more than the number of
    IPDs given at initialization, then generate the cover sequence
    using the distribution defined by the initialization IPDs.
    '''

    def __init__(self, IPDseq, RNGSource, q, seed = None):
        '''
            Construct the distribution upon creating the object

        '''
        self.randHQ = RNGSource(seed)
    
        self.initIPDs = np.array([x for x in IPDseq])
        self.dist = lambda: self.randHQ.choice(self.initIPDs)

        q.sort()

        self.alphabet_size=len(q)+1
        self.chunk_size=int(np.log2(self.alphabet_size))

        #gives a ndarray of IPDs that are the bin dividers 
        self.IPDcutoffs=np.percentile(self.initIPDs,q, interpolation='midpoint')
        #print('IPDinit:', self.initIPDs)
        #print('IPDcutoffs:', self.IPDcutoffs)
        #print('\n')
        
        #gives bin membership to each IPD
        self.initIPD_bin_memberships=np.digitize(self.initIPDs, self.IPDcutoffs)

        #group intial IPDs by membership
        self.binned_initIPDs=[self.initIPDs[self.initIPD_bin_memberships==i] for i in range(self.alphabet_size)] 
        if(0 in [len(group) for group in self.binned_initIPDs]): 
            raise ValueError('Need more distinct IPD cover values!')

    def generateIPDseq(self, msg, prob,numt=-1):
        '''
        Construct the stego traffic
        '''

        if prob <0 or prob>1.0:
            excStr='Probability must be between 0 and 1'
            print(excStr)
            raise EncodingError('Invalid probability used')

        if numt < 0:
            numt = int(len(msg)/self.chunk_size)

        stegoIPDs = []

        if numt <= len(self.initIPDs):
            stegoIPDs = self.initIPDs[0:numt]
            # use IPD from construction of object to initialize stego
        else:
            stegoIPDs = [self.randHQ.choice(self.initIPDs) for i in range(numt)] 
            #sample from distribution to intialize stego.

        msgsegnum=0
        ##convert 0th segment of message (which is a binary string) to integer 
        msgletter=int(msg[msgsegnum*self.chunk_size : (msgsegnum+1)*self.chunk_size],2)
        # embedd_word is the word exchanged by Alice and Bob before
        # embedding to decide which cover bits message will be
        # embedded into
        embedd_word=[]

        for i in range(numt):
            randnum=self.randHQ.random()
            if(randnum <prob):
                embedd_word.append('1')
                if(bisect_left(self.IPDcutoffs, stegoIPDs[i]) != msgletter):
                    stegoIPDs[i] = self.randHQ.choice(self.binned_initIPDs[msgletter])
                msgsegnum+=1
                #set msgletter to next message letter to encode 
                msgletter=int(msg[msgsegnum*self.chunk_size : (msgsegnum+1)*self.chunk_size],2)
            else:
                embedd_word.append('0')
        #return length of message (in bits) actually embedded. In larger alphabet, length of message embedded is msgsegnum
        return np.ndarray.tolist(stegoIPDs), msgsegnum*self.alphabet_size