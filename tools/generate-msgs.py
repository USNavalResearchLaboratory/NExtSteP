#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from os.path import isfile
import json
import os

def genRandom(nbits=4000000):
    '''
    Generate random data for embedding using os.urandom
    
    Save the data as a hex string.  Confirm that the correct number of bits was produced and that the file can be decoded to match what was generated.
    '''
    
    # put the file one level up from current directory
    ofile = '../urand-b' + str(nbits) + '.txt'

    if isfile(ofile):
        print('FILE: Output file exists: ' + ofile)
        print('Exiting without overwriting the file.')
    else:
        bytes_ = nbits//8
        if bytes_ * 8 < nbits:
            bytes_ += 1
    
        urandStr = os.urandom(bytes_)
        urandHexStr = urandStr.hex()
        urandBinStr = "".join(["{0:04b}".format(int(str(c),16)) for c in urandHexStr])
    
        print('Generated', len(urandBinStr), 'bits using os.urandom.')

        with open(ofile, 'w') as ostream:
            ostream.write(urandHexStr)
    
        with open(ofile, 'r') as istream:
            line = istream.readline()

        print('Line written to file has length:', len(line))
        confirmBinStr = "".join(["{0:04b}".format(int(c,16)) for c in line])
        print('Line written to file produces', len(confirmBinStr), 'bits.')
        print('Line written to file has', confirmBinStr.count('1'), '1s.')
        
        OK = True
        for i in range(len(urandBinStr)):
            if urandBinStr[i] != confirmBinStr[i]:
                OK = False
        if OK:
            print('Binary strings are the same.')
        else:
            print('Binary strings differ!')


##################

if __name__ == "__main__":
    
    genRandom(4000000)

