#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#

# Exception classes for use in this code
class NoDistributionError(Exception):
    '''
    Exception raised when a distribution is expected but not present in an object.
    '''
    def __init__(self,value):
        self.value = value

class EncodingError(Exception):
    '''
    Exception raised when problems arise encoding a covert message.
    '''
    def __init__(self,value):
        self.value = value

class DecodingError(Exception):
    '''
    Exception raised when problems arise decoding a covert message.
    '''
    def __init__(self,value):
        self.value = value
