'''
This file will be for tracking utility stuff, wall times, error handling etc


'''

class Error(Exception):
    pass


#analyze_lig_set debug
class SimMatrixError(Error):
    
    def __init__(self, message = '\nError computing the similiariy matrix\n'):
        
        self.message = message

class PCAError(Error):
    
    def __init__(self, message = '\nError performing PCA on similarity matrix'):
        
        self.message = message