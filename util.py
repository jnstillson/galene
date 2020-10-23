'''
This file will be for tracking utility stuff, wall times, error handling etc


'''


class Error(Exception):
    pass


# need screen_library debug


# analyze_lig_set debug
class LigSetAnalysisError(Error):

    def __init__(self, message='\nError while precomputing values for unsupervised ligand analysis\n'):
        self.message = message
