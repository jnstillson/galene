'''
This file will be for tracking utility stuff, wall times, error handling etc


'''

def get_lig_set_properties(lig_set):

    mmass = []
    logp = []
    mr = []
    h_don = []
    h_acc = []

    for lig in lig_set['ligand']:
        mmass.append(lig.mmass)
        logp.append(lig.logp)
        mr.append(lig.mr)
        h_don.append(lig.h_don)
        h_acc.append(lig.h_acc)

    lig_set_prop = {'size' : len(lig_set['ligand']),
                    'av mmass' : sum(mmass)/len(mmass),
                    'av logp' : sum(logp)/len(logp),
                    'av mr' : sum(mr)/len(mr),
                    'av h_don' : sum(h_don)/len(h_don),
                    'av h_acc' : sum(h_acc)/len(h_acc)
                    }
    return lig_set_prop

def get_rec_set_properties(rec_set):

    rec_set_prop = {'size' : len(rec_set),}

    return rec_set_prop


class Error(Exception):
    pass


# need screen_library debug


# analyze_lig_set debug
class LigSetAnalysisError(Error):

    def __init__(self, message='\nError while precomputing values for unsupervised ligand analysis\n'):
        self.message = message
