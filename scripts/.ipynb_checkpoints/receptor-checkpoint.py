#!/usr/bin/env python
# coding: utf-8
from rdkit.Chem import AllChem


class Receptor():

    def __init__(self,
                 name,
                 rig_file,
                 pos,
                 dim,
                 ):
        self.name = name

        self.rig_file = rig_file

        # self.mol = AllChem.MolFromMol2File(mol2_file, sanitize = False)

        self.pos = pos
        self.dim = dim

        # in the event of receptor analysis:
        '''
        self.sites = pd.DataFrame({'name' : None,
                                   'pdb' : None,
                                   'pos' : None,
                                   'dim' : None,
                                   })
                                   
        '''
