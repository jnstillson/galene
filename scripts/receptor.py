#!/usr/bin/env python
# coding: utf-8
from rdkit.Chem import AllChem


class Receptor():

    def __init__(self,
                 name,
                 rig_file,
                 pos,
                 dim,
                 pdb_file,
                 ):
        self.name = name

        self.rig_file = rig_file
        self.pdb_file = pdb_file
        self.mol_block = AllChem.MolToPDBBlock(AllChem.MolFromPDBFile(pdb_file, sanitize=False)).encode(encoding='UTF-8')

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
