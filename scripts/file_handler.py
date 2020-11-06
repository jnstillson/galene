'''
Module: file_handler

Purpose: Import, Export, File Handling for data sets
'''
from .ligand import *
from .receptor import *
from .screen_library import *
import pandas as pd


def lig_set_from_smi(smi):
    # need similar function to create lig_set from a smiles code

    # intake smiles file
    smifile = open(smi, 'r')
    lines = smifile.readlines()

    ligs = []
    names = []

    i = 0
    for line in lines:
        tabs = line.split()
        if len(tabs) == 2:
            lig = Ligand(smiles=str(tabs[0]), name=str(tabs[1]))
            lig.prep_default_conf()
            ligs.append(lig)
            names.append(str(tabs[1]))
        elif len(tabs) == 1:
            lig = Ligand(smiles=str(tabs[0]), name=f'lig{i}')
            lig.prep_default_conf()
            ligs.append(lig)
            names.append(f'lig{i}')
        i += 1

    lig_set = pd.DataFrame({'name': names,
                            'ligand': ligs,
                            'ligid': [i for i in range(len(ligs))],
                            })
    smifile.close()

    return lig_set


def lig_set_from_gypsum(sdf):
    # intake sdf file
    molsdf = AllChem.SDMolSupplier(sdf)

    mols = [mol for mol in molsdf][1:]

    # find unique names
    names = list(dict.fromkeys([mol.GetProp('_Name') for mol in mols]))

    lig_set = pd.DataFrame({'name': names,
                            'ligand': [Ligand() for i in range(len(names))],
                            'ligid': [i for i in range(len(names))],
                            'conf count': None,
                            'frags': None})

    for i in range(len(lig_set)):

        for mol in mols:
            if mol.GetProp('_Name') == lig_set.loc[i, 'name']:
                lig_set.loc[i, 'ligand'].conformer_set.append(mol)

        # count the conformers
        lig_set.loc[i, 'conf count'] = len(lig_set.loc[i, 'ligand'].conformer_set)

        # set the first conformer to be the representative of the molecule in 2D
        # (this could be problematic if input has ambiguous sterocentres)
        lig_set.loc[i, 'ligand'].smiles = AllChem.MolToSmiles(lig_set.loc[i, 'ligand'].conformer_set[0])
        lig_set.loc[i, 'ligand'].mol = AllChem.MolFromSmiles(lig_set.loc[i, 'ligand'].smiles)
        lig_set.loc[i, 'ligand'].name = names[i]

        lig_set.loc[i, 'ligand'].cheminformatic_analysis()

    return lig_set


def rec_set_from_file(txt):

    recs = []
    i = 0

    txtfile = open(txt, 'r')
    lines = txtfile.readlines()

    for line in lines:
        tabs = line.split()
        rec = Receptor(name=str(tabs[0]),
                       rig_file=str(tabs[7]),
                       pos=[tabs[1], tabs[2], tabs[3]],
                       dim=[tabs[4], tabs[5], tabs[6]],
                       pdb_file=tabs[8],
                       )
        recs.append(rec)

    rec_set = pd.DataFrame({'name' : [rec.name for rec in recs],
                            'receptor' : recs,
                            'recid' : [i for i in range(len(recs))],

               })

    txtfile.close()
    
    return rec_set


