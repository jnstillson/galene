'''
Module: file_handler

Purpose: Import, Export, File Handling for data sets
'''
from .ligand import *
import pandas as pd

def lig_set_from_smi(smi):
    
    #need similar function to create lig_set from a smiles code
    
    #intake smiles file
    smifile = open(smi, 'r')
    lines = smifile.readlines()
    
    ligs = []
    names = []
    
    for line in lines:
        tabs = line.split()
        
        ligs.append(Ligand(smiles = str(tabs[0])))
        names.append(str(tabs[1]))     
           

        
    lig_set = pd.DataFrame({'name' : names,
                            'ligand' : ligs,
                            'ligid' : [i for i in range(len(ligs))],
                            'frags' : None})
    
    
    return lig_set


def lig_set_from_gypsum(sdf):
    
    #intake sdf file 
    molsdf = AllChem.SDMolSupplier(sdf)
    
    mols = [mol for mol in molsdf][1:]
    
    #find unique names
    names = list(dict.fromkeys([mol.GetProp('_Name') for mol in mols]))
    
    
    lig_set = pd.DataFrame({'name' : names,
                            'ligand' : [Ligand() for i in range(len(names))],
                            'ligid' : [i for i in range(len(names))],
                            'conf count' : None,
                            'frags' : None})
    
    
    for i in range(len(lig_set)):
            
        for mol in mols:
            if mol.GetProp('_Name') == lig_set.loc[i, 'name']:
                lig_set.loc[i, 'ligand'].conformer_set.append(mol)
                    
        cols = [('Con'+str(i)) for i in range(len(lig_set.loc[i, 'ligand'].conformer_set))]
        lig_set.loc[i, 'ligand'].vina_result_matrix = pd.DataFrame(columns = cols) 
        
        #count the conformers
        lig_set.loc[i, 'conf count'] = len(lig_set.loc[i, 'ligand'].conformer_set)
                
        #set the first conformer to be the representative of the molecule in 2D 
        #(this could be problematic if input has ambiguous sterocentres)
        lig_set.loc[i, 'ligand'].smiles = AllChem.MolToSmiles(lig_set.loc[i, 'ligand'].conformer_set[0])
        lig_set.loc[i, 'ligand'].mol = AllChem.MolFromSmiles(lig_set.loc[i, 'ligand'].smiles)   
    
    return lig_set