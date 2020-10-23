#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import rdkit

from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps as maps
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen #J. Chem. Inf. Comput. Sci. 1999, 39, 868-873



class Ligand():
    
    def __init__(self, name = None, smiles = None, mol = None):
        
        #basic ligand data
        self.name = name
        self.smiles = smiles
        self.mol = mol
        
        if smiles is not None:
            self.mol = AllChem.MolFromSmiles(smiles)
            
        if mol is not None:
            self.smiles = AllChem.MolToSmiles(mol)
        
        self.frag_set = None
        self.conformer_set = []
        self.vina_result_matrix = pd.DataFrame() 
        
        
        #degub data
        
        self.flag = False
        
        if self.mol is not None:
        
            self.cheminformatic_analysis()
            
        
    ## Debug ##
    
    def flag_issue(self):
        
        self.flag = True
    
    def cheminformatic_analysis(self):
        
        self.mmass = Descriptors.MolWt(self.mol)
        self.h_don = Lipinski.NumHDonors(self.mol)
        self.h_acc = Lipinski.NumHAcceptors(self.mol)
        self.logp = Crippen.MolLogP(self.mol)
        self.mr = Crippen.MolMR(self.mol)
        self.ecfp = AllChem.GetMorganFingerprintAsBitVect(self.mol,2) #obviously these needs to be improved

    
    ## Ligand Preparation ##
        
        
    #prepares the 3D structure of the molecule to yeild a useful pdb block
    
    def prep_struct(self):
    
        # Assign Structure
        self.mol = AllChem.MolFromSmiles(self.smiles)
        
        # Add Hydrogens
        AllChem.AddHs(self.mol) #account for pKa!!
    
        #Generate 3D Coordinates
        AllChem.EmbedMolecule(self.mol, randomSeed = 0xf00d) #random seed???
        
        #Minimize coordinates
        AllChem.MMFFOptimizeMolecule(self.mol)

        ### NOTE if you want any precomputed 3D dependent coordinates this would be the place to compute them
        
           

    
    ## Miscellaneous ##
    
    
    #compute partial charges and map them onto an embedded image
        
    def map_charge_dist(self):
        
        #recreated instead of calling self.mol so that a 2D version can be constructed
        mol = AllChem.MolFromSmiles(self.smiles)
        
        #Determine Partial Charges 
        AllChem.ComputeGasteigerCharges(mol)
    
        #gasteiger charges
        full_charge_gasteiger = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
    
        #extended Huekel theory (eHT) charges?

        #flattening to 2D was problematic?
    
        maps.GetSimilarityMapFromWeights(mol, full_charge_gasteiger)
