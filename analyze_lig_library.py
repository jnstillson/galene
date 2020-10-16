from .screen_library import * 
import pandas as pd
import numpy as np
from sklearn import decomposition
from sklearn import cluster
from sklearn import metrics
from rdkit.Chem import Recap
from itertools import compress

class analyze_lig_library():
    
    def __init__(self, lig_set = None):
        
        self.lig_set = lig_set
        self.lig_ids = [self.lig_set.index(lig) for lig in self.lig_set]
        self.similarity_matrix = pd.DataFrame(columns = self.lig_ids, index = self.lig_ids)
        
        self.lig_space2 = None
        self.frag_set = pd.DataFrame(columns = ['LigandID','FragID','SMILES','Frag'])
        
        
        self.lig_grouped = None
    
    
    ### Lipinski Type Values ###    
        
        
        
    def logp_spread(self):
        
        spread = pd.Series([lig.logp for lig in self.lig_set])
        spread.hist()
        
    def mmass_spread(self):
        
        spread = pd.Series([lig.mmass for lig in self.lig_set])
        spread.hist()
        
        
        
    ### PCA Analysis of similarity matrix and k-means clustering ###
    
    
        
    def compute_similarity_index(self, fp1 = None, fp2 = None):
        
        return rdkit.DataStructs.FingerprintSimilarity(fp1, fp2) 
    
    
    def compute_similarity_matrix(self):
        
        for lig1 in self.lig_set:
            
            for lig2 in self.lig_set:
                
                try:
                    
                    i = self.compute_similarity_index(lig1.ecfp, lig2.ecfp)
                    
                except: 
                    
                    i = 0

                self.similarity_matrix.loc[self.lig_set.index(lig1), self.lig_set.index(lig2)] =  i
                
                
    def pca_similarity_matrix(self):
        
        pca = decomposition.PCA(n_components = 2)
            
        pca.fit(self.similarity_matrix)
        
        self.lig_space2 = pd.DataFrame(pca.transform(self.similarity_matrix))
        
        self.lig_space2.plot.scatter(x = 0,y = 1)
        
    def k_means_plot(self, clusters = 5, best_cluster = False):
    
        kmeans = cluster.KMeans(init = 'k-means++', n_clusters = clusters, n_init = 10, max_iter = 300)
        
        kmeans.fit(self.lig_space2.values)
        
        self.lig_space2.plot.scatter(x = 0, y = 1, c = kmeans.labels_, colormap = 'viridis')
        
        if best_cluster:
            
            self.lig_grouped = pd.Series(self.lig_set).groupby(kmeans.labels_)
            
        
    def k_means_analysis(self, k_min = 2, k_max = 8):
        
        sil_scores = []
        
        for k in range(k_min, k_max):
            
            kmeans = cluster.KMeans(init = 'k-means++', n_clusters = k, n_init = 10, max_iter = 300)
            
            kmeans.fit(self.lig_space2.values)
            
            sil_score = metrics.silhouette_score(self.lig_space2.values, kmeans.labels_)
            
            sil_scores.append(sil_score)
            
        best_score = max(sil_scores)
            
        best_score_index = sil_scores.index(best_score) + k_min
        
        print('\nBest Fit (by silhouette score) found at: '+ str(best_score_index) + ': ' + str(best_score) + '\n')
        
        self.k_means_plot(best_score_index, True)
        
        
    #need functions that encompass a series of analytical functions (compute, pca, k-means) and
    #produce a single graphic at the end.
    
    
    
    ### Extraction of Cores / Scaffold Analysis ###
    
    
    ## Adapted from the code provided by Medina-Franco and Naveja (2019) ##
    
    
    def single_find_cores(self, mol, c = 2/3):
        
        smi = AllChem.MolToSmiles(mol)
        mol_h = mol.GetNumHeavyAtoms()
        
        cores_full = list(Recap.RecapDecompose(mol).GetAllChildren().keys()) + [smi]
        cores_comp = list(compress(cores_full,[AllChem.MolFromSmiles(core).GetNumHeavyAtoms() >= c*mol_h for core in cores_full]))
        
        frags = [AllChem.MolFromSmiles(core) for core in cores_comp]
        
        return frags
    
    def find_cores(self, c = 2/3):
        
        ligID = 0
        
        for lig in self.lig_set:
            
            fragID = 0
            
            for frag in self.single_find_cores(lig.mol, c):
            
                self.frag_set = self.frag_set.append({'LigandID' : ligID,
                                      'FragID' : fragID,
                                      'SMILES' : AllChem.MolToSmiles(frag),
                                      'Frag' : frag}, ignore_index = True)
                
                fragID += 1
                
            ligID += 1
            
    def sort_cores(self):
        
        sorted_cores = self.frag_set.groupby("SMILES").LigandID.apply(lambda x: str(sorted(set(x))))
        
        df = pd.DataFrame({'Fragment' : sorted_cores.index, 'LigandID' : sorted_cores.values})
        
        return df
        
        
            
        
        
        
        
        
        
            
            
        
        
        
        
        
        
        
