from .screen_library import *
import pandas as pd
import numpy as np
from sklearn import decomposition
from sklearn import cluster
from sklearn import metrics
from rdkit.Chem import Recap
from itertools import compress
import networkx
import matplotlib.pyplot as plt

class analyze_lig_set():
    
    def __init__(self, lig_set = None):
        
        self.lig_set = lig_set
        
        '''
        # New Formats #
        lig_set ['ligand','ligid','frags'] 
        frag_set ['fragment','fragid','ligs']
        where  ligand / fragment : ligand()
               ligid /fragid : a unique ligand identifier
               frags / ligs : a set of ids defining the two surjective maps
                                   m: lig_set -> frag_set 
                                   n: frag_set -> lig_set
        
        Note that the surjectivity of both maps is assured becuase fragmenting in this context has an identity
        component (i.e. the molecule will always be its own fragment)
                        
        '''
        self.similarity_matrix = pd.DataFrame(columns = self.lig_set['ligid'], index = self.lig_set['ligid'])
        self.lig_space2 = None
        
        
        self.frag_set = pd.DataFrame(columns = ['fragment','fragid','ligs'])
        
        self.fsimilarity_matrix = None
        self.frag_space2 = None
        
        
        
        #will likeley be replaced by neighbourhood and cluster searching
        self.lig_grouped = None
        
    
    
    ### Lipinski Type Values ###    
        
        
        
    def logp_spread(self):
        
        spread = pd.Series([lig.logp for lig in self.lig_set['ligand']])
        spread.hist()
        
    def mmass_spread(self):
        
        spread = pd.Series([lig.mmass for lig in self.lig_set['ligand']])
        spread.hist()
        
        
        
    ### PCA Analysis of similarity matrix and k-means clustering ###
    
    
        
    def compute_similarity_index(self, fp1 = None, fp2 = None):
        
        return rdkit.DataStructs.FingerprintSimilarity(fp1, fp2) 
    
    
    def compute_similarity_matrix(self):
        
        index_0 = 0
        
        for lig1 in self.lig_set['ligand']:
            
            index_1 = 0
            
            for lig2 in self.lig_set['ligand']:
                
                try:
                    
                    i = self.compute_similarity_index(lig1.ecfp, lig2.ecfp)
                    
                except: 
                    
                    i = 0
                
                self.similarity_matrix.iloc[index_0, index_1] =  i
                
                index_1 += 1
                
            index_0 += 1
            
            
            
    def compute_fsimilarity_matrix(self):
        
        self.fsimilarity_matrix = pd.DataFrame(columns = self.frag_set['fragid'], index = self.frag_set['fragid'])
        
        index_0 = 0
        
        for frag1 in self.frag_set['fragment']:
            
            index_1 = 0
            
            for frag2 in self.frag_set['fragment']:
                
                try:
                    
                    i = self.compute_similarity_index(frag1.ecfp, frag2.ecfp)
                    
                except:
                    
                    i = 0
                
                self.fsimilarity_matrix.iloc[index_0, index_1] = i
                
                index_1 += 1
                
            index_0 += 1
                
                
    def pca_similarity_matrix(self, similarity_matrix):
        
        pca = decomposition.PCA(n_components = 2)
            
        pca.fit(similarity_matrix)
        
        return pd.DataFrame(pca.transform(similarity_matrix))
        
        
        
    def k_means_plot(self, space, clusters = 5, best_cluster = False):
    
        kmeans = cluster.KMeans(init = 'k-means++', n_clusters = clusters, n_init = 10, max_iter = 300)
        
        kmeans.fit(space.values)
        
        return space.plot.scatter(x = 0, y = 1, c = kmeans.labels_, colormap = 'viridis')
        
        #need to assign best cluster ids to the dictionary
        
        #if best_cluster:
            
        #    self.lig_grouped = pd.Series(self.lig_set['ligand']).groupby(kmeans.labels_)
            
        
    def k_means_analysis(self, space, k_min = 2, k_max = 8):
        
        sil_scores = []
        
        for k in range(k_min, k_max):
            
            kmeans = cluster.KMeans(init = 'k-means++', n_clusters = k, n_init = 10, max_iter = 300)
            
            kmeans.fit(space.values)
            
            sil_score = metrics.silhouette_score(space.values, kmeans.labels_)
            
            sil_scores.append(sil_score)
            
        best_score = max(sil_scores)
            
        best_score_index = sil_scores.index(best_score) + k_min
        
        print('\nBest Fit (by silhouette score) found at: '+ str(best_score_index) + ': ' + str(best_score) + '\n')
        
        return self.k_means_plot(space, best_score_index, True)
        
        
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
        
        tmp_frag_set = pd.DataFrame(columns = ['fragment','fragid','ligid','SMILES','parentSMILES'])
        
        lig_count = 0
        
        for lig in self.lig_set['ligand']:
            
            frag_count = 0
            
            for frag in self.single_find_cores(lig.mol, c):
            
                tmp_frag_set = tmp_frag_set.append({'fragment' : AllChem.MolToSmiles(frag),
                                      'ligid' : lig_count, #JUST PASSES INDEX NOT ID
                                      'SMILES' : AllChem.MolToSmiles(frag),
                                      'parentSMILES' : AllChem.MolToSmiles(lig.mol)}, ignore_index = True) #need ligs
                
                frag_count += 1
                
            lig_count += 1
            
        return tmp_frag_set
            
    def sort_cores(self, amb_cores):
        
        unique_cores = amb_cores.groupby("SMILES")
        unique_ligs = amb_cores.groupby("ligid")
       
        #correspondence between unqiue fragid and abmiguous coreids 
        sorted_cores = list(unique_cores.groups.values())

        #correspondence between unique ligid and ambiguous coreids
        sorted_ligs = list(unique_ligs.groups.values())

        f = [ligand(smiles = smi, non_empty_init = True) for smi in list(unique_cores.groups)]
        
        frag_set = pd.DataFrame({'fragment' : f, 'ligs' : None})
        
        frag_set.insert(1,'fragid', [i for i in frag_set.index])
        
        frag_id = 0    
    
        for u_core in sorted_cores:
    
            assoc_ligs = [amb_cores.loc[amb_core, 'ligid'] for amb_core in u_core]  
    
            frag_set.at[frag_id, 'ligs'] = assoc_ligs
    
            for amb_core in u_core:
                amb_cores.at[amb_core, 'fragid'] = frag_id 
        
            frag_id += 1
    
        lig_id = 0

        for u_lig in sorted_ligs:
    
            assoc_cores = [amb_cores.loc[amb_core, 'fragid'] for amb_core in u_lig]
    
            self.lig_set.at[lig_id, 'frags'] = assoc_cores
    
            lig_id += 1  
        
        self.frag_set = frag_set
        
        
## since the lig and frag sets have similar analyses applied to them a number of functions 
## should be move outside the class object and an initiation set of functions to precompute certain values
## needs to be established. After this, neighbourhood visualization can be addressed
        
def create_ligand_graph(lig_set, frag_set):
    
    lig_graph = networkx.Graph()
    
    lig_graph.add_nodes_from(lig_set['ligid'])
    
    for lig_group in frag_set['ligs']:
        
            k = networkx.complete_graph(lig_group)
        
            lig_graph.add_edges_from(k.edges)
        
    return lig_graph


def create_scaffold_graph(lig_set, frag_set):
    
    frag_graph = networkx.Graph()
    
    frag_graph.add_nodes_from(frag_set['fragid'])
    
    for frag_group in lig_set['frags']:
        
        k = networkx.complete_graph(frag_group)
        
        frag_graph.add_edges_from(k.edges)   
    
    return frag_graph