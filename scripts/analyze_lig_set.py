from .screen_library import *
from .util import *
import pandas as pd
from sklearn import decomposition
from sklearn import cluster
from sklearn import metrics
from rdkit.Chem import Recap
from itertools import compress
import networkx
import numpy as np


class AnalyzeLigSet():

    def __init__(self, lig_set=None, full_init=False, pca_plot=1, frag_plot=1):

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
        self.similarity_matrix = pd.DataFrame(columns=self.lig_set['ligid'], index=self.lig_set['ligid'])
        self.lig_space2 = None

        self.frag_set = pd.DataFrame(columns=['fragment', 'fragid', 'ligs'])

        self.fsimilarity_matrix = None
        self.frag_space2 = None

        self.lig_graph = None
        self.frag_graph = None

        # will likeley be replaced by neighbourhood and cluster searching
        self.lig_grouped = None

        self.frag_grouped = None

        if full_init:
            self.full_init()

        if pca_plot is not None:
            if int(pca_plot) == 0:
                pca_plot = 1
            try:
                self.compute_similarity_matrix()
                self.lig_space2 = self.pca_similarity_matrix(self.similarity_matrix)
                clusters = self.k_means_lables(self.lig_space2, clusters=int(pca_plot))
                self.pca_plot = self.scatter_to_js(self.lig_space2, clusters=clusters)
            except:
                self.pca_plot = '{x:0,y:0}'

        if frag_plot is not None:
            if int(frag_plot) == 0:
                frag_plot = 1

            self.sort_cores(self.find_cores(c=2 / 3))
            self.compute_fsimilarity_matrix()
            self.frag_space2 = self.pca_similarity_matrix(self.fsimilarity_matrix)
            clusters = self.k_means_lables(self.frag_space2, clusters=int(frag_plot))
            self.fpca_plot = self.scatter_to_js(self.frag_space2, clusters=clusters)

            try:
                pass
            except:
                self.fpca_plot = '{x:0,y:0}'

    def full_init(self):

        try:
            # 1
            self.compute_similarity_matrix()

            # 2
            self.lig_space2 = self.pca_similarity_matrix(self.similarity_matrix)

            # 3
            self.k_means_analysis(self.lig_space2)

            # 4
            self.sort_cores(self.find_cores(c=2 / 3))

            # 5
            self.compute_fsimilarity_matrix()

            # 6
            self.frag_space2 = self.pca_similarity_matrix(self.fsimilarity_matrix)

            # 7
            self.k_means_analysis(self.frag_space2)

            # 8
            self.lig_graph = self.create_ligand_graph(self.lig_set, self.frag_set)
            self.frag_graph = self.create_scaffold_graph(self.lig_set, self.frag_set)

        except:
            raise LigSetAnalysisError()

        '''
        Now you can set a number of things to be pre-computed (decided by full_init):
        1. The ligand similarity matrix
        2. The 2-D reduced similarity matrix
        3. The k-means grouping
        4. The unique scaffolds (fragments)
        5. The the fsimilarity matrix
        6. The 2-D reduced fsimilarity matrix
        7. The k-means grouping
        8. The ligand space and fragments space graphs
        '''

    ### Lipinski Type Values ###

    def logp_spread(self):

        spread = pd.Series([lig.logp for lig in self.lig_set['ligand']])
        spread.hist()

    def mmass_spread(self):

        spread = pd.Series([lig.mmass for lig in self.lig_set['ligand']])
        spread.hist()

    ### PCA Analysis of similarity matrix and k-means clustering ###

    def compute_similarity_index(self, fp1=None, fp2=None):

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

                self.similarity_matrix.iloc[index_0, index_1] = i

                index_1 += 1

            index_0 += 1

    def compute_fsimilarity_matrix(self):

        self.fsimilarity_matrix = pd.DataFrame(columns=self.frag_set['fragid'], index=self.frag_set['fragid'])

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

        pca = decomposition.PCA(n_components=2)

        pca.fit(similarity_matrix)

        return pd.DataFrame(pca.transform(similarity_matrix))

    def k_means_plot(self, space, clusters=5, best_cluster=False):

        kmeans = cluster.KMeans(init='k-means++', n_clusters=clusters, n_init=10, max_iter=300)

        kmeans.fit(space.values)

        return space.plot.scatter(x=0, y=1, c=kmeans.labels_, colormap='viridis')

        # need to assign best cluster ids to the dictionary

        # if best_cluster:

        #    self.lig_grouped = pd.Series(self.lig_set['ligand']).groupby(kmeans.labels_)

    def k_means_lables(self, space, clusters=5):
        kmeans = cluster.KMeans(init='k-means++', n_clusters=clusters, n_init=10, max_iter=300)

        kmeans.fit(space.values)

        return kmeans.labels_

    def k_means_analysis(self, space, k_min=2, k_max=8):

        sil_scores = []

        for k in range(k_min, k_max):
            kmeans = cluster.KMeans(init='k-means++', n_clusters=k, n_init=10, max_iter=300)

            kmeans.fit(space.values)

            sil_score = metrics.silhouette_score(space.values, kmeans.labels_)

            sil_scores.append(sil_score)

        best_score = max(sil_scores)

        best_score_index = sil_scores.index(best_score) + k_min

        print('\nBest Fit (by silhouette score) found at: ' + str(best_score_index) + ': ' + str(best_score) + '\n')

        return self.k_means_plot(space, best_score_index, True)

    # need functions that encompass a series of analytical functions (compute, pca, k-means) and
    # produce a single graphic at the end.

    ### Extraction of Cores / Scaffold Analysis ###

    ## Adapted from the code provided by Medina-Franco and Naveja (2019) ##

    def single_find_cores(self, mol, c=2 / 3):

        smi = AllChem.MolToSmiles(mol)
        mol_h = mol.GetNumHeavyAtoms()

        cores_full = list(Recap.RecapDecompose(mol).GetAllChildren().keys()) + [smi]
        cores_comp = list(
            compress(cores_full, [AllChem.MolFromSmiles(core).GetNumHeavyAtoms() >= c * mol_h for core in cores_full]))

        frags = [AllChem.MolFromSmiles(core) for core in cores_comp]

        return frags

    def find_cores(self, c=2 / 3):

        tmp_frag_set = pd.DataFrame(columns=['fragment', 'fragid', 'ligid', 'SMILES', 'parentSMILES'])

        lig_count = 0

        for lig in self.lig_set['ligand']:

            frag_count = 0

            for frag in self.single_find_cores(lig.mol, c):
                tmp_frag_set = tmp_frag_set.append({'fragment': AllChem.MolToSmiles(frag),
                                                    'ligid': lig_count,  # JUST PASSES INDEX NOT ID
                                                    'SMILES': AllChem.MolToSmiles(frag),
                                                    'parentSMILES': AllChem.MolToSmiles(lig.mol)},
                                                   ignore_index=True)  # need ligs

                frag_count += 1

            lig_count += 1

        return tmp_frag_set

    def sort_cores(self, amb_cores):

        unique_cores = amb_cores.groupby("SMILES")
        unique_ligs = amb_cores.groupby("ligid")

        # correspondence between unqiue fragid and abmiguous coreids
        sorted_cores = list(unique_cores.groups.values())

        # correspondence between unique ligid and ambiguous coreids
        sorted_ligs = list(unique_ligs.groups.values())

        f = [Ligand(smiles=smi) for smi in list(unique_cores.groups)]

        frag_set = pd.DataFrame({'fragment': f, 'ligs': None})

        frag_set.insert(1, 'fragid', [i for i in frag_set.index])

        frag_id = 0

        for u_core in sorted_cores:

            assoc_ligs = [amb_cores.loc[amb_core, 'ligid'] for amb_core in u_core]

            # frag_set.loc[frag_id, 'ligs'] = assoc_ligs

            for amb_core in u_core:
                amb_cores.loc[amb_core, 'fragid'] = frag_id

            frag_id += 1

        lig_id = 0

        for u_lig in sorted_ligs:
            assoc_cores = [amb_cores.loc[amb_core, 'fragid'] for amb_core in u_lig]

            # self.lig_set.loc[lig_id, 'frags'] = assoc_cores

            lig_id += 1

        self.frag_set = frag_set

    ## since the lig and frag sets have similar analyses applied to them a number of functions
    ## should be move outside the class object and an initiation set of functions to precompute certain values
    ## needs to be established. After this, neighbourhood visualization can be addressed

    def create_ligand_graph(self, lig_set, frag_set):

        lig_graph = networkx.Graph()

        lig_graph.add_nodes_from(lig_set['ligid'])

        for lig_group in frag_set['ligs']:
            k = networkx.complete_graph(lig_group)

            lig_graph.add_edges_from(k.edges)

        return lig_graph

    def create_scaffold_graph(self, lig_set, frag_set):

        frag_graph = networkx.Graph()

        frag_graph.add_nodes_from(frag_set['fragid'])

        for frag_group in lig_set['frags']:
            k = networkx.complete_graph(frag_group)

            frag_graph.add_edges_from(k.edges)

        return frag_graph

    # helpful functions (might move to util file later, but this currently feels most applicable under the current class

    def scatter_to_js(self, scatter, clusters):
        # assuming the structure of the plot is [[x vals],[y vals]]
        cluster_colors = ['#CFFBE2', '#AFDAF7', '#EFE583', '#B5AE4F',
                          '#F9CD9B', '#D3CDFF', '#E1D776', '#F2A694', ]
        scatter = np.array(scatter)
        scatter_js = '{ datasets: [ '
        for c in range(max(clusters) + 1):
            scatter_js += '{'
            scatter_js += f' label: "Cluster {c}", data: ['
            for i in range(len(scatter)):
                if clusters[i] == c:
                    scatter_js += str('{ x: '
                                      + str(scatter[i][0])
                                      + ', y: '
                                      + str(scatter[i][1])
                                      + '},\n')
            scatter_js += f'], backgroundColor : "{cluster_colors[c]}"'
            scatter_js += '},'
        scatter_js += '],}'
        return scatter_js
'''
graph 1: in ligand space
Nodes: unique ligands
Edges: when two ligands share a scaffold

Development: make edges when two ligands share a scaffold of similarity index s
(where s = 1 is the original graph)

graph 2: in scaffold space
Nodes: unique scaffolds
Edges: when two scaffolds share a ligand

Development: make edges when two scaffolds share a ligand of similarity index s
(again s = 1 is the original graph)


Notice that for each set (no matter the size) of ligs or frags in the ligand_set and fragment_set 
the connectivity of a given node is determined by selecting those nodes as a subgraph then making that
the complete subgraph with n nodes.


For the developemt ideas:
1. Select a fragment
2. Identify assocaited ligids
3. Identify all memebers of the corresponding row of the frag or lig similiarity matrix with value >= s
4. Create a complete subgraph of all these nodes

'''
'''

for better frag space analysis:
1. Eliminate fragments with only 1 ligand
2. Eliminate fragments with the same associated ligands but smaller

'''
