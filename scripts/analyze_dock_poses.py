'''
Module: AnalyzeDockPoses

Purpose: Functions for analyzing docking results and simplifying the results
'''
from .ligand import *
import numpy as np

def single_min_score(lig, id):
    scrs = []
    res = lig.vina_result_scores[id]
    if True:
        res = res.fillna(0)
        val = res.to_numpy()
        min = 0
        for va in val:
            for v in va:
                if v < min:
                    min = v

        scrs.append(min)

    return scrs

def min_score(lig_set, id):

    scores = []

    for lig in lig_set['ligand']:
        scores.append(single_min_score(lig, id))

    return scores


def compute_rmsd_matrix(lig,
                        rec_id,
                        ):
    poses = lig.vina_result_poses
#Compute RMSD Matrix (static alignment)
#cluster results and chose representative poses

#compute e3fp: pca and k-means?
#compute plif: same?
#identifying interacting residues and collecting those results for distributions, heat maps