'''
Module: AnalyzeDockPoses

Purpose: Functions for analyzing docking results and simplifying the results
'''


def compute_rmsd_matrix(lig,
                        rec_id,
                        ):
    poses = lig.vina_result_poses
#Compute RMSD Matrix (static alignment)
#cluster results and chose representative poses

#compute e3fp: pca and k-means?
#compute plif: same?
#identifying interacting residues and collecting those results for distributions, heat maps