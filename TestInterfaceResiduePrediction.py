# This script runs the IF algorithm for outlier detection to remove false positive patches and maps the predicted LSPs on the underlying residues.
# The results are compared to other predictor software packages.
# Please remember to set the path variable to the current location of the test set.

import numpy as np
from sklearn.neighbors.kd_tree import KDTree
from glob import glob
from math import copysign
from sklearn.ensemble import IsolationForest
from os.path import basename
from os import path, makedirs

from Bio.PDB.PDBParser import PDBParser
p = PDBParser(QUIET=True, PERMISSIVE=True)
from Bio.PDB.Polypeptide import three_to_one


def convert3to1(s):
    try :
        return three_to_one(s)
    except KeyError : 
        return "X"

import re
_hydrogen = re.compile("[123 ]*H.*")

def isHydrogen(atm):
    return _hydrogen.match(atm.get_id())

def isHETATM(atm):
    return atm.get_parent().get_id()[0] != " "


#######################
# import pickle
#######################

outlier_fraction = 0.18
threshold = 0.6232013
n_iterations = 100
mapping_distance = 6.0 

def compute_average_scores(testing_set_path, prediction_path):
    files = glob("%s/*_ab.pdb" % (testing_set_path))
    
    for pdb_filename in sorted(files) :
        file_id = basename(pdb_filename)[:-7]
        
        pdb_patch_coord = ("%s/%s_ab_patch_centers.txt" % (testing_set_path, file_id))
        pdb_patch_score = ("%s/%s_ab_patch_score.txt" % (testing_set_path, file_id))
      
        with open(pdb_patch_coord) as coord, open(pdb_patch_score) as score:
            patch_coord = [[float(x) for x in a.strip().split()] for a in coord.readlines()]
            patch_score = [float(x) - threshold for x in score.readlines()]
 
        min_v = min(patch_score)
        max_v = max(patch_score)
         
        patch_score_scaled = [(lambda x: -(x / min_v) if x < 0 else (x / max_v))(x) for x in patch_score]
     
        X = np.array([a[0] for a in zip(patch_coord, patch_score_scaled) if a[1] >= 0])
        X_weights = np.array([x for x in patch_score_scaled if x >= 0])
        
        pdb_structure = p.get_structure(file_id, pdb_filename)
        atoms = np.array([atm.get_coord() for atm in pdb_structure.get_atoms() if not isHydrogen(atm)])
        atoms_tree = KDTree(atoms)
        
        residues_coord = {}
        for residue in pdb_structure.get_residues() :
            for atm in residue :
                residues_coord[tuple(atm.get_coord())] = residue
    
        average_residues_scores = {residue : 0 for residue in pdb_structure.get_residues()}

        # since the isollation forest algorithm is random, we run it several times to assess the average performance of the method
        
        
        for iteration in xrange(n_iterations) :
            print "Running iteration %d of %d" % (iteration + 1, n_iterations)
            forest = IsolationForest(contamination=outlier_fraction, n_jobs=-1)
            forest.fit(X, sample_weight=X_weights)
     
            prediction_isolation_forest = forest.predict(patch_coord)
            patch_pred_no_outliers = [copysign(1, x) for x in prediction_isolation_forest]

            # here we map the patch predictions on the underlying residues
            for i in xrange(len(patch_coord)) : # for each patch
                # if it was predicted as non-interface continue to the next
                if patch_pred_no_outliers[i] < 0 : continue 
                # multiple residues can be underneath a given patch, we do not want to consider the same residue more than once
                marked_residues = set() 
                # get all atoms within mapping_distance from the given patch center
                indexes = atoms_tree.query_radius([patch_coord[i]], r=mapping_distance, count_only = False, return_distance=True, sort_results = True)
                for ind in zip(indexes[0][0], indexes[1][0]) :
                    # which residue does the current atom belong to?
                    current_res = residues_coord[tuple(atoms[ind[0]])] 
                    # if already considered continue to the next
                    if current_res in marked_residues : continue 
                    # increase the score of the current residue
                    average_residues_scores[current_res] += 1 / (1.0 + ind[1]) # patch_pred_no_outliers[i] / (1.0 + ind[1])
                    # mark as seen for the current patch
                    marked_residues.add(current_res)
             
        average_residues_scores.update((x, y / n_iterations) for x, y in average_residues_scores.items())
        
        residues_with_scores = [(lambda x, y, z : (convert3to1(z), x[2], x[3][1], x[3][2], y))(residue.get_full_id(), score, residue.get_resname()) for residue, score in average_residues_scores.items()]
        residues_with_scores.sort(key=lambda x : x[2])
        residues_with_scores.sort(key=lambda x : x[1])

        if not path.exists(prediction_path) : makedirs(prediction_path)
        print file_id
        with open("%s/%s_ab_residue_prediction.txt" % (prediction_path, file_id), "wb") as output_residue_scores :
            for r in residues_with_scores :
                output_residue_scores.write("%s;%s;%d;%s;%s\n" %(r[0], r[1], r[2], r[3], str(r[4])))

compute_average_scores("./our_dataset/testing_set/LH_NonProtein/structures/", "./method_comparison/our_method/predictions/LH_NonProtein/")
compute_average_scores("./our_dataset/testing_set/LH_Protein/structures/", "./method_comparison/our_method/predictions/LH_Protein/")
 
compute_average_scores("./our_dataset/homology/LH_NonProtein/structures/", "./method_comparison/our_method_homology/predictions/LH_NonProtein/")
compute_average_scores("./our_dataset/homology_90/LH_Protein/structures/", "./method_comparison/our_method_homology/predictions/LH_Protein/")


