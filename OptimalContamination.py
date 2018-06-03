# This script computes the best contamination value for the Isolation Forest algorithm that maximises the F1 score.
# Please set the validation_set_path variable to the current location of the validation samples before running the script.

import glob, os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import IsolationForest
from sklearn.metrics.classification import precision_recall_fscore_support
from multiprocessing import Manager
from joblib import Parallel, delayed

n_iterations = 100

THRESHOLD = 0.6232013

validation_set_path = "./structures/"
c = len(validation_set_path)
files = glob.glob(validation_set_path + "*_ab.pdb")

manager = Manager()


ab_pred = []
ab_truth = []
ab_coord = []
ab_X = []
ab_X_weights = []
pdb_ids = []
           
for f in sorted(files) :
    ab_filename = f
    ab_id = os.path.basename(f)[:-7]
    
    pdb_ids.append(ab_id)
    
    ab_patch_coord = "%s%s_ab_patch_centers.txt" % (validation_set_path, ab_id)
    ab_patch_score = "%s%s_ab_patch_score.txt" % (validation_set_path, ab_id)
    ab_patch_truth = "%s%s_ab_patch_truth.txt" % (validation_set_path, ab_id)

    with open(ab_patch_coord) as coord, open(ab_patch_score) as pred, open(ab_patch_truth) as truth :
        patch_pred   = [(float(x) - THRESHOLD) for x in pred.readlines()]
        patch_truth  = [int(x)   for x in truth.readlines()]
        patch_coord  = [[float(x) for x in a.split()] for a in coord.readlines()]
    
    min_v = min(patch_pred)
    max_v = max(patch_pred)
    
    patch_pred_scaled = [(lambda x: -(x / min_v) if x < 0 else (x / max_v))(x) for x in patch_pred]

    X = np.array([a[0] for a in zip(patch_coord, patch_pred) if a[1] >= 0])
    X_weights = np.array([x for x in patch_pred_scaled if x >= 0])
    
    ab_X.append(X)
    ab_X_weights.append(X_weights)

    ab_pred.append(patch_pred)
    ab_truth.append(patch_truth)
    ab_coord.append(patch_coord)

outlier_fractions = list(np.arange(0.01, 0.51, 0.01))

precision = manager.list([0 for _ in outlier_fractions])
recall    = manager.list([0 for _ in outlier_fractions])

def compute_scores(o, n_iterations, pdb_ids, ab_truth, ab_coord, ab_X, ab_X_weights, precision, recall):
    print outlier_fractions[o]
    forest = IsolationForest(contamination=outlier_fractions[o], n_jobs=4)
    for i in xrange(len(pdb_ids)) :
        print pdb_ids[i]
        current_precision = 0
        current_recall = 0
        for _ in xrange(n_iterations) :
            forest.fit(ab_X[i], sample_weight=ab_X_weights[i])
            patch_pred_no_outliers = forest.predict(ab_coord[i])
            p, r, _, _ = precision_recall_fscore_support(ab_truth[i], patch_pred_no_outliers, average='binary')
            current_precision += p
            current_recall += r
        current_precision /= n_iterations
        current_recall /= n_iterations
        precision[o] += current_precision
        recall[o] += current_recall
    precision[o] /= len(pdb_ids)
    recall[o] /= len(pdb_ids)

Parallel(n_jobs=12, verbose=5)(delayed(compute_scores)(o, n_iterations, pdb_ids, ab_truth, ab_coord, ab_X, ab_X_weights, precision, recall) for o in xrange(len(outlier_fractions)))


       
f1_mean = [2 * precision[o] * recall[o] / (precision[o] + recall[o])  for o in xrange(len(outlier_fractions))]

outlier_fractions.insert(0, 0)

f1_mean.insert(0, 0.5802242055419393)

print "outlier_fractions = %s" % outlier_fractions
print "f1_mean = %s" % f1_mean 

best_pair = max(zip(outlier_fractions, f1_mean), key=lambda x:x[1])

plt.figure(figsize=(10, 10), dpi=1200)
plt.xlim([0.0, 0.5])
plt.ylim([0.0, 1.0])
plt.xlabel('Outlier fraction')
plt.ylabel('Average F1 score') 
plt.title('The effect of the outlier fraction parameter on the average F1 score \nafter applying the Isolation Forest algorithm')
plt.plot(outlier_fractions, f1_mean, color='navy' , linestyle='-', linewidth=1)

plt.scatter(best_pair[0], best_pair[1], marker='x', color='red', s=40)
plt.plot([best_pair[0], best_pair[0]], [0, best_pair[1]], linestyle="dotted", linewidth=1, color='red')
plt.plot([0, best_pair[0]], [best_pair[1], best_pair[1]], linestyle="dotted", linewidth=1, color='red')
plt.annotate("(%.2f, %.4f)" % best_pair, xy=best_pair, xytext=(-140, 30),
    textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.3"))

plt.savefig("best_outlier_f1_score_all_th.pdf", dpi=1200, bbox_inches='tight')
plt.clf()
plt.close()
