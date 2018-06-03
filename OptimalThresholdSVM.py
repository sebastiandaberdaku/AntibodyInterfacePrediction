# This script computes the best SVM threshold value that maximises the F1 score.
# Please set the validation_set_path variable to the current location of the validation samples before running the script.

import numpy as np
import glob, os
import matplotlib.pyplot as plt
from math import copysign
from sklearn.metrics.classification import precision_recall_fscore_support
from multiprocessing import Process, Manager
import math, scipy


def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    
    return out

validation_set_path = "./structures/"

c = len(validation_set_path)
files = glob.glob(validation_set_path + "*_ab.pdb")

all_scores = []
all_truth = []
thresholds =[]
pdb_ids = []

manager = Manager()

for f in sorted(files) :
    interface_residues = set()

    ab_filename = f
    ab_id = os.path.basename(f)[:-7]
    
    ab_patch_score = "%s%s_ab_patch_score.txt" % (validation_set_path, ab_id)
    ab_patch_truth = "%s%s_ab_patch_truth.txt" % (validation_set_path, ab_id)

    with open(ab_patch_score) as pred, open(ab_patch_truth) as truth :
        patch_pred  = [float(x) for x in pred.readlines()]
        patch_truth = [int(x)   for x in truth.readlines()]
    
    pdb_ids.append(ab_id)

    all_scores.append(patch_pred)
    all_truth.append(patch_truth)

all_thresholds = np.unique(np.concatenate(all_scores))
print len(all_thresholds)
R = 20
pad_size = math.ceil(float(all_thresholds.size)/R)*R - all_thresholds.size
all_thresholds = np.append(all_thresholds, np.zeros(int(pad_size))*np.NaN)

all_thresholds = scipy.nanmean(all_thresholds.reshape(-1,R), axis=1)

print len(all_thresholds)

precision = manager.list([0 for _ in all_thresholds])
recall    = manager.list([0 for _ in all_thresholds])

def compute(indices, all_truth ,all_scores, pdb_ids, all_thresholds, precision, recall):
    print all_thresholds[indices[0]]
    for t in indices :
        for i in xrange(len(pdb_ids)) :
            p, r, _, _ = precision_recall_fscore_support(all_truth[i], [copysign(1, x - all_thresholds[t]) for x in all_scores[i]], average='binary')
            precision[t] += p
            recall[t] += r
        precision[t] /= len(pdb_ids)
        recall[t] /= len(pdb_ids)


# Parallel(n_jobs=12)(delayed(compute)(t, all_truth ,all_scores, pdb_ids, all_thresholds, precision, recall) for t in xrange(len(all_thresholds)))
L = chunkIt(range(len(all_thresholds)), 100)
job = [Process(target=compute, args=(indices, all_truth ,all_scores, pdb_ids, all_thresholds, precision, recall)) for indices in L]
_ = [p.start() for p in job]
_ = [p.join() for p in job]


thresholds_f1scores = [(all_thresholds[t], 2 * precision[t] * recall[t] / (precision[t] + recall[t])) for t in xrange(len(all_thresholds))]


best_pair = max(thresholds_f1scores, key=lambda x:x[1])
print ("Maximum F1 obtained for threshold: %s" % str(best_pair))

plt.figure(2, figsize=(10, 10), dpi=1200)
plt.xlim([all_thresholds[0], all_thresholds[-1]])
plt.ylim([0.0, 1.05])
plt.xlabel('Threshold values')
plt.ylabel('F1 score')
plt.title('Threshold versus F1 scores')
plt.plot(all_thresholds, [a[1] for a in thresholds_f1scores], color='navy', linestyle='solid', linewidth=2)
plt.scatter(best_pair[0], best_pair[1], marker='x', color='red', s=40)
plt.plot([best_pair[0], best_pair[0]], [0, best_pair[1]], linestyle="dotted", linewidth=1, color='red')
plt.plot([all_thresholds[0], best_pair[0]], [best_pair[1], best_pair[1]], linestyle="dotted", linewidth=1, color='red')
plt.annotate("(%.4f, %.4f)" % (best_pair[0], best_pair[1]), xy=(best_pair[0], best_pair[1]), xytext=(-140, 30),
    textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=90,rad=10"))
# plt.legend()
plt.savefig("threshold_for_best_F1_score.pdf", dpi=1200, bbox_inches='tight')
plt.close(2)
plt.clf()

