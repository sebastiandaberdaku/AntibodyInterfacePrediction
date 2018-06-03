from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVC
from sklearn.ensemble import BaggingClassifier

from sklearn.externals import joblib

from glob import glob
from os.path import basename
from numpy import append
from scipy.sparse import vstack
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_svmlight_file

def loadSamples(folder_path, file_wildcard, suffix_length, size=None):
    ''' Method that loads all the samples from the files in folder_path whose filename corresponds to file_wildcard.
        The size parameter can be used to load a uniformly sampled subset of the overall samples in the current folder.
        Parameters
        ----------
        folder_path : string
                    path to the folder that contains the sample files
        file_wildcard : string
                    file name wildcard
        size : float, int, or None, default None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.
        
        Returns
        -------
        samples : tuple (samples, class labels, file labels)
                Returns the loaded samples with relative class labels and 
                the file labels (PDB ids).        
    '''
    X = y = None
    labels = []
    files = glob(folder_path + file_wildcard)
    for f in sorted(files):
        print "Loading %s..." % basename(f) 
        X_f, y_f = load_svmlight_file(f, zero_based=False)
        if size is not None :
            X_f, _, y_f, _ = train_test_split(X_f, y_f, train_size=size, stratify = y_f)
        if X is None:
            X = X_f
            y = y_f
        else :
            X = vstack([X, X_f], "csr")
            y = append(y, y_f)
        current_label = basename(f)[:-suffix_length]
        labels += [current_label] * y_f.size
    return (X, y, labels)


training_set_path = "./our_dataset/training_set/structures/"

print ("Importing descriptors from the training set.")
X_train, y_train, labels_train = loadSamples(training_set_path, "*_ab_train_descriptors_N5.txt", len("_ab_train_descriptors_N5.txt"))
print ("Number of features: %d." % X_train.shape[-1])
 
print ("Scaling data.")
scaler = MinMaxScaler()
X_train_scale = scaler.fit_transform(X_train.todense())
 
print "Saving scaler model to file."
joblib.dump(scaler, "%s/scaler_object" % training_set_path)
# scaler = joblib.load("%s/scaler_object" % training_set_path)

kernel = "rbf"
C = 540.2102985534485
gamma = 0.007983367409838158

model = BaggingClassifier(SVC(kernel=kernel, C=C, gamma=gamma, cache_size = 500), max_samples=1.0 / 12, n_estimators=12, n_jobs=-1)
model.fit(X_train_scale, y_train)
 
print "Saving SVM model to file."
joblib.dump(model, "%s/svm_model" % training_set_path)
# model = joblib.load("%s/svm_model" % training_set_path)

for test_set_path in ["./our_dataset/testing_set/LH_Protein/structures/", "./our_dataset/testing_set/LH_NonProtein/structures/", "./our_dataset/validation_set/structures/", "./our_dataset/homology/LH_Protein/structures/", "./our_dataset/homology/LH_NonProtein/structures/"] : 
    print ("Importing descriptors from the testing set %s." % test_set_path)
    X_test, y_test, labels_test = loadSamples(test_set_path, "*_ab_test_descriptors_N5.txt", len("_ab_test_descriptors_N5.txt"))
    print ("Number of features: %d." % X_test.shape[-1])
    X_test_scale = scaler.transform(X_test.todense())
        
    print "Predicting the testing set %s." % test_set_path
    y_score = model.decision_function(X_test_scale)
    
    get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
    
    pdb_ids = sorted(set(labels_test))
    for file_id in pdb_ids :
        pdb_id_indices = get_indexes(file_id, labels_test)
        with open("%s/%s_ab_patch_score.txt" % (test_set_path, file_id), "w") as out_scores :
            for p in y_score[pdb_id_indices] :
                out_scores.write("%f\n" % p)

