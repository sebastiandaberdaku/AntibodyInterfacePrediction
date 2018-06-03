
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint
from sklearn.svm import SVC, LinearSVC
from sklearn.externals import joblib

from sklearn.model_selection import GroupKFold
from scipy.stats import uniform, rv_continuous

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


class exp_uniform_gen(rv_continuous):
    def rvs(self, loc=0, scale=1, size=1, random_state=None):
        u = uniform(loc=loc, scale=scale)
        return 2**u.rvs(size=size, random_state=random_state)
exp_uniform=exp_uniform_gen(a=0, b=1, name='exp_uniform')

def hyperparameterSearch(training_set_path):
    print ("Importing descriptors from the training set.")
    X, y, labels = loadSamples(training_set_path, "*_ab_train_descriptors_N5.txt", len("_ab_train_descriptors_N5.txt"), 0.05)
    print ("Number of features: %d." % X.shape[-1])
 
    print ("Scaling data.")
    scaler = MinMaxScaler()
    X_scale = scaler.fit_transform(X.todense())
    
    print ("Running randomized hyper-parameter search with 10 fold CV for the RBF kernel.")
    param_dist_rbf = {'kernel': ['rbf'], 'C': exp_uniform(loc=-5, scale=21), 'gamma': exp_uniform(loc=-15, scale=19)}
    random_sv_rbf = RandomizedSearchCV(SVC(), param_distributions=param_dist_rbf, n_iter=200, scoring='roc_auc', cv=GroupKFold(n_splits=10), n_jobs=-1, error_score=0, iid=False, refit=False, verbose=10)
    random_sv_rbf.fit(X_scale, y, groups=labels)
    print "Saving random_sv_rbf model to file."
    joblib.dump(random_sv_rbf, training_set_path + "random_sv_rbf")
     
    random_sv_rbf = joblib.load(training_set_path + "random_sv_rbf")
     
    print ("Running randomized hyper-parameter search with 10 fold CV for the sigmoid kernel.")
    param_dist_sigmoid = {'kernel': ['sigmoid'], 'C': exp_uniform(loc=-5, scale=21), 'coef0': uniform(loc=-2, scale=4), 'gamma': exp_uniform(loc=-15, scale=19)}
    random_sv_sigmoid = RandomizedSearchCV(SVC(), param_distributions=param_dist_sigmoid, n_iter=200, scoring='roc_auc', cv=GroupKFold(n_splits=10), n_jobs=-1, error_score=0, iid=False, refit=False, verbose=10)
    random_sv_sigmoid.fit(X_scale, y, groups=labels)
    print "Saving random_sv_sigmoid model to file."
    joblib.dump(random_sv_sigmoid, training_set_path + "random_sv_sigmoid")
      
#     random_sv_sigmoid = joblib.load(training_set_path + "random_sv_sigmoid")
     
          
    print ("Running randomized hyper-parameter search with 10 fold CV for the linear kernel.")
    param_dist_linear = {'C': exp_uniform(loc=-5, scale=21)}
    random_sv_linear = RandomizedSearchCV(LinearSVC(), param_distributions=param_dist_linear, n_iter=200, scoring='roc_auc', cv=GroupKFold(n_splits=10), n_jobs=-1, error_score=0, iid=False, refit=False, verbose=10)
    random_sv_linear.fit(X_scale, y, groups=labels)
    print "Saving random_sv_linear model to file."
    joblib.dump(random_sv_linear, training_set_path + "random_sv_linear")
      
#     random_sv_linear = joblib.load(training_set_path + "random_sv_linear")
            
    print ("Running randomized hyper-parameter search with 10 fold CV for the polynomial kernel.")
    param_dist_poly = {'kernel': ['poly'], 'C': exp_uniform(loc=-5, scale=21), 'degree': randint(2, 11), 'coef0': uniform(loc=-2, scale=4), 'gamma': exp_uniform(loc=-15, scale=19)}
    random_sv_poly = RandomizedSearchCV(SVC(), param_distributions=param_dist_poly, n_iter=100, scoring='roc_auc', cv=GroupKFold(n_splits=10), n_jobs=-1, error_score=0, iid=False, refit=False, verbose=10)
    random_sv_poly.fit(X_scale, y, groups=labels)
    print "Saving random_sv_poly model to file."
    joblib.dump(random_sv_poly, training_set_path + "random_sv_poly")
    
    with open("best_rbf2_parameters.txt", "w") as best_params :

        print("Best parameters found on training set with the RBF kernel:\n%s %s" % (random_sv_rbf.best_params_, random_sv_rbf.best_score_))
        best_params.write("Best parameters found on training set with the RBF kernel:\n%s %s\n" % (random_sv_rbf.best_params_, random_sv_rbf.best_score_))
        print("kernel = \"%s\"" % (random_sv_rbf.best_params_["kernel"]))
        best_params.write("\nkernel = \"%s\"\n" % (random_sv_rbf.best_params_["kernel"]))
        print("C = %f" % (random_sv_rbf.best_params_["C"]))
        best_params.write("C = %f\n" % (random_sv_rbf.best_params_["C"]))
        print("gamma = %f" % (random_sv_rbf.best_params_["gamma"]))    
        best_params.write("gamma = %f\n" % (random_sv_rbf.best_params_["gamma"]))
        print("Random 10 fold CV scores on development set:")
        best_params.write("Random 10 fold CV scores on development set:\n")
        means = random_sv_rbf.cv_results_['mean_test_score']
        stds = random_sv_rbf.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, random_sv_rbf.cv_results_['params']):
            print("%0.5f (stdev %0.5f) for %r" % (mean, std, params))
            best_params.write("%0.5f (stdev %0.5f) for %r\n" % (mean, std, params))
         
        print("Best parameters found on training set with the linear kernel:\n%s %s" % (random_sv_linear.best_params_, random_sv_linear.best_score_))
        best_params.write("Best parameters found on training set with the linear kernel:\n%s %s\n" % (random_sv_linear.best_params_, random_sv_linear.best_score_))
        print("kernel = \"%s\"" % ('linear'))
        best_params.write("\nkernel = \"%s\"\n" % ('linear'))
        print("C = %f" % (random_sv_linear.best_params_["C"]))
        best_params.write("C = %f\n" % (random_sv_linear.best_params_["C"]))
        print("Random 10 fold CV scores on development set:")
        best_params.write("Random 10 fold CV scores on development set:\n")
        means = random_sv_linear.cv_results_['mean_test_score']
        stds = random_sv_linear.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, random_sv_linear.cv_results_['params']):
            print("%0.5f (stdev %0.5f) for %r" % (mean, std, params))
            best_params.write("%0.5f (stdev %0.5f) for %r\n" % (mean, std, params))
 
        print("Best parameters found on training set with the polynomial kernel:\n%s %s" % (random_sv_poly.best_params_, random_sv_poly.best_score_))
        best_params.write("Best parameters found on training set with the polynomial kernel:\n%s %s\n" % (random_sv_poly.best_params_, random_sv_poly.best_score_))
        print("kernel = \"%s\"" % (random_sv_poly.best_params_["kernel"]))
        best_params.write("\nkernel = \"%s\"\n" % (random_sv_poly.best_params_["kernel"]))
        print("C = %f" % (random_sv_poly.best_params_["C"]))
        best_params.write("C = %f\n" % (random_sv_poly.best_params_["C"]))
        print("gamma = %f" % (random_sv_poly.best_params_["gamma"]))
        best_params.write("gamma = %f\n" % (random_sv_poly.best_params_["gamma"]))
        print("degree = %d" % (random_sv_poly.best_params_["degree"]))
        best_params.write("degree = %d\n" % (random_sv_poly.best_params_["degree"]))
        print("coef0 = %f" % (random_sv_poly.best_params_["coef0"]))
        best_params.write("coef0 = %f\n" % (random_sv_poly.best_params_["coef0"]))
        print("Random 10 fold CV scores on development set:")
        best_params.write("Random 10 fold CV scores on development set:\n")
        means = random_sv_poly.cv_results_['mean_test_score']
        stds = random_sv_poly.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, random_sv_poly.cv_results_['params']):
            print("%0.5f (stdev %0.5f) for %r" % (mean, std, params))
            best_params.write("%0.5f (stdev %0.5f) for %r\n" % (mean, std, params))
 
        print("Best parameters found on training set with the sigmoid kernel:\n%s %s" % (random_sv_sigmoid.best_params_, random_sv_sigmoid.best_score_))
        best_params.write("Best parameters found on training set with the sigmoid kernel:\n%s %s\n" % (random_sv_sigmoid.best_params_, random_sv_sigmoid.best_score_))
        print("kernel = \"%s\"" % (random_sv_sigmoid.best_params_["kernel"]))
        best_params.write("\nkernel = \"%s\"\n" % (random_sv_sigmoid.best_params_["kernel"]))
        print("C = %f" % (random_sv_sigmoid.best_params_["C"]))
        best_params.write("C = %f\n" % (random_sv_sigmoid.best_params_["C"]))
        print("gamma = %f" % (random_sv_sigmoid.best_params_["gamma"]))
        best_params.write("gamma = %f\n" % (random_sv_sigmoid.best_params_["gamma"]))
        print("coef0 = %f" % (random_sv_sigmoid.best_params_["coef0"]))
        best_params.write("coef0 = %f\n" % (random_sv_sigmoid.best_params_["coef0"]))
        print("Random 10 fold CV scores on development set:")
        best_params.write("Random 10 fold CV scores on development set:\n")
        means = random_sv_sigmoid.cv_results_['mean_test_score']
        stds = random_sv_sigmoid.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, random_sv_sigmoid.cv_results_['params']):
            print("%0.5f (stdev %0.5f) for %r" % (mean, std, params))
            best_params.write("%0.5f (stdev %0.5f) for %r\n" % (mean, std, params))
    
print ("Running hyper-parameter search")
training_set_path = "./our_dataset/training_set/structures/"
hyperparameterSearch(training_set_path)


