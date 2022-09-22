import numpy as np
from random import randrange
from sklearn import datasets
from sklearn import metrics
from sklearn.model_selection import cross_val_score
import datatable as dt
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier
import sys
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
import seaborn as sns
from sklearn.model_selection import GridSearchCV,RandomizedSearchCV
from sklearn.metrics import classification_report
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble import BaggingClassifier
import time

class Model:
    def __init__(self, datafile = None,
                 seed=1,scoring = None,n_jobs=None, model_params=None,crossval_n_folds=10,liste_model_name=['ext','knn','rdf','gnb', 'dt', 'nlp', 'svc']
                 ):
        self.scoring = scoring
        self.crossval_n_folds = crossval_n_folds
        self.seed = seed
        self.df = dt.fread(datafile)

        #del self.df[:,"index"]
        #del self.df[:,"seqIdInFam"]

        #self.df = pd.read_csv(datafile)
        #self.df.rename(columns={"familyId": "Classe"})
        #self.df.drop(self.df.columns[self.df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
        #self.df=self.df.astype(int)
        #self.col=self.df.drop(['Classe'], axis=1).columns
        self.n_jobs = n_jobs
        self.liste_model_name=liste_model_name
        self.list_model=[ExtraTreesClassifier(),KNeighborsClassifier(),RandomForestClassifier(),
                         GaussianNB(),
                         DecisionTreeClassifier(),
                         MLPClassifier(),
                         SVC()
                         ]

    #SVC()
    #MLPClassifier()
    #DecisionTreeClassifier()
    #GaussianNB()

    #   'gnb', 'dt', 'nlp', 'svc'


    def split(self, test_size):

        y=np.ravel(self.df[:,"familyId"])
        del self.df[:,"familyId"]
        X=self.df
        X = X.to_numpy()
        #X = np.array(self.df.drop(['Classe'], axis=1))
        #y = np.array(self.df['Classe'])
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size = test_size, random_state = 42)
        #self.X_train = pd.DataFrame(self.X_train)
        #self.X_train = dt.Frame(self.X_train)






    def cv_model(self, last_model_results_=None):


        last_model_results_ = pd.DataFrame(columns = ['model', 'cv_mean', 'cv_std','execution_time'])
        i=0
        for model in self.list_model:
            start_time = time.time()
            kfold = KFold( n_splits=self.crossval_n_folds, random_state=self.seed, shuffle=True)
            cv_scoress = cross_val_score(model, self.X_train, self.y_train, cv = kfold, scoring=self.scoring, n_jobs =self.n_jobs)

            if last_model_results_ is not None:
                last_model_results_ = last_model_results_.append(pd.DataFrame({'model': self.liste_model_name[i],
                                                                               'cv_mean': cv_scoress.mean(),
                                                                               'cv_std': cv_scoress.std(),
                                                                               'execution_time':(time.time() - start_time)},
                                                                              index = [0]),
                                                                 ignore_index = True)
            i=i+1
        return last_model_results_
    def model_test(self, last_model_results_=None):
        start_time = time.time()
        test_result=[]
        for model in self.list_model:
            model.fit(self.X_train, self.y_train)
            test_result.append(model.score(self.X_test,self.y_test))
        return test_result

if __name__ == '__main__':
    #path_data_set=R"C:\Users\ibra\OneDrive - USherbrooke\Project\MotifsExtractionSelection\cmake-build-debug\nbF050_nbSeqs20_40__min_1_max_6_beta_0.csv"
    path_data_set=sys.argv[1]

    print(path_data_set)

    model_instance = Model(datafile=path_data_set,n_jobs=-1)
    model_instance.split(0.3)

    print("train with cross_validate")
    traing_result=model_instance.cv_model()

    plt.figure(figsize=(10,10))
    #la barre pour voir la variation des scores selon les resultat de cross_validation un resultat bon avec minimum de variance est le meilleur
    plt.errorbar(list(traing_result['model']), list(traing_result['cv_mean']),
                 list(traing_result['cv_std']), fmt='o', lw=3,capsize=20,markersize=20)

    print("resultat____test__________________________________")

    test_result=model_instance.model_test()
    #print("les resultats de test de ['ext','knn','rdf', 'gnb', 'dt', 'nlp', 'svc']",test_result)

    traing_result['Score Test'] = test_result

    print(traing_result)

    #plt.bar(['ext','knn','rdf','gnb', 'dt', 'nlp', 'svc'],test_result)

