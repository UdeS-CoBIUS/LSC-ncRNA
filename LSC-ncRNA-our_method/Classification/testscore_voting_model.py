


import numpy as np
from random import randrange
from sklearn import datasets
from sklearn import metrics
from sklearn.model_selection import cross_val_score
import datatable as dt
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier,VotingClassifier
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
import joblib


class Model:
    def __init__(self, datafile = None,
                 seed=1,scoring = None,n_jobs=None, model_params=None,crossval_n_folds=10,liste_model_name=['ext','nlp','rdf']
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
        self.list_model=[
            ExtraTreesClassifier(),
            MLPClassifier(max_iter=800),
            RandomForestClassifier()#,
            #SVC(probability=True)
        ]
    #SVC()
    #MLPClassifier()
    #DecisionTreeClassifier()
    #GaussianNB()

    #'dt','gnb','nlp','svc'


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


    def model_voting(self, last_model_results_=None):
        start_time = time.time()

        model1 = self.list_model[0]
        model2 = self.list_model[1]
        model3 = self.list_model[2]
        #model4 = self.list_model[3]


        voted_model = VotingClassifier(estimators=[
            ('ext', model1),
            ('nlp', model2),
            ('rdf', model3)#,
            #('svc', model4)
        ], voting='soft',weights=[1,1,1])
        #eclf1 = eclf1.fit(X, y)
        voted_model.fit(self.X_train,self.y_train)
        #print("la prediction des données:",voted_model.predict(self.X_test))
        print("le score de voted_model",voted_model.score(self.X_test,self.y_test))
        print("total time : ",(time.time() - start_time))
if __name__ == '__main__':
    #data_set=R"C:\Users\ibra\OneDrive - USherbrooke\Project\MotifsExtractionSelection\cmake-build-debug\nbF050_nbSeqs20_40_allF__min_5_max_6_sameFamilyPct_40.csv"
    data_set=sys.argv[1]

    model_instance = Model(datafile=data_set,n_jobs=-1)
    model_instance.split(0.3)

    model_instance.model_voting()
