import os
import time
## import fcntl # this is only on linux and mac os, not on windows. we use portalocker instead
import portalocker # instead of fcntl, we use portalocker which is cross-platform compatible
import pandas as pd
import datatable as dt
import numpy as np
from Bio import SeqIO

from sklearn.ensemble import ExtraTreesClassifier, VotingClassifier, RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, classification_report
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import confusion_matrix

import construct_test_matrix as ctm

class Model:
    def __init__(self, n_job=1, train_csv_path=None):

        self.train_csv_path = train_csv_path
        self.model_name = None

        local_n_estimators = 700
        self.models = {
            'EXT': ExtraTreesClassifier(n_estimators=local_n_estimators, n_jobs=n_job),
            'MLP': MLPClassifier(max_iter=800),
            'RDF': RandomForestClassifier(n_estimators=local_n_estimators, n_jobs=n_job),
            'VOT': VotingClassifier(estimators=[
                ('EXT', ExtraTreesClassifier(n_estimators=local_n_estimators, n_jobs=n_job)),
                ('MLP', MLPClassifier(max_iter=800)),
                ('RDF', RandomForestClassifier(n_estimators=local_n_estimators, n_jobs=n_job))
            ], voting='soft', weights=[1, 1, 1], n_jobs=n_job)
        }

        self.dir_test_files = None # directory of the test files
        self.test_file_ext = None # extension of the test files

        self.full_test_name: str = None # full name of the test, will be set later in train()
        self.time_read_csv_datatable = 0 # time to read the csv file with datatable
        self.time_fit_model = 0 # time to fit the model
        self.total_time_train = 0 # total time train : read csv + fit model + save motifs
        self.score_train = 0
        self.train_pred_score = 0
        self.total_time_test = 0 # total time test = time read seqs + time matrix construction + time prediction
        self.accuracy = 0 # accuracy score test
        self.precision = 0 # precision score test
        self.recall = 0 # recall score test
        self.f1 = 0 # f1-score score test
        self.number_test_seqs = 0 # number of total sequences in the test set
        self.time_test_matrix_construction = 0 # during the test, the matrix (vectors) is constructed from the test input sequences: nbOccurrences of each motif in each sequence.
        self.time_test_prediction_only = 0 # during the test, the prediction is made by the model
        self.test_result_pred = None # during the test, the prediction is made by the model
        self.list_all_classes_ids = None # list of all classes ids in the test set
        self.list_all_classes_names = None # list of all classes names in the test set

    def _get_model(self, model_name):
            if model_name not in self.models:
                raise ValueError(f"Unknown model name: {model_name}")
            return self.models[model_name]

    def train(self, model_name):
            if model_name not in self.models:
                raise ValueError(f"Unknown model name: {model_name}. Available models are: {', '.join(self.models.keys())}")
            self.model_name = model_name
            self.full_test_name = f"{self.model_name}_{os.path.basename(self.train_csv_path)}"
            self._train_with_dt(self._get_model(model_name))

            
    def test(self, test_dir, file_ext):
        if not self.model_name:
            raise ValueError("Model has not been trained yet. Call train() first.")
        
        if self.model_name not in self.models:
            raise ValueError(f"Unknown model name: {self.model_name}. Available models are: {', '.join(self.models.keys())}")
        
        self.dir_test_files = test_dir # to be used later in _test(), and in print_detailed_results()
        self.test_file_ext = file_ext # to be used later in _test(), and in print_detailed_results()

        model = self._get_model(self.model_name)
        
        # check if the model is trained
        # this can be checked by sklearn.utils.validation.check_is_fitted
        # or see : https://stackoverflow.com/questions/39884009/whats-the-best-way-to-test-whether-an-sklearn-model-has-been-fitted
        if not hasattr(model, 'is_fitted_') or not model.is_fitted_:
            raise ValueError("Model is not trained. Call train() first.")
        
        self._test(model)

    def _train_with_dt(self, model):
        start_0_time = time.time()
        data_arn = dt.fread(self.train_csv_path)
        self.time_read_csv_datatable = time.time() - start_0_time

        data_classe = np.ravel(data_arn[:, "familyId"])
        del data_arn[:, "familyId"]

        start_time_train = time.time()
        model.fit(data_arn, data_classe)
        self.time_fit_model = time.time() - start_time_train
     
        self.list_motifs = data_arn.names  # save list of used motifs in classification
        
        # Total time train : read csv + train + save motifs
        self.total_time_train = time.time() - start_0_time
        
        # statistics
        self.score_train = model.score(data_arn, data_classe) # Calculate accuracy score on training data
        pred_train = model.predict(data_arn) # Make predictions on training data for further evaluation
        self.train_pred_score = accuracy_score(data_classe, pred_train) # Calculate accuracy score on training data predictions




    def train_test_from_one_CSV_file_pd(self, csv_file_path, test_size):
        start_time = time.time()
        data_arn = pd.read_csv(csv_file_path)
        data_arn = data_arn.rename(columns={"familyId": "Classe"})
        data_arn = data_arn.drop(['seqIdInFam', 'index'], axis=1)
        data_arn.drop(data_arn.filter(regex="Unname"), axis=1, inplace=True)
        data_arn = data_arn.astype(int)

        data_classe = data_arn["Classe"]
        data_motif = data_arn.drop(['Classe'], axis=1)

        X_train, X_test, y_train, y_test = train_test_split(data_motif, data_classe, test_size=test_size)

        print("len(y_test = ", len(y_test), "\n")

        print(" y_test  = \n", y_test)

        # ExTrCl = ExtraTreesClassifier(n_estimators=10, criterion='gini', max_depth=10)
        ExTrCl = ExtraTreesClassifier()

        ExTrCl.fit(X_train, y_train)

        end_time = time.time()

        print(" Train time = ", end_time - start_time)

        score_train = ExTrCl.score(data_motif, data_classe)
        pred_train = ExTrCl.predict(X_train)

        start_time = time.time()
        pred_test = ExTrCl.predict(X_test)
        end_time = time.time()
        print(" Test time = ", end_time - start_time)

        train_score, test_score = accuracy_score(y_train, pred_train), accuracy_score(y_test, pred_test)

        print(" score_train = ", score_train, " | Train score = ", train_score, " | Test score = ", test_score)
        print(" Len list motifs : ", len(list(data_motif)))

    def test_create_dt_df(self):
        lf = ['f1', 'f2', 'f3', 'f4', 'f5']
        list_lists = [
            [0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2],
        ]

        matrix = np.array(list_lists)

        print(matrix)

        dt_df = dt.Frame(matrix, names=lf)
        print(dt_df)

        pd_df = pd.DataFrame(list_lists, columns=lf)

        print(pd_df)

        pd_df_from_dt_df = dt_df.to_pandas()

        print(" dose pd_df == dt_df ? : ", pd_df == pd_df_from_dt_df)

    def train_pd(self, csv_file_path):
        start_time = time.time()
        data_arn = pd.read_csv(csv_file_path)  # On importe la dataset
        end_time = time.time()
        print(" time Read_csv file : ", end_time - start_time, " s")
        data_arn = data_arn.rename(columns={"familyId": "Classe"})
        data_arn = data_arn.drop(['seqIdInFam', 'index'], axis=1)
        data_arn.drop(data_arn.filter(regex="Unname"), axis=1, inplace=True)
        data_arn = data_arn.astype(int)  # On change le type to int

        data_classe = data_arn["Classe"]
        data_motif = data_arn.drop(['Classe'], axis=1)

        start_time_train = time.time()
        self.ExTrCl.fit(data_motif, data_classe)
        end_time = time.time()

        print(" train only time : ", end_time - start_time_train, " s")

        self.list_motifs = list(data_motif)  # save list of use motifs in classification

        # statistics
        score_train = self.ExTrCl.score(data_motif, data_classe)
        pred_train = self.ExTrCl.predict(data_motif)
        train_pred_score = accuracy_score(data_classe, pred_train)

        print("all Time Train = ", end_time - start_time, " s")
        print(" score_train = ", score_train, " | Train pred score = ", train_pred_score)
        print(" Len list motifs : ", len(list(data_motif)))
        # print(" data_motif.info() = ",data_motif.info())

    def test_one_pd(self, test_vec):
        dftest = pd.DataFrame([test_vec], columns=self.list_motifs)

        print(dftest)

        result = self.ExTrCl.predict(dftest)
        print(result)

    def test_groupe_pd(self, test_matrix):
        dftest = pd.DataFrame(test_matrix, columns=self.list_motifs)

        print(dftest)

        result = self.ExTrCl.predict(dftest)
        print(result)

    def test_groupe_score_pd(self, test_matrix, list_classes):
        start_time_0 = time.time()
        dftest = pd.DataFrame(test_matrix, columns=self.list_motifs)
        end_time = time.time()
        print(" time creating DataFrame = ", end_time - start_time_0)

        start_time = time.time()
        result = self.ExTrCl.predict(dftest)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")
        print(" Time create + pred = ", end_time - start_time_0, " s")

        score_test = accuracy_score(list_classes, result)

        print(" score test = ", score_test)

    def train_test_from_one_CSV_file_dt(self, csv_file_path, test_size):
        big_start_time = time.time()
        data_arn = dt.fread(csv_file_path)

        del data_arn[:, "seqIdInFam"]
        del data_arn[:, 'index']
        del data_arn[:, 'C0']

        data_classe = np.ravel(data_arn[:, "familyId"])

        del data_arn[:, "familyId"]

        data_motif = data_arn.to_numpy()

        X_train, X_test, y_train, y_test = train_test_split(data_motif, data_classe, test_size=test_size)

        ExTrCl = ExtraTreesClassifier()

        # X_train = dt.Frame(X_train)
        # print(" Len list motifs : ", len(X_train.names))
        # return

        start_time = time.time()
        ExTrCl.fit(X_train, y_train)
        end_time = time.time()

        print(" Only Train time = ", end_time - start_time)
        print(" all time till after traing = ", end_time - big_start_time)

        score_train = ExTrCl.score(data_motif, data_classe)
        pred_train = ExTrCl.predict(X_train)

        start_time = time.time()
        pred_test = ExTrCl.predict(X_test)
        end_time = time.time()
        print(" Test time = ", end_time - start_time)

        train_score, test_score = accuracy_score(y_train, pred_train), accuracy_score(y_test, pred_test)

        print(" score_train = ", score_train, " | Train score = ", train_score, " | Test score = ", test_score)
        # print(" Len list motifs : ",len(data_motif.names))

    def ext_train_with_dt(self, csv_file_path):
        start_0_time = time.time()
        data_arn = dt.fread(csv_file_path)
        end_time = time.time()
        print(" time Read_csv file : ", end_time - start_0_time, " s")

        # start_time = time.time()
        data_classe = np.ravel(data_arn[:, "familyId"])
        # end_time = time.time()
        # print(" time get class array = ", end_time - start_time)

        del data_arn[:, "familyId"]

        start_time_train = time.time()
        self.ExTrCl.fit(data_arn, data_classe)
        end_time = time.time()

        print(" train only time : ", end_time - start_time_train, " s")

        self.list_motifs = data_arn.names  # save list of use motifs in classification

        # statistics
        self.score_train = self.ExTrCl.score(data_arn, np.ravel(data_classe))
        pred_train = self.ExTrCl.predict(data_arn)
        self.train_pred_score = accuracy_score(np.ravel(data_classe), pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", self.score_train, " | Train pred score = ", self.train_pred_score)
        print(" Len list motifs : ", len(list(data_arn)))
        print("list motifs len : {}".format(len(self.list_motifs)))

        self.total_time_train = end_time - start_0_time

    def nlp_train_with_dt(self, csv_file_path):
        start_0_time = time.time()
        data_arn = dt.fread(csv_file_path)
        end_time = time.time()
        print(" time Read_csv file : ", end_time - start_0_time, " s")

        # start_time = time.time()
        data_classe = np.ravel(data_arn[:, "familyId"])
        # end_time = time.time()
        # print(" time get class array = ", end_time - start_time)

        del data_arn[:, "familyId"]

        start_time_train = time.time()
        self.nlp.fit(data_arn, data_classe)
        end_time = time.time()

        print(" train only time : ", end_time - start_time_train, " s")

        self.list_motifs = data_arn.names  # save list of use motifs in classification

        # statistics
        self.score_train = self.nlp.score(data_arn, np.ravel(data_classe))
        pred_train = self.nlp.predict(data_arn)
        self.train_pred_score = accuracy_score(np.ravel(data_classe), pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", self.score_train, " | Train pred score = ", self.train_pred_score)
        print(" Len list motifs : ", len(list(data_arn)))
        print("list motifs len : {}".format(len(self.list_motifs)))

        self.total_time_train = end_time - start_0_time


    def rdf_train_with_dt(self, csv_file_path):
        start_0_time = time.time()
        data_arn = dt.fread(csv_file_path)
        end_time = time.time()
        print(" time Read_csv file : ", end_time - start_0_time, " s")

        # start_time = time.time()
        data_classe = np.ravel(data_arn[:, "familyId"])
        # end_time = time.time()
        # print(" time get class array = ", end_time - start_time)

        del data_arn[:, "familyId"]

        start_time_train = time.time()
        self.rdf.fit(data_arn, data_classe)
        end_time = time.time()

        print(" train only time : ", end_time - start_time_train, " s")

        self.list_motifs = data_arn.names  # save list of use motifs in classification

        # statistics
        self.score_train = self.rdf.score(data_arn, np.ravel(data_classe))
        pred_train = self.rdf.predict(data_arn)
        self.train_pred_score = accuracy_score(np.ravel(data_classe), pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", self.score_train, " | Train pred score = ", self.train_pred_score)
        print(" Len list motifs : ", len(list(data_arn)))
        print("list motifs len : {}".format(len(self.list_motifs)))

        self.total_time_train = end_time - start_0_time


    def ext_prediction_pd(self, test_matrix):
        start_0_time = time.time()
        dt_dftest = pd.DataFrame(test_matrix, columns=self.list_motifs)
        end_time = time.time()

        print(" time creating DataFrame = ", end_time - start_0_time)

        start_time = time.time()
        result = self.ExTrCl.predict(dt_dftest)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")
        print(" Time pred + creat = ", end_time - start_0_time, " s")

        return result

    def ext_prediction_dt(self, test_matrix):
        start_0_time = time.time()
        dt_dftest = dt.Frame(np.array(test_matrix), names=self.list_motifs)
        end_time = time.time()

        print(" time creatind Fram dt = ", end_time - start_0_time)

        start_time = time.time()
        result = self.ExTrCl.predict(dt_dftest)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")
        print(" Time pred + creat = ", end_time - start_0_time, " s")

        return result

    def ext_prediction_numpy(self, test_matrix):
        start_time = time.time()
        result = self.ExTrCl.predict(test_matrix)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")

        return result
    
    def _prediction(self, test_matrix, model):
        start_time = time.time()
        result = model.predict(test_matrix)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")

        return result

    def nlp_prediction_numpy(self, test_matrix):
        start_time = time.time()
        result = self.nlp.predict(test_matrix)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")

        return result


    def rdf_prediction_numpy(self, test_matrix):
        start_time = time.time()
        result = self.rdf.predict(test_matrix)
        end_time = time.time()

        print(" Time pred only = ", end_time - start_time, " s")

        return result

    
    def test_groupe_score_dt(self, test_matrix, list_classes):
        start_0_time = time.time()
        dt_dftest = dt.Frame(np.array(test_matrix), names=self.list_motifs)
        end_time = time.time()

        print(" time creatind Fram dt = ", end_time - start_0_time)

        start_time = time.time()
        result = self.ExTrCl.predict(dt_dftest)
        end_time = time.time()

        score_test = accuracy_score(list_classes, result)

        print(" Time pred only = ", end_time - start_time, " s")
        print(" Time pred + creat = ", end_time - start_0_time, " s")
        print(" score test = ", score_test)

    def compute_classification_report_global_average(self, y_true, y_pred):
        
        print("{ ------------------- classification report : ")

        report = classification_report(y_true, y_pred)

        # Get accuracy, macro avg, and weighted avg from report
        lines = report.split('\n')
        print(" precision    recall  f1-score   support")
        print (lines[-4]) # accuracy
        print (lines[-3]) # macro avg
        print (lines[-2]) # weighted avg

        print(" classification report ------------------- }")




    def test_groupe_score_numpy(self, test_matrix, list_classes):
        start_0_time = time.time()
        start_time = time.time()
        # result = self.ExTrCl.predict(dt_dftest)
        result = self.ExTrCl.predict(test_matrix)
        end_time = time.time()

        score_test = accuracy_score(list_classes, result)

        print(" Time pred only = ", end_time - start_time, " s")
        print(" Time pred + creat = ", end_time - start_0_time, " s")
        print(" score test = ", score_test)

    def ext_test_group_score(self, dir_in_files, file_ext, method='np'):
        
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        end_time_ctm = time.time()

        if method == 'pd':
            result_pred = self.ext_prediction_pd(matrix_test)
        elif method == 'dt':
            result_pred = self.ext_prediction_dt(matrix_test)
        else:
            result_pred = self.ext_prediction_numpy(matrix_test)

        end_time_test = time.time()

        # check the score accuracy
        list_all_classes = self.get_all_classes_ids(dir_in_files, file_ext)
        score_test = accuracy_score(list_all_classes, result_pred)
        result_scores = precision_recall_fscore_support(list_all_classes, result_pred, average='macro')

        # add classification_report only global average score
        self.compute_classification_report_global_average(list_all_classes, result_pred)
        

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.total_time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.precision = result_scores[0]
        self.recall = result_scores[1]
        self.f1 = result_scores[2]

    def nlp_test_group_score(self, dir_in_files, file_ext, method='np'):
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        end_time_ctm = time.time()

        result_pred = self.nlp_prediction_numpy(matrix_test)

        end_time_test = time.time()

        # check the score accuracy
        list_all_classes_ids = self.get_all_classes_ids(dir_in_files, file_ext)
        score_test = accuracy_score(list_all_classes_ids, result_pred)
        result_scores = precision_recall_fscore_support(list_all_classes_ids, result_pred, average='macro')

        # add classification_report only global average score
        self.compute_classification_report_global_average(list_all_classes_ids, result_pred)


        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.total_time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.precision = result_scores[0]
        self.recall = result_scores[1]
        self.f1 = result_scores[2]

        print("-----------------------------------------------------")

        list_all_classes_names = self.get_all_classes_names(dir_in_files, file_ext)
        res_cla_reprot = classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4, output_dict=True)

        per_class_accuracies = self.get_accuracy_for_individual_class(list_all_classes_ids, list_all_classes_names,
                                                                      result_pred)

        #print(classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4))
        #print(res_cla_reprot)
        id_cls_sh = "RF00906_328_seqs_shuffled.fasta"
        id_cls_rd = "RF00906_328_random_seqs.fasta"
        id_cls_di = "RF00906_328_seqs_dinucleotide.fasta"

        if id_cls_sh in per_class_accuracies:
            print("{} = Accuracy: {} | {} ".format(id_cls_sh, per_class_accuracies[id_cls_sh], res_cla_reprot[id_cls_sh]))
        if id_cls_rd in per_class_accuracies:
            print("{} = Accuracy: {} | {} ".format(id_cls_rd, per_class_accuracies[id_cls_rd], res_cla_reprot[id_cls_rd]))
        if id_cls_di in per_class_accuracies:
            print("{} = Accuracy: {} | {} ".format(id_cls_di, per_class_accuracies[id_cls_di], res_cla_reprot[id_cls_di]))

    def _test_group_score(self, dir_in_files, file_ext, model):
            start_time_test = time.time()
            list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

            start_time_ctm = time.time()
            matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
            end_time_ctm = time.time()

            result_pred = self._prediction(matrix_test, model)

            end_time_test = time.time()

            # check the score accuracy
            list_all_classes_ids = self.get_all_classes_ids(dir_in_files, file_ext)
            score_test = accuracy_score(list_all_classes_ids, result_pred)
            result_scores = precision_recall_fscore_support(list_all_classes_ids, result_pred, average='macro')

            # add classification_report only global average score
            self.compute_classification_report_global_average(list_all_classes_ids, result_pred)


            # printing:
            print("list all seqs len: ", len(list_all_seqs))
            print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
            print(" --> time test (all):", end_time_test - start_time_test, " s")
            print(" Accuracy score test = ", score_test)
            print(" Other test score = ", result_scores)

            self.total_time_test = end_time_test - start_time_test
            self.accuracy = score_test
            self.precision = result_scores[0]
            self.recall = result_scores[1]
            self.f1 = result_scores[2]

            print("-----------------------------------------------------")

            list_all_classes_names = self.get_all_classes_names(dir_in_files, file_ext)
            res_cla_reprot = classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4, output_dict=True)

            per_class_accuracies = self.get_accuracy_for_individual_class(list_all_classes_ids, list_all_classes_names,
                                                                        result_pred)

            #print(classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4))
            #print(res_cla_reprot)
            id_cls_sh = "RF00906_328_seqs_shuffled.fasta"
            id_cls_rd = "RF00906_328_random_seqs.fasta"
            id_cls_di = "RF00906_328_seqs_dinucleotide.fasta"

            if id_cls_sh in per_class_accuracies:
                print("{} = Accuracy: {} | {} ".format(id_cls_sh, per_class_accuracies[id_cls_sh], res_cla_reprot[id_cls_sh]))
            if id_cls_rd in per_class_accuracies:
                print("{} = Accuracy: {} | {} ".format(id_cls_rd, per_class_accuracies[id_cls_rd], res_cla_reprot[id_cls_rd]))
            if id_cls_di in per_class_accuracies:
                print("{} = Accuracy: {} | {} ".format(id_cls_di, per_class_accuracies[id_cls_di], res_cla_reprot[id_cls_di]))


    def rdf_test_group_score(self, dir_in_files, file_ext):
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        end_time_ctm = time.time()

        result_pred = self.rdf_prediction_numpy(matrix_test)

        end_time_test = time.time()

        # check the score accuracy
        list_all_classes_ids = self.get_all_classes_ids(dir_in_files, file_ext)
        score_test = accuracy_score(list_all_classes_ids, result_pred)
        result_scores = precision_recall_fscore_support(list_all_classes_ids, result_pred, average='macro')

        # add classification_report only global average score
        self.compute_classification_report_global_average(list_all_classes_ids, result_pred)

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.total_time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.precision = result_scores[0]
        self.recall = result_scores[1]
        self.f1 = result_scores[2]

        print("-----------------------------------------------------")

        list_all_classes_names = self.get_all_classes_names(dir_in_files, file_ext)
        res_cla_reprot = classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4, output_dict=True)

        per_class_accuracies = self.get_accuracy_for_individual_class(list_all_classes_ids, list_all_classes_names,
                                                                      result_pred)

        #print(classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4))
        #print(res_cla_reprot)
        id_cls_sh = "RF00906_328_seqs_shuffled.fasta"
        id_cls_rd = "RF00906_328_random_seqs.fasta"
        id_cls_di = "RF00906_328_seqs_dinucleotide.fasta"

        if id_cls_sh in per_class_accuracies:
            print("{} = Accuracy: {} | {} ".format(id_cls_sh, per_class_accuracies[id_cls_sh], res_cla_reprot[id_cls_sh]))
        if id_cls_rd in per_class_accuracies:
            print("{} = Accuracy: {} | {} ".format(id_cls_rd, per_class_accuracies[id_cls_rd], res_cla_reprot[id_cls_rd]))
        if id_cls_di in per_class_accuracies:
            print("{} = Accuracy: {} | {} ".format(id_cls_di, per_class_accuracies[id_cls_di], res_cla_reprot[id_cls_di]))


    
    def get_accuracy_for_individual_class(self, list_class_ids, list_class_names, result_pred):

        # Get the confusion matrix
        cm = confusion_matrix(list_class_ids, result_pred)

        # We will store the results in a dictionary for easy access later
        per_class_accuracies = {}

        # Calculate the accuracy for each one of our classes
        for idx, cls in enumerate(list_class_names):
            # True negatives are all the samples that are not our current GT class (not the current row)
            # and were not predicted as the current class (not the current column)
            true_negatives = np.sum(np.delete(np.delete(cm, idx, axis=0), idx, axis=1))

            # True positives are all the samples of our current GT class that were predicted as such
            true_positives = cm[idx, idx]

            # The accuracy for the current class is ratio between correct predictions to all predictions
            per_class_accuracies[cls] = (true_positives + true_negatives) / np.sum(cm)

        return per_class_accuracies

    def train_motifs_oneChars(self, dir_in_files, file_ext):
        start_0_time = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        self.list_motifs = ["A", "C", "U", "G"]

        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        list_all_classes = self.get_all_classes_ids(dir_in_files, file_ext)

        dt_frame = dt.Frame(np.array(matrix_test), names=self.list_motifs)

        self.ExTrCl.fit(dt_frame, list_all_classes)

        end_time = time.time()

        # statistics
        score_train = self.ExTrCl.score(dt_frame, list_all_classes)
        pred_train = self.ExTrCl.predict(dt_frame)
        train_pred_score = accuracy_score(list_all_classes, pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", score_train, " | Train pred score = ", train_pred_score)
        print(" Len list motifs : ", len(list(dt_frame)))
        print("list motifs len : {}".format(len(self.list_motifs)))

    def train_motifs_oneChars_MLP(self, dir_in_files, file_ext):
        start_0_time = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        self.list_motifs = ["A", "C", "U", "G"]

        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        list_all_classes = self.get_all_classes_ids(dir_in_files, file_ext)

        dt_frame = dt.Frame(np.array(matrix_test), names=self.list_motifs)

        self.nlp.fit(dt_frame, list_all_classes)

        end_time = time.time()

        # statistics
        score_train = self.nlp.score(dt_frame, list_all_classes)
        pred_train = self.nlp.predict(dt_frame)
        train_pred_score = accuracy_score(list_all_classes, pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", score_train, " | Train pred score = ", train_pred_score)
        print(" Len list motifs : ", len(list(dt_frame)))
        print("list motifs len : {}".format(len(self.list_motifs)))

        self.total_time_train = end_time - start_0_time
        self.score_train = score_train
        self.train_pred_score = train_pred_score

    def train_voting(self, csv_file_path):
        start_0_time = time.time()
        data_arn = dt.fread(csv_file_path)
        end_time = time.time()
        print(" time Read_csv file : ", end_time - start_0_time, " s")

        # start_time = time.time()
        data_classe = np.ravel(data_arn[:, "familyId"])
        # end_time = time.time()
        # print(" time get class array = ", end_time - start_time)

        del data_arn[:, "familyId"]

        start_time_train = time.time()

        self.voted_model.fit(data_arn, data_classe)

        end_time = time.time()

        print(" train only time : ", end_time - start_time_train, " s")

        self.list_motifs = data_arn.names  # save list of use motifs in classification

        # statistics
        self.score_train = self.voted_model.score(data_arn, np.ravel(data_classe))
        pred_train = self.voted_model.predict(data_arn)
        self.train_pred_score = accuracy_score(np.ravel(data_classe), pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", self.score_train, " | Train pred score = ", self.train_pred_score)
        print(" Len list motifs : ", len(list(data_arn)))
        print("list motifs len : {}".format(len(self.list_motifs)))

        self.total_time_train = end_time - start_0_time

    def test_voting(self, dir_in_files, file_ext):
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        end_time_ctm = time.time()

        start_time = time.time()
        result_pred = self.voted_model.predict(matrix_test)
        end_time = time.time()
        print(" Time pred only = ", end_time - start_time, " s")

        end_time_test = time.time()

        # check the score accuracy
        list_all_classes = self.get_all_classes_ids(dir_in_files, file_ext)
        score_test = accuracy_score(list_all_classes, result_pred)
        result_scores = precision_recall_fscore_support(list_all_classes, result_pred, average='macro')

        # add classification_report only global average score
        self.compute_classification_report_global_average(list_all_classes, result_pred)

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.total_time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.precision = result_scores[0]
        self.recall = result_scores[1]
        self.f1 = result_scores[2]

    def _test(self, model):
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(self.dir_test_files, self.test_file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        end_time_ctm = time.time()
        
        start_time_pred_only = time.time()
        self.test_result_pred = model.predict(matrix_test)
        end_time_pred_only = time.time()
        

        end_time_test = time.time()
        
        # Check the score accuracy
        self.list_all_classes_ids = self.get_all_classes_ids(self.dir_test_files, self.test_file_ext)
        
        # Compute and store basic metrics
        self.accuracy = accuracy_score(self.list_all_classes_ids, self.test_result_pred)
        self.precision, self.recall, self.f1, _ = precision_recall_fscore_support(self.list_all_classes_ids, self.test_result_pred, average='macro')

        # Store results
        self.number_test_seqs = len(list_all_seqs)
        self.time_test_matrix_construction = end_time_ctm - start_time_ctm
        self.time_test_prediction_only = end_time_pred_only - start_time_pred_only
        self.total_time_test = end_time_test - start_time_test # total time test = time read seqs + time matrix construction + time prediction


    def print_detailed_results(self):
        list_all_classes_names = self.get_all_classes_names(self.dir_test_files, self.test_file_ext)
        
        # Generate the classification report as a string for printing
        report_str = classification_report(self.list_all_classes_ids, self.test_result_pred, 
                                           target_names=list_all_classes_names, digits=4)
        
        # Generate the classification report as a dictionary for further processing
        classification_report = classification_report(self.list_all_classes_ids, self.test_result_pred, 
                                                           target_names=list_all_classes_names, 
                                                           digits=4, output_dict=True)

        print("\nDetailed Classification Report:")
        print(report_str)

        print("\nConfusion Matrix:")
        cm = confusion_matrix(self.list_all_classes_ids, self.test_result_pred)
        print(cm)

        per_class_accuracies = self.get_accuracy_for_individual_class(self.list_all_classes_ids, 
                                                                      list_all_classes_names, 
                                                                      self.test_result_pred)
        print("\nPer-class Accuracies:")
        for cls, acc in per_class_accuracies.items():
            print(f"{cls}: {acc:.4f}")

        # Special classes reporting
        special_classes = ["RF00906_328_seqs_shuffled.fasta", 
                           "RF00906_328_random_seqs.fasta", 
                           "RF00906_328_seqs_dinucleotide.fasta"]
        print("\nSpecial Classes Report:")
        for cls in special_classes:
            if cls in per_class_accuracies and cls in self.classification_report:
                print(f"{cls} = Accuracy: {per_class_accuracies[cls]:.4f} | {self.classification_report[cls]}")

        


    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    @staticmethod
    def get_all_fasta_files_in_the_dir(dir_in, file_ext):

        mylist = []

        for file in os.listdir(dir_in):
            if file.endswith(file_ext):
                mylist.append(os.path.join(dir_in, file))

        mylist.sort()
        return mylist

    @staticmethod
    def get_list_seqs(infile):
        list_seqs = []
        for record in SeqIO.parse(infile, 'fasta'):
            list_seqs.append(str(record.seq))

        return list_seqs

    @staticmethod
    def get_all_seqs(dir_in, file_ext):

        list_files = Model.get_all_fasta_files_in_the_dir(dir_in, file_ext)
        list_seqs = []

        for file in list_files:
            list_seqs += Model.get_list_seqs(file)

        return list_seqs

    @staticmethod
    def get_all_classes_names(dir_in, file_ext):

        list_files = Model.get_all_fasta_files_in_the_dir(dir_in, file_ext)

        list_classes_names = []

        for file in list_files:
            list_classes_names.append(os.path.basename(file))

        return list_classes_names

    @staticmethod
    def get_all_classes_ids(dir_in, file_ext):

        list_files = Model.get_all_fasta_files_in_the_dir(dir_in, file_ext)

        list_seqs = None
        list_classes_ids = []
        idCls = 0

        for file in list_files:
            list_seqs = Model.get_list_seqs(file)
            for i in range(0, len(list_seqs)):
                list_classes_ids.append(idCls)
            idCls += 1

        return list_classes_ids

    @staticmethod
    def get_all_seqs_and_thier_classes(dir_in, file_ext):

        list_files = Model.get_all_fasta_files_in_the_dir(dir_in, file_ext)

        list_seqs = None
        list_classes = []
        idCls = 0

        for file in list_files:
            print(" file : ", file, " +++++++ Classe = ", idCls)
            list_seqs = Model.get_list_seqs(file)
            for i in range(0, len(list_seqs)):
                list_classes.append(idCls)
            idCls += 1

        return list_classes

    def write_results_to_csv_file(self, file_out):
        
        separator = ','
        line_result = self.full_test_name + separator + str(self.total_time_train) + separator + str(
            self.score_train) + separator + str(self.train_pred_score) + separator + str(
            self.total_time_test) + separator + str(self.accuracy) + separator + str(self.precision) + separator + str(
            self.recall) + separator + str(self.f1) + separator + "\n"

        lock = portalocker.RedisLock('csv_write_lock')

        with lock:
            with open(file_out, "a+") as my_file:
                my_file.write(line_result)

    def print_results(self, is_detailed_report=True):
        print("EXPERIMENT_RESULTS_START")
        print(f"Full_test_name: {self.full_test_name}")
        print(f"Training_Time: {self.total_time_train}")
        print(f"Training_Score: {self.score_train}")
        print(f"Training_Prediction_Score: {self.train_pred_score}")
        print(f"Classification_Time: {self.total_time_test}")
        print(f"Accuracy: {self.accuracy}")
        print(f"Precision: {self.precision}")
        print(f"Recall: {self.recall}")
        print(f"F1_Score: {self.f1}")
        print("EXPERIMENT_RESULTS_END")

        if is_detailed_report:
            self.print_detailed_results()
