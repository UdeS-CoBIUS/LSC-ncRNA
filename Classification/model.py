import os
import time
import fcntl

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

from xgboost import XGBClassifier


class Model:

    def __init__(self, n_job=1):
        local_n_estimators = 700
        self.ExTrCl = ExtraTreesClassifier(n_estimators=local_n_estimators, n_jobs=n_job)
        self.nlp = MLPClassifier(max_iter=800)
        #self.rdf = RandomForestClassifier(n_jobs=n_job)
        self.rdf = RandomForestClassifier(n_estimators=local_n_estimators, n_jobs=n_job)
        self.xgb = XGBClassifier(n_estimators=local_n_estimators, n_jobs=n_job)
        self.list_motifs = None
        self.voted_model = VotingClassifier(estimators=[
            ('ext', ExtraTreesClassifier(n_estimators=local_n_estimators, n_jobs=n_job)),
            ('nlp', MLPClassifier(max_iter=800)),
            ('rdf', ExtraTreesClassifier(n_estimators=local_n_estimators, n_jobs=n_job))
        ], voting='soft', weights=[1, 1, 1], n_jobs=n_job)

        self.time_train = 0
        self.score_train = 0
        self.train_pred_score = 0
        self.time_test = 0
        self.accuracy = 0
        self.Precision = 0
        self.Recall = 0
        self.Fbs = 0

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

    def train_with_dt(self, csv_file_path):
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

        self.time_train = end_time - start_0_time

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

        self.time_train = end_time - start_0_time


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

        self.time_train = end_time - start_0_time


    def xgb_train_with_dt(self, csv_file_path):
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
        self.xgb.fit(data_arn, data_classe)
        end_time = time.time()

        print(" train only time : ", end_time - start_time_train, " s")

        self.list_motifs = data_arn.names  # save list of use motifs in classification

        # statistics
        self.score_train = self.xgb.score(data_arn, np.ravel(data_classe))
        pred_train = self.xgb.predict(data_arn)
        self.train_pred_score = accuracy_score(np.ravel(data_classe), pred_train)

        print("all Time Train = ", end_time - start_0_time, " s")
        print(" score_train = ", self.score_train, " | Train pred score = ", self.train_pred_score)
        print(" Len list motifs : ", len(list(data_arn)))
        print("list motifs len : {}".format(len(self.list_motifs)))

        self.time_train = end_time - start_0_time



    def prediction_pd(self, test_matrix):
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

    def prediction_dt(self, test_matrix):
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

    def prediction_numpy(self, test_matrix):
        start_time = time.time()
        result = self.ExTrCl.predict(test_matrix)
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

    def xgb_prediction_numpy(self, test_matrix):

        start_time = time.time()
        result = self.xgb.predict(test_matrix)
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
        
        report = classification_report(y_true, y_pred)

        # Get accuracy, macro avg, and weighted avg from report
        lines = report.split('\n')

        print ("accuracy = ",lines[-3])
        print ("macro avg = ", lines[-2])
        print ("weighted avg = ", lines[-1])

        # precision at idx 1
        # recall    at idx 2
        # f1-score  at idx 3
        # support   at idx 4
        
        print (" computing  : ")
        accuracy = float(lines[-3].split()[3])
        macro_avg = [float(num) for num in lines[-2].split()[-4:]]
        weighted_avg = [float(num) for num in lines[-1].split()[-4:]]

        # Print to console
        print(" precision    recall  f1-score   support")
        print(f"Accuracy: {accuracy}")
        print(f"Macro Avg: {macro_avg}")
        print(f"Weighted Avg: {weighted_avg}")




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

    def test_group_score(self, dir_in_files, file_ext, method='np'):
        
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)
        end_time_ctm = time.time()

        if method == 'pd':
            result_pred = self.prediction_pd(matrix_test)
        elif method == 'dt':
            result_pred = self.prediction_dt(matrix_test)
        else:
            result_pred = self.prediction_numpy(matrix_test)

        end_time_test = time.time()

        # check the score accuracy
        list_all_classes = self.get_all_classes_ids(dir_in_files, file_ext)
        score_test = accuracy_score(list_all_classes, result_pred)
        result_scores = precision_recall_fscore_support(list_all_classes, result_pred, average='macro')

        # add classification_report only global average score
        compute_classification_report_global_average(list_all_classes, result_pred)
        

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.Precision = result_scores[0]
        self.Recall = result_scores[1]
        self.Fbs = result_scores[2]

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
        compute_classification_report_global_average(list_all_classes, result_pred)


        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.Precision = result_scores[0]
        self.Recall = result_scores[1]
        self.Fbs = result_scores[2]

        print("-----------------------------------------------------")

        list_all_classes_names = self.get_all_classes_names(dir_in_files, file_ext)
        res_cla_reprot = classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4, output_dict=True)

        per_class_accuracies = self.get_accuracy_for_individeul_class(list_all_classes_ids, list_all_classes_names,
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
        compute_classification_report_global_average(list_all_classes, result_pred)

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.Precision = result_scores[0]
        self.Recall = result_scores[1]
        self.Fbs = result_scores[2]

        print("-----------------------------------------------------")

        list_all_classes_names = self.get_all_classes_names(dir_in_files, file_ext)
        res_cla_reprot = classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4, output_dict=True)

        per_class_accuracies = self.get_accuracy_for_individeul_class(list_all_classes_ids, list_all_classes_names,
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


    def xgb_test_group_score(self, dir_in_files, file_ext):
        start_time_test = time.time()
        list_all_seqs = self.get_all_seqs(dir_in_files, file_ext)

        start_time_ctm = time.time()
        matrix_test = ctm.get_matrix_nbOcrrs_listStr_AhoCorasick(self.list_motifs, list_all_seqs)

        # added to slove xgb pbm
        np_matrix_test = np.array(matrix_test)
        dt_df_test = dt.Frame(np_matrix_test, names=self.list_motifs)

        end_time_ctm = time.time()

        #result_pred = self.xgb_prediction_numpy(matrix_test)
        result_pred = self.xgb_prediction_numpy(dt_df_test)

        end_time_test = time.time()

        # check the score accuracy
        list_all_classes_ids = self.get_all_classes_ids(dir_in_files, file_ext)
        score_test = accuracy_score(list_all_classes_ids, result_pred)
        result_scores = precision_recall_fscore_support(list_all_classes_ids, result_pred, average='macro')

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.Precision = result_scores[0]
        self.Recall = result_scores[1]
        self.Fbs = result_scores[2]

        print("-----------------------------------------------------")

        list_all_classes_names = self.get_all_classes_names(dir_in_files, file_ext)
        res_cla_reprot = classification_report(list_all_classes_ids, result_pred, target_names=list_all_classes_names, digits=4, output_dict=True)

        per_class_accuracies = self.get_accuracy_for_individeul_class(list_all_classes_ids, list_all_classes_names,
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

    def get_accuracy_for_individeul_class(self, list_class_ids, list_class_names, result_pred):

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

        self.time_train = end_time - start_0_time
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

        self.time_train = end_time - start_0_time

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
        compute_classification_report_global_average(list_all_classes, result_pred)

        # printing:
        print("list all seqs len: ", len(list_all_seqs))
        print(" time matrix construction By AC :", end_time_ctm - start_time_ctm, " s")
        print(" --> time test (all):", end_time_test - start_time_test, " s")
        print(" Accuracy score test = ", score_test)
        print(" Other test score = ", result_scores)

        self.time_test = end_time_test - start_time_test
        self.accuracy = score_test
        self.Precision = result_scores[0]
        self.Recall = result_scores[1]
        self.Fbs = result_scores[2]

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

    def write_results_to_csv_file(self, file_out, family_name):
        my_file = open(file_out, "a+")
        separator = ','
        line_result = family_name + separator + str(self.time_train) + separator + str(
            self.score_train) + separator + str(self.train_pred_score) + separator + str(
            self.time_test) + separator + str(self.accuracy) + separator + str(self.Precision) + separator + str(
            self.Recall) + separator + str(self.Fbs) + separator + "\n"

        fcntl.flock(my_file, fcntl.LOCK_EX) # because i use multi processing to write to the same file.
        my_file.write(line_result)
        fcntl.flock(my_file, fcntl.LOCK_UN)
        my_file.close()
