import os
import sys

# Generalized Suffix-Tree example.
# a = ["xxxabcxxx", "adsaabc", "ytysabcrew", "qqqabcqw", "aaabc"]
# st = STree.STree(a)
# print(st.lcs()) # "abc"

from model import Model


def main():

    choose_mlm = sys.argv[1]  # 'EXT', 'NLP', 'VOT'
    dir_in_train_csv_matrix = sys.argv[2]
    file_ext = ".fasta.txt" #todo: remove the dependcy on the extention, because some times I change the exetesion, and I get error only for that, and it take time to descover that error come from here.
    #file_ext = ".fa"
    dir_in_test_files = "/data/chei2402/ibra/test_infernal/nbF_all_nbSeqs_min_3/Test"
    dir_in_test_files = sys.argv[3]

    # n_job = -1 # -1 use multi-processing which is availabale only in EXT and rdf not in NLP
    # n_job = int(sys.argv[2])
    n_job = -1

    clm = Model(n_job)

    # clm.test_create_dt_df()
    # return

    # print(" train_test_from_one_CSV_file : ---------------------- ")
    # clm.train_test_from_one_CSV_file(dir_in_train_test_csv_matrix,test_size)
    # clm.train_test_from_one_CSV_file_dt(dir_in_train_test_csv_matrix,test_size)

    print("\n Train and then Test from real seqs : -----------------------")
    print(" Train using matrix : -----------------------")
    # clm.train(dir_in_train_csv_matrix)

    if choose_mlm == 'EXT':
        print("train for EXT ")
        clm.train_with_dt(dir_in_train_csv_matrix)
    elif choose_mlm == 'NLP':
        print("train for nlp ")
        clm.nlp_train_with_dt(dir_in_train_csv_matrix)
        # clm.train_motifs_oneChars_MLP(dir_in_train_csv_matrix,file_ext)
    elif choose_mlm == 'RDF':
        print("train for RDF ")
        clm.rdf_train_with_dt(dir_in_train_csv_matrix)
    elif choose_mlm == 'XGB':
        print("train for XGB ")
        clm.xgb_train_with_dt(dir_in_train_csv_matrix)
    else:
        print("train for voting model")
        clm.train_voting(dir_in_train_csv_matrix)

    # clm.train_motifs_oneChars(dir_in_train_files, file_ext)

    # print("\n Pred Train by seqs in : ----------------------- ")
    # clm.test_group_score(dir_in_train_files, file_ext)
    # return 0

    print("\n Pred Test by seqs in : ----------------------- ")

    if choose_mlm == 'EXT':
        print("test for EXT")
        clm.test_group_score(dir_in_test_files, file_ext)
    elif choose_mlm == 'NLP':
        print("test for nlp")
        clm.nlp_test_group_score(dir_in_test_files, file_ext)
    elif choose_mlm == 'RDF':
        print("test for RDF")
        clm.rdf_test_group_score(dir_in_test_files, file_ext)
    elif choose_mlm == "XGB":
        print("test for Xgboost")
        clm.xgb_test_group_score(dir_in_test_files, file_ext)
    else:
        print("test for voting model")
        clm.test_voting(dir_in_test_files, file_ext)

    file_out = "Secondary_noSecondary_trainTest_results.csv"
    test_name = "{}_{}".format(choose_mlm,os.path.basename(dir_in_train_csv_matrix))
    clm.write_results_to_csv_file(file_out, test_name)

    return 0

if __name__ == '__main__':
    main()
