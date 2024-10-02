import os
import sys

# Generalized Suffix-Tree example.
# a = ["xxxabcxxx", "adsaabc", "ytysabcrew", "qqqabcqw", "aaabc"]
# st = STree.STree(a)
# print(st.lcs()) # "abc"

from model import Model

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="RNA Classification")
    parser.add_argument("-h", "--help", action="store_true", help="Show this help message and exit")
    parser.add_argument("-m", "--mlm", choices=['EXT', 'MLP', 'VOT', 'RDF', 'XGB'], help="Choose machine learning model")
    parser.add_argument("-t", "--train-csv", help="Path to training CSV matrix")
    parser.add_argument("-d", "--test-dir", help="Path to test files directory")
    parser.add_argument("-e", "--file-ext", default=".fasta.txt", help="File extension for test files")
    parser.add_argument("-j", "--n-job", type=int, default=-1, help="Number of jobs to run in parallel. -1 means using all processors.")
    return parser.parse_args()

def print_args_info():
    print("RNA Classification - Usage Information")
    print("======================================")
    print("This program classifies RNA sequences using various machine learning models.")
    print("\nUsage: python Main.py -m <model> -t <train_csv> -d <test_dir> [-e <file_ext>] [-j <n_job>]")
    print("\nArguments:")
    print("  -h, --help            : Show this help message and exit")
    print("  -m, --mlm <model>     : Machine learning model to use. Choose from:")
    print("                          EXT (Extremely Randomized Trees)")
    print("                          MLP (Multi-Layer Perceptron)")
    print("                          VOT (Voting Classifier)")
    print("                          RDF (Random Forest)")
    print("                          XGB (XGBoost)")
    print("  -t, --train-csv <file>: Path to the training CSV matrix file")
    print("  -d, --test-dir <dir>  : Path to the directory containing test files")
    print("  -e, --file-ext <ext>  : File extension for test files (default: .fasta.txt)")
    print("  -j, --n-job <int>     : Number of jobs to run in parallel (-1 uses all processors, default: -1)")
    print("\nExample:")
    print("  python Main.py -m EXT -t path/to/train.csv -d path/to/test/dir -e .fasta -j 4")


def print_input_args(args):
    print(f"Machine learning model: {args.mlm}")
    print(f"Training CSV file: {args.train_csv}")
    print(f"Test directory: {args.test_dir}")
    print(f"File extension: {args.file_ext}")
    print(f"Number of jobs: {args.n_job}")

def main():
    args = parse_arguments()
    
    if args.help or len(sys.argv) == 1:
        print_args_info()
        sys.exit(0)
    
    if not all([args.mlm, args.train_csv, args.test_dir]):
        print("Error: Missing required arguments.")
        print_args_info()
        sys.exit(1)

    choose_mlm = args.mlm
    dir_in_train_csv_matrix = args.train_csv
    dir_in_test_files = args.test_dir
    file_ext = args.file_ext
    n_job = args.n_job

    print_input_args(args)
    # exit for debug
    exit()
    
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
    
    test_name = f"{args.mlm}_{os.path.basename(args.train_csv)}"
    clm.write_results_to_csv_file(file_out, test_name)

    return 0

if __name__ == '__main__':
    main()
