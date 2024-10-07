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
    parser.add_argument("-m", "--mlm", choices=['EXT', 'MLP', 'VOT', 'RDF'], help="Choose machine learning model")
    parser.add_argument("-t", "--train-csv", help="Path to training CSV matrix")
    parser.add_argument("-d", "--test-dir", help="Path to test files directory")
    parser.add_argument("-e", "--file-ext", default=".fasta.txt", help="File extension for test files")
    parser.add_argument("-j", "--n-job", type=int, default=1, help="Number of jobs to run in parallel. -1 means using all processors.")
    parser.add_argument("-o", "--output-csv", help="Path to output CSV file for results")
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
    print("  -t, --train-csv <file>: Path to the training CSV matrix file")
    print("  -d, --test-dir <dir>  : Path to the directory containing test files")
    print("  -e, --file-ext <ext>  : File extension for test files (default: .fasta.txt)")
    print("  -j, --n-job <int>     : Number of jobs to run in parallel (-1 uses all processors, default: 1)")
    print("  -o, --output-csv <file>: Path to the output CSV file for results")
    print("\nExample:")
    print("  python Main.py -m EXT -t path/to/train.csv -d path/to/test/dir -e .fasta -j 4")


def print_input_args(args):
    print(f"Machine learning model: {args.mlm}")
    print(f"Training CSV file: {args.train_csv}")
    print(f"Test directory: {args.test_dir}")
    print(f"File extension: {args.file_ext}")
    print(f"Number of jobs: {args.n_job}")
    print(f"Output CSV file: {args.output_csv if args.output_csv else 'Not specified'}")


def main():
    args = parse_arguments()
    
    if len(sys.argv) == 1:
        print_args_info()
        sys.exit(0)
    
    if not all([args.mlm, args.train_csv, args.test_dir]):
        print("Error: Missing required arguments.")
        print_args_info()
        sys.exit(1)

    print_input_args(args)

    clm = Model(args.n_job, args.train_csv)

    print("\nTraining model...")
    clm.train(args.mlm)

    print("\nTesting model...")
    clm.test(args.test_dir, args.file_ext)

    if args.output_csv:
        print(f"\nWriting results to {args.output_csv}...")
        clm.write_results_to_csv_file(args.output_csv)

    clm.print_results(is_detailed_report=False)

    return 0

if __name__ == '__main__':
    main()
