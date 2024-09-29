# this file is used as a pipeline to run all experiments of the paper,
#  it have 3 main steps:
# 1. prepare the dataset (unzip, split, etc.)
# 2. run experiments: 
#   a- extract and select motifs using C++ code --> generate a matrix saved in csv file. 
#   b- run the experiments with the matrix and the models
#   c- plot the results
# 3. for each step, save the parameters like time and space complexity, memory and disk usage, etc. 

import os
import zipfile
import subprocess
import re
import csv
import re
import csv
import matplotlib.pyplot as plt
import pandas as pd

from pathlib import Path

# define a global variable to use debug datasets
is_debug_datasets_global_var: bool = True
debug_datasets_size_global_var: list[int] = [5, 10, 15, 20, 25, 30]
dir_main_path_debug_datasets_global_var: str = "datasets/data/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test"

def find_project_root(project_name: str = "LSC-ncRNA") -> Path:
    """
    Find the project root directory by looking for a specific directory name.
    :param project_name: The name of the project directory
    :return: The absolute path to the project root directory as a Path object
    """
    current_dir: Path = Path(__file__).resolve().parent

    print(f"current_dir: {current_dir}")

    while current_dir != current_dir.parent:  # Stop at the root directory
        if current_dir.name == project_name:
            return current_dir
        current_dir = current_dir.parent

    raise FileNotFoundError(f"Project root directory '{project_name}' not found.")


# datasets: ------------------------------
# Datasets are pre-prepared and available in the `datasets/data/` directory.
# there are 3 datasets:
# datasets/data/Clans_ncRNA_from_Rfam_14.8.zip
# datasets/data/deep_ncrna_datasets.zip
# datasets/data/Rfam_14.1_dataset.zip


def prepare_dataset():
    # Get the project root directory
    project_root: Path = find_project_root()

    # Define the directory containing the datasets
    data_dir: Path = project_root / "datasets/data"

    # List of dataset zip files
    datasets: list[str] = [
        "Clans_ncRNA_from_Rfam_14.8.zip",
        "deep_ncrna_datasets.zip",
        "Rfam_14.1_dataset.zip"
    ]

    for dataset in datasets:
        zip_path: Path = os.path.join(data_dir, dataset)
        extract_dir: Path = os.path.join(data_dir, dataset[:-4])  # Remove .zip extension

        print(f"Unzipping {dataset}...")

        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
            print(f"Successfully unzipped {dataset} to {extract_dir}")
        except zipfile.BadZipFile:
            print(f"Error: {dataset} is not a valid zip file.")
        except FileNotFoundError:
            print(f"Error: {dataset} not found in {data_dir}")
        except Exception as e:
            print(f"An error occurred while unzipping {dataset}: {str(e)}")


def prepare_dataset_from_scratch():
    # ... existing code ...
    print("Preparing dataset from scratch...")
    # Implement dataset preparation steps here
    # e.g., unzipping, splitting, etc.


# delete sub motifs: ------------------------------
# We use the medium dataset including 600 ncRNA families for which we considered subsamples of increasing sizes \{100,200,300,400,500,600\}.
# The initial set of motifs is obtained by naively computing all common linear motifs (sub-strings) of length between $minl = 2$ and $maxl = 20$ between any two sequences of a family.
# there is no filtering parameters is used.
# beta: 0 , is defined as percentage_same_family, a cm is accepted if :   ((double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo) >= percentage_same_family )
# alpha: -1, no filtering on the number of occurences of the cm
# gamma nbOccrs_allowed: 1, no filtering on the number of occurences of the cm, at least we have 1 (the default)

# is_delete_subMotifs: True, to delete sub motifs
#   true: with filtering
#   false: without filtering

# there is no filtering parameters is used.
# beta: 0 , is defined as percentage_same_family, a cm is accepted if :   ((double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo) >= percentage_same_family )
# alpha: -1, no filtering on the number of occurences of the cm
# gamma nbOccrs_allowed: 1, no filtering on the number of occurences of the cm, at least we have 1 (the default)

# is_delete_subMotifs: True, to delete sub motifs
#   true: with filtering
#   false: without filtering

# We generate 2 figures:
# Evolution of the data size (Left figure) and processing time (Right figure) for the generation of datasetby filtering out exact submotifs\textbf{(F)}, and without any  filtering \textbf{(NF)}.
# for that:
# test on \{100,200,300,400,500,600\}
# $minl = 2$ and $maxl = 20$
# -d : <integer> (0: false, 1 or other: true), is delete sub-motifs

# for that:
# test on \{100,200,300,400,500,600\}
# $minl = 2$ and $maxl = 20$
# -d : <integer> (0: false, 1 or other: true), is delete sub-motifs


def get_rfam14_sample_train_test_dir_path(size: int, is_debug_datasets: bool = False) -> Path:
    """
    used:  datasets/data/Rfam_14.1_dataset/Rfam14.1_Sample: 50, 100, 150, 200, 250, 300, 350, 400, 500, 600
    debug: datasets/data/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test: 5, 10, 15, 20, 25, 30
    Get the directory path for the given dataset size.
    :param size: The size of the dataset
    :param is_debug_datasets: Whether to use debug datasets
    :return: The directory path as a Path object
    """
    # Get the project root directory
    project_root: Path = find_project_root()

    if is_debug_datasets:
        # Check for double subfolder (when unzipped by Python)
        double_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test/Train"
        # Check for single subfolder (when unzipped by other apps)
        single_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test/Train"
    else:
        # Check for double subfolder (when unzipped by Python)
        double_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test/Train"
        # Check for single subfolder (when unzipped by other apps)
        single_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test/Train"

    if double_subfolder_path.exists():
        return double_subfolder_path
    elif single_subfolder_path.exists():
        return single_subfolder_path
    else:
        raise FileNotFoundError(
            f"Directory not found for size {size}. Checked paths:\n{double_subfolder_path}\n{single_subfolder_path}")


def deletion_sub_motifs(is_debug_datasets: bool = False) -> None:
    # Get the project root directory
    project_root: Path = find_project_root()


    dataset_sizes: list[int] = [100, 200, 300, 400, 500, 600]
    
    if is_debug_datasets:
        dataset_sizes = debug_datasets_size_global_var
    
    test_name: str = f"test_dnd"  # test delete no delete
    min_length: int = 2
    max_length: int = 20
    beta: int = 0  # don't use beta (0 is the smallest value)
    alpha: int = -1  # don't use alpha
    gamma: int = 1  # 1 is smallest value

    results: list[dict[str, int | float | str]] = []

    # Step 2: Iterate through dataset sizes and submotif deletion options
    for size in dataset_sizes:
        # Step 3: Set up input file paths
        dir_path = get_rfam14_sample_train_test_dir_path(size, is_debug_datasets)

        for is_delete_submotifs in [0, 1]:
            
            try:
                # Step 4: Run the C++ program using the new function
                result = run_cpp_motif_extraction_and_selection(
                    dir_path=str(dir_path),  # Convert Path to string
                    test_name=test_name,
                    dataset_size=size,
                    min_length=min_length,
                    max_length=max_length,
                    is_delete_submotifs=is_delete_submotifs,
                    beta=beta,
                    alpha=alpha,
                    gamma=gamma
                )

                if result is None:
                    print(f"Error processing dataset size {size} with submotif deletion {'enabled' if is_delete_submotifs else 'disabled'}")
                    continue

                # Step 8: Store the results
                results.append({
                    "dataset_size": size,
                    "is_delete_submotifs": bool(is_delete_submotifs),
                    "execution_time": result['execution_time'],
                    "file_size_gb": result['file_size_gb']
                })

                print(
                    f"Processed dataset size {size} with submotif deletion {'enabled' if is_delete_submotifs else 'disabled'}")

                # Step 9: Save results to a CSV file
                results_dir: Path = Path("results")
                os.makedirs(results_dir, exist_ok=True)  # Create the directory if it doesn't exist
                csv_path: Path = os.path.join(results_dir, "deletion_sub_motifs_results.csv")

                with open(csv_path, "w", newline="") as csvfile:
                    fieldnames: list[str] = ["dataset_size", "is_delete_submotifs", "execution_time", "file_size_gb"]
                    writer: csv.DictWriter = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for row in results:
                        writer.writerow(row)
                print(f"Results saved to {csv_path}")

                # Step 10: Clean up the generated CSV file, since it will have big size
                ## os.remove(result['output_csv_file'])
                ## print(f"Cleaned up {result['output_csv_file']}")

            except Exception as e:
                print(
                    f"Error processing dataset size {size} with submotif deletion {'enabled' if is_delete_submotifs else 'disabled'}: {str(e)}")

                # generate the plot of the results in 2 figures:
    # we compare 2 methods:
    # A- filtering out exact submotifs \textbf{(F)}, 
    # B- Without any  filtering \textbf{(NF)}
    # 1. generate Evolution of the data size for F vs NF:  
    # 2. generate evolution of processing time in minute  for F vs NF

    # Call the function to generate plots
    generate_result_sub_motifs_F_vs_NF_plots(csv_path)


def generate_result_sub_motifs_F_vs_NF_plots(csv_file_path: Path) -> None:
    # Read the CSV file into a DataFrame
    df: pd.DataFrame = pd.read_csv(csv_file_path)

    # Plot 1: Evolution of the data size for F vs NF
    plt.figure(figsize=(10, 6))
    for is_delete, group in df.groupby('is_delete_submotifs'):
        label: str = 'F' if is_delete else 'NF'
        plt.plot(group['dataset_size'], group['file_size_gb'], marker='o', label=label)

    plt.xlabel('Dataset Size')
    plt.ylabel('File Size (GB)')
    plt.title('Evolution of Data Size')
    plt.legend()
    plt.grid(True)
    plt.savefig('results/evolution_of_data_size.png')
    plt.close()
    print("Plot saved as results/evolution_of_data_size.png")

    # Plot 2: Evolution of processing time for F vs NF
    plt.figure(figsize=(10, 6))
    for is_delete, group in df.groupby('is_delete_submotifs'):
        label: str = 'F' if is_delete else 'NF'
        plt.plot(group['dataset_size'], group['execution_time'] / 60, marker='o', label=label)

    plt.xlabel('Dataset Size')
    plt.ylabel('Processing Time (minutes)')
    plt.title('Evolution of Processing Time')
    plt.legend()
    plt.grid(True)
    plt.savefig('results/evolution_of_processing_time.png')
    plt.close()
    print("Plot saved as results/evolution_of_processing_time.png")


# -------------------------------------------
# -------------------------------------------
# Experiments on: number of generated motifs (fixed-length and combined-length) vs classification accuracy (EXT and MLP): ------------------------------
# datasets used: 50, 150, 250, 350
#
# figures fixed-length:
# - Growth_NB_motifs.PNG: Number of generated motifs for fixed length motifs minl=maxl from 2 to 20. for dataset : 50, 150, 250, 350.
# - k-mer_EXT.PNG: Accuracy EXT for fixed length motifs for minl=maxl from 1 to 20. for dataset : 50, 150, 250, 350.
# - k-mer_MLP.PNG: Accuracy MLP for fixed length motifs for minl=maxl from 1 to 20. for dataset : 50, 150, 250, 350.

# Figure combined-length:
# - Growth_NB_motifs_For_350_min_2_maxTo_20.PNG: Number of generated motifs for combined lenght motifs minl=2, maxl from 2 to 20. for dataset : 350.

# - EXT_NLP_nbF350_min2_max_variable.PNG: Accuracy EXT and MLP for minl motif = 2. maxl motifs from 2 to 20 ([2-20]).
# - EXT_NLP_nbF350_min2_max_variable_time.PNG: Time taken by whole program (matrix generation + classification):  for minl motif = 2. maxl motifs from 2 to 10 ([2-10]).

# \begin{figure*}[h]
#    %\centering
#    \includegraphics[width=.5\textwidth]{images/Growth_NB_motifs.PNG}\hfill
#    \includegraphics[width=.5\textwidth]{images/Growth_NB_motifs_For_350_min_2_maxTo_20.PNG}

#    \includegraphics[width=.5\textwidth]{images/k-mer_EXT.PNG}\hfill
#    \includegraphics[width=.5\textwidth]{images/EXT_NLP_nbF350_min2_max_variable.PNG}\hfill

#    %\includegraphics[width=.5\textwidth]{images/k-mer_NLP.PNG}\hfill
#    %\includegraphics[width=.5\textwidth]{}\hfill%images/EXT_NLP_nbF350_min2_max_variable_time.PNG}\hfill

#    \caption{Relation between the number of generated motifs and : (Top-Left figure) the size of fixed-length motifs; (Top-Right figure) the size of combined-length motifs. Classification accuracy for: (Bottom-Left figure) different sizes of fixed-length motifs using EXT; (Bottom-Right figure) different sizes of combined-length motifs using EXT and MLP.}
#    \label{fig:fixed_combined_ext_nlp}
#    \vspace{-.5cm}
#   \end{figure*}

# generate the common motifs csv file:
# we use the following parameters:
# fixed-length motifs: minl=maxl from 1 to 20
# combined-length motifs: minl=2, maxl from 2 to 20
# there is no filtering parameters is used.
# -d : <integer> (0: false), there is no delete sub-motifs
# -b: beta: 0 , is defined as percentage_same_family, a cm is accepted if :   ((double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo) >= percentage_same_family )
# -a : alpha: -1, no filtering on the number of occurences of the cm
# -g : gamma nbOccrs_allowed: 1, no filtering on the number of occurences of the cm, at least we have 1 (the default)
# => get the number of generated motifs
# => get the processing time for generating the csv file.
#
# classification:
# use EXT and MLP for the classification of the motifs
# use the csv file generated in the previous step
# => get the accuracy for the classification
# => get the processing time for the classification
#
# save the results:
# number of generated motifs
# total processing time: generating the csv file + classification
# accuracy for the classification

# plot the results: as explained before for 6 figures.


import os
import subprocess
import re
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.tree import ExtraTreeClassifier
from sklearn.metrics import accuracy_score


def run_motif_length_experiments(is_debug_datasets: bool = False) -> None:
    # Step 1: Define parameters
    dataset_sizes: list[int] = [50, 150, 250, 350]
    
    if is_debug_datasets:
        dataset_sizes = debug_datasets_size_global_var

    fixed_lengths: range = range(1, 21)  # 1 to 20
    combined_lengths: list[tuple[int, int]] = [(2, max_len) for max_len in range(2, 21)]  # (2, 2) to (2, 20)

    results: dict[str, dict[int, dict[int | tuple[int, int], dict[str, int | float | str]]]] = {
        'fixed_length': {size: {} for size in dataset_sizes},
        'combined_length': {dataset_sizes[-1]: {}} # take only last size, so 350 or 30 in debug datasets
    }

    # Step 2: Run experiments for fixed-length motifs
    for size in dataset_sizes:
        for cm_length in fixed_lengths:
            results['fixed_length'][size][cm_length] = run_single_experiment(size, cm_length, cm_length)

    # Step 3: Run experiments for combined-length motifs (only for size 350)
    for min_len, max_len in combined_lengths:
        results['combined_length'][350][(min_len, max_len)] = run_single_experiment(350, min_len, max_len)

    # Step 4: Generate plots
    generate_plots(results)


def run_single_experiment(dataset_size, min_length, max_length):
    # Get the project root directory
    project_root = find_project_root()

    # Set up paths and parameters
    dir_path = get_rfam14_sample_train_test_dir_path(dataset_size)
    test_name = "cm_len"

    # Run the C++ program using run_cpp_motif_extraction_and_selection
    result = run_cpp_motif_extraction_and_selection(
        dir_path=str(dir_path),
        test_name=test_name,
        dataset_size=dataset_size,
        min_length=min_length,
        max_length=max_length,
        is_delete_submotifs=0,  # no deletion of sub motifs
        beta=0,  # don't use beta
        alpha=-1,  # don't use alpha
        gamma=1  # gamma = 1, no filtering on the number of occurrences of the cm, at least we have 1 (the default)
    )

    if result is None:
        print(f"Error processing dataset size {dataset_size} with min_length {min_length} and max_length {max_length}")
        return None

    # Perform classification
    ext_accuracy, mlp_accuracy, classification_time = perform_classification(result['output_csv_file'])

    # Clean up the generated CSV file
    if os.path.exists(result['output_csv_file']):
        os.remove(result['output_csv_file'])
        print(f"Cleaned up {result['output_csv_file']}")
    else:
        print(f"Warning: Output file {result['output_csv_file']} not found")

    return {
        'num_motifs': result['num_motifs'],
        'execution_time': result['execution_time'],
        'ext_accuracy': ext_accuracy,
        'mlp_accuracy': mlp_accuracy,
        'classification_time': classification_time,
        'total_time': result['execution_time'] + classification_time
    }


def perform_classification(csv_file):
    # Load data from CSV
    data = pd.read_csv(csv_file)
    X = data.iloc[:, :-1]  # All columns except the last one
    y = data.iloc[:, -1]  # Last column as target

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Scale the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # ExtraTreeClassifier
    ext = ExtraTreeClassifier(random_state=42)
    ext.fit(X_train_scaled, y_train)
    ext_accuracy = accuracy_score(y_test, ext.predict(X_test_scaled))

    # MLPClassifier
    mlp = MLPClassifier(random_state=42, max_iter=1000)
    mlp.fit(X_train_scaled, y_train)
    mlp_accuracy = accuracy_score(y_test, mlp.predict(X_test_scaled))

    # Measure classification time (you might want to use a more precise timing method in practice)
    classification_time = 0  # Placeholder, implement actual timing if needed

    return ext_accuracy, mlp_accuracy, classification_time


def generate_plots(results):
    # Plot 1: Growth_NB_motifs.PNG
    plt.figure(figsize=(10, 6))
    for size in results['fixed_length']:
        lengths = list(results['fixed_length'][size].keys())
        num_motifs = [results['fixed_length'][size][l]['num_motifs'] for l in lengths]
        plt.plot(lengths, num_motifs, marker='o', label=f'Dataset size {size}')
    plt.xlabel('Motif Length')
    plt.ylabel('Number of Motifs')
    plt.title('Number of Generated Motifs (Fixed Length)')
    plt.legend()
    plt.savefig('results/Growth_NB_motifs.png')
    plt.close()

    # Plot 2: k-mer_EXT.PNG
    plt.figure(figsize=(10, 6))
    for size in results['fixed_length']:
        lengths = list(results['fixed_length'][size].keys())
        accuracies = [results['fixed_length'][size][l]['ext_accuracy'] for l in lengths]
        plt.plot(lengths, accuracies, marker='o', label=f'Dataset size {size}')
    plt.xlabel('Motif Length')
    plt.ylabel('EXT Accuracy')
    plt.title('EXT Accuracy for Fixed Length Motifs')
    plt.legend()
    plt.savefig('results/k-mer_EXT.png')
    plt.close()

    # Plot 3: k-mer_MLP.PNG
    plt.figure(figsize=(10, 6))
    for size in results['fixed_length']:
        lengths = list(results['fixed_length'][size].keys())
        accuracies = [results['fixed_length'][size][l]['mlp_accuracy'] for l in lengths]
        plt.plot(lengths, accuracies, marker='o', label=f'Dataset size {size}')
    plt.xlabel('Motif Length')
    plt.ylabel('MLP Accuracy')
    plt.title('MLP Accuracy for Fixed Length Motifs')
    plt.legend()
    plt.savefig('results/k-mer_MLP.png')
    plt.close()

    # Plot 4: Growth_NB_motifs_For_350_min_2_maxTo_20.PNG
    plt.figure(figsize=(10, 6))
    max_lengths = [max_len for _, max_len in results['combined_length'][350].keys()]
    num_motifs = [results['combined_length'][350][(2, max_len)]['num_motifs'] for max_len in max_lengths]
    plt.plot(max_lengths, num_motifs, marker='o')
    plt.xlabel('Max Motif Length')
    plt.ylabel('Number of Motifs')
    plt.title('Number of Generated Motifs (Combined Length, Dataset Size 350)')
    plt.savefig('results/Growth_NB_motifs_For_350_min_2_maxTo_20.png')
    plt.close()

    # Plot 5: EXT_NLP_nbF350_min2_max_variable.PNG
    plt.figure(figsize=(10, 6))
    max_lengths = [max_len for _, max_len in results['combined_length'][350].keys()]
    ext_accuracies = [results['combined_length'][350][(2, max_len)]['ext_accuracy'] for max_len in max_lengths]
    mlp_accuracies = [results['combined_length'][350][(2, max_len)]['mlp_accuracy'] for max_len in max_lengths]
    plt.plot(max_lengths, ext_accuracies, marker='o', label='EXT')
    plt.plot(max_lengths, mlp_accuracies, marker='o', label='MLP')
    plt.xlabel('Max Motif Length')
    plt.ylabel('Accuracy')
    plt.title('EXT and MLP Accuracy (Combined Length, Dataset Size 350)')
    plt.legend()
    plt.savefig('results/EXT_NLP_nbF350_min2_max_variable.png')
    plt.close()

    # Plot 6: EXT_NLP_nbF350_min2_max_variable_time.PNG
    plt.figure(figsize=(10, 6))
    max_lengths = [max_len for _, max_len in results['combined_length'][350].keys() if max_len <= 10]
    total_times = [results['combined_length'][350][(2, max_len)]['total_time'] for max_len in max_lengths]
    plt.plot(max_lengths, total_times, marker='o')
    plt.xlabel('Max Motif Length')
    plt.ylabel('Total Time (seconds)')
    plt.title('Total Processing Time (Combined Length, Dataset Size 350)')
    plt.savefig('results/EXT_NLP_nbF350_min2_max_variable_time.png')
    plt.close()


# -------------------------------------------
# -------------------------------------------

# compile code_MotifsExtractionSelection: ------------------------------
# Struct to store the parsed command-line arguments
# struct Args {

#    string dir_name;  // -in
#    string test_name = "test"; // -tn
#    int nb_families = 10; // -nf
#    int min_nb_seqs_allowed = 4; // -mins
#    int max_nb_seqs_allowed = 1000; // -maxs
#    int min_length_motif = 2; // -minl
#    int max_length_motif = 10; // -maxl
#    int is_delete_subMotifs = 0; // -d (false, true)
#    int Beta = 40; // -b
#    int Alpha = -1; // -a   // -1, mean don't use the Alpha paramter ==> whatever the variance we accepte it.
#    int nbOccrs_allowed = 1; // -g (gamma) , lower bound, 1 is defult
# };
#
# then the output csv file name is constructed like:
#     string output_csv_file = "del_No_nbF_";
#
#    if(args.is_delete_subMotifs) output_csv_file = "del_Yes_nbF_";
#
#    output_csv_file += args.test_name;
#   output_csv_file += "_min_" + util::to_string(args.min_length_motif);
#    output_csv_file += "_max_" + util::to_string(args.max_length_motif);
#    output_csv_file += "_beta_" + util::to_string(args.Beta);
#    output_csv_file += "_alpha_" + util::to_string(args.Alpha);
#    output_csv_file += "_nbOccrs_" + util::to_string(args.nbOccrs_allowed);
# 
# The output csv file name is as follows: del_[No/Yes:-d]_nbF_[test_name:-tn]_min_[-minl]_max_[-maxl]_beta_[-b]_alpha_[-a]_nbOccrs_[-g].csv.

from typing import Dict, Union


def generate_csv_output_filename(args: Dict[str, Union[str, int, bool]]) -> str:
    """
    Generate the output CSV filename based on the given arguments.
    
    :param args: A dictionary containing the argument values
    :return: The generated output CSV filename
    """
    output_csv_file = (
        f"{args['-tn']}"  # test_name
        f"_nbF_{args['-nf']}"  # nb_families (if provided)
        f"_is_del_{'yes' if args['-d'] else 'no'}"  # is_delete_subMotifs
        f"_min_{args['-minl']}"  # min_length_motif
        f"_max_{args['-maxl']}"  # max_length_motif
        f"_beta_{args['-b']}"  # Beta
        f"_alpha_{args['-a']}"  # Alpha
        f"_nbOccrs_{args['-g']}"  # nbOccrs_allowed (gamma)
    )

    return output_csv_file + ".csv"


def compile_code_MotifsExtractionSelection():
    # Get the project root directory
    project_root = find_project_root()

    cpp_dir = project_root / "LSC-ncRNA-our_method/MotifsExtractionSelection"

    executable_name = "MotifsExtractionSelection"
    command = [
        "g++",
        "-std=c++14", # c++17
        "-o", os.path.join(cpp_dir, executable_name),
        os.path.join(cpp_dir, "Main.cpp"),
        os.path.join(cpp_dir, "SuffixTree_QuadraticTime.cpp"),
        os.path.join(cpp_dir, "Old_CommonMotifs.cpp"),
        os.path.join(cpp_dir, "FastaFilesReader.cpp"),
        os.path.join(cpp_dir, "CommonMotifs.cpp")
    ]

    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Successfully compiled {executable_name}")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed. Error: {e.stderr.decode()}")
    except Exception as e:
        print(f"An error occurred during compilation: {str(e)}")

# -------------------------------------------
# -------------------------------------------
# run C++ code for motifs extraction and selection

def run_cpp_motif_extraction_and_selection(dir_path, test_name, dataset_size, min_length, max_length, is_delete_submotifs=0, beta=0, alpha=-1, gamma=1):
    """
    Run the motif extraction and selection process and return the results.

    :param dir_path: Input directory path
    :param test_name: Test name
    :param dataset_size: number of families (Size of the dataset)
    :param min_length: Minimum common motif length
    :param max_length: Maximum common motif length
    :param is_delete_submotifs: Whether to delete submotifs (0 or 1)
    :param beta: Beta parameter
    :param alpha: Alpha parameter
    :param gamma: Gamma parameter
    :return: A dictionary containing the results
    """

    # Get the project root directory
    project_root = find_project_root()

    cpp_dir = project_root / "LSC-ncRNA-our_method/MotifsExtractionSelection"

    executable = os.path.join(cpp_dir, "MotifsExtractionSelection")

    # check if the executable exists
    if not os.path.exists(executable):
        print(f"Error executable does not exists: {executable} not found in {cpp_dir}")
        return

    command = [
        executable,
        "-in", dir_path,
        "-tn", test_name,
        "-nf", str(dataset_size),
        "-minl", str(min_length),
        "-maxl", str(max_length),
        "-d", str(is_delete_submotifs),
        "-b", str(beta),
        "-a", str(alpha),
        "-g", str(gamma)
    ]

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    print(f"Standard error output: {stderr}")
    print(f"Standard output: {stdout}")

    if process.returncode != 0:
        print(f"Error running C++ program. Return code: {process.returncode}")
        return None

    # Extract execution time and number of motifs
    time_match = re.search(r"Time taken by whole program is : (\d+\.\d+) sec", stdout)
    execution_time = float(time_match.group(1)) if time_match else None

    motifs_match = re.search(r"Total nb_motifs = (\d+)", stdout)
    num_motifs = int(motifs_match.group(1)) if motifs_match else None

    # Generate CSV output filename
    args = {
        "-tn": test_name,
        "-nf": dataset_size,
        "-minl": min_length,
        "-maxl": max_length,
        "-d": is_delete_submotifs,
        "-b": beta,
        "-a": alpha,
        "-g": gamma
    }
    output_csv_file = generate_csv_output_filename(args)

    # Get the size of the generated CSV file
    file_size = os.path.getsize(output_csv_file) / (1024 * 1024 * 1024)  # Convert to GB

    return {
        'num_motifs': num_motifs,
        'execution_time': execution_time,
        'file_size_gb': file_size,
        'output_csv_file': output_csv_file
    }




def main():

    is_debug_datasets: bool = is_debug_datasets_global_var

    # unzip already pre-prepared datasets
    prepare_dataset()

    # compile the c++ code for motifs extraction and selection
    compile_code_MotifsExtractionSelection()
    
    # run the sub motifs deletion experiments
    deletion_sub_motifs(is_debug_datasets)
    


if __name__ == "__main__":
    main()
