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
from typing import Optional, Dict, Any

import time 

from dic_values_exp_example_debug import dict_len_motifs_exp_example_debug

# define a global variable to use debug datasets
is_debug_datasets_global_var: bool = True
debug_datasets_size_global_var: list[int] = [5, 10, 15, 20, 25, 30]
#debug_datasets_size_global_var: list[int] = [30]
dir_main_path_debug_datasets_global_var: str = "datasets/data/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test"

def find_project_root(project_name: str = "LSC-ncRNA") -> Path:
    """
    Find the project root directory by looking for a specific directory name.
    :param project_name: The name of the project directory
    :return: The absolute path to the project root directory as a Path object
    """
    current_dir: Path = Path(__file__).resolve().parent

    # print(f"current_dir: {current_dir}") # just for debug...

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


def get_dir_path_rfam14_sample(size: int, data_type: str, is_debug_datasets: bool = False) -> Path:
    """
    used:  datasets/data/Rfam_14.1_dataset/Rfam14.1_Sample: 50, 100, 150, 200, 250, 300, 350, 400, 500, 600
    debug: datasets/data/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test: 5, 10, 15, 20, 25, 30
    Get the directory path for the given dataset size.
    
     Args:
        size: The size of the dataset
        data_type: Either 'Train' or 'Test'
        is_debug_datasets: Whether to use debug datasets
    
    Returns:
        The directory path as a Path object
    
    Raises:
        ValueError: If data_type is not 'Train' or 'Test'
        FileNotFoundError: If the directory doesn't exist
    """
    if data_type not in ['Train', 'Test']:
        raise ValueError("data_type must be either 'Train' or 'Test'")
    
    # Get the project root directory
    project_root: Path = find_project_root()

    if is_debug_datasets:
        # Check for double subfolder (when unzipped by Python)
        double_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test" / data_type
        # Check for single subfolder (when unzipped by other apps)
        single_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/debug_small_Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test" / data_type
    else:
        # Check for double subfolder (when unzipped by Python)
        double_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test" / data_type
        # Check for single subfolder (when unzipped by other apps)
        single_subfolder_path: Path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test" / data_type

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
        dir_path = get_dir_path_rfam14_sample(size=size, data_type='Train', is_debug_datasets=is_debug_datasets)

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
                os.remove(result['output_csv_file'])
                print(f"Cleaned up {result['output_csv_file']}")

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


def run_motif_length_experiments(is_debug_datasets: bool = False) -> None:
    # Step 1: Define parameters
    dataset_sizes: list[int] = [50, 150, 250, 350]
    
    if is_debug_datasets:
        dataset_sizes = debug_datasets_size_global_var

    beta: int= 0 # percentage of the cm in the family (called alpha in the paper)
    alpha: int= -1 # variance of nb cm between seqs in family , -1 we don't care
    gamma: int= 1 # min nb of cm to consider 
    is_delete_submotifs: int= 0 # delete of note sub motif, default: 0 so no.
    

    #fixed_lengths: range = range(1, 21)  # for 1 we should use an intern python method to generate the dataset, the method exist already.
    fixed_lengths: range = range(2, 21)  # 2 to 20, fro now from 2. to be changed later, we treat size min max = 1 seperatly.
    # fixed_lengths: list[tuple[int, int]] = [(f-1, f) for f in range(2, 21)]  # (1, 2), (2,3) to (19, 20). we do this since in c++ program it implemented (min_len, max_len) like min_len<max_len begin min, stop before max. so 1,2 mean only one I have to check...
    combined_lengths: list[tuple[int, int]] = [(2, max_len) for max_len in range(3, 21)]  # (2, 3) to (2, 20)

    # note: when we have min_len = i, max_len = i+1 ex:(2,3) this mean we have fixed motif len of len 2. (this is in my c++ code work like this, maybe I should take time to change it...)

    results: dict[str, dict[int, dict[int | tuple[int, int], dict[str, int | float | str]]]] = {
        'fixed_length': {size: {} for size in dataset_sizes},
        'combined_length': {dataset_sizes[-1]: {}} # take only last size, so 350 or 30 in debug datasets
    }


    # Step 2: Run experiments for fixed-length motifs
    for size in dataset_sizes:
        print(f"Running experiments for dataset size {size}")
        for cm_length in fixed_lengths:
            results['fixed_length'][size][cm_length] = run_single_experiment(
                dataset_size=size, 
                min_length=cm_length-1, 
                max_length=cm_length,
                beta_percentage=beta,
                alpha_nb_occ_variation=alpha,
                gamma_nb_occrs=gamma,
                is_delete_submotifs=is_delete_submotifs,
                test_name='cm_len', 
                is_debug_datasets=is_debug_datasets)

    # Step 3: Run experiments for combined-length motifs (only for size 350)
    for min_len, max_len in combined_lengths:
        size = 350
        if is_debug_datasets:
            size = 30
        results['combined_length'][size][(min_len, max_len)] = run_single_experiment(
            dataset_size=size, 
            min_length=min_len, 
            max_length=max_len,
            beta_percentage=beta,
            alpha_nb_occ_variation=alpha,
            gamma_nb_occrs=gamma,
            is_delete_submotifs=is_delete_submotifs, 
            test_name='cm_len',
            is_debug_datasets=is_debug_datasets)

    # Step 4: Generate plots
    plot_all_len_motifs_exp(results)


def run_single_experiment(dataset_size: int = 500, 
                          min_length: int = 2, 
                          max_length: int = 10,
                          beta_percentage: int = 0, # default: no filtering
                          alpha_nb_occ_variation: int = -1, # default: no filtering
                          gamma_nb_occrs: int = 1, # default: no filtering
                          is_delete_submotifs: int = 0, # default: false , no deletion 
                          test_name: str = 'test_no_f',
                          is_debug_datasets: bool = False) -> Optional[Dict[str, Any]]:
    """
    Runs a single experiment including motif extraction and classification.

    :param dataset_size: Size of the dataset (number of families)
    :param min_length: Minimum motif length
    :param max_length: Maximum motif length
    :return: Dictionary containing experiment results or None if failed
    """
    project_root = find_project_root()

    # Set up paths and parameters
    train_dir_path = get_dir_path_rfam14_sample(size=dataset_size, data_type='Train', is_debug_datasets=is_debug_datasets)
    test_dir_path = get_dir_path_rfam14_sample(size=dataset_size, data_type='Test', is_debug_datasets=is_debug_datasets)
    

    #print(' inside run single experiment ...')
    #print('Train dir:', train_dir_path)
    #print('Test dir:', test_dir_path)
    #print(' run cpp program...')

    # Run the C++ program using run_cpp_motif_extraction_and_selection
    result = run_cpp_motif_extraction_and_selection(
        dir_path=str(train_dir_path),
        test_name=test_name,
        dataset_size=dataset_size,
        min_length=min_length,
        max_length=max_length,
        is_delete_submotifs= is_delete_submotifs,  # no deletion of sub motifs
        beta= beta_percentage,  # don't use beta
        alpha=alpha_nb_occ_variation,  # don't use alpha
        gamma=gamma_nb_occrs  # gamma = 1, no filtering on the number of occurrences of the cm, at least we have 1 (the default)
    )

    #print(' end cpp, results:')
    #print(result)
    #print(' go to classification setp...')
    #time.sleep(5)

    if result is None:
        print(f"Error processing dataset size {dataset_size} with min_length {min_length} and max_length {max_length}")
        return None

    # Initialize classification results
    classification_results = {}

    # Define models to run
    models = ['EXT', 'MLP']

    for model in models:
        print(f"\nRunning classification with model: {model}")

        # Define output CSV for classification results
        classification_output_csv = os.path.join(
            "results",
            f"classification_results_{model}_dataset_{dataset_size}_min{min_length}_max{max_length}.csv"
        )

        #print('classification_output_csv: ', classification_output_csv)
        #print('befor run classification ...')
        #time.sleep(3)
        # Run classification experiment
        classification_output = run_classification_experiment(
            model=model,
            train_csv=result['output_csv_file'],
            test_dir=str(test_dir_path),  # Adjust based on your test data directory structure
            file_ext=".fasta.txt",
            n_job=1,  # Adjust as needed
            output_csv=classification_output_csv
        )

        #print('end classification')
        #print('classification_output:')
        #print(classification_output)
        #time.sleep(3)

        if "error" in classification_output:
            print(f"Classification error for model {model}: {classification_output['error']}")
            classification_results[model] = {
                "accuracy": None,
                "precision": None,
                "recall": None,
                "f1_score": None,
                "training_time_sec": None,
                "testing_time_sec": None
            }
            #print(' ibra error ibra error')
            continue

        # Store the extracted metrics
        classification_results[model] = {
            "accuracy": classification_output.get("accuracy"),
            "precision": classification_output.get("precision"),
            "recall": classification_output.get("recall"),
            "f1_score": classification_output.get("f1_score"),
            "training_time_sec": classification_output.get("training_time_sec"),
            "testing_time_sec": classification_output.get("testing_time_sec")
        }

        # Clean up the CSV file that store the classification result. 
        os.remove(classification_output_csv)
        print(f"Cleaned up {classification_output_csv}")

        print(f"Model: {model} | Accuracy: {classification_results[model]['accuracy']} | "
              f"Training Time: {classification_results[model]['training_time_sec']} sec | "
              f"Testing Time: {classification_results[model]['testing_time_sec']} sec")

    # Combine C++ results with classification results
    experiment_results = {
        "dataset_size": dataset_size,
        "min_length": min_length,
        "max_length": max_length,
        "num_motifs": result.get('num_motifs'),
        "cpp_execution_time_sec": result.get('execution_time'),
        "file_size_gb": result.get('file_size_gb'),
    }

    # Merge classification_results into experiment_results
    experiment_results.update(classification_results)

    # Clean up the generated CSV file (train csv file), since we no longer need it & it will have big size
    os.remove(result['output_csv_file'])
    print(f"Cleaned up {result['output_csv_file']}")
    

    print(f"Completed experiment for dataset size {dataset_size}, min_length {min_length}, max_length {max_length}")

    return experiment_results 
    

def run_single_experiment_models_choices(dataset_size: int = 500, 
                          min_length: int = 2, 
                          max_length: int = 10,
                          beta_percentage: int = 0, # default: no filtering
                          alpha_nb_occ_variation: int = -1, # default: no filtering
                          gamma_nb_occrs: int = 1, # default: no filtering
                          is_delete_submotifs: int = 0, # default: false , no deletion 
                          test_name: str = 'test_no_f',
                          is_debug_datasets: bool = False) -> Optional[Dict[str, Any]]:
    """
    Runs a single experiment including motif extraction and classification.

    :param dataset_size: Size of the dataset (number of families)
    :param min_length: Minimum motif length
    :param max_length: Maximum motif length
    :return: Dictionary containing experiment results or None if failed
    """
    project_root = find_project_root()

    # Set up paths and parameters
    train_dir_path = get_dir_path_rfam14_sample(size=dataset_size, data_type='Train', is_debug_datasets=is_debug_datasets)
    
    # Run the C++ program using run_cpp_motif_extraction_and_selection
    result = run_cpp_motif_extraction_and_selection(
        dir_path=str(train_dir_path),
        test_name=test_name,
        dataset_size=dataset_size,
        min_length=min_length,
        max_length=max_length,
        is_delete_submotifs=is_delete_submotifs,
        beta=beta_percentage,
        alpha=alpha_nb_occ_variation,
        gamma=gamma_nb_occrs
    )

    if result is None:
        print(f"Error processing dataset size {dataset_size} with min_length {min_length} and max_length {max_length}")
        return None

    print(f"\nRunning classification with models choices...")

    # Define output CSV for classification results
    classification_output_csv = os.path.join(
        "results",
        f"classification_results_models_choices_dataset_{dataset_size}_min{min_length}_max{max_length}.csv"
    )

    # Run classification experiment
    classification_output = run_classification_experiment_models_choices(
        train_csv=result['output_csv_file'],
        output_csv=classification_output_csv
    )

    if "error" in classification_output:
        print(f"Classification error: {classification_output['error']}")
        return None

    # Store classification results
    classification_results = {}
    for i, model_name in enumerate(classification_output["models"]):
        classification_results[model_name] = {
            "processing_time": classification_output["processing_times"][i],
            "accuracy": classification_output["accuracies"][i]
        }
        print(f"Model: {model_name} | Accuracy: {classification_output['accuracies'][i]} | "
              f"Processing Time: {classification_output['processing_times'][i]} sec")

    # Clean up the classification results CSV
    if os.path.exists(classification_output_csv):
        os.remove(classification_output_csv)
        print(f"Cleaned up {classification_output_csv}")

    # Combine C++ results with classification results
    experiment_results = {
        "dataset_size": dataset_size,
        "min_length": min_length,
        "max_length": max_length,
        "num_motifs": result.get('num_motifs'),
        "cpp_execution_time_sec": result.get('execution_time'),
        "file_size_gb": result.get('file_size_gb'),
        "classification_results": classification_results
    }

    # Clean up the generated CSV file
    if os.path.exists(result['output_csv_file']):
        os.remove(result['output_csv_file'])
        print(f"Cleaned up {result['output_csv_file']}")

    print(f"Completed experiment for dataset size {dataset_size}, min_length {min_length}, max_length {max_length}")

    return experiment_results


# -------------------------------------------
# -------------------------------------------


def plot_growth_nb_motifs_fixed_len(results, save_as_file=True, filename='Growth_NB_motifs.png'):
    """
    Plot the growth in the number of motifs for fixed-length motifs for various dataset sizes.

    Args:
    results (Dict[str, Any]): A dictionary containing experiment results.
    save_as_file (bool): If True, save the plot as a file. If False, display the plot.
    filename (str): The filename to save the plot if save_as_file is True.
    """
    plt.figure(figsize=(10, 6))
    
    # Extract dataset sizes from the results dictionary
    dataset_sizes = results['fixed_length'].keys()
    
    for dataset in dataset_sizes:
        lengths = []
        num_motifs = []
        
        # Extract motif lengths and corresponding number of motifs for each dataset
        for length, data in results['fixed_length'].get(dataset, {}).items():
            lengths.append(length)
            num_motifs.append(data.get('num_motifs', 0))  # Replace missing values with 0
        
        # Sort lengths and number of motifs by motif length
        sorted_data = sorted(zip(lengths, num_motifs))
        sorted_lengths, sorted_num_motifs = zip(*sorted_data)
        
        # Plot the growth of number of motifs for the current dataset size
        plt.plot(sorted_lengths, sorted_num_motifs, marker='o', label=f'Dataset size {dataset}')
    
    # Customize the plot
    plt.xlabel('Motif Length')
    plt.ylabel('Number of Motifs')
    plt.title('Number of Generated Motifs (Fixed Length)')
    # Set x-axis ticks dynamically based on the sorted motif lengths
    plt.xticks(sorted_lengths)

    plt.grid(True)
    plt.legend()
    
    # Save or display the plot
    if save_as_file:
        plt.savefig(filename)
        plt.close()  # Close the plot to avoid display when saving
    else:
        plt.show()  # Display the plot if not saving


def plot_growth_nb_motifs_combined_len(results, save_as_file=True, filename='Growth_NB_motifs_combined.png'):
    """
    Plot the growth in the number of motifs for combined-length motifs for various dataset sizes.

    Args:
    results (Dict[str, Any]): A dictionary containing experiment results.
    save_as_file (bool): If True, save the plot as a file. If False, display the plot.
    filename (str): The filename to save the plot if save_as_file is True.
    """
    plt.figure(figsize=(10, 6))
    
    # Extract dataset sizes from the results dictionary for 'combined_length'
    combined_length_results = results.get('combined_length', {})
    
    for dataset_size, lengths_data in combined_length_results.items():
        lengths = []
        num_motifs = []
        
        # Extract motif (min_len, max_len) and corresponding number of motifs for each dataset
        for (min_len, max_len), data in lengths_data.items():
            lengths.append((min_len, max_len))  # Store as a tuple for sorting
            num_motifs.append(data.get('num_motifs', 0))  # Replace missing values with 0
        
        # Sort the lengths and number of motifs by max_len, since min_len is always the same
        sorted_data = sorted(zip(lengths, num_motifs), key=lambda x: x[0][1])  # Sort by max_len
        sorted_lengths, sorted_num_motifs = zip(*sorted_data)
        
        # Create formatted x-axis labels (min_len, max_len) after sorting
        sorted_labels = [f'({min_len},{max_len})' for min_len, max_len in sorted_lengths]
        
        # Plot the growth of number of motifs for the current dataset size
        plt.plot(sorted_labels, sorted_num_motifs, marker='o', label=f'Dataset size {dataset_size}')
    
    # Customize the plot
    plt.xlabel('Motif Length (min_len, max_len)')
    plt.ylabel('Number of Motifs')
    plt.title('Number of Generated Motifs (Combined Length)')
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.grid(True)
    plt.legend()
    
    # Save or display the plot
    if save_as_file:
        plt.savefig(filename)
        plt.close()  # Close the plot to avoid display when saving
    else:
        plt.show()  # Display the plot if not saving


def plot_accuracy_fixed_len(results, save_as_file=True):
    """
    Plot the accuracy of models for fixed-length motifs for various dataset sizes.

    Args:
    results (Dict[str, Any]): A dictionary containing experiment results.
    save_as_file (bool): If True, save the plots as files. If False, display the plots.
    """
    # List of models to plot
    models = ['EXT', 'MLP']  # Example models, adjust as needed

    for model in models:
        plt.figure(figsize=(10, 6))
        
        # Loop over dataset sizes
        for dataset, dataset_data in results['fixed_length'].items():
            lengths = []
            accuracies = []

            # Extract motif lengths and corresponding accuracies for the model
            for length, data in dataset_data.items():
                model_data = data.get(model)
                
                if model_data is not None:  # Ensure model data exists
                    accuracy = model_data.get('accuracy', 0)  # Default to 0 if 'accuracy' is missing
                    lengths.append(length)
                    accuracies.append(accuracy)

            # Sort lengths and accuracies by motif length
            if lengths and accuracies:  # Ensure we have data to plot
                sorted_lengths, sorted_accuracies = zip(*sorted(zip(lengths, accuracies)))
                plt.plot(sorted_lengths, sorted_accuracies, marker='o', label=f'Dataset size {dataset}')

        # Customize the plot
        plt.xlabel('Motif Length')
        plt.ylabel('Accuracy')
        plt.title(f'Accuracy for Model {model} (Fixed Length)')
        if lengths:  # Use sorted_lengths for dynamic x-axis if available
            plt.xticks(sorted_lengths)
        
        plt.ylim(0, 1)  # Set y-axis limits for accuracy

        plt.grid(True)
        plt.legend()

        # Save or display the plot
        filename = f'accuracy_fixed_len_{model}.png'
        if save_as_file:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()



def plot_accuracy_combined_len(results, save_as_file=True, filename='Accuracy_combined_len_EXT_MLP.png'):
    """
    Plot the accuracy of EXT and MLP models for combined-length motifs for various dataset sizes.

    Args:
    results (Dict[str, Any]): A dictionary containing experiment results.
    save_as_file (bool): If True, save the plot as a file. If False, display the plot.
    filename (str): The filename to save the plot if save_as_file is True.
    """
    plt.figure(figsize=(12, 6))
    
    # Extract dataset sizes from the results dictionary for 'combined_length'
    combined_length_results = results.get('combined_length', {})
    
    models = ['EXT', 'MLP']
    markers = ['o', 's']  # Different markers for EXT and MLP
    colors = ['blue', 'red']  # Different colors for EXT and MLP
    
    for dataset_size, lengths_data in combined_length_results.items():
        lengths = []
        ext_accuracies = []
        mlp_accuracies = []
        
        # Extract motif (min_len, max_len) and corresponding accuracies for each dataset
        for (min_len, max_len), data in lengths_data.items():
            lengths.append((min_len, max_len))  # Store as a tuple for sorting
            ext_accuracies.append(data.get('EXT', {}).get('accuracy', 0) if isinstance(data.get('EXT'), dict) else 0)
            mlp_accuracies.append(data.get('MLP', {}).get('accuracy', 0) if isinstance(data.get('MLP'), dict) else 0)
        
        # Sort the lengths and accuracies by max_len, since min_len is always the same
        sorted_data = sorted(zip(lengths, ext_accuracies, mlp_accuracies), key=lambda x: x[0][1])  # Sort by max_len
        sorted_lengths, sorted_ext_accuracies, sorted_mlp_accuracies = zip(*sorted_data)
        
        # Create formatted x-axis labels (min_len, max_len) after sorting
        sorted_labels = [f'({min_len},{max_len})' for min_len, max_len in sorted_lengths]
        
        # Plot the accuracies for both models
        plt.plot(sorted_labels, sorted_ext_accuracies, marker=markers[0], color=colors[0], linestyle='-', 
                 label=f'EXT (Dataset size {dataset_size})')
        plt.plot(sorted_labels, sorted_mlp_accuracies, marker=markers[1], color=colors[1], linestyle='--', 
                 label=f'MLP (Dataset size {dataset_size})')
    
    # Customize the plot
    plt.xlabel('Motif Length (min_len, max_len)')
    plt.ylabel('Accuracy')
    plt.title('Accuracy of EXT and MLP Models (Combined Length)')
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.ylim(0, 1)  # Set y-axis limits for accuracy
    plt.grid(True)
    plt.legend()
    
    # Adjust layout to prevent cutoff of labels
    plt.tight_layout()
    
    # Save or display the plot
    if save_as_file:
        plt.savefig(filename)
        plt.close()  # Close the plot to avoid display when saving
    else:
        plt.show()  # Display the plot if not saving

import matplotlib.pyplot as plt
from typing import Dict, Any

def plot_total_time_combined_len(results, save_as_file=True, filename='Total_time_combined_len_EXT_MLP.png'):
    """
    Plot the total execution time of EXT and MLP models for combined-length motifs for various dataset sizes.

    Args:
    results (Dict[str, Any]): A dictionary containing experiment results.
    save_as_file (bool): If True, save the plot as a file. If False, display the plot.
    filename (str): The filename to save the plot if save_as_file is True.
    """
    plt.figure(figsize=(12, 6))
    
    # Extract dataset sizes from the results dictionary for 'combined_length'
    combined_length_results = results.get('combined_length', {})
    
    models = ['EXT', 'MLP']
    markers = ['o', 's']  # Different markers for EXT and MLP
    colors = ['blue', 'red']  # Different colors for EXT and MLP
    
    for dataset_size, lengths_data in combined_length_results.items():
        lengths = []
        ext_total_times = []
        mlp_total_times = []
        
        # Extract motif (min_len, max_len) and corresponding total execution times for each dataset
        for (min_len, max_len), data in lengths_data.items():
            lengths.append((min_len, max_len))  # Store as a tuple for sorting
            
            # Calculate total execution time for EXT and MLP
            #print('Calculate total execution time for EXT and MLP')
            #print('data[EXT][training_time_sec] = ', data['EXT']['training_time_sec'])
            cpp_time = data.get('cpp_execution_time_sec', 0)
            ext_training_time = data.get('EXT', {}).get('training_time_sec', 0)
            ext_testing_time = data.get('EXT', {}).get('testing_time_sec', 0)
            mlp_training_time = data.get('MLP', {}).get('training_time_sec', 0)
            mlp_testing_time = data.get('MLP', {}).get('testing_time_sec', 0)
            #print(f'cpp_time : {cpp_time} type : {type(cpp_time)}')
            #print(f'ext_training_time : {ext_training_time} type : {type(ext_training_time)}')
            #print(f'ext_testing_time : {ext_testing_time} type : {type(ext_testing_time)}')


            ext_total_time = cpp_time + ext_training_time + ext_testing_time
            mlp_total_time = cpp_time + mlp_training_time + mlp_testing_time
            
            #print(f' ext_total_time = {ext_total_time}...')

            ext_total_times.append(ext_total_time)
            mlp_total_times.append(mlp_total_time)
        
        # Sort the lengths and total times by max_len, since min_len is always the same
        sorted_data = sorted(zip(lengths, ext_total_times, mlp_total_times), key=lambda x: x[0][1])  # Sort by max_len
        sorted_lengths, sorted_ext_total_times, sorted_mlp_total_times = zip(*sorted_data)
        
        # Create formatted x-axis labels (min_len, max_len) after sorting
        sorted_labels = [f'({min_len},{max_len})' for min_len, max_len in sorted_lengths]
        
        # Plot the total execution times for both models
        plt.plot(sorted_labels, sorted_ext_total_times, marker=markers[0], color=colors[0], linestyle='-', 
                 label=f'EXT (Dataset size {dataset_size})')
        plt.plot(sorted_labels, sorted_mlp_total_times, marker=markers[1], color=colors[1], linestyle='--', 
                 label=f'MLP (Dataset size {dataset_size})')
    
    # Customize the plot
    plt.xlabel('Motif Length (min_len, max_len)')
    plt.ylabel('Total Execution Time (seconds)')
    plt.title('Total Execution Time of EXT and MLP Models (Combined Length)')
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.grid(True)
    plt.legend()
    
    # Adjust layout to prevent cutoff of labels
    plt.tight_layout()
    
    # Save or display the plot
    if save_as_file:
        plt.savefig(filename)
        plt.close()  # Close the plot to avoid display when saving
    else:
        plt.show()  # Display the plot if not saving


def plot_all_len_motifs_exp(results):
    
    print('plot_growth_nb_motifs_fixed_len: ...')
    plot_growth_nb_motifs_fixed_len(results, save_as_file=True)
    print('plot_growth_nb_motifs_combined_len: ...')
    plot_growth_nb_motifs_combined_len(results, save_as_file=True)

    print(' Accuracy for models ...')
    plot_accuracy_fixed_len(results)

    print(' Accuracy combined len ...')
    plot_accuracy_combined_len(results)

    print(" Total time for combined len...")
    plot_total_time_combined_len(results)

    print(' Finish exp len motifs :) ')


# -------------------------------------------
# -------------------------------------------
# Experiments on: Beta filtering
# Filtering according to motif conservation in their family (parameter $\alpha$)
# Percentage of given cm in the family (until now we call it Beta), but now in previous phrase (from paper we call it alpha)
# because in the beginning we call it beta (and also in code), then later in the paper we changed to alpha.
# So, as parameter we still call it beta.
# in the figure we call it alpha.
#
# 
# In this experiment, we compared different settings including the parameter $\alpha$. Figure \ref{fig:beta_max_10} and Supplementary Figure \ref{fig:beta_max_10_s} show the evolution of the dataset size, the classification accuracy, the training time and the data generation time for different values of $\alpha$.
# \begin{figure}[b]
#    \centering
#    %\includegraphics[width=.49\textwidth]{images/Beta_0_Total_motifs.PNG}\hfill
#    %\includegraphics[width=.49\textwidth]{images/Beta_0_Accuracy.PNG}
#    %~\\
#    %\includegraphics[width=.49\textwidth]{images/Beta_gamma_Total_motifs.PNG}\hfill
#    \includegraphics[width=.49\textwidth]{images/Beta_gamma_Accuracy.PNG}
#    
#    \caption{Relation between $\alpha$ and the classification accuracy for EXT and MLP with  $\gamma = 1$ or $\gamma = 2$.}
#    
#    \label{fig:beta_gamma_test}
#\end{figure}
# need to show : Beta_gamma_Accuracy.PNG
# datasets used: 
# beta [0,50,100] (cm percentage in family), in the paper called alpha.
# alpha: -1 (in paper called beta)
# gamma = [1, 2]
# model: [EXT, MLP]
# min_len = 2
# max_len = 10
# data-size (nb families): 500
# No deletion

def run_beta_experiments(is_debug_datasets: bool = False) -> None:
    """
    Run experiments to analyze the effect of beta (percentage of common motifs) and gamma (minimum occurrences) parameters.

    :param is_debug_datasets: If True, use debug dataset size
    """
    # Step 1: Define parameters
    dataset_size: int = 500 if not is_debug_datasets else debug_datasets_size_global_var[-1]
    min_len: int = 2
    max_len: int = 10
    list_betas: list[int] = [0, 50, 100]  # percentage of common motifs in the family
    alpha: int = -1  # no filtering
    list_gammas: list[int] = [1, 2]
    is_sub_motif_delete: int = 0  # false, no deletion of submotifs

    # Initialize results dictionary with nested structure
    results: dict[int, dict[int, Optional[dict[str, Any]]]] = {
        beta: {gamma: None for gamma in list_gammas} for beta in list_betas
    }

    # Step 2: Run experiments
    for beta in list_betas:
        print(f"Running experiments for beta {beta}")
        for gamma in list_gammas:
            print(f"  Running experiment for gamma {gamma}")
            results[beta][gamma] = run_single_experiment(
                dataset_size=dataset_size,
                min_length=min_len,
                max_length=max_len,
                beta_percentage=beta,
                alpha_nb_occ_variation=alpha,
                gamma_nb_occrs=gamma,
                is_delete_submotifs=is_sub_motif_delete,
                test_name='beta',
                is_debug_datasets=is_debug_datasets
            )

    # Step 3: Generate plots
    plot_accuracy_beta_gamma(results, 'beta_gamma_accuracy_plot.png')
    plot_motifs_beta_gamma(results, 'beta_gamma_motifs_plot.png')


    # Step 4: Save results (optional, but recommended)
    save_results_to_file(results, "beta_experiments_results.json")

def save_results_to_file(results: dict, filename: str) -> None:
    """
    Save the results dictionary to a JSON file.

    :param results: The results dictionary to save
    :param filename: The name of the file to save the results to
    """
    import json
    with open(filename, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {filename}")


# Example of combined results dictionary
results_debug_beta_gamma = {
    0: {
        1: {
            'EXT': {'accuracy': 0.85},
            'MLP': {'accuracy': 0.82},
            'num_motifs': 150
        },
        2: {
            'EXT': {'accuracy': 0.87},
            'MLP': {'accuracy': 0.84},
            'num_motifs': 120
        }
    },
    50: {
        1: {
            'EXT': {'accuracy': 0.89},
            'MLP': {'accuracy': 0.86},
            'num_motifs': 200
        },
        2: {
            'EXT': {'accuracy': 0.91},
            'MLP': {'accuracy': 0.88},
            'num_motifs': 180
        }
    },
    100: {
        1: {
            'EXT': {'accuracy': 0.92},
            'MLP': {'accuracy': 0.90},
            'num_motifs': 250
        },
        2: {
            'EXT': {'accuracy': 0.94},
            'MLP': {'accuracy': 0.92},
            'num_motifs': 220
        }
    }
}

# Now let's update both plotting functions to work with this new structure
import numpy as np

def plot_accuracy_beta_gamma(results, output_filename='accuracy_plot.png'):
    betas = list(results.keys())
    gammas = list(results[betas[0]].keys())
    models = ['EXT', 'MLP']
    
    x = np.arange(len(betas))
    width = 0.1
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for i, model in enumerate(models):
        for j, gamma in enumerate(gammas):
            accuracies = [results[beta][gamma][model]['accuracy'] for beta in betas]
            offset = width * (i * len(gammas) + j - (len(models) * len(gammas) - 1) / 2)
            rects = ax.bar(x + offset, accuracies, width, label=f'{model}, γ={gamma}')
            ax.bar_label(rects, fmt='{:.2f}', padding=3, rotation=90, fontsize=8)

    ax.set_xlabel('Beta (%)')
    ax.set_ylabel('Accuracy')
    ax.set_title('Accuracy for EXT and MLP Models with Different Beta and Gamma Values')
    ax.set_xticks(x)
    ax.set_xticklabels(betas)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")

def plot_motifs_beta_gamma(results, output_filename='motifs_plot.png'):
    betas = list(results.keys())
    gammas = list(results[betas[0]].keys())
    
    x = np.arange(len(betas))
    width = 0.2
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for j, gamma in enumerate(gammas):
        num_motifs = [results[beta][gamma]['num_motifs'] for beta in betas]
        offset = width * (j - len(gammas)/2 + 0.5)
        rects = ax.bar(x + offset, num_motifs, width, label=f'γ={gamma}')
        ax.bar_label(rects, fmt='{:.0f}', padding=3, rotation=90, fontsize=8)

    ax.set_xlabel('Beta (%)')
    ax.set_ylabel('Number of Motifs')
    ax.set_title('Number of Motifs for Different Beta and Gamma Values')
    ax.set_xticks(x)
    ax.set_xticklabels(betas)
    ax.legend()
    
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")


# -------------------------------------------
# -------------------------------------------
# Experiments on: Alpha filtering (in the paper called beta :p)
# D. Filtering according to the number of occurrences variation (parameter $\beta$)
# the number of motifs for $\beta = 0$ and $\beta = +\infty$ for different values of $\alpha$.
# so, in the code
# beta (paper alpha) = 0, 50, 100
# alpha (paper beta) = 0 , +infty 
# gamma = 1 # no filtering
# model: [EXT, MLP]
# min_len = 2
# max_len = 10
# data-size (nb families): 500
# No deletion


def run_alpha_variance_experiments(is_debug_datasets: bool = False) -> None:
    """
    Run experiments to analyze the effect of beta (percentage of common motifs) and gamma (minimum occurrences) parameters.

    :param is_debug_datasets: If True, use debug dataset size
    """
    # Step 1: Define parameters
    dataset_size: int = 500 if not is_debug_datasets else debug_datasets_size_global_var[-1]
    min_len: int = 2
    max_len: int = 10
    list_betas: list[int] = [0, 50, 100]  # percentage of common motifs in the family
    list_alphas: list[int] = [1,100]  # 1 and 100 (as +infiny :p) since all test from 10, 20, ...and on give the same results. # We can Max positive value for 32-bit signed integer, but it give pbm, I don't know why !!!, but since we don't need it, so we don't care :p
    gamma: int = 1 # no filtering
    is_sub_motif_delete: int = 0  # false, no deletion of submotifs

    # Initialize results dictionary with nested structure
    results: dict[int, dict[int, Optional[dict[str, Any]]]] = {
        beta: {alpha: None for alpha in list_alphas} for beta in list_betas
    }

    # Step 2: Run experiments
    for beta in list_betas:
        print(f"Running experiments for beta {beta}")
        for alpha in list_alphas:
            print(f"  Running experiment for alpha {alpha}")
            results[beta][alpha] = run_single_experiment(
                dataset_size=dataset_size,
                min_length=min_len,
                max_length=max_len,
                beta_percentage=beta,
                alpha_nb_occ_variation=alpha,
                gamma_nb_occrs=gamma,
                is_delete_submotifs=is_sub_motif_delete,
                test_name='alpha',
                is_debug_datasets=is_debug_datasets
            )

    # Step 3: Generate plots
    plot_accuracy_beta_alpha(results, 'beta_alpha_accuracy_plot.png')
    plot_motifs_beta_alpha(results, 'beta_alpha_motifs_plot.png')


    # Step 4: Save results (optional, but recommended)
    save_results_to_file(results, "alpha_experiments_results.json")


def plot_accuracy_beta_alpha(results, output_filename='accuracy_plot.png'):
    betas = list(results.keys())
    alphas = list(results[betas[0]].keys())
    models = ['EXT', 'MLP']
    
    x = np.arange(len(betas))
    width = 0.1
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for i, model in enumerate(models):
        for j, alpha in enumerate(alphas):
            accuracies = [results[beta][alpha][model]['accuracy'] for beta in betas]
            offset = width * (i * len(alphas) + j - (len(models) * len(alphas) - 1) / 2)
            alpha_label = "1" if alpha == 1 else "∞"
            rects = ax.bar(x + offset, accuracies, width, label=f'{model}, α={alpha_label}')
            ax.bar_label(rects, fmt='{:.2f}', padding=3, rotation=90, fontsize=8)

    ax.set_xlabel('Beta (%)')
    ax.set_ylabel('Accuracy')
    ax.set_title('Accuracy for EXT and MLP Models with Different Beta and Alpha Values')
    ax.set_xticks(x)
    ax.set_xticklabels(betas)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")

def plot_motifs_beta_alpha(results, output_filename='motifs_plot.png'):
    betas = list(results.keys())
    alphas = list(results[betas[0]].keys())
    
    x = np.arange(len(betas))
    width = 0.2
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for j, alpha in enumerate(alphas):
        num_motifs = [results[beta][alpha]['num_motifs'] for beta in betas]
        offset = width * (j - len(alphas)/2 + 0.5)
        alpha_label = "1" if alpha == 1 else "∞"
        rects = ax.bar(x + offset, num_motifs, width, label=f'α={alpha_label}')
        ax.bar_label(rects, fmt='{:.0f}', padding=3, rotation=90, fontsize=8)

    ax.set_xlabel('Beta (%)')
    ax.set_ylabel('Number of Motifs')
    ax.set_title('Number of Motifs for Different Beta and Alpha Values')
    ax.set_xticks(x)
    ax.set_xticklabels(betas)
    ax.legend()
    
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")

# -------------------------------------------
# -------------------------------------------
# Experiments on: cla algs choices
# Fixed len: 5, 7
#     gamma = 1, 2
#     => classification time
#     => accuracy
# combined len: (2,5), (2,8)
#     gamma = 1, 2
#     => classification time
#     => accuracy
# size = 500
# beta (paper alpha): 50
# alpha variance : -1 # no filtering
# no deletion

def run_algs_choice_experiments(is_debug_datasets: bool = False) -> None:
    # Step 1: Define parameters
    dataset_size: int = 500
    
    if is_debug_datasets:
        dataset_sizes = debug_datasets_size_global_var[-1] 

    beta: int= 50 # percentage of the cm in the family (called alpha in the paper)
    alpha: int= -1 # variance of nb cm between seqs in family , -1 we don't care
    list_gammas: list[int]= [1, 2] # min nb of cm to consider 
    is_sub_motif_delete: int= 0 # delete of note sub motif, default: 0 so no.

    fixed_lengths: list[int] = [5, 7]
    combined_lengths: list[tuple[int, int]] = [(2, 5), (2, 8)]

    # note: when we have min_len = i, max_len = i+1 ex:(2,3) this mean we have fixed motif len of len 2. (this is in my c++ code work like this, maybe I should take time to change it...)

    # Initialize results dictionary with nested structure
    results_fixed_len: dict[int, dict[int, Optional[dict[str, Any]]]] = {
        cm_len: {gamma: None for gamma in list_gammas} for cm_len in fixed_lengths
    }

    results_combiened_len: dict[int, dict[int, Optional[dict[str, Any]]]] = {
        (cm_min, cm_max): {gamma: None for gamma in list_gammas} for (cm_min, cm_max) in combined_lengths
    }


    # Step 2: Run experiments for fixed-length motifs
    for cm_len in fixed_lengths:
        print(f"Running experiments for fixed len {cm_len}")
        for gamma in list_gammas:
            print(f"  Running experiment for gamma {gamma}")
            results_fixed_len[cm_len][gamma] = run_single_experiment_models_choices(
                dataset_size=dataset_size,
                min_length=cm_len-1,
                max_length=cm_len,
                beta_percentage=beta,
                alpha_nb_occ_variation=alpha,
                gamma_nb_occrs=gamma,
                is_delete_submotifs=is_sub_motif_delete,
                test_name='algs',
                is_debug_datasets=is_debug_datasets
            )


    # Step 3: Run experiments for combined-length motifs (only for size 350)
    for cm_min, cm_max in combined_lengths:
        print(f"Running experiments for Combined len {cm_min}, {cm_max}")
        for gamma in list_gammas:
            print(f"  Running experiment for gamma {gamma}")
            results_combiened_len[(cm_min, cm_max)][gamma] = run_single_experiment_models_choices(
                dataset_size=dataset_size,
                min_length=cm_min,
                max_length=cm_max,
                beta_percentage=beta,
                alpha_nb_occ_variation=alpha,
                gamma_nb_occrs=gamma,
                is_delete_submotifs=is_sub_motif_delete,
                test_name='algs',
                is_debug_datasets=is_debug_datasets
            )
    

    # Step 4: Generate plots
    


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
        os.path.join(cpp_dir, "CommonMotifs.cpp"),
        os.path.join(cpp_dir, "SequenceIdManager.cpp"),
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

    print(command)
    print(" ".join(command))

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



# -------------------------------------------
# -------------------------------------------
# run the classification experiments

def run_classification_experiment(
    model: str,
    train_csv: str,
    test_dir: str,
    file_ext: str = ".fasta.txt",
    n_job: int = 1,
    output_csv: Optional[str] = None
) -> Dict[str, Any]:
    """
    Runs a classification experiment using Main.py and retrieves results from the output CSV.

    :param model: Machine learning model name (e.g., 'EXT', 'MLP', 'VOT', 'RDF')
    :param train_csv: Path to the training CSV matrix file
    :param test_dir: Path to the test files directory
    :param file_ext: File extension for test files
    :param n_job: Number of parallel jobs
    :param output_csv: Path to the output CSV file for results. If None, a default path is used.
    :return: Dictionary containing the experiment results
    """
    # Define project root and Main.py path
    project_root = find_project_root()
    classification_script = project_root / "LSC-ncRNA-our_method/Classification/Main.py"
    
    # Define default output_csv if not provided
    if output_csv is None:
        output_csv = Path("results") / f"classification_results_{model}_{os.getpid()}.csv"
    
    # Ensure output directory exists
    output_csv_path = Path(output_csv)
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    #print(' inside run classification ...')
    #print('train_csv: ', train_csv)
    #print('output_csv_path: ', output_csv_path)
    #time.sleep(3)

    # Command to run Main.py
    command = [
        "python",
        str(classification_script),
        "-m", model,
        "-t", train_csv,
        "-d", test_dir,
        "-e", file_ext,
        "-j", str(n_job),
        "-o", str(output_csv)
    ]
    
    print(f"Running classification experiment for model: {model}")
    
    try:
        # Execute the command
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print(f"Experiment for model {model} completed successfully.")
        
    except subprocess.CalledProcessError as e:
        print(f"Error running classification experiment for model {model}: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return {"error": f"Subprocess failed: {e.stderr}"}
    
    #time.sleep(3)
    #print('Read the output CSV and extract results...')
    # Read the output CSV and extract results
    if output_csv_path.exists():
        try:
            df = pd.read_csv(output_csv_path)
            
            # Extract necessary metrics
            # Adjust the column names based on your CSV structure
            # Assuming columns like 'Training_Time_sec', 'Testing_Time_sec', 'Accuracy', etc.
            training_time = df.get('Training_Time_sec', [None])[0]
            testing_time = df.get('Testing_Time_sec', [None])[0]
            accuracy = df.get('Accuracy', [None])[0]
            precision = df.get('Precision', [None])[0]
            recall = df.get('Recall', [None])[0]
            f1_score = df.get('F1_Score', [None])[0]
            num_test_sequences = df.get('Number_Test_Sequences', [None])[0]
            num_total_motifs = df.get('Number_Total_Motifs', [None])[0]
            num_tested_classes = df.get('Number_Tested_Classes', [None])[0]
            
            return {
                "model": model,
                "training_time_sec": training_time,
                "testing_time_sec": testing_time,
                "accuracy": accuracy,
                "precision": precision,
                "recall": recall,
                "f1_score": f1_score,
                "num_test_sequences": num_test_sequences,
                "num_total_motifs": num_total_motifs,
                "num_tested_classes": num_tested_classes
            }
        
        except Exception as e:
            print(f"Error reading classification output CSV for model {model}: {e}")
            return {"error": f"CSV read failed: {str(e)}"}
    else:
        print(f"Output CSV file not found: {output_csv_path}")
        return {"error": "Output CSV file not found."}


def run_classification_experiment_models_choices(
    train_csv: str,
    output_csv: Optional[str] = None
) -> Dict[str, Any]:
    """
    Runs a classification experiment using modelstest.py and retrieves results from the output CSV.
    Returns processing time and accuracy for each model.

    Args:
        train_csv: Path to the training CSV matrix file
        output_csv: Path to the output CSV file for results. If None, a default path is used.
    
    Returns:
        Dictionary containing the experiment results for each model
    """
    # Define project root and modelstest.py path
    project_root = find_project_root()
    classification_script = project_root / "LSC-ncRNA-our_method/Classification/modelstest.py"
    
    # Define default output_csv if not provided
    if output_csv is None:
        output_csv = Path("results") / f"classification_results_models_choices_{os.getpid()}.csv"
    
    # Ensure output directory exists
    output_csv_path = Path(output_csv)
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Command to run modelstest.py
    command = [
        "python",
        str(classification_script),
        train_csv,
        str(output_csv)
    ]
    
    try:
        # Execute the command
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
    except subprocess.CalledProcessError as e:
        print(f"Error running classification experiment: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return {"error": f"Subprocess failed: {e.stderr}"}
    
    # Read the output CSV and extract results
    if output_csv_path.exists():
        try:
            df = pd.read_csv(output_csv_path)
            
            # Initialize results dictionary
            results = {
                "models": [],
                "processing_times": [],
                "accuracies": []
            }
            
            # Extract results for each model
            for _, row in df.iterrows():
                model_name = row['model']
                processing_time = row['execution_time']  # from the original CSV
                accuracy = row['Score Test']  # from the original CSV
                
                results["models"].append(model_name)
                results["processing_times"].append(processing_time)
                results["accuracies"].append(accuracy)
            
            return results
        
        except Exception as e:
            print(f"Error reading classification output CSV: {e}")
            return {"error": f"CSV read failed: {str(e)}"}
    else:
        print(f"Output CSV file not found: {output_csv_path}")
        return {"error": "Output CSV file not found."}



import csv
from typing import Dict, Any, Union

def write_results_to_csv(results: Dict[str, Any], csv_filename: str) -> None:
    """
    Write experiment results to a CSV file.

    Args:
    results (Dict[str, Any]): A dictionary containing experiment results.
    csv_filename (str): The name of the CSV file to write to.

    Raises:
    IOError: If there's an issue writing to the file.
    """
    try:
        with open(csv_filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            write_header(writer)
            write_fixed_length_results(writer, results.get('fixed_length', {}))
            write_combined_length_results(writer, results.get('combined_length', {}))
    except IOError as e:
        print(f"Error writing to file {csv_filename}: {e}")

def write_header(writer: csv.writer) -> None:
    writer.writerow([
        'Motif Type', 'Dataset Size', 'CM Min', 'CM Max', 'Number motifs',
        'Time cpp', 'File size GB', 'Model', 'Accuracy', 'Precision',
        'Recall', 'F1 Score', 'Training Time (sec)', 'Testing Time (sec)'
    ])

def write_fixed_length_results(writer: csv.writer, fixed_length_results: Dict[str, Any]) -> None:
    for size, cm_lengths in fixed_length_results.items():
        for cm_length, experiment_results in cm_lengths.items():
            write_experiment_results(writer, 'fixed_length', size, cm_length, cm_length, experiment_results)

def write_combined_length_results(writer: csv.writer, combined_length_results: Dict[str, Any]) -> None:
    for size, length_ranges in combined_length_results.items():
        for (min_len, max_len), experiment_results in length_ranges.items():
            write_experiment_results(writer, 'combined_length', size, min_len, max_len, experiment_results)

def write_experiment_results(writer: csv.writer, motif_type: str, size: Union[int, str], 
                             min_len: Union[int, str], max_len: Union[int, str], 
                             experiment_results: Dict[str, Any]) -> None:
    for model, metrics in experiment_results.items():
        if isinstance(metrics, dict) and model not in [
            'dataset_size', 'min_length', 'max_length',
            'num_motifs', 'cpp_execution_time_sec', 'file_size_gb']:
            writer.writerow([
                motif_type, size, min_len, max_len,
                experiment_results.get('num_motifs', 'N/A'),
                experiment_results.get('cpp_execution_time_sec', 'N/A'),
                experiment_results.get('file_size_gb', 'N/A'),
                model,
                metrics.get('accuracy', 'N/A'),
                metrics.get('precision', 'N/A'),
                metrics.get('recall', 'N/A'),
                metrics.get('f1_score', 'N/A'),
                metrics.get('training_time_sec', 'N/A'),
                metrics.get('testing_time_sec', 'N/A')
            ])
# Example usage:
# write_results_to_csv(experiment_results, 'experiment_results.csv')


def debug_cpp_fixed_len():
    
    # ibra todo:
    # I checked but didn't find the pbm
    # I know that the input min_len will be converted to min_len+1
    # so, after the main
    # we can put min_len = arg.min_len -1
    # like this after it will be chnaged to min_len+1 again, and work as desired.
    # Now, I don't find where this is donne, maybe in loop in AddStringMinMax...
    # but anyway, this will solve the pbm.
    

        # Get the project root directory
    project_root: Path = find_project_root()

    is_debug_datasets = True
    dataset_sizes = [5]


    test_name: str = f"fixed"  # test delete no delete
    min_length: int = 2
    max_length: int = 4
    beta: int = 0  # don't use beta (0 is the smallest value)
    alpha: int = -1  # don't use alpha
    gamma: int = 1  # 1 is smallest value

    results: list[dict[str, int | float | str]] = []

    # Step 2: Iterate through dataset sizes and submotif deletion options
    for size in dataset_sizes:
        # Step 3: Set up input file paths
        dir_path = get_dir_path_rfam14_sample(size=size, data_type='Train', is_debug_datasets=is_debug_datasets)

        for is_delete_submotifs in [0]:
            
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



def main():

    is_debug_datasets: bool = is_debug_datasets_global_var

    # unzip already pre-prepared datasets
    ## prepare_dataset()

    # compile the c++ code for motifs extraction and selection
    compile_code_MotifsExtractionSelection()
    
    # debug min max length
    debug_cpp_fixed_len()

    # run the sub motifs deletion experiments
    ## deletion_sub_motifs(is_debug_datasets)
    
    # run the motifs length experiments
    ## print("run the motifs length experiments...")
    ## run_motif_length_experiments(is_debug_datasets)

    #
    ## print("run beta gamma experiments...")
    ## run_beta_experiments(is_debug_datasets)
    # debug by internal dict 
    ## plot_all_len_motifs_exp(dict_len_motifs_exp_example_debug)

    ## print("run beta alpha experiments...")
    ## run_alpha_variance_experiments(is_debug_datasets)

    ## print('run algs models choices experiments...')


if __name__ == "__main__":
    main()
    