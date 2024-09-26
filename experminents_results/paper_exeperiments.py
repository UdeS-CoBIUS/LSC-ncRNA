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

def find_project_root(project_name="LSC-ncRNA") -> Path:
    """
    Find the project root directory by looking for a specific directory name.
    :param project_name: The name of the project directory
    :return: The absolute path to the project root directory as a Path object
    """
    current_dir = Path(__file__).resolve().parent
    
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
    project_root = find_project_root()
    
    # Define the directory containing the datasets
    data_dir = project_root / "datasets/data"

    # List of dataset zip files
    datasets = [
        "Clans_ncRNA_from_Rfam_14.8.zip",
        "deep_ncrna_datasets.zip",
        "Rfam_14.1_dataset.zip"
    ]
    
    for dataset in datasets:
        zip_path = os.path.join(data_dir, dataset)
        extract_dir = os.path.join(data_dir, dataset[:-4])  # Remove .zip extension
        
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
 

def get_delete_sub_motifs_dir_path(size):
    # Get the project root directory
    project_root = find_project_root()
    
    # Check for double subfolder (when unzipped by Python)
    double_subfolder_path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test/Train"
    
    # Check for single subfolder (when unzipped by other apps)
    single_subfolder_path = project_root / f"datasets/data/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_{size}_Train_Test/Train"
    
    if double_subfolder_path.exists():
        return double_subfolder_path
    elif single_subfolder_path.exists():
        return single_subfolder_path
    else:
        raise FileNotFoundError(f"Directory not found for size {size}. Checked paths:\n{double_subfolder_path}\n{single_subfolder_path}")


def deletion_sub_motifs():
    # Get the project root directory
    project_root = find_project_root()
    
    # Step 1: Define parameters and paths
    cpp_dir = project_root / "LSC-ncRNA-our_method/MotifsExtractionSelection"
    executable = os.path.join(cpp_dir, "MotifsExtractionSelection")
    
    # check if the executable exists
    if not os.path.exists(executable):
        print(f"Error executable does not exists: {executable} not found in {cpp_dir}")
        return

    dataset_sizes = [100, 200, 300, 400, 500, 600]
    min_length, max_length = 2, 20
    beta = 0 # don't use beta (0 is the smallest value)
    alpha = -1 # don't use alpha
    gamma = 1 # 1 is smallest value
    
    results = []

    # Step 2: Iterate through dataset sizes and submotif deletion options
    for size in dataset_sizes:
        for is_delete_submotifs in [0, 1]:
            # Step 3: Set up input file paths
            dir_path = get_delete_sub_motifs_dir_path(size)
            test_name = f"test_dnd" # test delete no delete
            
            # Step 4: Prepare the command for running MotifsExtractionSelection
            command = [
                executable,
                "-in", dir_path, # dir_name
                "-tn", test_name,
                "-nf", str(size),
                "-minl", str(min_length),
                "-maxl", str(max_length),
                "-d", str(is_delete_submotifs),
                "-b", str(beta),
                "-a", str(alpha),
                "-g", str(gamma)
            ]

            # get csv output file name
            args = {
                "-tn": test_name,
                "-nf": size,
                "-minl": min_length,
                "-maxl": max_length,
                "-d": is_delete_submotifs,
                "-b": beta,
                "-a": alpha,
                "-g": gamma
            }

            output_csv_file = generate_csv_output_filename(args)

            try:
                print(f'python file name: {output_csv_file}')
                # Step 5: Run the C++ program and capture its output
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                stdout, stderr = process.communicate()
                print(f'stdout: {stdout}')

                # Step 6: Extract execution time from stdout : in the program we have ( cout << "Time taken by whole program is : " << fixed << time_taken << setprecision(9) << " sec" << endl;)
                time_match = re.search(r"Time taken by whole program is : (\d+\.\d+) sec", stdout)
                execution_time = float(time_match.group(1)) if time_match else None

                # Step 7: Get the size of the generated CSV file
                file_size = os.path.getsize(output_csv_file) / (1024 * 1024 * 1024)  # Convert to GB

                # Step 8: Store the results
                results.append({
                    "dataset_size": size,
                    "is_delete_submotifs": bool(is_delete_submotifs),
                    "execution_time": execution_time,
                    "file_size_gb": file_size
                })

                print(f"Processed dataset size {size} with submotif deletion {'enabled' if is_delete_submotifs else 'disabled'}")

                # Step 9: Save results to a CSV file
                results_dir = "results"
                os.makedirs(results_dir, exist_ok=True)  # Create the directory if it doesn't exist
                csv_path = os.path.join(results_dir, "deletion_sub_motifs_results.csv")
                
                with open(csv_path, "w", newline="") as csvfile:
                    fieldnames = ["dataset_size", "is_delete_submotifs", "execution_time", "file_size_gb"]
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for row in results:
                        writer.writerow(row)
                print(f"Results saved to {csv_path}")

                # Step 10: Clean up the generated CSV file, since it will have big size
                os.remove(output_csv_file)
                print(f"Cleaned up {output_csv_file}")
            
            except Exception as e:
                print(f"Error processing dataset size {size} with submotif deletion {'enabled' if is_delete_submotifs else 'disabled'}: {str(e)}")    
    
    # generate the plot of the results in 2 figures:
    # we compare 2 methods:
    # A- filtering out exact submotifs \textbf{(F)}, 
    # B- Without any  filtering \textbf{(NF)}
    # 1. generate Evolution of the data size for F vs NF:  
    # 2. generate evolution of processing time in minute  for F vs NF
    
    # Call the function to generate plots
    generate_result_sub_motifs_F_vs_NF_plots(csv_path)

def generate_result_sub_motifs_F_vs_NF_plots(csv_file_path):
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file_path)

    # Plot 1: Evolution of the data size for F vs NF
    plt.figure(figsize=(10, 6))
    for is_delete, group in df.groupby('is_delete_submotifs'):
        label = 'F' if is_delete else 'NF'
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
        label = 'F' if is_delete else 'NF'
        plt.plot(group['dataset_size'], group['execution_time'] / 60, marker='o', label=label)

    plt.xlabel('Dataset Size')
    plt.ylabel('Processing Time (minutes)')
    plt.title('Evolution of Processing Time')
    plt.legend()
    plt.grid(True)
    plt.savefig('results/evolution_of_processing_time.png')
    plt.close()
    print("Plot saved as results/evolution_of_processing_time.png")



# compile code_MotifsExtractionSelection: ------------------------------
# Struct to store the parsed command-line arguments
#struct Args {

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
#};
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
        "-std=c++17",
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



def main():
    prepare_dataset() # unzip already pre-prepared datasets
    
    # compile the c++ code for motifs extraction and selection
    compile_code_MotifsExtractionSelection()

    deletion_sub_motifs()

if __name__ == "__main__":
    main()

