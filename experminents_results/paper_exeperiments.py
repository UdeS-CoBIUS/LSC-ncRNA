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


# datasets: ------------------------------
# Datasets are pre-prepared and available in the `datasets/data/` directory.
# there are 3 datasets:
# datasets/data/Clans_ncRNA_from_Rfam_14.8.zip
# datasets/data/deep_ncrna_datasets.zip
# datasets/data/Rfam_14.1_dataset.zip


def prepare_dataset():
    # Define the directory containing the datasets
    data_dir = "datasets/data/"
    
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
# We generate 2 figures:
# Evolution of the data size (Left figure) and processing time (Right figure) for the generation of datasetby filtering out exact submotifs\textbf{(F)}, and without any  filtering \textbf{(NF)}.

def deletion_sub_motifs():
    pass


def compile_code_MotifsExtractionSelection():
    cpp_dir = "LSC-ncRNA-our_method/MotifsExtractionSelection"
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
    # prepare_dataset() # unzip already pre-prepared datasets
    
    # compile the c++ code for motifs extraction and selection
    compile_code_MotifsExtractionSelection()

if __name__ == "__main__":
    main()

