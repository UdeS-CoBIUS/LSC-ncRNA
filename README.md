# LSC-ncRNA: Large-Scale Classification of non-coding RNA

LSC-ncRNA is a sequence-based method for classifying ncRNA families. It builds a vector representation of sequences from common motifs (computed/selected by C++), then trains supervised models (Python) to predict families at scale.

This repository contains all source code plus a single Python pipeline that can reproduce the paper’s experiments end-to-end (on a small “debug” sample locally; full datasets require a server with large RAM).

- Paper pipeline script: `experminents_results/paper_exeperiments.py`
- C++ motif extractor: `LSC-ncRNA-our_method/MotifsExtractionSelection`
- Python classifiers: `LSC-ncRNA-our_method/Classification`


## Table of contents

1. Quickstart (small demo on a laptop)
2. Full reproduction (paper experiments)
3. Requirements and installation
4. Build the C++ motif extractor
5. Datasets (prepared and from-scratch options)
6. Running experiments and reproducing figures
7. Standalone classification CLI
8. Compared methods (BLAST, Infernal)
9. Project structure
10. Troubleshooting and tips
11. License, authors, citation


## 1) Quickstart: small demo on a laptop

This runs a tiny “debug” experiment to verify your setup. It compiles the C++ extractor, generates a small motif matrix, and runs classification, producing sample plots in `experminents_results/`.

### Prerequisites
- C++14 toolchain (g++)
- Python 3.9+ recommended (3.8+ works); pip

### Steps
1. Create and activate a virtual environment
     - macOS/Linux:
         ```sh
         python3 -m venv .venv
         source .venv/bin/activate
         ```
     - Windows (PowerShell):
         ```ps1
         py -3 -m venv .venv
         .venv\Scripts\Activate.ps1
         ```

2. Install Python packages (see notes in section 3 for platform tips)
     ```sh
     pip install --upgrade pip
     pip install biopython datatable scikit-learn pandas numpy matplotlib portalocker pyahocorasick suffix-trees
     ```

3. Build the C++ motif extractor (details in section 4)

4. Run the pipeline in debug mode
     ```sh
     python experminents_results/paper_exeperiments.py
     ```
     By default, the script compiles the C++ code and runs a tiny fixed-length experiment (`debug_cpp_fixed_len`). Results and example figures are saved under `experminents_results/` and `results/`.


## 2) Full reproduction: paper experiments

Running the full Rfam-scale experiments requires a machine with substantial memory and disk (tens of GB RAM recommended). The pipeline functions are all implemented in `experminents_results/paper_exeperiments.py`. You can enable experiments by editing the `main()` function (uncomment the calls you need) or by importing and calling functions directly from a Python session.

### Key experiments and where they’re implemented
- Submotif deletion impact: `deletion_sub_motifs()`
- Number of motifs vs accuracy (fixed- and combined-length): `run_motif_length_experiments()`
- Conservation and occurrence filtering sweeps:
    - `run_same_family_threshold_experiments()` (called `run_beta_experiments` in older versions)
    - `run_occurrence_variation_experiments()` (called `run_alpha_variance_experiments` in older versions)
- Algorithm choices and timings: `run_algs_choice_experiments()`

### Outputs
- CSV summaries in `results/` or `experminents_results/results/`
- Figures in `experminents_results/` (for example: `Growth_NB_motifs.png`, `Accuracy_combined_len_EXT_MLP.png`, `beta_gamma_accuracy_plot.png`)


## 3) Requirements and installation

### Core requirements
- C++14 compiler (g++)
- Python 3.8+ (3.9+ recommended)

### Python packages (install via pip)
- `biopython`
- `datatable`
- `scikit-learn`
- `pandas`
- `numpy`
- `matplotlib`
- `portalocker`
- `pyahocorasick`
- `suffix-trees`
- Optional: `xgboost` (not required for the paper’s results). On macOS you may need `brew install libomp` first.

### Platform notes
- **macOS**: if `pip install datatable` fails, install from source:
    ```sh
    pip install git+https://github.com/h2oai/datatable
    ```
- **Windows**: we recommend using WSL2 or a recent MSVC build tools setup for C++14.


## 4) Build the C++ motif extractor

The C++ program computes and selects sequence motifs, creating a vectorial representation of ncRNA sequences. The result is a 2D matrix (saved as a CSV file) where each row is a sequence and each column is the occurrence count of a selected motif.

The Python pipeline can compile the C++ source for you, or you can build it manually.

### Manual Compilation
```sh
cd LSC-ncRNA-our_method/MotifsExtractionSelection
g++ -std=c++14 
    -o MotifsExtractionSelection 
    Main.cpp SuffixTree_QuadraticTime.cpp 
    FastaFilesReader.cpp CommonMotifs.cpp SequenceIdManager.cpp
```
This produces the executable `LSC-ncRNA-our_method/MotifsExtractionSelection/MotifsExtractionSelection`, which is used by the Python pipeline.

### Command-Line Usage
The program can be used like this:
```sh
./MotifsExtractionSelection -in <string> [-nf <integer> -mins <integer> -maxs <integer> -minl <integer> -maxl <integer> -d <integer> -b <integer> -a <integer> -g <integer> -tn <string>]
```

**Parameters:**
- `-in`: `<string>` **(Required)** Path to the directory of FASTA files.
- `-tn`: `<string>` Experiment name. Give a specific name to your experiment (default: "test").
- `-nf`: `<integer>` Number of families (default: 10).
- `-mins`: `<integer>` Minimum number of sequences per family (default: 4).
- `-maxs`: `<integer>` Maximum number of sequences per family (default: 1000).
- `-minl`: `<integer>` Minimum length of motifs (default: 2).
- `-maxl`: `<integer>` Maximum length of motifs (default: 10).
- `-d`: `<integer>` Delete sub-motifs (0: false, 1: true) (default: 0).
- `-b`: `<integer>` `same_family_percentage_threshold` (legacy: beta). A motif must be present in at least this percentage of sequences within the same family. Range: [0, 100] (default: 40).
- `-a`: `<integer>` `occurrence_variation_tolerance` (legacy: alpha). Tolerance for variation in motif occurrences across different families. A value of -1 disables this filter (default: -1).
- `-g`: `<integer>` `min_occurrence_count` (legacy: gamma). Minimum number of occurrences for a motif to be considered (default: 1).

**Example:**
```shell
nohup ./MotifsExtractionSelection 
    -in "/data/ibra/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_600_Train_Test/Train" 
    -minl 2 -maxl 8 -b 50 -g 1 -tn F_600 > out_F_600 &
```

The output CSV file name is constructed from the parameters, for example: `del_No_nbF_F_600_min_2_max_8_sameFamilyPct_50_occVarTol_-1_nbOccrs_1.csv`.


## 5) Datasets

### Prepared datasets (recommended)
The repository includes pre-prepared datasets under `datasets/data/`. Unzip the archives you need (Rfam 14.1, deep ncRNA, Rfam 14.8 clans). Most paper experiments use the Rfam 14.1 set with predefined train/test splits.

### From scratch (optional)
If you wish to prepare the datasets from raw data, use the scripts located in `datasets/preparation/`.

**1. Splitting Stockholm files to FASTA per-family:**
Use the C++ utility in `datasets/preparation/RNAFamilies_Stockholm_SeedAlignment_To_PlainFastaFiles/`.
```sh
# Compile the utility
g++ main.cpp -o extract -std=c++14
# Run it
./extract -in path_to/Rfam.seed -out path_to/Rfam_out_files
```

**2. Preparing Clan Datasets:**
The script `datasets/preparation/clans_dataset.py` downloads and prepares the clan-based datasets.

**3. Splitting Train/Test files and getting stats:**
The script `datasets/preparation/constructTrainTestFiles/` splits families into training and testing sets and computes statistics.
```sh
constructTrainTestFiles -in <string> [-out <string> -nf <integer> -mins <integer> -maxs <integer> -pt <integer> -m <string>]
```
- `-in`: Path to the directory of FASTA files.
- `-out`: Output directory for results.
- `-pt`: Percentage of sequences for the test set (default: 30).
- `-m`: Mode. Use `stt` to split all files in the input folder into train/test sets. Other modes (`i`, `s`, `sttmm`, `sttm`) are available for more specific sampling and statistics generation.

**Note on file encoding:** Ensure your FASTA files use UTF-8 encoding to prevent errors from extra spaces when reading lines. See [this Stack Overflow answer](https://stackoverflow.com/a/73952980/3429103) for details.


## 6) Running experiments and reproducing figures

The pipeline script `experminents_results/paper_exeperiments.py` centralizes all experiment logic.

### Common workflow
1. **Compile extractor once:** `compile_code_MotifsExtractionSelection()`
2. **Run C++ extractor with chosen parameters:** `run_cpp_motif_extraction_and_selection(...)`
3. **Run classification and save metrics:** `run_classification_experiment(...)` or `run_classification_experiment_models_choices(...)`
4. **Plot and write CSV summaries** with the provided helper functions.

### Experiment helpers
You can call these functions from the `main()` block of the script or import and use them in an interactive session.
- **Submotif deletion F vs NF:** `deletion_sub_motifs(is_debug_datasets=False)`
    - Generates a CSV and two figures: evolution of data size and processing time.
- **Motif length sweeps:** `run_motif_length_experiments(is_debug_datasets=False)`
    - Produces multiple plots for motif growth, accuracy, and timing.
- **Parameter sweeps:** `run_same_family_threshold_experiments()` and `run_occurrence_variation_experiments()`
    - Plots accuracy and motif counts across parameter grids.
- **Algorithm choices:** `run_algs_choice_experiments()`
    - Compares EXT vs MLP timing and accuracy.

### Debug vs full datasets
- The script uses a global flag `is_debug_datasets_global_var` (default `True`) to run quickly on a small sample.
- For full reproduction, set `is_debug_datasets_global_var=False` near the top of `paper_exeperiments.py` and ensure the full, unzipped datasets are available on disk.


## 7) Standalone classification CLI

You can run classification independently of the main pipeline using the CLI tools under `LSC-ncRNA-our_method/Classification`.

### Model selection via cross-validation
The script `modelstest.py` evaluates multiple classifiers using 10-fold cross-validation.
```sh
python LSC-ncRNA-our_method/Classification/modelstest.py path/to/motifs.csv
```
Tests: `ext`, `knn`, `rdf`, `gnb`, `dt`, `nlp`, `svc`.

### Train and predict with a chosen model
The script `Main.py` trains a model on the full training set and evaluates it on a test set.
```sh
python LSC-ncRNA-our_method/Classification/Main.py <EXT|RDF|NLP|VOT> path/to/motifs.csv path/to/Test
```
- **Model:** `EXT` (ExtraTrees), `RDF` (RandomForest), `NLP` (MLP), or `VOT` (Voting).
- **motifs.csv:** Path to the motif occurrence matrix generated by the C++ program.
- **Test:** Path to the directory containing test set FASTA files.


## 8) Compared methods

- **BLAST:** `compared_methods/blast_classification.py`
    - Creates a BLAST database from the training set, searches test sequences against it, and post-processes the best hits to classify.
    - Example usage is documented in the script's comments.

- **Infernal:** `compared_methods/infernal/`
    - Contains example inputs and result parsers for reproducing Infernal-based baselines.


## 9) Project structure

```
LSC-ncRNA/
├── LSC-ncRNA-our_method/
│   ├── MotifsExtractionSelection/         # C++ motif extraction/selection (build here)
│   └── Classification/                    # Python classifiers (Main.py, modelstest.py, etc.)
├── compared_methods/                      # BLAST and Infernal baselines
├── datasets/
│   ├── data/                              # Pre-prepared datasets (zip/unzip here)
│   ├── information/                       # Dataset stats
│   └── preparation/                       # From-scratch preparation scripts
├── experminents_results/                  # Pipeline, plots, CSVs for paper experiments
└── results/                               # Aggregated CSV outputs
```


## 10) Troubleshooting and tips

- **Memory and time:** Full Rfam runs are resource-intensive. Use the debug mode first, then scale up on a server with large RAM and an SSD.
- **macOS datatable install:** Use the GitHub source install if wheels are unavailable: `pip install git+https://github.com/h2oai/datatable`.
- **xgboost on macOS (optional):** Run `brew install libomp` before `pip install xgboost`.
- **Paths:** The pipeline constructs output CSV names automatically. Results are saved under `results/` and `experminents_results/`.


## 11) License, authors, citation

**License:** See `LICENSE`.

**Authors:** Ibrahim Chegrane, Nabil Bendjaa, Aida Ouangraoua — CoBIUS LAB, Department of Computer Science, Université de Sherbrooke, Canada.

If you use this code or results in your research, please cite our upcoming paper (details will be added when available).

**Contact:** Aida.Ouangraoua@usherbrooke.ca
