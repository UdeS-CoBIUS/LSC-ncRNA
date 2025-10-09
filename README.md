# LSC-ncRNA: Large-Scale Classification of non-coding RNA

LSC-ncRNA is a sequence-based method for classifying ncRNA families. It builds a vector representation of sequences from common motifs (computed/selected by C++), then trains supervised models (Python) to predict families at scale.

This repository bundles the full workflow—from motif discovery to model evaluation—along with a Python pipeline that can reproduce the paper’s experiments end-to-end. Run the bundled “debug” preset locally to validate your environment, then scale to the full datasets on a high-memory machine.

| Component | Path | Purpose |
| --- | --- | --- |
| Pipeline orchestrator | `experminents_results/paper_exeperiments.py` | Compiles the C++ extractor, launches experiments, aggregates results, and produces plots. |
| Motif extractor (C++) | `LSC-ncRNA-our_method/MotifsExtractionSelection/` | Builds motif occurrence matrices from FASTA inputs. |
| Classifiers (Python) | `LSC-ncRNA-our_method/Classification/` | Trains and scores supervised models using the generated motif matrices. |


## Table of contents

- [LSC-ncRNA: Large-Scale Classification of non-coding RNA](#lsc-ncrna-large-scale-classification-of-non-coding-rna)
  - [Table of contents](#table-of-contents)
  - [1. Quickstart: small demo on a laptop](#1-quickstart-small-demo-on-a-laptop)
    - [Prerequisites](#prerequisites)
    - [Steps](#steps)
  - [2. Full reproduction: paper experiments](#2-full-reproduction-paper-experiments)
    - [Key experiments and helper names](#key-experiments-and-helper-names)
    - [Outputs](#outputs)
  - [3. Requirements and installation](#3-requirements-and-installation)
    - [Core stack](#core-stack)
    - [Python dependencies](#python-dependencies)
  - [4. Build the C++ motif extractor](#4-build-the-c-motif-extractor)
    - [Manual Compilation](#manual-compilation)
    - [Command-Line Usage](#command-line-usage)
  - [5. Datasets (prepared and from-scratch options)](#5-datasets-prepared-and-from-scratch-options)
    - [Prepared archives (recommended starting point)](#prepared-archives-recommended-starting-point)
    - [From scratch (optional)](#from-scratch-optional)
  - [6. Running experiments and reproducing figures](#6-running-experiments-and-reproducing-figures)
    - [Common workflow](#common-workflow)
    - [Experiment helpers](#experiment-helpers)
    - [Debug vs full datasets](#debug-vs-full-datasets)
  - [7. Standalone classification CLI](#7-standalone-classification-cli)
    - [Model selection via cross-validation](#model-selection-via-cross-validation)
    - [Train and predict with a chosen model](#train-and-predict-with-a-chosen-model)
  - [8. Compared methods](#8-compared-methods)
  - [9. Project structure](#9-project-structure)
  - [10. Troubleshooting and tips](#10-troubleshooting-and-tips)
  - [11. License, authors, citation](#11-license-authors-citation)


## 1. Quickstart: small demo on a laptop

Run the lightweight “debug” preset to make sure your environment is ready. The pipeline compiles the C++ extractor, generates a small motif matrix, trains classifiers, and writes sample figures under `experminents_results/` and `results/`.

### Prerequisites
- **C++14 toolchain** (e.g., `g++` on macOS/Linux, MSVC/WSL on Windows)
- **Python 3.9+** (3.8 works) with `pip` available
- **Git LFS** if you plan to pull the larger datasets (optional for the debug run)

### Steps
1. From the repository root, create and activate a virtual environment.
    - macOS/Linux:
      ```bash
      python3 -m venv .venv
      source .venv/bin/activate
      ```
    - Windows (PowerShell):
      ```powershell
      py -3 -m venv .venv
      .venv\Scripts\Activate.ps1
      ```

2. Install the Python dependencies (see [Section 3](#3-requirements-and-installation) for platform-specific notes).
    ```bash
    pip install --upgrade pip
    pip install -r requirements.txt
    ```

3. Build the C++ motif extractor, or let the pipeline compile it automatically (details in [Section 4](#4-build-the-c-motif-extractor)).

4. Run the pipeline in debug mode (ensure your working directory is the repository root).
    ```bash
    python experminents_results/paper_exeperiments.py
    ```

5. Inspect the outputs.
    - Figures such as `experminents_results/accuracy_fixed_len_EXT.png`
    - CSV summaries in `experminents_results/results/`
    A successful run indicates that both the C++ and Python components are wired correctly.


## 2. Full reproduction: paper experiments

Running the full Rfam-scale experiments requires a machine with substantial memory (tens of GB RAM) and fast storage. All experiment entry points live in `experminents_results/paper_exeperiments.py`. Enable experiments by editing the `main()` function (uncomment the calls you need) or by importing the module and calling the helpers from a Python session:

```python
from experminents_results.paper_exeperiments import (
    run_same_family_threshold_experiments,
    run_occurrence_variation_experiments,
)

run_same_family_threshold_experiments(
    same_family_pct_values=[30, 40, 50],
    min_occurrence_counts=[1, 2, 3],
)
```

### Key experiments and helper names
- **Submotif deletion impact**: `deletion_sub_motifs()`
- **Motif length sweeps (fixed & combined)**: `run_motif_length_experiments()`
- **Conservation vs occurrence sweeps**:
  - `run_same_family_threshold_experiments()` (legacy name `run_beta_experiments`)
  - `run_occurrence_variation_experiments()` (legacy name `run_alpha_variance_experiments`)
  - Each helper takes the “new” descriptive parameters (`sameFamilyPct`, `occVarTol`, `min_occurrence_count`) while still accepting the legacy keywords for backwards compatibility.
- **Algorithm choices & timings**: `run_algs_choice_experiments()`

### Outputs
- CSV summaries in `experminents_results/results/`
- Figures in `experminents_results/` (e.g., `Growth_NB_motifs.png`, `Accuracy_combined_len_EXT_MLP.png`, `beta_gamma_accuracy_plot.png`)
- Logs written to the console and to `experminents_results/logs/` (created on demand)


## 3. Requirements and installation

### Core stack
- **C++14 compiler** (e.g., `g++` on macOS/Linux, MSVC/WSL on Windows)
- **Python 3.8+** (Python 3.9+ recommended) with `pip`
- **Git LFS** for downloading the larger dataset archives (optional for debug mode)

### Python dependencies
Install the baseline dependencies via `requirements.txt`:

```bash
pip install -r requirements.txt
```

The file lists:
- `biopython`, `datatable`, `scikit-learn`, `pandas`, `numpy`, `matplotlib`
- `portalocker`, `pyahocorasick`, `suffix-trees`
- Optional extras (commented) such as `xgboost`

> **macOS tip**: If the wheel for `datatable` is unavailable, install from source:
> ```bash
> pip install git+https://github.com/h2oai/datatable
> ```
>
> **Windows tip**: Use WSL2 or ensure the MSVC build tools are available to compile the C++ binaries.

After installation, you can validate the environment with:

```bash
python -m pip check
```

This command should report no dependency conflicts before you proceed to the full experiments.


## 4. Build the C++ motif extractor

The C++ program computes and selects sequence motifs, creating a vectorial representation of ncRNA sequences. The result is a 2D matrix (saved as a CSV file) where each row is a sequence and each column is the occurrence count of a selected motif.

The Python pipeline can compile the C++ source for you, or you can build it manually.

### Manual Compilation
```bash
cd LSC-ncRNA-our_method/MotifsExtractionSelection
g++ -std=c++14 
    -o MotifsExtractionSelection 
    Main.cpp SuffixTree_QuadraticTime.cpp 
    FastaFilesReader.cpp CommonMotifs.cpp SequenceIdManager.cpp
```
This produces the executable `LSC-ncRNA-our_method/MotifsExtractionSelection/MotifsExtractionSelection`, which is used by the Python pipeline.

### Command-Line Usage
The program can be used like this:
```bash
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
```bash
nohup ./MotifsExtractionSelection 
    -in "/data/ibra/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_600_Train_Test/Train" 
    -minl 2 -maxl 8 -b 50 -g 1 -tn F_600 > out_F_600 &
```

The output CSV file name is constructed from the parameters, for example: `del_No_nbF_F_600_min_2_max_8_sameFamilyPct_50_occVarTol_-1_nbOccrs_1.csv`.


## 5. Datasets (prepared and from-scratch options)

### Prepared archives (recommended starting point)
Unzip the dataset archives in `datasets/data/` before running full experiments.

| Dataset | Archive | Description |
| --- | --- | --- |
| Rfam 14.1 families | [`datasets/data/Rfam_14.1_dataset.zip`](datasets/data/Rfam_14.1_dataset.zip) | Default train/test splits used throughout the paper. |
| Deep ncRNA benchmark | [`datasets/data/deep_ncrna_datasets.zip`](datasets/data/deep_ncrna_datasets.zip) | 88-family dataset with noise variants from the DeepNC RNA paper. |
| Rfam 14.8 clans | [`datasets/data/Clans_ncRNA_from_Rfam_14.8.zip`](datasets/data/Clans_ncRNA_from_Rfam_14.8.zip) | Clan-level and family-level splits for low/high similarity scenarios. |

After extraction, you should see ready-to-use `Train/` and `Test/` folders matching the experiment scripts.

### From scratch (optional)
The `datasets/preparation/` folder contains scripts to recreate the datasets:

1. **Split Stockholm seeds to per-family FASTA** — `RNAFamilies_Stockholm_SeedAlignment_To_PlainFastaFiles/`
    ```bash
    g++ main.cpp -o extract -std=c++14
    ./extract -in path/to/Rfam.seed -out path/to/Rfam_out_files
    ```

2. **Build clan datasets** — `clans_dataset.py` downloads and assembles the clan-level corpora.

3. **Construct train/test splits & summary statistics** — `constructTrainTestFiles/`
    ```bash
    constructTrainTestFiles -in path/to/fasta_dir -out path/to/output -pt 30 -m stt
    ```
    - `-pt` controls the percentage of sequences assigned to the test split (default 30%).
    - Alternate modes (`i`, `s`, `sttmm`, `sttm`) support sampling and statistics-only workflows.

> **Encoding reminder**: Keep FASTA files in UTF-8 to avoid trailing-space issues when parsing lines. See [this Stack Overflow answer](https://stackoverflow.com/a/73952980/3429103) for the background.


## 6. Running experiments and reproducing figures

All experiment orchestration lives in `experminents_results/paper_exeperiments.py`. The pipeline follows the sequence below:

```
compile_code_MotifsExtractionSelection()
        ↓
run_cpp_motif_extraction_and_selection(...)
        ↓
run_classification_experiment(...)
        ↓
plot_*(...)
```

### Common workflow
1. **Compile the extractor** once per session using `compile_code_MotifsExtractionSelection()`.
2. **Generate motif matrices** with `run_cpp_motif_extraction_and_selection(...)`.
3. **Train and evaluate models** by calling `run_classification_experiment(...)` or `run_classification_experiment_models_choices(...)`.
4. **Visualise and export results** via the plotting helpers (`plot_accuracy_*`, `plot_motifs_*`, etc.).

### Experiment helpers
Use these from the `main()` function or from an interactive Python session.
- **Submotif deletion F vs NF** — `deletion_sub_motifs(is_debug_datasets=False)` generates CSVs and plots for data size vs processing time.
- **Motif length sweeps** — `run_motif_length_experiments(is_debug_datasets=False)` produces growth, accuracy, and timing figures.
- **Parameter sweeps** — `run_same_family_threshold_experiments()` and `run_occurrence_variation_experiments()` chart accuracy and motif counts across grids.
- **Algorithm choices** — `run_algs_choice_experiments()` compares EXT vs MLP timing and accuracy.

### Debug vs full datasets
The module-level variable `is_debug_datasets_global_var` (near the top of `paper_exeperiments.py`) controls dataset size.

```python
is_debug_datasets_global_var = True  # change to False for full Rfam runs
```

- `True` (default) runs a curated mini sample suitable for laptops.
- `False` switches to the full datasets—ensure the archives from [Section 5](#5-datasets-prepared-and-from-scratch-options) are extracted and that your machine has sufficient RAM.


## 7. Standalone classification CLI

Run classifiers independently of the main pipeline via the scripts in `LSC-ncRNA-our_method/Classification/`.

### Model selection via cross-validation
`modelstest.py` scores several algorithms with 10-fold CV (70/30 split) and writes metrics to stdout.

```bash
python LSC-ncRNA-our_method/Classification/modelstest.py results/del_No_nbF_F_600_min_2_max_8_sameFamilyPct_50_occVarTol_-1_nbOccrs_1.csv
```

Algorithms tested: `ext`, `knn`, `rdf`, `gnb`, `dt`, `nlp`, `svc`. The script reports accuracy, precision, recall, F1, and execution time for each model.

### Train and predict with a chosen model
`Main.py` trains the selected classifier on the full motif matrix and evaluates it against a FASTA test folder.

```bash
python LSC-ncRNA-our_method/Classification/Main.py EXT \
    results/del_No_nbF_F_600_min_2_max_8_sameFamilyPct_50_occVarTol_-1_nbOccrs_1.csv \
    datasets/data/Rfam_14.1_dataset/Rfam_600_Train_Test/Test
```

- **Model (arg 1)** — `EXT`, `RDF`, `NLP`, or `VOT`.
- **Motif matrix (arg 2)** — CSV produced by the C++ extractor.
- **Test folder (arg 3)** — Directory containing FASTA files to classify.

Outputs are streamed to stdout and saved to `Classification/logs/` (created if absent). You can use shell redirection to persist the metrics, e.g. `... > logs/ext_run.txt`.


## 8. Compared methods

- **BLAST** — `compared_methods/blast_classification.py`
    - Builds a BLAST database from the training FASTA files, queries the test set, and post-processes hits to produce predictions.
    - Example run (after installing BLAST+):
        ```bash
        python compared_methods/blast_classification.py \
                datasets/data/Rfam_14.1_dataset/Rfam_600_Train_Test/Train \
                datasets/data/Rfam_14.1_dataset/Rfam_600_Train_Test/Test \
                Rfam_600_demo > results/blast_Rfam_600.log
        ```
    - Outputs classification metrics and timing to stdout (captured above) and writes intermediate BLAST files under `compared_methods/blast_tmp/`.

- **Infernal** — `compared_methods/infernal/`
    - Provides example configs, test inputs, and parsers to reproduce Infernal baselines.
    - Follow the README within the folder to run `cmbuild`, `cmcalibrate`, and `cmsearch`, then parse outputs using the supplied Python scripts in `ParseResultInfernal/`.


## 9. Project structure

```
LSC-ncRNA/
├── LSC-ncRNA-our_method/               # Core method implementation
│   ├── MotifsExtractionSelection/      # C++ sources + build artifacts for motif extraction
│   └── Classification/                 # Python classifiers (Main.py, modelstest.py, utilities)
├── compared_methods/                   # BLAST and Infernal baselines + helpers
├── datasets/
│   ├── data/                           # Prepared archives (see Section 5)
│   ├── information/                    # Dataset statistics spreadsheets and CSVs
│   └── preparation/                    # Scripts to recreate datasets from raw sources
├── experminents_results/               # Experiment pipeline, plots, CSVs, logs
└── results/                            # Aggregated outputs copied from experiment runs
```


## 10. Troubleshooting and tips

- **Memory and time:** Full Rfam runs are resource-intensive. Use the debug mode first, then scale up on a server with large RAM and an SSD.
- **macOS datatable install:** Use the GitHub source install if wheels are unavailable: `pip install git+https://github.com/h2oai/datatable`.
- **xgboost on macOS (optional):** Run `brew install libomp` before `pip install xgboost`.
- **Paths:** The pipeline constructs output CSV names automatically. Results are saved under `results/` and `experminents_results/`.


## 11. License, authors, citation

**License:** See `LICENSE`.

**Authors:** Ibrahim Chegrane, Nabil Bendjaa, Aida Ouangraoua — CoBIUS LAB, Department of Computer Science, Université de Sherbrooke, Canada.

If you use this code or results in your research, please cite our upcoming paper (details will be added when available).

**Contact:** Aida.Ouangraoua@usherbrooke.ca
