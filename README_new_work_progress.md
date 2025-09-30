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

Prerequisites
- C++14 toolchain (g++)
- Python 3.9+ recommended (3.8+ works); pip

Steps
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

Key experiments and where they’re implemented
- Submotif deletion impact: `deletion_sub_motifs()`
- Number of motifs vs accuracy (fixed- and combined-length): `run_motif_length_experiments()`
- Conservation (alpha/beta) and occurrence (gamma) filtering sweeps:
    - `run_beta_experiments()` (alpha in paper terminology)
    - `run_alpha_variance_experiments()` (beta in paper terminology)
- Algorithm choices and timings: `run_algs_choice_experiments()`

Outputs
- CSV summaries in `results/` or `experminents_results/results/`
- Figures in `experminents_results/` (for example: `Growth_NB_motifs.png`, `Accuracy_combined_len_EXT_MLP.png`, `beta_gamma_accuracy_plot.png`)


## 3) Requirements and installation

Core requirements
- C++14 compiler (g++)
- Python 3.8+ (3.9+ recommended)

Python packages (install via pip)
- biopython, datatable, scikit-learn, pandas, numpy, matplotlib, portalocker, pyahocorasick, suffix-trees
- Optional: xgboost (not required for the paper’s results). On macOS you may need `brew install libomp` first.

Platform notes
- macOS: if `pip install datatable` fails, install from source:
    ```sh
    pip install git+https://github.com/h2oai/datatable
    ```
- Windows: we recommend using WSL2 or a recent MSVC build tools setup for C++14.


## 4) Build the C++ motif extractor

The Python pipeline compiles exactly these sources:
`LSC-ncRNA-our_method/MotifsExtractionSelection/{Main.cpp, SuffixTree_QuadraticTime.cpp, FastaFilesReader.cpp, CommonMotifs.cpp, SequenceIdManager.cpp}`

You can either let the pipeline compile it for you (it calls `compile_code_MotifsExtractionSelection()`), or build it manually:

```sh
cd LSC-ncRNA-our_method/MotifsExtractionSelection
g++ -std=c++14 \
    -o MotifsExtractionSelection \
    Main.cpp SuffixTree_QuadraticTime.cpp \
    FastaFilesReader.cpp CommonMotifs.cpp SequenceIdManager.cpp
```

This produces the executable `LSC-ncRNA-our_method/MotifsExtractionSelection/MotifsExtractionSelection` used by the pipeline.


## 5) Datasets

Prepared datasets (recommended)
- The repository includes pre-prepared datasets under `datasets/data/`. Unzip what you need (Rfam 14.1, deep ncRNA, Rfam 14.8 clans). Most paper experiments use the Rfam 14.1 set with predefined train/test splits.

From scratch (optional)
- Seed splitting to FASTA per-family: `datasets/preparation/RNAFamilies_Stockholm_SeedAlignment_To_PlainFastaFiles/`
- Clan preparation script: `datasets/preparation/clans_dataset.py`
- Train/test splitting and dataset stats: `datasets/preparation/constructTrainTestFiles/`

Train/test splitting modes are documented in the preparation scripts. Beware of file encoding (use UTF-8) to avoid extra spaces when reading lines (see https://stackoverflow.com/a/73952980/3429103).


## 6) Running experiments and reproducing figures

The pipeline script centralizes experiment logic: `experminents_results/paper_exeperiments.py`.

Common steps
1) Compile extractor once: `compile_code_MotifsExtractionSelection()`
2) Run C++ extractor with chosen parameters: `run_cpp_motif_extraction_and_selection(...)`
3) Run classification for EXT/MLP and save metrics: `run_classification_experiment(...)` or `run_classification_experiment_models_choices(...)`
4) Plot and write CSV summaries with the provided helpers

Experiment helpers (call from `main()` or import and use)
- Submotif deletion F vs NF: `deletion_sub_motifs(is_debug_datasets=False)`
    - Generates CSV and two figures: evolution of data size and processing time.
- Motif length sweeps: `run_motif_length_experiments(is_debug_datasets=False)`
    - Produces: `Growth_NB_motifs.png`, `Growth_NB_motifs_combined.png`, `accuracy_fixed_len_EXT.png`, `Accuracy_combined_len_EXT_MLP.png`, `Total_time_combined_len_EXT_MLP.png`.
- Alpha/Beta/Gamma sweeps: `run_beta_experiments()` and `run_alpha_variance_experiments()`
    - Plots accuracy and motif counts across parameter grids.
- Algorithm choices: `run_algs_choice_experiments()`
    - Compares EXT vs MLP timing and accuracy for selected (min,max) and gamma.

Debug vs full datasets
- The script has `is_debug_datasets_global_var` (default True) and a small preset of sizes to run quickly on a laptop.
- For full reproduction, set `is_debug_datasets_global_var=False` near the top of `paper_exeperiments.py`, and ensure the unzipped datasets are available on disk.


## 7) Standalone classification CLI

You can run classification independently of the pipeline using the CLI under `LSC-ncRNA-our_method/Classification`.

- Model selection via cross-validation
    ```sh
    python LSC-ncRNA-our_method/Classification/modelstest.py path/to/motifs.csv
    ```
    Tests: ext, knn, rdf, gnb, dt, nlp, svc using 10-fold CV (70/30 split).

- Train and predict with a chosen model
    ```sh
    python LSC-ncRNA-our_method/Classification/Main.py <EXT|RDF|NLP|VOT> path/to/motifs.csv path/to/Test
    ```


## 8) Compared methods

- BLAST: `compared_methods/blast_classification.py`
    - Creates a DB from Train, searches Test, then post-processes the best hits to produce classification and timing.
    - Example usage is embedded in the script and the historical README (see comments in the file).

- Infernal: `compared_methods/infernal/`
    - Contains example inputs and parsers for reproducing Infernal-based baselines.


## 9) Project structure (short)

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

- Memory and time: full Rfam runs are heavy. Use the debug mode first, then scale up on a server with large RAM and SSD.
- macOS datatable install: use the GitHub source install if wheels are unavailable: `pip install git+https://github.com/h2oai/datatable`.
- xgboost on macOS (optional): `brew install libomp` before `pip install xgboost`.
- Encoding: ensure input FASTA files are UTF-8 encoded.
- Paths: the pipeline constructs output CSV names automatically; results are under `results/` and `experminents_results/`.


## 11) License, authors, citation

License: see `LICENSE`.

Authors: Ibrahim Chegrane, Nabil Bendjaa, Aida Ouangraoua — CoBIUS LAB, Department of Computer Science, Université de Sherbrooke, Canada.

If you use this code or results in your research, please cite our upcoming paper (details will be added when available).

Contact: Aida.Ouangraoua@usherbrooke.ca
