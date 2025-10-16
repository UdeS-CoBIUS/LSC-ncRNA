## Experiments workspace

This directory groups everything needed to reproduce the paper results and to keep reference outputs under version control. It is the main entry point to re-run experiments.

### Directory layout

- `scripts/` – runnable Python helpers such as `paper_exeperiments.py`, lightweight utilities, and any future notebooks. These consume datasets from the project root and write new material into `outputs/` by default.
- `outputs/` – scratch space for regenerated results (plots, CSVs, logs). It is intentionally empty and should remain git-ignored apart from placeholder files. The automation currently creates two top-level folders:
	- `debug_small_dataset/` – results from fast sanity runs that use the tiny debug datasets (`is_debug_datasets_global_var = True`).
	- `full_paper_dataset/` – results from the full-size configuration that mirrors the published experiments (`is_debug_datasets_global_var = False`).
	Each folder may contain subdirectories such as `plots/`, `tables/`, or `logs/`, created automatically by `paper_exeperiments.py`.
- `artifacts_paper_results/` – read-only copies of the official paper artifacts produced on the 48‑CPU/600 GB server. Do **not** modify or overwrite these files; instead regenerate into `outputs/` and compare.

### Running experiments

The main entry point is `scripts/paper_exeperiments.py`, which now exposes a lightweight command-line interface so you can queue exactly the steps you want without commenting code in the file.

- Show the available steps in their execution order:

	```bash
	python experiments/scripts/paper_exeperiments.py --list
	```

- Run a single step (for example only the motif-length sweep) while staying in debug mode:

	```bash
	python experiments/scripts/paper_exeperiments.py --experiments motif-length --mode debug
	```

- Chain multiple steps explicitly (order follows the command line):

	```bash
	python experiments/scripts/paper_exeperiments.py --experiments compile-cpp deletion-sub-motifs motif-length
	```

- Skip any costly stage while keeping the rest (here we keep the default experiments but avoid recompiling C++):

	```bash
	python experiments/scripts/paper_exeperiments.py --skip compile-cpp
	```

- Run everything in the curated order:

	```bash
	python experiments/scripts/paper_exeperiments.py --experiments all --mode full
	```

By default the script still compiles the C++ binary and launches the deletion/no-deletion experiment, matching the previous behaviour. Use `--mode full` to switch from the default debug datasets to the full paper configuration.

### Provenance

The contents of `artifacts_paper_results/` were generated on a high-memory (≈600 GB RAM) machine with 48 CPU cores. Keep the folder immutable so the archived numbers continue to match the published paper.
