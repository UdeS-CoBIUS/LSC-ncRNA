## Experiments workspace

This directory groups everything needed to reproduce the paper results and to keep reference outputs under version control. It is the main entry point to re-run experiments.

### Directory layout

- `scripts/` – runnable Python helpers such as `paper_exeperiments.py`, lightweight utilities, and any future notebooks. These consume datasets from the project root and write new material into `outputs/` by default.
- `outputs/` – scratch space for regenerated results (plots, CSVs, logs). It is intentionally empty and should remain git-ignored apart from placeholder files. The automation currently creates two top-level folders:
	- `debug_small_dataset/` – results from fast sanity runs that use the tiny debug datasets (`is_debug_datasets_global_var = True`).
	- `full_paper_dataset/` – results from the full-size configuration that mirrors the published experiments (`is_debug_datasets_global_var = False`).
	Each folder may contain subdirectories such as `plots/`, `tables/`, or `logs/`, created automatically by `paper_exeperiments.py`.
- `artifacts_paper_results/` – read-only copies of the official paper artifacts produced on the 48‑CPU/600 GB server. Do **not** modify or overwrite these files; instead regenerate into `outputs/` and compare.

### Provenance

The contents of `artifacts_paper_results/` were generated on a high-memory (≈600 GB RAM) machine with 48 CPU cores. Keep the folder immutable so the archived numbers continue to match the published paper.
