## Experiments workspace

This directory groups everything needed to reproduce the paper results and to keep reference outputs under version control. It is the main entry point to re-run experiments.

### Directory layout

- `scripts/` – runnable Python helpers such as `paper_exeperiments.py`, lightweight utilities, and any future notebooks. These consume datasets from the project root and write new material into `outputs/` by default.
- `outputs/` – scratch space for regenerated results (plots, CSVs, logs). It is intentionally empty and should remain git-ignored apart from placeholder files. Feel free to create subfolders like `reproducible/` or `debug/` depending on the hardware profile you used for a run.
- `artifacts paper results/` – read-only copies of the official paper artifacts produced on the 48‑CPU/600 GB server. Do **not** modify or overwrite these files; instead regenerate into `outputs/` and compare.

### Provenance

The contents of `artifacts paper results/` were generated on a high-memory (≈600 GB RAM) machine with 48 CPU cores. Keep the folder immutable so the archived numbers continue to match the published paper.
