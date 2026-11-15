# Troubleshooting

## Installation cannot find a compatible wheel

Use Python 3.11 or later with a recent pip. Binary availability varies across platforms and GPU plugins. A Python 3.11 or 3.12 environment is a practical choice when an upstream package has limited wheels.

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install ".[cpu]"
```

FreeSASA may build from source and requires a C compiler and Python development headers. Follow its upstream installation instructions for your operating system. [FreeSASA Python installation](https://freesasa.github.io/python/intro.html#installation)

## The requested accelerator is unavailable

Install the matching extra and check device construction and execution directly using the test in the [installation guide](installation.md#verify-the-backend). Device creation in QProtSim can silently fall back to `default.qubit`. A `lightning.gpu` execution error containing GPU-related terms triggers an additional CPU retry. The report records the requested device name, so independently verify the backend before benchmarking.

For CPU execution, set the device explicitly.

```yaml
pl_device: default.qubit
```

## No residues match the selection

Check the PDB chain identifier and use standard amino-acid residues with `protein_only: true`. For all chains, set `chain: null`. Multi-chain structures and sequence gaps require careful review of the coupling model's adjacency assumptions.

Nonstandard residues can lack entries in the energy, hydrophobicity, or maximum-SASA tables. `protein_only: false` retains additional residues during selection, and those residues can raise a lookup error downstream.

## PDB download or structure images fail

Check the four-character identifier, network connectivity, and download-directory permissions. Use `pdb_path` to analyse a cached local structure. Download validation compares remote checksum metadata when available and also checks basic file integrity.

Structure images use RCSB downloads, optional PyMOL, and a Matplotlib fallback. Missing PyMOL can produce a warning while analysis continues. Install `.[visualization]` when a compatible PyMOL distribution is available.

## Information metrics contain NaN values

The residue-information estimator can encounter empty conditional histograms when a residue is observed in only one state. Column normalization can also divide by zero when all divergences vanish. Inspect state counts and repeat the analysis with an appropriate shot count. Finite sampling can leave very rare states unobserved.

The information functions infer ±1 encoding from the presence of a negative value. An all-`+1` array is ambiguous and can be interpreted as binary excited states. Convert known Pauli-Z samples explicitly before direct information-function calls.

```python
states = (z_samples == -1).astype(int)
```

See the [API reference](api.md) for each function's accepted encoding and the [scientific model](theory.md) for the estimators.

## Basin command rejects exported samples

The `basin` command expects a CSV whose first column is named `shot`, followed by residue columns containing ±1 states. The report's `<molecule>_z_samples.csv` uses this format. Custom `--E0` and `--E1` inputs are NumPy `.npy` arrays and must be supplied together in the same residue order.
