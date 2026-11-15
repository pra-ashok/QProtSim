# First analysis

## Install and check the command

Use Python 3.11 or later in an isolated environment. Start with the [CPU installation](installation.md#cpu), then check the CLI.

```bash
qprotsim --help
qprotsim analyze --help
```

`QProtSim`, `QPS`, and `qps` are additional command names. The importable Python package is `qprotsim`.

## Create a configuration

Save the following as `config.yaml` in your working directory. The 20-residue Trp-cage structure is a practical starting point for a statevector simulator.

```yaml
pdb_code: 1L2Y
molecule_name: Trp Cage
download_dir: ./pdb_cache
chain: A
protein_only: true
shots: 10000
pl_device: default.qubit
bins: 100
output_dir: ./outputs/1l2y
report_file_name: trp_cage_report.pdf
```

The first run downloads the PDB file. For a local structure, replace `pdb_code` and `download_dir` with an existing file path.

```yaml
pdb_path: ./input/1L2Y.pdb
```

Relative PDB and output paths are interpreted from the shell's working directory. Explicitly select `chain` for reproducible residue selection. The [configuration reference](configuration.md) covers defaults, path precedence, and batch settings.

## Run the report workflow

```bash
qprotsim analyze --config config.yaml
```

The configured directory receives `trp_cage_report.pdf`, CSV tables, sampled states, and figures. The PDF path is printed at completion. Figure numbering depends on which optional structure and motif plots succeed. See the [output reference](outputs.md) for exact naming and interpretation.

To select a different output directory or enable verbose debugging, use the CLI options.

```bash
qprotsim analyze -c config.yaml --output-dir ./outputs/repeat
qprotsim analyze -c config.yaml --debug
```

Use a fresh directory for each analysis. Repeated runs reuse artifact filenames, including when `run_stamp` changes.

## Run from Python

```python
from pathlib import Path

from qprotsim.report_analysis import run_analysis_from_yaml

report_path = run_analysis_from_yaml(Path("config.yaml"))
print(report_path)
```

For direct access to data frames, state samples, and matrices, call `run_analysis`. Its return value is an `AnalysisOutputs` object. The [Python API](api.md#run_analysis) includes a complete example that generates the PDF separately.

## Plan memory and sampling

A dense complex128 statevector alone uses approximately $16\times2^N$ bytes for $N$ residues. At 20 residues that is 16 MiB, at 25 it is 512 MiB, and at 30 it is 16 GiB. Simulator workspaces, sample arrays, plots, and pairwise analysis add to this cost. A GPU also needs enough device memory for the statevector.

Increasing `shots` improves sampling precision while increasing sample storage and analysis time. Rare conditional events can remain poorly sampled. Compare repeated runs and inspect the [metric limitations](theory.md#assumptions-and-interpretation) before interpreting small differences.
