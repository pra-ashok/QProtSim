# Configuration

A configuration file tells QProtSim which structure to analyse, how many samples to take, and where to save the results. Save it as a plain-text `.yaml` file, then run:

```bash
qprotsim analyze --config config.yaml
```

If you are new to QProtSim, start with the complete example in [your first analysis](getting-started.md).

## A single protein

```yaml
pdb_code: 1L2Y
molecule_name: Trp Cage
download_dir: ./pdb_cache
chain: A
protein_only: true
shots: 10000
pl_device: default.qubit
bins: 100
output_dir: ./outputs/trp_cage
report_file_name: report.pdf
```

You must provide a `molecule_name` and either a `pdb_code` to download or a `pdb_path` pointing to an existing file. To use a local file, replace `pdb_code` and `download_dir` with:

```yaml
pdb_path: ./input/your_protein.pdb
```

If both input settings are present, QProtSim uses `pdb_path`.

### Structure and sampling

Defaults apply when a setting is omitted. The first-run example uses fewer samples and plot bins than the defaults to keep the initial analysis manageable.

| Setting | Default | Meaning |
| --- | --- | --- |
| `pdb_path` | No default | Path to a local PDB structure file; required unless you supply `pdb_code`. |
| `pdb_code` | No default | Four-character Protein Data Bank identifier; required unless you supply `pdb_path`. |
| `molecule_name` | Required | Name used in the report and output filenames. Use a simple, descriptive name. |
| `chain` | All chains for a single analysis | Chain identifier such as `A`. Use `null` for all chains. An empty string or `"--"` also selects all chains. |
| `protein_only` | `true` | Restricts the analysis to recognised protein residues. Keep this enabled for ordinary use. |
| `shots` | `100000` | Number of samples drawn from the circuit. More samples improve precision but need more time and storage. |
| `backend` | `pennylane` | Circuit framework: `pennylane` or `qiskit`. |
| `pl_device` | `default.qubit` | PennyLane device, used when `backend: pennylane`. See [installation](installation.md). |
| `qiskit_backend` | `aer_simulator` | Local Aer backend, used when `backend: qiskit`. Requires the `qiskit` installation extra. |
| `backend_options` | `{}` | Mapping of simulator options. See the Qiskit example below. |
| `seed` | `null` | Optional nonnegative integer for reproducible sampling. |
| `bins` | `1000` | Number of intervals in the energy plots. This does not change the binning used for residue-sensitivity scores. |

Use positive whole numbers for `shots` and `bins`. Write `true` and `false` without quotation marks. A setting such as `protein_only: "false"` is text rather than a YAML boolean and can give unexpected results.

Set `chain` explicitly when comparing analyses. Selecting all chains can greatly increase memory use. Setting `protein_only: false` may include residues for which the model has no reference values, causing an error.

### Select Qiskit

Install with `python -m pip install ".[qiskit]"`, then add these settings to a single-protein or batch configuration:

```yaml
backend: qiskit
qiskit_backend: aer_simulator
backend_options:
  method: statevector
  device: CPU
seed: 42
```

`backend_options` must be a YAML mapping. Aer also supports other simulation methods, such as `matrix_product_state`; valid options depend on the selected backend and installed hardware. See the [Aer simulator options](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html). Invalid backend names or options stop the run with an error. Use `pl_device` to choose a PennyLane device when `backend: pennylane`.

Both frameworks return the same sample convention: computational state `|0>` is the ground model state with Pauli-Z value `+1`; `|1>` is the excited model state with Pauli-Z value `-1`. Qiskit measurements use binary values `b = 0, 1`, which QProtSim converts using `z = 1 - 2*b`. In a displayed Qiskit bit string, the rightmost bit belongs to qubit 0; exported sample columns follow residue order, starting with qubit 0. See [IBM's bit-ordering guide](https://quantum.cloud.ibm.com/docs/en/qiskit/bit-ordering). The kets `|+>` and `|->` describe Pauli-X eigenstates and are distinct from the numerical Pauli-Z outcomes `+1` and `-1`.

### Output files

| Setting | Default | Meaning |
| --- | --- | --- |
| `output_dir` | `qprotsim_outputs` beside the configuration file | Folder for the PDF, tables, and figures. |
| `report_file_name` | `report.pdf` | PDF filename inside that folder. Use a filename such as `trp_cage_report.pdf`, without a folder path. |

Choose a separate output folder for each analysis you want to keep. Reusing a folder can overwrite previous results.

### Downloads

These settings apply when you use `pdb_code`:

| Setting | Default | Meaning |
| --- | --- | --- |
| `download_dir` | The folder where you run the command | Folder for downloaded PDB files. Existing files can be reused. |
| `force_download` | `false` | Download the structure again even when a saved copy exists. |
| `verify_checksum` | `true` | Check downloaded file integrity using remote checksum information when available. |

### Paths and command-line overrides

An explicit relative path, such as `./input/protein.pdb`, starts from the folder where you run the command. It does not start from the configuration file's folder. If you omit `output_dir`, the default output folder is placed beside the configuration file.

You can change the output location without editing the YAML:

```bash
qprotsim analyze -c config.yaml --output-dir ./results/trp_cage
```

To choose a PDF path for a single analysis:

```bash
qprotsim analyze -c config.yaml --report ./results/trp_cage/custom_report.pdf
```

`--output-dir` takes priority over the YAML output directory. If you supply `--report` alone, its parent folder becomes the output directory when the path includes a folder. If you supply both options, the PDF is placed in `--output-dir` using the filename from `--report`.

### Supported settings

Use the settings listed on this page. Adding keys such as `temperature`, `beta`, `gamma`, or `weights` does not change the standard YAML analysis. The [science guide](theory.md) explains the model's fixed assumptions, and the [Python guide](api.md) covers direct numerical work.

## Several proteins in one run

Batch mode downloads and analyses a list of PDB entries using shared settings. Start with small proteins; a batch does not reduce the memory needed for each individual analysis.

Save this example as `batch.yaml`:

```yaml
batch_mode: true
pdb_codes:
  - 1L2Y
  - 1UAO
molecule_names:
  1L2Y: Trp Cage
  1UAO: Chignolin
download_dir: ./pdb_cache
chain: A
protein_only: true
shots: 10000
pl_device: default.qubit
bins: 100
output_dir: ./outputs/batch
report_file_name: "{pdb_code}_report.pdf"
```

Run it with:

```bash
qprotsim analyze -c batch.yaml
```

Each protein gets a separate output folder, such as `outputs/batch/1l2y/1l2y_report.pdf`.

| Batch setting | Default | Meaning |
| --- | --- | --- |
| `batch_mode` | `false` | Set to `true` to analyse a list of entries. |
| `pdb_codes` | Required for a batch | List of PDB identifiers to download. |
| `molecule_names` | PDB code used as the name | Optional display names, keyed by uppercase PDB code. |
| `chain` | `A` | Shared chain selection. Use `null` for all chains. This default differs from a single analysis. |
| `output_dir` | `qprotsim_batch_outputs` beside the configuration file | Parent folder for the individual result folders. |
| `report_file_name` | `"{pdb_code}_report.pdf"` | Filename pattern used for each report. |

The filename pattern accepts `{pdb_code}` and `{molecule_name}`. Both substitutions use lowercase text; spaces in molecule names become underscores. Other shared settings use the single-analysis defaults above.

Use `--output-dir` to change the parent output folder. For batch reports, set the filename pattern in YAML; `--report` applies only to a single analysis. Batch mode accepts PDB codes rather than a list of local file paths.

All requested downloads must succeed before analysis begins. If one protein fails during analysis, QProtSim continues with the others. Check the final success count and confirm that each expected PDF exists.
