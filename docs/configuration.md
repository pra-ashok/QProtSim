# Configuration

The `qprotsim analyze` command reads a YAML mapping. Single analyses require `molecule_name` and either `pdb_path` or `pdb_code`.

## A single protein

```yaml
pdb_code: 1L2Y
molecule_name: Trp Cage
download_dir: ./input
chain: A
protein_only: true
shots: 10000
pl_device: default.qubit
bins: 100
output_dir: ./outputs/trp_cage
report_file_name: report.pdf
```

```bash
qprotsim analyze -c config.yaml
```

To use a local structure, replace `pdb_code` and the download settings with `pdb_path: ./input/1L2Y.pdb`. If both input keys appear, `pdb_path` takes precedence.

| Key | Default | Meaning |
| --- | --- | --- |
| `pdb_path` | Required unless `pdb_code` is supplied | Local PDB file. |
| `pdb_code` | Required unless `pdb_path` is supplied | Four-character RCSB identifier to download. |
| `molecule_name` | Required | Report title and CSV filename prefix. |
| `chain` | `null` for a single run | Chain ID. `null`, an empty string, or `"--"` selects all chains. |
| `protein_only` | `true` | Select protein residues recognized by the package. |
| `shots` | `100000` | Number of sampled circuit configurations. |
| `pl_device` | `default.qubit` | PennyLane device name. See [installation](installation.md). |
| `bins` | `1000` | Number of bins for energy plots. |
| `output_dir` | `<config_dir>/qprotsim_outputs` | Destination for figures, CSV files, and the PDF. |
| `report_file_name` | `report.pdf` | PDF filename inside the output directory. |
| `download_dir` | Current working directory | Cache for downloaded PDB files. |
| `force_download` | `false` | Fetch a fresh PDB file even if a cached file exists. |
| `verify_checksum` | `true` | Check file integrity when checksum metadata is available. |

Use YAML booleans `true` and `false` without quotes. Choose positive integer values for `shots` and `bins`. Residues must have reference energies supported by the package. Setting `protein_only: false` includes additional residues in selection, which can fail SASA normalization or reference-energy lookup.

### Paths and overrides

Explicit relative paths in the YAML are interpreted from the shell's current working directory. The default output path is constructed from the configuration file's parent directory.

```bash
qprotsim analyze -c config.yaml --output-dir results/trp_cage
qprotsim analyze -c config.yaml --report results/trp_cage/custom_report.pdf
```

`--output-dir` overrides the YAML output directory. `--report` supplies the report filename and, when `--output-dir` is absent, its parent becomes the output directory. If both options are supplied, the report basename is written into `--output-dir`.

The Python wrapper follows the same priority for its `output_dir` and `report_filename` arguments. Use a basename such as `report.pdf` for `report_file_name` and `report_filename`.

### Advanced numerical settings

The YAML reader passes `pdb_path`, `protein_only`, `shots`, and `pl_device` to `Config`. It reads chain, plot bins, download settings, and output settings separately. Additional keys such as `temperature`, `beta`, `gamma`, `phi_cap`, or `weights` have no effect in this workflow.

The Python functions `compute_directed_transfer_schedule`, `prepare_state_propensity`, and `build_Kij_pairs` expose the corresponding lower-level controls. The high-level report uses their defaults. The `bins` setting reaches the energy plots. `run_analysis` calls `residue_information_gain` with that function's own default of 10,000 bins.

## A batch of PDB entries

```yaml
batch_mode: true
pdb_codes:
  - 1L2Y
  - 1UBQ
molecule_names:
  1L2Y: Trp Cage
  1UBQ: Ubiquitin
download_dir: ./input
chain: A
protein_only: true
shots: 10000
pl_device: default.qubit
bins: 100
output_dir: ./outputs/batch
report_file_name: "{pdb_code}_report.pdf"
```

```bash
qprotsim analyze -c batch.yaml
```

Batch mode accepts a `pdb_codes` list. After every download succeeds, each entry receives a separate analysis using the shared settings. `molecule_names` optionally maps the uppercase PDB codes returned by the downloader to display names.

| Batch key | Default | Meaning |
| --- | --- | --- |
| `batch_mode` | `false` | Enable the batch dispatcher. |
| `pdb_codes` | Required in batch mode | PDB identifiers to download and analyze. |
| `molecule_names` | PDB code as display name | Optional mapping from PDB code to display name. |
| `chain` | `A` | Shared chain selection. Use `null` for all chains. |
| `output_dir` | `<config_dir>/qprotsim_batch_outputs` | Base output directory. |
| `report_file_name` | `{pdb_code}_report.pdf` | Report filename pattern. |

Shared `protein_only`, `shots`, `pl_device`, `bins`, and download settings use the single-run defaults above. Batch mode creates lowercase PDB subdirectories, such as `outputs/batch/1l2y/1l2y_report.pdf`. The placeholders `{pdb_code}` and `{molecule_name}` are lowercased, and spaces in the molecule name become underscores.

Use `--output-dir` to override the batch output base directory. Batch report filenames come from `report_file_name`. The CLI's `--report` option applies to single analyses. A list of local-file configuration objects is currently unsupported by the batch dispatcher.

The downloader attempts all entries and raises an error if any download failed, stopping the batch before analysis begins. During the analysis stage, an individual protein's failure is printed and processing continues with the next entry. The return mapping contains completed reports. Check the final success count or the dictionary returned by `run_batch_analysis_from_yaml` before treating a batch as complete.

The configuration behavior is implemented in [`report_analysis.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/report_analysis.py) and [`__main__.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/__main__.py).
