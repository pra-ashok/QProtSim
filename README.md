# QProtSim

`QProtSim` is a Python package and command-line workflow for residue-level
thermodynamic analysis of protein structures using quantum-circuit sampling.
It maps a PDB structure to a residue graph, estimates solvent-accessible
surface area (SASA), builds residue-residue coupling terms, samples a
PennyLane circuit, and produces publication-oriented PDF, CSV, and figure
outputs for downstream analysis.

The default pipeline is designed for single-chain protein structures, but the
CLI also supports batch analysis from PDB identifiers.

## Scientific Abstract

Proteins encode stability, conformational response, and long-range residue
communication through a heterogeneous network of solvent exposure, residue
chemistry, and structural contacts. `QProtSim` models this network by assigning
each residue a two-state thermodynamic variable: a ground-like solvation state
and an excited solvation state. Solvent-accessible surface area provides the
local environmental signal, residue-specific transfer free energies provide the
energetic scale, and a structure-derived Kij coupling model links residues
through backbone connectivity, hydrogen bonding, salt bridges, van der Waals
contacts, water-mediated interactions, and mass compatibility.

The resulting residue graph is encoded as a quantum circuit in which local
rotations initialize residue excitation propensities and entangling layers
redistribute excitation according to Kij-weighted thermodynamic driving terms.
Sampling the circuit yields an ensemble of residue state configurations from
which `QProtSim` estimates total energy distributions, residue-level stability
and sensitivity metrics, pairwise mutual information, directional excitation
information, and joint-state frequencies. These outputs are intended to help
identify residues that are energetically stable, conformationally sensitive, or
informative for allosteric communication.

The framework is not a replacement for molecular dynamics or quantum chemistry.
Instead, it is a reduced thermodynamic and information-theoretic model for
screening protein residue networks, generating hypotheses about coupled sites,
and producing reproducible figures and tables suitable for manuscript-scale
analysis.

## Key Capabilities

- Residue-level SASA calculation with FreeSASA.
- Thermodynamic excitation propensity estimation from residue exposure and
  reference amino-acid energetics.
- Residue-residue coupling construction from backbone, hydrogen-bond, salt
  bridge, van der Waals, water-mediated, and mass-mobility terms.
- Quantum circuit sampling through PennyLane backends, including
  `lightning.gpu` when the CUDA environment is available.
- Energy ensemble analysis, residue KL metrics, mutual information, directional
  excitation information, joint-state frequencies, and motif allostery outputs.
- Automated PDF report generation with structure figures, plots, and summary
  tables.
- Batch workflows for multiple PDB codes with automatic RCSB PDB downloads.

## Installation

### Requirements

- Python 3.11 or newer.
- A working C/C++ compiler if Python wheels are unavailable for compiled
  dependencies.
- FreeSASA runtime support.
- CUDA libraries only if using GPU-backed PennyLane devices.

On Debian/Ubuntu systems, install FreeSASA and common build tools first:

```bash
sudo apt-get update
sudo apt-get install -y freesasa build-essential
```

Create and activate an environment:

```bash
conda create -n QProtSim python=3.11 -y
conda activate QProtSim
```

Install the package from the repository root:

```bash
cd /path/to/QProtSim
pip install .
```

For development and tests:

```bash
pip install -e ".[dev]"
```

### GPU Setup

`QProtSim` defaults to `lightning.gpu`. If the GPU plugin or CUDA runtime is not
available, the device creation code falls back to `default.qubit` for circuit
execution. For reproducible GPU execution, install CuPy and CUDA through conda:

```bash
conda install cupy cuda-version=13.1 -c conda-forge
```

The repository also includes `run_with_cuda.sh`, a local wrapper that sets CUDA
environment variables and then runs `QProtSim`. Edit the paths in that script if
your conda environment is located elsewhere.

## Quick Start

Create a YAML config file:

```yaml
pdb_path: /absolute/path/to/1L2Y.pdb
molecule_name: Trp Cage
chain: A
protein_only: true
shots: 100000
pl_device: lightning.gpu
bins: 10000
output_dir: /absolute/path/to/QProtSim_outputs
report_file_name: trp_cage_report.pdf
```

Run the analysis:

```bash
QProtSim analyze -c config.yaml
```

The command writes the report PDF and all companion artifacts into
`output_dir`.

You can also use a PDB code instead of a local file:

```yaml
pdb_code: 1L2Y
molecule_name: Trp Cage
download_dir: /absolute/path/to/pdb_cache
force_download: false
verify_checksum: true
chain: A
protein_only: true
shots: 100000
pl_device: lightning.gpu
bins: 10000
output_dir: /absolute/path/to/QProtSim_outputs
report_file_name: trp_cage_report.pdf
```

## Command-Line Usage

Show top-level help: QProtSim with alises as qprotsim, QPS, qps works to access package.

```bash
QProtSim --help 
qprotsim --help
QPS --help
qps --help
```

Run a single analysis:

```bash
qps analyze -c config.yaml
```

Override the output directory:

```bash
qps analyze -c config.yaml --output-dir results/trp_cage
```

Override the report path:

```bash
qps analyze -c config.yaml --report results/trp_cage/report.pdf
```

Enable verbose debug output from the information-theory and circuit routines:

```bash
qps analyze -c config.yaml --debug
```

Accepted aliases for `analyze` are `analyse`, `analysis`, and
`analyze-protein`.

### Preprocess Command

The CLI includes a preprocessing entry point for extracting molecular graph
features from a PDB file:

```bash
qps preprocess \
  --pdb input.pdb \
  --selection protein \
  --chain A \
  --output preprocessed.pkl
```

At present this command parses and reports the requested preprocessing options.
The full report-producing workflow is `QProtSim analyze`.

### Basin Command

Run basin analysis from sampled states and optional reference energy arrays:

```bash
qps basin \
  --z_samples outputs/trp_cage_z_samples.csv \
  --basins -60 -40 -40 -20 \
  --out basin_results.csv
```

Optional `--E0` and `--E1` arguments can point to `.npy` files. If omitted,
the package derives reference values from the internal residue energy table.

## Configuration Reference

Single-run configurations must include either `pdb_path` or `pdb_code`, plus
`molecule_name`.

| Key | Required | Default | Description |
| --- | --- | --- | --- |
| `pdb_path` | one of `pdb_path`/`pdb_code` | none | Local PDB file to analyze. |
| `pdb_code` | one of `pdb_path`/`pdb_code` | none | Four-character RCSB PDB ID to download before analysis. |
| `molecule_name` | yes | none | Human-readable name used in reports and artifact prefixes. |
| `chain` | no | `A` | Chain ID to analyze. Use `null`, empty string, or `--` to include all chains. |
| `protein_only` | no | `true` | Restrict analysis to standard protein residues. |
| `shots` | no | `100000` | Number of circuit samples. Increase for smoother ensemble statistics. |
| `pl_device` | no | `lightning.gpu` | PennyLane device name. Falls back to `default.qubit` if unavailable. |
| `bins` | no | `1000` | Number of bins for energy histograms and KL calculations. |
| `output_dir` | no | `<config_dir>/QProtSim_outputs` | Directory for report assets and CSV files. |
| `report_file_name` | no | `report.pdf` | PDF filename written at the root of `output_dir`. |
| `download_dir` | no | current working directory | Directory for downloaded PDB files when `pdb_code` is used. |
| `force_download` | no | `false` | Re-download the PDB file even if a cached file exists. |
| `verify_checksum` | no | `true` | Validate downloaded files when RCSB checksum metadata is available. |

## Batch Analysis

Batch mode downloads and analyzes multiple PDB entries from one YAML file:

```yaml
batch_mode: true
pdb_codes:
  - 1L2Y
  - 1AKG
  - 1B03

download_dir: /absolute/path/to/pdb_cache
output_dir: /absolute/path/to/batch_outputs
force_download: false
verify_checksum: true

molecule_names:
  1L2Y: Trp Cage
  1AKG: Protein 1AKG
  1B03: Protein 1B03

chain: A
protein_only: true
shots: 100000
pl_device: lightning.gpu
bins: 10000
report_file_name: "{pdb_code}_report.pdf"
```

Run:

```bash
qps analyze -c batch_config.yaml
```

Each protein is written to its own subdirectory under `output_dir`, for example
`batch_outputs/1l2y/1l2y_report.pdf`.

## Output Artifacts

A typical single analysis produces:

- `report.pdf` or the configured report filename.
- `<molecule>_residue_stats.csv` with residue index, residue name, normalized
  SASA, excitation propensity, and sampled ground-state probability.
- `<molecule>_residue_information_gain.csv` with normalized `KL_ground` and
  `KL_exc` residue metrics.
- `<molecule>_z_samples.csv` with per-shot residue states.
- `<molecule>_entanglement_parameters.csv` with the directed transfer schedule,
  gate type, residue labels, coupling strengths, transfer probabilities, and
  CRY/YY rotation angles used by the circuit.
- `<molecule>_kij_interaction_matrix.csv` and per-term Kij matrix CSV files.
- `<molecule>_directional_excitation_information.csv`.
- `<molecule>_mutual_information.csv`.
- Energy distribution and density plots.
- Residue KL, stability/sensitivity, coupling, information, joint-state, and
  ground-state comparison figures.
- Optional structure images from RCSB and PyMOL when image generation succeeds.
- Motif allostery files when the structure and sampled states satisfy the
  requirements for that analysis, including peptide dimer/trimer definitions,
  triangulated peptide/non-covalent triads, and separate dimer and trimer
  motif-region affinity matrix heatmaps.

Figure files are named with stable `figNN_...` prefixes so they can be cited
directly from the generated PDF or reused in manuscripts.

## Method Overview

The main analysis flow is:

1. Load a PDB file with MDAnalysis.
2. Compute normalized per-residue SASA using FreeSASA and Tien et al. maximum
   ASA reference values.
3. Convert SASA and residue hydrophobicity into residue excitation
   propensities.
4. Build a residue coupling graph and Kij terms from structural interactions.
5. Initialize a PennyLane circuit with residue excitation probabilities.
6. Apply directed transfer layers based on Kij-weighted thermodynamic driving
   terms.
7. Sample residue states and compute total energy ensembles.
8. Export residue statistics, information-theory matrices, plots, and a PDF
   report.

The sampled-state convention is:

- `+1` is the ground state.
- `-1` is the excited state.

Several exported information metrics convert this internally to a binary
encoding where `0` is ground and `1` is excited.

## Python API

The CLI is the recommended interface for full reports. Advanced users can call
the core functions directly:

```python
from pathlib import Path

from QProtSim.circuit import Config
from QProtSim.report_analysis import run_analysis, generate_report

cfg = Config(
    pdb_path="/absolute/path/to/1L2Y.pdb",
    protein_only=True,
    shots=100000,
    pl_device="lightning.gpu",
)

outputs = run_analysis(
    cfg=cfg,
    molecule_name="Trp Cage",
    chain="A",
    output_dir=Path("QProtSim_outputs"),
    bins=10000,
)

generate_report(
    output_path=Path("QProtSim_outputs/trp_cage_report.pdf"),
    molecule_name="Trp Cage",
    cfg=cfg,
    chain="A",
    outputs=outputs,
)
```

For a YAML-driven run from Python:

```python
from pathlib import Path
from QProtSim.report_analysis import run_analysis_from_yaml

report_path = run_analysis_from_yaml(Path("config.yaml"))
print(report_path)
```

## Function Parameter Reference

This section documents the main functions a user or developer is expected to
call directly. Functions whose names begin with `_` are internal implementation
helpers and should be treated as unstable.

### Core Configuration

`QProtSim.circuit.Config`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_path` | `str` | Path to the input PDB file used for structure loading, SASA calculation, residue graph construction, and reporting. |
| `protein_only` | `bool` | When `True`, restricts calculations to standard amino-acid residues. Set `False` only when nonstandard residues should remain in the model. |
| `shots` | `int` | Number of quantum-circuit samples. Larger values reduce sampling noise but increase runtime and memory use. |
| `pl_device` | `str` | PennyLane device name, for example `lightning.gpu` or `default.qubit`. GPU devices require the matching CUDA/PennyLane environment. |

### Report and Workflow Functions

`QProtSim.report_analysis.load_yaml_config(path)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `path` | `Path` | YAML config file to load. The file must contain a mapping of option names to values. |

`QProtSim.report_analysis.run_analysis_from_yaml(config_path, output_dir=None, report_filename=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `config_path` | `Path` | YAML file describing one analysis run. Must include `molecule_name` and either `pdb_path` or `pdb_code`. |
| `output_dir` | `Path` or `None` | Optional override for the artifact directory. Takes priority over `output_dir` in the YAML file. |
| `report_filename` | `str` or `None` | Optional report PDF filename. The report is written at the root of the resolved output directory. |

`QProtSim.report_analysis.run_batch_analysis_from_yaml(config_path, output_dir=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `config_path` | `Path` | YAML file with `batch_mode: true` and a `pdb_codes` list. Common analysis settings are applied to each protein. |
| `output_dir` | `Path` or `None` | Optional base directory for all batch outputs. Each PDB code receives its own subdirectory. |

`QProtSim.report_analysis.run_analysis(cfg, molecule_name, chain, output_dir, bins=1000, run_stamp=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `cfg` | `Config` | Runtime configuration containing PDB path, residue filtering, shot count, and PennyLane device. |
| `molecule_name` | `str` | Human-readable name used in plot titles, report text, and artifact prefixes. |
| `chain` | `str` or `None` | Chain identifier to analyze. Use `None` for all chains. |
| `output_dir` | `Path` | Directory where CSV files, figures, and intermediate analysis outputs are written. |
| `bins` | `int` | Number of bins used for energy histograms and KL calculations. |
| `run_stamp` | `str` or `None` | Optional timestamp string used for reproducible figure naming. If omitted, the function creates one. |

`QProtSim.report_analysis.generate_report(output_path, molecule_name, cfg, chain, outputs)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `output_path` | `Path` | Final PDF report path. |
| `molecule_name` | `str` | Name shown in the report title and captions. |
| `cfg` | `Config` | Analysis settings to summarize in the report. |
| `chain` | `str` or `None` | Chain label reported in the run summary. |
| `outputs` | `AnalysisOutputs` | Object returned by `run_analysis`; contains data frames, matrices, figures, and run metadata. |

`QProtSim.pipeline.run_pipeline(cfg, molecule_name, chain="A", plot_energy=True, bins=1, plot_path=None, save_scatter_plot=True, scatter_plot_path=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `cfg` | `Config` | Core run settings. |
| `molecule_name` | `str` | Name used in generated energy plots. |
| `chain` | `str` or `None` | Chain filter for SASA, residue order, and coupling construction. |
| `plot_energy` | `bool` | When `True`, writes energy distribution plots from sampled energies. |
| `bins` | `int` | Number of bins for energy plots. |
| `plot_path` | `str` or `None` | Optional path for the energy histogram. |
| `save_scatter_plot` | `bool` | When `True`, also saves the energy density plot. |
| `scatter_plot_path` | `str` or `None` | Optional path for the energy density plot. If omitted, a default path is derived from `molecule_name` or `plot_path`. |

`QProtSim.pipeline.export_residue_statistics(results, out_csv="residue_propensity.csv")`

| Parameter | Type | Usage |
| --- | --- | --- |
| `results` | `dict` | Output dictionary from `run_pipeline`; must include residue names, SASA values, excitation propensities, and sampled ground-state probabilities. |
| `out_csv` | `str` | Destination CSV path for per-residue statistics. |

`QProtSim.pipeline.export_entanglement_parameters(results, out_csv="entanglement_parameters.csv", yy_cap=np.pi/6, yy_probability_threshold=0.02, apply_yy_diffusion=True)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `results` | `dict` | Output dictionary from `run_pipeline`; must include `entanglement_schedule`, residue names, propensities, and residue energies. |
| `out_csv` | `str` | Destination CSV path for the circuit entanglement schedule. |
| `yy_cap` | `float` or `None` | Cap applied to derived IsingYY weak-diffusion angles. |
| `yy_probability_threshold` | `float` | Threshold used to label a schedule entry as `IsingYY` rather than `CRY_pair`. |
| `apply_yy_diffusion` | `bool` | Whether weak transfers are exported as YY diffusion gates. |

### PDB, Download, and SASA Helpers

`QProtSim.helper.compute_file_checksum(file_path, algorithm="sha256")`

| Parameter | Type | Usage |
| --- | --- | --- |
| `file_path` | `str` | File whose checksum should be computed. |
| `algorithm` | `str` | Hash algorithm accepted by `hashlib`, such as `md5`, `sha1`, or `sha256`. |

`QProtSim.helper.get_pdb_checksum(pdb_code)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_code` | `str` | Four-character RCSB PDB identifier. Used to query remote ETag/checksum metadata. |

`QProtSim.helper.verify_pdb_file(file_path, pdb_code, expected_checksum=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `file_path` | `str` | Local PDB file to validate. |
| `pdb_code` | `str` | PDB code used in warning and error messages. |
| `expected_checksum` | `str` or `None` | Optional expected checksum. When omitted, only existence and basic file-size checks are performed. |

`QProtSim.helper.fetch_pdb(pdb_code, output_dir=None, force_download=False, verify_checksum=True)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_code` | `str` | Four-character RCSB identifier to download. |
| `output_dir` | `str` or `None` | Directory for cached PDB files. Defaults to the current working directory. |
| `force_download` | `bool` | Re-download even if a local file already exists. |
| `verify_checksum` | `bool` | Validate downloaded or cached files when checksum metadata is available. |

`QProtSim.helper.batch_fetch_pdb(pdb_codes, output_dir=None, force_download=False, verify_checksum=True, delay=0.5)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_codes` | `list[str]` | PDB identifiers to download sequentially. |
| `output_dir` | `str` or `None` | Shared cache directory for all downloaded files. |
| `force_download` | `bool` | Re-download all requested structures even if cached. |
| `verify_checksum` | `bool` | Validate files when remote checksum metadata is available. |
| `delay` | `float` | Delay in seconds between downloads to avoid overloading RCSB services. |

`QProtSim.helper.per_residue_sasa(pdb_path, chain="A", protein_only=True)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_path` | `str` | PDB file used for FreeSASA calculation. |
| `chain` | `str` or `None` | Chain filter. Use `None` to include all chains. |
| `protein_only` | `bool` | When `True`, keep only standard amino-acid residues. |

`QProtSim.helper.map_sasa_to_mda_resorder(u, sasa_dict, protein_only=True, chain="A")`

| Parameter | Type | Usage |
| --- | --- | --- |
| `u` | `MDAnalysis.Universe` | Loaded structure whose residue order should define the output order. |
| `sasa_dict` | `dict` | SASA dictionary returned by `per_residue_sasa`. |
| `protein_only` | `bool` | Match the residue filter used during SASA calculation. |
| `chain` | `str` or `None` | Chain filter matching the analysis target. |

`QProtSim.helper.per_residue_ke(u, T=300.0, protein_only=True, seed=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `u` | `MDAnalysis.Universe` | Structure or trajectory used for residue masses and optional velocities. |
| `T` | `float` | Temperature in Kelvin for Maxwell-Boltzmann velocity sampling when velocities are absent. |
| `protein_only` | `bool` | Restrict kinetic-energy estimates to protein residues. |
| `seed` | `int` or `None` | Random seed for reproducible synthetic velocities. |

### Coupling and Circuit Functions

`QProtSim.coupling_factor.build_Kij_pairs(u, protein_only=True, chain="A", exclude_bonded=False, min_seq_separation=1, use_pbc=False, weights=None, omega=None, sigma_omega=None, estimate_omega_from_local_stiffness=True, normalize=True, return_terms=True)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `u` | `MDAnalysis.Universe` | Structure used to select residues and compute pairwise contacts. |
| `protein_only` | `bool` | Keep only standard amino-acid residues. |
| `chain` | `str` or `None` | Chain to analyze; `None` includes all selected residues. |
| `exclude_bonded` | `bool` | Skip neighboring residues when `True`. Keep `False` if backbone coupling should contribute. |
| `min_seq_separation` | `int` | Sequence-distance cutoff used only when `exclude_bonded=True`. |
| `use_pbc` | `bool` | Use periodic boundary dimensions from the MDAnalysis universe. |
| `weights` | `Mapping[str, float]` or `None` | Optional weights for coupling terms: `BB`, `HB`, `salt`, `vdW`, and `wat`. |
| `omega` | `Sequence`, `Mapping`, or `None` | Optional residue frequency values. If omitted, frequencies can be estimated from local stiffness. |
| `sigma_omega` | `float` or `None` | Width of the frequency-matching Gaussian. If omitted, an automatic value is estimated. |
| `estimate_omega_from_local_stiffness` | `bool` | Estimate relative residue frequencies from local coupling stiffness and mass when `omega` is omitted. |
| `normalize` | `bool` | Normalize final Kij values by the maximum Kij so couplings are on a comparable scale. |
| `return_terms` | `bool` | Return decomposed term values for reporting and heatmaps. |

`QProtSim.circuit.energies_ref(energies_ref="bulky")`

| Parameter | Type | Usage |
| --- | --- | --- |
| `energies_ref` | `str` | Select reference energy table. Use `"bulky"` for octanol-like transfer values or `"alkyl"` for alkyl-chain values. |

`QProtSim.circuit.compute_directed_transfer_schedule(exc_propensities, pairs, kijs, E0, E1, beta=1.0, eps=1.0e-12, gamma=0.95, max_layers=None, min_deltaP=0.0, phi_scale=1.0, phi_cap=np.pi/2, discount_start_power=1)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `exc_propensities` | `Sequence[float]` | Initial excitation probability for each residue/qubit. |
| `pairs` | `Sequence[tuple[int, int]]` | Undirected residue pairs that can exchange excitation. |
| `kijs` | `Mapping[tuple[int, int], float]` | Coupling strengths for residue pairs. |
| `E0` | `Sequence[float]` | Ground-state residue energies. |
| `E1` | `Sequence[float]` | Excited-state residue energies. |
| `beta` | `float` | Inverse-energy scaling factor used in thermodynamic driving terms. |
| `eps` | `float` | Numerical guard for divisions and logarithms. |
| `gamma` | `float` | Discount factor applied to later transfer layers. Must be in `(0, 1]`. |
| `max_layers` | `int` or `None` | Optional maximum number of directed transfers to keep. |
| `min_deltaP` | `float` | Ignore transfers at or below this probability threshold. |
| `phi_scale` | `float` | Multiplicative scale for CRY angles. |
| `phi_cap` | `float` or `None` | Maximum absolute CRY angle. Use `None` for no cap. |
| `discount_start_power` | `int` | Power used for the first gamma discount. Use `0` for undiscounted first layer or `1` for `gamma` on layer one. |

`QProtSim.circuit.prepare_state_propensity(exc_propensities, kij_pairs, kijs, E0, E1, yy_cap=np.pi/6, apply_yy_diffusion=True, return_schedule=True, yy_probability_threshold=0.02)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `exc_propensities` | `Sequence[float]` | Initial residue excitation probabilities. |
| `kij_pairs` | `Sequence[tuple[int, int]]` | Residue pairs included in entangling layers. |
| `kijs` | `Mapping[tuple[int, int], float]` | Kij coupling values for the pairs. |
| `E0` | `Sequence[float]` | Ground-state residue energies. |
| `E1` | `Sequence[float]` | Excited-state residue energies. |
| `yy_cap` | `float` or `None` | Maximum IsingYY angle for weak-diffusion corrections. |
| `apply_yy_diffusion` | `bool` | Use weak YY diffusion for small transfer drives instead of only CRY gates. |
| `return_schedule` | `bool` | Return the computed transfer schedule after applying gates. |
| `yy_probability_threshold` | `float` | Threshold below which a transfer is treated as weak diffusion. |

`QProtSim.circuit.sample_bits(exc_propensities, kij_pairs, kijs, E0, E1, cfg, return_schedule=False)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `exc_propensities` | `np.ndarray` | Per-residue excitation probabilities. |
| `kij_pairs` | `Sequence[tuple[int, int]]` | Residue pairs used to construct entangling layers. |
| `kijs` | `dict` | Coupling values for each pair. |
| `E0` | array-like | Ground-state energies. |
| `E1` | array-like | Excited-state energies. |
| `cfg` | `Config` | Supplies shot count and PennyLane device. |
| `return_schedule` | `bool` | When `True`, returns `(samples, entanglement_schedule)` instead of only the sample matrix. |

`QProtSim.circuit.energy_from_samples(Z_samples, E0, E1)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `Z_samples` | `np.ndarray` | Sample matrix of shape `(shots, n_res)`, where `+1` is ground and `-1` is excited. |
| `E0` | `np.ndarray` | Ground-state residue energies. |
| `E1` | `np.ndarray` | Excited-state residue energies. |

`QProtSim.circuit.per_residue_energies(resnames)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `resnames` | `Sequence[str]` | Three-letter residue names in analysis order. |

Additional public circuit helpers:

| Function | Parameters | Usage |
| --- | --- | --- |
| `residue_coms(u, protein_only=True, chain="A")` | `u`: MDAnalysis universe; `protein_only`: standard residues only; `chain`: chain filter. | Returns residue centers of mass in analysis order. |
| `build_peptide_pairs(n_res)` | `n_res`: number of residues. | Returns adjacent sequence pairs `(i, i + 1)`. |
| `build_disulfide_pairs(u, protein_only=True, chain="A", max_sg_distance=2.2)` | `u`: structure; `protein_only`: residue filter; `chain`: chain filter; `max_sg_distance`: SG-SG cutoff. | Detects cysteine pairs consistent with disulfide bonds. |
| `build_vdW_pairs(u, protein_only=True, chain="A", exclude_bonded=True, min_seq_separation=1, use_pbc=False)` | Structure and filtering parameters plus contact/exclusion settings. | Legacy contact-weight pair builder; the main pipeline now uses `build_Kij_pairs`. |
| `sigmoid(x, bias=-2.4)` | `x`: numeric input; `bias`: horizontal shift. | Converts a score to a bounded probability-like value. |
| `norm01(x)` | `x`: numeric array. | Normalizes an array by its maximum value. |
| `kd_to_scores(resnames)` | `resnames`: residue names. | Converts Kyte-Doolittle hydrophobicity to hydrophobic/hydrophilic score arrays. |
| `is_hydrophobic(resnames)` | `resnames`: residue names. | Returns a boolean hydrophobicity mask. |
| `angles_from_probability(p)` | `p`: probability array. | Converts probabilities to RY rotation angles. |
| `boltzmann_expectation(delta_E, T=298.0, k_B=8.314462618e-3)` | `delta_E`: energy gap; `T`: temperature; `k_B`: Boltzmann constant in matching units. | Computes two-state Boltzmann probabilities and expected energy. |
| `excitation_propensities(SASA, etype=1)` | `SASA`: per-residue normalized SASA dictionary; `etype`: reserved energy-type selector. | Converts residue exposure and hydrophobicity into excitation probabilities. |
| `make_device(n_qubits, name)` | `n_qubits`: wire count; `name`: PennyLane device name. | Creates a PennyLane device and falls back to `default.qubit` if the requested device fails. |

### Information-Theory and Allostery Functions

`QProtSim.info_theory.residue_information_gain(Z_samples, E_total, bins=10000)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `Z_samples` | `np.ndarray` | Sample matrix of residue states. Accepts `+1/-1` or `0/1` encodings. |
| `E_total` | `np.ndarray` | Total energy for each sampled configuration. Must align with the rows of `Z_samples`. |
| `bins` | `int` | Number of energy bins for KL divergence estimates. |

`QProtSim.info_theory.directional_excitation_information(Z_samples, threshold=0.0)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `Z_samples` | `np.ndarray` | Sample matrix of shape `(shots, n_res)`. |
| `threshold` | `float` | Optional noise floor subtracted from each directional score after KL calculation. |

`QProtSim.info_theory.mutual_information_matrix(Z_samples, eps=1e-12)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `Z_samples` | `np.ndarray` | Sample matrix used to compute pairwise residue mutual information. |
| `eps` | `float` | Small value added for numerical stability in probability ratios. |

`QProtSim.info_theory.joint_state_matrices(Z_samples)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `Z_samples` | `np.ndarray` | Sample matrix used to count `(0,0)`, `(0,1)`, `(1,0)`, and `(1,1)` pair states. |

`QProtSim.info_theory.compute_static_mi_from_counts(counts_00, counts_01, counts_10, counts_11, eps=1e-12)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `counts_00` | `np.ndarray` | Pairwise count/probability matrix for both residues ground. |
| `counts_01` | `np.ndarray` | Pairwise count/probability matrix for first ground, second excited. |
| `counts_10` | `np.ndarray` | Pairwise count/probability matrix for first excited, second ground. |
| `counts_11` | `np.ndarray` | Pairwise count/probability matrix for both residues excited. |
| `eps` | `float` | Numerical guard for logarithms and divisions. |

`QProtSim.info_theory.analyse_motif_allostery(...)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `Z_samples` | `np.ndarray` | Required sample matrix of residue states. |
| `residue_labels` | `Sequence[str]` or `None` | Optional labels for residues. Defaults to `R1`, `R2`, ... |
| `save_dir` | `str` or `None` | Directory for CSV and text outputs. If omitted, results are returned but not saved. |
| `save_prefix` | `str` | Prefix for motif allostery output filenames, default `motif_allostery`. |
| `random_state` | `int` | Seed for the k-means step inside spectral clustering. |
| `pdb_path` | `str` | PDB path used to derive peptide and non-covalent motif regions from the Kij decomposition. |
| `chain` | `str` or `None` | Chain filter for motif construction. |
| `protein_only` | `bool` | Keep only standard amino-acid residues during motif construction. |
| `min_support` | `int` | Legacy alias used when `min_support_res` is omitted. |
| `min_support_res` | `int` or `None` | Minimum number of excited observations required for a residue to act as a supported source. |
| `motif_agg_mode` | `str` | Aggregation mode for residue-to-motif DEI blocks: `mean`, `sum_abs`, or `topk_mean`. |
| `motif_agg_topk` | `int` | Number of largest absolute block values used by `topk_mean`. |
| `affinity_scale_percentile` | `float` | Percentile used to robustly scale the symmetric motif affinity matrix. |
| `n_clusters` | `int` or `None` | Number of motif clusters for spectral clustering. If omitted, a bounded automatic value is used. |
| `peptide_threshold` | `float` | Minimum `kappa_BB` value required for a peptide edge. |
| `noncovalent_threshold` | `float` | Minimum non-covalent edge strength from `kappa_HB`, `kappa_salt`, or `kappa_vdW`. |
| `include_dimers` | `bool` | Include peptide-bonded residue pairs as dimer motif regions. |
| `include_peptide_trimers` | `bool` | Include consecutive peptide-bonded triples as trimer motif regions. |
| `include_triangulated_triads` | `bool` | Include triads with at least one peptide edge and at least one non-covalent edge. |
| `save_affinity_plot` | `bool` | Save motif-region affinity heatmaps as PNG and PDF when `save_dir` is provided. |
| `eps` | `float` | Numerical guard for probabilities and normalizations. |

### Plotting and Structure Rendering Functions

Most plotting functions return a Matplotlib figure and axes object and save
both PNG and PDF versions when `save_path` ends in `.png`.

Common plotting parameters:

| Parameter | Usage |
| --- | --- |
| `info_df` | Residue information DataFrame, usually from `residue_information_gain`. |
| `results_df` | Per-residue statistics DataFrame from `export_residue_statistics`. |
| `results` | Pipeline result dictionary containing arrays such as `p_ground`, `p_exc`, and `E`. |
| `molecule_name` | Name used in plot titles. |
| `figsize` / `fig_size` | Figure dimensions in inches. |
| `save_path` | Destination figure path. When provided, the figure is saved. |
| `fig_num` | Optional figure number used by report-generation code. |
| `labels` | Residue labels for heatmaps and pairwise matrices. |
| `matrix` / `DEI` / `matrices` | Numeric matrix or matrices to visualize. |
| `color1`, `color2`, `point_color`, `connector_color` | Plot color controls. |
| `dpi` | Output resolution for saved figures. |

Individual plotting functions:

| Function | Parameters and Usage |
| --- | --- |
| `plot_directional_excitation_information_clustermap(csv_path, save_png=None, save_pdf=None, figsize=(7, 7), cmap=HEATMAP_CMAP)` | `csv_path` loads a DEI matrix CSV; `save_png`/`save_pdf` write files; `figsize` controls dimensions; `cmap` controls colors. |
| `set_academic_style()` | Takes no parameters; applies shared seaborn/matplotlib styling for report figures. |
| `plot_energy_distribution(E_total, molecule_name, bins=1000, figsize=(8, 4), save_path=None, fig_num=None)` | `E_total` is the energy vector; `molecule_name` labels the plot; `bins` controls histogram resolution; `save_path` writes the figure. |
| `plot_energy_distribution_histogram(E_total, molecule_name, bins=1000, save_path=None, fig_num=None)` | Compact histogram variant used by report generation. |
| `plot_energy_density(E_total, molecule_name, bins=1000, figsize=(8, 6), save_path=None, fig_num=None)` | Estimates and plots the sampled energy density. |
| `plot_energy_sample_scatter(E_total, molecule_name, bins=1000, figsize=(8, 6), save_path=None, fig_num=None)` | Compatibility wrapper that currently calls `plot_energy_density`. |
| `plot_dual_axis_comparison(info_df, results_df, results, molecule_name="Trp Cage", metric1_col="p_ground", metric2_col="KL_exc", metric1_label=..., metric2_label=..., color1=..., color2=..., figsize=(10, 5), save_path=..., dpi=600, fig_num=None)` | Compares one metric from `results` with one metric from `info_df` on two y-axes. Metric column and label parameters select what is plotted. |
| `plot_residue_information_gain(info_df, results_df, molecule_name, figsize=(10, 4), save_path=None, fig_num=None)` | Plots normalized per-residue `KL_ground` and `KL_exc`. |
| `plot_kl_divergence_scatter(info_df, results_df, molecule_name, figsize=(6.5, 6), save_path=None, fig_num=None)` | Scatter plot comparing `KL_ground` with `KL_exc` for each residue. |
| `plot_kl_divergence_line_comparison(info_df, results_df, molecule_name, color1=..., color2=..., connector_color=..., connector_alpha=0.6, connector_linewidth=0.9, figsize=(10, 5), save_path=None, fig_num=None)` | Line plot showing the difference between ground and excited KL contribution per residue. |
| `plot_stability_vs_information(info_df, results_df, results, molecule_name, figsize=(10, 4), save_path=None, fig_num=None)` | Bar plot comparing sampled ground-state probability with excited-state information contribution. |
| `plot_directional_excitation_information_heatmap(DEI, labels, figsize=(6, 5), save_path=None, fig_num=None)` | Heatmap for a directional excitation information matrix. |
| `plot_pairwise_matrix_heatmap(matrix, labels, title, colorbar_label, cmap=HEATMAP_CMAP, fig_size=None, save_path=None, fig_num=None, font_size=10, annotate=False)` | General heatmap for Kij, mutual information, or other residue-residue matrices. `annotate=True` writes values into cells. |
| `plot_joint_state_frequencies(matrices, titles, labels, figsize=(10, 8), save_path=None, fig_num=None)` | Plots the four joint-state frequency matrices in a 2x2 panel. |
| `plot_ground_state_comparison(results_df, results, molecule_name, figsize=(10, 4), save_path=None, fig_num=None)` | Compares per-residue ground-state probability before and after entangling layers. |
| `plot_dual_axis_line_comparison(info_df, results_df, results, molecule_name="Trp Cage", metric1_col="p_ground", metric2_col="KL_exc", metric1_label=..., metric2_label=..., color1=..., color2=..., connector_color=..., connector_alpha=0.6, connector_linewidth=0.9, figsize=(10, 5), save_path=..., dpi=600, fig_num=None)` | Dual-axis line version of the stability/sensitivity comparison. |
| `plot_stability_vs_sensitivity_quadrants(info_df, results_df, results, molecule_name="Trp Cage", metric1_col="p_ground", metric2_col="KL_exc", metric1_label=..., metric2_label=..., quadrant_method="mean", x_threshold=None, y_threshold=None, point_color=..., quadrant_line_color=..., annotate_residues=True, figsize=(7, 7), save_path=..., dpi=600, fig_num=None)` | Quadrant scatter plot. `quadrant_method`, `x_threshold`, and `y_threshold` control the reference lines. |

`QProtSim.structure_viz.get_structure_image(pdb_path, output_dir, pdb_code=None, chain="A", prefer_download=True, output_path=None)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_path` | `Path` | Local PDB file to render if download is unavailable or disabled. |
| `output_dir` | `Path` | Directory for generated structure images. |
| `pdb_code` | `str` or `None` | Optional PDB identifier used to download an RCSB structure image. |
| `chain` | `str` or `None` | Chain filter for fallback matplotlib rendering. |
| `prefer_download` | `bool` | Try RCSB image download before local rendering. |
| `output_path` | `Path` or `None` | Optional exact image path. If omitted, a default filename is used. |

`QProtSim.structure_viz.render_structure_pymol(pdb_path, output_path, width=800, height=600, style="cartoon", color_scheme="spectrum", ray_trace=True)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `pdb_path` | `Path` | Input PDB file. |
| `output_path` | `Path` | Output image path; PNG is written regardless of suffix. |
| `width` | `int` | Image width in pixels. |
| `height` | `int` | Image height in pixels. |
| `style` | `str` | PyMOL representation such as `cartoon`, `sticks`, `spheres`, `ribbon`, or `surface`. |
| `color_scheme` | `str` | Coloring mode: `spectrum`, `chain`, `secondary`, or `element`. |
| `ray_trace` | `bool` | Use PyMOL ray tracing for higher quality output. |

`QProtSim.structure_viz.optimize_image_for_report(image_path, max_width=800, max_height=600, quality=85)`

| Parameter | Type | Usage |
| --- | --- | --- |
| `image_path` | `Path` | Image to resize and optimize. |
| `max_width` | `int` | Maximum output width in pixels. |
| `max_height` | `int` | Maximum output height in pixels. |
| `quality` | `int` | Output image quality setting, from 1 to 100. |

## Development

Install development dependencies:

```bash
pip install -e ".[dev]"
```

Run tests when a `tests/` directory is present:

```bash
pytest
```

Or use the bundled runner:

```bash
python run_tests.py --fast
python run_tests.py --coverage
```

The Makefile exposes common development commands:

```bash
make test
make test-coverage
make format-check
make lint
make type-check
```

## Troubleshooting

### FreeSASA installation fails

Install system FreeSASA and compiler tooling first:

```bash
sudo apt-get install -y freesasa build-essential
```

Then reinstall:

```bash
pip install -e .
```

### GPU device is unavailable

Confirm that CUDA, CuPy, and `pennylane_lightning_gpu` are installed in the
active environment. For CPU execution, set:

```yaml
pl_device: default.qubit
```

### PDB download fails

Check the four-character PDB code, network access, and `download_dir`
permissions. To reuse a local file instead of downloading, provide `pdb_path`
in the config.

### No residues are found

Check that the requested `chain` matches the PDB file. To analyze all chains,
set:

```yaml
chain: null
```

### PyMOL structure rendering is skipped

The analysis continues if PyMOL rendering fails. Install or repair
`pymol-open-source` if PyMOL figures are required.

## License

This project is distributed under the MIT License. See `LICENSE`.

## Citation

If you use `QProtSim` in your research/publication, cite the repository and include the exact
software version or commit hash used to generate results.
Repo: https://github.com/pra-ashok/QProtSim 

Publication citation: 
Pratik Patil, Bhushan Bonde and Bhaskar Choubey, Quantum Circuit–Based Protein Modelling Reveals Distal Energetic Couplings in Protein, (2026), Preprint (To be communicated).
