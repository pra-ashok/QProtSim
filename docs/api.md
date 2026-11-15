# Python API

The import namespace is `qprotsim`. The command aliases `QProtSim`, `QPS`, and `qps` all call the same CLI.

## Generate an analysis and a PDF

`run_analysis` is the main entry point for working with the results in Python. It samples the circuit, writes CSV files and figures, and returns an `AnalysisOutputs` object. Call `generate_report` with that object to write the PDF.

```python
from pathlib import Path

from qprotsim.circuit import Config
from qprotsim.report_analysis import generate_report, run_analysis

cfg = Config(
    pdb_path="input/1L2Y.pdb",
    protein_only=True,
    shots=10_000,
    pl_device="default.qubit",
)
output_dir = Path("outputs/trp_cage")

outputs = run_analysis(
    cfg=cfg,
    molecule_name="Trp Cage",
    chain="A",
    output_dir=output_dir,
    bins=100,
)

generate_report(
    output_path=output_dir / "report.pdf",
    molecule_name="Trp Cage",
    cfg=cfg,
    chain="A",
    outputs=outputs,
)

print(outputs.results_df.head())
print(outputs.info_df.head())
print(outputs.results["Z_samples"].shape)
```

The PDB file must exist before this example runs. Use `fetch_pdb` below or a YAML configuration with `pdb_code` to download it.

### `Config`

```text
Config(
    pdb_path: str,
    protein_only: bool = True,
    shots: int = 100_000,
    pl_device: str = "default.qubit",
)
```

`pdb_path` identifies the structure. `protein_only` filters the residue selection. `shots` controls the number of sampled configurations. `pl_device` chooses the PennyLane device. Chain selection belongs to `run_analysis` or `run_pipeline`.

### `run_analysis`

```text
run_analysis(
    cfg: Config,
    molecule_name: str,
    chain: str | None,
    output_dir: Path,
    bins: int = 1000,
    run_stamp: str | None = None,
) -> AnalysisOutputs
```

| Argument | Meaning |
| --- | --- |
| `cfg` | Structure, residue filter, shot count, and device settings. |
| `molecule_name` | Display name and artifact prefix. Use a simple name containing letters, numbers, spaces, or underscores. |
| `chain` | Chain ID, or `None` to select all chains. |
| `output_dir` | `pathlib.Path` for the output directory, created as needed. |
| `bins` | Resolution of energy plots. The residue information calculation uses its separate default of 10,000 bins. |
| `run_stamp` | Optional report timestamp string, normally `YYYYMMDD_HHMMSS`. |

`run_analysis` first attempts structure images, then runs the numerical pipeline, exports residue and circuit data, calculates information matrices, and generates plots. Motif allostery runs with `min_support=20`. A failed motif analysis emits a warning and the remaining report artifacts continue.

### `AnalysisOutputs`

For `S` shots and `N` selected residues, the returned dataclass has these fields.

| Field | Contents |
| --- | --- |
| `results` | Raw pipeline dictionary described below. |
| `results_df` | `N` rows with `res_index`, `resname`, `SASA_norm`, `p_exc`, and `p_ground`. |
| `info_df` | `N` rows with `res_index`, `KL_ground`, `KL_exc`, and the merged `resname`. |
| `dei_matrix` | `N × N` directional excitation information in bits, indexed as `[target, source]`. |
| `mi_matrix` | `N × N` symmetric mutual information in bits. |
| `joint_counts` | Four `N × N` arrays containing joint-state probabilities, ordered `00`, `01`, `10`, `11`. |
| `plot_paths` | Mapping from plot keys to `Path` objects. Some keys refer to the same figure. |
| `fig_numbers` | Mapping from plot keys to report figure numbers. |
| `run_stamp` | Timestamp used in the report summary. |

Both information matrices have zero diagonals. Joint-state arrays also have zero diagonals. Residue indices in the data frames are sequential, one-based positions in the selected analysis order. Recover original PDB residue numbers and chain IDs from the input structure when mapping results back to atoms.

### `generate_report`

```text
generate_report(
    output_path: Path,
    molecule_name: str,
    cfg: Config,
    chain: str | None,
    outputs: AnalysisOutputs,
) -> None
```

This function assembles the existing results and image files into a ReportLab PDF. Ensure the parent of `output_path` exists. The `run_analysis` call above already creates `output_dir`.

## Run a YAML workflow

```python
from pathlib import Path
from qprotsim.report_analysis import run_analysis_from_yaml

report_path = run_analysis_from_yaml(
    config_path=Path("config.yaml"),
    output_dir=Path("outputs/trp_cage"),
    report_filename="report.pdf",
)
print(report_path)
```

```text
run_analysis_from_yaml(
    config_path: Path,
    output_dir: Path | None = None,
    report_filename: str | None = None,
) -> Path

run_batch_analysis_from_yaml(
    config_path: Path,
    output_dir: Path | None = None,
) -> dict[str, Path]

load_yaml_config(path: Path) -> dict
```

`run_analysis_from_yaml` returns the PDF path after generating all outputs. Function arguments override YAML output settings. `run_batch_analysis_from_yaml` returns a mapping from successful PDB codes to their PDF paths. When `batch_mode` is false, its return value is `{"single": report_path}`. See the [configuration reference](configuration.md) for supported keys, path resolution, and batch behavior.

## Use the numerical pipeline

```python
from qprotsim.pipeline import run_pipeline

results = run_pipeline(
    cfg=cfg,
    molecule_name="Trp Cage",
    chain="A",
    plot_energy=False,
)

samples = results["Z_samples"]
energies = results["E"]
print(samples.shape, energies.shape)
```

```text
run_pipeline(
    cfg: Config,
    molecule_name: str,
    chain: str | None = "A",
    plot_energy: bool = True,
    bins: int = 1,
    plot_path: str | None = None,
    save_scatter_plot: bool = True,
    scatter_plot_path: str | None = None,
) -> dict
```

`plot_energy=True` enables an energy histogram. `save_scatter_plot=True` additionally generates the energy density plot. `plot_path` and `scatter_plot_path` supply the corresponding output paths. Pass an explicit `bins` value such as `100` when enabling plots.

| Result key | Shape or type | Meaning |
| --- | --- | --- |
| `Z_samples` | `(S, N)` array | Pauli-Z outcomes. `+1` means ground and `-1` means excited. |
| `E` | `(S,)` array | Total energy of each sampled configuration, in kJ mol⁻¹. |
| `E0`, `E1` | `(N,)` arrays | Ground and excited reference energies, in kJ mol⁻¹. |
| `p_exc` | Length-`N` pandas Series | Initial excitation propensity from exposure and residue chemistry. |
| `p_ground` | `(N,)` array | Ground-state frequency measured after the coupling circuit. |
| `resnames` | Length-`N` list | Three-letter names in analysis order. |
| `alpha_SASA` | `(N,)` array | Relative SASA, capped at `0.999`. |
| `pairs` | List of `(i, j)` tuples | Coupled zero-based residue positions. |
| `kijs` | Mapping from pair to float | Normalized structural coupling values. |
| `kij_terms` | Mapping from pair to dictionary | Coupling component decomposition and intermediate factors. |
| `entanglement_schedule` | List of dictionaries | Ordered directed transfers with drive, probability, and rotation parameters. |

`p_exc` describes initialization, while `1 - p_ground` describes the sampled excitation frequency after the coupling gates.

### Export tables

```python
from qprotsim.pipeline import (
    export_entanglement_parameters,
    export_residue_statistics,
)

residue_table = export_residue_statistics(
    results, out_csv="residue_statistics.csv"
)
gate_table = export_entanglement_parameters(
    results, out_csv="entanglement_parameters.csv"
)
```

Each function writes a CSV and returns its pandas data frame. The residue export defaults to `residue_propensity.csv`. The circuit export defaults to `entanglement_parameters.csv` and uses `yy_cap=pi/6`, `yy_probability_threshold=0.02`, and `apply_yy_diffusion=True`. These values match the normal sampling circuit. When supplying a custom circuit, pass its settings to the exporter as well.

## Calculate information from samples

```python
from qprotsim.info_theory import (
    directional_excitation_information,
    joint_state_matrices,
    mutual_information_matrix,
    residue_information_gain,
)

mi = mutual_information_matrix(samples)
dei = directional_excitation_information(samples, threshold=0.0)
p00, p01, p10, p11 = joint_state_matrices(samples)
info = residue_information_gain(samples, energies, bins=100)
```

| Function | Return value and conventions |
| --- | --- |
| `mutual_information_matrix(Z_samples, eps=1e-12)` | Symmetric `N × N` mutual information in bits. Accepts Pauli-Z or binary excitation samples. |
| `directional_excitation_information(Z_samples, threshold=0.0)` | `N × N` conditional KL scores in bits, with target rows and excited-source columns. The threshold is subtracted and results are clamped at zero. |
| `joint_state_matrices(Z_samples)` | Tuple of four joint probability matrices. Requires Pauli-Z samples. The `0` and `1` labels mean ground and excited. |
| `residue_information_gain(Z_samples, E_total, bins=10000)` | Data frame containing `res_index`, `KL_ground`, and `KL_exc`. Each KL column is normalized across residues. |

For the first, second, and fourth functions, binary input uses `0=ground`, `1=excited`. Encoding is inferred from whether any negative value occurs. Explicit binary conversion avoids ambiguity for a sample matrix containing only `+1` Pauli-Z outcomes.

```python
binary_excited = (samples == -1).astype("uint8")
mi = mutual_information_matrix(binary_excited)
dei = directional_excitation_information(binary_excited)
```

The residue information implementation requires observed samples in both states for each residue and a nonzero total in each KL column to produce finite normalized values. Sparse or deterministic samples can produce `NaN`. Check `info[["KL_ground", "KL_exc"]].isna()` before interpreting rankings. See [the theory](theory.md) for the equations and units.

### Motif allostery

```python
from qprotsim.info_theory import analyse_motif_allostery

labels = [f"{i + 1}-{name}" for i, name in enumerate(results["resnames"])]
motifs = analyse_motif_allostery(
    Z_samples=samples,
    residue_labels=labels,
    pdb_path=cfg.pdb_path,
    chain="A",
    protein_only=cfg.protein_only,
    min_support=20,
    save_dir="outputs/motifs",
)
```

The structure supplies peptide dimers, peptide trimers, and noncovalent triads. This workflow calculates the signed excitation contribution `q_exc * ln(q_exc / p_exc)` in nats before aggregating regions. Default aggregation uses `motif_agg_mode="mean"`, `affinity_scale_percentile=95.0`, and automatic cluster count. `min_support` sets the minimum excited-source sample count, with `min_support_res` available as an explicit override. Returned arrays include `residue_allostery_matrix`, `region_to_residue_matrix`, `membership_matrix`, `directional_region_coupling_matrix`, `affinity_matrix`, and `cluster_labels`. The dictionary also provides data-frame views and `saved_paths` for files written when `save_dir` is supplied.

## Fetch and inspect structures

```python
from qprotsim.helper import fetch_pdb, per_residue_sasa

pdb_path = fetch_pdb("1L2Y", output_dir="input")
relative_sasa = per_residue_sasa(str(pdb_path), chain="A")
```

`fetch_pdb(pdb_code, output_dir=None, force_download=False, verify_checksum=True)` returns a `Path`, reusing a cached file when possible. `batch_fetch_pdb(pdb_codes, output_dir=None, force_download=False, verify_checksum=True, delay=0.5)` returns a mapping from uppercase PDB codes to paths when all downloads succeed. It attempts each requested download and raises `RuntimeError` if any download failed.

`per_residue_sasa(pdb_path, chain="A", protein_only=True)` returns a dictionary keyed by `(chain_label, resid, resname)`. Values are normalized relative SASA. `map_sasa_to_mda_resorder(u, sasa_dict, protein_only=True, chain="A")` aligns those values to an MDAnalysis universe and returns `(values, residue_names)`.

## Inspect coupling and energy calculations

```python
import MDAnalysis as mda
from qprotsim.coupling_factor import build_Kij_pairs

universe = mda.Universe(cfg.pdb_path)
pairs, couplings, components = build_Kij_pairs(
    universe,
    chain="A",
    protein_only=True,
    normalize=True,
    return_terms=True,
)
```

`build_Kij_pairs` returns zero-based residue pairs and their coupling values. With `return_terms=True`, the third result contains interaction components and intermediate factors. The function accepts custom `weights`, `omega`, and `sigma_omega`. Its defaults include `exclude_bonded=False`, `min_seq_separation=1`, `use_pbc=False`, and `estimate_omega_from_local_stiffness=True`. Keep the selection and chain settings aligned with those used to create the sample columns.

```python
from qprotsim.circuit import energy_from_samples, per_residue_energies

ground_energies, excited_energies = per_residue_energies(results["resnames"])
configuration_energies = energy_from_samples(
    samples, ground_energies, excited_energies
)
```

`per_residue_energies(resnames)` selects the default octanol reference table in the supplied residue order. `energies_ref(energies_ref="bulky")` returns that three-letter-name lookup table, and `energies_ref("alkyl")` returns the alternate alkyl-chain table. `energy_from_samples` requires Pauli-Z samples and sums ground or excited residue energies for each configuration.

The implementation of these interfaces is in [`report_analysis.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/report_analysis.py), [`pipeline.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/pipeline.py), [`info_theory.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/info_theory.py), and [`helper.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/helper.py).
