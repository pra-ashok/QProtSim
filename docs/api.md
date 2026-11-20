# Using QProtSim from Python

Use Python when you want to inspect tables, work with sample arrays, or combine QProtSim with your own analysis. The package name in imports is `qprotsim`.

If you prefer running commands, start with the [first-analysis guide](getting-started.md). You can also run the same configuration file from Python, as shown below. The later sections cover direct access to results and information measures.

## Run an existing configuration

This is the shortest route to a complete report. Create `config.yaml` using the [configuration guide](configuration.md), then run:

```python
from pathlib import Path
from qprotsim.report_analysis import run_analysis_from_yaml

report_path = run_analysis_from_yaml(Path("config.yaml"))
print(report_path)
```

The function writes the PDF, CSV tables, and figures, then returns the PDF's path. To override the output settings in the file:

```python
report_path = run_analysis_from_yaml(
    config_path=Path("config.yaml"),
    output_dir=Path("outputs/repeat"),
    report_filename="report.pdf",
)
```

Use a new output directory when you want to keep earlier results. For a batch configuration, use `run_batch_analysis_from_yaml(Path("batch.yaml"))` from the same module. It returns a dictionary mapping successful PDB codes to report paths.

## Run an analysis and inspect its results

This example downloads a small structure, samples the model, and creates a report. An internet connection is needed for the download. For your own local file, replace the `fetch_pdb` call with `Path("input/my_protein.pdb")`.

```python
from pathlib import Path

from qprotsim.circuit import Config
from qprotsim.helper import fetch_pdb
from qprotsim.report_analysis import generate_report, run_analysis

pdb_path = fetch_pdb("1L2Y", output_dir="input")
cfg = Config(
    pdb_path=str(pdb_path),
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

`run_analysis` returns the data and writes tables and figures. `generate_report` uses those results to create the PDF. The output directory is created automatically by `run_analysis`.

### `Config`

`Config` holds the structure and sampling settings.

| Setting | Default | Meaning |
| --- | --- | --- |
| `pdb_path` | Required | Path to an existing PDB file, passed as a string. |
| `protein_only` | `True` | Select protein residues. Standard amino acids have the reference values needed by the model. |
| `shots` | `100_000` | Number of sampled configurations. |
| `backend` | `"pennylane"` | Circuit framework: `"pennylane"` or `"qiskit"`. |
| `pl_device` | `"default.qubit"` | Device used by PennyLane. See [installation](installation.md). |
| `qiskit_backend` | `"aer_simulator"` | Local Aer backend used by Qiskit. |
| `backend_options` | `{}` | Mapping of simulator options. |
| `seed` | `None` | Optional nonnegative integer for reproducible sampling. |

Pass the chain separately to `run_analysis` or `run_pipeline`. Use `chain=None` to select all chains; review the selected size before sampling.

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

`molecule_name` appears in the report and exported filenames. Use a simple name containing letters, numbers, spaces, or underscores. `output_dir` should be a `pathlib.Path`. The optional `run_stamp` labels the report; it does not give each exported file a unique name.

`bins` controls the energy plots. The residue information scores use a separate default of 10,000 histogram bins. Motif analysis requires at least 20 excited-source observations by default. If that optional stage fails, the remaining analysis can still finish.

### Read `AnalysisOutputs`

For `S` shots and `N` selected residues:

| Field | What it contains |
| --- | --- |
| `results_df` | One row per residue: `res_index`, `resname`, `SASA_norm`, `p_exc`, and `p_ground`. |
| `info_df` | Residue names and normalized `KL_ground` and `KL_exc` sensitivity scores. |
| `mi_matrix` | Symmetric `N × N` mutual information matrix, in bits. |
| `dei_matrix` | `N × N` directional excitation information matrix, in bits; rows are targets and columns are sources. |
| `joint_counts` | Four `N × N` joint probability arrays, ordered `00`, `01`, `10`, `11`. Despite the field name, these contain probabilities. |
| `results` | Samples, energies, and model parameters; the main keys are described below. |
| `plot_paths` | Dictionary of generated figure paths. |
| `fig_numbers` | Dictionary of figure numbers used by the report. |
| `run_stamp` | Report timestamp. |

The MI, DEI, and joint probability matrices have zero diagonals by convention. In `joint_counts`, `0` means ground and `1` means excited.

Array positions follow the selected residue order and start at zero in Python. Exported `res_index` values start at one. These positions may differ from PDB residue numbers. Use the original structure to map results back to chain identifiers and residue numbers.

```python
samples = outputs.results["Z_samples"]
energies = outputs.results["E"]

# Each row is a shot; each column is a selected residue.
print(samples.shape, energies.shape)

# Sampled excitation frequency for the first selected residue.
print((samples[:, 0] == -1).mean())

# Excited source at position 1, target at position 0.
print(outputs.dei_matrix[0, 1])
```

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

Pass the same configuration, molecule name, and chain used in `run_analysis`. The parent directory of `output_path` must exist. The example above writes into the directory already created for the analysis.

## Work with samples and energies

For numerical results without the full set of report tables and figures, use `run_pipeline`:

```python
from qprotsim.pipeline import run_pipeline

# Reuse cfg from the earlier example.
results = run_pipeline(
    cfg=cfg,
    molecule_name="Trp Cage",
    chain="A",
    plot_energy=False,
)
samples = results["Z_samples"]
energies = results["E"]
```

| Result key | Shape | Meaning |
| --- | --- | --- |
| `Z_samples` | `(S, N)` | Residue states: `+1` is ground, `-1` is excited. |
| `E` | `(S,)` | Additive model energy for each shot, in kJ mol⁻¹. |
| `E0`, `E1` | `(N,)` | Ground and excited reference energies, in kJ mol⁻¹. |
| `p_exc` | Length `N` | Starting excitation probability for each residue. |
| `p_ground` | `(N,)` | Ground-state frequency after coupling and sampling. |
| `resnames` | Length `N` | Three-letter residue names in sample-column order. |
| `alpha_SASA` | `(N,)` | Relative solvent exposure, capped at `0.999`. |

The final excitation frequency is `1 - results["p_ground"]`. It can differ from `results["p_exc"]` because the latter describes the starting state.

For a detailed examination of the model, `pairs`, `kijs`, `kij_terms`, and `entanglement_schedule` contain the residue pairs, coupling strengths, coupling components, and ordered circuit parameters. See the [model definitions](theory.md#mathematical-definitions) for their scientific meaning.

## Build and inspect a Qiskit circuit

Install the optional dependencies with `python -m pip install ".[qiskit]"`. To use Qiskit for the workflows above, create the configuration with:

```python
cfg = Config(
    pdb_path="input/1L2Y.pdb",
    shots=10_000,
    backend="qiskit",
    qiskit_backend="aer_simulator",
    backend_options={"method": "statevector", "device": "CPU"},
    seed=42,
)
```

To inspect the native `QuantumCircuit`, reuse the numerical inputs from `run_pipeline`:

```python
from qprotsim.circuit import build_qiskit_circuit

circuit = build_qiskit_circuit(
    results["p_exc"],
    results["pairs"],
    results["kijs"],
    results["E0"],
    results["E1"],
    schedule=results["entanglement_schedule"],
    measure=False,
)
print(circuit.draw(output="text"))
```

The builder uses the same ordered RY, CRY, and YY rotation gates as the PennyLane circuit. Omit `measure=False` to include computational-basis measurements. It returns a circuit without executing it.

In both execution paths, `sample_bits` returns an array of shape `(shots, n_residues)` containing Pauli-Z outcomes: `|0>` maps to `+1` and `|1>` maps to `-1`. Qiskit's binary outcomes are converted with `z = 1 - 2*b`. A Qiskit bit string displays qubit 0 on the right; QProtSim places that qubit in sample column 0. See [IBM's bit-ordering guide](https://quantum.cloud.ibm.com/docs/en/qiskit/bit-ordering). The labels `|+>` and `|->` are Pauli-X eigenstates, rather than names for these Pauli-Z outcomes.

## Calculate information measures

You can calculate these measures from the samples returned by either workflow. Convert the Pauli-Z samples to explicit binary excitation states for MI, DEI, and residue information. Keep the original Pauli-Z samples for the joint probability function.

```python
from qprotsim.info_theory import (
    directional_excitation_information,
    joint_state_matrices,
    mutual_information_matrix,
    residue_information_gain,
)

binary_excited = (samples == -1).astype("uint8")
mi = mutual_information_matrix(binary_excited)
dei = directional_excitation_information(binary_excited)
info = residue_information_gain(binary_excited, energies, bins=100)
p00, p01, p10, p11 = joint_state_matrices(samples)
```

| Function | Result |
| --- | --- |
| `mutual_information_matrix(samples)` | Symmetric `N × N` dependence scores in bits. |
| `directional_excitation_information(samples, threshold=0.0)` | `N × N` conditional KL scores in bits, indexed `[target, source]`. A nonnegative threshold is subtracted, then values are clamped at zero. |
| `residue_information_gain(samples, energies, bins=10000)` | Data frame with `res_index`, `KL_ground`, and `KL_exc`. Each score column is normalized separately across residues. |
| `joint_state_matrices(z_samples)` | Four joint probability matrices for `00`, `01`, `10`, and `11`; requires Pauli-Z `+1`/`-1` samples. |

The first three functions also accept Pauli-Z samples, but an array containing only `+1` values has an ambiguous encoding. The explicit binary conversion above avoids that ambiguity. Binary states use `0=ground` and `1=excited`.

Check for missing scores before interpreting rankings:

```python
print(info[["KL_ground", "KL_exc"]].isna().any())
```

`NaN` values can occur when a residue has no observations in one state or when the total raw divergence is zero. A DEI source column is zero when the source is never observed excited. Review state counts as well as scores, especially for rare excitations. The [theory guide](theory.md) explains the definitions and limitations.

## Examine groups of connected residues

For advanced analysis, `analyse_motif_allostery` groups associations over small structural motifs, including residue pairs and triads. Reuse the structure, chain, and residue selection that produced the samples:

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

`min_support` is the minimum number of excited-source observations. The returned dictionary includes `residue_allostery_matrix`, `directional_region_coupling_matrix`, `affinity_matrix`, `cluster_labels`, and `saved_paths`.

Motif scores use a signed excitation contribution in nats, while MI and DEI use their own definitions in bits. Read motif scores using the [motif association definition](theory.md#motif-association). These associations help form hypotheses about connected residues; establishing biological allostery requires independent evidence.
