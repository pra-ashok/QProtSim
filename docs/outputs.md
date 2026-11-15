# Outputs and interpretation

The CLI and `run_analysis_from_yaml` produce a PDF plus numerical tables and figures. `run_analysis` produces the tables and figures and returns them through `AnalysisOutputs`. Call `generate_report` to build its PDF.

## Numerical files

The filename prefix below is `molecule_name.lower().replace(" ", "_")`. For `Trp Cage`, `<molecule>` is `trp_cage`.

| Filename | Contents |
| --- | --- |
| `report.pdf` | Default PDF report, including run settings, figures, and summary tables. |
| `<molecule>_residue_stats.csv` | `res_index`, `resname`, normalized SASA, initial `p_exc`, and sampled `p_ground`. |
| `<molecule>_entanglement_parameters.csv` | One row per directed transfer layer, including residue labels, gate type, coupling, drive, probabilities, discount, and gate angles. |
| `<molecule>_residue_information_gain.csv` | `res_index`, `KL_ground`, and `KL_exc`. |
| `<molecule>_z_samples.csv` | Pauli-Z sample matrix with a leading `shot` index column. |
| `<molecule>_kij_interaction_matrix.csv` | Symmetric structural coupling matrix. |
| `<molecule>_kij_term_<term>_matrix.csv` | Symmetric matrices for available scalar terms of the coupling decomposition. |
| `<molecule>_directional_excitation_information.csv` | DEI matrix in bits, with target residues in rows and source residues in columns. |
| `<molecule>_mutual_information.csv` | Symmetric pairwise mutual information matrix in bits. |

Residue labels such as `1-ASN` and `2-LEU` use consecutive, one-based positions in the selected structure. `source_qubit` and `target_qubit` in the circuit table use zero-based positions. Input PDB residue numbering can contain gaps, insertion codes, or chain restarts, so retain the selected structure order when mapping results back to a PDB.

### Read sampled states

```python
import pandas as pd

sample_table = pd.read_csv(
    "outputs/trp_cage/trp_cage_z_samples.csv", index_col="shot"
)
z_samples = sample_table.to_numpy()
excited = (z_samples == -1)
sampled_excitation_frequency = excited.mean(axis=0)
```

Each row is one measured configuration. `+1` is ground and `-1` is excited. `p_ground` is the measured fraction of `+1` outcomes. The residue statistics column `p_exc` contains the initial propensity used to prepare the circuit. Coupling gates can change the measured excitation frequency.

### Read matrices

```python
dei = pd.read_csv(
    "outputs/trp_cage/trp_cage_directional_excitation_information.csv",
    index_col=0,
)
target_label = dei.index[0]
source_label = dei.columns[1]
score_bits = dei.loc[target_label, source_label]
```

The entry `[i, j]` of DEI measures the shift in residue `i`'s state distribution conditional on residue `j` being excited. A source never observed excited has a zero DEI column. MI uses both ground and excited source states and is symmetric. Both matrices have zero diagonals by convention.

`KL_ground` and `KL_exc` measure changes in the ensemble energy histogram conditional on a residue's state. Each output column is divided by its sum over residues, yielding dimensionless relative contributions. A value of `0.1` therefore represents ten percent of that column's total conditional KL contribution. The raw histogram divergences use natural logarithms. Sparse or deterministic state samples can yield `NaN` in these columns. The [theory page](theory.md) gives the definitions.

## Figures

Numbered figures use names such as `fig03_energy_distribution.png`. Successful structure images are allocated the first numbers, and available coupling components add further figures. Figure numbers therefore depend on which artifacts were generated.

The report includes structure views when available, energy histograms and density, residue information plots, stability and sensitivity comparisons, coupling heatmaps, DEI heatmaps and clustering, MI, joint-state probabilities, and sampled ground-state comparisons. Most plot exports include a companion PDF.

`AnalysisOutputs.plot_paths` provides the exact generated figure paths. The keys `stability_vs_info` and `stability_vs_sensitivity_quadrants` refer to the same plot. `joint_counts` contains normalized joint probabilities in memory. The main report draws them as four heatmaps.

The RCSB image represents the deposited entry and the PyMOL view renders the loaded structure. Use the run's chain selection and residue labels to establish which residues contribute to the numerical analysis.

## Motif allostery files

The analysis attempts an additional motif workflow using peptide dimers, peptide trimers, and noncovalent triads. The following files are present when that stage succeeds.

| Filename | Contents |
| --- | --- |
| `motif_allostery_residue_allostery_matrix.csv` | Signed residue excitation contribution in nats, with a minimum excited-source support filter. |
| `motif_allostery_region_to_residue_matrix.csv` | Directional scores aggregated from motif source regions to target residues. |
| `motif_allostery_membership_matrix.csv` | Residue membership weights for each motif region. |
| `motif_allostery_region_summary.csv` | Region definitions, support, activity, and cluster assignment. |
| `motif_allostery_residue_memberships.csv` | Residue-to-region membership records. |
| `motif_allostery_directional_region_coupling.csv` | Directional aggregation between motif regions. |
| `motif_allostery_affinity_matrix.csv` | Symmetric, scaled motif affinity used for clustering. |
| `motif_allostery_cluster_summary.csv` | Region counts and support summaries by cluster. |
| `motif_allostery_affinity_matrix_heatmap.png` | Combined motif affinity heatmap, with a companion PDF. |
| `allostery_summary.txt` | Textual summary of motifs, support, and prominent connections. |

Motif residue scores use `q_exc * ln(q_exc / p_exc)`, where `q_exc` is the target excitation probability conditional on source excitation. Values can be signed. Region aggregation and symmetric affinity scaling follow this calculation. The [theory page](theory.md) explains the relationship to the full Bernoulli DEI used in the main report.

Available dimer and trimer groups also receive `motif_allostery_dimer_*` and `motif_allostery_trimer_*` affinity matrices, directional coupling matrices, region summaries, and heatmaps. The report includes their affinity plots when available. A failed motif stage emits a warning while the core analysis continues.

## Repeated runs

Artifacts use stable filenames within the chosen output directory. The timestamp is recorded inside the report, while repeated runs can overwrite earlier CSV and figure files. Use a distinct output directory for each run when retaining results, and keep the input structure, configuration, package version, and selected backend alongside the outputs.

The numerical exports and artifact names are defined in [`report_analysis.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/report_analysis.py), [`pipeline.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/pipeline.py), and [`info_theory.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/info_theory.py).
