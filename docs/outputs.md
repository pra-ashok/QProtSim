# Outputs and interpretation

After an analysis, open the PDF in your chosen output folder. It brings together the run settings, plots, and summary tables. CSV files contain the numerical data so you can inspect them in a spreadsheet or use them in Python.

If the terminology is new, start with [key concepts](concepts.md). For equations and units, see the [scientific model](theory.md).

## Read your first report

1. **Check the run settings.** Confirm the input structure, selected chain, and number of samples.
2. **Look at the energy distribution.** The histogram shows how often different model energies were sampled. Its width describes variation within this model and run.
3. **Read the residue statistics.** Compare solvent exposure with the sampled ground-state frequencies.
4. **Explore residue information scores.** These show how selecting a residue state changes the sampled energy distribution.
5. **Inspect pair and motif plots.** Use these to identify associations that may deserve further study.

Read plot labels such as “stability” and “sensitivity” in terms of the quantities being plotted. They summarise model state frequencies and information scores; experimental stability or a biological response needs independent validation.

## Find the main files

For `molecule_name: Trp Cage`, the filename prefix is `trp_cage`. In the table below, replace `<molecule>` with that prefix. The default report name is `report.pdf`; the first-analysis example changes it to `trp_cage_report.pdf`.

| File | What it contains |
| --- | --- |
| `report.pdf` | Run settings, figures, and summary tables. |
| `<molecule>_residue_stats.csv` | Residue labels, normalised solvent exposure, starting excitation probability, and sampled ground-state frequency. |
| `<molecule>_residue_information_gain.csv` | Relative energy-information scores, `KL_ground` and `KL_exc`, for each residue. |
| `<molecule>_z_samples.csv` | Every sampled configuration, with one residue per column and a leading `shot` column. |
| `<molecule>_kij_interaction_matrix.csv` | Structural coupling scores between residue pairs. |
| `<molecule>_kij_term_<term>_matrix.csv` | Separate contributions to the structural coupling scores. |
| `<molecule>_mutual_information.csv` | Pairwise statistical dependence, measured in bits. |
| `<molecule>_directional_excitation_information.csv` | Changes in target state distributions when source residues are excited, measured in bits. |
| `<molecule>_entanglement_parameters.csv` | Circuit operations and their parameters, for detailed analysis of how the model was sampled. |

## Understand the residue statistics

| Column | Meaning |
| --- | --- |
| `res_index` | Consecutive residue position in the selected analysis, starting at 1. |
| `resname` | Three-letter amino-acid name, such as `ALA`. |
| `SASA_norm` | Solvent-accessible area divided by a reference maximum for that residue type, capped at 0.999. |
| `p_exc` | Excitation probability used to prepare the circuit, before coupling gates act. |
| `p_ground` | Fraction of measured samples in the ground-like state, after coupling gates act. |

The final sampled excitation frequency is **`1 - p_ground`**. For example, `p_ground = 0.8` means the residue was ground-like in 80% of the samples and excited in 20%. It need not equal the starting `p_exc` because the coupling gates can change state probabilities.

**Check residue numbering before mapping results to a structure.** Labels such as `1-ASN` and `2-LEU` use analysis order. They can differ from the residue numbers in the PDB file, which may have gaps or restart between chains. The circuit table's `source_qubit` and `target_qubit` positions start at zero.

## Understand the information scores

### Residue energy-information scores

`KL_ground` and `KL_exc` compare the overall energy distribution with the distribution obtained by selecting samples where a residue is ground-like or excited. A larger value means that selecting that state makes a larger relative difference to the energy distribution.

Each column is normalised separately across residues. A score of `0.1` represents 10% of that column's total conditional information contribution. These are dimensionless relative scores, so a ranking is most directly interpreted within the same analysis. A high score alone does not predict the effect of a mutation.

A missing value shown as `NaN` means the score could not be estimated. This can happen when a residue is only observed in one state or when the column cannot be normalised. Check the sample counts before interpreting rankings; the [troubleshooting guide](troubleshooting.md) explains what to do.

### Pairwise matrices

A matrix is a table with residues on both axes. Read the labels to identify the pair:

| Matrix | How to read it |
| --- | --- |
| Structural coupling, `Kij` | Larger values indicate a stronger model connection derived from the structure. |
| Mutual information, MI | Larger values indicate greater dependence between the two sampled residue states. The table is symmetric: swapping the pair gives the same score. |
| Directional excitation information, DEI | Read **target rows** and **source columns**. An entry describes the target's state distribution when the source is excited, compared with the target's overall distribution. Swapping the pair can give a different score. |

MI and DEI use bits, a unit of information. Their diagonals are set to zero by convention. If a source is never observed excited, its DEI column is also zero, so check the source's sample count before treating a zero as evidence of weak association.

A directional score describes conditional association. Establishing cause and effect or biological communication needs additional evidence. A structural coupling score and an information score also describe different parts of the model: input connections and sampled relationships, respectively.

## Use the data in Python

Each row of the sample CSV is one configuration. **`+1` means ground-like and `-1` means excited.** Rows are repeated circuit measurements, not time steps in a trajectory.

```python
import pandas as pd

samples = pd.read_csv(
    "outputs/1l2y/trp_cage_z_samples.csv", index_col="shot"
)
excitation_frequency = (samples == -1).mean(axis=0)
print(excitation_frequency.head())
```

Read a DEI entry using its residue labels:

```python
dei = pd.read_csv(
    "outputs/1l2y/trp_cage_directional_excitation_information.csv",
    index_col=0,
)
target = dei.index[0]
source = dei.columns[1]
print(dei.loc[target, source])
```

The [Python API](api.md) explains how to access arrays and tables directly during an analysis.

## Figures and motif analysis

Figure names include numbers, such as `fig03_energy_distribution.png`. Available structure images and optional analyses affect the numbering, so use the actual report captions and filenames. Most plots also have a PDF version. A structure illustration may show more of the deposited structure than the chain selected for numerical analysis; confirm the selection using the run settings and residue labels.

The joint-state plots show how often each pair is ground–ground, ground–excited, excited–ground, or excited–excited. Each pair's four off-diagonal probabilities sum to one; diagonal entries are set to zero.

A **motif** is a small structural group, such as a pair or a three-residue group. When motif analysis succeeds, additional files describe associations between these groups:

| File or group | What it contains |
| --- | --- |
| `motif_allostery_residue_allostery_matrix.csv` | Signed excitation-association scores between residues, in nats. |
| `motif_allostery_region_to_residue_matrix.csv` | Associations from source groups to target residues. |
| `motif_allostery_membership_matrix.csv` and `motif_allostery_residue_memberships.csv` | Which residues belong to each group and their membership weights. |
| `motif_allostery_region_summary.csv` | Group definitions, sample support, activity, and cluster assignment. |
| `motif_allostery_directional_region_coupling.csv` | Directional scores between groups. |
| `motif_allostery_affinity_matrix.csv` | Symmetric, scaled group-association scores used for clustering. |
| `motif_allostery_cluster_summary.csv` | Group counts and sample support by cluster. |
| `motif_allostery_affinity_matrix_heatmap.png` | Group-association heatmap, with a companion PDF. |
| `allostery_summary.txt` | A text summary of groups and prominent connections. |

Separate pair and three-residue summaries use `motif_allostery_dimer_*` and `motif_allostery_trimer_*` filenames. Some groups overlap. A motif-stage warning can leave these optional outputs unavailable while the main analysis continues.

Motif residue scores can be positive or negative and use **nats**, an information unit based on natural logarithms. Their definition differs from the main DEI score in bits. See [motif association](theory.md#motif-association) before comparing these outputs. The filenames use “allostery” for exploratory associations; confirming an allosteric effect requires independent evidence.

## Keep results from repeated runs

Use a different output folder for each run. Reusing a folder can overwrite earlier tables and figures. Keep the input structure, configuration, software version, and backend information with the results so you can reproduce and compare the analyses.
