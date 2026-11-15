# Key concepts

You can run your first QProtSim analysis without knowing quantum computing or information theory. This page introduces the terms you will meet in the guides and reports. The [scientific model](theory.md) provides the equations when you want more depth.

## What does QProtSim model?

A protein is a chain of amino acids folded into a three-dimensional structure. An amino acid within that chain is called a **residue**. A protein may contain several chains, often labelled A, B, and so on.

QProtSim starts from a structure saved in a **PDB file**. You can supply this file yourself or use a Protein Data Bank identifier, such as `1L2Y`, to download it. You choose which chain to analyse.

For each selected residue, the model uses two possible solvation states: **ground-like** and **excited**. Solvation describes how a molecule interacts with its surrounding solvent, usually water. These are simplified model labels, rather than direct observations of a residue folding, moving, or becoming unstable.

## How does the analysis work?

1. **Measure exposure.** The structure shows how much of each residue is accessible to solvent.
2. **Set starting probabilities.** Exposure and residue chemistry set how likely each residue is to start in the excited model state.
3. **Connect residues.** Structural contacts provide coupling scores that influence how residue states depend on one another.
4. **Sample combinations of states.** A simulated quantum circuit produces many configurations, with one state for each residue in every sample.
5. **Summarise the samples.** QProtSim calculates model energies and statistical relationships, then writes the report, tables, and figures.

The input coordinates stay fixed. The sampled configurations describe combinations of model states, rather than a sequence of protein movements through time.

## Why a quantum circuit?

The model represents each residue with one **qubit**, the two-state unit used in a quantum circuit. **Gates** are operations that prepare the probabilities and introduce dependencies between residues. QProtSim uses PennyLane, a Python library, to simulate these operations on a CPU or a supported GPU.

A **shot** is one measurement of the circuit. For example, `shots: 10000` requests 10,000 sampled configurations. The collection is called an **ensemble**. More shots can improve estimates of probabilities, although rare states may still be poorly sampled.

Each extra residue makes a standard statevector simulation much larger. Start with the small protein in the [first-analysis guide](getting-started.md), even if you plan to use a GPU later.

## Terms in the settings and results

| Term | Plain-language meaning | How to read it |
| --- | --- | --- |
| Solvent-accessible surface area (SASA) | How much of a residue's surface the solvent can reach. | `SASA_norm` expresses exposure relative to a reference maximum for that residue type. |
| Excitation propensity, `p_exc` | The starting probability of the excited state. | This is set before the coupling gates act. |
| Ground-state frequency, `p_ground` | The fraction of samples where a residue is ground-like. | The final excited-state frequency is `1 - p_ground`. |
| Coupling, `Kij` | A model score for a structural connection between two residues. | Stronger coupling influences the circuit; it does not by itself prove biological communication. |
| Energy distribution | How often different model energies occur in the samples. | A histogram groups nearby energies into bins. |
| Residue information score, `KL_ground` or `KL_exc` | How much the energy distribution differs when you select samples with a particular residue state. | Scores are normalised separately within each column. They are relative contributions within that analysis. |
| Mutual information (MI) | How much knowing one residue's sampled state tells you about another's. | It measures statistical dependence and is symmetric between the pair. |
| Directional excitation information (DEI) | How a target's state distribution changes when a source is excited. | Rows are targets; columns are sources. Directional conditioning does not establish cause and effect. |
| Motif | A small group of residues defined from the structure. | Motif analysis summarises associations across groups, which can overlap. |
| Allostery | A biological effect where a change at one site influences another site. | QProtSim's association scores can suggest sites to investigate; confirming allostery needs additional evidence. |

## A small reading example

Suppose a residue has `p_exc = 0.20` and `p_ground = 0.65`. Its starting excitation probability was 20%, while it was excited in 35% of the measured samples after the coupling gates acted. Those values describe different stages of the analysis, so they need not add up to one.

If that residue also has a high information score, its state is strongly associated with changes in the model's energy distribution relative to other residues in the same run. This gives you a candidate for closer study. It does not establish the residue's experimental importance.

Continue with [installation](installation.md) and [your first analysis](getting-started.md), or explore [how to interpret the outputs](outputs.md).
