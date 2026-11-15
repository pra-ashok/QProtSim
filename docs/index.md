---
title: QProtSim
sidebar_title: Overview
---

<div class="hero" markdown>

# Protein ensembles through quantum circuits

QProtSim connects solvent exposure, structural contacts, and a two-state residue model through a PennyLane quantum circuit. Use it to sample residue configurations, examine energy distributions and information-sharing patterns, and generate a PDF report with reusable tables and figures.

[Run your first analysis](getting-started.md){ .button }
[Read the paper](paper.md){ .button .outline }

</div>

**Authors:** [Pratik Patil](https://pra-ashok.github.io/CV/) · [Bhushan Bonde](https://scholar.google.com/citations?hl=en&user=O3hkOEUAAAAJ) · [Bhaskar Choubey](https://scholar.google.com/citations?user=b88L0aEAAAAJ&hl=en)

## Start with a structure

```bash
python -m pip install ".[cpu]"
qprotsim analyze --config config.yaml
```

Run these commands from a checkout of the repository. Python imports use the lowercase package name `qprotsim`.

<div class="grid cards" markdown>

- **Run an analysis**

    Install a backend and analyse a small protein from a local PDB file or an RCSB identifier.

    [Installation](installation.md) · [First analysis](getting-started.md)

- **Understand the model**

    Follow the equations from solvent exposure through coupling and circuit gates to energy and information metrics.

    [Theory and diagrams](theory.md)

- **Use Python directly**

    Work with `run_analysis`, access sampled arrays and matrices, or run the numerical pipeline independently.

    [Python API](api.md)

- **Read the results**

    Identify CSV schemas, state encodings, matrix orientation, and the quantities shown in reports.

    [Outputs and interpretation](outputs.md)

</div>

## From a PDB structure to a report

Each selected residue becomes one qubit. FreeSASA exposure estimates set the initial excitation probabilities. Structural interactions supply the coupling graph, and an ordered circuit transforms the initial product state. Measurements provide the ensemble used by the report workflow.

<figure class="paper-figure" markdown>

[![Original paper workflow illustrating the mapping from Chignolin residues to qubits, controlled rotations, measurement, and energy analysis.](assets/paper/figure-1.png){ loading=lazy }](assets/paper/figure-1.png)

<figcaption>Figure 2.1 from Patil, Bonde, and Choubey (2026), illustrated with Chignolin (1UAO). Reproduced unchanged. Select the image for full resolution. <a href="https://arxiv.org/abs/2608.05491">Paper</a> · <a href="https://creativecommons.org/licenses/by-nc-nd/4.0/">CC BY-NC-ND 4.0</a></figcaption>
</figure>

The numerical pipeline returns state samples with shape `(shots, residues)` and one additive energy per shot. `run_analysis` adds information metrics, CSV exports, and figures. `generate_report` assembles those outputs into a PDF. The YAML wrapper and CLI perform all three stages.

## Scientific scope

The model uses a fixed input structure, empirical residue energy tables, and a prescribed gate schedule. Its outputs support hypotheses about residue coupling and sensitivity. Calibration and comparison with experimental or simulation data are needed to establish physical predictive accuracy. The [scientific model](theory.md) records the implemented equations and assumptions.

## Software and attribution

<div class="publication" markdown>

**A Quantum Circuit Framework for Protein Ensemble-Level Energetics**

Pratik Patil, Bhushan Bonde, and Bhaskar Choubey. 2026. arXiv:2608.05491.

[Read on arXiv](https://arxiv.org/abs/2608.05491) · [Research overview and paper figures](paper.md) · [Download BibTeX](assets/citations/qprotsim.bib)

</div>

When using QProtSim in research, cite the paper and record the exact software version or commit. These pages describe QProtSim 0.6.1 and the accompanying source tree. The research page includes a summary of the alphaXiv AI overview and original manuscript figures with attribution. [alphaXiv overview](https://www.alphaxiv.org/abs/2608.05491)

QProtSim software is distributed under the MIT License. [Source repository](https://github.com/pra-ashok/QProtSim) · [Software license](https://github.com/pra-ashok/QProtSim/blob/main/LICENSE)
