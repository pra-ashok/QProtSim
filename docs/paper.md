# The research behind QProtSim

This page links the research paper, shows its original figures, and explains how to cite it. For a gentle introduction to the software, read [key concepts](concepts.md). For the model's definitions and assumptions, read the [scientific model](theory.md).

<div class="publication" markdown>

## A Quantum Circuit Framework for Protein Ensemble-Level Energetics

**Authors:** [Pratik Patil](https://pra-ashok.github.io/CV/) · [Bhushan Bonde](https://scholar.google.com/citations?hl=en&user=O3hkOEUAAAAJ) · [Bhaskar Choubey](https://scholar.google.com/citations?user=b88L0aEAAAAJ&hl=en)

arXiv:2608.05491 · 2026 · Version 1, submitted 6 August 2026

[Read the paper](https://arxiv.org/abs/2608.05491){ .button }
[Explore on alphaXiv](https://www.alphaxiv.org/abs/2608.05491){ .button .outline }
[Download BibTeX](assets/citations/qprotsim.bib){ .button .outline download }

</div>

<figure class="qps-cover" markdown>

[![Original conference overview of quantum encoding of proteins, with residue interaction graphs, circuit implementation, an energy landscape, and entanglement illustrations.](assets/branding/qprotsim-cover.png){ .qps-cover-image width=989 height=551 }](assets/branding/qprotsim-cover.png)

<figcaption>“Quantum Encoding of Protein for exploring Energy Landscape,” Pratik Patil et al. Presented at the MGMS conference in Oxford, UK, on 29 November 2025. Reproduced unchanged from the supplied figure. Select the image to view it at full size.</figcaption>
</figure>

## What is the main idea?

The framework represents each amino-acid residue with one qubit and two possible solvation states. Solvent exposure and residue chemistry set the starting probabilities. Structural contacts determine circuit operations that introduce dependencies between residues. Repeated measurements produce many combinations of residue states, called an ensemble.

QProtSim uses that ensemble to explore model energy distributions, residue sensitivity, and statistical associations between residues. The paper illustrates the framework with protein examples, including Trp-cage (1L2Y) and 9GDL. The figures below introduce the workflow and show example energy distributions.

The practical use is to identify patterns and candidate residues for further study. Predicting experimental stability, mutation effects, or biological communication requires validation beyond these model scores. Read the paper alongside the [model assumptions](theory.md#assumptions-and-interpretation) when assessing those possibilities.

## Figures from the paper

The three figures below are reproduced from Pratik Patil, Bhushan Bonde, and Bhaskar Choubey, arXiv:2608.05491v1, under the publication's CC BY-NC-ND 4.0 licence. [Original paper](https://arxiv.org/html/2608.05491v1) · [Figure licence](https://creativecommons.org/licenses/by-nc-nd/4.0/)

Select a figure to open its original resolution. [Figure attribution and source record](assets/paper/provenance.json)

<figure class="paper-figure" markdown>

[![The paper's residue-to-qubit workflow showing propensity calculation, quantum state preparation, structural interactions, and measurement.](assets/paper/figure-1.png){ loading=lazy }](assets/paper/figure-1.png)

<figcaption>Paper Figure 2.1. Residue-to-qubit workflow, illustrated with Chignolin (1UAO). <a href="https://arxiv.org/html/2608.05491v1">Source</a></figcaption>
</figure>

<figure class="paper-figure" markdown>

[![Original paper figure illustrating the circuit's ensemble sampling and its relationship to protein energetic states.](assets/paper/figure-2.png){ loading=lazy }](assets/paper/figure-2.png)

<figcaption>Paper Figure 4.2. A conceptual view of repeated circuit measurement and ensemble sampling. <a href="https://arxiv.org/html/2608.05491v1">Source</a></figcaption>
</figure>

<figure class="paper-figure" markdown>

[![Paper energy histograms comparing circuit samples for the 1L2Y and 9GDL protein structures.](assets/paper/figure-3.png){ loading=lazy }](assets/paper/figure-3.png)

<figcaption>Paper Figure 3.2. Energy histograms for 1L2Y and 9GDL from one million samples each. The histograms use different bin widths. <a href="https://arxiv.org/html/2608.05491v1">Source</a></figcaption>
</figure>

These histograms describe sampled model energies. Each protein has its own residue composition, and the plots use different bin widths, so compare their scales and captions before interpreting their shapes. The [scientific model](theory.md) explains how energies and information scores are calculated. [Paper methods and results](https://arxiv.org/html/2608.05491v1)

## Cite the paper

Pratik Patil, Bhushan Bonde, and Bhaskar Choubey. *A Quantum Circuit Framework for Protein Ensemble-Level Energetics*. arXiv:2608.05491, 2026. DOI 10.48550/arXiv.2608.05491. [Paper and citation record](https://arxiv.org/abs/2608.05491)

Include the QProtSim version or commit and the input configuration when reporting software results. [Download the BibTeX entry](assets/citations/qprotsim.bib)

## Credits

**Authors:** [Pratik Patil](https://pra-ashok.github.io/CV/) · [Bhushan Bonde](https://scholar.google.com/citations?hl=en&user=O3hkOEUAAAAJ) · [Bhaskar Choubey](https://scholar.google.com/citations?user=b88L0aEAAAAJ&hl=en)

Paper figures retain their publication attribution and licence, recorded separately from the package's MIT licence. [Publication record](https://arxiv.org/abs/2608.05491)
