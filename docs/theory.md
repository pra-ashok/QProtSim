# Theory and mathematical model

QProtSim maps the residues in a PDB structure to two-state variables, prepares a quantum circuit from solvent exposure and structural contacts, and analyzes the sampled ensemble. The equations below describe the executable implementation in version 0.6.1. The circuit defines a model distribution over residue states. Its interpretation depends on the chosen energy table, contact weights, and circuit schedule.

The accompanying research page provides the paper citation, an alphaXiv overview summary, and original manuscript figures. [Research paper and figures](paper.md)

<figure markdown>

![Workflow from protein coordinates through solvent exposure and contact couplings to a sampled circuit ensemble and analysis outputs.](assets/diagrams/workflow.svg)

<figcaption>Figure 1. The structure supplies local exposure values and a residue contact graph. These become circuit parameters. Information measures summarize the resulting samples.</figcaption>

</figure>

## State convention and notation

For residue \(i\), the computational state \(|0\rangle\) represents ground-like solvation and \(|1\rangle\) represents excited solvation. Measurement returns a Pauli-Z eigenvalue \(Z_{ti}\), where \(t\) indexes a shot.

\[
R_{ti}=\frac{1-Z_{ti}}{2},\qquad
Z_{ti}\in\{+1,-1\},\qquad R_{ti}\in\{0,1\}.
\]

Thus \(Z=+1\) means ground and \(Z=-1\) means excited. The sampled array has shape \((S,N)\) for \(S\) shots and \(N\) selected residues. Residue indices in matrices follow the selected MDAnalysis residue order. CSV `res_index` values count from one. Input PDB residue numbering can contain offsets, gaps, and insertion codes.

| Symbol | Meaning | Numerical convention |
| --- | --- | --- |
| \(\alpha_i\) | Relative solvent exposure | Dimensionless |
| \(p_i\) | Excitation probability used to initialize the circuit | Dimensionless |
| \(\widehat p_i\) | Excitation probability estimated after circuit sampling | Dimensionless |
| \(\Delta G_i=E_{1i}-E_{0i}\) | Reference energy gap | Treated as kJ mol\(^{-1}\) |
| \(K_{ij}\) | Normalized residue contact coupling | Dimensionless |
| \(\theta_i,\phi_\ell\) | Circuit rotation angles | Radians |

## Solvent exposure and local excitation

FreeSASA calculates the residue solvent-accessible area \(A_i\) from the input coordinates. QProtSim divides this area by the residue-specific theoretical maximum from Tien et al. and caps the ratio. The returned dictionary from `per_residue_sasa` contains this relative exposure, which is also exported as `SASA_norm`. FreeSASA evaluates the full input structure before QProtSim filters the returned residues by chain. [FreeSASA Python calculation guide](https://freesasa.github.io/python/intro.html), [Tien et al., 2013](https://doi.org/10.1371/journal.pone.0080635)

\[
\alpha_i=\min\!\left(0.999,\frac{A_i}{A_{\max}(a_i)}\right).
\]

The code selects the `energies_bulk_octabol` reference table by default. All entries in this table use \(E_{0i}<0\) and \(E_{1i}=0\). The Kyte–Doolittle value \(h_i\) selects an exposure-dependent effective gap \(g_i\). The exact branching in `excitation_propensities` is

\[
g_i=\begin{cases}
\Delta G_i(1-\alpha_i), & h_i>0.5,\\[2pt]
\Delta G_i[1-\max(\alpha_i,0.1)], & -0.5<h_i<0.5,\\[2pt]
\Delta G_i\alpha_i, & \text{otherwise}.
\end{cases}
\]

The two-level Boltzmann calculation then supplies the initial excitation probability

\[
p_i=\frac{e^{-g_i/(RT)}}{1+e^{-g_i/(RT)}}
   =\frac{1}{1+e^{g_i/(RT)}},
\qquad T=298\ \mathrm{K},\qquad
R=8.314462618\times10^{-3}\ \mathrm{kJ\,mol^{-1}\,K^{-1}}.
\]

The gas constant is named `k_B` in `boltzmann_expectation`. These numerical units make it a molar energy calculation. Exposure increases excitation propensity in the hydrophobic branch and decreases it in the hydrophilic branch. The neutral branch imposes a minimum exposure of 0.1. With nonnegative effective gaps, initial probabilities satisfy \(0<p_i\leq0.5\). The current YAML interface uses the fixed temperature above. [Implementation in `circuit.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/circuit.py)

<figure markdown>

![Effective energy gaps and initial excitation probabilities as functions of relative solvent exposure for alanine, glycine, and lysine at 298 kelvin.](assets/diagrams/excitation.svg)

<figcaption>Figure 2. Curves use the implemented default reference values for ALA, GLY, and LYS. They illustrate the initialization rule. Sampled probabilities can change when coupling gates act.</figcaption>

</figure>

## Structure-derived couplings

`build_Kij_pairs` builds an undirected graph of selected residues. Its structural score combines backbone connectivity, disulfide bonds, approximate hydrogen bonding, salt bridges, and heavy-atom packing.

\[
S_{ij}=1.00\kappa^{\mathrm{BB}}_{ij}
      +1.00\kappa^{\mathrm{SS}}_{ij}
      +0.50\kappa^{\mathrm{HB}}_{ij}
      +0.40\kappa^{\mathrm{salt}}_{ij}
      +0.05\kappa^{\mathrm{vdW}}_{ij}.
\]

The backbone term is one for adjacent selected residues with a backbone C–N distance at most 1.8 Å. When those atoms are absent, sequence adjacency supplies the fallback. The disulfide term is one for cysteine-like SG pairs within 2.3 Å. Other terms use bounded distance-based sums.

\[
\kappa^x_{ij}=1-e^{-u^x_{ij}},\qquad
u^x_{ij}=\frac{1}{c_x}
\sum_{(a,b)\in\mathcal C^x_{ij}}
\left(\frac{r_x}{d_{ab}}\right)^2.
\]

| Term \(x\) | Atom-pair selection \(\mathcal C^x_{ij}\) | Contact cutoff | \(r_x\) | \(c_x\) |
| --- | --- | --- | --- | --- |
| HB | N/O/S pairs | 3.5 Å | 3.5 Å | 1 |
| Salt | Oppositely charged side-chain groups | 4.0 Å | 4.0 Å | 1 |
| vdW | Heavy-atom pairs | 4.0 Å | 3.5 Å | 8 |

The HB score uses a heavy-atom distance approximation. Salt and HB scores can both contribute to the same residue pair. Water-mediated terms exist as helpers, while the active coupling calculation uses the five terms in \(S_{ij}\).

Residue masses \(M_i\) and estimated local stiffnesses add frequency and mass compatibility factors.

\[
k_i=\sum_j S_{ij},\qquad
\omega_i=\sqrt{\frac{\max(k_i,10^{-12})}{\max(M_i,10^{-12})}},
\]

\[
\Phi_{ij}=\exp\!\left[-\frac{(\omega_i-\omega_j)^2}{2\sigma_\omega^2}\right],
\qquad
\Lambda_{ij}=\frac{2\sqrt{M_iM_j}}{M_i+M_j},
\]

\[
\widetilde K_{ij}=S_{ij}\Phi_{ij}\Lambda_{ij},\qquad
K_{ij}=\frac{\widetilde K_{ij}}{\max_{(a,b)}\widetilde K_{ab}}.
\]

The default \(\sigma_\omega\) is the standard deviation of estimated residue frequencies, with a fallback value of one. These frequencies are relative model quantities. The exported `kij_terms` records individual contributions and factors. [Implementation in `coupling_factor.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/coupling_factor.py)

## Circuit initialization and gate schedule

Each qubit starts in \(|0\rangle\). A local rotation encodes its excitation propensity exactly before coupling gates act.

\[
\theta_i=2\arcsin\sqrt{p_i},\qquad
R_Y(\theta_i)|0\rangle
=\sqrt{1-p_i}|0\rangle+\sqrt{p_i}|1\rangle.
\]

The initialized state is \(\bigotimes_i R_Y(\theta_i)|0\rangle\). PennyLane defines these rotations in the computational basis. [PennyLane `RY` reference](https://docs.pennylane.ai/en/stable/code/api/pennylane.RY.html)

The executable transfer schedule clips \(p_i\) to \([\epsilon,1-\epsilon]\), with \(\epsilon=10^{-12}\), then computes

\[
\mu_i=\beta\Delta G_i p_i,\qquad
d_i=\sum_{k\in\mathcal N(i)}K_{ik},\qquad \beta=1.0,
\]

\[
a_{i\rightarrow j}=\max\!\left(0,
\frac{\mu_i-\mu_j}{\mu_i+\mu_j+10^{-5}}\right),
\qquad
\delta_{i\rightarrow j}=\operatorname{clip}_{[0,1]}
\left(\frac{2K_{ij}a_{i\rightarrow j}}{d_i+d_j+\epsilon}\right).
\]

Here \(\beta\) is the schedule's numerical energy scale. The local excitation calculation above uses its own \(1/(RT)\). Schedule entries are sorted by descending \(\delta\). Entry \(\ell\), counting from one, receives

\[
\phi_\ell=\min\!\left(\frac{\pi}{2},
2\arcsin\sqrt{\gamma^\ell\delta_\ell}\right),
\qquad\gamma=0.95.
\]

For \(a_{i\rightarrow j}>0.02\), the circuit applies `CRY(phi, wires=[i, j])` followed by `CRY(-phi, wires=[j, i])`. For weaker scheduled entries it applies an `IsingYY` gate with \(\vartheta=\min(\pi/6,2\arcsin\sqrt{a_{i\rightarrow j}})\). Each schedule entry contributes a sequential gate operation. [PennyLane `CRY` reference](https://docs.pennylane.ai/en/stable/code/api/pennylane.CRY.html), [PennyLane `IsingYY` reference](https://docs.pennylane.ai/en/stable/code/api/pennylane.IsingYY.html)

The exported `deltaP` records a classical schedule score. Actual marginal changes follow the complete ordered unitary circuit. The initial probabilities stay fixed during schedule construction. The `phi_scale` argument accepted by `compute_directed_transfer_schedule` is currently unused. [Implementation in `circuit.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/circuit.py)

<figure markdown>

![A three-residue circuit schematic and a two-residue joint probability example showing symmetric mutual information and directional excitation information.](assets/diagrams/circuit-information.svg)

<figcaption>Figure 3. The circuit panel shows one strong scheduled pair. The probability panel uses a synthetic joint distribution to illustrate matrix orientation and directional conditioning. The example values are mathematical illustrations.</figcaption>

</figure>

## Energies and ground-state probabilities

The sampled configuration energy is the sum of the selected residue state energies.

\[
E_t=\sum_{i=1}^{N}
\left[(1-R_{ti})E_{0i}+R_{ti}E_{1i}\right],\qquad
\widehat p_{\mathrm{ground},i}=\frac{1}{S}\sum_t\mathbf1[Z_{ti}=+1].
\]

With the default reference table, \(\sum_iE_{0i}\leq E_t\leq0\). Couplings influence the sampled configurations through circuit gates. The energy accumulation uses the residue state terms in this equation. The exported `p_exc` stores the initialization propensity \(p_i\), while `p_ground` is measured after coupling gates. The final excitation estimate is \(1-\widehat p_{\mathrm{ground},i}\). [Implementation in `pipeline.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/pipeline.py)

## Residue sensitivity from energy histograms

`residue_information_gain` compares the total-energy histogram with histograms conditioned on each residue's state. Let \(P_b\) be the probability mass in energy bin \(b\), and \(Q^{(s)}_{ib}=P(E\in b\mid R_i=s)\). Its raw divergences use natural logarithms.

\[
D_i^{(s)}=\sum_b Q^{(s)}_{ib}
\ln\frac{Q^{(s)}_{ib}}{P_b},\qquad
\mathrm{KL}_{i,s}=\frac{D_i^{(s)}}{\sum_jD_j^{(s)}}.
\]

The implementation adds \(10^{-12}\) to the global histogram and \(10^{-5}\) to each conditional histogram. SciPy normalizes these vectors during divergence calculation. The returned `KL_ground` and `KL_exc` are dimensionless contributions normalized independently across residues. Their magnitudes depend on the residue set, histogram discretization, and sampling. [SciPy `entropy` reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html)

The report workflow currently calls this function with its default 10,000 energy bins. The `run_analysis` `bins` parameter controls the displayed energy plots. A residue state with zero sampled support produces an undefined conditional histogram. A zero sum of raw divergences also leaves the normalization undefined, so these situations can yield `NaN` values.

## Mutual information and directional excitation information

Mutual information summarizes pairwise dependence over all four joint states. The implementation estimates joint probabilities from samples, adds a numerical stabilizer, and reports bits.

\[
I(R_i;R_j)=\sum_{a,b\in\{0,1\}}
P(R_i=a,R_j=b)\log_2
\frac{P(R_i=a,R_j=b)}{P(R_i=a)P(R_j=b)}.
\]

The MI matrix is symmetric and has a zero diagonal by convention. `joint_state_matrices` returns the four empirical joint-probability matrices \(P_{00},P_{01},P_{10},P_{11}\), also with zero diagonals. Off diagonal, they sum to one for each pair.

Directional excitation information conditions on the source residue being excited. Matrix entry \(D_{ij}\) uses **source column \(j\)** and **target row \(i\)**.

\[
D_{ij}=\sum_{a\in\{0,1\}}
P(R_i=a\mid R_j=1)\log_2
\frac{P(R_i=a\mid R_j=1)}{P(R_i=a)}.
\]

Both target states contribute to this Bernoulli KL divergence. It is nonnegative and can differ between directions. The function uses clipped probabilities with \(\epsilon=10^{-12}\), sets the diagonal to zero, and sets a source column to zero when its excited state has no observations. A positive `threshold` returns \(\max(D_{ij}-\texttt{threshold},0)\). Scores can be large when conditioning on a rarely observed source excitation, so source support matters when comparing residues. [Implementation in `info_theory.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/info_theory.py)

The ideal distribution-level identity connects the two measures through an average over source states,

\[
I(R_i;R_j)=\sum_{b\in\{0,1\}}P(R_j=b)
D_{\mathrm{KL},2}\!\left(P(R_i\mid R_j=b)\,\Vert\,P(R_i)\right).
\]

DEI takes the excited-source conditional divergence directly. Mutual information weights both source conditions by their occurrence probabilities.

## Motif allostery

The report also aggregates residue associations over structural motifs such as peptide dimers and triads. The motif implementation uses an excitation-only signed contribution in natural-log units,

\[
B_{ij}=q_{ij}\ln\frac{q_{ij}+\epsilon}{\widehat p_i+\epsilon},
\qquad q_{ij}=P(R_i=1\mid R_j=1).
\]

This contribution is positive for excitation enrichment and negative for excitation depletion. The default source-support threshold is 20 observations. Overlapping motif membership weights aggregate residue blocks into directional motif scores \(C_{VW}\). A symmetric affinity uses \(|C_{VW}|+|C_{WV}|\), divides by its positive-value 95th percentile, and clips to \([0,1]\) before clustering. Read these motif outputs with their own signed-score definition and units. [Motif implementation in `info_theory.py`](https://github.com/pra-ashok/QProtSim/blob/main/src/qprotsim/info_theory.py)

## Assumptions and interpretation

This is a residue solvation and contact model with a circuit-defined ensemble. The implementation supplies no calibration that identifies the full coupled sample distribution with a canonical equilibrium distribution. MI, DEI, and motif association describe sampled statistical relationships. Interpreting them as biological communication requires experimental or independently modeled evidence.

Record the PDB input, chain, package version, requested and actual backend, shot count, energy table, and coupling schedule with each analysis. Finite-shot variation changes estimated probabilities and information scores. Standard state-vector simulation stores \(2^N\) complex amplitudes, so residue count strongly affects memory requirements. The plotting script in `docs/assets/diagrams/generate.py` regenerates the explanatory figures with NumPy and Matplotlib.
