# First analysis

This guide takes you from a protein structure to a PDF report and data tables. You can follow it without knowing quantum computing or writing Python code.

The example uses **Trp-cage**, a small protein with 20 amino-acid residues. A residue is one amino-acid unit in a protein. QProtSim represents each selected residue with one qubit, a two-state unit in a quantum circuit simulated on your computer.

## 1. Check the installation

Follow the [CPU installation](installation.md#cpu), activate your environment, and run:

```bash
qprotsim --help
```

If the command is unavailable, see [troubleshooting](troubleshooting.md#the-qprotsim-command-is-not-found).

## 2. Save the analysis settings

Create a plain-text file named `config.yaml` in your working folder. YAML is a settings format that uses `name: value` pairs. Copy this example into the file:

```yaml
pdb_code: 1L2Y
molecule_name: Trp Cage
download_dir: ./pdb_cache
chain: A
protein_only: true
shots: 10000
pl_device: default.qubit
bins: 100
output_dir: ./outputs/1l2y
report_file_name: trp_cage_report.pdf
```

The main settings mean:

| Setting | What it does |
| --- | --- |
| `pdb_code: 1L2Y` | Downloads the structure with this identifier from the Protein Data Bank. Internet access is needed for the first download. |
| `chain: A` | Selects protein chain A. A chain is one connected amino-acid sequence in the structure. |
| `shots: 10000` | Takes 10,000 samples from the model. Each sample records a state for every selected residue. |
| `pl_device: default.qubit` | Runs the circuit on your CPU. |
| `bins: 100` | Divides the displayed energy range into 100 histogram intervals. |
| `output_dir` | Chooses the folder where the report and tables will be saved. |

Keep the spelling and indentation shown. Run the next command from the folder containing `config.yaml`; paths beginning with `./` are measured from that folder. The [configuration guide](configuration.md) explains all settings.

## 3. Run the analysis

```bash
qprotsim analyze --config config.yaml
```

QProtSim downloads the structure if needed, samples the circuit, calculates the results, and saves the report. It prints the PDF path when finished.

## 4. Read the results

Open `outputs/1l2y/trp_cage_report.pdf`. The same folder contains CSV tables, which you can open in a spreadsheet, and figures for closer inspection.

Start with three questions:

1. **What energies does the model sample?** Look at the energy distribution.
2. **Which residues stand out?** Compare their sampled ground-state probabilities and energy-sensitivity scores.
3. **Which residues vary together?** Examine the pairwise information plots.

These are results from a simplified residue-state model. They help form hypotheses; biological conclusions need independent evidence. The [results guide](outputs.md) explains how to read each measure, and the [science guide](theory.md) explains the assumptions.

## Use your own structure

To analyse an existing PDB file, replace `pdb_code` and `download_dir` with:

```yaml
pdb_path: ./input/your_protein.pdb
```

The file must exist at that path. Update `molecule_name`, `chain`, and `output_dir` to match your analysis. Begin with a small protein or a small chain: the memory needed to simulate the circuit grows rapidly with the number of residues.

## Repeat an analysis

Choose a new output folder to keep the earlier results:

```bash
qprotsim analyze -c config.yaml --output-dir ./outputs/repeat
```

Files in a reused output folder can be overwritten. Sampling is random, so repeated runs can give slightly different probabilities and scores. Increase `shots` when you need more precise estimates, and check whether the patterns persist across runs. More samples cannot correct the model's scientific assumptions.

## Go further

Use the [configuration guide](configuration.md) to select other inputs or analyse several proteins. The [Python guide](api.md) shows how to work with the results in scripts and notebooks.

For planning larger analyses, a dense simulator needs approximately $16\times2^N$ bytes just to store the state for $N$ residues: 16 MiB at 20 residues, 512 MiB at 25, and 16 GiB at 30. Sampling, plots, and other calculations need additional memory. GPU execution also needs enough memory on the GPU.
