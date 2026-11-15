# Troubleshooting

Find the symptom below and try the suggested checks. For a first run, use the small CPU example in [First analysis](getting-started.md).

## Installation fails

Check that Python is version 3.11 or later, activate your environment, and update the installer before trying again from the QProtSim folder:

```bash
python --version
python -m pip install --upgrade pip setuptools wheel
python -m pip install ".[cpu]"
```

A “no matching distribution” error can mean that a dependency has no compatible package for your Python version or operating system. Check the package named in the error. CPU and GPU packages can have different platform requirements.

FreeSASA may need a C compiler and Python development headers if a ready-made package is unavailable. Follow the [FreeSASA installation instructions](https://freesasa.github.io/python/intro.html#installation) for your platform.

## The qprotsim command is not found

Activate the same environment you used for installation, then try:

```bash
python -m qprotsim --help
```

If Python also cannot find `qprotsim`, install it in that environment using the [installation guide](installation.md#cpu).

## The PDB file cannot be found or downloaded

For a local file, check that `pdb_path` points to an existing `.pdb` file. Relative paths start from the folder where you run the command. For example:

```yaml
pdb_path: ./input/your_protein.pdb
```

For a download, check the four-character PDB identifier, internet connection, and permission to write to `download_dir`. You can use a previously downloaded file with `pdb_path` instead.

## No residues match the selection

Check the chain identifier in your structure. To include all chains, use:

```yaml
chain: null
```

Keep `protein_only: true` for ordinary protein analysis. Modified or nonstandard residues may lack the reference values needed by the model. Including them with `protein_only: false` does not add those missing values.

## The analysis is slow or runs out of memory

Start with a smaller protein or select a smaller chain. Each selected residue adds a qubit, and the memory used by a dense circuit simulator roughly doubles with each additional qubit. See [memory planning](getting-started.md#go-further) for example sizes.

Reducing `shots` lowers sample storage and some analysis work. It does not reduce the size of the simulated quantum state. Selecting a GPU also requires enough memory on that GPU.

## The requested GPU is unavailable

Install the package matching your GPU and check its hardware and driver requirements in the [installation guide](installation.md#optional-gpu-execution). Run the [backend check](installation.md#verify-the-backend) before a large analysis.

QProtSim can continue on the CPU when the requested backend is unavailable. The device name in the report describes the requested setting, so verify GPU execution separately before comparing performance.

To run explicitly on the CPU, use:

```yaml
pl_device: default.qubit
```

## Some structure pictures are missing

Structure pictures can depend on internet access or the optional PyMOL package. Check whether the PDF and numerical tables were still produced. A missing picture can occur even when the numerical analysis succeeds.

To enable PyMOL where a compatible distribution is available, follow [Optional structure pictures](installation.md#optional-structure-pictures).

## A score is NaN or changes between runs

`NaN` means that a numerical value is undefined. A residue-sensitivity score can be undefined when a residue appears in only one state, or when there is no total divergence to normalise.

Inspect the sampled-state table to see whether both states were observed. More `shots` may help with rare states, but cannot ensure that a state with zero model probability will appear. Do not interpret an undefined score as zero sensitivity.

Small changes between runs are expected because the model is sampled randomly. Compare repeated runs before treating a small difference as meaningful. The [results guide](outputs.md) and [science guide](theory.md) explain the measures and their limitations.

If you use the Python information functions with your own arrays, check the state encoding. For exported Pauli-Z samples, `+1` means ground and `-1` means excited. Convert them explicitly to binary excitation states when the function accepts binary input:

```python
binary_excited = (z_samples == -1).astype(int)
```

This avoids ambiguity when every sampled value is `+1`. See the [Python guide](api.md) for each function's input convention.

## Earlier results have changed or disappeared

QProtSim reuses filenames within an output folder. Save each analysis in a different folder, for example:

```bash
qprotsim analyze -c config.yaml --output-dir ./outputs/repeat
```

Keep the input structure and configuration with the results so that you can repeat the analysis later.

## The basin command rejects a sample file

Basin analysis groups sampled configurations into energy ranges. Its `--z_samples` input must be a CSV with a leading `shot` column and residue columns containing `+1` or `-1`. Use the exported `<molecule>_z_samples.csv` file as the starting point.

If you supply custom reference energies, provide both `--E0` and `--E1` as NumPy `.npy` arrays in the same residue order as the sample columns.
