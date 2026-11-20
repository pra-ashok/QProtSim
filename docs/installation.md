# Installation

Start with the CPU installation below. It runs on your computer's processor and does not require a quantum computer or a graphics card. GPU options are available for users with compatible hardware.

## Get QProtSim

You need **Python 3.11 or later** and a copy of the [QProtSim repository](https://github.com/pra-ashok/QProtSim). Download and extract the repository ZIP, then open a terminal in the extracted folder. If you use Git, you can get the same files with:

```bash
git clone https://github.com/pra-ashok/QProtSim.git
cd QProtSim
```

The commands on this page run from that folder. Check your Python version:

```bash
python --version
```

On systems where the command is `python3`, use `python3` to create the environment in the next step.

## Create an environment

A virtual environment keeps QProtSim's packages separate from your other Python tools.

```bash
python -m venv .venv
```

Activate it using the command for your terminal:

| Terminal | Activation command |
| --- | --- |
| Linux or macOS | `source .venv/bin/activate` |
| Windows PowerShell | `.venv\Scripts\Activate.ps1` |
| Windows Command Prompt | `.venv\Scripts\activate.bat` |

Then update the package installer:

```bash
python -m pip install --upgrade pip
```

Activate this environment again whenever you open a new terminal to use QProtSim.

## CPU

Install QProtSim and its CPU simulator:

```bash
python -m pip install ".[cpu]"
qprotsim --help
```

The dot in `.[cpu]` means “install from this folder.” If the help text appears, continue to [your first analysis](getting-started.md).

Use this setting in your analysis configuration:

```yaml
pl_device: default.qubit
```

`pl_device` names the simulator that runs the quantum circuit. The CPU installation also includes `lightning.qubit`, an alternative CPU simulator. You can select it by changing that value.

If installation reports a missing compiler or a failed FreeSASA build, see [installation troubleshooting](troubleshooting.md#installation-fails). FreeSASA is the tool QProtSim uses to estimate how much of each residue is exposed to solvent.

## Optional Qiskit execution

Install the Qiskit SDK and Aer simulator through the extra declared in `pyproject.toml`:

```bash
python -m pip install ".[qiskit]"
```

Select Qiskit in your YAML configuration:

```yaml
backend: qiskit
qiskit_backend: aer_simulator
backend_options:
  method: statevector
  device: CPU
seed: 42
```

QProtSim builds a native Qiskit circuit with the same initial rotations and coupling gates as its PennyLane circuit. `qiskit_backend` selects a local Aer backend; `backend_options` supplies its simulator options. See the [Aer simulator documentation](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html) for available methods and devices. This installation uses CPU simulation; Aer GPU execution requires a compatible Aer GPU installation and hardware.

Run the supplied example with `qprotsim analyze -c input/1l2y_qiskit_config.yaml`. It uses the local `input/1L2Y.pdb` structure. Missing Qiskit dependencies or an invalid Qiskit backend produce an error. Existing configurations without `backend` continue to use PennyLane.

## Optional GPU execution

A GPU is a graphics processor that can also run numerical calculations. Choose the option matching your hardware and system software:

| Hardware | Install from the QProtSim folder | Configuration setting |
| --- | --- | --- |
| Supported NVIDIA GPU with CUDA | `python -m pip install ".[cuda]"` | `pl_device: lightning.gpu` |
| Supported AMD GPU with ROCm | `python -m pip install ".[rocm]"` | `pl_device: lightning.amdgpu` |

These packages use PennyLane Lightning 0.44.0. Installing a GPU package does not install your graphics driver or make every GPU compatible. Check the requirements for your chosen simulator before installing: [NVIDIA requirements](https://github.com/PennyLaneAI/pennylane-lightning/blob/v0.44.0/doc/lightning_gpu/device.rst) · [AMD requirements](https://github.com/PennyLaneAI/pennylane-lightning/blob/v0.44.0/doc/lightning_amdgpu/device.rst).

### Verify the backend

A backend is the simulator used to execute the circuit. The PennyLane execution path can continue on the CPU when a requested device is unavailable. The report lists the requested device, so that label alone does not confirm GPU execution. The Qiskit execution path reports backend failures directly.

Before measuring GPU performance, run this small check in Python. Replace `lightning.qubit` with your chosen device name. A successful run prints the device name and a sample shape of `(100, 2)`.

```python
import pennylane as qml

device = qml.device("lightning.qubit", wires=2)

@qml.set_shots(100)
@qml.qnode(device)
def check_backend():
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    return qml.sample(wires=[0, 1])

samples = check_backend()
print(device.name, samples.shape)
```

If this check fails, resolve the backend error or use `default.qubit` for CPU execution.

## Optional structure pictures

QProtSim can include pictures of the protein in its report. To enable PyMOL rendering where a compatible package is available, install:

```bash
python -m pip install ".[visualization]"
```

PyMOL is optional. Other image sources and simpler structure drawings are available when it is absent.
