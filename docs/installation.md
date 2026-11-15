# Installation

## Choose an environment

QProtSim requires Python 3.11 or later. The examples below install the checked-out source tree from the repository root. Use a virtual environment and a recent pip.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

On Windows, activate the environment with `.venv\Scripts\activate`. Accelerator wheels have additional platform requirements described below.

FreeSASA and other scientific dependencies contain compiled components. If a wheel is unavailable, install your platform's compiler and Python development headers before building from source. [FreeSASA installation](https://freesasa.github.io/python/intro.html#introduction)

| Extra | PennyLane device | Intended execution |
| --- | --- | --- |
| `cpu` | `default.qubit` or `lightning.qubit` | CPU |
| `cuda` | `lightning.gpu` | Supported NVIDIA GPU on Linux |
| `rocm` | `lightning.amdgpu` | Supported AMD GPU with ROCm |

Installing an extra makes its backend available. Select the desired device through YAML `pl_device` or `Config(pl_device=...)`. The default is `default.qubit`. PennyLane and the Lightning plugins are kept at version 0.44.0 to align their interfaces.

## CPU

```bash
python -m pip install ".[cpu]"
```

The base installation, `python -m pip install .`, also supports CPU execution. The explicit `cpu` extra provides a consistent installation command for backend selection.

```yaml
pl_device: default.qubit
```

For the compiled Lightning CPU simulator, select `lightning.qubit`.

```yaml
pl_device: lightning.qubit
```

## NVIDIA CUDA

```bash
python -m pip install ".[cuda]"
```

```yaml
pl_device: lightning.gpu
```

This installs `pennylane-lightning-gpu==0.44.0`. Its wheel dependencies supply cuStateVec and CUDA 12 user-space libraries. A compatible NVIDIA driver and supported GPU are required separately. The pinned release requires compute capability SM 7.0 or later. Follow the upstream release's platform and driver requirements. [Lightning-GPU 0.44 device requirements](https://github.com/PennyLaneAI/pennylane-lightning/blob/v0.44.0/doc/lightning_gpu/device.rst)

CuPy is optional in QProtSim. The quantum simulator uses the Lightning plugin, and the helper module supports NumPy when CuPy is absent. Installing CuPy separately requires a build matching your CUDA environment.

## AMD ROCm

```bash
python -m pip install ".[rocm]"
```

```yaml
pl_device: lightning.amdgpu
```

This installs `pennylane-lightning-amdgpu==0.44.0`. Upstream prebuilt wheels target Linux x86_64, AMD MI300-series GPUs, and ROCm 7.0 or later. Other AMD architectures or ROCm versions require a matching Lightning source build with HIP and the appropriate Kokkos architecture flags. The extra installs the Python backend, while the system ROCm stack must be available. [Lightning-AMDGPU 0.44 device requirements](https://github.com/PennyLaneAI/pennylane-lightning/blob/v0.44.0/doc/lightning_amdgpu/device.rst) · [Source-build instructions](https://github.com/PennyLaneAI/pennylane-lightning/blob/v0.44.0/doc/lightning_amdgpu/installation.rst)

## Verify the backend

Execute a small circuit directly before running a large analysis. Set `backend` to the device you intend to use. This call reaches PennyLane directly and reports backend errors.

```python
import pennylane as qml

backend = "lightning.qubit"
device = qml.device(backend, wires=2)

@qml.set_shots(100)
@qml.qnode(device)
def check_backend():
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    return qml.sample(wires=[0, 1])

samples = check_backend()
assert samples.shape == (100, 2)
print(device.name)
```

QProtSim's device factory falls back to `default.qubit` when a requested device fails to initialize. Its CUDA execution path can also retry on CPU after GPU-related errors. The report records the requested device, so the direct check above helps verify accelerator availability.

## Optional tools

```bash
python -m pip install -e ".[cpu,docs]"
python -m pip install -e ".[cpu,dev]"
python -m pip install ".[visualization]"
```

The `docs` extra installs MkDocs, the Shadcn theme, PyMdown Extensions, and Pygments. The `dev` extra includes the documentation tools and development dependencies. The `visualization` extra requests PyMOL for local structure rendering when a compatible distribution is available. The report workflow can also use RCSB images and Matplotlib rendering.

Continue with a [first analysis](getting-started.md) or [build the documentation](development.md).
