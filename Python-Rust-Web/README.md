# MagneticTB Python / Rust / Web

English | [简体中文](README.zh-CN.md)

MagneticTB constructs symmetry-constrained tight-binding Hamiltonians from a lattice, Wyckoff positions, magnetic space-group data, and local orbitals.

- **Rust** is the sole mathematical core. It implements groups, representations, bond constraints, exact null spaces, Hamiltonians, and physical-property calculations.
- **Python** provides the user-facing API without reimplementing the Rust mathematics.
- **Web** is installed with the Python package and provides a local browser interface.

The current package version is `0.1.0` and requires CPython 3.9 or newer. See the [Python user guide](docs/user-guide/web/index.html) for the web documentation.

## Install a wheel (recommended)

Most users only need Python and a matching wheel. Rust, Cargo, and compiler tools are not required.

`PLATFORM` below is a placeholder. Replace it with the platform tag in the actual wheel filename.

Installing the wheel installs both the Python API and the Web runtime dependencies. No separate installation option is needed.

macOS / Linux / WSL:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install ./magnetictb-0.1.0-cp39-abi3-PLATFORM.whl
```

Windows PowerShell:

```powershell
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install .\magnetictb-0.1.0-cp39-abi3-PLATFORM.whl
```

Verify the installation on macOS / Linux / WSL:

```bash
.venv/bin/python -c "import magnetictb; print(magnetictb.__version__)"
```

On Windows PowerShell:

```powershell
.venv\Scripts\python.exe -c "import magnetictb; print(magnetictb.__version__)"
```

### Install the wheel with Conda

Use an existing Conda installation, or install [Miniforge](https://github.com/conda-forge/miniforge#install) for your operating system and CPU. On Windows, open Miniforge Prompt or Anaconda Prompt; on macOS / Linux, use a Conda-enabled terminal.

The following commands replace the `.venv` steps above. Run the installation command in the wheel directory and replace `PLATFORM` with the tag in your actual wheel filename:

```bash
# Skip creation if this environment already exists; activate it instead.
conda create -n magnetictb --override-channels -c conda-forge python=3.12 pip
conda activate magnetictb

# Install the Python API and Web dependencies together; no Rust is needed.
python -m pip install "./magnetictb-0.1.0-cp39-abi3-PLATFORM.whl"
python -c "import magnetictb; print(magnetictb.__version__)"

# Start the local Web interface.
magnetictb-web
```

Open <http://127.0.0.1:8000/>. In a new terminal, activate `magnetictb` before running `python` or `magnetictb-web`. Do not create another `.venv` inside this environment. Conda manages Python; pip installs MagneticTB from the local wheel.

## Start the Web interface

The wheel installs the Web runtime dependencies. Start the service in the same environment, from the directory containing `.venv`:

Conda users run `conda activate magnetictb` followed by `magnetictb-web`. The paths below are for `.venv` users.

macOS / Linux / WSL:

```bash
.venv/bin/magnetictb-web
```

Windows PowerShell:

```powershell
.venv\Scripts\magnetictb-web.exe
```

Open the following local URLs in a browser:

- Web interface: <http://127.0.0.1:8000/>
- API documentation: <http://127.0.0.1:8000/api/docs>
- Health check: <http://127.0.0.1:8000/api/health>

The service listens only on `127.0.0.1:8000` by default. Press `Ctrl+C` to stop it.

## Basic examples

### Query magnetic space-group operations

```python
from magnetictb import gray, msgop

operations = msgop(gray[191])
print(len(operations))
print(operations[0])
```

`msgop` returns a complete, ordered ordinary Python list. Each operation has the form:

```text
[label, rotation, translation, "F" or "T"]
```

### Construct a Hamiltonian

The following is a minimal three-dimensional model with one site and one `s` orbital:

```python
from magnetictb import SymmetryOperation, init, symham

identity = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)

init(
    lattice=identity,
    wyckoffposition=(((0, 0, 0), (0, 0, 0)),),
    symminformation=(SymmetryOperation("E", identity),),
    basis_functions=("s",),
    initial_bond_shells=2,
)

onsite = symham(1)
nearest = symham(2)

print(onsite)
print(nearest)
```

On success, `init(...)` returns `None` and installs the current model. `symham(shell)` returns the complete two-dimensional symbolic matrix for one physical bond shell. Shell numbering starts at 1.

See the [Python user guide](docs/user-guide/web/index.html) for Graphene, the three-band MoS₂ model, exact inputs, multi-shell combinations, and explicit representations. The Web home page also includes runnable Graphene, MoS₂, and additional crystal-system examples.

## Build from source (optional)

Most users can skip this section: installing a wheel needs no Rust toolchain.

For detailed platform instructions, see [Build from source](docs/user-guide/web/guide/BuildFromSource.html).

### 1. Prerequisites

Building from source requires:

- CPython 3.9 or newer;
- Rust 1.85 or newer;
- Cargo, which is installed together with Rust by `rustup`;
- `maturin==1.14.1`;
- MSVC C++ Build Tools on Windows, Xcode Command Line Tools on macOS, or GCC/Clang and a linker on Linux.

Check Python, Rust, and Cargo first:

```text
python3 --version
rustc --version
cargo --version
```

On Windows, use `py -3 --version` instead of `python3 --version`. Platform-specific toolchain instructions are available in the [source-build guide](docs/user-guide/web/guide/BuildFromSource.html).

### 2. macOS / Linux / WSL

Enter the repository's `python/` directory:

```bash
cd python
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install "maturin==1.14.1"
.venv/bin/maturin develop --release --locked
```

Verify the installation:

```bash
.venv/bin/python -c "import magnetictb; print(magnetictb.__version__)"
```

The expected output is:

```text
0.1.0
```

### 3. Windows PowerShell

Enter the repository's `python` directory:

```powershell
cd python
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install "maturin==1.14.1"
.venv\Scripts\maturin.exe develop --release --locked
```

Verify the installation:

```powershell
.venv\Scripts\python.exe -c "import magnetictb; print(magnetictb.__version__)"
```

`maturin develop` compiles the Rust core and installs the Python package into the current `.venv`. It does not create a remote release or upload any files.

### Build in Conda

For a source installation instead, install the platform build tools described above, activate this Conda environment, and run from the repository root:

```bash
# Build into the active Conda environment.
cd python
python -m pip install "maturin==1.14.1"
python -m maturin develop --release --locked
```

### Build a wheel

From the `python/` directory, using the virtual environment created above:

#### macOS / Linux / WSL

```bash
.venv/bin/maturin build --release --locked --out dist
```

#### Windows PowerShell

```powershell
.venv\Scripts\maturin.exe build --release --locked --out dist
```

The wheel is written to `python/dist/`. Its filename contains operating-system, CPU-architecture, and ABI tags, for example:

```text
magnetictb-0.1.0-cp39-abi3-PLATFORM.whl
```

The wheel contains the compiled Rust crate code and the Web static assets. A wheel user does not need Rust, Cargo, or maturin, but the wheel must match the user's operating system and CPU architecture.

## Citation

If you use MagneticTB in research, please cite:

> Zeying Zhang, Zhi-Ming Yu, Gui-Bin Liu, Yugui Yao, “MagneticTB: A package for tight-binding model of magnetic and non-magnetic materials,” *Computer Physics Communications* **270**, 108153 (2022). <https://doi.org/10.1016/j.cpc.2021.108153>

```bibtex
@article{ZHANG2022108153,
  title   = {MagneticTB: A package for tight-binding model of
             magnetic and non-magnetic materials},
  journal = {Computer Physics Communications},
  volume  = {270},
  pages   = {108153},
  year    = {2022},
  doi     = {10.1016/j.cpc.2021.108153},
  author  = {Zeying Zhang and Zhi-Ming Yu and
             Gui-Bin Liu and Yugui Yao}
}
```
