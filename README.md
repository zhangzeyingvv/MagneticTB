# MagneticTB

English | [简体中文](README.zh-CN.md)

MagneticTB constructs symmetry-constrained tight-binding Hamiltonians for
magnetic space groups, nonmagnetic space groups, and spin-space groups. Given
the symmetry and orbital information at the selected Wyckoff positions, it
generates the symmetry-allowed Hamiltonian and provides tools for symmetry
operations, band structures, and related tight-binding calculations.

This repository contains two independently usable implementations:

- [`Mathematica/`](Mathematica/): the Wolfram Language package, magnetic
  symmetry data, English and Simplified Chinese documentation, and example
  notebooks.
- [`Python-Rust-Web/`](Python-Rust-Web/): the Rust computational core, Python
  API, Web interface, runtime data, and user documentation.

The two implementations do not require one another at runtime.

## Installation

### Mathematica

The current Mathematica paclet version is `2.0.10` and requires Wolfram
Language 12.1 or later.

Download `MagneticTB-2.0.10.paclet` from the release assets, then install it in
a Mathematica kernel using its absolute path:

```wl
PacletInstall["/absolute/path/to/MagneticTB-2.0.10.paclet"]
```

After installation, quit and restart Mathematica. Then verify and load the
package in a new kernel:

```wl
PacletFind["MagneticTB"]
Needs["MagneticTB`"]
```

To open the installed documentation, choose **Help > Wolfram Documentation**
in Mathematica and search for `MagneticTB`. Open the **MagneticTB** guide from
the search results to browse the function pages, tutorials, and examples. You
can also search for a function name such as `init`, `initfromrep`, or `symham`
to open its reference page directly.

See the [Mathematica README](Mathematica/README.md) for updating, uninstalling,
documentation, and examples.

### Python + Rust

The Python interface requires CPython 3.9 or later. Download the `.whl` file
matching your Python version, operating system, and CPU architecture from the
release assets. The wheel includes the compiled Rust core; users do not need
to install Rust or Cargo separately.

Using `venv` on macOS, Linux, or WSL:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install /absolute/path/to/downloaded-wheel.whl
.venv/bin/python -c "import magnetictb; print(magnetictb.__version__)"
.venv/bin/magnetictb-web
```

Using `venv` in Windows PowerShell:

```powershell
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install C:\absolute\path\to\downloaded-wheel.whl
.venv\Scripts\python.exe -c "import magnetictb; print(magnetictb.__version__)"
.venv\Scripts\magnetictb-web.exe
```

Using Conda is recommended. Open a Conda-enabled terminal (Anaconda Prompt or
Miniforge Prompt on Windows), create an environment, and install the wheel
with `pip` inside that environment:

```sh
conda create -n magnetictb --override-channels -c conda-forge python=3.12 pip
conda activate magnetictb
python -m pip install /absolute/path/to/downloaded-wheel.whl
python -c "import magnetictb; print(magnetictb.__version__)"
magnetictb-web
```

The current package version prints `0.1.0`.

After starting `magnetictb-web`, open <http://127.0.0.1:8000/> in a browser.
Interactive API documentation is available at <http://127.0.0.1:8000/api/docs>.
The service listens only on the local machine by default; press `Ctrl+C` to
stop it.

See the [Python/Rust/Web README](Python-Rust-Web/README.md) and the
[Web help center](Python-Rust-Web/docs/user-guide/web/index.html) for the
modeling workflow, API, exact inputs, and examples.

## Capabilities

- Construct symmetry-constrained tight-binding Hamiltonians.
- Work with magnetic, nonmagnetic, and spin-space-group symmetry data.
- Obtain matrix representations of symmetry operations.
- Generate real-space and momentum-space Hamiltonians by bond shell.
- Manipulate and analyze band structures and related model properties.
- Use either the Mathematica interface or the Python API backed by the Rust
  computational core.

## Examples and documentation

- Mathematica examples: [`Mathematica/Examples/`](Mathematica/Examples/)
- Mathematica bilingual help: [`Mathematica/Documentation/`](Mathematica/Documentation/)
- Python/Rust Web help: [`Python-Rust-Web/docs/user-guide/web/`](Python-Rust-Web/docs/user-guide/web/index.html)

## Release Notes

Releases are listed from newest to oldest.

### Python/Rust 0.1.0 (2026-08-18)

- Added the Rust computational core, Python API, and local Web interface.

### Mathematica 2.0.10 (2026-08-17)

- Expanded the English and Simplified Chinese documentation, tutorials, and
  quick-start material.
- Added real-space Hamiltonian, slab, surface Green-function, Berry geometry,
  Wilson-loop, crystal, Brillouin-zone, and k-path tools, Improved band visualization.

### Mathematica 2.0.0  (2026-03-15)

- Rewrote the core algorithms using linear algebra and group
  representation theory.
- Added the induced-representation mode, enabling construction of minimal
  tight-binding models.
- Added `initfromrep`, enabling construction of a tight-binding model from
  site symmetry group representation input without basis-functions.
- Added full support for spin-space groups (SSGs), including collinear,
  coplanar, and non-coplanar cases.
- Added the cyclotomic exact null-space kernel as a `KernelMethod` available to
  `symham`.

### Mathematica 1.06 (2025-12-11)

- Added beta support for spin-space groups.
- Added magnetic-space-group Wyckoff-position display and symmetry-operation
  data for magnetic layer and rod groups.

### Mathematica 1.05 (2024-12-04)

- Fixed a bug in `hop`.
- Added `readHR` for importing `wannier90_hr.dat` files.

### Mathematica 1.04 (2024-06-21)

- Added the `CartesianCoordinates` option to `symham`.
- Added `banddata` for generating `band.dat` files.

### Mathematica 1.03 (2023-02-17)

- Added a greedy algorithm for automatically finding space-group generators,
  significantly improving computational efficiency.

### Mathematica 1.02b (2023-02-14)

- Added an example showing how to obtain tight-binding parameters by hand.
- Added the English manual.

### Mathematica 1.02 (2022-12-01)

- Added tight-binding band co-representation calculations using the optional
  `SpaceGroupIrep` and `MSGCorep` packages.
- Added `getMSGElemFromMSGCorep` and `getTBBandCorep`.

### Mathematica 1.01 (2022-07-22)

- Fixed a rare basis-function ordering change caused by automatic
  unitarization.

### Mathematica 1.00c

- Added the Chinese manual.
- Fixed monoclinic lattice-vector display.

### Mathematica 1.00b

- Added the charge-4 Weyl-point example for double magnetic space groups.

## Citation

If MagneticTB is useful in your research, please cite:

Z. Zhang, Z.-M. Yu, G.-B. Liu, and Y. Yao, “MagneticTB: A package for
tight-binding model of magnetic and nonmagnetic materials,” *Computer Physics
Communications* **270**, 108153 (2022).

- [Journal article](https://www.sciencedirect.com/science/article/abs/pii/S0010465521002654)
- [arXiv:2105.09504](https://arxiv.org/abs/2105.09504)

## License

MagneticTB is licensed under the GNU General Public License version 3 only
(`GPL-3.0-only`). See [LICENSE](LICENSE).

Copyright (C) 2021-2026 Zhang Zeying.
