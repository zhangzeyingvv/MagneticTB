# MagneticTB

[English](README.md) | [简体中文](README.zh-CN.md)

MagneticTB is a Wolfram Language package for constructing symmetry-constrained
tight-binding Hamiltonians. It supports magnetic and nonmagnetic crystal
symmetries, local orbital bases, symbolic hopping Hamiltonians, band structures,
and related analysis and visualization tools.

The current paclet requires Wolfram Language 12.1 or later. English and
Simplified Chinese documentation are included.

QQ discussion group: **625192239**.

## Install, update, or uninstall

MagneticTB is installed from a local `.paclet` archive. This repository does
not provide an automatic online update channel.

In Mathematica, set the absolute path of your local `.paclet` archive and
install it:

```wl
pacletArchive =
  "/absolute/path/to/MagneticTB/build/MagneticTB-2.0.10.paclet";

PacletInstall[pacletArchive]
```

After installation, save your work, fully quit Mathematica, and reopen it
(not just its kernel). Then verify and load the package:

```wl
PacletFind["MagneticTB"]
Needs["MagneticTB`"]
```

To update, install a newer local archive with the same paclet name, then fully
quit and reopen Mathematica:

```wl
PacletInstall["/absolute/path/to/MagneticTB-<new-version>.paclet"]
```

`PacletInstall::samevers` means that the same version is already installed. To
reinstall that exact version, uninstall it first:

```wl
PacletUninstall["MagneticTB"]
PacletInstall["/absolute/path/to/MagneticTB-<version>.paclet"]
```

To remove MagneticTB without reinstalling it:

```wl
PacletUninstall["MagneticTB"]
PacletFind["MagneticTB"]
```

After an update, reinstall, or uninstall, save your work and fully quit and
reopen Mathematica. This avoids continuing with previously loaded definitions
or an old help page left open in the front end.

## Open the Mathematica documentation

After installing the paclet and restarting Mathematica as described above,
search for `MagneticTB` in the Documentation Center or open the help homepage
directly:

```wl
SystemOpen["paclet:MagneticTB/guide/MagneticTB"]
```

A particular function page can be opened in the same way:

```wl
SystemOpen["paclet:MagneticTB/ref/init"]
```

In a notebook, placing the cursor on a MagneticTB function and pressing the
standard function-help key also opens its reference page. Mathematica selects
the installed English or Simplified Chinese documentation according to the
front end language.

## Basic example: graphene Hamiltonian and bands

The following example constructs a spinless graphene model from the complete
ordered symmetry-operation list returned by `msgop`, obtains the onsite through
next-nearest-neighbour Hamiltonian, and plots its bands.

```wl
Needs["MagneticTB`"];

grapheneOperations = msgop[gray[191]];

init[
  lattice -> {
    {a, 0, 0},
    {-a/2, Sqrt[3] a/2, 0},
    {0, 0, c}
  },
  lattpar -> {a -> 1, c -> 3},
  wyckoffposition -> {
    {{1/3, 2/3, 0}, {0, 0, 0}}
  },
  symminformation -> grapheneOperations,
  basisFunctions -> {{"pz"}}
];

grapheneHamiltonian = Total[symham /@ Range[3]];
MatrixForm[grapheneHamiltonian]
```

The path below is given in reciprocal fractional coordinates. Fixing the model
parameters produces a reproducible band plot:

```wl
graphenePath = {
  {{{0, 0, 0}, {0, 1/2, 0}}, {"G", "M"}},
  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
  {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}}
};

bandplot[
  graphenePath,
  60,
  grapheneHamiltonian,
  {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}
]
```

For longer workflows, open the `MagneticTB` help homepage and the
`Getting Started with MagneticTB` tutorial included with the paclet.

## Build the paclet from source

Run the release builder from the repository root:

```bash
sh Developer/Paclet/BuildRelease.sh
```

On macOS the script uses
`/Applications/Mathematica.app/Contents/MacOS/WolframKernel` by default. If the
kernel is installed elsewhere, set its absolute path explicitly:

```bash
WOLFRAM_KERNEL="/absolute/path/to/WolframKernel" \
  sh Developer/Paclet/BuildRelease.sh
```

The builder stages the runtime files and bilingual documentation, then writes
the installable archive to:

```text
build/MagneticTB-<version>.paclet
```

Generated `build/` and `tmp/` directories are local build output and are not
part of the maintained source tree.

## Citation

If MagneticTB contributes to published work, cite:

Zeying Zhang, Zhi-Ming Yu, Gui-Bin Liu, and Yugui Yao,
“MagneticTB: A package for tight-binding model of magnetic and non-magnetic
materials,” *Computer Physics Communications* **270**, 108153 (2022),
[doi:10.1016/j.cpc.2021.108153](https://doi.org/10.1016/j.cpc.2021.108153).

```bibtex
@article{Zhang2022MagneticTB,
  author  = {Zhang, Zeying and Yu, Zhi-Ming and Liu, Gui-Bin and Yao, Yugui},
  title   = {MagneticTB: A package for tight-binding model of magnetic and non-magnetic materials},
  journal = {Computer Physics Communications},
  volume  = {270},
  pages   = {108153},
  year    = {2022},
  doi     = {10.1016/j.cpc.2021.108153}
}
```
