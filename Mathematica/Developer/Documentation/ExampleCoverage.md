# MagneticTB 2.0 example coverage

This file is the release gate required by `AGENTS.md`.  A row is complete only
when the help page contains the full input, an output produced by that input,
and the stated physical validation has been rerun in a fresh kernel.  English
and Simplified Chinese pages use identical Wolfram Language input and output.

## `Examples/GeneralExamples.nb`

The source notebook currently contains 52 `Input` cells and 35 `Output` cells.
The workflows are reorganized into readable sections in the official tutorial;
they are not replaced by smaller toy models.

| Source workflow | Help URI | Executed user-visible output | Physical or numerical validation |
|---|---|---|---|
| Graphene `pz` model | `MagneticTB/tutorial/GeneralExamples` and `MagneticTB/tutorial/GettingStarted` | orbital table, bond table, shell Hamiltonians, total Hamiltonian, interactive and static bands | exact Hermiticity; unitary and antiunitary covariance; non-Gamma checks are covered by the Wyckoff regression |
| Spinful graphene and SOC workflow | `MagneticTB/tutorial/GeneralExamples` | four-band symbolic Hamiltonian, `bandManipulate`, fixed-parameter bands, and exported real-space hopping text | the ordered spinor basis and antiunitary representation are evaluated before export |
| MoS2 three-band model | `MagneticTB/tutorial/GeneralExamples` | symbolic Hamiltonian, Cartesian-coordinate Hamiltonian, band plot, and the contents of exported `band.dat` | the two momentum-coordinate conventions give the same sampled bands |
| Spinful MoS2 in two local basis orders | `MagneticTB/tutorial/GeneralExamples` | both orbital tables and the corresponding onsite Hamiltonians | the row/column permutation is read from `orbitalTable` and retained in each rebuilt Hamiltonian |
| MoS2 lowest-band constant-energy contours | `MagneticTB/tutorial/GeneralExamples` | contour graphics over several Brillouin zones and the reciprocal primitive-cell polygon | the plotted cell is computed by `init`, not entered as a decorative polygon |
| Explicit broken-symmetry/subgroup reinitialization | `MagneticTB/tutorial/GeneralExamples` | complete rules from `brokenSymmetryInitRules`, a new initialized session, and the resulting Hamiltonian | the removed `symmetryset` workflow is replaced by an explicit initialization boundary |
| Magnetic C3 Weyl point | `MagneticTB/tutorial/GeneralExamples` | Weyl Hamiltonian, band plot, Berry-curvature/Chern workflow result | Chern-number phase flow is evaluated from the displayed model |
| Magnetic cubic nodal line | `MagneticTB/tutorial/GeneralExamples` | symbolic Hamiltonian and band plot along the chosen path | nodal degeneracy is visible on the symmetry line |
| C4 topological-insulator model | `MagneticTB/tutorial/GeneralExamples` | symbolic Hamiltonian and topological-model band plot | C4-constrained spectrum is evaluated from the documented Hamiltonian |
| MSG 198.11 spinless model | `MagneticTB/tutorial/GeneralExamples` | symbolic Hamiltonian and bands | generated representation matrices and Hamiltonian covariance are evaluated |
| MSG 198.11 spinful C4 Weyl model | `MagneticTB/tutorial/GeneralExamples` | spinful Hamiltonian and bands | double-valued representation and symmetry-enforced degeneracy are evaluated |
| Magnetic Wyckoff-position display | `MagneticTB/tutorial/GeneralExamples` and `MagneticTB/ref/showMSGWyckoff` | the actual Wyckoff grid | result is read from the packaged database, not handwritten |
| Magnetic layer-group operations | `MagneticTB/tutorial/GeneralExamples` and `MagneticTB/ref/mlgop` | the actual operation list | operation records are read from the packaged database |
| Magnetic rod-group operations | `MagneticTB/tutorial/GeneralExamples` and `MagneticTB/ref/mrgop` | the actual operation list | operation records are read from the packaged database |

## Models used by the 2.0 regression suite

| Model | Help URI | Executed user-visible output | Physical or numerical validation |
|---|---|---|---|
| `SimpleCubicS` | `MagneticTB/tutorial/ValidatedModels` | orbital table, bond table, shell Hamiltonian, bands | exact Hermiticity and shell-boundary behavior |
| `GraphenePz` | `MagneticTB/tutorial/GettingStarted` | two-site orbital table, shell Hamiltonians, symmetry matrices, bands | strict corepresentation multiplication and full k-space covariance |
| `CsClRectangular1x3` | `MagneticTB/tutorial/GettingStarted` and `MagneticTB/tutorial/ValidatedModels` | explicit `1 x 3` cross block, full Hamiltonian, hopping-source table, bands | rectangular-basis rank and exact Hermiticity |
| `SpinfulSquarePxPy` | `MagneticTB/tutorial/ValidatedModels` | onsite block, spinful representation matrices, bands | double-valued multiplication and Hamiltonian covariance |
| `P4Direct` / `P4Induced` | `MagneticTB/tutorial/ValidatedModels` | both Hamiltonians, both sets of full representation matrices, both band plots | equal real parameter subspaces in both directions and aligned bands |
| `FeSeTwoWyckoffRectangular3x2` | `MagneticTB/tutorial/ValidatedModels` | orbital table, `3 x 2` block, full Hamiltonian, bands | multiple Wyckoff orbits, rectangular hopping, exact Hermiticity |
| `CollinearDiscreteSSG` | `MagneticTB/tutorial/ValidatedModels` | Hamiltonian, spin-space representation matrices, bands | strict multiplication and covariance for all 8 operations |
| `NoncollinearOctahedralSSG` | `MagneticTB/tutorial/ValidatedModels` | Hamiltonian, representative spin matrices, bands | strict multiplication and covariance for all 24 operations |
| MoS2 11-band and FeSe 10-band applications | `MagneticTB/tutorial/GeneralExamples` and `MagneticTB/tutorial/ValidatedModels` | physically selected Hamiltonian blocks and band plots | multiple orbitals and multiple Wyckoff-orbit assembly |
| Pure-internal C4 and half-translation/time-reversal SSG models | `MagneticTB/tutorial/ValidatedModels` | full matrices, Hamiltonians, bands | unitary and semilinear antiunitary covariance |
| Two-site `E + PT` induced model | `MagneticTB/tutorial/ValidatedModels` | induced swap matrix, Hamiltonian, bands | antiunitary coset representative, strict corepresentation multiplication, semilinear covariance |
| `C3 x Z2^T / C3` complex local character | `MagneticTB/tutorial/ValidatedModels` | all six full matrices, induced Hamiltonian, bands | the target antiunitary representative uses the conjugated local block; strict multiplication |
| One-generator `C_infinity` model | `MagneticTB/tutorial/ContinuousSymmetry` | `D(theta)`, orbital table, constrained Hamiltonian | exact differentiation to the Hermitian generator followed by finite constraints |
| Spin-layer group `14.1.2.3.L.1`, `000` point | `MagneticTB/tutorial/ContinuousSymmetry` | bundled exact finite matrices, four-band DirectProduct model, two-band Induced model, bands | immutable matrices generated from SpinLayerCorepresentations, continuous constraint, and finite strict multiplication |

## Seven crystal systems from `wyckoffMSG.mx`

Every row is evaluated twice, with `RepresentationMode -> "DirectProduct"`
and `RepresentationMode -> "Induced"`.  The page displays the database seed,
both concrete Hamiltonians, representative matrices, interactive band controls,
and fixed-parameter bands.

| Crystal system | BNS / Wyckoff | Help URI | Executed output and validation |
|---|---|---|---|
| Triclinic | 1.3 `a` | `MagneticTB/tutorial/WyckoffCrystalSystems` | database seed, DirectProduct/Induced Hamiltonians and bands; full comparison below |
| Triclinic | 2.6 `i` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Monoclinic | 3.3 `e` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Monoclinic | 3.4 `a` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Orthorhombic | 16.3 `q` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Orthorhombic | 16.3 `r` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Tetragonal | 75.3 `c` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Trigonal | 143.3 `a` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Hexagonal | 168.111 `b` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |
| Cubic | 195.3 `a` | `MagneticTB/tutorial/WyckoffCrystalSystems` | same |

For all ten cases the release regression checks, independently of the visible
page, both directions of the DirectProduct/Induced real parameter-subspace
comparison, aligned bands at five momenta, strict representation or
corepresentation multiplication, six generic-k covariance reports, and all
non-Gamma little groups on the `{0, Pi/2, Pi}^3` grid with the reciprocal
sewing matrix.  A Gamma-only result is never used as the space-group proof.

## Archived backend

| Source workflow | Help URI | Executed user-visible output | Validation |
|---|---|---|---|
| Complete legacy Graphene model | `MagneticTB/tutorial/LegacyBackend`, `MagneticTB/ref/initold`, `MagneticTB/ref/symhamold`, and `MagneticTB/ref/bandManipulateold` | legacy orbital/bond preparation, shell Hamiltonians, total Hamiltonian, interactive and static bands | evaluated in a dedicated old-only child kernel; no `MagneticTB`` symbols are loaded |

The new and old backends are never loaded in the same kernel.  Cross-version
comparison passes only neutral Hamiltonian data between separate child kernels.
