# MagneticTB

A Mathematica program package MagneticTB, which can generate
the tight-binding model for arbitrary magnetic space group. The only
input parameters in MagneticTB are the (magnetic) space group number
and the orbital information in each Wyckoff positions. Some useful
functions including getting the matrix expression for symmetry operators,
manipulating the energy band structure by parameters and interfacing
with other software are also developed.

## Installation

 Unzip the "MagneticTB-main.zip" file and copy the MagneticTB directory to any of the following four paths:

* ```FileNameJoin[{$UserBaseDirectory, "Applications"}]```
* ```FileNameJoin[{$BaseDirectory, "Applications"}]```
* ```FileNameJoin[{$InstallationDirectory, "AddOns", "Packages"}]```
* ```FileNameJoin[{$InstallationDirectory, "AddOns", "Applications"}]```


Then one can use the package after running ```Needs["MagneticTB`"]```.
The version of Mathematica should higher or equal to 11.0.

## Capabilities of MagneticTB

* Construct the tight-binding model for arbitrary magnetic space group
* Get the matrix expression for symmetry operators
* Interface with other software
* Manipulate the energy band structure by parameters
* Calculate the band co-representations of tight-binding model

See [Comput. Phys. Commun. 270, 108153 (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0010465521002654) [(arXiv:2012.08871)](https://arxiv.org/abs/2105.09504) for detail (please cite this paper if you use our code for your research).

## Examples

See examples.nb.

## Release Notes

v1.00b

* Add a new example for Charge-4 Weyl point in double magnetic space group [PRB **105**, 104426 (2022)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.104426)[(arXiv:2112.10479)](https://arxiv.org/abs/2112.10479).

v1.00c
* Add Chinese manual.
* Fix a bug for displaying the lattice vector of monoclinic lattice.

v1.01 2022/07/22
* Fixed an issue where in very rare cases the order of basis functions would change due to the automatic unitarization.

v1.02 2022/12/1
* MagneticTB can calculate the band co-representations of tight-binding model!

   * Before using this capability, [```SpaceGroupIrep```](https://github.com/goodluck1982/SpaceGroupIrep) and [```MSGCorep```](https://github.com/goodluck1982/MSGCorep)  package needed to be installed.

   * Add two functions ```getMSGElemFromMSGCorep``` and ```getTBBandCorep```, both two funcitons are depending on the ```MSGCorep``` package.

   * ```getMSGElemFromMSGCorep[{N1, N2}]``` gives the magnetic space group element from ```MSGCorep``` package, where N1.N2 is the BNS magnetic space group number.

   * ```getTBBandCorep[BNSNo, Hamiltonian, paramaters, kset]```, give the co-representations of tight-binding model, where ```BNSNo``` is the BNS magnetic space group number,  ```Hamiltonian``` is the Hamiltonian generated by MagneticTB, ```parameters``` is the parameter in the tight-binding model, ``kset`` is the list contains several k points.

   * For Orthorhombic and Monoclinic lattices, if you want to calculate the co-representations of tight-binding model please use ```getMSGElemFromMSGCorep```  to get the magnetic space group elements rather than ```msgop```. For other lattices both ```msgop``` and ```getMSGElemFromMSGCorep``` should OK.

   * See  examples.nb for concrete example.

   * Please also consider to cite ```MSGCorep``` package \([arXiv:2211.10740](https://arxiv.org/abs/2211.10740)\) if you are using this capability.


