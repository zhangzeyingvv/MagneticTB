(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "BasisCatalog.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "SymmetryAlgebra.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "SpinRepresentation.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "FunctionBasisRepresentation.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "MagneticGroupInput.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "SpinSpaceGroupInput.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "SpinSpaceGroupBasisAction.wl"}]];
]
