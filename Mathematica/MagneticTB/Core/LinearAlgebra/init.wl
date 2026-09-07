(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "MatrixPredicates.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "OperatorSpaceBasis.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "LinearConstraintKernel.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "OperatorSpaceConstraints.wl"}]];
]
