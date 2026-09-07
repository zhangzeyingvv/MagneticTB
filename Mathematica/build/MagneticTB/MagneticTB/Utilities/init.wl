(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "MatrixUtilities.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "HamiltonianUtilities.wl"}]];
]
