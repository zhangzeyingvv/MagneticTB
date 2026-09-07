(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "DataLoader.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "SymmetryData.wl"}]];
]
