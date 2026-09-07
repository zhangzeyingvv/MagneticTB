(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "Common.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "BandFitting.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "VaspEigenvalues.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "CompareBand.wl"}]];
]
