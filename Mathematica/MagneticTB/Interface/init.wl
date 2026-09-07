(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "SessionState.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "OrbitalTable.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "PublicAPI.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "BrokenSymmetryInit.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "CompatibilityFunctions.wl"}]];
]
