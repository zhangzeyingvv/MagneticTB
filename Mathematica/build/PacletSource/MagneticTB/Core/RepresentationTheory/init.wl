(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "RepresentationData.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "DirectProductRepresentation.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "InducedRepresentation.wl"}]];
]
