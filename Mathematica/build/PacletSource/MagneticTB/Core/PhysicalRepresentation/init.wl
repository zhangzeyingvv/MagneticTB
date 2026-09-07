(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "SiteActionData.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "PhysicalRepresentationCompiler.wl"}]];
]
