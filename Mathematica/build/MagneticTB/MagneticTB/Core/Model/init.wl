(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{moduleDirectory, "ModelSchema.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "ExplicitRepresentationInput.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "InducedRepresentationInput.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "ModelPreparation.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "ModelInputCompiler.wl"}]];
]
