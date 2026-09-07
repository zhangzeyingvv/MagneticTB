(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Quiet[
    Get[FileNameJoin[{moduleDirectory, "GroupEnumeration.wl"}]],
    {General::shdw}
  ];
  Get[FileNameJoin[{moduleDirectory, "GroupAlgebra.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "GeneratorsAndSubgroups.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "CosetsAndSchreier.wl"}]];
  Get[FileNameJoin[{moduleDirectory, "GroupActions.wl"}]];
]
