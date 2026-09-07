(* ::Package:: *)

With[{moduleDirectory = DirectoryName[$InputFileName]},
  Block[{$ContextPath = DeleteCases[$ContextPath, "MagneticTB`"]},
    Get[FileNameJoin[{moduleDirectory, "NullSpace", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "LinearAlgebra", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "AbstractGroup", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "RepresentationTheory", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "Symmetry", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "CrystalGeometry", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "PhysicalRepresentation", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "TightBinding", "init.wl"}]];
    Get[FileNameJoin[{moduleDirectory, "Model", "init.wl"}]];
  ]
]
