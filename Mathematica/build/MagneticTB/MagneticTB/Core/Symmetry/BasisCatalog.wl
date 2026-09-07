(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

basisdict = <|
  "s" -> 1,
  "px" -> x,
  "py" -> y,
  "pz" -> z,
  "px+ipy" -> x + I y,
  "px-ipy" -> x - I y,
  "dx2-y2" -> x^2 - y^2,
  "dz2" -> 2 z^2 - x^2 - y^2,
  "dxy" -> 2 x y,
  "dyz" -> 2 y z,
  "dxz" -> 2 x z,

  "sup" -> {1, 0},
  "pxup" -> {x, 0},
  "pyup" -> {y, 0},
  "pzup" -> {z, 0},
  "px+ipy up" -> {x + I y, 0},
  "px-ipy up" -> {x - I y, 0},
  "dx2-y2up" -> {x^2 - y^2, 0},
  "dz2up" -> {2 z^2 - x^2 - y^2, 0},
  "dxyup" -> {2 x y, 0},
  "dyzup" -> {2 y z, 0},
  "dxzup" -> {2 x z, 0},

  "sdn" -> {0, 1},
  "pxdn" -> {0, x},
  "pydn" -> {0, y},
  "pzdn" -> {0, z},
  "px+ipy dn" -> {0, x + I y},
  "px-ipy dn" -> {0, x - I y},
  "dx2-y2dn" -> {0, x^2 - y^2},
  "dz2dn" -> {0, 2 z^2 - x^2 - y^2},
  "dxydn" -> {0, 2 x y},
  "dyzdn" -> {0, 2 y z},
  "dxzdn" -> {0, 2 x z},

  "ptest3" -> {0, (x + I y)/Sqrt[2]},
  "ptest4" -> {(x - I y)/Sqrt[2], 0}
|>;

resolveMagneticTBBasisEntry[entry_] := Module[{resolved},
  resolved = basisdict[entry];
  If[MissingQ[resolved], entry, resolved]
];

End[]

EndPackage[]
