(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

texOutput[mat_] := Module[{ptex},
  ptex = ToString[TeXForm@mat];
  StringReplace[
    ToString[ptex],
    {
      RegularExpression["\\\\text\\{(\\D)(\\d+)\\}"] -> "$1" <> "_" <> "{$2}",
      RegularExpression["\\\\text\\{(k)(.)\\}"] -> "$1" <> "_" <> "$2",
      RegularExpression["\\\\text\\{(\\D)(\\d+)(\\D)(\\d+)\\}"] ->
        "$1" <> "_" <> "{$4}" <> "^" <> "{$2}"
    }
  ]
];

End[]
EndPackage[]
