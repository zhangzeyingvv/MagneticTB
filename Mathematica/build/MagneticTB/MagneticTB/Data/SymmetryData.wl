(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

braLatt = <|
  "CubicP" -> {
    {{a, 0, 0}, {0, a, 0}, {0, 0, a}},
    {{a, 0, 0}, {0, a, 0}, {0, 0, a}}
  },
  "CubicF" -> {
    {{a, 0, 0}, {0, a, 0}, {0, 0, a}},
    {{0, a/2, a/2}, {a/2, 0, a/2}, {a/2, a/2, 0}}
  },
  "CubicI" -> {
    {{a, 0, 0}, {0, a, 0}, {0, 0, a}},
    {{-a/2, a/2, a/2}, {a/2, -a/2, a/2}, {a/2, a/2, -a/2}}
  },
  "TetragonalP" -> {
    {{a, 0, 0}, {0, a, 0}, {0, 0, c}},
    {{a, 0, 0}, {0, a, 0}, {0, 0, c}}
  },
  "TetragonalI" -> {
    {{a, 0, 0}, {0, a, 0}, {0, 0, c}},
    {{-a/2, a/2, c/2}, {a/2, -a/2, c/2}, {a/2, a/2, -c/2}}
  },
  "OrthorhombicP" -> {
    {{a, 0, 0}, {0, b, 0}, {0, 0, c}},
    {{a, 0, 0}, {0, b, 0}, {0, 0, c}}
  },
  "OrthorhombicF" -> {
    {{a, 0, 0}, {0, b, 0}, {0, 0, c}},
    {{0, b/2, c/2}, {a/2, 0, c/2}, {a/2, b/2, 0}}
  },
  "OrthorhombicI" -> {
    {{a, 0, 0}, {0, b, 0}, {0, 0, c}},
    {{-a/2, b/2, c/2}, {a/2, -b/2, c/2}, {a/2, b/2, -c/2}}
  },
  "OrthorhombicC" -> {
    {{a, 0, 0}, {0, b, 0}, {0, 0, c}},
    {{a/2, b/2, 0}, {-a/2, b/2, 0}, {0, 0, c}}
  },
  "OrthorhombicA" -> {
    {{a, 0, 0}, {0, b, 0}, {0, 0, c}},
    {{a, 0, 0}, {0, -b/2, c/2}, {0, b/2, c/2}}
  },
  "HexagonalP" -> {
    {{a, 0, 0}, {-a/2, (Sqrt[3] a)/2, 0}, {0, 0, c}},
    {{a, 0, 0}, {-a/2, (Sqrt[3] a)/2, 0}, {0, 0, c}}
  },
  "TrigonalR" -> {
    {{Sqrt[3] a, 0, 0}, {-(Sqrt[3] a)/2, (3 a)/2, 0}, {0, 0, 3 c}},
    {{(Sqrt[3] a)/2, a/2, c}, {-(Sqrt[3] a)/2, a/2, c}, {0, -a, c}}
  },
  "MonoclinicP" -> {
    {{a, 0, 0}, {0, 0, b}, {c Cos[\[Beta]], c Sin[\[Beta]], 0}},
    {{a, 0, 0}, {0, 0, b}, {c Cos[\[Beta]], c Sin[\[Beta]], 0}}
  },
  "MonoclinicB" -> {
    {{a, 0, 0}, {0, 0, b}, {c Cos[\[Beta]], c Sin[\[Beta]], 0}},
    {{a/2, 0, b/2}, {-a/2, 0, b/2}, {c Cos[\[Beta]], c Sin[\[Beta]], 0}}
  },
  "TriclinicP" -> {
    {
      {a, 0, 0},
      {b Cos[\[Gamma]], b Sin[\[Gamma]], 0},
      {
        c Cos[\[Beta]],
        c (Cos[\[Alpha]] - Cos[\[Beta]] Cos[\[Gamma]]) Csc[\[Gamma]],
        c Sqrt[
          1 - Cos[\[Alpha]]^2 - Cos[\[Beta]]^2 +
            2 Cos[\[Alpha]] Cos[\[Beta]] Cos[\[Gamma]] - Cos[\[Gamma]]^2
        ] Csc[\[Gamma]]
      }
    },
    {
      {a, 0, 0},
      {b Cos[\[Gamma]], b Sin[\[Gamma]], 0},
      {
        c Cos[\[Beta]],
        c (Cos[\[Alpha]] - Cos[\[Beta]] Cos[\[Gamma]]) Csc[\[Gamma]],
        c Sqrt[
          1 - Cos[\[Alpha]]^2 - Cos[\[Beta]]^2 +
            2 Cos[\[Alpha]] Cos[\[Beta]] Cos[\[Gamma]] - Cos[\[Gamma]]^2
        ] Csc[\[Gamma]]
      }
    }
  }
|>;

msgop[number_] := Module[{data},
  data = MSGOP[number];
  Print["Magnetic space group (BNS): ", data["MSG"]];
  Print["Lattice: ", data["BRAV"]];
  Print["Primitive Lattice Vactor: ", braLatt[data["BRAV"]][[2]]];
  Print["Conventional Lattice Vactor: ", braLatt[data["BRAV"]][[1]]];
  data["SymmetryOperation"]
];

mlgop[mlg_] := Module[{},
  Print["Magnetic layer group:", {StringRiffle[mlg, "."], lgOgToSymbol[mlg]}];
  Print["Lattice:", lgToBrav@First[mlg]];
  Print["Primitive Lattice Vactor:", BasicVectorsMLG[lgToBrav@First[mlg]]];
  layerop[mlg]
];

mrgop[mlg_] := Module[{},
  Print["Magnetic rod group:", {StringRiffle[mlg, "."], rgOgToSymbol[mlg]}];
  Print["Lattice:", rgToBrav@First[mlg]];
  Print["Primitive Lattice Vactor:", BasicVectorsMRG[rgToBrav@First[mlg]]];
  rodop[mlg]
];

End[]
EndPackage[]
