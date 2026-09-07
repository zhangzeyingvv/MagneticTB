(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  approximateScalarEqualQ,
  inferStandardBravaisData,
  inferStandardBravaisType,
  resolveStandardBravaisData,
  standardKPathDefinition,
  standardKPathDefinitionRecord,
  standardKPathBCLabelOverrides,
  standardKPathDisplayLabels,
  pathFromLabelSequence,
  pathFromLabelSequences,
  standardKPathData
];

standardKPath::lattice =
  "The current model does not contain a valid standardized numerical lattice from which a conventional k path can be selected.";
standardKPath::bravais =
  "The Bravais/BZ type `1` is not supported for the current standardized lattice.";

approximateScalarEqualQ[x_, y_, tolerance_] :=
  Abs[N[x - y]] <= tolerance Max[1., Abs[N[x]], Abs[N[y]]];

standardKPathDefinitionRecord[points_Association, sequences_List] := <|
  "Points" -> points,
  "Sequences" -> sequences
|>;

pathFromLabelSequence[
    points_Association,
    sequence_List,
    displayLabels_Association
  ] := Map[
  {
    {points[#[[1]]], points[#[[2]]]},
    {displayLabels[#[[1]]], displayLabels[#[[2]]]}
  } &,
  Partition[sequence, 2, 1]
];

pathFromLabelSequences[
    points_Association,
    sequences_List,
    displayLabels_Association
  ] := Flatten[
  pathFromLabelSequence[points, #, displayLabels] & /@ sequences,
  1
];

$standardBZTypeFamilies = <|
  "CUB" -> "CubicPrimitive",
  "FCC" -> "CubicFaceCentered",
  "BCC" -> "CubicBodyCentered",
  "TET" -> "TetragonalPrimitive",
  "BCT1" -> "TetragonalBodyCentered",
  "BCT2" -> "TetragonalBodyCentered",
  "ORC" -> "OrthorhombicPrimitive",
  "ORCF1" -> "OrthorhombicFaceCentered",
  "ORCF2" -> "OrthorhombicFaceCentered",
  "ORCF3" -> "OrthorhombicFaceCentered",
  "ORCI" -> "OrthorhombicBodyCentered",
  "ORCC" -> "OrthorhombicBaseCentered",
  "HEX" -> "HexagonalPrimitive",
  "RHL1" -> "RhombohedralPrimitive",
  "RHL2" -> "RhombohedralPrimitive",
  "MCL" -> "MonoclinicPrimitive",
  "MCLC1" -> "MonoclinicBaseCentered",
  "MCLC2" -> "MonoclinicBaseCentered",
  "MCLC3" -> "MonoclinicBaseCentered",
  "MCLC4" -> "MonoclinicBaseCentered",
  "MCLC5" -> "MonoclinicBaseCentered",
  "TRI1a" -> "TriclinicPrimitive",
  "TRI1b" -> "TriclinicPrimitive",
  "TRI2a" -> "TriclinicPrimitive",
  "TRI2b" -> "TriclinicPrimitive"
|>;

(* Point coordinates and path topology follow the Setyawan-Curtarolo
   definitions below. Display names follow Bradley-Cracknell whenever BC
   assigns the corresponding point, or an unambiguous representative on a
   named symmetry line. A missing or ambiguous BC point name keeps the source
   definition's label. Equivalent arms of one BC k-vector star may therefore
   intentionally share a display label. *)
standardKPathBCLabelOverrides["BCT1"] := <|
  "M" -> "Z", "Z" -> "\[CapitalLambda]", "Z1" -> "V"
|>;
standardKPathBCLabelOverrides["BCT2"] := <|"Y1" -> "U"|>;
standardKPathBCLabelOverrides["ORC"] := <|
  "T" -> "U", "U" -> "T", "X" -> "Y", "Y" -> "X"
|>;
standardKPathBCLabelOverrides["ORCF1"] := <|
  "A" -> "B", "A1" -> "B", "T" -> "Y",
  "X" -> "\[CapitalDelta]", "X1" -> "\[CapitalDelta]", "Y" -> "X"
|>;
standardKPathBCLabelOverrides["ORCF2"] := <|
  "C" -> "B", "C1" -> "D", "D" -> "C", "D1" -> "A",
  "H" -> "G", "H1" -> "H", "X" -> "Y", "Y" -> "X"
|>;
standardKPathBCLabelOverrides["ORCI"] := <|
  "R" -> "S", "S" -> "R", "X" -> "U",
  "X1" -> "\[CapitalDelta]", "Y" -> "\[CapitalSigma]",
  "Y1" -> "F", "Z" -> "X"
|>;
standardKPathBCLabelOverrides["ORCC"] := <|
  "A1" -> "E", "X" -> "\[CapitalSigma]", "X1" -> "C"
|>;
standardKPathBCLabelOverrides["MCL"] := <|
  "X" -> "Y", "Y" -> "Z", "Y1" -> "Z", "Z" -> "B"
|>;
standardKPathBCLabelOverrides[type : ("MCLC1" | "MCLC2")] := <|
  "N" -> "A", "N1" -> "V", "L" -> "M", "M" -> "L",
  "Y" -> "L", "Y1" -> "L", "Z" -> "V"
|>;
standardKPathBCLabelOverrides[type : ("MCLC3" | "MCLC4")] := <|
  "I" -> "M", "M" -> "L", "N" -> "A", "N1" -> "V",
  "X" -> "L", "Z" -> "V"
|>;
standardKPathBCLabelOverrides["MCLC5"] := <|
  "L" -> "M", "M" -> "L", "N" -> "A", "N1" -> "V",
  "X" -> "L", "Z" -> "V"
|>;
standardKPathBCLabelOverrides[type : ("TRI1a" | "TRI2a")] := <|
  "X" -> "B", "Y" -> "F", "Z" -> "G"
|>;
standardKPathBCLabelOverrides[type : ("TRI1b" | "TRI2b")] := <|
  "M" -> "G", "X" -> "F", "Y" -> "B"
|>;
standardKPathBCLabelOverrides[_] := <||>;

standardKPathDisplayLabels[bzType_String, points_Association] := Module[
  {overrides = standardKPathBCLabelOverrides[bzType]},
  AssociationMap[Lookup[overrides, #, #] &, Keys[points]]
];

inferStandardBravaisData[
    lattice_?MatrixQ,
    tolerance_?NumericQ
  ] := Module[
  {
    numericalLattice, gram, diagonal, offDiagonal, scale, closeQ,
    zeroQ, equalDiagonalQ, reciprocal, reciprocalGram,
    reciprocalCosines, lengthA, lengthB, lengthC, cosineAlpha,
    sineAlpha, criterion, cosineGamma, bzType, family, parameters,
    squaredA, squaredB, squaredC
  },
  numericalLattice = N[lattice];
  If[
    Dimensions[numericalLattice] =!= {3, 3} ||
      !MatrixQ[numericalLattice, NumericQ] ||
      Abs[Det[numericalLattice]] <= tolerance,
    Return[$Failed]
  ];
  gram = Chop[numericalLattice.Transpose[numericalLattice]];
  diagonal = Diagonal[gram];
  offDiagonal = {gram[[1, 2]], gram[[1, 3]], gram[[2, 3]]};
  scale = Max[1., Max[Abs[gram]]];
  closeQ[x_, y_] := Abs[N[x - y]] <= tolerance scale;
  zeroQ[x_] := closeQ[x, 0];
  equalDiagonalQ = closeQ[diagonal[[1]], diagonal[[2]]] &&
    closeQ[diagonal[[2]], diagonal[[3]]];
  reciprocal = nativeReciprocalLattice[numericalLattice];
  If[reciprocal === $Failed, Return[$Failed]];
  reciprocalGram = reciprocal.Transpose[reciprocal];
  reciprocalCosines = {
    reciprocalGram[[2, 3]]/Sqrt[reciprocalGram[[2, 2]] reciprocalGram[[3, 3]]],
    reciprocalGram[[1, 3]]/Sqrt[reciprocalGram[[1, 1]] reciprocalGram[[3, 3]]],
    reciprocalGram[[1, 2]]/Sqrt[reciprocalGram[[1, 1]] reciprocalGram[[2, 2]]]
  };

  Which[
    And @@ (zeroQ /@ offDiagonal),
      parameters = <|
        "A" -> Sqrt[diagonal[[1]]],
        "B" -> Sqrt[diagonal[[2]]],
        "C" -> Sqrt[diagonal[[3]]]
      |>;
      bzType = Which[
        equalDiagonalQ, "CUB",
        closeQ[diagonal[[1]], diagonal[[2]]], "TET",
        True, "ORC"
      ],

    equalDiagonalQ && And @@ (
      closeQ[#, diagonal[[1]]/2] & /@ offDiagonal
    ),
      bzType = "FCC";
      parameters = <|"A" -> Sqrt[2 diagonal[[1]]]|>,

    equalDiagonalQ && And @@ (
      closeQ[#, -diagonal[[1]]/3] & /@ offDiagonal
    ),
      bzType = "BCC";
      parameters = <|"A" -> Sqrt[4 diagonal[[1]]/3]|>,

    closeQ[diagonal[[1]], diagonal[[2]]] &&
      closeQ[offDiagonal[[1]], -diagonal[[1]]/2] &&
      zeroQ[offDiagonal[[2]]] && zeroQ[offDiagonal[[3]]],
      bzType = "HEX";
      parameters = <|
        "A" -> Sqrt[diagonal[[1]]],
        "C" -> Sqrt[diagonal[[3]]]
      |>,

    equalDiagonalQ &&
      closeQ[offDiagonal[[1]], offDiagonal[[2]]] &&
      closeQ[offDiagonal[[2]], offDiagonal[[3]]],
      cosineAlpha = offDiagonal[[1]]/diagonal[[1]];
      bzType = If[cosineAlpha > 0, "RHL1", "RHL2"];
      parameters = <|
        "A" -> Sqrt[diagonal[[1]]],
        "CosAlpha" -> cosineAlpha,
        "Alpha" -> ArcCos[cosineAlpha]
      |>,

    equalDiagonalQ && closeQ[offDiagonal[[2]], offDiagonal[[3]]] &&
      !closeQ[offDiagonal[[1]], offDiagonal[[2]]] &&
      offDiagonal[[2]] < -tolerance scale,
      squaredC = -4 offDiagonal[[2]];
      squaredA = squaredC/2 - 2 offDiagonal[[1]];
      If[Min[squaredA, squaredC] <= 0, Return[$Failed]];
      lengthA = Sqrt[squaredA]; lengthC = Sqrt[squaredC];
      bzType = If[lengthC < lengthA, "BCT1", "BCT2"];
      parameters = <|"A" -> lengthA, "C" -> lengthC|>,

    equalDiagonalQ,
      squaredA = -2 (offDiagonal[[1]] + offDiagonal[[2]]);
      squaredB = -2 (offDiagonal[[1]] + offDiagonal[[3]]);
      squaredC = -2 (offDiagonal[[2]] + offDiagonal[[3]]);
      If[Min[squaredA, squaredB, squaredC] <= 0, Return[$Failed]];
      bzType = "ORCI";
      parameters = <|
        "A" -> Sqrt[squaredA], "B" -> Sqrt[squaredB],
        "C" -> Sqrt[squaredC]
      |>,

    And @@ (# > tolerance scale & /@ offDiagonal) &&
      closeQ[diagonal[[1]], offDiagonal[[1]] + offDiagonal[[2]]] &&
      closeQ[diagonal[[2]], offDiagonal[[1]] + offDiagonal[[3]]] &&
      closeQ[diagonal[[3]], offDiagonal[[2]] + offDiagonal[[3]]],
      squaredA = 4 offDiagonal[[3]];
      squaredB = 4 offDiagonal[[2]];
      squaredC = 4 offDiagonal[[1]];
      lengthA = Sqrt[squaredA]; lengthB = Sqrt[squaredB];
      lengthC = Sqrt[squaredC];
      criterion = 1 - squaredA/squaredB - squaredA/squaredC;
      bzType = Which[
        Abs[criterion] <= tolerance, "ORCF3",
        criterion > 0, "ORCF1",
        True, "ORCF2"
      ];
      parameters = <|"A" -> lengthA, "B" -> lengthB, "C" -> lengthC|>,

    closeQ[diagonal[[1]], diagonal[[2]]] &&
      closeQ[offDiagonal[[2]], offDiagonal[[3]]],
      If[zeroQ[offDiagonal[[2]]],
        squaredA = 2 (diagonal[[1]] + offDiagonal[[1]]);
        squaredB = 2 (diagonal[[1]] - offDiagonal[[1]]);
        squaredC = diagonal[[3]];
        If[Min[squaredA, squaredB, squaredC] <= 0, Return[$Failed]];
        bzType = "ORCC";
        parameters = <|
          "A" -> Sqrt[squaredA], "B" -> Sqrt[squaredB],
          "C" -> Sqrt[squaredC]
        |>,
        squaredA = 2 (diagonal[[1]] - offDiagonal[[1]]);
        squaredB = 2 (diagonal[[1]] + offDiagonal[[1]]);
        squaredC = diagonal[[3]];
        If[Min[squaredA, squaredB, squaredC] <= 0, Return[$Failed]];
        lengthA = Sqrt[squaredA]; lengthB = Sqrt[squaredB];
        lengthC = Sqrt[squaredC];
        cosineAlpha = 2 offDiagonal[[2]]/(lengthB lengthC);
        sineAlpha = Sqrt[Max[0, 1 - cosineAlpha^2]];
        cosineGamma = reciprocalCosines[[3]];
        criterion = lengthB cosineAlpha/lengthC +
          lengthB^2 sineAlpha^2/lengthA^2;
        bzType = Which[
          cosineGamma < -tolerance, "MCLC1",
          Abs[cosineGamma] <= tolerance, "MCLC2",
          Abs[criterion - 1] <= tolerance, "MCLC4",
          criterion < 1, "MCLC3",
          True, "MCLC5"
        ];
        parameters = <|
          "A" -> lengthA, "B" -> lengthB, "C" -> lengthC,
          "CosAlpha" -> cosineAlpha, "SinAlpha" -> sineAlpha,
          "Alpha" -> ArcCos[cosineAlpha]
        |>
      ],

    zeroQ[offDiagonal[[1]]] && zeroQ[offDiagonal[[2]]] &&
      !zeroQ[offDiagonal[[3]]],
      lengthA = Sqrt[diagonal[[1]]]; lengthB = Sqrt[diagonal[[2]]];
      lengthC = Sqrt[diagonal[[3]]];
      cosineAlpha = offDiagonal[[3]]/(lengthB lengthC);
      bzType = "MCL";
      parameters = <|
        "A" -> lengthA, "B" -> lengthB, "C" -> lengthC,
        "CosAlpha" -> cosineAlpha,
        "SinAlpha" -> Sqrt[Max[0, 1 - cosineAlpha^2]],
        "Alpha" -> ArcCos[cosineAlpha]
      |>,

    True,
      bzType = Which[
        AllTrue[reciprocalCosines, # < -tolerance &] &&
          reciprocalCosines[[3]] >= Max[reciprocalCosines[[1 ;; 2]]] - tolerance,
          "TRI1a",
        AllTrue[reciprocalCosines, # > tolerance &] &&
          reciprocalCosines[[3]] <= Min[reciprocalCosines[[1 ;; 2]]] + tolerance,
          "TRI1b",
        reciprocalCosines[[1]] < -tolerance &&
          reciprocalCosines[[2]] < -tolerance &&
          Abs[reciprocalCosines[[3]]] <= tolerance,
          "TRI2a",
        reciprocalCosines[[1]] > tolerance &&
          reciprocalCosines[[2]] > tolerance &&
          Abs[reciprocalCosines[[3]]] <= tolerance,
          "TRI2b",
        True,
          "TRI"
      ];
      parameters = <|"ReciprocalCosines" -> reciprocalCosines|>
  ];
  family = Lookup[$standardBZTypeFamilies, bzType, "TriclinicPrimitive"];
  <|
    "BravaisType" -> family,
    "BZType" -> bzType,
    "Parameters" -> parameters,
    "GramMatrix" -> gram,
    "ReciprocalCosines" -> reciprocalCosines
  |>
];

inferStandardBravaisType[lattice_?MatrixQ, tolerance_?NumericQ] := Module[
  {data = inferStandardBravaisData[lattice, tolerance]},
  If[AssociationQ[data], data["BravaisType"], $Failed]
];

resolveStandardBravaisData[lattice_, Automatic, tolerance_] :=
  inferStandardBravaisData[lattice, tolerance];

resolveStandardBravaisData[lattice_, requested_String, tolerance_] := Module[
  {data = inferStandardBravaisData[lattice, tolerance]},
  If[!AssociationQ[data], Return[$Failed]];
  Which[
    requested === data["BZType"], data,
    requested === data["BravaisType"], data,
    True, $Failed
  ]
];

resolveStandardBravaisData[_, _, _] := $Failed;

standardKPathDefinition["CUB", _] := standardKPathDefinitionRecord[
  <|
    "\[CapitalGamma]" -> {0, 0, 0}, "X" -> {0, 1/2, 0},
    "M" -> {1/2, 1/2, 0}, "R" -> {1/2, 1/2, 1/2}
  |>,
  {{"\[CapitalGamma]", "X", "M", "\[CapitalGamma]", "R", "X"}, {"M", "R"}}
];

standardKPathDefinition["FCC", _] := standardKPathDefinitionRecord[
  <|
    "\[CapitalGamma]" -> {0, 0, 0}, "K" -> {3/8, 3/8, 3/4},
    "L" -> {1/2, 1/2, 1/2}, "U" -> {5/8, 1/4, 5/8},
    "W" -> {1/2, 1/4, 3/4}, "X" -> {1/2, 0, 1/2}
  |>,
  {{"\[CapitalGamma]", "X", "W", "K", "\[CapitalGamma]", "L", "U", "W", "L", "K"}, {"U", "X"}}
];

standardKPathDefinition["BCC", _] := standardKPathDefinitionRecord[
  <|
    "\[CapitalGamma]" -> {0, 0, 0}, "H" -> {1/2, -1/2, 1/2},
    "P" -> {1/4, 1/4, 1/4}, "N" -> {0, 0, 1/2}
  |>,
  {{"\[CapitalGamma]", "H", "N", "\[CapitalGamma]", "P", "H"}, {"P", "N"}}
];

standardKPathDefinition["TET", _] := standardKPathDefinitionRecord[
  <|
    "\[CapitalGamma]" -> {0, 0, 0}, "A" -> {1/2, 1/2, 1/2},
    "M" -> {1/2, 1/2, 0}, "R" -> {0, 1/2, 1/2},
    "X" -> {0, 1/2, 0}, "Z" -> {0, 0, 1/2}
  |>,
  {{"\[CapitalGamma]", "X", "M", "\[CapitalGamma]", "Z", "R", "A", "Z"}, {"X", "R"}, {"M", "A"}}
];

standardKPathDefinition["BCT1", parameters_] := Module[{a, c, eta},
  a = parameters["A"]; c = parameters["C"];
  eta = (1 + c^2/a^2)/4;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "M" -> {-1/2, 1/2, 1/2},
      "N" -> {0, 1/2, 0}, "P" -> {1/4, 1/4, 1/4},
      "X" -> {0, 0, 1/2}, "Z" -> {eta, eta, -eta},
      "Z1" -> {-eta, 1 - eta, eta}
    |>,
    {{"\[CapitalGamma]", "X", "M", "\[CapitalGamma]", "Z", "P", "N", "Z1", "M"}, {"X", "P"}}
  ]
];

standardKPathDefinition["BCT2", parameters_] := Module[
  {a, c, eta, zeta},
  a = parameters["A"]; c = parameters["C"];
  eta = (1 + a^2/c^2)/4;
  zeta = a^2/(2 c^2);
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "N" -> {0, 1/2, 0},
      "P" -> {1/4, 1/4, 1/4},
      "\[CapitalSigma]" -> {-eta, eta, eta},
      "\[CapitalSigma]1" -> {eta, 1 - eta, -eta},
      "X" -> {0, 0, 1/2}, "Y" -> {-zeta, zeta, 1/2},
      "Y1" -> {1/2, 1/2, -zeta}, "Z" -> {1/2, 1/2, -1/2}
    |>,
    {{"\[CapitalGamma]", "X", "Y", "\[CapitalSigma]", "\[CapitalGamma]", "Z", "\[CapitalSigma]1", "N", "P", "Y1", "Z"}, {"X", "P"}}
  ]
];

standardKPathDefinition["ORC", _] := standardKPathDefinitionRecord[
  <|
    "\[CapitalGamma]" -> {0, 0, 0}, "R" -> {1/2, 1/2, 1/2},
    "S" -> {1/2, 1/2, 0}, "T" -> {0, 1/2, 1/2},
    "U" -> {1/2, 0, 1/2}, "X" -> {1/2, 0, 0},
    "Y" -> {0, 1/2, 0}, "Z" -> {0, 0, 1/2}
  |>,
  {{"\[CapitalGamma]", "X", "S", "Y", "\[CapitalGamma]", "Z", "U", "R", "T", "Z"}, {"Y", "T"}, {"U", "X"}, {"S", "R"}}
];

standardKPathDefinition[type : ("ORCF1" | "ORCF3"), parameters_] := Module[
  {a, b, c, zeta, eta, sequences},
  a = parameters["A"]; b = parameters["B"]; c = parameters["C"];
  zeta = (1 + a^2/b^2 - a^2/c^2)/4;
  eta = (1 + a^2/b^2 + a^2/c^2)/4;
  sequences = If[
    type === "ORCF1",
    {{"\[CapitalGamma]", "Y", "T", "Z", "\[CapitalGamma]", "X", "A1", "Y"}, {"T", "X1"}, {"X", "A", "Z"}, {"L", "\[CapitalGamma]"}},
    {{"\[CapitalGamma]", "Y", "T", "Z", "\[CapitalGamma]", "X", "A1", "Y"}, {"X", "A", "Z"}, {"L", "\[CapitalGamma]"}}
  ];
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0},
      "A" -> {1/2, 1/2 + zeta, zeta},
      "A1" -> {1/2, 1/2 - zeta, 1 - zeta},
      "L" -> {1/2, 1/2, 1/2}, "T" -> {1, 1/2, 1/2},
      "X" -> {0, eta, eta}, "X1" -> {1, 1 - eta, 1 - eta},
      "Y" -> {1/2, 0, 1/2}, "Z" -> {1/2, 1/2, 0}
    |>,
    sequences
  ]
];

standardKPathDefinition["ORCF2", parameters_] := Module[
  {a, b, c, eta, delta, phi},
  a = parameters["A"]; b = parameters["B"]; c = parameters["C"];
  eta = (1 + a^2/b^2 - a^2/c^2)/4;
  delta = (1 + b^2/a^2 - b^2/c^2)/4;
  phi = (1 + c^2/b^2 - c^2/a^2)/4;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0},
      "C" -> {1/2, 1/2 - eta, 1 - eta},
      "C1" -> {1/2, 1/2 + eta, eta},
      "D" -> {1/2 - delta, 1/2, 1 - delta},
      "D1" -> {1/2 + delta, 1/2, delta},
      "L" -> {1/2, 1/2, 1/2},
      "H" -> {1 - phi, 1/2 - phi, 1/2},
      "H1" -> {phi, 1/2 + phi, 1/2},
      "X" -> {0, 1/2, 1/2}, "Y" -> {1/2, 0, 1/2},
      "Z" -> {1/2, 1/2, 0}
    |>,
    {{"\[CapitalGamma]", "Y", "C", "D", "X", "\[CapitalGamma]", "Z", "D1", "H", "C"}, {"C1", "Z"}, {"X", "H1"}, {"H", "Y"}, {"L", "\[CapitalGamma]"}}
  ]
];

standardKPathDefinition["ORCI", parameters_] := Module[
  {a, b, c, zeta, delta, eta, mu},
  a = parameters["A"]; b = parameters["B"]; c = parameters["C"];
  zeta = (1 + a^2/c^2)/4;
  delta = (b^2 - a^2)/(4 c^2);
  eta = (1 + b^2/c^2)/4;
  mu = (a^2 + b^2)/(4 c^2);
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0},
      "L" -> {-mu, mu, 1/2 - delta},
      "L1" -> {mu, -mu, 1/2 + delta},
      "L2" -> {1/2 - delta, 1/2 + delta, -mu},
      "R" -> {0, 1/2, 0}, "S" -> {1/2, 0, 0},
      "T" -> {0, 0, 1/2}, "W" -> {1/4, 1/4, 1/4},
      "X" -> {-zeta, zeta, zeta},
      "X1" -> {zeta, 1 - zeta, -zeta},
      "Y" -> {eta, -eta, eta},
      "Y1" -> {1 - eta, eta, -eta},
      "Z" -> {1/2, 1/2, -1/2}
    |>,
    {{"\[CapitalGamma]", "X", "L", "T", "W", "R", "X1", "Z", "\[CapitalGamma]", "Y", "S", "W"}, {"L1", "Y"}, {"Y1", "Z"}}
  ]
];

standardKPathDefinition["ORCC", parameters_] := Module[{a, b, zeta},
  a = parameters["A"]; b = parameters["B"];
  zeta = (1 + a^2/b^2)/4;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "A" -> {zeta, zeta, 1/2},
      "A1" -> {-zeta, 1 - zeta, 1/2}, "R" -> {0, 1/2, 1/2},
      "S" -> {0, 1/2, 0}, "T" -> {-1/2, 1/2, 1/2},
      "X" -> {zeta, zeta, 0}, "X1" -> {-zeta, 1 - zeta, 0},
      "Y" -> {-1/2, 1/2, 0}, "Z" -> {0, 0, 1/2}
    |>,
    {{"\[CapitalGamma]", "X", "S", "R", "A", "Z", "\[CapitalGamma]", "Y", "X1", "A1", "T", "Y"}, {"Z", "T"}}
  ]
];

standardKPathDefinition["HEX", _] := standardKPathDefinitionRecord[
  <|
    "\[CapitalGamma]" -> {0, 0, 0}, "A" -> {0, 0, 1/2},
    "H" -> {1/3, 1/3, 1/2}, "K" -> {1/3, 1/3, 0},
    "L" -> {1/2, 0, 1/2}, "M" -> {1/2, 0, 0}
  |>,
  {{"\[CapitalGamma]", "M", "K", "\[CapitalGamma]", "A", "L", "H", "A"}, {"L", "M"}, {"K", "H"}}
];

standardKPathDefinition["RHL1", parameters_] := Module[
  {cosineAlpha, eta, nu},
  cosineAlpha = parameters["CosAlpha"];
  eta = (1 + 4 cosineAlpha)/(2 + 4 cosineAlpha);
  nu = 3/4 - eta/2;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0},
      "B" -> {eta, 1/2, 1 - eta},
      "B1" -> {1/2, 1 - eta, eta - 1},
      "F" -> {1/2, 1/2, 0}, "L" -> {1/2, 0, 0},
      "L1" -> {0, 0, -1/2}, "P" -> {eta, nu, nu},
      "P1" -> {1 - nu, 1 - nu, 1 - eta},
      "P2" -> {nu, nu, eta - 1}, "Q" -> {1 - nu, nu, 0},
      "X" -> {nu, 0, -nu}, "Z" -> {1/2, 1/2, 1/2}
    |>,
    {{"\[CapitalGamma]", "L", "B1"}, {"B", "Z", "\[CapitalGamma]", "X"}, {"Q", "F", "P1", "Z"}, {"L", "P"}}
  ]
];

standardKPathDefinition["RHL2", parameters_] := Module[
  {alpha, eta, nu},
  alpha = parameters["Alpha"];
  eta = 1/(2 Tan[alpha/2]^2);
  nu = 3/4 - eta/2;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "F" -> {1/2, -1/2, 0},
      "L" -> {1/2, 0, 0}, "P" -> {1 - nu, -nu, 1 - nu},
      "P1" -> {nu, nu - 1, nu - 1}, "Q" -> {eta, eta, eta},
      "Q1" -> {1 - eta, -eta, -eta}, "Z" -> {1/2, -1/2, 1/2}
    |>,
    {{"\[CapitalGamma]", "P", "Z", "Q", "\[CapitalGamma]", "F", "P1", "Q1", "L", "Z"}}
  ]
];

standardKPathDefinition["MCL", parameters_] := Module[
  {b, c, cosineAlpha, sineAlpha, eta, nu},
  b = parameters["B"]; c = parameters["C"];
  cosineAlpha = parameters["CosAlpha"];
  sineAlpha = parameters["SinAlpha"];
  eta = (1 - b cosineAlpha/c)/(2 sineAlpha^2);
  nu = 1/2 - eta c cosineAlpha/b;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "A" -> {1/2, 1/2, 0},
      "C" -> {0, 1/2, 1/2}, "D" -> {1/2, 0, 1/2},
      "D1" -> {1/2, 0, -1/2}, "E" -> {1/2, 1/2, 1/2},
      "H" -> {0, eta, 1 - nu}, "H1" -> {0, 1 - eta, nu},
      "H2" -> {0, eta, -nu}, "M" -> {1/2, eta, 1 - nu},
      "M1" -> {1/2, 1 - eta, nu}, "M2" -> {1/2, eta, -nu},
      "X" -> {0, 1/2, 0}, "Y" -> {0, 0, 1/2},
      "Y1" -> {0, 0, -1/2}, "Z" -> {1/2, 0, 0}
    |>,
    {{"\[CapitalGamma]", "Y", "H", "C", "E", "M1", "A", "X", "H1"}, {"M", "D", "Z"}, {"Y", "D"}}
  ]
];

standardKPathDefinition[type : ("MCLC1" | "MCLC2"), parameters_] := Module[
  {a, b, c, cosineAlpha, sineAlpha, zeta, eta, psi, phi, sequences},
  a = parameters["A"]; b = parameters["B"]; c = parameters["C"];
  cosineAlpha = parameters["CosAlpha"];
  sineAlpha = parameters["SinAlpha"];
  zeta = (2 - b cosineAlpha/c)/(4 sineAlpha^2);
  eta = 1/2 + 2 zeta c cosineAlpha/b;
  psi = 3/4 - a^2/(4 b^2 sineAlpha^2);
  phi = psi + (3/4 - psi) b cosineAlpha/c;
  sequences = If[
    type === "MCLC1",
    {{"\[CapitalGamma]", "Y", "F", "L", "I"}, {"I1", "Z", "F1"}, {"Y", "X1"}, {"X", "\[CapitalGamma]", "N"}, {"M", "\[CapitalGamma]"}},
    {{"\[CapitalGamma]", "Y", "F", "L", "I"}, {"I1", "Z", "F1"}, {"N", "\[CapitalGamma]", "M"}}
  ];
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "N" -> {1/2, 0, 0},
      "N1" -> {0, -1/2, 0}, "F" -> {1 - zeta, 1 - zeta, 1 - eta},
      "F1" -> {zeta, zeta, eta}, "F2" -> {-zeta, -zeta, 1 - eta},
      "F3" -> {1 - zeta, -zeta, 1 - eta},
      "I" -> {phi, 1 - phi, 1/2}, "I1" -> {1 - phi, phi - 1, 1/2},
      "L" -> {1/2, 1/2, 1/2}, "M" -> {1/2, 0, 1/2},
      "X" -> {1 - psi, psi - 1, 0}, "X1" -> {psi, 1 - psi, 0},
      "X2" -> {psi - 1, -psi, 0}, "Y" -> {1/2, 1/2, 0},
      "Y1" -> {-1/2, -1/2, 0}, "Z" -> {0, 0, 1/2}
    |>,
    sequences
  ]
];

standardKPathDefinition[type : ("MCLC3" | "MCLC4"), parameters_] := Module[
  {a, b, c, cosineAlpha, sineAlpha, mu, delta, zeta, eta, phi, psi, sequences},
  a = parameters["A"]; b = parameters["B"]; c = parameters["C"];
  cosineAlpha = parameters["CosAlpha"];
  sineAlpha = parameters["SinAlpha"];
  mu = (1 + b^2/a^2)/4;
  delta = b c cosineAlpha/(2 a^2);
  zeta = mu - 1/4 + (1 - b cosineAlpha/c)/(4 sineAlpha^2);
  eta = 1/2 + 2 zeta c cosineAlpha/b;
  phi = 1 + zeta - 2 mu;
  psi = eta - 2 delta;
  sequences = If[
    type === "MCLC3",
    {{"\[CapitalGamma]", "Y", "F", "H", "Z", "I", "F1"}, {"H1", "Y1", "X", "\[CapitalGamma]", "N"}, {"M", "\[CapitalGamma]"}},
    {{"\[CapitalGamma]", "Y", "F", "H", "Z", "I"}, {"H1", "Y1", "X", "\[CapitalGamma]", "N"}, {"M", "\[CapitalGamma]"}}
  ];
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "F" -> {1 - phi, 1 - phi, 1 - psi},
      "F1" -> {phi, phi, psi}, "F2" -> {1 - phi, -phi, 1 - psi},
      "H" -> {zeta, zeta, eta}, "H1" -> {1 - zeta, -zeta, 1 - eta},
      "H2" -> {-zeta, -zeta, 1 - eta}, "I" -> {1/2, -1/2, 1/2},
      "M" -> {1/2, 0, 1/2}, "N" -> {1/2, 0, 0},
      "N1" -> {0, -1/2, 0}, "X" -> {1/2, -1/2, 0},
      "Y" -> {mu, mu, delta}, "Y1" -> {1 - mu, -mu, -delta},
      "Y2" -> {-mu, -mu, -delta}, "Y3" -> {mu, mu - 1, delta},
      "Z" -> {0, 0, 1/2}
    |>,
    sequences
  ]
];

standardKPathDefinition["MCLC5", parameters_] := Module[
  {
    a, b, c, cosineAlpha, sineAlpha, zeta, eta, mu, nu,
    omega, delta, rho
  },
  a = parameters["A"]; b = parameters["B"]; c = parameters["C"];
  cosineAlpha = parameters["CosAlpha"];
  sineAlpha = parameters["SinAlpha"];
  zeta = (b^2/a^2 + (1 - b cosineAlpha/c)/sineAlpha^2)/4;
  eta = 1/2 + 2 zeta c cosineAlpha/b;
  mu = eta/2 + b^2/(4 a^2) - b c cosineAlpha/(2 a^2);
  nu = 2 mu - zeta;
  omega = (4 nu - 1 - b^2 sineAlpha^2/a^2) c/(2 b cosineAlpha);
  delta = zeta c cosineAlpha/b + omega/2 - 1/4;
  rho = 1 - zeta a^2/b^2;
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "F" -> {nu, nu, omega},
      "F1" -> {1 - nu, 1 - nu, 1 - omega},
      "F2" -> {nu, nu - 1, omega}, "H" -> {zeta, zeta, eta},
      "H1" -> {1 - zeta, -zeta, 1 - eta},
      "H2" -> {-zeta, -zeta, 1 - eta},
      "I" -> {rho, 1 - rho, 1/2},
      "I1" -> {1 - rho, rho - 1, 1/2}, "L" -> {1/2, 1/2, 1/2},
      "M" -> {1/2, 0, 1/2}, "N" -> {1/2, 0, 0},
      "N1" -> {0, -1/2, 0}, "X" -> {1/2, -1/2, 0},
      "Y" -> {mu, mu, delta}, "Y1" -> {1 - mu, -mu, -delta},
      "Y2" -> {-mu, -mu, -delta}, "Y3" -> {mu, mu - 1, delta},
      "Z" -> {0, 0, 1/2}
    |>,
    {{"\[CapitalGamma]", "Y", "F", "L", "I"}, {"I1", "Z", "H", "F1"}, {"H1", "Y1", "X", "\[CapitalGamma]", "N"}, {"M", "\[CapitalGamma]"}}
  ]
];

standardKPathDefinition[type : ("TRI1a" | "TRI2a"), _] :=
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "L" -> {1/2, 1/2, 0},
      "M" -> {0, 1/2, 1/2}, "N" -> {1/2, 0, 1/2},
      "R" -> {1/2, 1/2, 1/2}, "X" -> {1/2, 0, 0},
      "Y" -> {0, 1/2, 0}, "Z" -> {0, 0, 1/2}
    |>,
    {{"X", "\[CapitalGamma]", "Y"}, {"L", "\[CapitalGamma]", "Z"}, {"N", "\[CapitalGamma]", "M"}, {"R", "\[CapitalGamma]"}}
  ];

standardKPathDefinition[type : ("TRI1b" | "TRI2b"), _] :=
  standardKPathDefinitionRecord[
    <|
      "\[CapitalGamma]" -> {0, 0, 0}, "L" -> {1/2, -1/2, 0},
      "M" -> {0, 0, 1/2}, "N" -> {-1/2, -1/2, 1/2},
      "R" -> {0, -1/2, 1/2}, "X" -> {0, -1/2, 0},
      "Y" -> {1/2, 0, 0}, "Z" -> {-1/2, 0, 1/2}
    |>,
    {{"X", "\[CapitalGamma]", "Y"}, {"L", "\[CapitalGamma]", "Z"}, {"N", "\[CapitalGamma]", "M"}, {"R", "\[CapitalGamma]"}}
  ];

standardKPathDefinition[_, _] := $Failed;

Options[standardKPath] = {
  "BravaisType" -> Automatic,
  "Tolerance" -> 10^-6
};

standardKPathData[OptionsPattern[standardKPath]] := Module[
  {
    session, model, metadata, lattice, tolerance, classification,
    definition, displayLabels, path
  },
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  model = Lookup[session, "ModelSpecification", $Failed];
  metadata = If[AssociationQ[model], Lookup[model, "Metadata", $Failed], $Failed];
  lattice = If[AssociationQ[metadata], Lookup[metadata, "LatticePlot", $Failed], $Failed];
  tolerance = OptionValue["Tolerance"];
  If[
    !MatrixQ[lattice, NumericQ] || Dimensions[lattice] =!= {3, 3} ||
      !NumericQ[tolerance] || tolerance <= 0,
    Message[standardKPath::lattice];
    Return[$Failed]
  ];
  classification = resolveStandardBravaisData[
    lattice,
    OptionValue["BravaisType"],
    tolerance
  ];
  If[!AssociationQ[classification],
    Message[standardKPath::bravais, OptionValue["BravaisType"]];
    Return[$Failed]
  ];
  definition = standardKPathDefinition[
    classification["BZType"],
    classification["Parameters"]
  ];
  If[!AssociationQ[definition],
    Message[standardKPath::bravais, classification["BZType"]];
    Return[$Failed]
  ];
  displayLabels = standardKPathDisplayLabels[
    classification["BZType"],
    definition["Points"]
  ];
  path = pathFromLabelSequences[
    definition["Points"],
    definition["Sequences"],
    displayLabels
  ];
  <|
    "BravaisType" -> classification["BravaisType"],
    "BZType" -> classification["BZType"],
    "Parameters" -> classification["Parameters"],
    "Points" -> definition["Points"],
    "PointLabels" -> displayLabels,
    "LabelConvention" ->
      "Bradley-CracknellWithSetyawan-CurtaroloFallback",
    "Sequences" -> definition["Sequences"],
    "Sequence" -> definition["Sequences"],
    "Path" -> path
  |>
];

standardKPath[opts : OptionsPattern[]] := Module[{data},
  data = standardKPathData[opts];
  If[AssociationQ[data], data["Path"], $Failed]
];

End[]
EndPackage[]
