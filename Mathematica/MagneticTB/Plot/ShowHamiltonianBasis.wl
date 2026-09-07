(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  proportionalBasisFunctionCoefficient,
  canonicalTransportedBasisLabel,
  displayedHamiltonianBasisState,
  hamiltonianBasisDisplayRow,
  renderHamiltonianBasisOrder,
  renderHamiltonianElementBasis
];

proportionalBasisFunctionCoefficient[transformed_, candidate_] := Module[
  {
    transformedVector, candidateVector, nonzeroPosition,
    coefficient
  },
  transformedVector = If[ListQ[transformed], transformed, {transformed}];
  candidateVector = If[ListQ[candidate], candidate, {candidate}];
  If[Length[transformedVector] =!= Length[candidateVector],
    Return[$Failed]
  ];
  nonzeroPosition = SelectFirst[
    Range[Length[candidateVector]],
    !MatrixPredicates`ExactZeroExpressionQ[
      candidateVector[[#]]
    ] &,
    Missing["NotFound"]
  ];
  If[MissingQ[nonzeroPosition], Return[$Failed]];
  coefficient = FullSimplify[
    Cancel@Together[
      transformedVector[[nonzeroPosition]]/
        candidateVector[[nonzeroPosition]]
    ]
  ];
  If[
    !FreeQ[coefficient, x | y | z] ||
      !And @@ (
        MatrixPredicates`ExactZeroExpressionQ /@
          FullSimplify[
            transformedVector - coefficient candidateVector
          ]
      ),
    $Failed,
    coefficient
  ]
];

canonicalTransportedBasisLabel[transformed_] := Module[
  {matches, label, coefficient, coefficientText},
  matches = Cases[
    Normal[basisdict],
    Rule[currentLabel_, candidate_] :>
      With[{
        currentCoefficient = proportionalBasisFunctionCoefficient[
          transformed,
          candidate
        ]
      },
        If[
          currentCoefficient === $Failed,
          Nothing,
          {currentLabel, currentCoefficient}
        ]
      ]
  ];
  If[matches === {}, Return[transformed]];
  {label, coefficient} = First[matches];
  Which[
    TrueQ[coefficient === 1],
      label,
    TrueQ[coefficient === -1],
      "-" <> label,
    TrueQ[coefficient === I],
      "I " <> label,
    TrueQ[coefficient === -I],
      "-I " <> label,
    True,
      coefficientText = ToString[coefficient, InputForm];
      coefficientText <> " " <> label
  ]
];

canonicalTransportedBasisLabel[transformed_Association] := Row[{
  Lookup[transformed, "AbstractOrbitalLabel", transformed],
  " transported by operation ",
  Lookup[transformed, "TransportOperationIndex", "?"]
}];

displayedHamiltonianBasisState[record_Association] := If[
  record["BasisConvention"] === "Induced",
  canonicalTransportedBasisLabel[record["TransportedBasisFunction"]],
  record["BasisFunction"]
];

showHamiltonianBasis::data =
  "The Hamiltonian basis order stored by init is missing or inconsistent.";
showHamiltonianBasis::index =
  "Hamiltonian row and column indices must lie between 1 and `1`; received (`2`,`3`).";

hamiltonianBasisDisplayRow[record_Association] := {
  record["HamiltonianIndex"],
  record["SiteIndex"],
  record["WyckoffOrbitIndex"],
  record["EquivalentAtomIndex"],
  record["FractionalPosition"],
  record["LocalOrbitalIndex"],
  displayedHamiltonianBasisState[record],
  If[
    record["BasisConvention"] === "Induced",
    record["BasisFunction"],
    None
  ],
  If[
    record["BasisConvention"] === "Induced",
    Row[{"Induced from atom ", record["ReferenceEquivalentAtomIndex"]}],
    "Direct product"
  ],
  If[
    record["BasisConvention"] === "Induced",
    Row[{
      record["TransportOperationIndex"],
      ": ",
      record["TransportOperationLabel"]
    }],
    None
  ]
};

renderHamiltonianBasisOrder[data_Association] := Grid[
  Join[
    {{
      Row[{
        "Hamiltonian dimension = ", data["Dimension"],
        "; H[[i,j]] = <basis i|H|basis j>"
      }],
      SpanFromLeft, SpanFromLeft, SpanFromLeft,
      SpanFromLeft, SpanFromLeft, SpanFromLeft,
      SpanFromLeft, SpanFromLeft, SpanFromLeft
    }},
    {{
      "H index",
      "Site index",
      "Wyckoff orbit",
      "Equivalent atom",
      "Fractional position",
      "Local orbital",
      "Actual basis state",
      "Induced reference input",
      "Basis convention",
      "Transport operation"
    }},
    hamiltonianBasisDisplayRow /@ data["Orbitals"]
  ],
  Frame -> All,
  Alignment -> Left
];

renderHamiltonianElementBasis[
    data_Association,
    row_Integer?Positive,
    column_Integer?Positive
  ] := Module[{rowRecord, columnRecord},
  rowRecord = data["Orbitals"][[row]];
  columnRecord = data["Orbitals"][[column]];
  Grid[
    {
      {
        Row[{"H[[", row, ",", column, "]] = <row|H|column>"}],
        SpanFromLeft, SpanFromLeft, SpanFromLeft,
        SpanFromLeft, SpanFromLeft, SpanFromLeft,
        SpanFromLeft, SpanFromLeft, SpanFromLeft,
        SpanFromLeft
      },
      {
        "Role", "H index", "Site index", "Wyckoff orbit",
        "Equivalent atom", "Fractional position", "Local orbital",
        "Actual basis state", "Induced reference input",
        "Basis convention",
        "Transport operation"
      },
      Prepend[hamiltonianBasisDisplayRow[rowRecord], "Row (bra)"],
      Prepend[
        hamiltonianBasisDisplayRow[columnRecord],
        "Column (ket)"
      ]
    },
    Frame -> All,
    Alignment -> Left
  ]
];

showHamiltonianBasis[] := Module[{data},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  data = buildHamiltonianBasisOrderData[$CurrentModelSession];
  If[data === $Failed,
    Message[showHamiltonianBasis::data];
    Return[$Failed]
  ];
  renderHamiltonianBasisOrder[data]
];

showHamiltonianBasis[
    row_Integer?Positive,
    column_Integer?Positive
  ] := Module[{data, dimension},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  data = buildHamiltonianBasisOrderData[$CurrentModelSession];
  If[data === $Failed,
    Message[showHamiltonianBasis::data];
    Return[$Failed]
  ];
  dimension = data["Dimension"];
  If[row > dimension || column > dimension,
    Message[showHamiltonianBasis::index, dimension, row, column];
    Return[$Failed]
  ];
  renderHamiltonianElementBasis[data, row, column]
];

showHamiltonianBasis[row_, column_] := Module[{data, dimension},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  data = buildHamiltonianBasisOrderData[$CurrentModelSession];
  If[data === $Failed,
    Message[showHamiltonianBasis::data];
    Return[$Failed]
  ];
  dimension = data["Dimension"];
  Message[showHamiltonianBasis::index, dimension, row, column];
  $Failed
];

End[]
EndPackage[]
