(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  nativeGenerateSpinSpaceGroupInformation,
  nativeSpinSpaceGroupCompatibilityInformation,
  nativeSpinSpaceGroupCartesianActions,
  compileSpinSpaceGroupInput
];

nativeGenerateSpinSpaceGroupInformation[
    symmetryInformation_Association
  ] := Module[
  {
    supplied, suppliedLabels, identity, operations, labels,
    position, generatedLabel
  },
  supplied = Values[symmetryInformation];
  suppliedLabels = Keys[symmetryInformation];
  identity = <|
    "space" -> {IdentityMatrix[3], {0, 0, 0}},
    "spin" -> {IdentityMatrix[3], 0}
  |>;
  operations = Quiet@Check[
    GroupEnumeration`GenerateGroup[
      supplied,
      identity,
      SymmetryAlgebra`SpinSpaceGroupElementProduct,
      SameTest -> SymmetryAlgebra`EquivalentSpinSpaceGroupElementQ
    ],
    $Failed
  ];
  If[operations === $Failed, Return[$Failed]];
  labels = MapIndexed[
    Function[{operation, index},
      position = FirstPosition[
        supplied,
        candidate_ /;
          SymmetryAlgebra`EquivalentSpinSpaceGroupElementQ[
            candidate,
            operation
          ]
      ];
      If[
        MissingQ[position],
        generatedLabel = "GeneratedSSG" <> ToString[First[index]];
        While[MemberQ[suppliedLabels, generatedLabel],
          generatedLabel = generatedLabel <> "$"
        ];
        generatedLabel,
        suppliedLabels[[First[position]]]
      ]
    ],
    operations
  ];
  If[!DuplicateFreeQ[labels], Return[$Failed]];
  AssociationThread[labels, operations]
];

nativeSpinSpaceGroupCompatibilityInformation[
    symmetryInformation_Association
  ] := MapThread[
  {
    #1,
    #2["space"][[1]],
    #2["space"][[2]],
    If[#2["spin"][[2]] === 1, "T", "F"]
  } &,
  {Keys[symmetryInformation], Values[symmetryInformation]}
];

nativeSpinSpaceGroupCartesianActions[
    elements_List,
    evaluatedLattice_?MatrixQ
  ] := Map[
  Function[element,
    FullSimplify[
      Transpose[evaluatedLattice] . element["spin"][[1]] .
        Inverse[Transpose[evaluatedLattice]]
    ]
  ],
  elements
];

compileSpinSpaceGroupInput[
    suppliedSymmetryInformation_Association,
    evaluatedLattice_?MatrixQ,
    generateGroup : (True | False)
  ] := Module[
  {spinSpaceGroupInformation, elements, groupData},
  spinSpaceGroupInformation = If[
    generateGroup,
    nativeGenerateSpinSpaceGroupInformation[suppliedSymmetryInformation],
    suppliedSymmetryInformation
  ];
  If[!AssociationQ[spinSpaceGroupInformation], Return[$Failed]];
  elements = Values[spinSpaceGroupInformation];
  groupData = SymmetryAlgebra`CompileSpinSpaceGroupData[elements];
  If[groupData === $Failed, Return[$Failed]];
  <|
    "SpinSpaceGroupInformation" -> spinSpaceGroupInformation,
    "SpinSpaceGroupElements" -> elements,
    "SymmetryInformation" ->
      nativeSpinSpaceGroupCompatibilityInformation[
        spinSpaceGroupInformation
      ],
    "GroupAlgebra" -> groupData["GroupAlgebra"],
    "SpatialActions" -> Map[
      <|
        "Rotation" -> #["space"][[1]],
        "Translation" -> #["space"][[2]]
      |> &,
      elements
    ],
    "SpinActions" -> nativeSpinSpaceGroupCartesianActions[
      elements,
      evaluatedLattice
    ],
    "AntiunitaryFlags" -> (#["spin"][[2]] === 1 & /@ elements)
  |>
];

compileSpinSpaceGroupInput[___] := $Failed;

End[]

EndPackage[]
