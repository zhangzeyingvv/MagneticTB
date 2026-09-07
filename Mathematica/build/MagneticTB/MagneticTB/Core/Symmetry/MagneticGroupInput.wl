(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  nativeSeitzOperations,
  nativeSpatialActions,
  nativeSeitzGroupProduct,
  nativeSeitzOperationSameQ,
  nativeGenerateSymmetryInformation,
  compileMagneticGroupInput
];

nativeSeitzOperations[symmetryInformation_List] := Map[
  <|
    "Rotation" -> #[[2]],
    "Translation" -> #[[3]],
    "Antiunitary" -> (#[[4]] === "T")
  |> &,
  symmetryInformation
];

nativeSpatialActions[seitzOperations_List] :=
  KeyTake[#, {"Rotation", "Translation"}] & /@ seitzOperations;

nativeSeitzGroupProduct[left_Association, right_Association] := Module[
  {product = SymmetryAlgebra`SeitzOperationProduct[left, right]},
  <|
    "Rotation" -> product["Rotation"],
    "Translation" -> Mod[product["Translation"], 1],
    "Antiunitary" -> product["Antiunitary"]
  |>
];

nativeSeitzOperationSameQ[first_Association, second_Association] :=
  SameQ[first, second] ||
    SymmetryAlgebra`EquivalentSeitzOperationQ[first, second];

nativeGenerateSymmetryInformation[symmetryInformation_List] := Module[
  {supplied, identity, operations, labels, position},
  supplied = nativeSeitzOperations[symmetryInformation];
  identity = <|
    "Rotation" -> IdentityMatrix[3],
    "Translation" -> {0, 0, 0},
    "Antiunitary" -> False
  |>;
  operations = Quiet@Check[
    GroupEnumeration`GenerateGroup[
      supplied,
      identity,
      nativeSeitzGroupProduct,
      SameTest -> nativeSeitzOperationSameQ
    ],
    $Failed
  ];
  If[operations === $Failed, Return[$Failed]];
  labels = MapIndexed[
    Function[{operation, index},
      position = FirstPosition[
        supplied,
        candidate_ /; nativeSeitzOperationSameQ[candidate, operation]
      ];
      If[
        MissingQ[position],
        "Generated" <> ToString[First[index]],
        symmetryInformation[[First[position], 1]]
      ]
    ],
    operations
  ];
  MapThread[
    {
      #1,
      #2["Rotation"],
      #2["Translation"],
      If[#2["Antiunitary"], "T", "F"]
    } &,
    {labels, operations}
  ]
];

compileMagneticGroupInput[
    suppliedSymmetryInformation_List,
    evaluatedLattice_?MatrixQ,
    generateGroup : (True | False)
  ] := Module[
  {symmetryInformation, seitzOperations, groupData, spatialActions},
  symmetryInformation = If[
    generateGroup,
    nativeGenerateSymmetryInformation[suppliedSymmetryInformation],
    suppliedSymmetryInformation
  ];
  If[symmetryInformation === $Failed, Return[$Failed]];
  seitzOperations = nativeSeitzOperations[symmetryInformation];
  groupData = SymmetryAlgebra`CompileSeitzGroupData[seitzOperations];
  If[groupData === $Failed, Return[$Failed]];
  spatialActions = nativeSpatialActions[seitzOperations];
  <|
    "SymmetryInformation" -> symmetryInformation,
    "GroupAlgebra" -> groupData["GroupAlgebra"],
    "SpatialActions" -> spatialActions,
    "SpinActions" -> nativeMSGSpinActions[
      spatialActions,
      evaluatedLattice
    ],
    "AntiunitaryFlags" ->
      (TrueQ[Lookup[#, "Antiunitary", False]] & /@ seitzOperations)
  |>
];

compileMagneticGroupInput[___] := $Failed;

End[]

EndPackage[]
