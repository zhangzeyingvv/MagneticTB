(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  explicitSquareMatrixQ,
  nativeSplitExplicitSymmetryInformation,
  nativeSplitExplicitRepresentationInformation,
  nativeContinuousExactZeroQ,
  nativeContinuousExactMatrixEqualQ,
  nativeContinuousRecordParameter,
  nativeContinuousMatrixGenerator,
  nativeExplicitRepresentationDimensions,
  nativeExplicitOrbitalLabels,
  nativeCompileExplicitContinuousRepresentation,
  nativeValidateExplicitContinuousCompatibility,
  nativeAbstractPointOperationRecords,
  nativeAbstractInducedBasisOrderingData
];

explicitSquareMatrixQ[matrix_] :=
  MatrixPredicates`ExactSquareMatrixQ[matrix];

nativeSplitExplicitSymmetryInformation[
    symmetryInformation_Association
  ] := Module[{continuousLabels, finiteInformation, label, record},
  continuousLabels = Select[
    Keys[symmetryInformation],
    AssociationQ[symmetryInformation[#]] &&
      KeyExistsQ[symmetryInformation[#], "continuous"] &
  ];
  If[Length[continuousLabels] > 1, Return[$Failed]];
  finiteInformation = KeyDrop[symmetryInformation, continuousLabels];
  If[finiteInformation === <||>, Return[$Failed]];
  If[continuousLabels === {},
    <|"Finite" -> finiteInformation, "Continuous" -> None|>,
    label = First[continuousLabels];
    record = symmetryInformation[label];
    <|
      "Finite" -> finiteInformation,
      "Continuous" -> <|"Label" -> label, "Record" -> record|>
    |>
  ]
];

nativeSplitExplicitSymmetryInformation[symmetryInformation_List] :=
  <|"Finite" -> symmetryInformation, "Continuous" -> None|>;

nativeSplitExplicitSymmetryInformation[_] := $Failed;

nativeSplitExplicitRepresentationInformation[
    representationInformation_List
  ] := <|
  "Discrete" -> representationInformation,
  "Continuous" -> None
|>;

nativeSplitExplicitRepresentationInformation[
    representationInformation_Association
  ] := Module[{keys, discrete, continuous},
  keys = Keys[representationInformation];
  If[
    !SubsetQ[{"Discrete", "Continuous"}, keys] ||
      Complement[keys, {"Discrete", "Continuous"}] =!= {},
    Return[$Failed]
  ];
  discrete = Lookup[representationInformation, "Discrete", $Failed];
  continuous = Lookup[representationInformation, "Continuous", $Failed];
  If[!ListQ[discrete] || !AssociationQ[continuous], Return[$Failed]];
  If[Keys[continuous] =!= {"SiteMatrices"}, Return[$Failed]];
  <|"Discrete" -> discrete, "Continuous" -> continuous|>
];

nativeSplitExplicitRepresentationInformation[_] := $Failed;

nativeContinuousExactZeroQ[expression_, assumptions_] :=
  TrueQ[FullSimplify[expression, Assumptions -> assumptions] === 0] ||
    TrueQ[
      PossibleZeroQ[
        FullSimplify[expression, Assumptions -> assumptions]
      ]
    ];

nativeContinuousExactMatrixEqualQ[
    first_?MatrixQ,
    second_?MatrixQ,
    assumptions_
  ] := Dimensions[first] === Dimensions[second] &&
  And @@ (
    nativeContinuousExactZeroQ[#, assumptions] & /@
      Flatten[Normal[first - second]]
  );

nativeContinuousRecordParameter[continuousInput_Association] := Module[
  {record, parameter, space, spin},
  record = Lookup[continuousInput, "Record", $Failed];
  If[!AssociationQ[record], Return[$Failed]];
  parameter = Lookup[record, "continuous", $Failed];
  space = Lookup[record, "space", $Failed];
  spin = Lookup[record, "spin", $Failed];
  If[
    !MatchQ[parameter, _Symbol] ||
      !MatchQ[space, {_?MatrixQ, _List}] ||
      Dimensions[space[[1]]] =!= {3, 3} || Length[space[[2]]] =!= 3 ||
      !MatrixPredicates`ExactMatrixEqualQ[
        space[[1]],
        IdentityMatrix[3]
      ] ||
      !And @@ (
        MatrixPredicates`ExactZeroExpressionQ /@ space[[2]]
      ) ||
      !MatchQ[spin, {_?MatrixQ, 0}] || Dimensions[spin[[1]]] =!= {3, 3},
    $Failed,
    parameter
  ]
];

nativeContinuousMatrixGenerator[
    matrix_?explicitSquareMatrixQ,
    parameter_Symbol
  ] := Module[
  {dimension, assumptions, atIdentity, unitaryResidual,
   firstParameter, secondParameter, groupResidual, generator,
   hermitianResidual},
  dimension = Length[matrix];
  assumptions = Element[parameter, Reals];
  atIdentity = matrix /. parameter -> 0;
  If[
    !nativeContinuousExactMatrixEqualQ[
      atIdentity,
      IdentityMatrix[dimension],
      True
    ],
    Return[$Failed]
  ];
  unitaryResidual = ConjugateTranspose[matrix].matrix -
    IdentityMatrix[dimension];
  If[
    !nativeContinuousExactMatrixEqualQ[
      unitaryResidual,
      ConstantArray[0, {dimension, dimension}],
      assumptions
    ],
    Return[$Failed]
  ];
  firstParameter = Unique["continuousParameter1$"];
  secondParameter = Unique["continuousParameter2$"];
  groupResidual =
    (matrix /. parameter -> firstParameter) .
      (matrix /. parameter -> secondParameter) -
      (matrix /. parameter -> firstParameter + secondParameter);
  If[
    !nativeContinuousExactMatrixEqualQ[
      groupResidual,
      ConstantArray[0, {dimension, dimension}],
      Element[{firstParameter, secondParameter}, Reals]
    ],
    Return[$Failed]
  ];
  generator = FullSimplify[
    I D[matrix, parameter] /. parameter -> 0
  ];
  hermitianResidual = generator - ConjugateTranspose[generator];
  If[
    !nativeContinuousExactMatrixEqualQ[
      hermitianResidual,
      ConstantArray[0, {dimension, dimension}],
      True
    ],
    Return[$Failed]
  ];
  generator
];

nativeContinuousMatrixGenerator[_, _] := $Failed;

nativeCompileExplicitContinuousRepresentation[
    None,
    None,
    _List,
    _List
  ] := None;

nativeCompileExplicitContinuousRepresentation[
    continuousInput_Association,
    continuousRepresentation_Association,
    siteOrbits_List,
    localDimensions_List
  ] := Module[
  {parameter, siteMatrices, generators},
  parameter = nativeContinuousRecordParameter[continuousInput];
  siteMatrices = Lookup[
    continuousRepresentation,
    "SiteMatrices",
    $Failed
  ];
  If[
    parameter === $Failed ||
      !ListQ[siteMatrices] ||
      Length[siteMatrices] =!= Length[siteOrbits] ||
      Length[localDimensions] =!= Length[siteOrbits] ||
      !And @@ MapThread[
        Function[{matrices, orbit, dimension},
          ListQ[matrices] && Length[matrices] == Length[orbit] &&
            And @@ (
              Dimensions[#] === {dimension, dimension} & /@ matrices
            )
        ],
        {siteMatrices, siteOrbits, localDimensions}
      ],
    Return[$Failed]
  ];
  generators = Map[
    nativeContinuousMatrixGenerator[#, parameter] &,
    siteMatrices,
    {2}
  ];
  If[MemberQ[generators, $Failed, Infinity], Return[$Failed]];
  <|
    "GeneratorCount" -> 1,
    "Label" -> continuousInput["Label"],
    "SiteGenerators" -> generators,
    "RepresentationSource" -> "ExplicitMatrices",
    "Verified" -> True
  |>
];

nativeCompileExplicitContinuousRepresentation[___] := $Failed;

nativeValidateExplicitContinuousCompatibility[
    _Association,
    None
  ] := <|"NormalizerSigns" -> {}, "Verified" -> True|>;

nativeValidateExplicitContinuousCompatibility[
    prepared_Association,
    continuousData_Association
  ] := Module[
  {actionData, siteGenerators, imageSiteIndices, localBlocks,
   antiunitaryFlags, operationCount, signs, mapped, sourceGenerator,
   targetGenerator, plusQ, minusQ, nontrivialSigns},
  actionData = Lookup[prepared, "RepresentationActionData", $Failed];
  siteGenerators = Lookup[continuousData, "SiteGenerators", $Failed];
  If[!AssociationQ[actionData] || !ListQ[siteGenerators], Return[$Failed]];
  imageSiteIndices = Lookup[actionData, "ImageSiteIndices", $Failed];
  localBlocks = Lookup[actionData, "LocalBlocks", $Failed];
  antiunitaryFlags = Lookup[actionData, "AntiunitaryFlags", $Failed];
  If[
    !ListQ[imageSiteIndices] || !ListQ[localBlocks] ||
      !ListQ[antiunitaryFlags],
    Return[$Failed]
  ];
  operationCount = Length[antiunitaryFlags];
  signs = Table[
    nontrivialSigns = DeleteCases[
      Flatten@Table[
        sourceGenerator = siteGenerators[[orbit, sourceSite]];
        targetGenerator = siteGenerators[[
          orbit,
          imageSiteIndices[[orbit, operation, sourceSite]]
        ]];
        mapped = FullSimplify[
          localBlocks[[orbit, operation, sourceSite]] .
            If[
              antiunitaryFlags[[operation]],
              Conjugate[sourceGenerator],
              sourceGenerator
            ] .
            ConjugateTranspose[
              localBlocks[[orbit, operation, sourceSite]]
            ]
        ];
        plusQ = MatrixPredicates`ExactMatrixEqualQ[
          mapped,
          targetGenerator
        ];
        minusQ = MatrixPredicates`ExactMatrixEqualQ[
          mapped,
          -targetGenerator
        ];
        Which[
          plusQ && minusQ, 0,
          plusQ, 1,
          minusQ, -1,
          True, $Failed
        ],
        {orbit, Length[siteGenerators]},
        {sourceSite, Length[siteGenerators[[orbit]]]}
      ],
      0
    ];
    If[
      MemberQ[nontrivialSigns, $Failed] ||
        Length[DeleteDuplicates[nontrivialSigns]] > 1,
      $Failed,
      If[nontrivialSigns === {}, 1, First[nontrivialSigns]]
    ],
    {operation, operationCount}
  ];
  If[MemberQ[signs, $Failed], Return[$Failed]];
  <|"NormalizerSigns" -> signs, "Verified" -> True|>
];

nativeValidateExplicitContinuousCompatibility[___] := $Failed;

nativeExplicitRepresentationDimensions[
    representationInformation_List,
    "DirectProduct",
    orbitCount_Integer?Positive,
    operationCount_Integer?Positive
  ] := Module[{dimensions},
  If[
    Length[representationInformation] =!= orbitCount ||
      !And @@ (
        ListQ[#] && Length[#] == operationCount &&
          And @@ (explicitSquareMatrixQ /@ #) & /@
          representationInformation
      ),
    Return[$Failed]
  ];
  dimensions = Map[Dimensions, representationInformation, {2}];
  If[!And @@ (SameQ @@ # & /@ dimensions), Return[$Failed]];
  First[First[#]] & /@ dimensions
];

nativeExplicitRepresentationDimensions[
    representationInformation_List,
    "Induced",
    orbitCount_Integer?Positive,
    operationCount_Integer?Positive
  ] := Module[{matrixLists, dimensions},
  If[
    Length[representationInformation] =!= orbitCount ||
      !And @@ (AssociationQ /@ representationInformation),
    Return[$Failed]
  ];
  matrixLists = Lookup[
    representationInformation,
    "SiteSymmetryMatrices",
    $Failed
  ];
  If[
    !ListQ[matrixLists] || Length[matrixLists] =!= orbitCount ||
      MemberQ[matrixLists, $Failed] ||
      !And @@ (
        ListQ[#] && # =!= {} &&
          And @@ (explicitSquareMatrixQ /@ #) & /@ matrixLists
      ),
    Return[$Failed]
  ];
  dimensions = Map[Dimensions, matrixLists, {2}];
  If[!And @@ (SameQ @@ # & /@ dimensions), Return[$Failed]];
  First[First[#]] & /@ dimensions
];

nativeExplicitRepresentationDimensions[___] := $Failed;

nativeExplicitOrbitalLabels[
    Automatic,
    dimensions_List
  ] := Map[
  Function[dimension,
    Table["orb" <> ToString[index], {index, dimension}]
  ],
  dimensions
];

nativeExplicitOrbitalLabels[
    suppliedLabels_List,
    dimensions_List
  ] := Module[{labels = suppliedLabels},
  If[
    Length[labels] =!= Length[dimensions] ||
      !And @@ MapThread[
        ListQ[#1] && Length[#1] == #2 && #1 =!= {} &,
        {labels, dimensions}
      ],
    $Failed,
    labels
  ]
];

nativeExplicitOrbitalLabels[_, _] := $Failed;

nativeAbstractPointOperationRecords[
    orbitCount_Integer?Positive,
    operationLabels_List,
    antiunitaryFlags_List
  ] := Table[
  <|
    "OperationIndex" -> operationIndex,
    "OperationLabel" -> operationLabels[[operationIndex]],
    "Antiunitary" -> antiunitaryFlags[[operationIndex]],
    "RepresentationSource" -> "Matrices"
  |>,
  {orbit, orbitCount},
  {operationIndex, Length[operationLabels]}
];

nativeAbstractInducedBasisOrderingData[
    actionCompilation_Association,
    orbitalLabels_List
  ] := Module[{orbitData, records},
  orbitData = Lookup[actionCompilation, "SiteSymmetryData", $Failed];
  If[
    !ListQ[orbitData] ||
      Length[orbitData] =!= Length[orbitalLabels] ||
      !And @@ (AssociationQ /@ orbitData),
    Return[$Failed]
  ];
  records = MapThread[
    Function[{data, labels},
      Module[{representatives, transportedBases},
        representatives = Lookup[
          data,
          "CosetRepresentativeIndices",
          $Failed
        ];
        If[
          !ListQ[representatives] ||
            !And @@ (IntegerQ[#] && Positive[#] & /@ representatives),
          Return[$Failed, Module]
        ];
        transportedBases = Map[
          Function[operationIndex,
            Map[
              <|
                "AbstractOrbitalLabel" -> #,
                "TransportOperationIndex" -> operationIndex
              |> &,
              labels
            ]
          ],
          representatives
        ];
        Join[
          KeyTake[
            data,
            {"ReferenceSiteIndex", "CosetRepresentativeIndices"}
          ],
          <|"TransportedBases" -> transportedBases|>
        ]
      ]
    ],
    {orbitData, orbitalLabels}
  ];
  If[
    MemberQ[records, $Failed, Infinity],
    $Failed,
    <|"Mode" -> "Induced", "OrbitData" -> records|>
  ]
];

End[]

EndPackage[]
