(* ::Package:: *)

BeginPackage["ModelSchema`"]

CreateModelSpecification::usage =
  "CreateModelSpecification[spec] validates and freezes a mode-independent tight-binding model. spec[\"Symmetry\"] contains GroupAlgebra, SpatialActions, SpinActions, and AntiunitaryFlags in one common group-element order.";
ModelSpecificationQ::usage =
  "ModelSpecificationQ[model] tests whether model is a validated ModelSpecification.";
CreatePreparedModel::usage =
  "CreatePreparedModel[model,actionCompilation] combines a ModelSpecification with verified representation action data and returns a clean PreparedModel.";
PreparedModelQ::usage =
  "PreparedModelQ[model] tests whether model is a validated PreparedModel.";
ReplaceModelBondClasses::usage =
  "ReplaceModelBondClasses[model,bondClasses] returns a revalidated ModelSpecification with new periodic bond classes.";

Begin["`Private`"]

ClearAll[
  validLocalBasesQ,
  validPointOperationRecordsQ,
  validContinuousSymmetryQ,
  validSpinActionsQ,
  validFiniteSymmetryLawQ,
  validBondClassesQ
];

ClearAll[modelSpecificationValidationCache, preparedModelValidationCache];
modelSpecificationValidationCache = <||>;
preparedModelValidationCache = <||>;

validLocalBasesQ[bases_, orbitCount_Integer] :=
  ListQ[bases] && Length[bases] == orbitCount &&
    And @@ (ListQ[#] && Length[#] > 0 & /@ bases);

validPointOperationRecordsQ[records_, orbitCount_Integer, operationCount_Integer] :=
  ListQ[records] && Length[records] == orbitCount &&
    And @@ (
      ListQ[#] && Length[#] == operationCount &&
        And @@ (AssociationQ /@ #) & /@ records
    );

validContinuousSymmetryQ[
    None,
    _List,
    _List
  ] := True;

validContinuousSymmetryQ[
    data_Association,
    siteOrbits_List,
    localDimensions_List
  ] := Module[{siteGenerators},
  siteGenerators = Lookup[data, "SiteGenerators", $Failed];
  TrueQ[Lookup[data, "GeneratorCount", 0] === 1] &&
    TrueQ[Lookup[data, "Verified", False]] &&
    ListQ[siteGenerators] &&
    Length[siteGenerators] === Length[siteOrbits] &&
    Length[localDimensions] === Length[siteOrbits] &&
    And @@ MapThread[
      Function[{generators, orbit, dimension},
        ListQ[generators] && Length[generators] === Length[orbit] &&
          And @@ (
            MatrixPredicates`ExactSquareMatrixQ[#] &&
              Dimensions[#] === {dimension, dimension} &&
              MatrixPredicates`ExactMatrixEqualQ[
                #,
                ConjugateTranspose[#]
              ] & /@ generators
          )
      ],
      {siteGenerators, siteOrbits, localDimensions}
    ]
];

validContinuousSymmetryQ[_, _, _] := False;

validSpinActionsQ[actions_, operationCount_Integer] :=
  ListQ[actions] && Length[actions] === operationCount &&
    And @@ (
      MatrixQ[#] && Dimensions[#] === {3, 3} &&
        MatrixPredicates`ExactUnitaryMatrixQ[#] &&
        MatrixPredicates`ExactZeroExpressionQ[
          FullSimplify[Det[#] - 1]
        ] & /@ actions
    );

validFiniteSymmetryLawQ[
    spatialActions_List,
    spinActions_List,
    antiunitaryFlags_List,
    groupAlgebra_Association
  ] := Module[
  {productTable, spatialProduct, spatialDifference, productIndex},
  productTable = groupAlgebra["ProductTable"];
  And @@ Flatten@Table[
    productIndex = productTable[[left, right]];
    spatialProduct = SymmetryAlgebra`SpatialOperationProduct[
      spatialActions[[left]],
      spatialActions[[right]]
    ];
    spatialDifference = SymmetryAlgebra`SpatialOperationDifference[
      spatialProduct,
      spatialActions[[productIndex]]
    ];
    spatialDifference =!= $Failed &&
      MatrixPredicates`IntegerVectorQ[spatialDifference] &&
      MatrixPredicates`ExactMatrixEqualQ[
        spinActions[[left]].spinActions[[right]],
        spinActions[[productIndex]]
      ] &&
      Xor[
        antiunitaryFlags[[left]],
        antiunitaryFlags[[right]]
      ] === antiunitaryFlags[[productIndex]],
    {left, groupAlgebra["Order"]},
    {right, groupAlgebra["Order"]}
  ]
];

validBondClassesQ[bondClasses_] := ListQ[bondClasses] &&
  And @@ Map[
    Function[shell,
      ListQ[shell] && shell =!= {} && And @@ Map[
        MatchQ[#, {_?NumericQ, _Integer?NonNegative,
            {{_List, _List} ...}}] && #[[2]] === Length[#[[3]]] &,
        shell
      ]
    ],
    bondClasses
  ];

CreateModelSpecification::spec =
  "The supplied Association is not a valid ModelSpecification input.";

CreateModelSpecification[spec_Association] := Module[
  {
    lattice, siteOrbits, localBases, variables,
    pointOperationRecords, symmetry, spatialActions, spinActions,
    groupAlgebra, antiunitaryFlags, continuousSymmetry, bondClasses,
    orbitCount, operationCount, dimension, metadata, candidate
  },
  lattice = Lookup[spec, "Lattice", $Failed];
  siteOrbits = Lookup[spec, "SiteOrbits", $Failed];
  localBases = Lookup[spec, "LocalBases", $Failed];
  variables = Lookup[spec, "Variables", $Failed];
  pointOperationRecords = Lookup[spec, "PointOperationRecords", $Failed];
  symmetry = Lookup[spec, "Symmetry", $Failed];
  If[!AssociationQ[symmetry],
    Message[CreateModelSpecification::spec];
    Return[$Failed]
  ];
  spatialActions = Lookup[symmetry, "SpatialActions", $Failed];
  spinActions = Lookup[symmetry, "SpinActions", $Failed];
  groupAlgebra = Lookup[symmetry, "GroupAlgebra", $Failed];
  antiunitaryFlags = Lookup[symmetry, "AntiunitaryFlags", $Failed];
  continuousSymmetry = Lookup[spec, "ContinuousSymmetry", None];
  bondClasses = Lookup[spec, "BondClasses", $Failed];
  metadata = Lookup[spec, "Metadata", <||>];

  If[
    !MatrixPredicates`ExactSquareMatrixQ[lattice] ||
      !SitePermutationCompiler`ValidSiteOrbitsQ[siteOrbits] ||
      !ListQ[variables] || Length[variables] =!= Length[lattice] ||
      !ListQ[spatialActions] || spatialActions === {} ||
      !GroupAlgebra`GroupAlgebraQ[groupAlgebra] ||
      !ListQ[antiunitaryFlags] ||
      !AssociationQ[metadata],
    Message[CreateModelSpecification::spec];
    Return[$Failed]
  ];
  orbitCount = Length[siteOrbits];
  operationCount = Length[spatialActions];
  dimension = Length[lattice];
  If[
    !validLocalBasesQ[localBases, orbitCount] ||
      !validPointOperationRecordsQ[
        pointOperationRecords,
        orbitCount,
        operationCount
      ] ||
      !validContinuousSymmetryQ[
        continuousSymmetry,
        siteOrbits,
        Length /@ localBases
      ] ||
      !validSpinActionsQ[spinActions, operationCount] ||
      !And @@ (
        SymmetryAlgebra`ValidSpatialOperationQ[#, dimension] & /@
          spatialActions
      ) ||
      groupAlgebra["Order"] =!= operationCount ||
      Length[antiunitaryFlags] =!= operationCount ||
      !And @@ (BooleanQ /@ antiunitaryFlags) ||
      !validFiniteSymmetryLawQ[
        spatialActions,
        spinActions,
        antiunitaryFlags,
        groupAlgebra
      ] ||
      !validBondClassesQ[bondClasses],
    Message[CreateModelSpecification::spec];
    Return[$Failed]
  ];

  candidate = <|
    "Schema" -> "ModelSpecification",
    "SchemaVersion" -> 1,
    "Lattice" -> lattice,
    "SiteOrbits" -> siteOrbits,
    "LocalBases" -> localBases,
    "LocalDimensions" -> (Length /@ localBases),
    "Variables" -> variables,
    "PointOperationRecords" -> pointOperationRecords,
    "Symmetry" -> <|
      "GroupAlgebra" -> groupAlgebra,
      "SpatialActions" -> spatialActions,
      "SpinActions" -> spinActions,
      "AntiunitaryFlags" -> antiunitaryFlags
    |>,
    "ContinuousSymmetry" -> continuousSymmetry,
    "BondClasses" -> bondClasses,
    "Metadata" -> metadata,
    "Verified" -> True
  |>;
  AssociateTo[modelSpecificationValidationCache, Hash[candidate] -> True];
  candidate
];
CreateModelSpecification[_] :=
  (Message[CreateModelSpecification::spec]; $Failed);

ModelSpecificationQ[model_] := Module[{cached, input, canonical, result},
  cached = Lookup[
    modelSpecificationValidationCache, Hash[model], Missing["NotCached"]];
  If[!MissingQ[cached], Return[cached]];
  result = If[
    !AssociationQ[model] || Lookup[model, "Schema", None] =!= "ModelSpecification" ||
      Lookup[model, "SchemaVersion", None] =!= 1 ||
      !TrueQ[Lookup[model, "Verified", False]],
    False,
    input = KeyDrop[model, {"Schema", "SchemaVersion", "Verified", "LocalDimensions"}];
    canonical = Quiet[CreateModelSpecification[input]];
    AssociationQ[canonical] &&
      And @@ (KeyExistsQ[model, #] && model[#] === canonical[#] & /@
        Keys[canonical])
  ];
  result = TrueQ[result];
  AssociateTo[modelSpecificationValidationCache, Hash[model] -> result];
  result
];

CreatePreparedModel::model =
  "The first argument must be a validated ModelSpecification.";
CreatePreparedModel::action =
  "The action compilation is incompatible with the supplied ModelSpecification.";

CreatePreparedModel[
    model_Association,
    actionCompilation_Association
  ] := Module[
  {actionData, localDimensions, localRepresentations,
   representationMatrices, orbitRepresentationMatrices, method,
   continuousSymmetry, candidate},
  If[!ModelSpecificationQ[model],
    Message[CreatePreparedModel::model];
    Return[$Failed]
  ];
  actionData = Lookup[actionCompilation, "ActionData", $Failed];
  localDimensions = Lookup[actionCompilation, "LocalDimensions", $Failed];
  localRepresentations = Lookup[
    actionCompilation,
    "LocalRepresentationMatrices",
    Missing["NotAvailable"]
  ];
  representationMatrices = Lookup[
    actionCompilation, "RepresentationMatrices", $Failed];
  orbitRepresentationMatrices = Lookup[
    actionCompilation, "OrbitRepresentationMatrices", $Failed];
  method = Lookup[actionCompilation, "Method", "Unspecified"];
  continuousSymmetry = Lookup[model, "ContinuousSymmetry", None];
  If[
    !RepresentationActionData`SiteActionDataQ[actionData] ||
      Lookup[actionData, "SiteOrbits", $Failed] =!= model["SiteOrbits"] ||
      Lookup[actionData, "SpatialActions", $Failed] =!=
        model["Symmetry", "SpatialActions"] ||
      Lookup[actionData, "GroupAlgebra", $Failed] =!=
        model["Symmetry", "GroupAlgebra"] ||
      Lookup[actionData, "AntiunitaryFlags", $Failed] =!=
        model["Symmetry", "AntiunitaryFlags"] ||
      localDimensions =!= model["LocalDimensions"] ||
      !ListQ[representationMatrices] ||
      Length[representationMatrices] =!= model["Symmetry", "GroupAlgebra", "Order"] ||
      !And @@ (MatrixPredicates`ExactSquareMatrixQ /@ representationMatrices),
    Message[CreatePreparedModel::action];
    Return[$Failed]
  ];
  candidate = <|
    "Schema" -> "PreparedModel",
    "SchemaVersion" -> 1,
    "ModelSpecification" -> model,
    "Lattice" -> model["Lattice"],
    "SiteOrbits" -> model["SiteOrbits"],
    "Symmetry" -> model["Symmetry"],
    "ContinuousSymmetry" -> continuousSymmetry,
    "BondClasses" -> model["BondClasses"],
    "LocalDimensions" -> localDimensions,
    "RepresentationMethod" -> method,
    "RepresentationActionData" -> actionData,
    "LocalRepresentationMatrices" -> localRepresentations,
    "RepresentationMatrices" -> representationMatrices,
    "OrbitRepresentationMatrices" -> orbitRepresentationMatrices,
    "Metadata" -> model["Metadata"],
    "Verified" -> True
  |>;
  AssociateTo[preparedModelValidationCache, Hash[candidate] -> True];
  candidate
];
CreatePreparedModel[_, _] :=
  (Message[CreatePreparedModel::model]; $Failed);

PreparedModelQ[model_] := Module[
  {cached, specification, actionCompilation, canonical, result},
  cached = Lookup[
    preparedModelValidationCache, Hash[model], Missing["NotCached"]];
  If[!MissingQ[cached], Return[cached]];
  result = If[
    !AssociationQ[model] || Lookup[model, "Schema", None] =!= "PreparedModel" ||
      Lookup[model, "SchemaVersion", None] =!= 1 ||
      !TrueQ[Lookup[model, "Verified", False]],
    False,
    specification = Lookup[model, "ModelSpecification", $Failed];
    actionCompilation = <|
      "ActionData" -> Lookup[model, "RepresentationActionData", $Failed],
      "LocalDimensions" -> Lookup[model, "LocalDimensions", $Failed],
      "LocalRepresentationMatrices" ->
        Lookup[model, "LocalRepresentationMatrices", Missing["NotAvailable"]],
      "RepresentationMatrices" -> Lookup[model, "RepresentationMatrices", $Failed],
      "OrbitRepresentationMatrices" ->
        Lookup[model, "OrbitRepresentationMatrices", $Failed],
      "Method" -> Lookup[model, "RepresentationMethod", "Unspecified"]
    |>;
    canonical = Quiet[CreatePreparedModel[specification, actionCompilation]];
    AssociationQ[canonical] &&
      And @@ (KeyExistsQ[model, #] && model[#] === canonical[#] & /@
        Keys[canonical])
  ];
  result = TrueQ[result];
  AssociateTo[preparedModelValidationCache, Hash[model] -> result];
  result
];

ReplaceModelBondClasses::model =
  "The first argument must be a validated ModelSpecification.";
ReplaceModelBondClasses[
    model_Association,
    bondClasses_List
  ] := Module[{input},
  If[!ModelSpecificationQ[model],
    Message[ReplaceModelBondClasses::model];
    Return[$Failed]
  ];
  input = KeyDrop[model, {"Schema", "SchemaVersion", "Verified"}];
  CreateModelSpecification[Join[input, <|"BondClasses" -> bondClasses|>]]
];

End[]

EndPackage[]
