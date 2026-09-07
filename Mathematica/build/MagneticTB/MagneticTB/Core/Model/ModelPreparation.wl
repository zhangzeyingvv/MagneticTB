(* ::Package:: *)

BeginPackage["ModelPreparation`"]

CompileModelPointRepresentations::usage =
  "CompileModelPointRepresentations[model] computes one exact local representation matrix per site orbit and point operation from a clean ModelSpecification.";
PrepareDirectProductModel::usage =
  "PrepareDirectProductModel[model] compiles local representations and direct-product action data once and returns a clean PreparedModel. PrepareDirectProductModel[model,localRepresentations] accepts already compiled exact local matrices from another representation frontend.";
AssemblePreparedRepresentation::usage =
  "AssemblePreparedRepresentation[prepared] assembles the full exact unitary matrices from a PreparedModel.";

Begin["`Private`"]

CompileModelPointRepresentations::model =
  "Expected a validated ModelSpecification.";
CompileModelPointRepresentations[model_Association] := Module[
  {localBases, records, variables, representations},
  If[!ModelSchema`ModelSpecificationQ[model],
    Message[CompileModelPointRepresentations::model];
    Return[$Failed]
  ];
  localBases = model["LocalBases"];
  records = model["PointOperationRecords"];
  variables = model["Variables"];
  representations = MapThread[
    Function[{basis, operationRecords},
      FunctionBasisRepresentation`PointOperationRepresentationMatrices[
        basis,
        operationRecords,
        variables
      ]
    ],
    {localBases, records}
  ];
  If[MemberQ[representations, $Failed], $Failed, representations]
];
CompileModelPointRepresentations[_] :=
  (Message[CompileModelPointRepresentations::model]; $Failed);

PrepareDirectProductModel::model =
  "Expected a validated ModelSpecification compatible with the direct-product representation compiler.";
PrepareDirectProductModel[model_Association] := Module[
  {permutationData},
  If[!ModelSchema`ModelSpecificationQ[model],
    Message[PrepareDirectProductModel::model];
    Return[$Failed]
  ];
  permutationData = SitePermutationCompiler`CompileSitePermutationData[
    model["SiteOrbits"], model["Symmetry", "SpatialActions"]];
  If[permutationData === $Failed, Return[$Failed]];
  PrepareDirectProductModel[model, permutationData]
];

PrepareDirectProductModel[
    model_Association, permutationData_Association
  ] := Module[{localRepresentations},
  If[!ModelSchema`ModelSpecificationQ[model] ||
      !SitePermutationCompiler`SitePermutationDataQ[permutationData],
    Message[PrepareDirectProductModel::model]; Return[$Failed]
  ];
  localRepresentations = CompileModelPointRepresentations[model];
  If[localRepresentations === $Failed, Return[$Failed]];
  PrepareDirectProductModel[model, localRepresentations, permutationData]
];

PrepareDirectProductModel[
    model_Association,
    localRepresentations_List
  ] := Module[{actionCompilation},
  If[!ModelSchema`ModelSpecificationQ[model],
    Message[PrepareDirectProductModel::model];
    Return[$Failed]
  ];
  actionCompilation = Module[{permutationData},
    permutationData = SitePermutationCompiler`CompileSitePermutationData[
      model["SiteOrbits"], model["Symmetry", "SpatialActions"]];
    If[permutationData === $Failed, Return[$Failed, Module]];
    PhysicalRepresentation`CompileDirectProductActionData[
      model["SiteOrbits"], localRepresentations,
      model["Symmetry", "SpatialActions"], model["Symmetry", "GroupAlgebra"],
      model["Symmetry", "AntiunitaryFlags"], permutationData]
  ];
  If[actionCompilation === $Failed, Return[$Failed]];
  ModelSchema`CreatePreparedModel[model, actionCompilation]
];
PrepareDirectProductModel[
    model_Association, localRepresentations_List,
    permutationData_Association
  ] := Module[{actionCompilation},
  If[!ModelSchema`ModelSpecificationQ[model] ||
      !SitePermutationCompiler`SitePermutationDataQ[permutationData] ||
      permutationData["SiteOrbits"] =!= model["SiteOrbits"] ||
      permutationData["SpatialActions"] =!= model["Symmetry", "SpatialActions"],
    Message[PrepareDirectProductModel::model]; Return[$Failed]
  ];
  actionCompilation = PhysicalRepresentation`CompileDirectProductActionData[
    model["SiteOrbits"], localRepresentations,
    model["Symmetry", "SpatialActions"], model["Symmetry", "GroupAlgebra"],
    model["Symmetry", "AntiunitaryFlags"], permutationData];
  If[actionCompilation === $Failed, Return[$Failed]];
  ModelSchema`CreatePreparedModel[model, actionCompilation]
];
PrepareDirectProductModel[_] :=
  (Message[PrepareDirectProductModel::model]; $Failed);
PrepareDirectProductModel[_, _] :=
  (Message[PrepareDirectProductModel::model]; $Failed);
PrepareDirectProductModel[_, _, _] :=
  (Message[PrepareDirectProductModel::model]; $Failed);

AssemblePreparedRepresentation::model =
  "Expected a validated PreparedModel.";
AssemblePreparedRepresentation[
    prepared_Association
  ] := Module[{compilation, result},
  If[!ModelSchema`PreparedModelQ[prepared],
    Message[AssemblePreparedRepresentation::model];
    Return[$Failed]
  ];
  compilation = <|
    "Method" -> prepared["RepresentationMethod"],
    "SpatialActions" -> prepared["Symmetry", "SpatialActions"],
    "SpinActions" -> prepared["Symmetry", "SpinActions"],
    "GroupAlgebra" -> prepared["Symmetry", "GroupAlgebra"],
    "AntiunitaryFlags" -> prepared["Symmetry", "AntiunitaryFlags"],
    "RepresentationMatrices" -> prepared["RepresentationMatrices"],
    "OrbitRepresentationMatrices" ->
      prepared["OrbitRepresentationMatrices"],
    "LocalRepresentationMatrices" ->
      prepared["LocalRepresentationMatrices"],
    "LocalDimensions" -> prepared["LocalDimensions"]
  |>;
  result = RepresentationData`CreateRepresentationData[
    <|
      "RepresentationMatrices" -> compilation["RepresentationMatrices"],
      "AntiunitaryFlags" -> compilation["AntiunitaryFlags"],
      "GroupAlgebra" -> compilation["GroupAlgebra"],
      "Convention" -> Lookup[
        compilation,
        "Convention",
        "columns are source states"
      ]
    |>
  ];
  If[result === $Failed, Return[$Failed]];
  Join[
    compilation,
    KeyDrop[result, {"Verified"}],
    <|"PreparedModel" -> prepared|>
  ]
];
AssemblePreparedRepresentation[___] :=
  (Message[AssemblePreparedRepresentation::model]; $Failed);

End[]

EndPackage[]
