(* ::Package:: *)

BeginPackage["PhysicalRepresentation`"]

CompileDirectProductActionData::usage =
  "CompileDirectProductActionData[sites,localMatrices,spatial,group,flags,permutations] combines verified site geometry and direct-product local matrices.";
CompileInducedActionData::usage =
  "CompileInducedActionData[sites,localData,spatial,group,flags,permutations] combines verified site geometry and induced local data.";
AutomaticCosetRepresentativeIndices::usage =
  "AutomaticCosetRepresentativeIndices[group,images,flags,reference] chooses one operation carrying the reference site to each target site, preferring a unitary operation when available.";

Begin["`Private`"]

validPermutationDataQ[data_, siteOrbits_, spatialActions_] :=
  SitePermutationCompiler`SitePermutationDataQ[data] &&
    Lookup[data, "SiteOrbits", $Failed] === siteOrbits &&
    Lookup[data, "SpatialActions", $Failed] === spatialActions;

CompileDirectProductActionData::data =
  "The site permutation data, local matrices, group, or antiunitary flags are incompatible.";

CompileDirectProductActionData[
    siteOrbits_List, localRepresentations_List, spatialActions_List,
    group_Association, antiunitaryFlags_List, permutationData_Association
  ] := Module[{mathematical, actionData},
  If[!validPermutationDataQ[permutationData, siteOrbits, spatialActions],
    Message[CompileDirectProductActionData::data]; Return[$Failed]
  ];
  mathematical = DirectProductRepresentation`CompileDirectProductRepresentation[
    permutationData["ImageSiteIndices"], localRepresentations,
    group, antiunitaryFlags];
  If[mathematical === $Failed,
    Message[CompileDirectProductActionData::data]; Return[$Failed]
  ];
  actionData = RepresentationActionData`CreateRepresentationActionData[
    <|
      "Method" -> "DirectProduct",
      "SiteOrbits" -> siteOrbits,
      "SpatialActions" -> spatialActions,
      "GroupAlgebra" -> group,
      "AntiunitaryFlags" -> antiunitaryFlags,
      "ImageSiteIndices" -> permutationData["ImageSiteIndices"],
      "CellTranslations" -> permutationData["CellTranslations"],
      "LocalBlocks" -> mathematical["LocalBlocks"],
      "LocalDimensions" -> mathematical["LocalDimensions"]
    |>];
  If[actionData === $Failed, Return[$Failed]];
  Join[
    KeyTake[mathematical, {
      "Method", "Convention", "GroupAlgebra", "AntiunitaryFlags",
      "LocalRepresentationMatrices", "LocalBlocks", "LocalDimensions",
      "OrbitRepresentationMatrices", "RepresentationMatrices",
      "BlockDimensions", "Dimension", "UnitaryVerified"}],
    <|"ActionData" -> actionData, "SitePermutationData" -> permutationData|>
  ]
];
CompileDirectProductActionData[___] :=
  (Message[CompileDirectProductActionData::data]; $Failed);

AutomaticCosetRepresentativeIndices::data =
  "The group, site action, antiunitary flags, or reference-site index is invalid.";
AutomaticCosetRepresentativeIndices[
    group_Association, images_List, flags_List, reference_Integer
  ] := Module[{representatives, matches},
  If[!GroupAlgebra`GroupActionTableQ[group, images] ||
      Length[flags] =!= group["Order"] ||
      !And @@ (BooleanQ /@ flags) ||
      !Between[reference, {1, Length[First[images]]}],
    Message[AutomaticCosetRepresentativeIndices::data];
    Return[$Failed]
  ];
  representatives = Table[
    matches = GroupAlgebra`TransporterIndices[
      group, images, reference, target
    ];
    If[
      matches === $Failed || matches === {},
      $Failed,
      SelectFirst[matches, !flags[[#]] &, First[matches]]
    ],
    {target, Length[First[images]]}
  ];
  If[MemberQ[representatives, $Failed], $Failed, representatives]
];
AutomaticCosetRepresentativeIndices[___] :=
  (Message[AutomaticCosetRepresentativeIndices::data]; $Failed);

CompileInducedActionData::data =
  "The site permutation data, local site-symmetry data, group, or flags are incompatible.";
CompileInducedActionData::orbitdata =
  "Orbit `1` has invalid subgroup, local-matrix, coset-representative, or explicit function-basis action data.";

CompileInducedActionData[
    siteOrbits_List, localData_List, spatialActions_List,
    group_Association, antiunitaryFlags_List, permutationData_Association
  ] := Module[
  {orbitSpecifications, orbitResults, reference, expectedSubgroup, suppliedSubgroup,
   matrices, representatives, mathematical, actionData, failure,
   explicitFunctionBasisQ, localBlocks, localDimensions, assembled},
  If[!validPermutationDataQ[permutationData, siteOrbits, spatialActions] ||
      !GroupAlgebra`GroupAlgebraQ[group] || group["Order"] =!= Length[spatialActions] ||
      Length[antiunitaryFlags] =!= group["Order"] ||
      !And @@ (BooleanQ /@ antiunitaryFlags) ||
      Length[localData] =!= Length[siteOrbits] ||
      !And @@ (AssociationQ /@ localData),
    Message[CompileInducedActionData::data]; Return[$Failed]
  ];
  orbitSpecifications = Table[
    reference = Lookup[localData[[orbit]], "ReferenceSiteIndex", 1];
    If[!IntegerQ[reference] || !Between[reference, {1, Length[siteOrbits[[orbit]]]}],
      Message[CompileInducedActionData::data]; Return[$Failed]
    ];
    expectedSubgroup = SitePermutationCompiler`SiteSymmetryOperationIndices[
      permutationData, orbit, reference];
    suppliedSubgroup = Lookup[localData[[orbit]],
      "SiteSymmetryOperationIndices", expectedSubgroup];
    If[expectedSubgroup === $Failed ||
        Sort[suppliedSubgroup] =!= Sort[expectedSubgroup],
      Message[CompileInducedActionData::data]; Return[$Failed]
    ];
    matrices = Lookup[localData[[orbit]], "SiteSymmetryMatrices", $Failed];
    representatives = Lookup[localData[[orbit]],
      "CosetRepresentativeIndices", Automatic];
    If[representatives === Automatic,
      representatives = AutomaticCosetRepresentativeIndices[
        group,
        permutationData["ImageSiteIndices"][[orbit]],
        antiunitaryFlags, reference]
    ];
    If[representatives === $Failed,
      Message[CompileInducedActionData::data];
      Return[$Failed]
    ];
    Join[
      <|
        "ReferenceSiteIndex" -> reference,
        "SiteSymmetryOperationIndices" -> suppliedSubgroup,
        "SiteSymmetryMatrices" -> matrices,
        "CosetRepresentativeIndices" -> representatives
      |>,
      KeyTake[localData[[orbit]], {"TransportedBases", "LocalBlocks"}]
    ],
    {orbit, Length[siteOrbits]}
  ];
  explicitFunctionBasisQ = And @@ (
    KeyExistsQ[#, "TransportedBases"] && KeyExistsQ[#, "LocalBlocks"] & /@
      orbitSpecifications
  );
  If[explicitFunctionBasisQ,
    localBlocks = Lookup[orbitSpecifications, "LocalBlocks"];
    localDimensions = Length[First[#]] & /@
      Lookup[orbitSpecifications, "SiteSymmetryMatrices"];
    assembled = RepresentationData`AssembleBlockMonomialRepresentation[
      permutationData["ImageSiteIndices"], localBlocks,
      localDimensions, antiunitaryFlags
    ];
    If[assembled === $Failed,
      Message[CompileInducedActionData::data]; Return[$Failed]
    ];
    mathematical = Join[
      <|
        "Method" -> "Induced",
        "Convention" -> "explicit function-basis orbit action",
        "GroupAlgebra" -> group,
        "AntiunitaryFlags" -> antiunitaryFlags,
        "LocalBlocks" -> localBlocks,
        "LocalDimensions" -> localDimensions
      |>,
      assembled
    ];
    orbitResults = Map[
      Join[#, <|"LocalDimension" -> Length[First[#["SiteSymmetryMatrices"]]],
        "Verified" -> True|>] &,
      orbitSpecifications
    ],
    mathematical = InducedRepresentation`CompileInducedRepresentation[
      orbitSpecifications, group, antiunitaryFlags
    ];
    If[mathematical === $Failed,
      failure = SelectFirst[
        Range[Length[orbitSpecifications]],
        InducedRepresentation`CompileInducedOrbitRepresentation[
          group,
          orbitSpecifications[[#, "SiteSymmetryOperationIndices"]],
          orbitSpecifications[[#, "SiteSymmetryMatrices"]],
          orbitSpecifications[[#, "CosetRepresentativeIndices"]],
          antiunitaryFlags
        ] === $Failed &,
        Missing["UndiagnosedOrbit"]
      ];
      If[MissingQ[failure],
        Message[CompileInducedActionData::data],
        Message[CompileInducedActionData::orbitdata, failure]
      ];
      Return[$Failed]
    ];
    If[Lookup[mathematical["SiteSymmetryData"], "ImageSiteIndices"] =!=
        permutationData["ImageSiteIndices"],
      Message[CompileInducedActionData::data]; Return[$Failed]
    ];
    orbitResults = MapThread[
      Join[#1, KeyTake[#2, {
        "SiteSymmetryOperationIndices", "CosetRepresentativeIndices",
        "SchreierRecords", "LocalBlocks", "LocalDimension", "Verified"}]] &,
      {orbitSpecifications, mathematical["SiteSymmetryData"]}
    ]
  ];
  actionData = RepresentationActionData`CreateRepresentationActionData[
    <|
      "Method" -> "Induced",
      "SiteOrbits" -> siteOrbits,
      "SpatialActions" -> spatialActions,
      "GroupAlgebra" -> group,
      "AntiunitaryFlags" -> antiunitaryFlags,
      "ImageSiteIndices" -> permutationData["ImageSiteIndices"],
      "CellTranslations" -> permutationData["CellTranslations"],
      "LocalBlocks" -> mathematical["LocalBlocks"],
      "LocalDimensions" -> mathematical["LocalDimensions"]
    |>];
  If[actionData === $Failed, Return[$Failed]];
  Join[
    KeyTake[mathematical, {
      "Method", "Convention", "GroupAlgebra", "AntiunitaryFlags",
      "LocalBlocks", "LocalDimensions", "OrbitRepresentationMatrices",
      "RepresentationMatrices", "BlockDimensions", "Dimension", "UnitaryVerified"}],
    <|
      "ActionData" -> actionData,
      "SitePermutationData" -> permutationData,
      "SiteSymmetryData" -> orbitResults
    |>
  ]
];
CompileInducedActionData[___] :=
  (Message[CompileInducedActionData::data]; $Failed);

End[]

EndPackage[]
