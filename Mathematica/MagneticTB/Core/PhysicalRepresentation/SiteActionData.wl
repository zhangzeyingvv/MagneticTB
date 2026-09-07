(* ::Package:: *)

BeginPackage["RepresentationActionData`"]

CreateRepresentationActionData::usage =
  "CreateRepresentationActionData[spec] validates the physical site/local-block action consumed by bond compilation.";
SiteActionDataQ::usage =
  "SiteActionDataQ[data] deeply validates the canonical physical site-action schema.";
CompileSiteActionAccessor::usage =
  "CompileSiteActionAccessor[data] deeply validates SiteActionData once and returns an opaque accessor for repeated internal site-action lookup.";
RepresentationSiteAction::usage =
  "RepresentationSiteAction[data,orbit,source,operation] returns image site, cell translation, local matrix, and antiunitary flag. It also accepts the opaque result of CompileSiteActionAccessor for repeated internal lookup.";
AssembleFullRepresentation::usage =
  "AssembleFullRepresentation[data] assembles full matrices through the pure RepresentationData module.";

Begin["`Private`"]

ClearAll[
  siteActionValidationCache,
  validSiteActionDataAssociationQ,
  verifiedSiteActionAccessor,
  compiledSiteActionAccessorQ
];
siteActionValidationCache = <||>;

validSiteActionDataAssociationQ[data_] := Module[
  {siteOrbits, spatialActions, group, images, shifts, blocks, dimensions,
   flags, orbitCount, operationCount, coordinateDimension},
  If[!AssociationQ[data] || Lookup[data, "Schema", None] =!= "SiteActionData" ||
      !TrueQ[Lookup[data, "Verified", False]], Return[False]];
  siteOrbits = Lookup[data, "SiteOrbits", $Failed];
  spatialActions = Lookup[data, "SpatialActions", $Failed];
  group = Lookup[data, "GroupAlgebra", $Failed];
  images = Lookup[data, "ImageSiteIndices", $Failed];
  shifts = Lookup[data, "CellTranslations", $Failed];
  blocks = Lookup[data, "LocalBlocks", $Failed];
  dimensions = Lookup[data, "LocalDimensions", $Failed];
  flags = Lookup[data, "AntiunitaryFlags", $Failed];
  If[!SitePermutationCompiler`ValidSiteOrbitsQ[siteOrbits] ||
      !ListQ[spatialActions] || spatialActions === {} ||
      !GroupAlgebra`GroupAlgebraQ[group] || !ListQ[images] ||
      !ListQ[shifts] || !ListQ[blocks] || !ListQ[dimensions] ||
      !ListQ[flags] || !And @@ (BooleanQ /@ flags), Return[False]];
  orbitCount = Length[siteOrbits];
  operationCount = Length[spatialActions];
  coordinateDimension = Length[siteOrbits[[1, 1]]];
  If[group["Order"] =!= operationCount || Length[flags] =!= operationCount ||
      Length[images] =!= orbitCount || Length[shifts] =!= orbitCount ||
      Length[blocks] =!= orbitCount || Length[dimensions] =!= orbitCount ||
      !And @@ (IntegerQ[#] && Positive[#] & /@ dimensions) ||
      !And @@ (SymmetryAlgebra`ValidSpatialOperationQ[#, coordinateDimension] & /@
        spatialActions) ||
      !And @@ (GroupAlgebra`GroupActionTableQ[group, #] & /@ images),
    Return[False]
  ];
  And @@ Flatten@Table[
    Length[images[[orbit]]] === operationCount &&
      Length[shifts[[orbit]]] === operationCount &&
      Length[blocks[[orbit]]] === operationCount &&
      Length[images[[orbit, operation]]] === Length[siteOrbits[[orbit]]] &&
      Sort[images[[orbit, operation]]] === Range[Length[siteOrbits[[orbit]]]] &&
      Length[shifts[[orbit, operation]]] === Length[siteOrbits[[orbit]]] &&
      And @@ (VectorQ[#] && Length[#] === coordinateDimension &&
          MatrixPredicates`IntegerVectorQ[#] & /@ shifts[[orbit, operation]]) &&
      Length[blocks[[orbit, operation]]] === Length[siteOrbits[[orbit]]] &&
      And @@ (Dimensions[#] === {dimensions[[orbit]], dimensions[[orbit]]} &&
          MatrixPredicates`ExactUnitaryMatrixQ[#] & /@ blocks[[orbit, operation]]),
    {orbit, orbitCount}, {operation, operationCount}
  ]
];

SiteActionDataQ[data_] := Module[{cached, result},
  cached = Lookup[siteActionValidationCache, Hash[data], Missing["NotCached"]];
  If[!MissingQ[cached], Return[cached]];
  result = TrueQ[validSiteActionDataAssociationQ[data]];
  AssociateTo[siteActionValidationCache, Hash[data] -> result];
  result
];

CreateRepresentationActionData::spec =
  "The site-action specification has incompatible or unverified geometry, local blocks, flags, or dimensions.";

CreateRepresentationActionData[spec_Association] := Module[
  {method, siteOrbits, spatialActions, group, images, shifts, blocks,
   dimensions, flags, orbitCount, operationCount, coordinateDimension,
   candidate},
  method = Lookup[spec, "Method", "Unspecified"];
  siteOrbits = Lookup[spec, "SiteOrbits", $Failed];
  spatialActions = Lookup[spec, "SpatialActions", $Failed];
  group = Lookup[spec, "GroupAlgebra", $Failed];
  images = Lookup[spec, "ImageSiteIndices", $Failed];
  shifts = Lookup[spec, "CellTranslations", $Failed];
  blocks = Lookup[spec, "LocalBlocks", $Failed];
  dimensions = Lookup[spec, "LocalDimensions", $Failed];
  flags = Lookup[spec, "AntiunitaryFlags", $Failed];
  If[!StringQ[method] || !SitePermutationCompiler`ValidSiteOrbitsQ[siteOrbits] ||
      !ListQ[spatialActions] || spatialActions === {} ||
      !GroupAlgebra`GroupAlgebraQ[group] || !ListQ[images] ||
      !ListQ[shifts] || !ListQ[blocks] || !ListQ[dimensions] ||
      !ListQ[flags] || !And @@ (BooleanQ /@ flags),
    Message[CreateRepresentationActionData::spec]; Return[$Failed]
  ];
  orbitCount = Length[siteOrbits]; operationCount = Length[spatialActions];
  coordinateDimension = Length[siteOrbits[[1, 1]]];
  If[group["Order"] =!= operationCount || Length[flags] =!= operationCount ||
      Length[images] =!= orbitCount || Length[shifts] =!= orbitCount ||
      Length[blocks] =!= orbitCount || Length[dimensions] =!= orbitCount ||
      !And @@ (IntegerQ[#] && Positive[#] & /@ dimensions) ||
      !And @@ (GroupAlgebra`GroupActionTableQ[group, #] & /@ images) ||
      !And @@ (SymmetryAlgebra`ValidSpatialOperationQ[#, coordinateDimension] & /@
        spatialActions),
    Message[CreateRepresentationActionData::spec]; Return[$Failed]
  ];
  Do[
    If[Length[images[[orbit]]] =!= operationCount ||
        Length[shifts[[orbit]]] =!= operationCount ||
        Length[blocks[[orbit]]] =!= operationCount,
      Message[CreateRepresentationActionData::spec]; Return[$Failed]
    ];
    Do[
      If[Length[images[[orbit, operation]]] =!= Length[siteOrbits[[orbit]]] ||
          Sort[images[[orbit, operation]]] =!= Range[Length[siteOrbits[[orbit]]]] ||
          Length[shifts[[orbit, operation]]] =!= Length[siteOrbits[[orbit]]] ||
          !And @@ (VectorQ[#] && Length[#] === coordinateDimension &&
              MatrixPredicates`IntegerVectorQ[#] & /@ shifts[[orbit, operation]]) ||
          Length[blocks[[orbit, operation]]] =!= Length[siteOrbits[[orbit]]] ||
          !And @@ (Dimensions[#] === {dimensions[[orbit]], dimensions[[orbit]]} &&
              MatrixPredicates`ExactUnitaryMatrixQ[#] & /@ blocks[[orbit, operation]]),
        Message[CreateRepresentationActionData::spec]; Return[$Failed]
      ],
      {operation, operationCount}
    ],
    {orbit, orbitCount}
  ];
  candidate = <|
    "Schema" -> "SiteActionData",
    "Method" -> method,
    "Convention" -> "orbit-local indices; columns are source local states",
    "SiteOrbits" -> siteOrbits,
    "SpatialActions" -> spatialActions,
    "GroupAlgebra" -> group,
    "AntiunitaryFlags" -> flags,
    "ImageSiteIndices" -> images,
    "CellTranslations" -> shifts,
    "LocalBlocks" -> blocks,
    "LocalDimensions" -> dimensions,
    "Verified" -> True
  |>;
  AssociateTo[siteActionValidationCache, Hash[candidate] -> True];
  candidate
];
CreateRepresentationActionData[_] :=
  (Message[CreateRepresentationActionData::spec]; $Failed);

CompileSiteActionAccessor::data = "Expected verified SiteActionData.";
CompileSiteActionAccessor[data_Association] := If[
  SiteActionDataQ[data],
  verifiedSiteActionAccessor[
    data["SiteOrbits"],
    data["ImageSiteIndices"],
    data["CellTranslations"],
    data["LocalBlocks"],
    data["AntiunitaryFlags"]
  ],
  Message[CompileSiteActionAccessor::data];
  $Failed
];
CompileSiteActionAccessor[_] :=
  (Message[CompileSiteActionAccessor::data]; $Failed);

compiledSiteActionAccessorQ[
    verifiedSiteActionAccessor[
      siteOrbits_List,
      images_List,
      shifts_List,
      blocks_List,
      flags_List
    ]
  ] := siteOrbits =!= {} &&
  Length[images] === Length[siteOrbits] &&
  Length[shifts] === Length[siteOrbits] &&
  Length[blocks] === Length[siteOrbits] &&
  flags =!= {};
compiledSiteActionAccessorQ[_] := False;

RepresentationSiteAction::data =
  "Expected verified SiteActionData and valid orbit/source/operation indices.";
RepresentationSiteAction[
    accessor_verifiedSiteActionAccessor,
    orbit_Integer?Positive,
    source_Integer?Positive, operation_Integer?Positive
  ] := Module[
  {siteOrbits, images, shifts, blocks, flags},
  {siteOrbits, images, shifts, blocks, flags} = List @@ accessor;
  If[!compiledSiteActionAccessorQ[accessor] ||
      orbit > Length[siteOrbits] || source > Length[siteOrbits[[orbit]]] ||
      orbit > Length[images] || orbit > Length[shifts] || orbit > Length[blocks] ||
      operation > Length[flags] || operation > Length[images[[orbit]]] ||
      operation > Length[shifts[[orbit]]] || operation > Length[blocks[[orbit]]] ||
      source > Length[images[[orbit, operation]]] ||
      source > Length[shifts[[orbit, operation]]] ||
      source > Length[blocks[[orbit, operation]]],
    Message[RepresentationSiteAction::data]; Return[$Failed]
  ];
  <|
    "OrbitIndex" -> orbit,
    "SourceSiteIndex" -> source,
    "ImageSiteIndex" -> images[[orbit, operation, source]],
    "CellTranslation" -> shifts[[orbit, operation, source]],
    "Matrix" -> blocks[[orbit, operation, source]],
    "Antiunitary" -> flags[[operation]]
  |>
];
RepresentationSiteAction[
    data_Association,
    orbit_Integer?Positive,
    source_Integer?Positive,
    operation_Integer?Positive
  ] := Module[{accessor},
  accessor = Quiet@CompileSiteActionAccessor[data];
  If[accessor === $Failed,
    Message[RepresentationSiteAction::data];
    $Failed,
    RepresentationSiteAction[accessor, orbit, source, operation]
  ]
];
RepresentationSiteAction[___] :=
  (Message[RepresentationSiteAction::data]; $Failed);

AssembleFullRepresentation::data = "Expected verified SiteActionData.";
AssembleFullRepresentation[data_Association] := Module[{assembled},
  If[!SiteActionDataQ[data],
    Message[AssembleFullRepresentation::data]; Return[$Failed]
  ];
  assembled = RepresentationData`AssembleBlockMonomialRepresentation[
    data["ImageSiteIndices"], data["LocalBlocks"],
    data["LocalDimensions"], data["AntiunitaryFlags"]];
  If[assembled === $Failed, Return[$Failed]];
  Join[
    <|
      "Method" -> data["Method"],
      "Convention" -> "orbit-major, site-major, local-index-minor",
      "GroupAlgebra" -> data["GroupAlgebra"],
      "AntiunitaryFlags" -> data["AntiunitaryFlags"]
    |>,
    assembled
  ]
];
AssembleFullRepresentation[_] :=
  (Message[AssembleFullRepresentation::data]; $Failed);

End[]

EndPackage[]
