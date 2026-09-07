(* ::Package:: *)

BeginPackage["InducedRepresentation`"]

CompileInducedOrbitRepresentation::usage =
  "CompileInducedOrbitRepresentation[group,subgroup,matrices,reps,flags] constructs one induced orbit from integer group data only.";
CompileInducedRepresentation::usage =
  "CompileInducedRepresentation[orbitSpecifications,group,flags] assembles several induced orbits without geometric coordinates.";

Begin["`Private`"]

validateLocalMatrices[group_, subgroup_List, matrices_List] := And[
  GroupAlgebra`FindGeneratorIndices[group, subgroup] =!= $Failed,
  Length[matrices] === Length[subgroup],
  matrices =!= {},
  And @@ (MatrixPredicates`ExactSquareMatrixQ /@ matrices),
  SameQ @@ (Dimensions /@ matrices),
  And @@ (MatrixPredicates`ExactUnitaryMatrixQ /@ matrices)
];

CompileInducedOrbitRepresentation::data =
  "The subgroup, local matrices, coset representatives, or Schreier decomposition is invalid.";

CompileInducedOrbitRepresentation[
    group_Association,
    subgroupIndices_List,
    siteSymmetryMatrices_List,
    cosetRepresentatives_List,
    antiunitaryFlags_List
  ] := Module[
  {subgroup, positions, operationCount, siteCount, records, images, blocks},
  subgroup = DeleteDuplicates[subgroupIndices];
  operationCount = Lookup[group, "Order", 0];
  siteCount = Length[cosetRepresentatives];
  If[!GroupAlgebra`GroupAlgebraQ[group] ||
      Length[antiunitaryFlags] =!= operationCount ||
      !And @@ (BooleanQ /@ antiunitaryFlags) || siteCount < 1 ||
      !And @@ (IntegerQ[#] && Between[#, {1, operationCount}] & /@ cosetRepresentatives) ||
      !validateLocalMatrices[group, subgroup, siteSymmetryMatrices],
    Message[CompileInducedOrbitRepresentation::data]; Return[$Failed]
  ];
  positions = AssociationThread[subgroup, Range[Length[subgroup]]];
  records = Table[
    GroupAlgebra`SchreierDecompositionIndex[
      group, subgroup, cosetRepresentatives, operation, source],
    {operation, operationCount}, {source, siteCount}
  ];
  If[MemberQ[records, $Failed, Infinity],
    Message[CompileInducedOrbitRepresentation::data]; Return[$Failed]
  ];
  images = records[[All, All, 1]];
  If[!And @@ (Sort[#] === Range[siteCount] & /@ images),
    Message[CompileInducedOrbitRepresentation::data]; Return[$Failed]
  ];
  blocks = Map[
    Function[record,
      With[{matrix = siteSymmetryMatrices[[positions[record[[2]]]]]},
        If[
          antiunitaryFlags[[cosetRepresentatives[[record[[1]]]]]],
          Conjugate[matrix],
          matrix
        ]
      ]
    ],
    records,
    {2}
  ];
  <|
    "SiteSymmetryOperationIndices" -> subgroup,
    "CosetRepresentativeIndices" -> cosetRepresentatives,
    "SchreierRecords" -> records,
    "ImageSiteIndices" -> images,
    "LocalBlocks" -> blocks,
    "LocalDimension" -> Length[First[siteSymmetryMatrices]],
    "Verified" -> True
  |>
];
CompileInducedOrbitRepresentation[___] :=
  (Message[CompileInducedOrbitRepresentation::data]; $Failed);

CompileInducedRepresentation::data =
  "Each orbit specification must contain subgroup indices, local matrices, and coset representatives.";
CompileInducedRepresentation[
    orbitSpecifications_List,
    group_Association,
    antiunitaryFlags_List
  ] := Module[{orbitResults, assembled},
  If[orbitSpecifications === {} || !And @@ (AssociationQ /@ orbitSpecifications),
    Message[CompileInducedRepresentation::data]; Return[$Failed]
  ];
  orbitResults = Map[
    CompileInducedOrbitRepresentation[
      group,
      Lookup[#, "SiteSymmetryOperationIndices", $Failed],
      Lookup[#, "SiteSymmetryMatrices", $Failed],
      Lookup[#, "CosetRepresentativeIndices", $Failed],
      antiunitaryFlags
    ] &,
    orbitSpecifications
  ];
  If[MemberQ[orbitResults, $Failed], Return[$Failed]];
  assembled = RepresentationData`AssembleBlockMonomialRepresentation[
    Lookup[orbitResults, "ImageSiteIndices"],
    Lookup[orbitResults, "LocalBlocks"],
    Lookup[orbitResults, "LocalDimension"],
    antiunitaryFlags
  ];
  If[assembled === $Failed, Return[$Failed]];
  Join[
    <|
      "Method" -> "Induced",
      "Convention" -> "orbit-major induced site action",
      "GroupAlgebra" -> group,
      "AntiunitaryFlags" -> antiunitaryFlags,
      "SiteSymmetryData" -> orbitResults,
      "LocalBlocks" -> Lookup[orbitResults, "LocalBlocks"],
      "LocalDimensions" -> Lookup[orbitResults, "LocalDimension"]
    |>,
    assembled
  ]
];
CompileInducedRepresentation[___] :=
  (Message[CompileInducedRepresentation::data]; $Failed);

End[]

EndPackage[]
