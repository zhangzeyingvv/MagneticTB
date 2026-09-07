(* ::Package:: *)

BeginPackage["SitePermutationCompiler`"]

ValidSiteOrbitsQ::usage = "ValidSiteOrbitsQ[siteOrbits] validates finite coordinate orbits.";
CompileSitePermutationData::usage =
  "CompileSitePermutationData[siteOrbits,spatialActions] returns operation-by-source site images and integer cell translations for every orbit.";
SitePermutationDataQ::usage =
  "SitePermutationDataQ[data] deeply validates operation-by-source site images and integer cell translations.";
SiteSymmetryOperationIndices::usage =
  "SiteSymmetryOperationIndices[sitePermutationData,orbit,site] returns operations fixing one site.";

Begin["`Private`"]

ClearAll[
  sitePermutationValidationCache,
  exactSiteImageEquationQ,
  validSitePermutationDataAssociationQ
];
sitePermutationValidationCache = <||>;

ValidSiteOrbitsQ[siteOrbits_List] := Module[{dimensions},
  If[siteOrbits === {} ||
      !And @@ (ListQ[#] && # =!= {} & /@ siteOrbits) ||
      !And @@ Map[VectorQ, siteOrbits, {2}], Return[False]];
  dimensions = Length /@ Flatten[siteOrbits, 1];
  dimensions =!= {} && First[dimensions] > 0 && SameQ @@ dimensions
];
ValidSiteOrbitsQ[_] := False;

exactSiteImageEquationQ[
    siteOrbit_List,
    operation_Association,
    source_Integer,
    target_Integer,
    cellShift_List
  ] := Module[{difference},
  If[!Between[source, {1, Length[siteOrbit]}] ||
      !Between[target, {1, Length[siteOrbit]}], Return[False]];
  difference = FullSimplify[
    operation["Rotation"].siteOrbit[[source]] + operation["Translation"] -
      siteOrbit[[target]] - cellShift
  ];
  VectorQ[difference] &&
    And @@ (MatrixPredicates`ExactZeroExpressionQ /@ difference)
];
exactSiteImageEquationQ[___] := False;

validSitePermutationDataAssociationQ[data_] := Module[
  {siteOrbits, spatialActions, records, images, shifts, orbitCount,
   operationCount, dimension},
  If[!AssociationQ[data] || Lookup[data, "Schema", None] =!= "SitePermutationData" ||
      Lookup[data, "SchemaVersion", None] =!= 1 ||
      !TrueQ[Lookup[data, "Verified", False]], Return[False]];
  siteOrbits = Lookup[data, "SiteOrbits", $Failed];
  spatialActions = Lookup[data, "SpatialActions", $Failed];
  records = Lookup[data, "ImageRecords", $Failed];
  images = Lookup[data, "ImageSiteIndices", $Failed];
  shifts = Lookup[data, "CellTranslations", $Failed];
  If[!ValidSiteOrbitsQ[siteOrbits] || !ListQ[spatialActions] ||
      spatialActions === {} || !ListQ[records] || !ListQ[images] ||
      !ListQ[shifts], Return[False]];
  orbitCount = Length[siteOrbits]; operationCount = Length[spatialActions];
  dimension = Length[siteOrbits[[1, 1]]];
  If[Length[records] =!= orbitCount || Length[images] =!= orbitCount ||
      Length[shifts] =!= orbitCount ||
      !And @@ (SymmetryAlgebra`ValidSpatialOperationQ[#, dimension] & /@
        spatialActions), Return[False]];
  And @@ Flatten@Table[
    Length[records[[orbit]]] === operationCount &&
      Length[images[[orbit]]] === operationCount &&
      Length[shifts[[orbit]]] === operationCount &&
      Length[records[[orbit, operation]]] === Length[siteOrbits[[orbit]]] &&
      Length[images[[orbit, operation]]] === Length[siteOrbits[[orbit]]] &&
      Sort[images[[orbit, operation]]] === Range[Length[siteOrbits[[orbit]]]] &&
      Length[shifts[[orbit, operation]]] === Length[siteOrbits[[orbit]]] &&
      And @@ Table[
        AssociationQ[records[[orbit, operation, source]]] &&
          records[[orbit, operation, source, "SourceSite"]] === source &&
          records[[orbit, operation, source, "TargetSite"]] ===
            images[[orbit, operation, source]] &&
          records[[orbit, operation, source, "CellTranslation"]] ===
            shifts[[orbit, operation, source]] &&
          VectorQ[shifts[[orbit, operation, source]]] &&
          Length[shifts[[orbit, operation, source]]] === dimension &&
          MatrixPredicates`IntegerVectorQ[shifts[[orbit, operation, source]]] &&
          exactSiteImageEquationQ[
            siteOrbits[[orbit]],
            spatialActions[[operation]],
            source,
            images[[orbit, operation, source]],
            shifts[[orbit, operation, source]]
          ],
        {source, Length[siteOrbits[[orbit]]]}],
    {orbit, orbitCount}, {operation, operationCount}
  ]
];

SitePermutationDataQ[data_] := Module[{cached, result},
  cached = Lookup[
    sitePermutationValidationCache, Hash[data], Missing["NotCached"]];
  If[!MissingQ[cached], Return[cached]];
  result = TrueQ[validSitePermutationDataAssociationQ[data]];
  AssociateTo[sitePermutationValidationCache, Hash[data] -> result];
  result
];

siteImageRecord[siteOrbit_List, operation_Association, source_Integer] := Module[
  {transformed, matches, target, shift},
  transformed = FullSimplify[
    operation["Rotation"].siteOrbit[[source]] + operation["Translation"]];
  matches = Select[Range[Length[siteOrbit]],
    MatrixPredicates`IntegerVectorQ[
      FullSimplify[transformed - siteOrbit[[#]]]] &];
  If[Length[matches] =!= 1, Return[$Failed]];
  target = First[matches];
  shift = FullSimplify[transformed - siteOrbit[[target]]];
  <|"SourceSite" -> source, "TargetSite" -> target,
    "CellTranslation" -> shift|>
];

CompileSitePermutationData::sites =
  "siteOrbits must be nonempty finite coordinate orbits with a common dimension.";
CompileSitePermutationData::spatial =
  "Every spatial action must have a valid Rotation and Translation.";
CompileSitePermutationData::image =
  "Operation `2` does not map orbit `1` bijectively to itself modulo integer translations.";

CompileSitePermutationData[siteOrbits_List, spatialActions_List] := Module[
  {dimension, records, images, shifts, failure, candidate},
  If[!ValidSiteOrbitsQ[siteOrbits],
    Message[CompileSitePermutationData::sites]; Return[$Failed]
  ];
  dimension = Length[siteOrbits[[1, 1]]];
  If[spatialActions === {} ||
      !And @@ (SymmetryAlgebra`ValidSpatialOperationQ[#, dimension] & /@
        spatialActions),
    Message[CompileSitePermutationData::spatial]; Return[$Failed]
  ];
  records = Table[
    siteImageRecord[siteOrbits[[orbit]], spatialActions[[operation]], source],
    {orbit, Length[siteOrbits]},
    {operation, Length[spatialActions]},
    {source, Length[siteOrbits[[orbit]]]}
  ];
  If[MemberQ[records, $Failed, Infinity],
    failure = FirstPosition[records, $Failed];
    Message[CompileSitePermutationData::image, failure[[1]], failure[[2]]];
    Return[$Failed]
  ];
  images = Map[Lookup[#, "TargetSite"] &, records, {2}];
  If[!And @@ Flatten@Map[Sort[#] === Range[Length[#]] &, images, {2}],
    Message[CompileSitePermutationData::image, "unknown", "unknown"];
    Return[$Failed]
  ];
  shifts = Map[Lookup[#, "CellTranslation"] &, records, {2}];
  candidate = <|
    "Schema" -> "SitePermutationData",
    "SchemaVersion" -> 1,
    "SiteOrbits" -> siteOrbits,
    "SpatialActions" -> spatialActions,
    "ImageRecords" -> records,
    "ImageSiteIndices" -> images,
    "CellTranslations" -> shifts,
    "Verified" -> True
  |>;
  AssociateTo[sitePermutationValidationCache, Hash[candidate] -> True];
  candidate
];
CompileSitePermutationData[___] :=
  (Message[CompileSitePermutationData::sites]; $Failed);

SiteSymmetryOperationIndices::data =
  "Expected verified SitePermutationData and valid orbit/site indices.";
SiteSymmetryOperationIndices[
    data_Association, orbit_Integer?Positive, site_Integer?Positive
  ] := Module[{images},
  images = Lookup[data, "ImageSiteIndices", $Failed];
  If[!SitePermutationDataQ[data] ||
      orbit > Length[images] || site > Length[images[[orbit, 1]]],
    Message[SiteSymmetryOperationIndices::data]; Return[$Failed]
  ];
  Pick[
    Range[Length[images[[orbit]]]],
    images[[orbit]][[All, site]],
    site
  ]
];
SiteSymmetryOperationIndices[___] :=
  (Message[SiteSymmetryOperationIndices::data]; $Failed);

End[]

EndPackage[]
