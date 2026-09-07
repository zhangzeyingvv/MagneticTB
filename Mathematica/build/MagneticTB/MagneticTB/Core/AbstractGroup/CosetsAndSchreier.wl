(* ::Package:: *)

BeginPackage["GroupAlgebra`"]

RightCosetIndices::usage =
  "RightCosetIndices[data,subgroup,representative] returns representative*subgroup.";
LeftCosetIndices::usage =
  "LeftCosetIndices[data,subgroup,representative] returns subgroup*representative.";
SchreierDecompositionIndex::usage =
  "SchreierDecompositionIndex[data,subgroup,reps,g,source] returns {target,h} with g reps[[source]] = reps[[target]] h.";

Begin["`Private`"]

closedSubgroupQ[data_, subgroup_List] :=
  FindGeneratorIndices[data, subgroup] =!= $Failed;

RightCosetIndices::data =
  "Expected a verified group algebra, closed subgroup, and valid representative.";
RightCosetIndices[data_Association, subgroupIndices_List, representative_Integer] := Module[
  {subgroup = Sort@DeleteDuplicates[subgroupIndices]},
  If[!GroupAlgebraQ[data] || !Between[representative, {1, data["Order"]}] ||
      !closedSubgroupQ[data, subgroup],
    Message[RightCosetIndices::data]; Return[$Failed]
  ];
  data["ProductTable"][[representative, subgroup]]
];
RightCosetIndices[___] := (Message[RightCosetIndices::data]; $Failed);

LeftCosetIndices::data = RightCosetIndices::data;
LeftCosetIndices[data_Association, subgroupIndices_List, representative_Integer] := Module[
  {subgroup = Sort@DeleteDuplicates[subgroupIndices]},
  If[!GroupAlgebraQ[data] || !Between[representative, {1, data["Order"]}] ||
      !closedSubgroupQ[data, subgroup],
    Message[LeftCosetIndices::data]; Return[$Failed]
  ];
  data["ProductTable"][[subgroup, representative]]
];
LeftCosetIndices[___] := (Message[LeftCosetIndices::data]; $Failed);

SchreierDecompositionIndex::data =
  "Expected valid subgroup, coset representatives, operation, and source indices.";
SchreierDecompositionIndex[
    data_Association, subgroupIndices_List, representatives_List,
    operation_Integer, source_Integer
  ] := Module[{table, leftProduct, matches},
  If[!GroupAlgebraQ[data] || !closedSubgroupQ[data, subgroupIndices] ||
      !And @@ (IntegerQ[#] && Between[#, {1, data["Order"]}] & /@ representatives) ||
      !Between[operation, {1, data["Order"]}] ||
      !Between[source, {1, Length[representatives]}],
    Message[SchreierDecompositionIndex::data]; Return[$Failed]
  ];
  table = data["ProductTable"];
  leftProduct = table[[operation, representatives[[source]]]];
  matches = Cases[
    Tuples[{Range[Length[representatives]], subgroupIndices}],
    {target_, h_} /; table[[representatives[[target]], h]] === leftProduct
  ];
  If[Length[matches] =!= 1,
    Message[SchreierDecompositionIndex::data]; Return[$Failed]
  ];
  First[matches]
];
SchreierDecompositionIndex[___] :=
  (Message[SchreierDecompositionIndex::data]; $Failed);

End[]

EndPackage[]
