(* ::Package:: *)

BeginPackage["GroupAlgebra`"]

GeneratedSubgroupIndices::usage =
  "GeneratedSubgroupIndices[data,generators] returns the generated subgroup indices.";
FindGeneratorIndices::usage =
  "FindGeneratorIndices[data] finds generators of the full group; FindGeneratorIndices[data,subgroup] finds generators of a closed subgroup.";

Begin["`Private`"]

GeneratedSubgroupIndices::data =
  "Expected a verified group algebra and valid generator indices.";
GeneratedSubgroupIndices[data_Association, generatorIndices_List] := Module[
  {order, table, generators, seen, queue, head = 1, current, products},
  If[!GroupAlgebraQ[data],
    Message[GeneratedSubgroupIndices::data]; Return[$Failed]
  ];
  order = data["Order"];
  generators = DeleteDuplicates[generatorIndices];
  If[!And @@ (IntegerQ[#] && Between[#, {1, order}] & /@ generators),
    Message[GeneratedSubgroupIndices::data]; Return[$Failed]
  ];
  table = data["ProductTable"];
  seen = Association[data["IdentityIndex"] -> True];
  queue = {data["IdentityIndex"]};
  While[head <= Length[queue],
    current = queue[[head++]];
    products = Join[
      table[[current, generators]],
      table[[generators, current]]
    ];
    Scan[
      Function[product,
        If[!KeyExistsQ[seen, product],
          AssociateTo[seen, product -> True]; AppendTo[queue, product]
        ]
      ],
      products
    ]
  ];
  Sort[Keys[seen]]
];
GeneratedSubgroupIndices[___] :=
  (Message[GeneratedSubgroupIndices::data]; $Failed);

FindGeneratorIndices::data =
  "Expected a verified group algebra and a closed subgroup containing the identity.";
FindGeneratorIndices[data_Association, subgroup_: All] := Module[
  {indices, table, identity, generators = {}, generated, candidateClosure},
  If[!GroupAlgebraQ[data],
    Message[FindGeneratorIndices::data]; Return[$Failed]
  ];
  indices = If[subgroup === All, Range[data["Order"]],
    Sort@DeleteDuplicates[subgroup]];
  table = data["ProductTable"];
  identity = data["IdentityIndex"];
  If[indices === {} || !MemberQ[indices, identity] ||
      !And @@ (IntegerQ[#] && Between[#, {1, data["Order"]}] & /@ indices) ||
      !SubsetQ[indices, Flatten[table[[indices, indices]]]],
    Message[FindGeneratorIndices::data]; Return[$Failed]
  ];
  If[indices === {identity}, Return[{identity}]];
  generated = {identity};
  Do[
    If[MemberQ[generated, candidate], Continue[]];
    candidateClosure = GeneratedSubgroupIndices[data, Append[generators, candidate]];
    If[candidateClosure === $Failed || !SubsetQ[indices, candidateClosure],
      Message[FindGeneratorIndices::data]; Return[$Failed]
    ];
    If[Length[candidateClosure] > Length[generated],
      AppendTo[generators, candidate]; generated = candidateClosure
    ];
    If[generated === indices, Break[]],
    {candidate, indices}
  ];
  If[generated === indices, generators,
    Message[FindGeneratorIndices::data]; $Failed]
];
FindGeneratorIndices[___] := (Message[FindGeneratorIndices::data]; $Failed);

End[]

EndPackage[]
