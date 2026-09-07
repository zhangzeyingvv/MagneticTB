(* ::Package:: *)

BeginPackage["GroupAlgebra`"]

CreateGroupAlgebra::usage =
  "CreateGroupAlgebra[productTable] validates an integer Cayley table and returns an abstract finite-group algebra.";
GroupAlgebraQ::usage = "GroupAlgebraQ[data] tests a verified group algebra.";
ProductIndex::usage = "ProductIndex[data,left,right] returns a product index.";

Begin["`Private`"]

ClearAll[groupAlgebraValidationCache, validGroupAlgebraAssociationQ];
groupAlgebraValidationCache = <||>;

validProductTableQ[table_List] := Module[{order, indices},
  order = Length[table];
  indices = Range[order];
  order > 0 && Dimensions[table] === {order, order} &&
    And @@ (IntegerQ[#] && Between[#, {1, order}] & /@ Flatten[table]) &&
    And @@ (Sort[#] === indices & /@ table) &&
    And @@ (Sort[#] === indices & /@ Transpose[table])
];
validProductTableQ[_] := False;

productTableAssociativeQ[table_List] := Module[{order = Length[table]},
  And @@ Flatten@Table[
    table[[table[[left, right]], All]] ===
      table[[left, table[[right, All]]]],
    {left, order}, {right, order}
  ]
];

CreateGroupAlgebra::table =
  "Expected a square integer Cayley table defining one associative finite group.";

CreateGroupAlgebra[productTable_List] := Module[
  {order, indices, identityIndices, identityIndex, inverseIndices,
   candidate},
  If[!validProductTableQ[productTable] ||
      !productTableAssociativeQ[productTable],
    Message[CreateGroupAlgebra::table]; Return[$Failed]
  ];
  order = Length[productTable];
  indices = Range[order];
  identityIndices = Select[
    indices,
    productTable[[#, All]] === indices &&
      productTable[[All, #]] === indices &
  ];
  If[Length[identityIndices] =!= 1,
    Message[CreateGroupAlgebra::table]; Return[$Failed]
  ];
  identityIndex = First[identityIndices];
  inverseIndices = Table[
    SelectFirst[
      indices,
      productTable[[element, #]] === identityIndex &&
        productTable[[#, element]] === identityIndex &,
      Missing["NotFound"]
    ],
    {element, order}
  ];
  If[AnyTrue[inverseIndices, MissingQ],
    Message[CreateGroupAlgebra::table]; Return[$Failed]
  ];
  candidate = <|
    "Schema" -> "GroupAlgebra",
    "Order" -> order,
    "ProductTable" -> productTable,
    "IdentityIndex" -> identityIndex,
    "InverseIndices" -> inverseIndices,
    "Verified" -> True
  |>;
  AssociateTo[groupAlgebraValidationCache, Hash[candidate] -> True];
  candidate
];
CreateGroupAlgebra[_] := (Message[CreateGroupAlgebra::table]; $Failed);

validGroupAlgebraAssociationQ[data_] := Module[
  {table, order, identity, inverses, indices},
  If[!AssociationQ[data] || Lookup[data, "Schema", None] =!= "GroupAlgebra" ||
      !TrueQ[Lookup[data, "Verified", False]], Return[False]];
  table = Lookup[data, "ProductTable", $Failed];
  order = Lookup[data, "Order", $Failed];
  identity = Lookup[data, "IdentityIndex", $Failed];
  inverses = Lookup[data, "InverseIndices", $Failed];
  If[!IntegerQ[order] || order < 1 || !ListQ[table] ||
      Length[table] =!= order || !validProductTableQ[table] ||
      !productTableAssociativeQ[table] || !IntegerQ[identity] ||
      !Between[identity, {1, order}] || !ListQ[inverses] ||
      Length[inverses] =!= order, Return[False]];
  indices = Range[order];
  table[[identity, All]] === indices && table[[All, identity]] === indices &&
    And @@ Table[
      IntegerQ[inverses[[element]]] &&
        Between[inverses[[element]], {1, order}] &&
        table[[element, inverses[[element]]]] === identity &&
        table[[inverses[[element]], element]] === identity,
      {element, order}
  ]
];

GroupAlgebraQ[data_] := Module[{cached, result},
  cached = Lookup[
    groupAlgebraValidationCache, Hash[data], Missing["NotCached"]];
  If[!MissingQ[cached], Return[cached]];
  result = TrueQ[validGroupAlgebraAssociationQ[data]];
  AssociateTo[groupAlgebraValidationCache, Hash[data] -> result];
  result
];

ProductIndex::data =
  "Expected a verified group algebra and valid left/right indices.";
ProductIndex[data_Association, left_Integer, right_Integer] := Module[
  {order = Lookup[data, "Order", 0]},
  If[!GroupAlgebraQ[data] || !Between[left, {1, order}] ||
      !Between[right, {1, order}],
    Message[ProductIndex::data]; Return[$Failed]
  ];
  data["ProductTable"][[left, right]]
];
ProductIndex[___] := (Message[ProductIndex::data]; $Failed);

End[]

EndPackage[]
