(* ::Package:: *)

BeginPackage["SymmetryAlgebra`"]

ValidSpatialOperationQ::usage = "ValidSpatialOperationQ[operation,dimension] validates an affine spatial action.";
SpatialOperationProduct::usage = "SpatialOperationProduct[left,right] composes affine spatial actions.";
SpatialOperationDifference::usage = "SpatialOperationDifference[first,second] returns their translation difference when rotations agree.";
SeitzOperationProduct::usage = "SeitzOperationProduct[left,right] composes Seitz records including antiunitary parity.";
EquivalentSeitzOperationQ::usage = "EquivalentSeitzOperationQ[first,second] compares Seitz records modulo integer translations.";
SeitzOperationDifference::usage = "SeitzOperationDifference[first,second] returns a compatible translation difference.";
CompileSeitzMultiplicationTable::usage = "CompileSeitzMultiplicationTable[operations] builds operation-index and lattice-shift records.";
CompileSeitzGroupData::usage = "CompileSeitzGroupData[operations] returns a verified abstract GroupAlgebra.";
ValidSpinSpaceGroupElementQ::usage = "ValidSpinSpaceGroupElementQ[element] validates the discrete SSG dictionary form.";
SpinSpaceGroupElementProduct::usage = "SpinSpaceGroupElementProduct[left,right] multiplies discrete SSG elements.";
EquivalentSpinSpaceGroupElementQ::usage = "EquivalentSpinSpaceGroupElementQ[first,second] compares SSG elements modulo integer translation.";
CompileSpinSpaceGroupData::usage = "CompileSpinSpaceGroupData[elements] returns a verified abstract GroupAlgebra.";

Begin["`Private`"]

ValidSpatialOperationQ[operation_Association, dimension_Integer?Positive] := Module[
  {rotation = Lookup[operation, "Rotation", $Failed],
   translation = Lookup[operation, "Translation", $Failed]},
  MatrixQ[rotation] && Dimensions[rotation] === {dimension, dimension} &&
    !MatrixPredicates`ExactZeroExpressionQ[Det[rotation]] &&
    VectorQ[translation] && Length[translation] === dimension
];
ValidSpatialOperationQ[_, _] := False;

SpatialOperationProduct[left_Association, right_Association] := <|
  "Rotation" -> FullSimplify[left["Rotation"].right["Rotation"]],
  "Translation" -> FullSimplify[
    left["Rotation"].right["Translation"] + left["Translation"]]
|>;

SpatialOperationDifference[first_Association, second_Association] := If[
  MatrixPredicates`ExactMatrixEqualQ[first["Rotation"], second["Rotation"]],
  FullSimplify[first["Translation"] - second["Translation"]],
  $Failed
];

SeitzOperationProduct[left_Association, right_Association] := Module[
  {leftFlag = Lookup[left, "Antiunitary", False],
   rightFlag = Lookup[right, "Antiunitary", False]},
  If[!BooleanQ[leftFlag] || !BooleanQ[rightFlag], Return[$Failed]];
  Join[SpatialOperationProduct[left, right],
    <|"Antiunitary" -> Xor[leftFlag, rightFlag]|>]
];
SeitzOperationProduct[___] := $Failed;

EquivalentSeitzOperationQ[first_Association, second_Association] :=
  BooleanQ[Lookup[first, "Antiunitary", False]] &&
  BooleanQ[Lookup[second, "Antiunitary", False]] &&
  SameQ[Lookup[first, "Antiunitary", False],
    Lookup[second, "Antiunitary", False]] &&
  MatrixPredicates`ExactMatrixEqualQ[first["Rotation"], second["Rotation"]] &&
  MatrixPredicates`IntegerVectorQ[
    FullSimplify[first["Translation"] - second["Translation"]]];
EquivalentSeitzOperationQ[_, _] := False;

SeitzOperationDifference[first_Association, second_Association] := If[
  BooleanQ[Lookup[first, "Antiunitary", False]] &&
    BooleanQ[Lookup[second, "Antiunitary", False]] &&
    SameQ[Lookup[first, "Antiunitary", False],
      Lookup[second, "Antiunitary", False]],
  SpatialOperationDifference[first, second],
  $Failed
];

canonicalSeitzKey[operation_Association] := ToString[
  {
    FullSimplify[operation["Rotation"]],
    FullSimplify[Mod[#, 1]] & /@ FullSimplify[operation["Translation"]],
    TrueQ[Lookup[operation, "Antiunitary", False]]
  },
  InputForm
];

CompileSeitzMultiplicationTable[operations_List] := Module[
  {keys, indexMap, records, product, index},
  If[operations === {} ||
      !And @@ (ValidSpatialOperationQ[
        #, Length[Lookup[#, "Translation", {}]]] &&
          BooleanQ[Lookup[#, "Antiunitary", False]] & /@ operations),
    Return[$Failed]
  ];
  keys = canonicalSeitzKey /@ operations;
  If[!DuplicateFreeQ[keys], Return[$Failed]];
  indexMap = AssociationThread[keys, Range[Length[keys]]];
  records = Table[
    product = SeitzOperationProduct[operations[[left]], operations[[right]]];
    index = Lookup[indexMap, canonicalSeitzKey[product], Missing["NotFound"]];
    If[MissingQ[index], $Failed,
      <|
        "OperationIndex" -> index,
        "LatticeTranslation" -> SeitzOperationDifference[product, operations[[index]]]
      |>
    ],
    {left, Length[operations]}, {right, Length[operations]}
  ];
  If[MemberQ[records, $Failed, Infinity], $Failed, records]
];
CompileSeitzMultiplicationTable[_] := $Failed;

CompileSeitzGroupData[operations_List] := Module[{records, table, group},
  records = CompileSeitzMultiplicationTable[operations];
  If[records === $Failed, Return[$Failed]];
  table = Map[Lookup[#, "OperationIndex", $Failed] &, records, {2}];
  group = GroupAlgebra`CreateGroupAlgebra[table];
  If[group === $Failed, $Failed, <|"GroupAlgebra" -> group|>]
];
CompileSeitzGroupData[_] := $Failed;

ValidSpinSpaceGroupElementQ[element_Association] := Module[
  {space = Lookup[element, "space", $Failed], spin = Lookup[element, "spin", $Failed]},
  MatchQ[space, {_?MatrixQ, _List}] && Dimensions[space[[1]]] === {3, 3} &&
    Length[space[[2]]] === 3 &&
    !MatrixPredicates`ExactZeroExpressionQ[Det[space[[1]]]] &&
  MatchQ[spin, {_?MatrixQ, 0 | 1}] && Dimensions[spin[[1]]] === {3, 3} &&
    !MatrixPredicates`ExactZeroExpressionQ[Det[spin[[1]]]]
];
ValidSpinSpaceGroupElementQ[_] := False;

SpinSpaceGroupElementProduct[left_Association, right_Association] := Module[
  {leftSpace, rightSpace, leftSpin, rightSpin},
  If[!ValidSpinSpaceGroupElementQ[left] || !ValidSpinSpaceGroupElementQ[right],
    Return[$Failed]
  ];
  leftSpace = left["space"]; rightSpace = right["space"];
  leftSpin = left["spin"]; rightSpin = right["spin"];
  <|
    "space" -> {
      FullSimplify[leftSpace[[1]].rightSpace[[1]]],
      FullSimplify[Mod[leftSpace[[1]].rightSpace[[2]] + leftSpace[[2]], 1]]
    },
    "spin" -> {
      FullSimplify[leftSpin[[1]].rightSpin[[1]]],
      Mod[leftSpin[[2]] + rightSpin[[2]], 2]
    }
  |>
];
SpinSpaceGroupElementProduct[___] := $Failed;

EquivalentSpinSpaceGroupElementQ[first_Association, second_Association] :=
  ValidSpinSpaceGroupElementQ[first] && ValidSpinSpaceGroupElementQ[second] &&
    MatrixPredicates`ExactMatrixEqualQ[first["space"][[1]], second["space"][[1]]] &&
    MatrixPredicates`IntegerVectorQ[
      FullSimplify[first["space"][[2]] - second["space"][[2]]]] &&
    MatrixPredicates`ExactMatrixEqualQ[first["spin"][[1]], second["spin"][[1]]] &&
    first["spin"][[2]] === second["spin"][[2]];
EquivalentSpinSpaceGroupElementQ[___] := False;

CompileSpinSpaceGroupData[elements_List] := Module[
  {order, productTable, product, matches, group},
  If[elements === {} || !And @@ (ValidSpinSpaceGroupElementQ /@ elements),
    Return[$Failed]
  ];
  order = Length[elements];
  productTable = Table[
    product = SpinSpaceGroupElementProduct[elements[[left]], elements[[right]]];
    matches = Select[Range[order], EquivalentSpinSpaceGroupElementQ[product, elements[[#]]] &];
    If[Length[matches] === 1, First[matches], $Failed],
    {left, order}, {right, order}
  ];
  If[MemberQ[productTable, $Failed, Infinity], Return[$Failed]];
  group = GroupAlgebra`CreateGroupAlgebra[productTable];
  If[group === $Failed, $Failed, <|"GroupAlgebra" -> group|>]
];
CompileSpinSpaceGroupData[_] := $Failed;

End[]

EndPackage[]
