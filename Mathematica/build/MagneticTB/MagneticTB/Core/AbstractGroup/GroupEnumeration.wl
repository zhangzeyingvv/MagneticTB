(* ::Package:: *)

BeginPackage["GroupEnumeration`"]

GenerateGroup::usage =
  "GenerateGroup[generators,identity,multiply] enumerates a finite group of opaque elements. SameTest supplies exact element equality.";
FindConcreteGeneratorElements::usage =
  "FindConcreteGeneratorElements[elements,identity,multiply] greedily selects concrete generators using the same exact multiplication and equality callbacks.";

Begin["`Private`"]

Options[GenerateGroup] = {SameTest -> SameQ};
GenerateGroup::input = "Expected a generator list, an identity element, and callable multiplication/equality functions.";
GenerateGroup::order =
  "A supplied generator did not return to the identity within 200 products.";

GenerateGroup[
    generators_List,
    identityElement_,
    multiply_,
    OptionsPattern[]
  ] := Module[
  {i, j, generatorCount, maxOrder = 200, orders, cyclicSubgroups,
   multiplicationList, firstLayer, nextLayer, newElements, sameTest,
   memberQ, unique, setDifference, identityPosition},
  sameTest = OptionValue[SameTest];
  memberQ[list_, element_] :=
    AnyTrue[list, TrueQ[sameTest[#, element]] &];
  unique[list_] := Fold[
    If[memberQ[#1, #2], #1, Append[#1, #2]] &,
    {},
    list
  ];
  setDifference[first_, second_] :=
    Select[first, !memberQ[second, #] &];
  If[generators === {}, Return[{identityElement}]];
  generatorCount = Length[generators];
  orders = ConstantArray[0, generatorCount];
  cyclicSubgroups = ConstantArray[{}, generatorCount];
  Do[
    multiplicationList = FoldList[
      multiply,
      Table[generators[[i]], maxOrder]
    ];
    identityPosition = SelectFirst[
      Range[Length[multiplicationList]],
      TrueQ[sameTest[multiplicationList[[#]], identityElement]] &,
      Missing["NotFound"]
    ];
    If[MissingQ[identityPosition],
      Message[GenerateGroup::order];
      Return[$Failed]
    ];
    orders[[i]] = identityPosition;
    cyclicSubgroups[[i]] = multiplicationList[[;; identityPosition]],
    {i, generatorCount}
  ];
  firstLayer = unique[Flatten[cyclicSubgroups, 1]];
  nextLayer = unique@Flatten@Table[
    multiply[firstLayer[[i]], firstLayer[[j]]],
    {i, Length[firstLayer]},
    {j, Length[firstLayer]}
  ];
  newElements = setDifference[nextLayer, firstLayer];
  firstLayer = nextLayer;
  While[newElements =!= {},
    nextLayer = unique@Flatten@Table[
      multiply[newElements[[i]], firstLayer[[j]]],
      {i, Length[newElements]},
      {j, Length[firstLayer]}
    ];
    newElements = setDifference[nextLayer, firstLayer];
    firstLayer = nextLayer
  ];
  SortBy[nextLayer, !TrueQ[sameTest[#, identityElement]] &]
];
GenerateGroup[___] := (Message[GenerateGroup::input]; $Failed);

Options[FindConcreteGeneratorElements] = Options[GenerateGroup];
FindConcreteGeneratorElements::input =
  "Expected a complete finite element list, identity, and multiplication/equality callbacks.";

FindConcreteGeneratorElements[
    elements_List,
    identityElement_,
    multiply_,
    OptionsPattern[]
  ] := Module[
  {sameTest, memberQ, setEqualQ, generators = {}, generated, candidateGenerated},
  sameTest = OptionValue[SameTest];
  memberQ[list_, element_] :=
    AnyTrue[list, TrueQ[sameTest[#, element]] &];
  setEqualQ[first_, second_] :=
    Length[first] === Length[second] &&
      And @@ (memberQ[second, #] & /@ first);
  If[elements === {} || !memberQ[elements, identityElement],
    Message[FindConcreteGeneratorElements::input];
    Return[$Failed]
  ];
  If[Length[elements] === 1, Return[{identityElement}]];
  generated = {identityElement};
  Do[
    If[memberQ[generated, candidate], Continue[]];
    candidateGenerated = GenerateGroup[
      Append[generators, candidate],
      identityElement,
      multiply,
      SameTest -> sameTest
    ];
    If[candidateGenerated === $Failed, Return[$Failed]];
    If[Length[candidateGenerated] > Length[generated],
      AppendTo[generators, candidate];
      generated = candidateGenerated
    ];
    If[setEqualQ[elements, generated], Break[]],
    {candidate, elements}
  ];
  If[setEqualQ[elements, generated], generators,
    Message[FindConcreteGeneratorElements::input];
    $Failed
  ]
];
FindConcreteGeneratorElements[___] :=
  (Message[FindConcreteGeneratorElements::input]; $Failed);

End[]

EndPackage[]
