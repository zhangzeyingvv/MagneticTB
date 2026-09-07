(* ::Package:: *)

BeginPackage["GroupAlgebra`"]

GroupActionTableQ::usage =
  "GroupActionTableQ[group,imageTable] checks an operation-by-source integer action table.";
CompileGroupAction::usage =
  "CompileGroupAction[group,imageTable] validates a finite group action once and returns an opaque action object consumed by orbit, stabilizer, and transporter functions.";
ActionOrbitIndices::usage = "ActionOrbitIndices[group,table,source] returns one orbit.";
ActionOrbits::usage = "ActionOrbits[group,table] partitions the acted-on finite set.";
StabilizerIndices::usage = "StabilizerIndices[group,table,source] returns fixing operation indices.";
TransporterIndices::usage = "TransporterIndices[group,table,source,target] returns mapping operation indices.";

Begin["`Private`"]

ClearAll[
  groupActionValidationCache,
  validGroupActionTableQ,
  verifiedGroupAction,
  compiledGroupActionQ,
  validCompiledActionIndexQ
];
groupActionValidationCache = <||>;

validGroupActionTableQ[group_Association, imageTable_List] := Module[
  {order, objectCount, indices, table},
  If[!GroupAlgebraQ[group] || imageTable === {}, Return[False]];
  order = group["Order"];
  objectCount = Length[First[imageTable]];
  indices = Range[objectCount];
  table = group["ProductTable"];
  Dimensions[imageTable] === {order, objectCount} && objectCount > 0 &&
    And @@ (Sort[#] === indices & /@ imageTable) &&
    imageTable[[group["IdentityIndex"]]] === indices &&
    And @@ Flatten@Table[
      imageTable[[table[[left, right]], source]] ===
        imageTable[[left, imageTable[[right, source]]]],
      {left, order}, {right, order}, {source, objectCount}
  ]
];
validGroupActionTableQ[_, _] := False;

GroupActionTableQ[group_Association, imageTable_List] := Module[
  {key, cached, result},
  key = Hash[{group, imageTable}];
  cached = Lookup[groupActionValidationCache, key, Missing["NotCached"]];
  If[!MissingQ[cached], Return[cached]];
  result = TrueQ[validGroupActionTableQ[group, imageTable]];
  AssociateTo[groupActionValidationCache, key -> result];
  result
];
GroupActionTableQ[_, _] := False;

CompileGroupAction::data =
  "Expected a verified group algebra and a valid operation-by-source action table.";

CompileGroupAction[group_Association, imageTable_List] := If[
  GroupActionTableQ[group, imageTable],
  verifiedGroupAction[group, imageTable, Length[First[imageTable]]],
  Message[CompileGroupAction::data];
  $Failed
];
CompileGroupAction[_, _] := (Message[CompileGroupAction::data]; $Failed);

compiledGroupActionQ[
    verifiedGroupAction[group_Association, table_List, objectCount_Integer]
  ] := objectCount > 0 &&
  Dimensions[table] === {Lookup[group, "Order", 0], objectCount};
compiledGroupActionQ[_] := False;

validCompiledActionIndexQ[action_, index_] :=
  compiledGroupActionQ[action] && IntegerQ[index] &&
    Between[index, {1, action[[3]]}];

ActionOrbitIndices::data =
  "Expected a compiled group action and valid source index, or a verified group with a valid action table.";
ActionOrbitIndices[action_verifiedGroupAction, source_Integer] :=
  If[validCompiledActionIndexQ[action, source],
    Sort@DeleteDuplicates[action[[2, All, source]]],
    Message[ActionOrbitIndices::data]; $Failed
  ];
ActionOrbitIndices[group_Association, table_List, source_Integer] := Module[
  {action = Quiet@CompileGroupAction[group, table]},
  If[action === $Failed,
    Message[ActionOrbitIndices::data]; $Failed,
    ActionOrbitIndices[action, source]
  ]
];
ActionOrbitIndices[___] := (Message[ActionOrbitIndices::data]; $Failed);

ActionOrbits::data = ActionOrbitIndices::data;
ActionOrbits[action_verifiedGroupAction] := Module[
  {unseen, orbits = {}, orbit},
  If[!compiledGroupActionQ[action],
    Message[ActionOrbits::data]; Return[$Failed]
  ];
  unseen = Range[action[[3]]];
  While[unseen =!= {},
    orbit = ActionOrbitIndices[action, First[unseen]];
    If[orbit === $Failed, Return[$Failed]];
    AppendTo[orbits, orbit]; unseen = Complement[unseen, orbit]
  ];
  orbits
];
ActionOrbits[group_Association, table_List] := Module[
  {action = Quiet@CompileGroupAction[group, table]},
  If[action === $Failed,
    Message[ActionOrbits::data]; $Failed,
    ActionOrbits[action]
  ]
];
ActionOrbits[___] := (Message[ActionOrbits::data]; $Failed);

StabilizerIndices::data = ActionOrbitIndices::data;
StabilizerIndices[action_verifiedGroupAction, source_Integer] :=
  If[validCompiledActionIndexQ[action, source],
    Pick[Range[action[[1]]["Order"]], action[[2, All, source]], source],
    Message[StabilizerIndices::data]; $Failed
  ];
StabilizerIndices[group_Association, table_List, source_Integer] := Module[
  {action = Quiet@CompileGroupAction[group, table]},
  If[action === $Failed,
    Message[StabilizerIndices::data]; $Failed,
    StabilizerIndices[action, source]
  ]
];
StabilizerIndices[___] := (Message[StabilizerIndices::data]; $Failed);

TransporterIndices::data =
  "Expected a compiled group action and valid source/target indices, or a verified group with a valid action table.";
TransporterIndices[
    action_verifiedGroupAction,
    source_Integer,
    target_Integer
  ] := If[
  validCompiledActionIndexQ[action, source] &&
    Between[target, {1, action[[3]]}],
  Pick[Range[action[[1]]["Order"]], action[[2, All, source]], target],
  Message[TransporterIndices::data]; $Failed
];
TransporterIndices[
    group_Association,
    table_List,
    source_Integer,
    target_Integer
  ] := Module[{action = Quiet@CompileGroupAction[group, table]},
  If[action === $Failed,
    Message[TransporterIndices::data]; $Failed,
    TransporterIndices[action, source, target]
  ]
];
TransporterIndices[___] := (Message[TransporterIndices::data]; $Failed);

End[]

EndPackage[]
