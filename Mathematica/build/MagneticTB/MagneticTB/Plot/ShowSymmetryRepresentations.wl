(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  resolveSymmetryRepresentationIndices,
  validSymmetryRepresentationDataQ,
  buildSymmetryRepresentationDisplayData,
  renderSymmetryRepresentationGrid
];

showSymmetryRepresentations::selection =
  "The operation selection must be Automatic, All, a positive integer, or a nonempty list of positive integers; received `1`.";
showSymmetryRepresentations::range =
  "Operation indices `1` are outside the prepared range 1 through `2`.";
showSymmetryRepresentations::data =
  "The cached full representation is missing or malformed.";

resolveSymmetryRepresentationIndices[
    selection_,
    operationCount_Integer?Positive,
    generatorIndices_List
  ] := Module[{indices},
  indices = Which[
    selection === Automatic,
      generatorIndices,
    selection === All,
      Range[operationCount],
    IntegerQ[selection] && selection > 0,
      {selection},
    ListQ[selection] && selection =!= {} &&
        And @@ (IntegerQ[#] && # > 0 & /@ selection),
      DeleteDuplicates[selection],
    True,
      Return[Failure["InvalidSelection", <|"Selection" -> selection|>]]
  ];
  If[
    indices === {} || !SubsetQ[Range[operationCount], indices],
    Failure[
      "SelectionOutOfRange",
      <|
        "Selection" -> indices,
        "OperationCount" -> operationCount
      |>
    ],
    indices
  ]
];

validSymmetryRepresentationDataQ[
    spatialActions_,
    spinActions_,
    matrices_,
    antiunitaryFlags_,
    labels_
  ] := Module[{operationCount, dimensions},
  If[
    !ListQ[spatialActions] || spatialActions === {} ||
      !ListQ[spinActions] ||
      !ListQ[matrices] || !ListQ[antiunitaryFlags] || !ListQ[labels],
    Return[False]
  ];
  operationCount = Length[spatialActions];
  If[
    Length[spinActions] =!= operationCount ||
      Length[matrices] =!= operationCount ||
      Length[antiunitaryFlags] =!= operationCount ||
      Length[labels] =!= operationCount ||
      !And @@ (AssociationQ /@ spatialActions) ||
      !And @@ (MatrixQ[#] && Dimensions[#] === {3, 3} & /@ spinActions) ||
      !And @@ (MatrixQ /@ matrices) ||
      !And @@ (MemberQ[{True, False}, #] & /@ antiunitaryFlags),
    Return[False]
  ];
  dimensions = Dimensions /@ matrices;
  SameQ @@ dimensions &&
    Length[First[dimensions]] == 2 &&
    SameQ @@ First[dimensions] &&
    And @@ (
      MatrixQ[Lookup[#, "Rotation", $Failed]] &&
        ListQ[Lookup[#, "Translation", $Failed]] & /@ spatialActions
    )
];

(* Pure data layer: this function only selects already assembled matrices from
   CurrentModelSession.  It performs no representation construction, group
   closure, simplification, or validation. *)
buildSymmetryRepresentationDisplayData[
    session_Association,
    selection_
  ] := Module[
  {
    representation, model, symmetry, metadata,
    spatialActions, spinActions, matrices,
    antiunitaryFlags, labels, generatorIndices, operationCount,
    selectedIndices, records
  },
  representation = Lookup[session, "FullRepresentation", $Failed];
  model = Lookup[session, "ModelSpecification", $Failed];
  symmetry = If[AssociationQ[model], Lookup[model, "Symmetry", $Failed], $Failed];
  metadata = If[AssociationQ[model], Lookup[model, "Metadata", $Failed], $Failed];
  generatorIndices = Lookup[session, "GeneratorIndices", $Failed];
  If[
    !AssociationQ[representation] ||
      !AssociationQ[symmetry] || !AssociationQ[metadata] ||
      !ListQ[generatorIndices],
    Return[$Failed]
  ];
  spatialActions = Lookup[
    symmetry,
    "SpatialActions",
    $Failed
  ];
  spinActions = Lookup[symmetry, "SpinActions", $Failed];
  matrices = Lookup[representation, "RepresentationMatrices", $Failed];
  antiunitaryFlags = Lookup[
    symmetry,
    "AntiunitaryFlags",
    $Failed
  ];
  labels = Lookup[metadata, "SymmetryInformation", $Failed];
  If[ListQ[labels], labels = labels[[All, 1]]];
  If[
    !validSymmetryRepresentationDataQ[
      spatialActions,
      spinActions,
      matrices,
      antiunitaryFlags,
      labels
    ],
    Return[$Failed]
  ];
  operationCount = Length[spatialActions];
  If[
    generatorIndices === {} ||
      !And @@ (IntegerQ[#] && Between[#, {1, operationCount}] & /@
        generatorIndices),
    Return[$Failed]
  ];
  selectedIndices = resolveSymmetryRepresentationIndices[
    selection,
    operationCount,
    generatorIndices
  ];
  If[FailureQ[selectedIndices], Return[selectedIndices]];

  records = Map[
    Function[index,
      <|
        "OperationIndex" -> index,
        "Label" -> labels[[index]],
        "Generator" -> MemberQ[generatorIndices, index],
        "Rotation" -> spatialActions[[index, "Rotation"]],
        "Translation" -> spatialActions[[index, "Translation"]],
        "SpinAction" -> spinActions[[index]],
        "Antiunitary" -> antiunitaryFlags[[index]],
        "Matrix" -> matrices[[index]]
      |>
    ],
    selectedIndices
  ];

  <|
    "Schema" -> "MagneticTBSymmetryRepresentationDisplayData",
    "SchemaVersion" -> 1,
    "RepresentationMode" -> Lookup[
      session,
      "RepresentationMode",
      Missing["RepresentationMode"]
    ],
    "Convention" -> Lookup[
      representation,
      "Convention",
      Missing["Convention"]
    ],
    "Dimension" -> First[Dimensions[First[matrices]]],
    "OperationCount" -> operationCount,
    "GeneratorIndices" -> generatorIndices,
    "SelectedIndices" -> selectedIndices,
    "UnitaryVerified" -> TrueQ[
      Lookup[representation, "UnitaryVerified", False]
    ],
    "Records" -> records
  |>
];

renderSymmetryRepresentationGrid[data_Association] := Module[
  {title, rows},
  title = Row[{
    "Representation mode: ", data["RepresentationMode"],
    "; dimension = ", data["Dimension"],
    "; showing ", Length[data["SelectedIndices"]],
    " of ", data["OperationCount"], " operations"
  }];
  rows = Map[
    Function[record,
      {
        record["OperationIndex"],
        record["Label"],
        If[TrueQ[record["Generator"]], "Yes", ""],
        MatrixForm[record["Rotation"]],
        record["Translation"],
        MatrixForm[record["SpinAction"]],
        If[TrueQ[record["Antiunitary"]], "Antiunitary", "Unitary"],
        MatrixForm[record["Matrix"]]
      }
    ],
    data["Records"]
  ];
  Grid[
    Join[
      {{
        title,
        SpanFromLeft,
        SpanFromLeft,
        SpanFromLeft,
        SpanFromLeft,
        SpanFromLeft,
        SpanFromLeft,
        SpanFromLeft
      }},
      {{
        "Index", "Label", "Generator", "R", "t", "SpinAction", "Type",
        "Matrix"
      }},
      rows
    ],
    Frame -> All,
    Alignment -> Left
  ]
];

showSymmetryRepresentations[selection_: Automatic] := Module[
  {session, data},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  data = buildSymmetryRepresentationDisplayData[session, selection];
  Which[
    MatchQ[data, Failure["InvalidSelection", _Association]],
      Message[showSymmetryRepresentations::selection, selection];
      $Failed,
    MatchQ[data, Failure["SelectionOutOfRange", _Association]],
      Message[
        showSymmetryRepresentations::range,
        data[[2, "Selection"]],
        data[[2, "OperationCount"]]
      ];
      $Failed,
    data === $Failed,
      Message[showSymmetryRepresentations::data];
      $Failed,
    True,
      renderSymmetryRepresentationGrid[data]
  ]
];

End[]
EndPackage[]
