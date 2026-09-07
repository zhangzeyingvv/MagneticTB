cqMatrixQ[matrix_] := AssociationQ[matrix] &&
  Lookup[matrix, "type", None] === "exact_matrix" &&
  cqContextQ[Lookup[matrix, "context", None]] &&
  IntegerQ[Lookup[matrix, "rows", None]] && Lookup[matrix, "rows", -1] >= 0 &&
  IntegerQ[Lookup[matrix, "columns", None]] && Lookup[matrix, "columns", -1] >= 0 &&
  ListQ[Lookup[matrix, "entries", None]] &&
  Length[matrix["entries"]] === matrix["rows"] matrix["columns"] &&
  And @@ (cqElementDataQ[matrix["context"], #] & /@ matrix["entries"]);

cqMatrixCreate[context_, rows_Integer, columns_Integer, entries_List] := Module[{},
  If[!cqContextQ[context], Return[cqFailure["ConductorMismatch", <||>]]];
  If[rows < 0 || columns < 0 || Length[entries] =!= rows columns,
    Return[cqFailure[
      "MalformedSerialization",
      <|"rows" -> rows, "columns" -> columns, "entry_count" -> Length[entries]|>
    ]]
  ];
  If[!And @@ (cqElementDataQ[context, #] & /@ entries),
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  <|
    "type" -> "exact_matrix",
    "context" -> context,
    "rows" -> rows,
    "columns" -> columns,
    "entries" -> entries
  |>
];

cqMatrixEntry[matrix_, row_Integer, column_Integer] :=
  matrix["entries"][[ (row - 1) matrix["columns"] + column ]];

cqMatrixIdentity[context_, dimension_Integer] := Module[{zero, one, entries, row, column},
  If[dimension < 0, Return[cqFailure["DimensionMismatch", <|"dimension" -> dimension|>]]];
  zero = cqElementZero[context];
  one = cqElementOne[context];
  entries = Flatten[
    Table[If[row === column, one, zero], {row, 1, dimension}, {column, 1, dimension}],
    1
  ];
  cqMatrixCreate[context, dimension, dimension, entries]
];

cqMatrixZero[context_, rows_Integer, columns_Integer] :=
  cqMatrixCreate[context, rows, columns, ConstantArray[cqElementZero[context], rows columns]];

cqMatrixTranspose[matrix_] := Module[{entries, row, column},
  If[!cqMatrixQ[matrix], Return[cqFailure["MalformedSerialization", <||>]]];
  entries = Flatten[
    Table[cqMatrixEntry[matrix, row, column],
      {column, 1, matrix["columns"]}, {row, 1, matrix["rows"]}],
    1
  ];
  cqMatrixCreate[matrix["context"], matrix["columns"], matrix["rows"], entries]
];

cqMatrixMultiply[left_, right_] := Module[
  {context, entries = {}, row, column, inner, sum, product},
  If[!cqMatrixQ[left] || !cqMatrixQ[right],
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  If[left["context"]["context_id"] =!= right["context"]["context_id"],
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  If[left["columns"] =!= right["rows"],
    Return[cqFailure[
      "DimensionMismatch",
      <|"left_columns" -> left["columns"], "right_rows" -> right["rows"]|>
    ]]
  ];
  context = left["context"];
  For[row = 1, row <= left["rows"], row++,
    For[column = 1, column <= right["columns"], column++,
      sum = cqElementZero[context];
      For[inner = 1, inner <= left["columns"], inner++,
        product = cqElementMultiply[
          context,
          cqMatrixEntry[left, row, inner],
          cqMatrixEntry[right, inner, column]
        ];
        If[FailureQ[product], Return[product]];
        sum = cqElementAdd[context, sum, product];
        If[FailureQ[sum], Return[sum]];
      ];
      AppendTo[entries, sum];
    ];
  ];
  cqMatrixCreate[context, left["rows"], right["columns"], entries]
];

cqMatrixZeroQ[matrix_] := cqMatrixQ[matrix] &&
  And @@ (cqElementZeroQ[matrix["context"], #] & /@ matrix["entries"]);

cqMatrixEqualQ[left_, right_] := cqMatrixQ[left] && cqMatrixQ[right] &&
  left["rows"] === right["rows"] && left["columns"] === right["columns"] &&
  left["context"]["context_id"] === right["context"]["context_id"] &&
  And @@ MapThread[cqElementEqualQ[left["context"], #1, #2] &, {left["entries"], right["entries"]}];

cqMatrixToNested[matrix_] := Module[{row},
  If[!cqMatrixQ[matrix], Return[cqFailure["MalformedSerialization", <||>]]];
  Table[
    Table[cqRestoreElement[matrix["context"], cqMatrixEntry[matrix, row, column]],
      {column, 1, matrix["columns"]}],
    {row, 1, matrix["rows"]}
  ]
];

cqNormalizeRawMatrix[raw_List, coordinateDimension_Integer] := Module[
  {rows, rowLengths, entries},
  rows = Length[raw];
  If[rows === 0,
    Return[<|"rows" -> 0, "columns" -> coordinateDimension, "entries" -> {}|>]
  ];
  If[!And @@ (ListQ /@ raw),
    Return[cqFailure["MalformedSerialization", <|"matrix" -> HoldForm[raw]|>]]
  ];
  rowLengths = Length /@ raw;
  If[!And @@ (# === coordinateDimension & /@ rowLengths),
    Return[cqFailure[
      "DimensionMismatch",
      <|"coordinate_dimension" -> coordinateDimension, "row_lengths" -> rowLengths|>
    ]]
  ];
  entries = Flatten[raw, 1];
  <|"rows" -> rows, "columns" -> coordinateDimension, "entries" -> entries|>
];

cqNormalizeRawMatrix[raw_Association, coordinateDimension_Integer] := Module[
  {rows, columns, entries},
  rows = Lookup[raw, "rows", Missing["rows"]];
  columns = Lookup[raw, "columns", Missing["columns"]];
  entries = Lookup[raw, "entries", Missing["entries"]];
  If[!IntegerQ[rows] || rows < 0 || !IntegerQ[columns] || columns < 0 || !ListQ[entries],
    Return[cqFailure["MalformedSerialization", <|"matrix" -> raw|>]]
  ];
  If[columns =!= coordinateDimension,
    Return[cqFailure[
      "DimensionMismatch",
      <|"coordinate_dimension" -> coordinateDimension, "columns" -> columns|>
    ]]
  ];
  If[Length[entries] =!= rows columns,
    Return[cqFailure[
      "MalformedSerialization",
      <|"rows" -> rows, "columns" -> columns, "entry_count" -> Length[entries]|>
    ]]
  ];
  <|"rows" -> rows, "columns" -> columns, "entries" -> entries|>
];

cqNormalizeRawMatrix[raw_, coordinateDimension_Integer] :=
  cqFailure["MalformedSerialization", <|"matrix" -> HoldForm[raw]|>];

CompileCyclotomicMatrices[input_Association] := Module[
  {coordinateDimension, constraintsRaw, maxDegree, normalizedMatrices,
   parsedMatrices, parsedEntryLists, orders, conductor, context,
   compiledMatrices, compiledEntryLists, cachedParse, cachedCompile},
  coordinateDimension = Lookup[input, "coordinate_dimension", Missing["coordinate_dimension"]];
  constraintsRaw = Lookup[input, "constraints", Missing["constraints"]];
  maxDegree = Lookup[input, "max_cyclotomic_degree", 128];
  If[!IntegerQ[coordinateDimension] || coordinateDimension < 0 || !ListQ[constraintsRaw],
    Return[cqFailure["MalformedSerialization", <|"input" -> input|>]]
  ];
  If[!IntegerQ[maxDegree] || maxDegree < 1,
    Return[cqFailure["ConductorDegreeLimitExceeded", <|"limit" -> maxDegree|>]]
  ];
  normalizedMatrices = cqNormalizeRawMatrix[#, coordinateDimension] & /@ constraintsRaw;
  If[AnyTrue[normalizedMatrices, FailureQ],
    Return[FirstCase[normalizedMatrices, _Failure]]
  ];
  cachedParse[value_] := cachedParse[value] = cqParseRootSum[value];
  parsedEntryLists = (cachedParse /@ # ["entries"]) & /@ normalizedMatrices;
  If[AnyTrue[Flatten[parsedEntryLists, 1], FailureQ],
    Return[FirstCase[Flatten[parsedEntryLists, 1], _Failure]]
  ];
  parsedMatrices = MapThread[
    <|
      "rows" -> #1["rows"],
      "columns" -> #1["columns"],
      "entries" -> #2
    |> &,
    {normalizedMatrices, parsedEntryLists}
  ];
  orders = Flatten[
    (# ["orders"] & /@ Flatten[(# ["entries"] & /@ parsedMatrices), 1]),
    1
  ];
  conductor = If[orders === {}, 1, cqIntegerLCMList[orders]];
  context = cqCreateContext[conductor, maxDegree];
  If[FailureQ[context], Return[context]];
  cachedCompile[value_] :=
    cachedCompile[value] = cqCompileParsedRootSum[value, context];
  compiledEntryLists = (cachedCompile /@ # ["entries"]) & /@ parsedMatrices;
  If[AnyTrue[Flatten[compiledEntryLists, 1], FailureQ],
    Return[FirstCase[Flatten[compiledEntryLists, 1], _Failure]]
  ];
  compiledMatrices = MapThread[
    cqMatrixCreate[context, #1["rows"], #1["columns"], #2] &,
    {parsedMatrices, compiledEntryLists}
  ];
  If[AnyTrue[compiledMatrices, FailureQ],
    Return[FirstCase[compiledMatrices, _Failure]]
  ];
  <|
    "type" -> "compiled_cyclotomic_problem",
    "context" -> context,
    "conductor" -> conductor,
    "coordinate_dimension" -> coordinateDimension,
    "constraint_count" -> Length[compiledMatrices],
    "constraints" -> compiledMatrices
  |>
];

CompileCyclotomicMatrices[input_] :=
  cqFailure["MalformedSerialization", <|"input" -> HoldForm[input]|>];
