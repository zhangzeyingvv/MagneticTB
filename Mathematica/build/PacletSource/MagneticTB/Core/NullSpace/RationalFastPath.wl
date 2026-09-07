(* Exact N=1 fast path.  It keeps native Mathematica rationals during linear
   algebra, while preserving the public cyclotomic matrix/result schema. *)

cqRationalRawMatrix[rows_Integer, columns_Integer, data_List] := Module[{},
  If[rows < 0 || columns < 0 || Length[data] =!= rows,
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  If[rows > 0 && !And @@ (ListQ[#] && Length[#] === columns & /@ data),
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  If[!And @@ (cqExactRationalQ /@ Flatten[data]),
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  <|"rows" -> rows, "columns" -> columns, "data" -> data|>
];

cqRationalRawMatrixQ[matrix_] := Module[{rows, columns, data},
  If[!AssociationQ[matrix], Return[False]];
  rows = Lookup[matrix, "rows", None];
  columns = Lookup[matrix, "columns", None];
  data = Lookup[matrix, "data", None];
  IntegerQ[rows] && rows >= 0 && IntegerQ[columns] && columns >= 0 &&
    ListQ[data] && Length[data] === rows &&
    (rows === 0 || And @@ (ListQ[#] && Length[#] === columns & /@ data)) &&
    And @@ (cqExactRationalQ /@ Flatten[data])
];

cqRationalRawMatrixStructureQ[matrix_] := Module[{rows, columns, data},
  If[!AssociationQ[matrix], Return[False]];
  rows = Lookup[matrix, "rows", None];
  columns = Lookup[matrix, "columns", None];
  data = Lookup[matrix, "data", None];
  If[!IntegerQ[rows] || rows < 0 || !IntegerQ[columns] || columns < 0 ||
     !ListQ[data] || Length[data] =!= rows,
    Return[False]
  ];
  Which[
    rows === 0, data === {},
    columns === 0, data === ConstantArray[{}, rows],
    True, Dimensions[data] === {rows, columns}
  ]
];

cqRationalRawFromExact[matrix_, validate_: True] := Module[
  {rows, columns, values, data},
  If[TrueQ[validate] &&
      (!cqMatrixQ[matrix] || matrix["context"]["conductor"] =!= 1),
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  rows = matrix["rows"];
  columns = matrix["columns"];
  values = First[# ["coefficients"]] & /@ matrix["entries"];
  data = Which[
    rows === 0, {},
    columns === 0, ConstantArray[{}, rows],
    True, Partition[values, columns]
  ];
  cqRationalRawMatrix[rows, columns, data]
];

cqRationalElementDirect[context_, value_] := <|
  "type" -> "cyclotomic_element",
  "context_id" -> context["context_id"],
  "coefficients" -> {value}
|>;

cqRationalRawToExact[context_, matrix_Association] := Module[{values, entries},
  If[!cqContextQ[context] || context["conductor"] =!= 1,
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  values = Flatten[matrix["data"]];
  entries = cqRationalElementDirect[context, #] & /@ values;
  cqMatrixCreate[context, matrix["rows"], matrix["columns"], entries]
];

cqRationalRawIdentity[dimension_Integer] :=
  cqRationalRawMatrix[
    dimension,
    dimension,
    If[dimension === 0, {}, IdentityMatrix[dimension]]
  ];

cqRationalRawTranspose[matrix_Association] := Module[{data},
  data = Which[
    matrix["columns"] === 0, {},
    matrix["rows"] === 0, ConstantArray[{}, matrix["columns"]],
    True, Transpose[matrix["data"]]
  ];
  cqRationalRawMatrix[matrix["columns"], matrix["rows"], data]
];

cqRationalRawMultiply[left_Association, right_Association] := Module[
  {data, uniqueLeftRows, rowIndices, uniqueProduct},
  If[left["columns"] =!= right["rows"],
    Return[cqFailure["DimensionMismatch", <||>]]
  ];
  data = Which[
    left["rows"] === 0, {},
    right["columns"] === 0, ConstantArray[{}, left["rows"]],
    left["columns"] === 0, ConstantArray[0, {left["rows"], right["columns"]}],
    True,
      uniqueLeftRows = DeleteDuplicates[left["data"]];
      rowIndices = (
        First@FirstPosition[uniqueLeftRows, #] & /@ left["data"]
      );
      uniqueProduct = uniqueLeftRows . right["data"];
      uniqueProduct[[rowIndices]]
  ];
  cqRationalRawMatrix[left["rows"], right["columns"], data]
];

cqRationalRawZeroQ[matrix_Association] := MatchQ[
  matrix["data"],
  {RepeatedNull[{RepeatedNull[0]}]}
];

cqRationalRawRREF[matrix_Association, preparedRows_: False] := Module[
  {originalRows, rows, columns, data, pivotRow = 1, pivotColumns = {}, column,
   foundRow, candidateRow, pivot, factor, zeroRow, allColumns, freeColumns},
  originalRows = matrix["rows"];
  columns = matrix["columns"];
  data = If[
    TrueQ[preparedRows],
    matrix["data"],
    DeleteDuplicates@Select[
      matrix["data"],
      !And @@ (# === 0 & /@ #) &
    ]
  ];
  rows = Length[data];
  For[column = 1, column <= columns && pivotRow <= rows, column++,
    foundRow = 0;
    For[candidateRow = pivotRow, candidateRow <= rows && foundRow === 0, candidateRow++,
      If[data[[candidateRow, column]] =!= 0, foundRow = candidateRow]
    ];
    If[foundRow =!= 0,
      If[foundRow =!= pivotRow,
        data[[{pivotRow, foundRow}]] = data[[{foundRow, pivotRow}]]
      ];
      pivot = data[[pivotRow, column]];
      If[pivot =!= 1,
        data[[pivotRow, column ;; columns]] =
          data[[pivotRow, column ;; columns]]/pivot
      ];
      For[candidateRow = 1, candidateRow <= rows, candidateRow++,
        If[candidateRow =!= pivotRow,
          factor = data[[candidateRow, column]];
          If[factor =!= 0,
            data[[candidateRow, column ;; columns]] =
              data[[candidateRow, column ;; columns]] -
                factor data[[pivotRow, column ;; columns]];
            data[[candidateRow, column]] = 0
          ]
        ]
      ];
      AppendTo[pivotColumns, column - 1];
      pivotRow++
    ]
  ];
  zeroRow = ConstantArray[0, columns];
  data = Join[data, ConstantArray[zeroRow, originalRows - rows]];
  allColumns = If[columns === 0, {}, Range[0, columns - 1]];
  freeColumns = Complement[allColumns, pivotColumns];
  <|
    "reduced" -> cqRationalRawMatrix[originalRows, columns, data],
    "pivot_columns" -> pivotColumns,
    "free_columns" -> freeColumns,
    "rank" -> Length[pivotColumns],
    "nullity" -> Length[freeColumns]
  |>
];

cqRationalRawNullSpace[
  matrix_Association, verifyResidual_: True, preparedRows_: False
] := Module[
  {rref, columns, pivotColumns, freeColumns, rank, basisRows, vector,
   freeColumn, pivotIndex, pivotColumn, nullspaceRows, basisMatrix, residual},
  rref = cqRationalRawRREF[matrix, preparedRows];
  If[FailureQ[rref], Return[rref]];
  columns = matrix["columns"];
  pivotColumns = rref["pivot_columns"];
  freeColumns = rref["free_columns"];
  rank = rref["rank"];
  basisRows = Table[
    freeColumn = freeColumns[[freeIndex]];
    vector = ConstantArray[0, columns];
    vector[[freeColumn + 1]] = 1;
    For[pivotIndex = 1, pivotIndex <= rank, pivotIndex++,
      pivotColumn = pivotColumns[[pivotIndex]];
      vector[[pivotColumn + 1]] =
        -rref["reduced"]["data"][[pivotIndex, freeColumn + 1]]
    ];
    vector,
    {freeIndex, Length[freeColumns]}
  ];
  nullspaceRows = cqRationalRawMatrix[Length[freeColumns], columns, basisRows];
  If[FailureQ[nullspaceRows], Return[nullspaceRows]];
  basisMatrix = cqRationalRawTranspose[nullspaceRows];
  If[FailureQ[basisMatrix], Return[basisMatrix]];
  If[TrueQ[verifyResidual],
    residual = cqRationalRawMultiply[matrix, basisMatrix];
    If[FailureQ[residual] || !cqRationalRawZeroQ[residual],
      Return[cqFailure["ResidualVerificationFailed", <||>]]
    ]
  ];
  Join[rref, <|
    "nullspace_rows" -> nullspaceRows,
    "basis_matrix" -> basisMatrix,
    "exact_residual_verified" -> True
  |>]
];

cqRationalRREFResult[matrix_] := Module[{raw, rref, reduced},
  raw = cqRationalRawFromExact[matrix];
  If[FailureQ[raw], Return[raw]];
  rref = cqRationalRawRREF[raw];
  If[FailureQ[rref], Return[rref]];
  reduced = cqRationalRawToExact[matrix["context"], rref["reduced"]];
  If[FailureQ[reduced], Return[reduced]];
  <|
    "type" -> "rref_result",
    "reduced_matrix" -> reduced,
    "pivot_columns" -> rref["pivot_columns"],
    "free_columns" -> rref["free_columns"],
    "rank" -> rref["rank"],
    "nullity" -> rref["nullity"]
  |>
];

cqRationalNullSpaceResult[matrix_] := Module[
  {raw, kernel, reduced, rows, basis},
  raw = cqRationalRawFromExact[matrix];
  If[FailureQ[raw], Return[raw]];
  kernel = cqRationalRawNullSpace[raw];
  If[FailureQ[kernel], Return[kernel]];
  reduced = cqRationalRawToExact[matrix["context"], kernel["reduced"]];
  rows = cqRationalRawToExact[matrix["context"], kernel["nullspace_rows"]];
  basis = cqRationalRawToExact[matrix["context"], kernel["basis_matrix"]];
  If[AnyTrue[{reduced, rows, basis}, FailureQ],
    Return[FirstCase[{reduced, rows, basis}, _Failure]]
  ];
  <|
    "type" -> "kernel_result",
    "context" -> matrix["context"],
    "conductor" -> 1,
    "reduced_matrix" -> reduced,
    "nullspace_rows" -> rows,
    "basis_matrix" -> basis,
    "pivot_columns" -> kernel["pivot_columns"],
    "free_columns" -> kernel["free_columns"],
    "rank" -> kernel["rank"],
    "nullity" -> kernel["nullity"],
    "exact_residual_verified" -> True
  |>
];
