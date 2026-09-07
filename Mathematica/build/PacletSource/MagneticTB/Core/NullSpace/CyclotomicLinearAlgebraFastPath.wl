(* Coefficient-vector linear algebra over an exact polynomial quotient field.
   Public cyclotomic contexts and private real-subfield contexts share these
   arithmetic primitives, but only public contexts may escape the backend. *)

$cqMultiplicationTensorCache = <||>;
$cqCoefficientInverseCache = <||>;
$cqCoefficientFieldCache = <||>;
$cqCoefficientInverseCacheMaximumEntries = 4096;
$cqCoefficientBatchEliminationMinimumDegree = 12;
$cqCoefficientBatchEliminationMinimumDimension = 20;

cqCoefficientFieldContextQ[context_] :=
  AssociationQ[context] &&
    Lookup[context, "type", None] === "coefficient_field_context" &&
    MemberQ[
      {"full_cyclotomic_field", "maximal_real_subfield"},
      Lookup[context, "field_kind", None]
    ] &&
    StringQ[Lookup[context, "field_id", None]] &&
    IntegerQ[Lookup[context, "degree", None]] && context["degree"] > 0 &&
    ListQ[Lookup[context, "defining_polynomial", None]] &&
    Length[context["defining_polynomial"]] === context["degree"] + 1 &&
    Last[context["defining_polynomial"]] === 1 &&
    And @@ (cqExactRationalQ /@ context["defining_polynomial"]);

cqCoefficientFieldFromCyclotomic[context_] := Module[{key, cached, field},
  If[!cqContextQ[context], Return[cqFailure["ConductorMismatch", <||>]]];
  key = context["context_id"];
  cached = Lookup[$cqCoefficientFieldCache, key, Missing["not_cached"]];
  If[!MissingQ[cached], Return[cached]];
  field = <|
    "type" -> "coefficient_field_context",
    "field_kind" -> "full_cyclotomic_field",
    "field_id" -> StringJoin["coefficient:", key],
    "public_context_id" -> key,
    "degree" -> context["degree"],
    "defining_polynomial" -> context["cyclotomic_polynomial"]
  |>;
  If[!cqCoefficientFieldContextQ[field],
    Return[cqFailure["MalformedSerialization", <|"context" -> context|>]]
  ];
  $cqCoefficientFieldCache[key] = field;
  field
];

cqCoefficientFieldID[context_] := context["field_id"];

cqCoefficientFieldModulus[context_] := context["defining_polynomial"];

cqCoefficientZero[context_?cqCoefficientFieldContextQ] :=
  ConstantArray[0, context["degree"]];
cqCoefficientOne[context_?cqCoefficientFieldContextQ] :=
  ReplacePart[cqCoefficientZero[context], 1 -> 1];
cqCoefficientZeroQ[value_List] := MatchQ[value, {RepeatedNull[0]}];

cqCoefficientRawMatrix[context_, rows_Integer, columns_Integer, data_List] := Module[{},
  If[!cqCoefficientFieldContextQ[context] || rows < 0 || columns < 0 ||
     Length[data] =!= rows,
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  If[rows > 0 && !And @@ (ListQ[#] && Length[#] === columns & /@ data),
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  <|
    "context" -> context,
    "rows" -> rows,
    "columns" -> columns,
    "data" -> data
  |>
];

cqCoefficientRawMatrixQ[matrix_] := Module[
  {context, rows, columns, data, values},
  If[!AssociationQ[matrix], Return[False]];
  context = Lookup[matrix, "context", None];
  rows = Lookup[matrix, "rows", None];
  columns = Lookup[matrix, "columns", None];
  data = Lookup[matrix, "data", None];
  If[!cqCoefficientFieldContextQ[context] ||
     !IntegerQ[rows] || rows < 0 || !IntegerQ[columns] || columns < 0 ||
     !ListQ[data] || Length[data] =!= rows ||
     (rows > 0 && !And @@ (ListQ[#] && Length[#] === columns & /@ data)),
    Return[False]
  ];
  values = Flatten[data, 1];
  And @@ (
    ListQ[#] && Length[#] === context["degree"] &&
      And @@ (cqExactRationalQ /@ #) & /@ values
  )
];

cqCoefficientRawMatrixStructureQ[matrix_] := Module[
  {context, rows, columns, data},
  If[!AssociationQ[matrix], Return[False]];
  context = Lookup[matrix, "context", None];
  rows = Lookup[matrix, "rows", None];
  columns = Lookup[matrix, "columns", None];
  data = Lookup[matrix, "data", None];
  If[!cqCoefficientFieldContextQ[context] ||
     !IntegerQ[rows] || rows < 0 || !IntegerQ[columns] || columns < 0 ||
     !ListQ[data] || Length[data] =!= rows,
    Return[False]
  ];
  Which[
    rows === 0, data === {},
    columns === 0, data === ConstantArray[{}, rows],
    True, Dimensions[data] === {rows, columns, context["degree"]}
  ]
];

cqCoefficientRawFromExact[matrix_, validate_: True] := Module[
  {rows, columns, coefficients, data, fieldContext},
  If[TrueQ[validate] && !cqMatrixQ[matrix],
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  rows = matrix["rows"];
  columns = matrix["columns"];
  fieldContext = cqCoefficientFieldFromCyclotomic[matrix["context"]];
  If[FailureQ[fieldContext], Return[fieldContext]];
  coefficients = # ["coefficients"] & /@ matrix["entries"];
  data = Which[
    rows === 0, {},
    columns === 0, ConstantArray[{}, rows],
    True, Partition[coefficients, columns]
  ];
  cqCoefficientRawMatrix[fieldContext, rows, columns, data]
];

cqCoefficientElementDirect[context_, coefficients_List] := <|
  "type" -> "cyclotomic_element",
  "context_id" -> context["context_id"],
  "coefficients" -> coefficients
|>;

cqCoefficientRawToExact[publicContext_, matrix_Association] := Module[
  {fieldContext, entries},
  If[!cqContextQ[publicContext],
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  fieldContext = cqCoefficientFieldFromCyclotomic[publicContext];
  If[FailureQ[fieldContext], Return[fieldContext]];
  If[!cqCoefficientFieldContextQ[Lookup[matrix, "context", None]] ||
     cqCoefficientFieldID[matrix["context"]] =!= cqCoefficientFieldID[fieldContext],
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  entries = cqCoefficientElementDirect[publicContext, #] & /@
    Flatten[matrix["data"], 1];
  cqMatrixCreate[
    publicContext, matrix["rows"], matrix["columns"], entries
  ]
];

cqCoefficientMultiplicationTensor[context_?cqCoefficientFieldContextQ] := Module[
  {key, cached, degree, modulus, powers, current, exponent, padded, tensor},
  key = cqCoefficientFieldID[context];
  cached = Lookup[$cqMultiplicationTensorCache, key, Missing["not_cached"]];
  If[!MissingQ[cached], Return[cached]];
  degree = context["degree"];
  modulus = cqCoefficientFieldModulus[context];
  current = {1};
  powers = {};
  For[exponent = 0, exponent <= 2 degree - 2, exponent++,
    padded = cqPolynomialPad[current, degree];
    If[FailureQ[padded], Return[padded]];
    AppendTo[powers, padded];
    current = cqPolynomialMod[Prepend[current, 0], modulus];
    If[FailureQ[current], Return[current]]
  ];
  (* Tensor order {left basis, result basis, right basis} makes a.tensor.b
     return the coefficient vector of the product. *)
  tensor = Table[
    powers[[leftIndex + rightIndex - 1, resultIndex]],
    {leftIndex, degree}, {resultIndex, degree}, {rightIndex, degree}
  ];
  $cqMultiplicationTensorCache[key] = tensor;
  tensor
];

cqCoefficientMultiply[context_?cqContextQ, left_List, right_List] := Module[
  {fieldContext},
  fieldContext = cqCoefficientFieldFromCyclotomic[context];
  If[FailureQ[fieldContext], Return[fieldContext]];
  cqCoefficientMultiply[fieldContext, left, right]
];

cqCoefficientMultiply[
  context_?cqCoefficientFieldContextQ, left_List, right_List
] := Module[{tensor},
  If[cqCoefficientZeroQ[left] || cqCoefficientZeroQ[right],
    Return[cqCoefficientZero[context]]
  ];
  tensor = cqCoefficientMultiplicationTensor[context];
  If[FailureQ[tensor], Return[tensor]];
  left . tensor . right
];

cqCoefficientInverse[
  context_?cqCoefficientFieldContextQ, value_List
] := Module[
  {key, cached, extended, padded},
  If[cqCoefficientZeroQ[value], Return[cqFailure["DivisionByZero", <||>]]];
  If[And @@ (# === 0 & /@ Rest[value]),
    Return[ReplacePart[cqCoefficientZero[context], 1 -> 1/First[value]]]
  ];
  key = {cqCoefficientFieldID[context], value};
  cached = Lookup[$cqCoefficientInverseCache, Key[key], Missing["not_cached"]];
  If[!MissingQ[cached], Return[cached]];
  extended = cqPolynomialExtendedGCD[value, cqCoefficientFieldModulus[context]];
  If[FailureQ[extended], Return[extended]];
  If[extended["gcd"] =!= {1},
    Return[cqFailure["SingularPolynomialElement", <||>]]
  ];
  padded = cqPolynomialPad[extended["left_coefficient"], context["degree"]];
  If[FailureQ[padded], Return[padded]];
  If[Length[$cqCoefficientInverseCache] < $cqCoefficientInverseCacheMaximumEntries,
    AssociateTo[$cqCoefficientInverseCache, key -> padded]
  ];
  padded
];

cqCoefficientRawIdentity[context_, dimension_Integer] := Module[{zero, one, data},
  zero = cqCoefficientZero[context];
  one = cqCoefficientOne[context];
  data = If[
    dimension === 0,
    {},
    Table[If[row === column, one, zero], {row, dimension}, {column, dimension}]
  ];
  cqCoefficientRawMatrix[context, dimension, dimension, data]
];

cqCoefficientRawTranspose[matrix_Association] := Module[{data},
  data = Which[
    matrix["columns"] === 0, {},
    matrix["rows"] === 0, ConstantArray[{}, matrix["columns"]],
    True, Transpose[matrix["data"]]
  ];
  cqCoefficientRawMatrix[
    matrix["context"], matrix["columns"], matrix["rows"], data
  ]
];

cqCoefficientRawMultiply[left_Association, right_Association] := Module[
  {context, zero, tensor, data, row, inner, leftValue, leftOperator,
   rightNonzeroColumns, activeColumns, contributions, resultRow,
   uniqueLeftRows, rowIndices, rowIndexMap, uniqueData},
  If[cqCoefficientFieldID[left["context"]] =!=
       cqCoefficientFieldID[right["context"]],
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  If[left["columns"] =!= right["rows"],
    Return[cqFailure["DimensionMismatch", <||>]]
  ];
  context = left["context"];
  zero = cqCoefficientZero[context];
  If[left["rows"] === 0,
    Return[cqCoefficientRawMatrix[context, 0, right["columns"], {}]]
  ];
  If[right["columns"] === 0,
    Return[cqCoefficientRawMatrix[
      context, left["rows"], 0, ConstantArray[{}, left["rows"]]
    ]]
  ];
  If[left["columns"] === 0,
    Return[cqCoefficientRawMatrix[
      context,
      left["rows"],
      right["columns"],
      ConstantArray[zero, {left["rows"], right["columns"]}]
    ]]
  ];
  tensor = cqCoefficientMultiplicationTensor[context];
  If[FailureQ[tensor], Return[tensor]];
  uniqueLeftRows = DeleteDuplicates[left["data"]];
  rowIndexMap = Association@MapIndexed[Rule[#1, First[#2]] &, uniqueLeftRows];
  rowIndices = rowIndexMap[#] & /@ left["data"];
  uniqueData = Table[zero, {Length[uniqueLeftRows]}, {right["columns"]}];
  rightNonzeroColumns = (
    Pick[
      Range[right["columns"]],
      cqCoefficientZeroQ /@ #,
      False
    ] & /@ right["data"]
  );
  For[row = 1, row <= Length[uniqueLeftRows], row++,
    resultRow = uniqueData[[row]];
    For[inner = 1, inner <= left["columns"], inner++,
      leftValue = uniqueLeftRows[[row, inner]];
      If[!cqCoefficientZeroQ[leftValue],
        leftOperator = leftValue . tensor;
        activeColumns = rightNonzeroColumns[[inner]];
        If[activeColumns =!= {},
          contributions =
            right["data"][[inner, activeColumns]] . Transpose[leftOperator];
          resultRow[[activeColumns]] =
            resultRow[[activeColumns]] + contributions
        ]
      ]
    ];
    uniqueData[[row]] = resultRow
  ];
  data = uniqueData[[rowIndices]];
  cqCoefficientRawMatrix[context, left["rows"], right["columns"], data]
];

cqCoefficientRawZeroQ[matrix_Association] := MatchQ[
  matrix["data"],
  {RepeatedNull[{RepeatedNull[{RepeatedNull[0]}]}]}
];

cqCoefficientRawRREF[matrix_Association, preparedRows_: False] := Module[
  {context, originalRows, rows, columns, data, pivotRow = 1, pivotColumns = {}, column,
   foundRow, candidateRow, pivot, inverse, inverseOperator, factor,
   factorOperator, tensor, zero, one, zeroRow, allColumns, freeColumns,
   pivotSlice, activeRows, factorOperators, corrections, batchEliminationQ},
  context = matrix["context"];
  originalRows = matrix["rows"];
  columns = matrix["columns"];
  data = If[
    TrueQ[preparedRows],
    matrix["data"],
    DeleteDuplicates@Select[
      matrix["data"],
      !And @@ (cqCoefficientZeroQ /@ #) &
    ]
  ];
  rows = Length[data];
  tensor = cqCoefficientMultiplicationTensor[context];
  If[FailureQ[tensor], Return[tensor]];
  zero = cqCoefficientZero[context];
  one = cqCoefficientOne[context];
  batchEliminationQ =
    context["degree"] >= $cqCoefficientBatchEliminationMinimumDegree &&
    rows >= $cqCoefficientBatchEliminationMinimumDimension &&
    columns >= $cqCoefficientBatchEliminationMinimumDimension;
  For[column = 1, column <= columns && pivotRow <= rows, column++,
    foundRow = 0;
    For[candidateRow = pivotRow, candidateRow <= rows && foundRow === 0, candidateRow++,
      If[!cqCoefficientZeroQ[data[[candidateRow, column]]], foundRow = candidateRow]
    ];
    If[foundRow =!= 0,
      If[foundRow =!= pivotRow,
        data[[{pivotRow, foundRow}]] = data[[{foundRow, pivotRow}]]
      ];
      pivot = data[[pivotRow, column]];
      inverse = cqCoefficientInverse[context, pivot];
      If[FailureQ[inverse], Return[inverse]];
      inverseOperator = inverse . tensor;
      pivotSlice =
        data[[pivotRow, column ;; columns]] . Transpose[inverseOperator];
      pivotSlice[[1]] = one;
      data[[pivotRow, column ;; columns]] = pivotSlice;
      If[batchEliminationQ,
        (* Batch all row eliminations for this pivot.  Both the complete
           cyclotomic field and real-subfield contexts use this same path. *)
        activeRows = Select[
          Range[rows],
          # =!= pivotRow && !cqCoefficientZeroQ[data[[#, column]]] &
        ];
        If[activeRows =!= {},
          factorOperators = data[[activeRows, column]] . tensor;
          corrections =
            (pivotSlice . Transpose[#] &) /@ factorOperators;
          data[[activeRows, column ;; columns]] =
            data[[activeRows, column ;; columns]] - corrections;
          data[[activeRows, column]] = ConstantArray[
            zero, Length[activeRows]
          ]
        ]
        ,
        For[candidateRow = 1, candidateRow <= rows, candidateRow++,
          If[candidateRow =!= pivotRow,
            factor = data[[candidateRow, column]];
            If[!cqCoefficientZeroQ[factor],
              factorOperator = factor . tensor;
              data[[candidateRow, column ;; columns]] =
                data[[candidateRow, column ;; columns]] -
                  pivotSlice . Transpose[factorOperator];
              data[[candidateRow, column]] = zero
            ]
          ]
        ]
      ];
      AppendTo[pivotColumns, column - 1];
      pivotRow++
    ]
  ];
  zeroRow = ConstantArray[zero, columns];
  data = Join[data, ConstantArray[zeroRow, originalRows - rows]];
  allColumns = If[columns === 0, {}, Range[0, columns - 1]];
  freeColumns = Complement[allColumns, pivotColumns];
  <|
    "reduced" -> cqCoefficientRawMatrix[context, originalRows, columns, data],
    "pivot_columns" -> pivotColumns,
    "free_columns" -> freeColumns,
    "rank" -> Length[pivotColumns],
    "nullity" -> Length[freeColumns]
  |>
];

cqCoefficientRawNullSpace[
  matrix_Association, verifyResidual_: True, preparedRows_: False
] := Module[
  {context, rref, columns, pivotColumns, freeColumns, rank, zero, one,
   basisRows, vector, freeColumn, pivotIndex, pivotColumn, nullspaceRows,
   basisMatrix, residual},
  context = matrix["context"];
  rref = cqCoefficientRawRREF[matrix, preparedRows];
  If[FailureQ[rref], Return[rref]];
  columns = matrix["columns"];
  pivotColumns = rref["pivot_columns"];
  freeColumns = rref["free_columns"];
  rank = rref["rank"];
  zero = cqCoefficientZero[context];
  one = cqCoefficientOne[context];
  basisRows = Table[
    freeColumn = freeColumns[[freeIndex]];
    vector = ConstantArray[zero, columns];
    vector[[freeColumn + 1]] = one;
    For[pivotIndex = 1, pivotIndex <= rank, pivotIndex++,
      pivotColumn = pivotColumns[[pivotIndex]];
      vector[[pivotColumn + 1]] =
        -rref["reduced"]["data"][[pivotIndex, freeColumn + 1]]
    ];
    vector,
    {freeIndex, Length[freeColumns]}
  ];
  nullspaceRows = cqCoefficientRawMatrix[
    context, Length[freeColumns], columns, basisRows
  ];
  If[FailureQ[nullspaceRows], Return[nullspaceRows]];
  basisMatrix = cqCoefficientRawTranspose[nullspaceRows];
  If[FailureQ[basisMatrix], Return[basisMatrix]];
  If[TrueQ[verifyResidual],
    residual = cqCoefficientRawMultiply[matrix, basisMatrix];
    If[FailureQ[residual] || !cqCoefficientRawZeroQ[residual],
      Return[cqFailure["ResidualVerificationFailed", <||>]]
    ]
  ];
  Join[rref, <|
    "nullspace_rows" -> nullspaceRows,
    "basis_matrix" -> basisMatrix,
    "exact_residual_verified" -> True
  |>]
];

cqCoefficientRREFResult[matrix_] := Module[{raw, rref, reduced},
  raw = cqCoefficientRawFromExact[matrix];
  If[FailureQ[raw], Return[raw]];
  rref = cqCoefficientRawRREF[raw];
  If[FailureQ[rref], Return[rref]];
  reduced = cqCoefficientRawToExact[matrix["context"], rref["reduced"]];
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

cqCoefficientNullSpaceResult[matrix_] := Module[
  {raw, kernel, reduced, rows, basis},
  raw = cqCoefficientRawFromExact[matrix];
  If[FailureQ[raw], Return[raw]];
  kernel = cqCoefficientRawNullSpace[raw];
  If[FailureQ[kernel], Return[kernel]];
  reduced = cqCoefficientRawToExact[matrix["context"], kernel["reduced"]];
  rows = cqCoefficientRawToExact[matrix["context"], kernel["nullspace_rows"]];
  basis = cqCoefficientRawToExact[matrix["context"], kernel["basis_matrix"]];
  If[AnyTrue[{reduced, rows, basis}, FailureQ],
    Return[FirstCase[{reduced, rows, basis}, _Failure]]
  ];
  <|
    "type" -> "kernel_result",
    "context" -> matrix["context"],
    "conductor" -> matrix["context"]["conductor"],
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
