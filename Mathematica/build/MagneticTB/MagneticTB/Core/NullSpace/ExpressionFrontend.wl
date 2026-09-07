cqFrontendRootSum[constant_, terms_List] := <|
  "kind" -> "root_sum",
  "constant" -> constant,
  "terms" -> terms
|>;

cqFrontendConstant[value_] := cqFrontendRootSum[value, {}];

cqFrontendRoot[order_Integer, power_Integer, coefficient_: 1] := Module[{},
  If[!cqExactRationalQ[coefficient],
    Return[cqFailure["NonExactInput", <|"value" -> HoldForm[coefficient]|>]]
  ];
  If[order < 1,
    Return[cqFailure["InvalidRootOfUnity", <|"order" -> order, "power" -> power|>]]
  ];
  If[coefficient === 0, Return[cqFrontendConstant[0]]];
  If[order === 1,
    Return[cqFrontendConstant[coefficient]]
  ];
  cqFrontendRootSum[0, {
    <|"coefficient" -> coefficient, "order" -> order, "power" -> Mod[power, order]|>
  }]
];

cqFrontendAdd[left_Association, right_Association] :=
  cqFrontendRootSum[
    left["constant"] + right["constant"],
    Join[left["terms"], right["terms"]]
  ];

cqFrontendNegate[value_Association] := cqFrontendRootSum[
  -value["constant"],
  Map[
    <|
      "coefficient" -> -# ["coefficient"],
      "order" -> # ["order"],
      "power" -> # ["power"]
    |> &,
    value["terms"]
  ]
];

cqFrontendConjugate[value_Association] := cqFrontendRootSum[
  value["constant"],
  Map[
    <|
      "coefficient" -> # ["coefficient"],
      "order" -> # ["order"],
      "power" -> Mod[-# ["power"], # ["order"]]
    |> &,
    value["terms"]
  ]
];

cqFrontendComponents[value_Association] := Module[{components = value["terms"]},
  If[value["constant"] =!= 0,
    components = Prepend[components, <|
      "coefficient" -> value["constant"],
      "order" -> 1,
      "power" -> 0
    |>]
  ];
  components
];

cqFrontendMultiply[left_Association, right_Association] := Module[
  {leftComponents, rightComponents, result = cqFrontendConstant[0], leftTerm,
   rightTerm, coefficient, order, power, leftIndex, rightIndex},
  leftComponents = cqFrontendComponents[left];
  rightComponents = cqFrontendComponents[right];
  For[leftIndex = 1, leftIndex <= Length[leftComponents], leftIndex++,
    leftTerm = leftComponents[[leftIndex]];
    For[rightIndex = 1, rightIndex <= Length[rightComponents], rightIndex++,
      rightTerm = rightComponents[[rightIndex]];
      coefficient = leftTerm["coefficient"] rightTerm["coefficient"];
      If[coefficient =!= 0,
        order = cqIntegerLCM[leftTerm["order"], rightTerm["order"]];
        power = Mod[
          leftTerm["power"] Quotient[order, leftTerm["order"]] +
          rightTerm["power"] Quotient[order, rightTerm["order"]],
          order
        ];
        If[order === 1,
          result["constant"] = result["constant"] + coefficient,
          AppendTo[result["terms"], <|
            "coefficient" -> coefficient,
            "order" -> order,
            "power" -> power
          |>]
        ];
      ];
    ];
  ];
  result
];

cqFrontendPower[value_Association, exponent_Integer] := Module[
  {power = exponent, factor = value, result = cqFrontendConstant[1]},
  If[exponent < 0,
    Return[cqFailure[
      "UnsupportedScalarInput",
      <|"reason" -> "negative_power_of_non_monomial_expression"|>
    ]]
  ];
  While[power > 0,
    If[OddQ[power], result = cqFrontendMultiply[result, factor]];
    power = Quotient[power, 2];
    If[power > 0, factor = cqFrontendMultiply[factor, factor]];
  ];
  result
];

cqFrontendSquareFreeDecomposition[value_Integer] := Module[
  {factors, outside = 1, squareFree = 1, factor, index},
  If[value < 1,
    Return[cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[value]|>]]
  ];
  factors = cqIntegerFactorization[value];
  If[FailureQ[factors], Return[factors]];
  For[index = 1, index <= Length[factors], index++,
    factor = factors[[index]];
    outside *= factor[[1]]^Quotient[factor[[2]], 2];
    If[OddQ[factor[[2]]], squareFree *= factor[[1]]];
  ];
  {outside, squareFree}
];

cqFrontendSqrtPositiveRational[value_] := Module[
  {combined, decomposition, outside, squareFree, discriminant, scale, terms = {},
   symbol, power},
  If[!cqExactRationalQ[value] || value <= 0,
    Return[cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[value]|>]]
  ];
  combined = Numerator[value] Denominator[value];
  decomposition = cqFrontendSquareFreeDecomposition[combined];
  If[FailureQ[decomposition], Return[decomposition]];
  outside = decomposition[[1]];
  squareFree = decomposition[[2]];
  scale = outside/Denominator[value];
  If[squareFree === 1, Return[cqFrontendConstant[scale]]];
  discriminant = If[Mod[squareFree, 4] === 1, squareFree, 4 squareFree];
  If[discriminant =!= squareFree, scale /= 2];
  For[power = 0, power < discriminant, power++,
    symbol = KroneckerSymbol[discriminant, power];
    If[symbol =!= 0,
      AppendTo[terms, <|
        "coefficient" -> scale symbol,
        "order" -> discriminant,
        "power" -> power
      |>]
    ];
  ];
  cqFrontendRootSum[0, terms]
];

cqFrontendSqrtRational[value_] := Module[{positive, imaginaryUnit},
  If[!cqExactRationalQ[value],
    Return[cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[value]|>]]
  ];
  If[value === 0, Return[cqFrontendConstant[0]]];
  If[value > 0, Return[cqFrontendSqrtPositiveRational[value]]];
  positive = cqFrontendSqrtPositiveRational[-value];
  If[FailureQ[positive], Return[positive]];
  imaginaryUnit = cqFrontendRoot[4, 1];
  cqFrontendMultiply[imaginaryUnit, positive]
];

cqFrontendHalfIntegerPower[base_, exponent_] := Module[{squareRoot, rationalFactor},
  If[!cqExactRationalQ[base] || Denominator[exponent] =!= 2,
    Return[cqFailure[
      "UnsupportedScalarInput",
      <|"base" -> HoldForm[base], "exponent" -> HoldForm[exponent]|>
    ]]
  ];
  squareRoot = cqFrontendSqrtRational[base];
  If[FailureQ[squareRoot], Return[squareRoot]];
  rationalFactor = base^Quotient[Numerator[exponent] - 1, 2];
  cqFrontendMultiply[cqFrontendConstant[rationalFactor], squareRoot]
];

cqFrontendRationalAngle[angle_, fullTurn_] := Module[{ratio},
  ratio = angle/fullTurn;
  If[cqExactRationalQ[ratio], ratio, Missing["not_rational_angle"]]
];

cqExactExpressionToRootSum[expression_] := Module[
  {head, arguments, result, converted, index, base, exponent, ratio, root, conjugate,
   difference, imaginaryScale, parsed},
  If[cqExactRationalQ[expression], Return[cqFrontendConstant[expression]]];
  If[AssociationQ[expression],
    parsed = cqParseRootSum[expression];
    If[FailureQ[parsed], Return[parsed]];
    Return[expression]
  ];
  If[Head[expression] === Real,
    Return[cqFailure["NonExactInput", <|"value" -> HoldForm[expression]|>]]
  ];
  If[expression === I, Return[cqFrontendRoot[4, 1]]];
  head = Head[expression];
  arguments = List @@ expression;
  If[head === Plus,
    result = cqFrontendConstant[0];
    For[index = 1, index <= Length[arguments], index++,
      converted = cqExactExpressionToRootSum[arguments[[index]]];
      If[FailureQ[converted], Return[converted]];
      result = cqFrontendAdd[result, converted];
    ];
    Return[result]
  ];
  If[head === Times,
    result = cqFrontendConstant[1];
    For[index = 1, index <= Length[arguments], index++,
      converted = cqExactExpressionToRootSum[arguments[[index]]];
      If[FailureQ[converted], Return[converted]];
      result = cqFrontendMultiply[result, converted];
    ];
    Return[result]
  ];
  If[head === Power,
    base = arguments[[1]];
    exponent = arguments[[2]];
    If[base === E,
      ratio = cqFrontendRationalAngle[exponent, 2 Pi I];
      If[MissingQ[ratio],
        Return[cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[expression]|>]]
      ];
      Return[cqFrontendRoot[Denominator[ratio], Numerator[ratio]]]
    ];
    If[base === -1 && cqExactRationalQ[exponent],
      ratio = exponent/2;
      Return[cqFrontendRoot[Denominator[ratio], Numerator[ratio]]]
    ];
    If[cqExactRationalQ[base] && cqExactRationalQ[exponent] && Denominator[exponent] === 2,
      Return[cqFrontendHalfIntegerPower[base, exponent]]
    ];
    If[IntegerQ[exponent],
      converted = cqExactExpressionToRootSum[base];
      If[FailureQ[converted], Return[converted]];
      Return[cqFrontendPower[converted, exponent]]
    ];
    Return[cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[expression]|>]]
  ];
  If[head === Cos || head === Sin,
    ratio = cqFrontendRationalAngle[arguments[[1]], 2 Pi];
    If[MissingQ[ratio],
      Return[cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[expression]|>]]
    ];
    root = cqFrontendRoot[Denominator[ratio], Numerator[ratio]];
    conjugate = cqFrontendConjugate[root];
    If[head === Cos,
      Return[cqFrontendMultiply[
        cqFrontendConstant[1/2],
        cqFrontendAdd[root, conjugate]
      ]]
    ];
    difference = cqFrontendAdd[root, cqFrontendNegate[conjugate]];
    imaginaryScale = cqFrontendRoot[4, 1, -1/2];
    Return[cqFrontendMultiply[imaginaryScale, difference]]
  ];
  If[head === Conjugate,
    converted = cqExactExpressionToRootSum[First[arguments]];
    If[FailureQ[converted], Return[converted]];
    Return[cqFrontendConjugate[converted]]
  ];
  cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[expression]|>]
];

cqInferDirectCoordinateDimension[matrices_List, requested_] := Module[
  {candidates = {}, matrix, rows, rowLengths, columns, matrixIndex},
  For[matrixIndex = 1, matrixIndex <= Length[matrices], matrixIndex++,
    matrix = matrices[[matrixIndex]];
    If[AssociationQ[matrix],
      columns = Lookup[matrix, "columns", Missing["columns"]];
      If[!IntegerQ[columns] || columns < 0,
        Return[cqFailure["MalformedSerialization", <|"matrix" -> matrix|>]]
      ];
      AppendTo[candidates, columns],
      If[!ListQ[matrix],
        Return[cqFailure["MalformedSerialization", <|"matrix" -> HoldForm[matrix]|>]]
      ];
      rows = Length[matrix];
      If[rows > 0,
        If[!And @@ (ListQ /@ matrix),
          Return[cqFailure["MalformedSerialization", <|"matrix" -> HoldForm[matrix]|>]]
        ];
        rowLengths = Length /@ matrix;
        If[Length[DeleteDuplicates[rowLengths]] =!= 1,
          Return[cqFailure["DimensionMismatch", <|"row_lengths" -> rowLengths|>]]
        ];
        AppendTo[candidates, First[rowLengths]];
      ];
    ];
  ];
  If[requested === Automatic,
    If[candidates === {},
      Return[cqFailure["CoordinateDimensionRequired", <||>]]
    ];
    columns = First[candidates],
    If[!IntegerQ[requested] || requested < 0,
      Return[cqFailure["DimensionMismatch", <|"coordinate_dimension" -> requested|>]]
    ];
    columns = requested
  ];
  If[!And @@ (# === columns & /@ candidates),
    Return[cqFailure[
      "DimensionMismatch",
      <|"coordinate_dimension" -> columns, "matrix_columns" -> candidates|>
    ]]
  ];
  columns
];

cqCompileDirectMatrices[matrices_List, coordinateDimension_, maxDegree_] := Module[
  {dimension, normalizedMatrices, rationalInputQ, context, compiledMatrices,
   convertedMatrices, convertedEntryLists, cachedConvert},
  dimension = cqInferDirectCoordinateDimension[matrices, coordinateDimension];
  If[FailureQ[dimension], Return[dimension]];
  normalizedMatrices = cqNormalizeRawMatrix[#, dimension] & /@ matrices;
  If[AnyTrue[normalizedMatrices, FailureQ],
    Return[FirstCase[normalizedMatrices, _Failure]]
  ];
  rationalInputQ = And @@ (
    And @@ (cqExactRationalQ /@ #["entries"]) & /@ normalizedMatrices
  );
  If[rationalInputQ,
    context = cqCreateContext[1, maxDegree];
    If[FailureQ[context], Return[context]];
    compiledMatrices = (
      cqMatrixCreate[
        context,
        #["rows"],
        #["columns"],
        (cqRationalElementDirect[context, #] & /@ #["entries"])
      ] & /@ normalizedMatrices
    );
    If[AnyTrue[compiledMatrices, FailureQ],
      Return[FirstCase[compiledMatrices, _Failure]]
    ];
    Return[<|
      "type" -> "compiled_cyclotomic_problem",
      "context" -> context,
      "conductor" -> 1,
      "coordinate_dimension" -> dimension,
      "constraint_count" -> Length[compiledMatrices],
      "constraints" -> compiledMatrices
    |>]
  ];
  cachedConvert[value_] := cachedConvert[value] = cqExactExpressionToRootSum[value];
  convertedEntryLists = (cachedConvert /@ # ["entries"]) & /@ normalizedMatrices;
  If[AnyTrue[Flatten[convertedEntryLists, 1], FailureQ],
    Return[FirstCase[Flatten[convertedEntryLists, 1], _Failure]]
  ];
  convertedMatrices = MapThread[
    <|
      "rows" -> #1["rows"],
      "columns" -> #1["columns"],
      "entries" -> #2
    |> &,
    {normalizedMatrices, convertedEntryLists}
  ];
  CompileCyclotomicMatrices[<|
    "coordinate_dimension" -> dimension,
    "constraints" -> convertedMatrices,
    "max_cyclotomic_degree" -> maxDegree
  |>]
];

Options[CyclotomicCommonNullSpace] = {
  "CoordinateDimension" -> Automatic,
  "MaxCyclotomicDegree" -> 128
};

$cqFrontendRadicalDegreeLimit = 32;

cqFrontendPrettyElement[context_, element_] := Module[
  {expression, candidates, conjugate, radicalCandidate},
  expression = cqRestoreElement[context, element];
  If[FailureQ[expression], Return[expression]];
  candidates = {expression, ExpToTrig[expression]};
  conjugate = cqElementConjugate[context, element];
  If[
    !FailureQ[conjugate] &&
    cqElementEqualQ[context, element, conjugate] &&
    context["degree"] <= $cqFrontendRadicalDegreeLimit,
    radicalCandidate = Quiet[Check[ToRadicals[RootReduce[expression]], expression]];
    If[FreeQ[radicalCandidate, _Root], AppendTo[candidates, radicalCandidate]];
  ];
  First[MinimalBy[DeleteDuplicates[candidates], LeafCount]]
];

cqFrontendPrettyMatrix[matrix_] := Module[{entries, context, pretty},
  If[!cqMatrixQ[matrix], Return[cqFailure["MalformedSerialization", <||>]]];
  context = matrix["context"];
  pretty[element_Association] := pretty[element] =
    cqFrontendPrettyElement[context, element];
  entries = pretty /@ matrix["entries"];
  If[AnyTrue[entries, FailureQ], Return[FirstCase[entries, _Failure]]];
  If[matrix["rows"] === 0, Return[{}]];
  If[matrix["columns"] === 0, Return[ConstantArray[{}, matrix["rows"]]]];
  Partition[entries, matrix["columns"]]
];

CyclotomicCommonNullSpace[matrices_List, OptionsPattern[]] := Module[
  {compiled, result, restoredRows},
  compiled = cqCompileDirectMatrices[
    matrices,
    OptionValue["CoordinateDimension"],
    OptionValue["MaxCyclotomicDegree"]
  ];
  If[FailureQ[compiled], Return[compiled]];
  result = CyclotomicCommonKernel[compiled];
  If[FailureQ[result], Return[result]];
  restoredRows = cqFrontendPrettyMatrix[result["nullspace_rows"]];
  If[FailureQ[restoredRows], Return[restoredRows]];
  restoredRows
];

CyclotomicCommonNullSpace[input_, OptionsPattern[]] :=
  cqFailure["MalformedSerialization", <|"input" -> HoldForm[input]|>];
