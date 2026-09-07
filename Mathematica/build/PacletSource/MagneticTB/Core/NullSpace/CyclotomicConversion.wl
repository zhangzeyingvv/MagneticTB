cqParseRootSum[value_] /; cqExactRationalQ[value] :=
  <|"constant" -> value, "terms" -> {}, "orders" -> {}|>;

cqParseRootSum[value_Real] := cqFailure["NonExactInput", <|"value" -> HoldForm[value]|>];

cqParseRootSum[data_Association] := Module[
  {kind, constantRaw, termsRaw, constant, accumulated = <||>, metadata = <||>, orders = {},
   coefficient, order, power, normalizedPower, key, normalizedTerms = {}, term,
   termIndex},
  kind = Lookup[data, "kind", Missing["kind"]];
  If[kind =!= "root_sum",
    Return[cqFailure["UnsupportedScalarInput", <|"value" -> data|>]]
  ];
  constantRaw = Lookup[data, "constant", Missing["constant"]];
  termsRaw = Lookup[data, "terms", Missing["terms"]];
  If[MissingQ[constantRaw] || MissingQ[termsRaw] || !ListQ[termsRaw],
    Return[cqFailure["MalformedSerialization", <|"value" -> data|>]]
  ];
  constant = cqParseRational[constantRaw];
  If[FailureQ[constant], Return[constant]];
  For[termIndex = 1, termIndex <= Length[termsRaw], termIndex++,
    term = termsRaw[[termIndex]];
    If[!AssociationQ[term],
      Return[cqFailure["MalformedSerialization", <|"term" -> HoldForm[term]|>]]
    ];
    coefficient = cqParseRational[Lookup[term, "coefficient", Missing["coefficient"]]];
    If[FailureQ[coefficient], Return[coefficient]];
    order = Lookup[term, "order", Missing["order"]];
    power = Lookup[term, "power", Missing["power"]];
    If[!IntegerQ[order] || order < 1 || !IntegerQ[power],
      Return[cqFailure[
        "InvalidRootOfUnity",
        <|"order" -> order, "power" -> power|>
      ]]
    ];
    If[coefficient =!= 0,
      AppendTo[orders, order];
      normalizedPower = Mod[power, order];
      key = StringJoin[ToString[order], ":", ToString[normalizedPower]];
      accumulated[key] = Lookup[accumulated, key, 0] + coefficient;
      metadata[key] = {order, normalizedPower};
    ];
  ];
  Do[
    If[accumulated[key] =!= 0,
      AppendTo[normalizedTerms, <|
        "coefficient" -> accumulated[key],
        "order" -> metadata[key][[1]],
        "power" -> metadata[key][[2]]
      |>]
    ];
  , {key, Keys[accumulated]}];
  <|"constant" -> constant, "terms" -> normalizedTerms, "orders" -> orders|>
];

cqParseRootSum[value_] :=
  cqFailure["UnsupportedScalarInput", <|"value" -> HoldForm[value]|>];

cqCompileParsedRootSum[parsed_Association, context_] := Module[
  {result, root, scaled, term, termIndex},
  result = cqElementFromPolynomial[context, {parsed["constant"]}];
  If[FailureQ[result], Return[result]];
  For[termIndex = 1, termIndex <= Length[parsed["terms"]], termIndex++,
    term = parsed["terms"][[termIndex]];
    root = cqRootOfUnityElement[context, term["order"], term["power"]];
    If[FailureQ[root], Return[root]];
    scaled = cqElementScale[context, term["coefficient"], root];
    If[FailureQ[scaled], Return[scaled]];
    result = cqElementAdd[context, result, scaled];
    If[FailureQ[result], Return[result]];
  ];
  result
];

cqRestoreElement[context_, element_] := Module[{terms, index},
  If[!cqElementQ[context, element], Return[cqFailure["ConductorMismatch", <||>]]];
  terms = Table[
    element["coefficients"][[index + 1]] Exp[2 Pi I index/context["conductor"]],
    {index, 0, context["degree"] - 1}
  ];
  Total[terms]
];

RestoreCyclotomicExpression[context_Association, element_Association] :=
  cqRestoreElement[context, element];

RestoreCyclotomicExpression[object_Association] := Module[{type, context, entries, restored},
  type = Lookup[object, "type", None];
  Switch[type,
    "exact_matrix",
      context = object["context"];
      entries = cqRestoreElement[context, #] & /@ object["entries"];
      If[AnyTrue[entries, FailureQ], Return[FirstCase[entries, _Failure]]];
      <|
        "type" -> "restored_matrix",
        "rows" -> object["rows"],
        "columns" -> object["columns"],
        "entries" -> entries
      |>,
    "kernel_result",
      restored = object;
      restored["reduced_matrix"] = RestoreCyclotomicExpression[object["reduced_matrix"]];
      restored["nullspace_rows"] = RestoreCyclotomicExpression[object["nullspace_rows"]];
      restored["basis_matrix"] = RestoreCyclotomicExpression[object["basis_matrix"]];
      restored,
    "common_kernel_result",
      restored = object;
      restored["basis_matrix"] = RestoreCyclotomicExpression[object["basis_matrix"]];
      restored["nullspace_rows"] = RestoreCyclotomicExpression[object["nullspace_rows"]];
      restored,
    _, cqFailure["MalformedSerialization", <|"object" -> object|>]
  ]
];

RestoreCyclotomicExpression[object_] :=
  cqFailure["MalformedSerialization", <|"object" -> HoldForm[object]|>];
