$cqContextCache = <||>;

cqContextQ[context_] := Module[{conductor, degree, polynomial},
  If[!AssociationQ[context] ||
     Lookup[context, "type", None] =!= "cyclotomic_context" ||
     !StringQ[Lookup[context, "context_id", None]],
    Return[False]
  ];
  conductor = Lookup[context, "conductor", None];
  degree = Lookup[context, "degree", None];
  polynomial = Lookup[context, "cyclotomic_polynomial", None];
  IntegerQ[conductor] && conductor >= 1 &&
    context["context_id"] === StringJoin["cyclotomic:", ToString[conductor]] &&
    IntegerQ[degree] && degree >= 1 &&
    ListQ[polynomial] && Length[polynomial] === degree + 1 &&
    Last[polynomial] === 1 &&
    And @@ (cqExactRationalQ /@ polynomial)
];

cqCreateContext[conductor_Integer, maxDegree_Integer : 128] := Module[
  {factors, degree, modulus, key, cached},
  If[conductor < 1, Return[cqFailure["InvalidRootOfUnity", <|"order" -> conductor|>]]];
  If[maxDegree < 1, Return[cqFailure["ConductorDegreeLimitExceeded", <|"limit" -> maxDegree|>]]];
  factors = cqIntegerFactorization[conductor];
  If[FailureQ[factors], Return[factors]];
  degree = cqEulerPhiFromFactors[conductor, factors];
  If[degree > maxDegree,
    Return[cqFailure[
      "ConductorDegreeLimitExceeded",
      <|"conductor" -> conductor, "degree" -> degree, "limit" -> maxDegree|>
    ]]
  ];
  key = ToString[conductor];
  cached = Lookup[$cqContextCache, key, Missing["not_cached"]];
  If[!MissingQ[cached], Return[cached]];
  modulus = cqCyclotomicPolynomial[conductor];
  If[FailureQ[modulus], Return[modulus]];
  If[Length[modulus] =!= degree + 1 || Last[modulus] =!= 1,
    Return[cqFailure[
      "SingularPolynomialElement",
      <|"conductor" -> conductor, "polynomial" -> modulus|>
    ]]
  ];
  cached = <|
    "type" -> "cyclotomic_context",
    "context_id" -> StringJoin["cyclotomic:", key],
    "conductor" -> conductor,
    "degree" -> degree,
    "cyclotomic_polynomial" -> modulus
  |>;
  $cqContextCache[key] = cached;
  cached
];

cqCreateContext[conductor_, maxDegree_: 128] :=
  cqFailure["InvalidRootOfUnity", <|"order" -> HoldForm[conductor], "limit" -> HoldForm[maxDegree]|>];

cqElementDataQ[context_, element_] := AssociationQ[element] &&
  Lookup[element, "type", None] === "cyclotomic_element" &&
  Lookup[element, "context_id", None] === context["context_id"] &&
  ListQ[Lookup[element, "coefficients", None]] &&
  Length[element["coefficients"]] === context["degree"] &&
  And @@ (cqExactRationalQ /@ element["coefficients"]);

cqElementQ[context_, element_] :=
  cqContextQ[context] && cqElementDataQ[context, element];

cqElementFromPolynomial[context_, polynomial_List] := Module[{reduced, padded},
  If[!cqContextQ[context], Return[cqFailure["ConductorMismatch", <||>]]];
  reduced = cqPolynomialMod[polynomial, context["cyclotomic_polynomial"]];
  If[FailureQ[reduced], Return[reduced]];
  padded = cqPolynomialPad[reduced, context["degree"]];
  If[FailureQ[padded], Return[padded]];
  <|
    "type" -> "cyclotomic_element",
    "context_id" -> context["context_id"],
    "coefficients" -> padded
  |>
];

cqElementZero[context_] := cqElementFromPolynomial[context, {0}];
cqElementOne[context_] := cqElementFromPolynomial[context, {1}];

cqCheckElementPair[context_, left_, right_] :=
  If[cqElementQ[context, left] && cqElementQ[context, right],
    True,
    cqFailure["ConductorMismatch", <||>]
  ];

cqElementZeroQ[context_, element_] :=
  If[cqElementQ[context, element], And @@ (# === 0 & /@ element["coefficients"]), False];

cqElementEqualQ[context_, left_, right_] :=
  cqElementQ[context, left] && cqElementQ[context, right] &&
    left["coefficients"] === right["coefficients"];

cqElementAdd[context_, left_, right_] := Module[{check},
  check = cqCheckElementPair[context, left, right];
  If[FailureQ[check], Return[check]];
  cqElementFromPolynomial[context, left["coefficients"] + right["coefficients"]]
];

cqElementNegate[context_, element_] :=
  If[cqElementQ[context, element],
    cqElementFromPolynomial[context, -element["coefficients"]],
    cqFailure["ConductorMismatch", <||>]
  ];

cqElementSubtract[context_, left_, right_] := Module[{negative},
  negative = cqElementNegate[context, right];
  If[FailureQ[negative], Return[negative]];
  cqElementAdd[context, left, negative]
];

cqElementScale[context_, scalar_, element_] := Module[{normalized},
  If[!cqElementQ[context, element], Return[cqFailure["ConductorMismatch", <||>]]];
  normalized = cqNormalizeRational[scalar];
  If[FailureQ[normalized], Return[normalized]];
  cqElementFromPolynomial[context, normalized element["coefficients"]]
];

cqElementMultiply[context_, left_, right_] := Module[{check, product},
  check = cqCheckElementPair[context, left, right];
  If[FailureQ[check], Return[check]];
  product = cqPolynomialMultiply[left["coefficients"], right["coefficients"]];
  If[FailureQ[product], Return[product]];
  cqElementFromPolynomial[context, product]
];

cqElementInverse[context_, element_] := Module[{extended, gcd, inverse},
  If[!cqElementQ[context, element], Return[cqFailure["ConductorMismatch", <||>]]];
  If[cqElementZeroQ[context, element], Return[cqFailure["DivisionByZero", <||>]]];
  extended = cqPolynomialExtendedGCD[
    element["coefficients"],
    context["cyclotomic_polynomial"]
  ];
  If[FailureQ[extended], Return[extended]];
  gcd = extended["gcd"];
  If[gcd =!= {1},
    Return[cqFailure["SingularPolynomialElement", <|"gcd" -> gcd|>]]
  ];
  inverse = cqElementFromPolynomial[context, extended["left_coefficient"]];
  If[FailureQ[inverse], Return[inverse]];
  inverse
];

cqElementDivide[context_, left_, right_] := Module[{inverse},
  inverse = cqElementInverse[context, right];
  If[FailureQ[inverse], Return[inverse]];
  cqElementMultiply[context, left, inverse]
];

cqElementPower[context_, element_, exponent_Integer] := Module[
  {power = exponent, factor = element, result},
  If[!cqElementQ[context, element], Return[cqFailure["ConductorMismatch", <||>]]];
  If[power < 0,
    factor = cqElementInverse[context, factor];
    If[FailureQ[factor], Return[factor]];
    power = -power;
  ];
  result = cqElementOne[context];
  While[power > 0,
    If[OddQ[power],
      result = cqElementMultiply[context, result, factor];
      If[FailureQ[result], Return[result]];
    ];
    power = Quotient[power, 2];
    If[power > 0,
      factor = cqElementMultiply[context, factor, factor];
      If[FailureQ[factor], Return[factor]];
    ];
  ];
  result
];

cqRootOfUnityElement[context_, order_Integer, power_Integer] := Module[
  {conductor, exponent, reducedPower, polynomial},
  If[!cqContextQ[context], Return[cqFailure["ConductorMismatch", <||>]]];
  If[order < 1, Return[cqFailure["InvalidRootOfUnity", <|"order" -> order|>]]];
  conductor = context["conductor"];
  If[Mod[conductor, order] =!= 0,
    Return[cqFailure[
      "ConductorMismatch",
      <|"source_order" -> order, "target_conductor" -> conductor|>
    ]]
  ];
  reducedPower = Mod[power, order];
  exponent = reducedPower Quotient[conductor, order];
  polynomial = cqPolynomialPowerMod[{0, 1}, exponent, context["cyclotomic_polynomial"]];
  If[FailureQ[polynomial], Return[polynomial]];
  cqElementFromPolynomial[context, polynomial]
];

cqRootOfUnityElement[context_, order_, power_] :=
  cqFailure["InvalidRootOfUnity", <|"order" -> HoldForm[order], "power" -> HoldForm[power]|>];

cqElementConjugate[context_, element_] := Module[
  {result, coefficient, basisImage, index, exponent},
  If[!cqElementQ[context, element], Return[cqFailure["ConductorMismatch", <||>]]];
  result = cqElementZero[context];
  For[index = 0, index < context["degree"], index++,
    coefficient = element["coefficients"][[index + 1]];
    If[coefficient =!= 0,
      exponent = Mod[-index, context["conductor"]];
      basisImage = cqPolynomialPowerMod[
        {0, 1}, exponent, context["cyclotomic_polynomial"]
      ];
      If[FailureQ[basisImage], Return[basisImage]];
      basisImage = cqElementFromPolynomial[context, basisImage];
      If[FailureQ[basisImage], Return[basisImage]];
      basisImage = cqElementScale[context, coefficient, basisImage];
      If[FailureQ[basisImage], Return[basisImage]];
      result = cqElementAdd[context, result, basisImage];
      If[FailureQ[result], Return[result]];
    ];
  ];
  result
];

cqEmbedElement[sourceContext_, targetContext_, element_] := Module[
  {sourceConductor, targetConductor, result, index, image, scaled},
  If[!cqElementQ[sourceContext, element] || !cqContextQ[targetContext],
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  sourceConductor = sourceContext["conductor"];
  targetConductor = targetContext["conductor"];
  If[Mod[targetConductor, sourceConductor] =!= 0,
    Return[cqFailure[
      "ConductorMismatch",
      <|"source_conductor" -> sourceConductor, "target_conductor" -> targetConductor|>
    ]]
  ];
  result = cqElementZero[targetContext];
  For[index = 0, index < sourceContext["degree"], index++,
    If[element["coefficients"][[index + 1]] =!= 0,
      image = cqRootOfUnityElement[targetContext, sourceConductor, index];
      If[FailureQ[image], Return[image]];
      scaled = cqElementScale[targetContext, element["coefficients"][[index + 1]], image];
      If[FailureQ[scaled], Return[scaled]];
      result = cqElementAdd[targetContext, result, scaled];
      If[FailureQ[result], Return[result]];
    ];
  ];
  result
];

cqSerializeElement[context_, element_] := Module[{serialized},
  If[!cqElementQ[context, element], Return[cqFailure["ConductorMismatch", <||>]]];
  serialized = cqSerializeRational /@ element["coefficients"];
  If[AnyTrue[serialized, FailureQ], Return[FirstCase[serialized, _Failure]]];
  <|"conductor" -> context["conductor"], "coefficients" -> serialized|>
];

cqDeserializeElement[data_Association, maxDegree_Integer : 128] := Module[
  {conductor, coefficientsRaw, context, coefficients, parsed, element, index},
  conductor = Lookup[data, "conductor", Missing["conductor"]];
  coefficientsRaw = Lookup[data, "coefficients", Missing["coefficients"]];
  If[!IntegerQ[conductor] || conductor < 1 || !ListQ[coefficientsRaw],
    Return[cqFailure["MalformedSerialization", <|"value" -> data|>]]
  ];
  context = cqCreateContext[conductor, maxDegree];
  If[FailureQ[context], Return[context]];
  If[Length[coefficientsRaw] =!= context["degree"],
    Return[cqFailure[
      "MalformedSerialization",
      <|"expected_coefficients" -> context["degree"], "actual_coefficients" -> Length[coefficientsRaw]|>
    ]]
  ];
  coefficients = {};
  For[index = 1, index <= Length[coefficientsRaw], index++,
    parsed = cqParseRational[coefficientsRaw[[index]]];
    If[FailureQ[parsed], Return[parsed]];
    AppendTo[coefficients, parsed];
  ];
  element = cqElementFromPolynomial[context, coefficients];
  If[FailureQ[element], Return[element]];
  <|"context" -> context, "element" -> element|>
];

cqDeserializeElement[data_, maxDegree_: 128] :=
  cqFailure["MalformedSerialization", <|"value" -> HoldForm[data], "limit" -> maxDegree|>];
