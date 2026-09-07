cqFailure[tag_String, data_: <||>] := Failure[tag, Join[<|"tag" -> tag|>, data]];

cqExactRationalQ[value_] := IntegerQ[value] || Head[value] === Rational;

cqNormalizeRational[value_] :=
  If[cqExactRationalQ[value], value, cqFailure["NonExactInput", <|"value" -> HoldForm[value]|>]];

cqParseInteger[value_Integer] := value;

cqParseInteger[value_String] := Module[
  {characters, sign = 1, digits},
  characters = Characters[value];
  If[characters === {}, Return[cqFailure["MalformedSerialization", <|"value" -> value|>]]];
  If[First[characters] === "-", sign = -1; characters = Rest[characters]];
  If[characters === {} || !And @@ (StringMatchQ[#, DigitCharacter] & /@ characters),
    Return[cqFailure["MalformedSerialization", <|"value" -> value|>]]
  ];
  digits = FromDigits[StringJoin[characters]];
  sign digits
];

cqParseInteger[value_] := cqFailure["MalformedSerialization", <|"value" -> HoldForm[value]|>];

cqParseRational[value_] /; cqExactRationalQ[value] := value;

cqParseRational[data_Association] := Module[
  {numeratorRaw, denominatorRaw, numerator, denominator},
  numeratorRaw = Lookup[data, "numerator", Missing["numerator"]];
  denominatorRaw = Lookup[data, "denominator", Missing["denominator"]];
  If[MissingQ[numeratorRaw] || MissingQ[denominatorRaw],
    Return[cqFailure["MalformedSerialization", <|"value" -> data|>]]
  ];
  numerator = cqParseInteger[numeratorRaw];
  If[FailureQ[numerator], Return[numerator]];
  denominator = cqParseInteger[denominatorRaw];
  If[FailureQ[denominator], Return[denominator]];
  If[denominator === 0, Return[cqFailure["DivisionByZero", <|"value" -> data|>]]];
  numerator/denominator
];

cqParseRational[value_] := cqFailure["NonExactInput", <|"value" -> HoldForm[value]|>];

cqSerializeRational[value_] := Module[{normalized},
  normalized = cqNormalizeRational[value];
  If[FailureQ[normalized], Return[normalized]];
  <|
    "numerator" -> ToString[Numerator[normalized], InputForm],
    "denominator" -> ToString[Denominator[normalized], InputForm]
  |>
];

cqIntegerGCD[a_Integer, b_Integer] := Module[{x = Abs[a], y = Abs[b], remainder},
  While[y =!= 0,
    remainder = Mod[x, y];
    x = y;
    y = remainder;
  ];
  x
];

cqIntegerLCM[a_Integer, b_Integer] :=
  If[a === 0 || b === 0, 0, Abs[Quotient[a, cqIntegerGCD[a, b]] b]];

cqIntegerLCMList[values_List] := Fold[cqIntegerLCM, 1, values];

cqIntegerFactorization[n_Integer] := Module[
  {remaining = n, divisor = 2, exponent, factors = {}},
  If[n < 1, Return[cqFailure["InvalidRootOfUnity", <|"order" -> n|>]]];
  While[divisor divisor <= remaining,
    exponent = 0;
    While[Mod[remaining, divisor] === 0,
      remaining = Quotient[remaining, divisor];
      exponent++;
    ];
    If[exponent > 0, AppendTo[factors, {divisor, exponent}]];
    divisor = If[divisor === 2, 3, divisor + 2];
  ];
  If[remaining > 1, AppendTo[factors, {remaining, 1}]];
  factors
];

cqEulerPhiFromFactors[n_Integer, factors_List] :=
  Fold[Quotient[#1, #2[[1]]] (#2[[1]] - 1) &, n, factors];

cqDivisorsFromFactors[factors_List] := Module[{divisors = {1}, prime, exponent},
  Do[
    prime = factor[[1]];
    exponent = factor[[2]];
    divisors = Flatten[Table[d prime^power, {d, divisors}, {power, 0, exponent}]];
  , {factor, factors}];
  Sort[divisors]
];

