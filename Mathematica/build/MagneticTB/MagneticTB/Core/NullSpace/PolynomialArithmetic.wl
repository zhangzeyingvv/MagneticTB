(* Coefficients are stored from lowest to highest degree. *)

cqPolynomialNormalize[coefficients_List] := Module[{result = coefficients},
  If[result === {}, Return[{0}]];
  If[!And @@ (cqExactRationalQ /@ result),
    Return[cqFailure["NonExactInput", <|"coefficients" -> HoldForm[coefficients]|>]]
  ];
  While[Length[result] > 1 && Last[result] === 0, result = Most[result]];
  result
];

cqPolynomialZeroQ[polynomial_List] := Module[{normalized = cqPolynomialNormalize[polynomial]},
  If[FailureQ[normalized], False, normalized === {0}]
];

cqPolynomialDegree[polynomial_List] := Module[{normalized = cqPolynomialNormalize[polynomial]},
  If[FailureQ[normalized], Return[normalized]];
  If[normalized === {0}, -1, Length[normalized] - 1]
];

cqPolynomialPad[polynomial_List, length_Integer] := Module[{normalized},
  normalized = cqPolynomialNormalize[polynomial];
  If[FailureQ[normalized], Return[normalized]];
  If[Length[normalized] > length,
    Return[cqFailure["MalformedSerialization", <|"polynomial" -> polynomial, "length" -> length|>]]
  ];
  PadRight[normalized, length, 0]
];

cqPolynomialAdd[left_List, right_List] := Module[{size, a, b},
  size = Max[Length[left], Length[right]];
  a = PadRight[left, size, 0];
  b = PadRight[right, size, 0];
  cqPolynomialNormalize[a + b]
];

cqPolynomialNegate[polynomial_List] := cqPolynomialNormalize[-polynomial];

cqPolynomialSubtract[left_List, right_List] :=
  cqPolynomialAdd[left, cqPolynomialNegate[right]];

cqPolynomialScale[scalar_, polynomial_List] := Module[{normalizedScalar},
  normalizedScalar = cqNormalizeRational[scalar];
  If[FailureQ[normalizedScalar], Return[normalizedScalar]];
  cqPolynomialNormalize[normalizedScalar polynomial]
];

cqPolynomialMultiply[left_List, right_List] := Module[
  {a, b, result, i, j},
  a = cqPolynomialNormalize[left];
  If[FailureQ[a], Return[a]];
  b = cqPolynomialNormalize[right];
  If[FailureQ[b], Return[b]];
  If[a === {0} || b === {0}, Return[{0}]];
  result = ConstantArray[0, Length[a] + Length[b] - 1];
  For[i = 1, i <= Length[a], i++,
    For[j = 1, j <= Length[b], j++,
      result[[i + j - 1]] += a[[i]] b[[j]];
    ];
  ];
  cqPolynomialNormalize[result]
];

cqPolynomialDivRem[dividend_List, divisor_List] := Module[
  {remainder, normalizedDivisor, quotient, divisorDegree, lead, shift, factor, term},
  remainder = cqPolynomialNormalize[dividend];
  If[FailureQ[remainder], Return[remainder]];
  normalizedDivisor = cqPolynomialNormalize[divisor];
  If[FailureQ[normalizedDivisor], Return[normalizedDivisor]];
  If[normalizedDivisor === {0}, Return[cqFailure["DivisionByZero", <||>]]];
  divisorDegree = Length[normalizedDivisor] - 1;
  quotient = ConstantArray[0, Max[1, Length[remainder] - divisorDegree]];
  lead = Last[normalizedDivisor];
  While[remainder =!= {0} && Length[remainder] - 1 >= divisorDegree,
    shift = (Length[remainder] - 1) - divisorDegree;
    factor = Last[remainder]/lead;
    quotient[[shift + 1]] += factor;
    term = Join[ConstantArray[0, shift], factor normalizedDivisor];
    remainder = cqPolynomialSubtract[remainder, term];
    If[FailureQ[remainder], Return[remainder]];
  ];
  <|
    "quotient" -> cqPolynomialNormalize[quotient],
    "remainder" -> cqPolynomialNormalize[remainder]
  |>
];

cqPolynomialExactQuotient[dividend_List, divisor_List] := Module[{division},
  division = cqPolynomialDivRem[dividend, divisor];
  If[FailureQ[division], Return[division]];
  If[division["remainder"] =!= {0},
    Return[cqFailure["SingularPolynomialElement", <|"remainder" -> division["remainder"]|>]]
  ];
  division["quotient"]
];

cqPolynomialMod[polynomial_List, modulus_List] := Module[{division},
  division = cqPolynomialDivRem[polynomial, modulus];
  If[FailureQ[division], division, division["remainder"]]
];

cqPolynomialMultiplyMod[left_List, right_List, modulus_List] := Module[{product},
  product = cqPolynomialMultiply[left, right];
  If[FailureQ[product], Return[product]];
  cqPolynomialMod[product, modulus]
];

cqPolynomialPowerMod[base_List, exponent_Integer, modulus_List] := Module[
  {power = exponent, factor, result = {1}},
  If[exponent < 0, Return[cqFailure["InvalidRootOfUnity", <|"power" -> exponent|>]]];
  factor = cqPolynomialMod[base, modulus];
  If[FailureQ[factor], Return[factor]];
  While[power > 0,
    If[OddQ[power],
      result = cqPolynomialMultiplyMod[result, factor, modulus];
      If[FailureQ[result], Return[result]];
    ];
    power = Quotient[power, 2];
    If[power > 0,
      factor = cqPolynomialMultiplyMod[factor, factor, modulus];
      If[FailureQ[factor], Return[factor]];
    ];
  ];
  cqPolynomialNormalize[result]
];

cqPolynomialExtendedGCD[left_List, right_List] := Module[
  {oldR, r, oldS = {1}, s = {0}, oldT = {0}, t = {1}, division, q,
   nextR, nextS, nextT, scale},
  oldR = cqPolynomialNormalize[left];
  If[FailureQ[oldR], Return[oldR]];
  r = cqPolynomialNormalize[right];
  If[FailureQ[r], Return[r]];
  While[r =!= {0},
    division = cqPolynomialDivRem[oldR, r];
    If[FailureQ[division], Return[division]];
    q = division["quotient"];
    nextR = division["remainder"];
    nextS = cqPolynomialSubtract[oldS, cqPolynomialMultiply[q, s]];
    If[FailureQ[nextS], Return[nextS]];
    nextT = cqPolynomialSubtract[oldT, cqPolynomialMultiply[q, t]];
    If[FailureQ[nextT], Return[nextT]];
    oldR = r; r = nextR;
    oldS = s; s = nextS;
    oldT = t; t = nextT;
  ];
  If[oldR === {0}, Return[<|"gcd" -> {0}, "left_coefficient" -> {0}, "right_coefficient" -> {0}|>]];
  scale = 1/Last[oldR];
  <|
    "gcd" -> cqPolynomialScale[scale, oldR],
    "left_coefficient" -> cqPolynomialScale[scale, oldS],
    "right_coefficient" -> cqPolynomialScale[scale, oldT]
  |>
];

cqCyclotomicPolynomial[n_Integer] := Module[
  {factors, divisors, computed = <||>, current, polynomial, properDivisor, quotient,
   properDivisors, divisorIndex, properIndex},
  If[n < 1, Return[cqFailure["InvalidRootOfUnity", <|"order" -> n|>]]];
  factors = cqIntegerFactorization[n];
  If[FailureQ[factors], Return[factors]];
  divisors = cqDivisorsFromFactors[factors];
  For[divisorIndex = 1, divisorIndex <= Length[divisors], divisorIndex++,
    current = divisors[[divisorIndex]];
    polynomial = Join[{-1}, ConstantArray[0, current - 1], {1}];
    properDivisors = Select[divisors, Divisible[current, #] && # < current &];
    For[properIndex = 1, properIndex <= Length[properDivisors], properIndex++,
      properDivisor = properDivisors[[properIndex]];
      quotient = cqPolynomialExactQuotient[polynomial, computed[ToString[properDivisor]]];
      If[FailureQ[quotient], Return[quotient]];
      polynomial = quotient;
    ];
    computed[ToString[current]] = cqPolynomialNormalize[polynomial];
  ];
  computed[ToString[n]]
];
