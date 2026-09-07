(* Deterministic sparsity-aware row blocking for iterative common kernels.
   Blocking changes only evaluation order; all rank decisions remain exact. *)

$cqAdaptiveInitialDenseRows = 20;
$cqAdaptiveMinimumDenseRows = 8;
$cqAdaptiveMaximumDenseRows = 40;
$cqAdaptiveMaximumBlockRows = 40;

cqAdaptiveRowNonzeroCount[row_List, zeroQ_] :=
  Count[zeroQ /@ row, False];

cqAdaptivePrepareRows[data_List, zeroQ_] := Module[
  {indexed},
  indexed = MapIndexed[
    {First[#2], #1, cqAdaptiveRowNonzeroCount[#1, zeroQ]} &,
    DeleteDuplicates[data]
  ];
  indexed = Select[indexed, #[[3]] > 0 &];
  #[[2]] & /@ SortBy[indexed, {#[[3]], #[[1]]} &]
];

cqAdaptiveTakeRowBlock[
  rows_List,
  start_Integer,
  basisColumns_Integer,
  originalColumns_Integer,
  targetDenseRows_Integer,
  zeroQ_
] := Module[
  {remaining, dimensionLimit, maxRows, effectiveDenseRows, targetWeight,
   end, rowCount = 0, nonzeroWeight = 0},
  remaining = Length[rows] - start + 1;
  If[remaining <= 0,
    Return[<|"data" -> {}, "next_start" -> start, "row_count" -> 0,
      "nonzero_weight" -> 0|>]
  ];
  dimensionLimit = Max[$cqAdaptiveMinimumDenseRows, 2 basisColumns];
  maxRows = Min[$cqAdaptiveMaximumBlockRows, dimensionLimit, remaining];
  effectiveDenseRows = Min[targetDenseRows, dimensionLimit];
  targetWeight = Max[1, effectiveDenseRows Max[1, originalColumns]];
  end = start - 1;
  While[
    end < Length[rows] && rowCount < maxRows &&
      (rowCount === 0 || nonzeroWeight < targetWeight),
    end++;
    rowCount++;
    nonzeroWeight += cqAdaptiveRowNonzeroCount[rows[[end]], zeroQ]
  ];
  <|
    "data" -> rows[[start ;; end]],
    "next_start" -> end + 1,
    "row_count" -> rowCount,
    "nonzero_weight" -> nonzeroWeight
  |>
];

cqAdaptiveNextDenseRows[
  current_Integer,
  rankGain_Integer,
  potentialRankGain_Integer
] := Which[
  potentialRankGain <= 0,
    current,
  rankGain === 0,
    Min[$cqAdaptiveMaximumDenseRows, 2 current],
  3 rankGain >= 2 potentialRankGain,
    Max[$cqAdaptiveMinimumDenseRows, Floor[3 current/4]],
  4 rankGain <= potentialRankGain,
    Min[$cqAdaptiveMaximumDenseRows, Ceiling[3 current/2]],
  True,
    current
];
