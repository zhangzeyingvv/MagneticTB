(* ::Package:: *)

BeginPackage["MagneticTBLinearAlgebra`"]

PeriodicBondClassesWithinRadius::usage =
  "PeriodicBondClassesWithinRadius[positions,lattice,radius] enumerates all directed periodic bonds inside a Cartesian ball and groups them into complete distance shells.";
FindPeriodicBondClasses::usage =
  "FindPeriodicBondClasses[positions,lattice,n] grows a Cartesian search ball until at least n complete distance shells have been found. Set \"ReturnMetadata\"->True to obtain the stopping radius and neighbour counts as well as the bond classes.";

Begin["`Private`"]

ClearAll[
  coordinateVectorQ,
  normalizeBondPositions,
  validBondLatticeQ,
  integerTranslationsInBall,
  splitDistanceShells,
  bondSearchDataWithinRadius,
  initialBondSearchRadius
];

coordinateVectorQ[vector_] :=
  VectorQ[vector, NumericQ] && And @@ (TrueQ[Im[N[#]] == 0.] & /@ vector);

normalizeBondPositions[positions_List] := Which[
  positions =!= {} && And @@ (coordinateVectorQ /@ positions),
    positions,
  positions =!= {} && And @@ (ListQ /@ positions) &&
      And @@ (coordinateVectorQ /@ Flatten[positions, 1]),
    Flatten[positions, 1],
  True,
    $Failed
];

validBondLatticeQ[lattice_, dimension_Integer] :=
  MatrixQ[lattice, NumericQ] &&
    Length[lattice] == dimension &&
    Length[Dimensions[lattice]] == 2 &&
    Last[Dimensions[lattice]] >= dimension &&
    MatrixRank[N[lattice]] == dimension &&
    Max[Abs[Im[N[lattice]]]] == 0.;

(* Cholesky/Fincke--Pohst enumeration of

     (translation + shift).G.(translation + shift)^T <= radius^2,

   where translation is integral.  This enumerates the ellipsoid itself;
   it does not first construct a potentially enormous integer box. *)
integerTranslationsInBall[
  shift_List,
  triangular_?MatrixQ,
  radius_?NumericQ,
  tolerance_?NumericQ
] := Module[
  {dimension, vector, shiftedVector, radiusSquared, absoluteTolerance,
   recurse, harvested},
  dimension = Length[shift];
  vector = ConstantArray[0, dimension];
  shiftedVector = ConstantArray[0., dimension];
  radiusSquared = radius^2;
  absoluteTolerance = tolerance Max[1., radiusSquared];

  recurse[level_Integer, remaining_?NumericQ] := Module[
    {tail, diagonal, centre, halfWidth, lower, upper, term, nextRemaining},
    If[level == 0,
      Sow[vector];
      Return[]
    ];
    tail = Sum[
      triangular[[level, column]] shiftedVector[[column]],
      {column, level + 1, dimension}
    ];
    diagonal = triangular[[level, level]];
    centre = -shift[[level]] - tail/diagonal;
    halfWidth = Sqrt[Max[0., remaining + absoluteTolerance]]/Abs[diagonal];
    lower = Ceiling[centre - halfWidth - 10 tolerance];
    upper = Floor[centre + halfWidth + 10 tolerance];
    Do[
      vector[[level]] = integer;
      shiftedVector[[level]] = integer + shift[[level]];
      term = diagonal shiftedVector[[level]] + tail;
      nextRemaining = remaining - term^2;
      If[nextRemaining >= -absoluteTolerance,
        recurse[level - 1, Max[0., nextRemaining]]
      ],
      {integer, lower, upper}
    ]
  ];

  harvested = Reap[recurse[dimension, radiusSquared]][[2]];
  If[harvested === {}, {}, First[harvested]]
];

splitDistanceShells[sortedRecords_List, tolerance_?NumericQ] := Module[
  {shells = {}, current = {}, reference = 0., closeQ},
  closeQ[first_, second_] :=
    Abs[first - second] <= tolerance Max[1., Abs[first], Abs[second]];
  Scan[
    Function[record,
      If[current === {} || closeQ[First[record], reference],
        If[current === {}, reference = First[record]];
        AppendTo[current, record],
        AppendTo[shells, current];
        current = {record};
        reference = First[record]
      ]
    ],
    sortedRecords
  ];
  If[current =!= {}, AppendTo[shells, current]];
  shells
];

bondSearchDataWithinRadius[
  sites_List,
  lattice_?MatrixQ,
  radius_?NumericQ,
  tolerance_?NumericQ
] := Module[
  {numericLattice, gram, triangular, radiusSquared, absoluteTolerance,
   records, sortedRecords, rawShells, bondClasses, nonzeroSourceIndices,
   neighbourCounts},
  numericLattice = Re[N[lattice]];
  gram = numericLattice.Transpose[numericLattice];
  triangular = CholeskyDecomposition[gram];
  radiusSquared = radius^2;
  absoluteTolerance = tolerance Max[1., radiusSquared];

  records = Flatten[
    Table[
      Module[{source, target, shift, translations, keptTranslations},
        source = sites[[sourceIndex]];
        target = sites[[targetIndex]];
        shift = N[target - source];
        translations = integerTranslationsInBall[
          shift,
          triangular,
          radius,
          tolerance
        ];
        keptTranslations = Select[
          translations,
          Total[(N[target + # - source].numericLattice)^2] <=
            radiusSquared + absoluteTolerance &
        ];
        Map[
          Function[translation,
            {
              Total[(N[target + translation - source].numericLattice)^2],
              sourceIndex,
              targetIndex,
              source,
              target + translation
            }
          ],
          keptTranslations
        ]
      ],
      {sourceIndex, Length[sites]},
      {targetIndex, Length[sites]}
    ],
    2
  ];
  sortedRecords = SortBy[
    records,
    {#[[1]] &, #[[2]] &, #[[3]] &, #[[5]] &}
  ];
  rawShells = splitDistanceShells[sortedRecords, tolerance];
  bondClasses = Map[
    Function[shell,
      Module[{distance, bySource},
        distance = Sqrt[Max[0., First[shell][[1]]]];
        bySource = KeySort@GroupBy[shell, #[[2]] &];
        KeyValueMap[
          Function[{sourceIndex, sourceRecords},
            With[{pairs = SortBy[sourceRecords[[All, {4, 5}]], Last]},
              {distance, Length[pairs], pairs}
            ]
          ],
          bySource
        ]
      ]
    ],
    rawShells
  ];

  nonzeroSourceIndices = Cases[
    records,
    record_ /; First[record] > absoluteTolerance :> record[[2]]
  ];
  neighbourCounts = Lookup[
    Counts[nonzeroSourceIndices],
    Range[Length[sites]],
    0
  ];
  <|
    "BondClasses" -> bondClasses,
    "Radius" -> radius,
    "ShellCount" -> Length[bondClasses],
    "NeighborCountsPerSite" -> neighbourCounts,
    "CandidateBondCount" -> Length[records]
  |>
];

initialBondSearchRadius[sites_List, lattice_?MatrixQ, tolerance_?NumericQ] := Module[
  {dimension, numericLattice, smallTranslations, translationDistances,
   sitePairDistances, positiveDistances},
  dimension = Length[First[sites]];
  numericLattice = Re[N[lattice]];
  smallTranslations = DeleteCases[
    Tuples[{-1, 0, 1}, dimension],
    ConstantArray[0, dimension]
  ];
  translationDistances = Norm[#.numericLattice] & /@ smallTranslations;
  sitePairDistances = Flatten@Table[
    Norm[N[target - source].numericLattice],
    {source, sites},
    {target, sites}
  ];
  positiveDistances = Select[
    Join[translationDistances, sitePairDistances],
    # > 10 tolerance &
  ];
  If[positiveDistances === {}, 1., Min[positiveDistances]]
];

PeriodicBondClassesWithinRadius::positions =
  "Positions must be a nonempty list of real fractional-coordinate vectors, optionally grouped by Wyckoff position.";
PeriodicBondClassesWithinRadius::lattice =
  "The lattice must contain one linearly independent real row vector per fractional-coordinate dimension.";
PeriodicBondClassesWithinRadius::radius =
  "The Cartesian cutoff radius must be a nonnegative real number; received `1`.";
PeriodicBondClassesWithinRadius::option =
  "DistanceTolerance must be positive and ReturnMetadata must be True or False.";

Options[PeriodicBondClassesWithinRadius] = {
  "DistanceTolerance" -> 10^-10,
  "ReturnMetadata" -> False
};

PeriodicBondClassesWithinRadius[
  positions_List,
  lattice_?MatrixQ,
  radius_,
  OptionsPattern[]
] := Module[{sites, tolerance, returnMetadata, data},
  sites = normalizeBondPositions[positions];
  If[sites === $Failed,
    Message[PeriodicBondClassesWithinRadius::positions];
    Return[$Failed]
  ];
  If[!validBondLatticeQ[lattice, Length[First[sites]]],
    Message[PeriodicBondClassesWithinRadius::lattice];
    Return[$Failed]
  ];
  If[!NumericQ[radius] || !TrueQ[radius >= 0],
    Message[PeriodicBondClassesWithinRadius::radius, radius];
    Return[$Failed]
  ];
  tolerance = OptionValue["DistanceTolerance"];
  returnMetadata = OptionValue["ReturnMetadata"];
  If[!NumericQ[tolerance] || !TrueQ[tolerance > 0] ||
      !BooleanQ[returnMetadata],
    Message[PeriodicBondClassesWithinRadius::option]; Return[$Failed]
  ];
  data = bondSearchDataWithinRadius[sites, lattice, N[radius], tolerance];
  If[returnMetadata, data, data["BondClasses"]]
];

FindPeriodicBondClasses::shells =
  "The requested number of complete distance shells must be a positive integer; received `1`.";
FindPeriodicBondClasses::option = "Invalid adaptive-search option value: `1`.";
FindPeriodicBondClasses::limit =
  "The adaptive search did not reach `1` complete shells and `2` neighbours per site after `3` iterations (last radius `4`).";

Options[FindPeriodicBondClasses] = {
  "MinimumNeighborsPerSite" -> 0,
  "InitialRadius" -> Automatic,
  "GrowthFactor" -> 5/4,
  "MaximumIterations" -> 24,
  "MaximumRadius" -> Infinity,
  "DistanceTolerance" -> 10^-10,
  "ReturnMetadata" -> False
};

FindPeriodicBondClasses[
  positions_List,
  lattice_?MatrixQ,
  requestedShells_,
  OptionsPattern[]
] := Module[
  {sites, minimumNeighbours, initialRadius, growthFactor, maximumIterations,
   maximumRadius, tolerance, returnMetadata, radius, data = <||>, success = False,
   iteration = 0, result},
  If[!IntegerQ[requestedShells] || requestedShells < 1,
    Message[FindPeriodicBondClasses::shells, requestedShells];
    Return[$Failed]
  ];
  sites = normalizeBondPositions[positions];
  If[sites === $Failed,
    Message[PeriodicBondClassesWithinRadius::positions];
    Return[$Failed]
  ];
  If[!validBondLatticeQ[lattice, Length[First[sites]]],
    Message[PeriodicBondClassesWithinRadius::lattice];
    Return[$Failed]
  ];

  minimumNeighbours = OptionValue["MinimumNeighborsPerSite"];
  initialRadius = OptionValue["InitialRadius"];
  growthFactor = OptionValue["GrowthFactor"];
  maximumIterations = OptionValue["MaximumIterations"];
  maximumRadius = OptionValue["MaximumRadius"];
  tolerance = OptionValue["DistanceTolerance"];
  returnMetadata = OptionValue["ReturnMetadata"];
  If[
    !IntegerQ[minimumNeighbours] || minimumNeighbours < 0 ||
    !(initialRadius === Automatic ||
      (NumericQ[initialRadius] && TrueQ[initialRadius >= 0])) ||
    !NumericQ[growthFactor] || !TrueQ[growthFactor > 1] ||
    !IntegerQ[maximumIterations] || maximumIterations < 1 ||
    !(maximumRadius === Infinity ||
      (NumericQ[maximumRadius] && TrueQ[maximumRadius > 0])) ||
    !NumericQ[tolerance] || !TrueQ[tolerance > 0] ||
    !BooleanQ[returnMetadata],
    Message[FindPeriodicBondClasses::option, {
      minimumNeighbours, initialRadius, growthFactor, maximumIterations,
      maximumRadius, tolerance, returnMetadata
    }];
    Return[$Failed]
  ];

  radius = If[
    initialRadius === Automatic,
    If[requestedShells == 1 && minimumNeighbours == 0,
      0.,
      initialBondSearchRadius[sites, lattice, tolerance]
    ],
    N[initialRadius]
  ];
  If[maximumRadius =!= Infinity, radius = Min[radius, N[maximumRadius]]];

  Do[
    iteration = currentIteration;
    data = bondSearchDataWithinRadius[sites, lattice, radius, tolerance];
    success = data["ShellCount"] >= requestedShells &&
      Min[data["NeighborCountsPerSite"]] >= minimumNeighbours;
    If[success, Break[]];
    If[maximumRadius =!= Infinity && radius >= N[maximumRadius], Break[]];
    radius = growthFactor Max[radius, 10 Sqrt[tolerance]];
    If[maximumRadius =!= Infinity, radius = Min[radius, N[maximumRadius]]],
    {currentIteration, maximumIterations}
  ];

  If[!success,
    Message[
      FindPeriodicBondClasses::limit,
      requestedShells,
      minimumNeighbours,
      iteration,
      radius
    ];
    Return[Failure[
      "SearchLimitReached",
      Join[data, <|"Iterations" -> iteration|>]
    ]]
  ];
  result = Join[
    data,
    <|
      "RequestedShells" -> requestedShells,
      "MinimumNeighborsPerSite" -> minimumNeighbours,
      "Iterations" -> iteration
    |>
  ];
  If[returnMetadata, result, result["BondClasses"]]
];

End[]

EndPackage[]
