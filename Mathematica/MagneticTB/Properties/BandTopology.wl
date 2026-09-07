(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  topologyFiniteNumberQ,
  topologyRealVectorQ,
  topologyOptionNames,
  topologyUnknownOptions,
  topologyHamiltonianRecord,
  topologyCubeTriangles,
  topologyOccupiedLink,
  topologyTriangleFlux,
  topologyBrillouinZoneQ,
  topologyCanonicalPoint,
  topologyPeriodicDistance,
  topologyGridCoordinates,
  pointChernNumber,
  findGaplessPoints
];

topologyFiniteNumberQ[value_] :=
  NumericQ[value] && Quiet@Check[
    FreeQ[N[value], Indeterminate | ComplexInfinity | DirectedInfinity],
    False
  ];

topologyRealVectorQ[vector_, dimension_Integer?Positive] :=
  ListQ[vector] && Length[vector] === dimension &&
    VectorQ[vector, topologyFiniteNumberQ] &&
    TrueQ[Max[Abs[Im[N[vector]]]] == 0];

topologyOptionNames[function_] := First /@ Options[function];
topologyUnknownOptions[function_, supplied_List] := Complement[
  DeleteDuplicates[
    First /@ Cases[supplied, _Rule | _RuleDelayed, {1}]
  ],
  topologyOptionNames[function]
];

topologyHamiltonianRecord[
    hamiltonian_, point_List, momentumSymbols_List,
    hermitianTolerance_, includeEigenvectors_: True
  ] := Module[
  {
    matrix, dimensions, hermitianResidual, values, vectors = None,
    order, record
  },
  matrix = Quiet@Check[
    If[
      MatrixQ[hamiltonian],
      N[hamiltonian /. Thread[
        Take[momentumSymbols, Length[point]] -> point
      ]],
      N[hamiltonian[point]]
    ],
    $Failed
  ];
  If[
    matrix === $Failed || !MatrixQ[matrix, topologyFiniteNumberQ],
    Return[Failure[
      "InvalidHamiltonian",
      <|"Point" -> point|>
    ]]
  ];
  dimensions = Dimensions[matrix];
  If[
    Length[dimensions] =!= 2 || First[dimensions] =!= Last[dimensions] ||
      First[dimensions] < 2,
    Return[Failure[
      "InvalidHamiltonian",
      <|"Point" -> point, "Dimensions" -> dimensions|>
    ]]
  ];
  matrix = Normal[matrix];
  hermitianResidual = Max@Abs@Flatten[
    matrix - ConjugateTranspose[matrix]
  ];
  If[!topologyFiniteNumberQ[hermitianResidual] ||
      TrueQ[hermitianResidual > hermitianTolerance],
    Return[Failure[
      "NonHermitianHamiltonian",
      <|
        "Point" -> point,
        "Residual" -> hermitianResidual
      |>
    ]]
  ];
  If[
    TrueQ[includeEigenvectors],
    {values, vectors} = Quiet@Check[
      Eigensystem[matrix],
      {$Failed, $Failed}
    ],
    values = Quiet@Check[Eigenvalues[matrix], $Failed]
  ];
  If[
    values === $Failed || !VectorQ[values, topologyFiniteNumberQ] ||
      (TrueQ[includeEigenvectors] &&
        !MatrixQ[vectors, topologyFiniteNumberQ]),
    Return[Failure[
      "EigensystemFailed",
      <|"Point" -> point|>
    ]]
  ];
  order = Ordering[Re[values]];
  record = <|
    "Point" -> point,
    "Dimension" -> First[dimensions],
    "Eigenvalues" -> values[[order]],
    "HermitianResidual" -> hermitianResidual
  |>;
  If[
    TrueQ[includeEigenvectors],
    AssociateTo[
      record,
      <|
        "Matrix" -> matrix,
        "Eigenvectors" -> vectors[[order]]
      |>
    ]
  ];
  record
];

(* Each square is split along its lower-left to upper-right diagonal. The
   orientation signs make the six cube faces point outward. *)
topologyCubeTriangles[center_List, radius_, subdivisions_Integer] := Module[
  {coordinates, faces, harvested, face, orientation, points, triangles},
  coordinates = Subdivide[-radius, radius, subdivisions];
  faces = {
    {Function[{u, v}, center + {radius, u, v}], 1},
    {Function[{u, v}, center + {-radius, u, v}], -1},
    {Function[{u, v}, center + {u, radius, v}], -1},
    {Function[{u, v}, center + {u, -radius, v}], 1},
    {Function[{u, v}, center + {u, v, radius}], 1},
    {Function[{u, v}, center + {u, v, -radius}], -1}
  };
  harvested = Reap[
    Do[
      face = faces[[faceIndex, 1]];
      orientation = faces[[faceIndex, 2]];
      points = {
        face[coordinates[[firstIndex]], coordinates[[secondIndex]]],
        face[coordinates[[firstIndex + 1]], coordinates[[secondIndex]]],
        face[coordinates[[firstIndex + 1]], coordinates[[secondIndex + 1]]],
        face[coordinates[[firstIndex]], coordinates[[secondIndex + 1]]]
      };
      triangles = {
        points[[{1, 2, 3}]],
        points[[{1, 3, 4}]]
      };
      If[orientation < 0, triangles = Reverse /@ triangles];
      Scan[Sow, triangles],
      {faceIndex, Length[faces]},
      {firstIndex, subdivisions},
      {secondIndex, subdivisions}
    ]
  ][[2]];
  If[harvested === {}, {}, First[harvested]]
];

topologyOccupiedLink[firstFrame_, secondFrame_, tolerance_] := Module[
  {overlap, singularValues, determinant},
  overlap = Conjugate[firstFrame] . Transpose[secondFrame];
  singularValues = Quiet@Check[SingularValueList[overlap], $Failed];
  If[
    singularValues === $Failed || singularValues === {} ||
      !VectorQ[singularValues, topologyFiniteNumberQ],
    Return[Failure["OverlapFailed", <||>]]
  ];
  If[TrueQ[Min[singularValues] <= tolerance],
    Return[Failure[
      "SingularOverlap",
      <|"MinimumSingularValue" -> Min[singularValues]|>
    ]]
  ];
  determinant = Det[overlap];
  If[
    !topologyFiniteNumberQ[determinant] ||
      TrueQ[Abs[determinant] <= tolerance],
    Return[Failure[
      "SingularOverlap",
      <|"MinimumSingularValue" -> Min[singularValues]|>
    ]]
  ];
  <|
    "Phase" -> determinant/Abs[determinant],
    "MinimumSingularValue" -> Min[singularValues]
  |>
];

topologyTriangleFlux[frames_List, tolerance_] := Module[
  {links, failure},
  links = {
    topologyOccupiedLink[frames[[1]], frames[[2]], tolerance],
    topologyOccupiedLink[frames[[2]], frames[[3]], tolerance],
    topologyOccupiedLink[frames[[3]], frames[[1]], tolerance]
  };
  failure = SelectFirst[links, FailureQ, Missing["NoFailure"]];
  If[!MissingQ[failure], Return[failure]];
  <|
    "Flux" -> Arg[Times @@ Lookup[links, "Phase"]],
    "MinimumSingularValue" -> Min[Lookup[links, "MinimumSingularValue"]]
  |>
];

pointChernNumber::point =
  "The point must be a finite real three-vector; received `1`.";
pointChernNumber::radius =
  "The enclosing cube radius must be a finite positive real number; received `1`.";
pointChernNumber::ham =
  "The Hamiltonian at `1` must be a finite numeric square Hermitian matrix of constant dimension; detail: `2`.";
pointChernNumber::occ =
  "The occupied-band count must be an integer from 1 through dimension-1; received `1` for dimension `2`.";
pointChernNumber::center =
  "The specified point is not gapless within tolerance: direct gap `1`, tolerance `2`. Set RequireGaplessCenter -> False only when intentionally measuring all charge inside the enclosing surface.";
pointChernNumber::surfacegap =
  "The occupied subspace is not isolated on the enclosing surface at `1`: direct gap `2` does not exceed tolerance `3`.";
pointChernNumber::overlap =
  "Adjacent occupied subspaces on the surface mesh have a singular overlap at triangle `1` (minimum singular value `2`); increase SurfaceSubdivisions or reduce the radius.";
pointChernNumber::quantization =
  "The computed surface flux `1` is not within `2` of an integer; increase SurfaceSubdivisions or check that the enclosing surface is gapped.";
pointChernNumber::option =
  "Invalid pointChernNumber option value(s): `1`.";

Options[pointChernNumber] = {
  "SurfaceSubdivisions" -> 8,
  "HermitianTolerance" -> 10^-10,
  "SurfaceGapTolerance" -> 10^-8,
  "CenterGapTolerance" -> 10^-6,
  "OverlapTolerance" -> 10^-10,
  "IntegerTolerance" -> 5 10^-3,
  "RequireGaplessCenter" -> True,
  "MomentumSymbols" -> {kx, ky, kz},
  "Output" -> "ChernNumber"
};

pointChernNumber[
    hamiltonian_, occupied_, point_, radius_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, subdivisions,
    hermitianTolerance, surfaceGapTolerance, centerGapTolerance,
    overlapTolerance, integerTolerance, requireGaplessCenter,
    momentumSymbols, output, center, numericRadius, centerRecord,
    dimension, centerGap, triangles, frameCache = <||>, frameAt,
    triangleResults, failurePosition, failure, surfaceGaps,
    minimumSurfaceGap, minimumSingularValue, rawChern, integerChern,
    quantizationError, data
  },
  unknown = topologyUnknownOptions[pointChernNumber, supplied];
  If[unknown =!= {},
    Message[pointChernNumber::option, unknown];
    Return[$Failed]
  ];
  subdivisions = OptionValue["SurfaceSubdivisions"];
  hermitianTolerance = OptionValue["HermitianTolerance"];
  surfaceGapTolerance = OptionValue["SurfaceGapTolerance"];
  centerGapTolerance = OptionValue["CenterGapTolerance"];
  overlapTolerance = OptionValue["OverlapTolerance"];
  integerTolerance = OptionValue["IntegerTolerance"];
  requireGaplessCenter = OptionValue["RequireGaplessCenter"];
  momentumSymbols = OptionValue["MomentumSymbols"];
  output = OptionValue["Output"];
  If[!topologyRealVectorQ[point, 3],
    Message[pointChernNumber::point, point];
    Return[$Failed]
  ];
  If[
    !topologyFiniteNumberQ[radius] || !TrueQ[Im[N[radius]] == 0] ||
      !TrueQ[radius > 0],
    Message[pointChernNumber::radius, radius];
    Return[$Failed]
  ];
  If[
    !IntegerQ[subdivisions] || subdivisions < 2 ||
      !And @@ (
        topologyFiniteNumberQ[#] && TrueQ[Im[N[#]] == 0] &&
          TrueQ[# >= 0] & /@ {
          hermitianTolerance, surfaceGapTolerance, centerGapTolerance,
          overlapTolerance, integerTolerance
        }
      ) ||
      !BooleanQ[requireGaplessCenter] ||
      !ListQ[momentumSymbols] || Length[momentumSymbols] < 3 ||
      !VectorQ[momentumSymbols, MatchQ[#, _Symbol] &] ||
      !DuplicateFreeQ[momentumSymbols] ||
      !MemberQ[{"ChernNumber", "Data"}, output],
    Message[pointChernNumber::option, {
      subdivisions, hermitianTolerance, surfaceGapTolerance,
      centerGapTolerance, overlapTolerance, integerTolerance,
      requireGaplessCenter, momentumSymbols, output
    }];
    Return[$Failed]
  ];
  center = N[point];
  numericRadius = N[radius];
  centerRecord = topologyHamiltonianRecord[
    hamiltonian,
    center,
    momentumSymbols,
    hermitianTolerance
  ];
  If[FailureQ[centerRecord],
    Message[
      pointChernNumber::ham,
      center,
      centerRecord[[1]]
    ];
    Return[$Failed]
  ];
  dimension = centerRecord["Dimension"];
  If[!IntegerQ[occupied] || !Between[occupied, {1, dimension - 1}],
    Message[pointChernNumber::occ, occupied, dimension];
    Return[$Failed]
  ];
  centerGap = Re[
    centerRecord["Eigenvalues"][[occupied + 1]] -
      centerRecord["Eigenvalues"][[occupied]]
  ];
  If[requireGaplessCenter && TrueQ[centerGap > centerGapTolerance],
    Message[
      pointChernNumber::center,
      centerGap,
      centerGapTolerance
    ];
    Return[$Failed]
  ];
  triangles = topologyCubeTriangles[
    center,
    numericRadius,
    subdivisions
  ];
  frameAt[coordinate_List] := Module[
    {key, cached, record, gap, result},
    key = ToString[N[coordinate, 16], InputForm];
    cached = Lookup[frameCache, key, Missing["NotCached"]];
    If[!MissingQ[cached], Return[cached]];
    record = topologyHamiltonianRecord[
      hamiltonian,
      coordinate,
      momentumSymbols,
      hermitianTolerance
    ];
    If[FailureQ[record],
      result = Failure[
        "HamiltonianFailure",
        <|"Point" -> coordinate, "Detail" -> record[[1]]|>
      ],
      If[record["Dimension"] =!= dimension,
        result = Failure[
          "HamiltonianFailure",
          <|"Point" -> coordinate, "Detail" -> "DimensionChanged"|>
        ],
        gap = Re[
          record["Eigenvalues"][[occupied + 1]] -
            record["Eigenvalues"][[occupied]]
        ];
        result = If[
          !topologyFiniteNumberQ[gap] ||
            TrueQ[gap <= surfaceGapTolerance],
          Failure[
            "SurfaceGapClosed",
            <|"Point" -> coordinate, "Gap" -> gap|>
          ],
          <|
            "Frame" -> record["Eigenvectors"][[;; occupied]],
            "Gap" -> gap,
            "HermitianResidual" -> record["HermitianResidual"]
          |>
        ]
      ]
    ];
    AssociateTo[frameCache, key -> result];
    result
  ];
  triangleResults = Table[
    Module[{vertexRecords, vertexFailure, fluxRecord},
      vertexRecords = frameAt /@ triangles[[triangleIndex]];
      vertexFailure = SelectFirst[
        vertexRecords,
        FailureQ,
        Missing["NoFailure"]
      ];
      If[!MissingQ[vertexFailure],
        vertexFailure,
        fluxRecord = topologyTriangleFlux[
          Lookup[vertexRecords, "Frame"],
          overlapTolerance
        ];
        If[
          FailureQ[fluxRecord],
          Failure[
            "SingularOverlap",
            <|
              "Triangle" -> triangleIndex,
              "MinimumSingularValue" -> Lookup[
                fluxRecord[[2]],
                "MinimumSingularValue",
                Missing["NotAvailable"]
              ]
            |>
          ],
          Join[
            fluxRecord,
            <|"MinimumGap" -> Min[Lookup[vertexRecords, "Gap"]]|>
          ]
        ]
      ]
    ],
    {triangleIndex, Length[triangles]}
  ];
  failurePosition = FirstPosition[triangleResults, _?FailureQ];
  If[!MissingQ[failurePosition],
    failure = Extract[triangleResults, failurePosition];
    Switch[
      failure[[1]],
      "HamiltonianFailure",
        Message[
          pointChernNumber::ham,
          Lookup[failure[[2]], "Point", Missing["NotAvailable"]],
          Lookup[failure[[2]], "Detail", failure[[1]]]
        ],
      "SurfaceGapClosed",
        Message[
          pointChernNumber::surfacegap,
          failure[[2, "Point"]],
          failure[[2, "Gap"]],
          surfaceGapTolerance
        ],
      "SingularOverlap",
        Message[
          pointChernNumber::overlap,
          failure[[2, "Triangle"]],
          failure[[2, "MinimumSingularValue"]]
        ]
    ];
    Return[$Failed]
  ];
  surfaceGaps = Lookup[triangleResults, "MinimumGap"];
  minimumSurfaceGap = Min[surfaceGaps];
  minimumSingularValue = Min[
    Lookup[triangleResults, "MinimumSingularValue"]
  ];
  rawChern = Total[Lookup[triangleResults, "Flux"]]/(2 Pi);
  integerChern = Round[rawChern];
  quantizationError = Abs[rawChern - integerChern];
  If[TrueQ[quantizationError > integerTolerance],
    Message[
      pointChernNumber::quantization,
      rawChern,
      integerTolerance
    ];
    Return[$Failed]
  ];
  data = <|
    "Schema" -> "MagneticTBPointChernNumber",
    "SchemaVersion" -> 1,
    "ChernNumber" -> integerChern,
    "RawChernNumber" -> rawChern,
    "QuantizationError" -> quantizationError,
    "Point" -> point,
    "Radius" -> radius,
    "OccupiedBands" -> occupied,
    "HamiltonianDimension" -> dimension,
    "CenterGap" -> centerGap,
    "RequireGaplessCenter" -> requireGaplessCenter,
    "MinimumSurfaceGap" -> minimumSurfaceGap,
    "MinimumOverlapSingularValue" -> minimumSingularValue,
    "SurfaceSubdivisions" -> subdivisions,
    "TriangleCount" -> Length[triangles],
    "UniqueVertexCount" -> Length[frameCache],
    "Method" -> "gauge-invariant occupied-subspace flux on an oriented cube mesh"
  |>;
  If[output === "Data", data, integerChern]
];

topologyBrillouinZoneQ[zone_] :=
  ListQ[zone] && Between[Length[zone], {1, 3}] &&
    And @@ Map[
      Function[interval,
        ListQ[interval] && Length[interval] === 2 &&
          topologyRealVectorQ[interval, 2] &&
          TrueQ[interval[[1]] < interval[[2]]]
      ],
      zone
    ];

topologyCanonicalPoint[point_List, zone_List] := MapThread[
  #2[[1]] + Mod[#1 - #2[[1]], #2[[2]] - #2[[1]]] &,
  {point, zone}
];

topologyPeriodicDistance[first_List, second_List, zone_List] := Module[
  {widths, differences},
  widths = #[[2]] - #[[1]] & /@ zone;
  differences = Abs[first - second];
  Norm[MapThread[Min[#1, #2 - #1] &, {differences, widths}]]
];

topologyGridCoordinates[zone_List, gridSize_List] := MapThread[
  Most@Subdivide[#1[[1]], #1[[2]], #2] &,
  {zone, gridSize}
];

findGaplessPoints::bz =
  "BrillouinZone must contain one to three finite real {minimum,maximum} intervals with minimum<maximum; received `1`.";
findGaplessPoints::grid =
  "GridSize must be an integer at least 2 or one such integer per Brillouin-zone dimension; received `1`.";
findGaplessPoints::ham = pointChernNumber::ham;
findGaplessPoints::occ = pointChernNumber::occ;
findGaplessPoints::refine =
  "Gap refinement encountered an invalid Hamiltonian at `1`; detail: `2`.";
findGaplessPoints::option =
  "Invalid findGaplessPoints option value(s): `1`.";

Options[findGaplessPoints] = {
  "BrillouinZone" -> ConstantArray[{-Pi, Pi}, 3],
  "GridSize" -> 15,
  "CandidateCount" -> 32,
  "GapTolerance" -> 10^-7,
  "MergeTolerance" -> 10^-4,
  "HermitianTolerance" -> 10^-10,
  "MaxIterations" -> 500,
  "RefinementMethod" -> "PrincipalAxis",
  "MomentumSymbols" -> {kx, ky, kz},
  "Output" -> "Points"
};

findGaplessPoints[
    hamiltonian_, occupied_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, zone, gridSetting,
    candidateCount, gapTolerance, mergeTolerance,
    hermitianTolerance, maxIterations, refinementMethod,
    momentumSymbols, output,
    dimension, gridSize, axisCoordinates, gridPoints,
    spectrumCache = <||>, spectrumAt, firstRecord, hamiltonianDimension,
    gapAt, gridGaps, failurePosition, failure, gapArray,
    gridIndices, neighborOffsets, localMinimumIndices,
    pointFromIndex, localRecords, gridRecords, localSeeds,
    extraGridRecords, seedRecords,
    variables, objective, refinementFailure = None,
    refinementRecords, refinedResult, refinedPoint, refinedRecord,
    qualifying, merged, finalRecords, minimumGap, failedRefinements,
    data
  },
  unknown = topologyUnknownOptions[findGaplessPoints, supplied];
  If[unknown =!= {},
    Message[findGaplessPoints::option, unknown];
    Return[$Failed]
  ];
  zone = OptionValue["BrillouinZone"];
  gridSetting = OptionValue["GridSize"];
  candidateCount = OptionValue["CandidateCount"];
  gapTolerance = OptionValue["GapTolerance"];
  mergeTolerance = OptionValue["MergeTolerance"];
  hermitianTolerance = OptionValue["HermitianTolerance"];
  maxIterations = OptionValue["MaxIterations"];
  refinementMethod = OptionValue["RefinementMethod"];
  momentumSymbols = OptionValue["MomentumSymbols"];
  output = OptionValue["Output"];
  If[!topologyBrillouinZoneQ[zone],
    Message[findGaplessPoints::bz, zone];
    Return[$Failed]
  ];
  dimension = Length[zone];
  gridSize = Which[
    IntegerQ[gridSetting], ConstantArray[gridSetting, dimension],
    ListQ[gridSetting] && Length[gridSetting] === dimension,
      gridSetting,
    True, $Failed
  ];
  If[
    gridSize === $Failed || !VectorQ[gridSize, IntegerQ] ||
      Min[gridSize] < 2,
    Message[findGaplessPoints::grid, gridSetting];
    Return[$Failed]
  ];
  If[
    !IntegerQ[candidateCount] || candidateCount < 1 ||
      !IntegerQ[maxIterations] || maxIterations < 1 ||
      !And @@ (
        topologyFiniteNumberQ[#] && TrueQ[Im[N[#]] == 0] &&
          TrueQ[# >= 0] & /@ {
          gapTolerance, mergeTolerance, hermitianTolerance
        }
      ) ||
      !ListQ[momentumSymbols] || Length[momentumSymbols] < dimension ||
      !VectorQ[momentumSymbols, MatchQ[#, _Symbol] &] ||
      !DuplicateFreeQ[momentumSymbols] ||
      !MemberQ[{"QuasiNewton", "PrincipalAxis"}, refinementMethod] ||
      !MemberQ[{"Points", "Data"}, output],
    Message[findGaplessPoints::option, {
      candidateCount, gapTolerance, mergeTolerance,
      hermitianTolerance, maxIterations, refinementMethod,
      momentumSymbols, output
    }];
    Return[$Failed]
  ];
  axisCoordinates = topologyGridCoordinates[N[zone], gridSize];
  gridPoints = Tuples[axisCoordinates];
  spectrumAt[point_List] := Module[{canonical, key, cached, record},
    canonical = topologyCanonicalPoint[N[point], N[zone]];
    key = ToString[N[canonical, 16], InputForm];
    cached = Lookup[spectrumCache, key, Missing["NotCached"]];
    If[!MissingQ[cached], Return[cached]];
    record = topologyHamiltonianRecord[
      hamiltonian,
      canonical,
      momentumSymbols,
      hermitianTolerance,
      False
    ];
    AssociateTo[spectrumCache, key -> record];
    record
  ];
  firstRecord = spectrumAt[First[gridPoints]];
  If[FailureQ[firstRecord],
    Message[
      findGaplessPoints::ham,
      First[gridPoints],
      firstRecord[[1]]
    ];
    Return[$Failed]
  ];
  hamiltonianDimension = firstRecord["Dimension"];
  If[
    !IntegerQ[occupied] ||
      !Between[occupied, {1, hamiltonianDimension - 1}],
    Message[
      findGaplessPoints::occ,
      occupied,
      hamiltonianDimension
    ];
    Return[$Failed]
  ];
  gapAt[point_List] := Module[{record},
    record = spectrumAt[point];
    If[
      FailureQ[record],
      record,
      If[record["Dimension"] =!= hamiltonianDimension,
        Failure[
          "DimensionChanged",
          <|"Point" -> point|>
        ],
        Re[
          record["Eigenvalues"][[occupied + 1]] -
            record["Eigenvalues"][[occupied]]
        ]
      ]
    ]
  ];
  gridGaps = gapAt /@ gridPoints;
  failurePosition = FirstPosition[gridGaps, _?FailureQ];
  If[!MissingQ[failurePosition],
    failure = Extract[gridGaps, failurePosition];
    Message[
      findGaplessPoints::ham,
      gridPoints[[First[failurePosition]]],
      failure[[1]]
    ];
    Return[$Failed]
  ];
  gapArray = ArrayReshape[gridGaps, gridSize];
  gridIndices = Tuples[Range /@ gridSize];
  neighborOffsets = DeleteCases[
    Tuples[{-1, 0, 1}, dimension],
    ConstantArray[0, dimension]
  ];
  localMinimumIndices = Select[
    gridIndices,
    Function[index,
      With[{value = Extract[gapArray, index]},
        And @@ (
          TrueQ[value <= Extract[
            gapArray,
            1 + Mod[index - 1 + #, gridSize]
          ]] & /@ neighborOffsets
        )
      ]
    ]
  ];
  pointFromIndex[index_List] := MapThread[
    #1[[#2]] &,
    {axisCoordinates, index}
  ];
  localRecords = SortBy[
    <|
      "Point" -> pointFromIndex[#],
      "Gap" -> Extract[gapArray, #],
      "Source" -> "LocalGridMinimum"
    |> & /@ localMinimumIndices,
    Lookup[#, "Gap"] &
  ];
  gridRecords = SortBy[
    MapThread[
      <|"Point" -> #1, "Gap" -> #2, "Source" -> "SmallGridGap"|> &,
      {gridPoints, gridGaps}
    ],
    Lookup[#, "Gap"] &
  ];
  localSeeds = Take[localRecords, UpTo[candidateCount]];
  extraGridRecords = Select[
    gridRecords,
    !MemberQ[Lookup[localSeeds, "Point", {}], Lookup[#, "Point"]] &
  ];
  seedRecords = DeleteDuplicatesBy[
    Join[
      localSeeds,
      Take[
        extraGridRecords,
        UpTo[Max[0, candidateCount - Length[localSeeds]]]
      ]
    ],
    ToString[N[Lookup[#, "Point"], 16], InputForm] &
  ];
  variables = Array[Unique["gaplessCoordinate"] &, dimension];
  objective[arguments__] /; VectorQ[{arguments}, NumericQ] := Module[
    {value},
    value = gapAt[topologyCanonicalPoint[{arguments}, N[zone]]];
    If[FailureQ[value],
      refinementFailure = value;
      10.^100,
      value^2
    ]
  ];
  refinementRecords = Catch[
    Table[
      If[
        TrueQ[seedRecord["Gap"] <= gapTolerance],
        Join[
          seedRecord,
          <|
            "GridSeed" -> seedRecord["Point"],
            "GridGap" -> seedRecord["Gap"],
            "Refined" -> False,
            "RefinementStatus" -> "NotNeeded"
          |>
        ],
        refinementFailure = None;
        refinedResult = Quiet@Check[
          FindMinimum[
            objective @@ variables,
            Evaluate[Thread[{variables, seedRecord["Point"]}]],
            Method -> refinementMethod,
            MaxIterations -> maxIterations
          ],
          $Failed
        ];
        If[FailureQ[refinementFailure],
          Throw[
            Failure[
              "RefinementHamiltonianFailure",
              <|
                "Point" -> Lookup[
                  refinementFailure[[2]],
                  "Point",
                  seedRecord["Point"]
                ],
                "Detail" -> refinementFailure[[1]]
              |>
            ],
            "GaplessRefinementFailure"
          ]
        ];
        If[
          refinedResult === $Failed,
          Join[
            seedRecord,
            <|
              "Refined" -> False,
              "RefinementStatus" -> "Failed"
            |>
          ],
          refinedPoint = topologyCanonicalPoint[
            N[variables /. Last[refinedResult]],
            N[zone]
          ];
          refinedRecord = spectrumAt[refinedPoint];
          If[FailureQ[refinedRecord],
            Throw[
              Failure[
                "RefinementHamiltonianFailure",
                <|
                  "Point" -> refinedPoint,
                  "Detail" -> refinedRecord[[1]]
                |>
              ],
              "GaplessRefinementFailure"
            ]
          ];
          <|
            "Point" -> refinedPoint,
            "Gap" -> Re[
              refinedRecord["Eigenvalues"][[occupied + 1]] -
                refinedRecord["Eigenvalues"][[occupied]]
            ],
            "Source" -> seedRecord["Source"],
            "GridSeed" -> seedRecord["Point"],
            "GridGap" -> seedRecord["Gap"],
            "Refined" -> True,
            "RefinementStatus" -> "Succeeded",
            "ObjectiveMinimum" -> First[refinedResult]
          |>
        ]
      ],
      {seedRecord, seedRecords}
    ],
    "GaplessRefinementFailure"
  ];
  If[FailureQ[refinementRecords],
    Message[
      findGaplessPoints::refine,
      refinementRecords[[2, "Point"]],
      refinementRecords[[2, "Detail"]]
    ];
    Return[$Failed]
  ];
  failedRefinements = Count[
    Lookup[refinementRecords, "RefinementStatus", "Failed"],
    "Failed"
  ];
  qualifying = Select[
    SortBy[refinementRecords, Lookup[#, "Gap"] &],
    topologyFiniteNumberQ[Lookup[#, "Gap", Infinity]] &&
      TrueQ[Lookup[#, "Gap"] <= gapTolerance] &
  ];
  merged = Fold[
    Function[{accepted, candidate},
      If[
        AnyTrue[
          accepted,
          topologyPeriodicDistance[
            Lookup[#, "Point"],
            Lookup[candidate, "Point"],
            N[zone]
          ] <= mergeTolerance &
        ],
        accepted,
        Append[accepted, candidate]
      ]
    ],
    {},
    qualifying
  ];
  finalRecords = Map[
    Function[record,
      With[{spectrum = spectrumAt[record["Point"]]},
        Join[
          record,
          <|"Eigenvalues" -> spectrum["Eigenvalues"]|>
        ]
      ]
    ],
    merged
  ];
  minimumGap = Min[
    Join[
      gridGaps,
      Lookup[refinementRecords, "Gap", {}]
    ]
  ];
  data = <|
    "Schema" -> "MagneticTBGaplessPointSearch",
    "SchemaVersion" -> 1,
    "Points" -> Lookup[finalRecords, "Point", {}],
    "PointRecords" -> finalRecords,
    "OccupiedBands" -> occupied,
    "HamiltonianDimension" -> hamiltonianDimension,
    "BrillouinZone" -> zone,
    "GridSize" -> gridSize,
    "GridPointCount" -> Length[gridPoints],
    "LocalMinimumCount" -> Length[localMinimumIndices],
    "SeedCount" -> Length[seedRecords],
    "FailedRefinements" -> failedRefinements,
    "GapTolerance" -> gapTolerance,
    "MergeTolerance" -> mergeTolerance,
    "RefinementMethod" -> refinementMethod,
    "MinimumGap" -> minimumGap,
    "Method" ->
      "periodic grid local minima plus configurable gap-squared refinement",
    "CompletenessGuaranteed" -> False
  |>;
  If[output === "Data", data, data["Points"]]
];

End[]

EndPackage[]
