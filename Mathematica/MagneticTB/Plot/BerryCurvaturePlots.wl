(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  berryPlotFiniteRealQ,
  berryPlotRangesQ,
  berryPlotGridSize,
  berryPlotCoordinateSteps,
  berryPlotColorFunction,
  plotBerryCurvature2D,
  plotBerryCurvature3D
];

berryPlotFiniteRealQ[value_] :=
  NumericQ[value] &&
    FreeQ[N[value], Indeterminate | ComplexInfinity | DirectedInfinity] &&
    TrueQ[Im[N[value]] == 0];

berryPlotRangesQ[ranges_, dimension_Integer?Positive] :=
  ListQ[ranges] && Length[ranges] === dimension &&
    And @@ Map[
      ListQ[#] && Length[#] === 2 &&
        And @@ (berryPlotFiniteRealQ /@ #) && TrueQ[#[[1]] < #[[2]]] &,
      ranges
    ];

berryPlotGridSize[setting_, dimension_Integer?Positive] := Which[
  IntegerQ[setting] && setting >= 2,
    ConstantArray[setting, dimension],
  ListQ[setting] && Length[setting] === dimension &&
    And @@ (IntegerQ[#] && # >= 2 & /@ setting),
    setting,
  True,
    $Failed
];

berryPlotCoordinateSteps[
    setting_, coordinateDimension_Integer?Positive,
    plottedDirections_List, plottedSpacings_List
  ] := Module[{steps},
  steps = Which[
    setting === Automatic,
      ConstantArray[Min[plottedSpacings]/2, coordinateDimension],
    berryPlotFiniteRealQ[setting] && setting > 0,
      ConstantArray[setting, coordinateDimension],
    ListQ[setting] && Length[setting] === coordinateDimension &&
      And @@ (berryPlotFiniteRealQ[#] && # > 0 & /@ setting),
      setting,
    True,
      Return[$Failed]
  ];
  If[setting === Automatic,
    Do[
      steps[[plottedDirections[[index]]]] = plottedSpacings[[index]]/2,
      {index, Length[plottedDirections]}
    ]
  ];
  steps
];

berryPlotColorFunction[Automatic] := ColorData["TemperatureMap"];
berryPlotColorFunction[name_String] := Quiet@Check[ColorData[name], $Failed];
berryPlotColorFunction[function_] := function;

plotBerryCurvature2D::input =
  "Wannier centers must be a nonempty list of equal-length real coordinate vectors with dimension at least two.";
plotBerryCurvature2D::range =
  "The plotting range must contain two ordered finite real intervals; received `1`.";
plotBerryCurvature2D::directions =
  "Directions must contain two distinct coordinate indices valid for dimension `1`; received `2`.";
plotBerryCurvature2D::fixed =
  "FixedCoordinates must be Automatic or one finite real coordinate per Hamiltonian momentum dimension; received `1`.";
plotBerryCurvature2D::option =
  "Invalid plotBerryCurvature2D option value(s): `1`.";
plotBerryCurvature2D::compute =
  "Berry-curvature evaluation failed at momentum `1`.";

Options[plotBerryCurvature2D] = {
  "Directions" -> {1, 2},
  "FixedCoordinates" -> Automatic,
  "GridSize" -> 31,
  "StepSize" -> Automatic,
  "HermitianTolerance" -> 10^-10,
  "GapTolerance" -> 10^-9,
  "CovarianceTolerance" -> 10^-8,
  "OverlapTolerance" -> 10^-10,
  "SymmetricColorRange" -> True,
  "Output" -> "Graphics",
  ColorFunction -> Automatic,
  ColorFunctionScaling -> True,
  PlotLegends -> Automatic,
  PlotRange -> Automatic,
  InterpolationOrder -> 1,
  Mesh -> None,
  FrameLabel -> Automatic,
  FontSize -> 18,
  FontFamily -> "Times",
  AspectRatio -> 1,
  ImageSize -> Large
};

plotBerryCurvature2D[
    hamiltonian_, centers_, occupied_, ranges_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, coordinateDimension,
    directions, fixedCoordinates, gridSize, gridValues, spacings,
    coordinateSteps, tolerances, output, symmetricColorRangeQ,
    fontSize, fontFamily, failedPoint = None, curvatures,
    densityData, maximumAbsolute, resolvedPlotRange, frameLabel,
    colorFunction, testColor, graphic, data
  },
  unknown = Complement[
    First /@ supplied,
    First /@ Options[plotBerryCurvature2D]
  ];
  If[unknown =!= {},
    Message[plotBerryCurvature2D::option, unknown];
    Return[$Failed]
  ];
  If[
    !ListQ[centers] || centers === {} || !ListQ[First[centers]] ||
      Length[First[centers]] < 2 ||
      !And @@ (ListQ[#] && Length[#] === Length[First[centers]] &&
        VectorQ[#, berryPlotFiniteRealQ] & /@ centers),
    Message[plotBerryCurvature2D::input];
    Return[$Failed]
  ];
  coordinateDimension = Length[First[centers]];
  If[!berryPlotRangesQ[ranges, 2],
    Message[plotBerryCurvature2D::range, ranges];
    Return[$Failed]
  ];
  directions = OptionValue["Directions"];
  If[
    !ListQ[directions] || Length[directions] =!= 2 ||
      !VectorQ[directions, IntegerQ] || !DuplicateFreeQ[directions] ||
      !And @@ (Between[#, {1, coordinateDimension}] & /@ directions),
    Message[
      plotBerryCurvature2D::directions,
      coordinateDimension,
      directions
    ];
    Return[$Failed]
  ];
  fixedCoordinates = Replace[
    OptionValue["FixedCoordinates"],
    Automatic -> ConstantArray[0., coordinateDimension]
  ];
  If[
    !ListQ[fixedCoordinates] ||
      Length[fixedCoordinates] =!= coordinateDimension ||
      !VectorQ[fixedCoordinates, berryPlotFiniteRealQ],
    Message[plotBerryCurvature2D::fixed, fixedCoordinates];
    Return[$Failed]
  ];
  gridSize = berryPlotGridSize[OptionValue["GridSize"], 2];
  If[gridSize === $Failed,
    Message[plotBerryCurvature2D::option, OptionValue["GridSize"]];
    Return[$Failed]
  ];
  gridValues = MapThread[
    Subdivide[#1[[1]], #1[[2]], #2 - 1] &,
    {N[ranges], gridSize}
  ];
  spacings = MapThread[(#1[[2]] - #1[[1]])/(#2 - 1) &, {N[ranges], gridSize}];
  coordinateSteps = berryPlotCoordinateSteps[
    OptionValue["StepSize"],
    coordinateDimension,
    directions,
    spacings
  ];
  tolerances = (OptionValue[#] &) /@ {
    "HermitianTolerance", "GapTolerance", "CovarianceTolerance",
    "OverlapTolerance"
  };
  output = OptionValue["Output"];
  symmetricColorRangeQ = OptionValue["SymmetricColorRange"];
  fontSize = OptionValue[FontSize];
  fontFamily = OptionValue[FontFamily];
  If[
    coordinateSteps === $Failed ||
      !And @@ (berryPlotFiniteRealQ[#] && # >= 0 & /@ tolerances) ||
      !MemberQ[{"Graphics", "Data"}, output] ||
      !MemberQ[{True, False}, symmetricColorRangeQ] ||
      !MemberQ[{0, 1}, OptionValue[InterpolationOrder]] ||
      !berryPlotFiniteRealQ[fontSize] || fontSize <= 0 ||
      !StringQ[fontFamily],
    Message[
      plotBerryCurvature2D::option,
      {
        OptionValue["StepSize"], tolerances, output,
        symmetricColorRangeQ, OptionValue[InterpolationOrder],
        fontSize, fontFamily
      }
    ];
    Return[$Failed]
  ];
  curvatures = Catch@Table[
    Module[{point, value},
      point = N[fixedCoordinates];
      point[[directions[[1]]]] = firstCoordinate;
      point[[directions[[2]]]] = secondCoordinate;
      value = berryCurvature[
        hamiltonian,
        centers,
        occupied,
        point,
        "Directions" -> directions,
        "StepSize" -> coordinateSteps[[directions]],
        "HermitianTolerance" -> tolerances[[1]],
        "GapTolerance" -> tolerances[[2]],
        "CovarianceTolerance" -> tolerances[[3]],
        "OverlapTolerance" -> tolerances[[4]]
      ];
      If[value === $Failed || !berryPlotFiniteRealQ[value],
        failedPoint = point;
        Throw[$Failed]
      ];
      value
    ],
    {firstCoordinate, gridValues[[1]]},
    {secondCoordinate, gridValues[[2]]}
  ];
  If[curvatures === $Failed,
    Message[plotBerryCurvature2D::compute, failedPoint];
    Return[$Failed]
  ];
  densityData = Flatten[
    Table[
      {
        gridValues[[1, firstIndex]],
        gridValues[[2, secondIndex]],
        curvatures[[firstIndex, secondIndex]]
      },
      {firstIndex, gridSize[[1]]},
      {secondIndex, gridSize[[2]]}
    ],
    1
  ];
  maximumAbsolute = Max[Abs[Flatten[N[curvatures]]]];
  resolvedPlotRange = Replace[
    OptionValue[PlotRange],
    Automatic -> If[
      symmetricColorRangeQ && maximumAbsolute > 0,
      {-maximumAbsolute, maximumAbsolute},
      All
    ]
  ];
  frameLabel = Replace[
    OptionValue[FrameLabel],
    Automatic -> ("k" <> ToString[#] & /@ directions)
  ];
  colorFunction = berryPlotColorFunction[OptionValue[ColorFunction]];
  testColor = Quiet@Check[colorFunction[0.5], $Failed];
  If[colorFunction === $Failed || !TrueQ[ColorQ[testColor]],
    Message[plotBerryCurvature2D::option, OptionValue[ColorFunction]];
    Return[$Failed]
  ];
  graphic = ListDensityPlot[
    densityData,
    ColorFunction -> colorFunction,
    ColorFunctionScaling -> OptionValue[ColorFunctionScaling],
    PlotLegends -> OptionValue[PlotLegends],
    PlotRange -> resolvedPlotRange,
    InterpolationOrder -> OptionValue[InterpolationOrder],
    Mesh -> OptionValue[Mesh],
    Frame -> True,
    FrameStyle -> Black,
    FrameLabel -> frameLabel,
    FrameTicksStyle -> Directive[
      Black,
      FontFamily -> fontFamily,
      FontSize -> fontSize
    ],
    BaseStyle -> Directive[
      FontFamily -> fontFamily,
      FontSize -> fontSize
    ],
    AspectRatio -> OptionValue[AspectRatio],
    ImageSize -> OptionValue[ImageSize]
  ];
  data = <|
    "Schema" -> "MagneticTBBerryCurvature2DPlot",
    "SchemaVersion" -> 1,
    "Directions" -> directions,
    "FixedCoordinates" -> fixedCoordinates,
    "GridValues" -> gridValues,
    "GridSize" -> gridSize,
    "StepSize" -> coordinateSteps[[directions]],
    "Curvature" -> curvatures,
    "MaximumAbsoluteCurvature" -> maximumAbsolute,
    "Graphics" -> graphic
  |>;
  If[output === "Data", data, graphic]
];

plotBerryCurvature3D::input =
  "Wannier centers must be a nonempty list of finite real three-vectors.";
plotBerryCurvature3D::range =
  "The plotting range must contain three ordered finite real intervals; received `1`.";
plotBerryCurvature3D::option =
  "Invalid plotBerryCurvature3D option value(s): `1`.";
plotBerryCurvature3D::compute =
  "Berry-curvature vector evaluation failed at momentum `1` for oriented component `2`.";

Options[plotBerryCurvature3D] = {
  "GridSize" -> 7,
  "StepSize" -> Automatic,
  "HermitianTolerance" -> 10^-10,
  "GapTolerance" -> 10^-9,
  "CovarianceTolerance" -> 10^-8,
  "OverlapTolerance" -> 10^-10,
  "VectorScale" -> Automatic,
  "MagnitudeThreshold" -> 0,
  "ArrowheadSize" -> 0.025,
  "ArrowThickness" -> 1.6,
  "Output" -> "Graphics",
  ColorFunction -> "SolarColors",
  Axes -> True,
  AxesLabel -> Automatic,
  FontSize -> 14,
  ViewPoint -> Automatic,
  PlotRange -> All,
  ImageSize -> Large
};

plotBerryCurvature3D[
    hamiltonian_, centers_, occupied_, ranges_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, gridSize, gridValues,
    spacings, coordinateSteps, tolerances, componentDirections,
    failedPoint = None, failedDirections = None, nestedRecords,
    records, vectors, magnitudes, maximumMagnitude, threshold,
    vectorScale, arrowheadSize, arrowThickness, output, colorFunction,
    visibleRecords, colors, displacements, axesLabel, graphic, data
  },
  unknown = Complement[
    First /@ supplied,
    First /@ Options[plotBerryCurvature3D]
  ];
  If[unknown =!= {},
    Message[plotBerryCurvature3D::option, unknown];
    Return[$Failed]
  ];
  If[
    !ListQ[centers] || centers === {} ||
      !And @@ (ListQ[#] && Length[#] === 3 &&
        VectorQ[#, berryPlotFiniteRealQ] & /@ centers),
    Message[plotBerryCurvature3D::input];
    Return[$Failed]
  ];
  If[!berryPlotRangesQ[ranges, 3],
    Message[plotBerryCurvature3D::range, ranges];
    Return[$Failed]
  ];
  gridSize = berryPlotGridSize[OptionValue["GridSize"], 3];
  If[gridSize === $Failed,
    Message[plotBerryCurvature3D::option, OptionValue["GridSize"]];
    Return[$Failed]
  ];
  gridValues = MapThread[
    Subdivide[#1[[1]], #1[[2]], #2 - 1] &,
    {N[ranges], gridSize}
  ];
  spacings = MapThread[(#1[[2]] - #1[[1]])/(#2 - 1) &, {N[ranges], gridSize}];
  coordinateSteps = berryPlotCoordinateSteps[
    OptionValue["StepSize"],
    3,
    {1, 2, 3},
    spacings
  ];
  tolerances = (OptionValue[#] &) /@ {
    "HermitianTolerance", "GapTolerance", "CovarianceTolerance",
    "OverlapTolerance"
  };
  threshold = OptionValue["MagnitudeThreshold"];
  vectorScale = OptionValue["VectorScale"];
  arrowheadSize = OptionValue["ArrowheadSize"];
  arrowThickness = OptionValue["ArrowThickness"];
  output = OptionValue["Output"];
  If[
    coordinateSteps === $Failed ||
      !And @@ (berryPlotFiniteRealQ[#] && # >= 0 & /@ tolerances) ||
      !berryPlotFiniteRealQ[threshold] || threshold < 0 ||
      !(vectorScale === Automatic ||
        (berryPlotFiniteRealQ[vectorScale] && vectorScale > 0)) ||
      !berryPlotFiniteRealQ[arrowheadSize] || arrowheadSize <= 0 ||
      !berryPlotFiniteRealQ[arrowThickness] || arrowThickness <= 0 ||
      !MemberQ[{"Graphics", "Data"}, output] ||
      !berryPlotFiniteRealQ[OptionValue[FontSize]] ||
      OptionValue[FontSize] <= 0,
    Message[
      plotBerryCurvature3D::option,
      {
        OptionValue["StepSize"], tolerances, threshold, vectorScale,
        arrowheadSize, arrowThickness, output, OptionValue[FontSize]
      }
    ];
    Return[$Failed]
  ];
  componentDirections = {{2, 3}, {3, 1}, {1, 2}};
  nestedRecords = Catch@Table[
    Module[{point, vector},
      point = {firstCoordinate, secondCoordinate, thirdCoordinate};
      vector = Table[
        Module[{value},
          value = berryCurvature[
            hamiltonian,
            centers,
            occupied,
            point,
            "Directions" -> directions,
            "StepSize" -> coordinateSteps[[directions]],
            "HermitianTolerance" -> tolerances[[1]],
            "GapTolerance" -> tolerances[[2]],
            "CovarianceTolerance" -> tolerances[[3]],
            "OverlapTolerance" -> tolerances[[4]]
          ];
          If[value === $Failed || !berryPlotFiniteRealQ[value],
            failedPoint = point;
            failedDirections = directions;
            Throw[$Failed]
          ];
          value
        ],
        {directions, componentDirections}
      ];
      <|
        "Point" -> point,
        "Vector" -> vector,
        "Magnitude" -> Norm[vector]
      |>
    ],
    {firstCoordinate, gridValues[[1]]},
    {secondCoordinate, gridValues[[2]]},
    {thirdCoordinate, gridValues[[3]]}
  ];
  If[nestedRecords === $Failed,
    Message[
      plotBerryCurvature3D::compute,
      failedPoint,
      failedDirections
    ];
    Return[$Failed]
  ];
  records = Flatten[nestedRecords, 2];
  vectors = ArrayReshape[
    Lookup[records, "Vector"],
    Append[gridSize, 3]
  ];
  magnitudes = ArrayReshape[Lookup[records, "Magnitude"], gridSize];
  maximumMagnitude = Max[Lookup[records, "Magnitude"]];
  vectorScale = Replace[vectorScale, Automatic -> 0.45 Min[spacings]];
  visibleRecords = Select[records, #1["Magnitude"] > threshold &];
  colorFunction = berryPlotColorFunction[OptionValue[ColorFunction]];
  If[colorFunction === $Failed,
    Message[plotBerryCurvature3D::option, OptionValue[ColorFunction]];
    Return[$Failed]
  ];
  colors = Quiet@Check[
    If[
      maximumMagnitude > 0,
      colorFunction[#1["Magnitude"]/maximumMagnitude] & /@ visibleRecords,
      {}
    ],
    $Failed
  ];
  If[
    colors === $Failed || !And @@ (TrueQ[ColorQ[#]] & /@ colors),
    Message[plotBerryCurvature3D::option, OptionValue[ColorFunction]];
    Return[$Failed]
  ];
  displacements = If[
    maximumMagnitude > 0,
    vectorScale Lookup[visibleRecords, "Vector"]/maximumMagnitude,
    {}
  ];
  axesLabel = Replace[
    OptionValue[AxesLabel],
    Automatic -> {"k1", "k2", "k3"}
  ];
  graphic = Graphics3D[
    {
      Arrowheads[arrowheadSize],
      AbsoluteThickness[arrowThickness],
      MapThread[
        {
          #1,
          Arrow[{#2["Point"], #2["Point"] + #3}]
        } &,
        {colors, visibleRecords, displacements}
      ]
    },
    Axes -> OptionValue[Axes],
    AxesLabel -> axesLabel,
    Boxed -> False,
    ViewProjection -> "Orthographic",
    ViewPoint -> OptionValue[ViewPoint],
    PlotRange -> OptionValue[PlotRange],
    SphericalRegion -> True,
    BaseStyle -> Directive[FontSize -> OptionValue[FontSize]],
    ImageSize -> OptionValue[ImageSize]
  ];
  data = <|
    "Schema" -> "MagneticTBBerryCurvature3DPlot",
    "SchemaVersion" -> 1,
    "ComponentDirections" -> componentDirections,
    "GridValues" -> gridValues,
    "GridSize" -> gridSize,
    "StepSize" -> coordinateSteps,
    "Vectors" -> vectors,
    "Magnitudes" -> magnitudes,
    "MaximumMagnitude" -> maximumMagnitude,
    "VisibleVectorCount" -> Length[visibleRecords],
    "VectorScale" -> vectorScale,
    "Graphics" -> graphic
  |>;
  If[output === "Data", data, graphic]
];

End[]
EndPackage[]
