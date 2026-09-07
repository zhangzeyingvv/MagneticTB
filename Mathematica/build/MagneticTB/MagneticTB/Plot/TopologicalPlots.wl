(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  topologicalPlotFiniteRealQ,
  topologicalPlotPointQ,
  topologicalPlotSegmentQ,
  topologicalPlotBoundaryLabels,
  sampleTopologicalPlotPath,
  styledTopologicalXTicks,
  plotWilsonLoop,
  plotSurfaceSpectrum
];

topologicalPlotFiniteRealQ[value_] :=
  NumericQ[value] &&
    FreeQ[N[value], Indeterminate | ComplexInfinity | DirectedInfinity] &&
    TrueQ[Im[N[value]] == 0];

topologicalPlotPointQ[point_, dimension_Integer?Positive] :=
  ListQ[point] && Length[point] === dimension &&
    VectorQ[point, topologicalPlotFiniteRealQ];

topologicalPlotSegmentQ[segment_, dimension_Integer?Positive] :=
  ListQ[segment] && Length[segment] === 2 &&
    ListQ[segment[[1]]] && Length[segment[[1]]] === 2 &&
    And @@ (topologicalPlotPointQ[#, dimension] & /@ segment[[1]]) &&
    MatchQ[segment[[2]], {_String, _String}];

topologicalPlotBoundaryLabels[path_List] := Module[{labels},
  labels = {path[[1, 2, 1]]};
  Do[
    AppendTo[
      labels,
      If[
        index < Length[path] &&
          path[[index, 2, 2]] =!= path[[index + 1, 2, 1]],
        path[[index, 2, 2]] <> "|" <> path[[index + 1, 2, 1]],
        path[[index, 2, 2]]
      ]
    ],
    {index, Length[path]}
  ];
  StringReplace[#, "\\Gamma" -> "\[CapitalGamma]"] & /@ labels
];

sampleTopologicalPlotPath[
    path_, subdivisions_Integer?Positive, dimension_Integer?Positive
  ] := Module[
  {sampledSegments, xSegments, points, coordinates, distances},
  Which[
    ListQ[path] && Length[path] >= 2 &&
      And @@ (topologicalPlotPointQ[#, dimension] & /@ path),
      distances = Norm /@ Differences[N[path]];
      If[!And @@ (# > 0 & /@ distances), Return[$Failed]];
      <|
        "Points" -> N[path],
        "Coordinates" -> FoldList[Plus, 0., distances],
        "BoundaryCoordinates" -> None,
        "BoundaryLabels" -> None
      |>,

    ListQ[path] && path =!= {} &&
      And @@ (topologicalPlotSegmentQ[#, dimension] & /@ path),
      If[
        !And @@ MapThread[
          Norm[N[#1[[1, 2]] - #2[[1, 1]]]] <= 10^-10 &,
          {Most[path], Rest[path]}
        ] ||
          !And @@ (Norm[N[#[[1, 2]] - #[[1, 1]]]] > 0 & /@ path),
        Return[$Failed]
      ];
      sampledSegments = Subdivide[#[[1]], #[[2]], subdivisions] & /@
        N[path[[All, 1]]];
      xSegments = Table[
        Subdivide[index - 1., index, subdivisions],
        {index, Length[path]}
      ];
      points = Flatten[
        MapIndexed[If[First[#2] === 1, #1, Rest[#1]] &, sampledSegments],
        1
      ];
      coordinates = Flatten[
        MapIndexed[If[First[#2] === 1, #1, Rest[#1]] &, xSegments],
        1
      ];
      <|
        "Points" -> points,
        "Coordinates" -> coordinates,
        "BoundaryCoordinates" -> N@Range[0, Length[path]],
        "BoundaryLabels" -> topologicalPlotBoundaryLabels[path]
      |>,

    True,
      $Failed
  ]
];

styledTopologicalXTicks[data_Association, fontFamily_, fontSize_] := If[
  data["BoundaryCoordinates"] === None,
  Automatic,
  MapThread[
    {
      #1,
      Style[
        #2,
        Black,
        FontFamily -> fontFamily,
        FontSize -> fontSize
      ]
    } &,
    {data["BoundaryCoordinates"], data["BoundaryLabels"]}
  ]
];

plotWilsonLoop::path =
  "The parameter path must be either at least two distinct real momentum offsets or a connected list of labeled momentum segments with the same dimension as the Wilson loop.";
plotWilsonLoop::option = "Invalid plotWilsonLoop option value(s): `1`.";
plotWilsonLoop::compute =
  "Wilson-loop evaluation failed at parameter-path point `1`.";

Options[plotWilsonLoop] = {
  "ParameterSubdivisions" -> 60,
  "LoopSubdivisions" -> 50,
  "HermitianTolerance" -> 10^-10,
  "GapTolerance" -> 10^-9,
  "CovarianceTolerance" -> 10^-8,
  "OverlapTolerance" -> 10^-10,
  "PhaseConvention" -> "PhaseOverPi",
  "Output" -> "Graphics",
  MagneticTB`yTicks -> Automatic,
  Joined -> False,
  PlotMarkers -> Automatic,
  PlotStyle -> Black,
  PlotRange -> Automatic,
  GridLines -> Automatic,
  FrameLabel -> Automatic,
  FontSize -> 18,
  FontFamily -> "Times",
  ImageSize -> Large
};

plotWilsonLoop[
    h_, centers_, occupied_, start_, end_, parameterPath_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    suppliedOptionList = {suppliedOptions}, unknownOptionNames,
    parameterSubdivisions, loopSubdivisions, tolerances,
    phaseConvention, output, fontSize, fontFamily, pathData,
    loopOptions, values, failedPosition, datasets, xTicks, styledYTicks,
    resolvedPlotRange, gridLines, frameLabel, graphic, data
  },
  unknownOptionNames = Complement[
    First /@ suppliedOptionList,
    First /@ Options[plotWilsonLoop]
  ];
  If[unknownOptionNames =!= {},
    Message[plotWilsonLoop::option, unknownOptionNames];
    Return[$Failed]
  ];
  parameterSubdivisions = OptionValue["ParameterSubdivisions"];
  loopSubdivisions = OptionValue["LoopSubdivisions"];
  tolerances = (OptionValue[#] &) /@ {
    "HermitianTolerance", "GapTolerance", "CovarianceTolerance",
    "OverlapTolerance"
  };
  phaseConvention = OptionValue["PhaseConvention"];
  output = OptionValue["Output"];
  fontSize = OptionValue[FontSize];
  fontFamily = OptionValue[FontFamily];
  If[
    !IntegerQ[parameterSubdivisions] || parameterSubdivisions < 1 ||
      !IntegerQ[loopSubdivisions] || loopSubdivisions < 1 ||
      !And @@ (topologicalPlotFiniteRealQ[#] && # >= 0 & /@ tolerances) ||
      !MemberQ[{"PhaseOverPi", "WannierCenters"}, phaseConvention] ||
      !MemberQ[{"Graphics", "Data"}, output] ||
      !MemberQ[{True, False}, OptionValue[Joined]] ||
      !topologicalPlotFiniteRealQ[fontSize] || fontSize <= 0 ||
      !StringQ[fontFamily],
    Message[
      plotWilsonLoop::option,
      {
        parameterSubdivisions, loopSubdivisions, tolerances,
        phaseConvention, output, OptionValue[Joined], fontSize,
        fontFamily
      }
    ];
    Return[$Failed]
  ];
  If[
    !ListQ[start] || !ListQ[end] || start === {} ||
      Length[start] =!= Length[end],
    Message[plotWilsonLoop::path];
    Return[$Failed]
  ];
  pathData = sampleTopologicalPlotPath[
    parameterPath,
    parameterSubdivisions,
    Length[start]
  ];
  If[pathData === $Failed,
    Message[plotWilsonLoop::path];
    Return[$Failed]
  ];
  loopOptions = {
    "PathSubdivisions" -> loopSubdivisions,
    "HermitianTolerance" -> tolerances[[1]],
    "GapTolerance" -> tolerances[[2]],
    "CovarianceTolerance" -> tolerances[[3]],
    "OverlapTolerance" -> tolerances[[4]],
    "Output" -> phaseConvention
  };
  values = Map[
    wilsonLoop[
      h,
      centers,
      occupied,
      N[start] + #,
      N[end] + #,
      Sequence @@ loopOptions
    ] &,
    pathData["Points"]
  ];
  failedPosition = FirstPosition[values, $Failed];
  If[!MissingQ[failedPosition],
    Message[
      plotWilsonLoop::compute,
      pathData["Points"][[First[failedPosition]]]
    ];
    Return[$Failed]
  ];
  If[
    !MatrixQ[values, topologicalPlotFiniteRealQ] ||
      values === {} || Length[First[values]] < 1,
    Message[plotWilsonLoop::compute, "non-numeric result"];
    Return[$Failed]
  ];
  datasets = Table[
    Transpose[{pathData["Coordinates"], values[[All, branch]]}],
    {branch, Length[First[values]]}
  ];
  xTicks = styledTopologicalXTicks[pathData, fontFamily, fontSize];
  styledYTicks = styledBandYTicks[
    OptionValue[MagneticTB`yTicks],
    fontFamily,
    fontSize
  ];
  If[styledYTicks === $Failed,
    Message[plotWilsonLoop::option, OptionValue[MagneticTB`yTicks]];
    Return[$Failed]
  ];
  resolvedPlotRange = Replace[
    OptionValue[PlotRange],
    Automatic -> {
      All,
      If[phaseConvention === "PhaseOverPi", {-1, 1}, {0, 1}]
    }
  ];
  gridLines = Replace[
    OptionValue[GridLines],
    Automatic -> {
      Replace[pathData["BoundaryCoordinates"], None -> None],
      If[phaseConvention === "PhaseOverPi", {0}, {0, 1/2, 1}]
    }
  ];
  frameLabel = Replace[
    OptionValue[FrameLabel],
    Automatic -> {
      "parameter path",
      If[
        phaseConvention === "PhaseOverPi",
        "Wilson phase / \[Pi]",
        "Wannier center"
      ]
    }
  ];
  graphic = ListPlot[
    datasets,
    Joined -> OptionValue[Joined],
    PlotMarkers -> OptionValue[PlotMarkers],
    PlotStyle -> OptionValue[PlotStyle],
    PlotRange -> resolvedPlotRange,
    GridLines -> gridLines,
    Frame -> True,
    FrameStyle -> Black,
    FrameTicks -> {{styledYTicks, None}, {xTicks, None}},
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
    ImageSize -> OptionValue[ImageSize]
  ];
  data = <|
    "Schema" -> "MagneticTBWilsonLoopPlot",
    "SchemaVersion" -> 1,
    "ParameterPath" -> pathData["Points"],
    "PathCoordinate" -> pathData["Coordinates"],
    "BoundaryCoordinates" -> pathData["BoundaryCoordinates"],
    "BoundaryLabels" -> pathData["BoundaryLabels"],
    "PhaseConvention" -> phaseConvention,
    "Values" -> values,
    "BranchCount" -> Length[First[values]],
    "Graphics" -> graphic
  |>;
  If[output === "Data", data, graphic]
];

plotSurfaceSpectrum::path =
  "The surface momentum path must be either at least two distinct real two-vectors or a connected list of labeled two-dimensional momentum segments.";
plotSurfaceSpectrum::energy =
  "The energy range must be an ordered pair {emin,emax} of finite real numbers with emin<emax.";
plotSurfaceSpectrum::option =
  "Invalid plotSurfaceSpectrum option value(s): `1`.";
plotSurfaceSpectrum::compute =
  "Surface Green-function evaluation failed at momentum `1` and energy `2`.";

Options[plotSurfaceSpectrum] = {
  "MomentumSubdivisions" -> 60,
  "EnergyPoints" -> 201,
  "CellMatrix" -> Automatic,
  "Broadening" -> 10^-3,
  "Tolerance" -> 10^-10,
  "MaxIterations" -> 200,
  "Surface" -> "Positive",
  "HermiticityTolerance" -> 10^-9,
  "Output" -> "Graphics",
  ColorFunction -> "SolarColors",
  ColorFunctionScaling -> True,
  PlotLegends -> Automatic,
  PlotRange -> All,
  InterpolationOrder -> 1,
  Mesh -> None,
  FrameLabel -> Automatic,
  FontSize -> 18,
  FontFamily -> "Times",
  AspectRatio -> 1/GoldenRatio,
  ImageSize -> Large
};

plotSurfaceSpectrum[
    input_, momentumPath_, energyRange_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    suppliedOptionList = {suppliedOptions}, unknownOptionNames,
    momentumSubdivisions, energyPoints, output, fontSize, fontFamily,
    pathData, energies, surfaceOptions, failedMomentum = None,
    failedEnergy = None, weights, densityData, xTicks, frameLabel,
    graphic, data
  },
  unknownOptionNames = Complement[
    First /@ suppliedOptionList,
    First /@ Options[plotSurfaceSpectrum]
  ];
  If[unknownOptionNames =!= {},
    Message[plotSurfaceSpectrum::option, unknownOptionNames];
    Return[$Failed]
  ];
  momentumSubdivisions = OptionValue["MomentumSubdivisions"];
  energyPoints = OptionValue["EnergyPoints"];
  output = OptionValue["Output"];
  fontSize = OptionValue[FontSize];
  fontFamily = OptionValue[FontFamily];
  If[
    !IntegerQ[momentumSubdivisions] || momentumSubdivisions < 1 ||
      !IntegerQ[energyPoints] || energyPoints < 2 ||
      !MemberQ[{"Graphics", "Data"}, output] ||
      !topologicalPlotFiniteRealQ[fontSize] || fontSize <= 0 ||
      !StringQ[fontFamily] ||
      !MemberQ[{0, 1}, OptionValue[InterpolationOrder]],
    Message[
      plotSurfaceSpectrum::option,
      {
        momentumSubdivisions, energyPoints, output, fontSize,
        fontFamily, OptionValue[InterpolationOrder]
      }
    ];
    Return[$Failed]
  ];
  If[
    !ListQ[energyRange] || Length[energyRange] =!= 2 ||
      !And @@ (topologicalPlotFiniteRealQ /@ energyRange) ||
      !TrueQ[energyRange[[1]] < energyRange[[2]]],
    Message[plotSurfaceSpectrum::energy];
    Return[$Failed]
  ];
  pathData = sampleTopologicalPlotPath[
    momentumPath,
    momentumSubdivisions,
    2
  ];
  If[pathData === $Failed,
    Message[plotSurfaceSpectrum::path];
    Return[$Failed]
  ];
  energies = N@Subdivide[
    energyRange[[1]],
    energyRange[[2]],
    energyPoints - 1
  ];
  surfaceOptions = {
    "CellMatrix" -> OptionValue["CellMatrix"],
    "Broadening" -> OptionValue["Broadening"],
    "Tolerance" -> OptionValue["Tolerance"],
    "MaxIterations" -> OptionValue["MaxIterations"],
    "Surface" -> OptionValue["Surface"],
    "HermiticityTolerance" -> OptionValue["HermiticityTolerance"]
  };
  weights = Catch@Table[
    Module[{value},
      value = surfaceSpectralFunction[
        input,
        momentum,
        energy,
        Sequence @@ surfaceOptions
      ];
      If[value === $Failed || !topologicalPlotFiniteRealQ[value],
        failedMomentum = momentum;
        failedEnergy = energy;
        Throw[$Failed]
      ];
      value
    ],
    {momentum, pathData["Points"]},
    {energy, energies}
  ];
  If[weights === $Failed,
    Message[
      plotSurfaceSpectrum::compute,
      failedMomentum,
      failedEnergy
    ];
    Return[$Failed]
  ];
  densityData = Flatten[
    MapThread[
      Function[{coordinate, row},
        MapThread[{coordinate, #1, #2} &, {energies, row}]
      ],
      {pathData["Coordinates"], weights}
    ],
    1
  ];
  xTicks = styledTopologicalXTicks[pathData, fontFamily, fontSize];
  frameLabel = Replace[
    OptionValue[FrameLabel],
    Automatic -> {"surface k path", "energy"}
  ];
  graphic = ListDensityPlot[
    densityData,
    ColorFunction -> OptionValue[ColorFunction],
    ColorFunctionScaling -> OptionValue[ColorFunctionScaling],
    PlotLegends -> OptionValue[PlotLegends],
    PlotRange -> OptionValue[PlotRange],
    InterpolationOrder -> OptionValue[InterpolationOrder],
    Mesh -> OptionValue[Mesh],
    Frame -> True,
    FrameStyle -> Black,
    FrameTicks -> {{Automatic, None}, {xTicks, None}},
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
    "Schema" -> "MagneticTBSurfaceSpectrumPlot",
    "SchemaVersion" -> 1,
    "SurfaceMomentumPath" -> pathData["Points"],
    "PathCoordinate" -> pathData["Coordinates"],
    "BoundaryCoordinates" -> pathData["BoundaryCoordinates"],
    "BoundaryLabels" -> pathData["BoundaryLabels"],
    "Energies" -> energies,
    "SpectralWeight" -> weights,
    "Broadening" -> OptionValue["Broadening"],
    "Surface" -> OptionValue["Surface"],
    "Graphics" -> graphic
  |>;
  If[output === "Data", data, graphic]
];

End[]
EndPackage[]
