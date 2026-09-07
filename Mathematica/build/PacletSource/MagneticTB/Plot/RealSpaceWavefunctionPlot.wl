(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  realSpaceWavefunctionFiniteNumberQ,
  realSpaceWavefunctionSiteRecords,
  plotRealSpaceWavefunction
];

realSpaceWavefunctionFiniteNumberQ[value_] :=
  NumericQ[value] &&
    FreeQ[N[value], Indeterminate | ComplexInfinity | DirectedInfinity];

realSpaceWavefunctionSiteRecords[
    basisRecords_List,
    amplitudes_List,
    aggregation_,
    tolerance_
  ] := Module[
  {positions, groups},
  positions = Lookup[basisRecords, "CartesianPosition", $Failed];
  If[
    !ListQ[positions] ||
      !And @@ (ListQ[#] && Length[#] === 3 &&
        VectorQ[#, realSpaceWavefunctionFiniteNumberQ] & /@ positions),
    Return[$Failed]
  ];
  groups = If[
    aggregation === "Atom",
    Gather[
      Range[Length[basisRecords]],
      Norm[N[positions[[#1]] - positions[[#2]]]] <= tolerance &
    ],
    List /@ Range[Length[basisRecords]]
  ];
  MapIndexed[
    Function[{basisIndices, sitePosition},
      Module[{dominantIndex, weight},
        dominantIndex = First@MaximalBy[
          basisIndices,
          Abs[amplitudes[[#]]] &
        ];
        weight = Total[Abs[amplitudes[[basisIndices]]]^2];
        <|
          "SiteIndex" -> First[sitePosition],
          "BasisIndices" -> basisIndices,
          "CellIndices" -> DeleteDuplicates@Lookup[
            basisRecords[[basisIndices]],
            "CellIndex"
          ],
          "Cells" -> DeleteDuplicates@Lookup[
            basisRecords[[basisIndices]],
            "Cell"
          ],
          "OrbitalIndices" -> Lookup[
            basisRecords[[basisIndices]],
            "OrbitalIndex"
          ],
          "CartesianPosition" -> Mean[positions[[basisIndices]]],
          "Weight" -> weight,
          "RepresentativeAmplitude" -> amplitudes[[dominantIndex]],
          "Phase" -> Arg[amplitudes[[dominantIndex]]]
        |>
      ]
    ],
    groups
  ]
];

plotRealSpaceWavefunction::data =
  "Expected Data output from buildRealSpaceHamiltonian with complete BasisRecords, Lattice, and WannierCenters position metadata.";
plotRealSpaceWavefunction::state =
  "The wavefunction must be a finite numeric vector of length `1`; received length `2`.";
plotRealSpaceWavefunction::zero =
  "The supplied wavefunction has zero norm and cannot be displayed.";
plotRealSpaceWavefunction::option =
  "Invalid plotRealSpaceWavefunction option value(s): `1`.";

Options[plotRealSpaceWavefunction] = {
  "Aggregation" -> "Atom",
  "Normalize" -> True,
  "WeightThreshold" -> 10^-8,
  "PositionTolerance" -> 10^-8,
  "PhaseColoring" -> False,
  "ShowAllSites" -> True,
  "ShowSiteLabels" -> False,
  "AtomRadius" -> Automatic,
  "WavefunctionScale" -> Automatic,
  "Output" -> "Graphics",
  ColorFunction -> Automatic,
  FontSize -> 12,
  ViewPoint -> Automatic,
  ImageSize -> Large
};

plotRealSpaceWavefunction[
    realSpaceData_, state_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, dimension, basisRecords,
    lattice, aggregation, normalizeQ, threshold, positionTolerance,
    phaseColoringQ, showAllSitesQ, showSiteLabelsQ, output,
    amplitudes, stateNorm, siteRecords, visibleRecords, maximumWeight,
    minimumLatticeLength, atomRadius, wavefunctionScale,
    colorSetting, colorFunction, waveColors, waveRadii, graphic, data
  },
  unknown = Complement[
    First /@ supplied,
    First /@ Options[plotRealSpaceWavefunction]
  ];
  If[unknown =!= {},
    Message[plotRealSpaceWavefunction::option, unknown];
    Return[$Failed]
  ];
  If[
    !AssociationQ[realSpaceData] ||
      Lookup[realSpaceData, "Schema", None] =!=
        "MagneticTBFiniteRealSpaceHamiltonian" ||
      !TrueQ[Lookup[realSpaceData, "PositionDataAvailable", False]],
    Message[plotRealSpaceWavefunction::data];
    Return[$Failed]
  ];
  dimension = Lookup[realSpaceData, "Dimension", $Failed];
  basisRecords = Lookup[realSpaceData, "BasisRecords", $Failed];
  lattice = Lookup[realSpaceData, "Lattice", $Failed];
  If[
    !IntegerQ[dimension] || dimension < 1 ||
      !ListQ[basisRecords] || Length[basisRecords] =!= dimension ||
      !MatrixQ[lattice, realSpaceWavefunctionFiniteNumberQ] ||
      Dimensions[lattice] =!= {3, 3},
    Message[plotRealSpaceWavefunction::data];
    Return[$Failed]
  ];
  If[
    !ListQ[state] || Length[state] =!= dimension ||
      !VectorQ[state, realSpaceWavefunctionFiniteNumberQ],
    Message[
      plotRealSpaceWavefunction::state,
      dimension,
      If[ListQ[state], Length[state], "not a list"]
    ];
    Return[$Failed]
  ];
  aggregation = OptionValue["Aggregation"];
  normalizeQ = OptionValue["Normalize"];
  threshold = OptionValue["WeightThreshold"];
  positionTolerance = OptionValue["PositionTolerance"];
  phaseColoringQ = OptionValue["PhaseColoring"];
  showAllSitesQ = OptionValue["ShowAllSites"];
  showSiteLabelsQ = OptionValue["ShowSiteLabels"];
  atomRadius = OptionValue["AtomRadius"];
  wavefunctionScale = OptionValue["WavefunctionScale"];
  output = OptionValue["Output"];
  colorSetting = OptionValue[ColorFunction];
  If[
    !MemberQ[{"Atom", "Orbital"}, aggregation] ||
      !MemberQ[{True, False}, normalizeQ] ||
      !realSpaceWavefunctionFiniteNumberQ[threshold] || threshold < 0 ||
      !realSpaceWavefunctionFiniteNumberQ[positionTolerance] ||
      positionTolerance <= 0 ||
      !MemberQ[{True, False}, phaseColoringQ] ||
      !MemberQ[{True, False}, showAllSitesQ] ||
      !MemberQ[{True, False}, showSiteLabelsQ] ||
      !MemberQ[{"Graphics", "Data"}, output] ||
      !(atomRadius === Automatic ||
        (realSpaceWavefunctionFiniteNumberQ[atomRadius] && atomRadius > 0)) ||
      !(wavefunctionScale === Automatic ||
        (realSpaceWavefunctionFiniteNumberQ[wavefunctionScale] &&
          wavefunctionScale > 0)) ||
      !realSpaceWavefunctionFiniteNumberQ[OptionValue[FontSize]] ||
      OptionValue[FontSize] <= 0,
    Message[
      plotRealSpaceWavefunction::option,
      {
        aggregation, normalizeQ, threshold, positionTolerance,
        phaseColoringQ, showAllSitesQ, showSiteLabelsQ, atomRadius,
        wavefunctionScale, output, OptionValue[FontSize]
      }
    ];
    Return[$Failed]
  ];
  amplitudes = N[state];
  stateNorm = Sqrt[Total[Abs[amplitudes]^2]];
  If[!realSpaceWavefunctionFiniteNumberQ[stateNorm] || stateNorm == 0,
    Message[plotRealSpaceWavefunction::zero];
    Return[$Failed]
  ];
  If[normalizeQ, amplitudes = amplitudes/stateNorm];
  siteRecords = realSpaceWavefunctionSiteRecords[
    basisRecords,
    amplitudes,
    aggregation,
    positionTolerance
  ];
  If[siteRecords === $Failed,
    Message[plotRealSpaceWavefunction::data];
    Return[$Failed]
  ];
  visibleRecords = Select[siteRecords, #1["Weight"] > threshold &];
  maximumWeight = Max[Lookup[siteRecords, "Weight"]];
  minimumLatticeLength = Min[Norm /@ N[lattice]];
  atomRadius = Replace[atomRadius, Automatic -> 0.055 minimumLatticeLength];
  wavefunctionScale = Replace[
    wavefunctionScale,
    Automatic -> 0.22 minimumLatticeLength
  ];
  colorFunction = Which[
    colorSetting === Automatic,
      ColorData["SolarColors"],
    StringQ[colorSetting],
      Quiet@Check[ColorData[colorSetting], $Failed],
    True,
      colorSetting
  ];
  waveColors = Quiet@Check[
    If[
      phaseColoringQ,
      Hue[Mod[#1["Phase"]/(2 Pi), 1], 0.85, 0.9] & /@ visibleRecords,
      colorFunction[#1["Weight"]/maximumWeight] & /@ visibleRecords
    ],
    $Failed
  ];
  If[
    waveColors === $Failed ||
      !And @@ (TrueQ[ColorQ[#]] & /@ waveColors),
    Message[plotRealSpaceWavefunction::option, colorSetting];
    Return[$Failed]
  ];
  waveRadii = wavefunctionScale (
    Lookup[visibleRecords, "Weight"]/maximumWeight
  )^(1/3);
  graphic = Graphics3D[
    {
      If[
        showAllSitesQ,
        {
          Directive[GrayLevel[0.72], Opacity[0.45]],
          Sphere[Lookup[siteRecords, "CartesianPosition"], atomRadius]
        },
        {}
      ],
      MapThread[
        {
          #1,
          Specularity[White, 25],
          Sphere[#2["CartesianPosition"], #3]
        } &,
        {waveColors, visibleRecords, waveRadii}
      ],
      If[
        showSiteLabelsQ,
        Map[
          Text[
            Style[
              #1["SiteIndex"],
              Black,
              FontSize -> OptionValue[FontSize]
            ],
            #1["CartesianPosition"] + {0, 0, 1.4 atomRadius}
          ] &,
          siteRecords
        ],
        {}
      ]
    },
    Axes -> False,
    Boxed -> False,
    ViewProjection -> "Orthographic",
    ViewPoint -> OptionValue[ViewPoint],
    Lighting -> "Neutral",
    PlotRange -> All,
    PlotRangePadding -> Scaled[0.08],
    SphericalRegion -> True,
    ImageSize -> OptionValue[ImageSize]
  ];
  data = <|
    "Schema" -> "MagneticTBRealSpaceWavefunctionPlot",
    "SchemaVersion" -> 1,
    "Aggregation" -> aggregation,
    "Normalized" -> normalizeQ,
    "InputNorm" -> stateNorm,
    "Amplitudes" -> amplitudes,
    "SiteRecords" -> siteRecords,
    "VisibleSiteCount" -> Length[visibleRecords],
    "MaximumWeight" -> maximumWeight,
    "PhaseColoring" -> phaseColoringQ,
    "Graphics" -> graphic
  |>;
  If[output === "Data", data, graphic]
];

End[]
EndPackage[]
