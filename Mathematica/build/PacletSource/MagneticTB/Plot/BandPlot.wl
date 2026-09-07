(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  validBandPathQ,
  resolveBandPath,
  bandPathBoundaryLabels,
  styledBandYTicks,
  bandEnergySegments,
  bandPlotData
];

bandplot::path =
  "The band path must contain one or more {{{k1},{k2}},{label1,label2}} segments, or be Automatic after init.";
bandplot::npoint =
  "The number of subdivisions per path segment must be a positive integer; received `1`.";
bandplot::yticks =
  "yTicks must be Automatic, None, a list of numeric positions, or a standard list of tick specifications.";

validBandPathQ[path_] := MatchQ[
  path,
  {
    {
      {{_, _, _}, {_, _, _}},
      {_String, _String}
    } ..
  }
];

resolveBandPath[Automatic] := Module[{path},
  path = standardKPath[];
  If[FailureQ[path] || path === $Failed, $Failed, path]
];
resolveBandPath[path_] := If[validBandPathQ[path], path, $Failed];

bandPathBoundaryLabels[path_List] := Module[{labels, segmentCount},
  segmentCount = Length[path];
  labels = ConstantArray["", segmentCount + 1];
  labels[[1]] = path[[1, 2, 1]];
  Do[
    labels[[index + 1]] = If[
      index < segmentCount &&
        path[[index, 2, 2]] =!= path[[index + 1, 2, 1]],
      path[[index, 2, 2]] <> "|" <> path[[index + 1, 2, 1]],
      path[[index, 2, 2]]
    ],
    {index, segmentCount}
  ];
  StringReplace[#, "\\Gamma" -> "\[CapitalGamma]"] & /@ labels
];

styledBandYTicks[Automatic, _, _] := Automatic;
styledBandYTicks[None, _, _] := None;
styledBandYTicks[ticks_List, fontFamily_, fontSize_] := Which[
  VectorQ[ticks, NumericQ],
    {
      #,
      Style[#, Black, FontFamily -> fontFamily, FontSize -> fontSize]
    } & /@ ticks,
  And @@ Map[
    ListQ[#] && Length[#] >= 2 && NumericQ[First[#]] &,
    ticks
  ],
    ReplacePart[
      #,
      2 -> Style[
        #[[2]],
        Black,
        FontFamily -> fontFamily,
        FontSize -> fontSize
      ]
    ] & /@ ticks,
  True,
    $Failed
];
styledBandYTicks[_, _, _] := $Failed;

bandEnergySegments[
    path_List,
    npoint_Integer?Positive,
    ham_,
    rules_List
  ] := Module[{segments},
  segments = Subdivide[#[[1]], #[[2]], npoint] & /@ path[[All, 1]];
  Map[
    Function[kpoints,
      Transpose@Table[
        Sort@Eigenvalues@Evaluate[
          N[ham] /. rules /. {
            kx -> 2 Pi kpoint[[1]],
            ky -> 2 Pi kpoint[[2]],
            kz -> 2 Pi kpoint[[3]]
          }
        ],
        {kpoint, kpoints}
      ]
    ],
    segments
  ]
];

bandPlotData[path_List, npoint_Integer?Positive, ham_, rules_List] := Module[
  {segmentEnergies},
  segmentEnergies = bandEnergySegments[path, npoint, ham, rules];
  Chop@Flatten[
    MapIndexed[
      Function[{segmentBands, segmentIndex},
        Map[
          MapIndexed[
            {
              npoint (First[segmentIndex] - 1) + First[#2] - 1,
              #1
            } &,
            #
          ] &,
          segmentBands
        ]
      ],
      segmentEnergies
    ],
    1
  ]
];

Options[bandplot] = {
  MagneticTB`plotRange -> All,
  MagneticTB`yTicks -> Automatic,
  FontSize -> 24,
  FontFamily -> "Times",
  ImageSize -> Automatic
};
Options[showband] = Options[bandplot];

bandplot[
    suppliedPath_,
    npoint_,
    ham_,
    rules_,
    OptionsPattern[]
  ] := Module[
  {
    path, labels, fontFamily, fontSize, xticks, yticks, data,
    segmentCount
  },
  If[!IntegerQ[npoint] || npoint < 1,
    Message[bandplot::npoint, npoint];
    Return[$Failed]
  ];
  If[!ListQ[rules], Message[bandplot::path]; Return[$Failed]];
  path = resolveBandPath[suppliedPath];
  If[path === $Failed,
    Message[bandplot::path];
    Return[$Failed]
  ];
  fontFamily = OptionValue[FontFamily];
  fontSize = OptionValue[FontSize];
  If[!StringQ[fontFamily] || !NumericQ[fontSize] || fontSize <= 0,
    Message[bandplot::yticks];
    Return[$Failed]
  ];
  yticks = styledBandYTicks[
    OptionValue[MagneticTB`yTicks],
    fontFamily,
    fontSize
  ];
  If[yticks === $Failed,
    Message[bandplot::yticks];
    Return[$Failed]
  ];

  segmentCount = Length[path];
  labels = bandPathBoundaryLabels[path];
  xticks = MapThread[
    {
      #1,
      Style[
        #2,
        Black,
        FontFamily -> fontFamily,
        FontSize -> fontSize
      ]
    } &,
    {npoint Range[0, segmentCount], labels}
  ];
  data = bandPlotData[path, npoint, ham, rules];

  ListLinePlot[
    data,
    PlotRange -> {
      {0, npoint segmentCount},
      OptionValue[MagneticTB`plotRange]
    },
    PlotStyle -> Black,
    GridLines -> {npoint Range[0, segmentCount], {0}},
    Frame -> True,
    FrameStyle -> Black,
    FrameTicks -> {{yticks, None}, {xticks, None}},
    FrameTicksStyle -> Directive[
      Black,
      FontFamily -> fontFamily,
      FontSize -> fontSize
    ],
    BaseStyle -> Directive[
      FontFamily -> fontFamily,
      FontSize -> fontSize
    ],
    GridLinesStyle -> Directive[Black],
    ImageSize -> OptionValue[ImageSize]
  ]
];

showband[
    npoint_,
    ham_,
    rules_,
    opts : OptionsPattern[]
  ] := bandplot[Automatic, npoint, ham, rules, opts];

showband[args___] := bandplot[args];

banddata[
    suppliedPath_,
    npoint_Integer?Positive,
    ham_,
    rules_List,
    save_
  ] := Module[{path, segmentEnergies, bands},
  path = resolveBandPath[suppliedPath];
  If[path === $Failed,
    Message[bandplot::path];
    Return[$Failed]
  ];
  segmentEnergies = bandEnergySegments[path, npoint, ham, rules];
  bands = Transpose[Flatten /@ Transpose[Chop@segmentEnergies]];
  Print@Export[save, bands, "CSV"]
];

End[]
EndPackage[]
