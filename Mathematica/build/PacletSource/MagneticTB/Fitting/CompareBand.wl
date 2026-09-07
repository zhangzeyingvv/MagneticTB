(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  validComparisonPathQ,
  comparisonPathLabels,
  comparisonReferencePlotData
];

validComparisonPathQ[path_] :=
  ListQ[path] && path =!= {} && And @@ Map[
    Function[segment,
      MatchQ[segment, {{_List, _List}, {_, _}}] &&
        And @@ (Length[#] === 3 && VectorQ[#, NumericQ] & /@
          segment[[1]])
    ],
    path
  ];

comparisonPathLabels[path_] := Module[{labels, boundaries},
  labels = Map[ToString, path[[All, 2]], {2}];
  boundaries = Join[
    {labels[[1, 1]]},
    Table[
      If[
        labels[[i, 2]] === labels[[i + 1, 1]],
        labels[[i, 2]],
        labels[[i, 2]] <> "|" <> labels[[i + 1, 1]]
      ],
      {i, Length[labels] - 1}
    ],
    {labels[[-1, 2]]}
  ];
  StringReplace[boundaries, "\\Gamma" -> "\[CapitalGamma]"]
];

comparisonReferencePlotData[segments_, npoint_Integer] := Flatten[
  MapIndexed[
    Function[{segment, index},
      With[
        {
          xvalues = Subdivide[
            npoint (index[[1]] - 1),
            npoint index[[1]],
            Length[segment] - 1
          ]
        },
        Map[
          MapThread[List, {xvalues, #}] &,
          Transpose[segment]
        ]
      ]
    ],
    segments
  ],
  1
];

compareBand::matrix =
  "The Hamiltonian must be a nonempty square matrix after display wrappers are removed.";
compareBand::path =
  "path must contain numeric three-component segment endpoints, and npoint must be a positive integer.";
compareBand::rules = "rules must be a list or Association of replacement rules.";
compareBand::data =
  "Reference bands must be valid numeric records with exactly npoint records per path segment.";
compareBand::nonnumeric =
  "The Hamiltonian remains nonnumeric after applying rules and momentum: `1`.";

Options[compareBand] = {plotRange -> All};

compareBand[pathstr_, npoint_, ham_, rules_, vaspband_, OptionsPattern[]] :=
  Module[
    {
      hamiltonian,
      parameterRules,
      parameterizedHamiltonian,
      segments,
      kSegments,
      matrices,
      invalid,
      bandSegments,
      data,
      referenceSegments,
      vaspplotdata,
      labels,
      xticks,
      yticks,
      maxAbsEnergy,
      font = "Times",
      tbplot
    },
    hamiltonian = normalizeFittingHamiltonian[ham];
    If[!validFittingHamiltonianQ[hamiltonian],
      Message[compareBand::matrix];
      Return[$Failed]
    ];
    If[!validComparisonPathQ[pathstr] || !IntegerQ[npoint] || npoint <= 0,
      Message[compareBand::path];
      Return[$Failed]
    ];
    parameterRules = fittingRulesList[rules];
    If[parameterRules === $Failed,
      Message[compareBand::rules];
      Return[$Failed]
    ];
    If[
      !validEigenvalueDataQ[vaspband] ||
        Length[vaspband] =!= Length[pathstr] npoint,
      Message[compareBand::data];
      Return[$Failed]
    ];

    segments = pathstr[[All, 1]];
    kSegments = N[2 Pi (Subdivide[#[[1]], #[[2]], npoint] & /@
      segments)];
    parameterizedHamiltonian = hamiltonian /. parameterRules;
    matrices = Map[
      Function[k,
        N[
          parameterizedHamiltonian /.
            Thread[{kx, ky, kz} -> k]
        ]
      ],
      kSegments,
      {2}
    ];
    invalid = SelectFirst[
      Flatten[matrices, 1],
      !MatrixQ[#, NumericQ] &,
      Missing[]
    ];
    If[!MissingQ[invalid],
      Message[compareBand::nonnumeric, invalid];
      Return[$Failed]
    ];
    bandSegments = Map[
      Transpose[fittingNumericEigenvalues /@ #] &,
      matrices
    ];
    data = Chop@Flatten[
      MapIndexed[
        {npoint (#2[[1]] - 1) + #2[[3]] - 1, #1} &,
        bandSegments,
        {3}
      ],
      1
    ];

    referenceSegments = Partition[vaspband[[All, 2]], npoint];
    vaspplotdata = comparisonReferencePlotData[
      referenceSegments,
      npoint
    ];

    labels = comparisonPathLabels[pathstr];
    xticks = Transpose@{
      npoint Range[0, Length[labels] - 1],
      Style[#, Black, FontFamily -> font, 24] & /@ labels
    };
    maxAbsEnergy = Max[Abs@Flatten[bandSegments]];
    yticks = {
      #,
      Style[Round[#, 0.01], Black, FontFamily -> font, 24]
    } & /@ Subdivide[-maxAbsEnergy, maxAbsEnergy, 4];

    tbplot = ListLinePlot[
      data,
      PlotRange -> {
        {0, npoint (Length[labels] - 1)},
        OptionValue[plotRange]
      },
      PlotStyle -> Purple,
      GridLines -> {npoint Range[Length[labels]], {0}},
      Frame -> {{True, True}, {True, True}},
      FrameStyle -> True,
      FrameTicks -> {{yticks, None}, {xticks, None}},
      FrameTicksStyle -> Directive[Black, 24],
      GridLinesStyle -> Directive[Black]
    ];
    Labeled[
      Show[
        tbplot,
        ListLinePlot[vaspplotdata, PlotStyle -> Lighter[Blue, 0.5]]
      ],
      {"Purple: TB, Blue: VASP"},
      {Top}
    ]
  ];

End[]

EndPackage[]
