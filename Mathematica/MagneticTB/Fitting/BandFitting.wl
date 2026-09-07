(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  fittingBandSeries,
  fittingComparisonPlot,
  fittingRealNumericQ,
  fittingBandIndices,
  fittingKNeighborhoodSpec,
  fittingKPointSelectedQ,
  fittingLinearHamiltonianModel,
  fittingLevenbergMarquardt
];

fittingBandSeries[indices_List, bands_List] :=
  fittingBandSeries[
    indices,
    bands,
    ConstantArray[True, Dimensions[bands]]
  ];

fittingBandSeries[
    indices_List,
    bands_List,
    selectionMask_List
  ] :=
  Select[
    MapThread[
      Pick[MapThread[List, {indices, #1}], #2, True] &,
      {Transpose[bands], Transpose[selectionMask]}
    ],
    # =!= {} &
  ];

fittingComparisonPlot[
    indices_List,
    referenceBands_List,
    initialBands_List,
    fittedBands_List
  ] :=
  fittingComparisonPlot[
    indices,
    referenceBands,
    initialBands,
    fittedBands,
    ConstantArray[True, Dimensions[referenceBands]]
  ];

fittingComparisonPlot[
    indices_List,
    referenceBands_List,
    initialBands_List,
    fittedBands_List,
    selectionMask_List
  ] := Legended[
  Show[
    ListPlot[
      fittingBandSeries[indices, referenceBands, selectionMask],
      PlotStyle -> Directive[Blue, AbsolutePointSize[5]],
      PlotMarkers -> Automatic
    ],
    ListLinePlot[
      fittingBandSeries[indices, initialBands, selectionMask],
      PlotStyle -> Directive[Gray, Dashed, Thin]
    ],
    ListLinePlot[
      fittingBandSeries[indices, fittedBands, selectionMask],
      PlotStyle -> Directive[Red, Thick]
    ],
    Frame -> True,
    FrameLabel -> {"k-point index", "Energy"},
    PlotRange -> All,
    ImageSize -> Large
  ],
  Placed[
    LineLegend[
      {
        Directive[Blue, Thick],
        Directive[Gray, Dashed],
        Directive[Red, Thick]
      },
      {"VASP", "Initial TB", "Fitted TB"}
    ],
    Above
  ]
];

fittingRealNumericQ[value_] :=
  NumericQ[value] && TrueQ[Chop[Im[N[value]]] === 0];

fittingBandIndices[selection_, bandCount_Integer] := Module[
  {indices, arguments},
  indices = Which[
    selection === All,
      Range[bandCount],
    IntegerQ[selection],
      {selection},
    ListQ[selection],
      selection,
    Head[selection] === Span,
      arguments = List @@ selection;
      Which[
        MatchQ[arguments, {_Integer, _Integer}],
          Range @@ arguments,
        MatchQ[arguments, {_Integer, _Integer, _Integer}] &&
          Last[arguments] > 0,
          Range @@ arguments,
        True,
          $Failed
      ],
    True,
      $Failed
  ];
  If[
    indices === $Failed || indices === {} ||
      !VectorQ[indices, IntegerQ] ||
      !And @@ (Between[#, {1, bandCount}] & /@ indices) ||
      Length[DeleteDuplicates[indices]] =!= Length[indices] ||
      !And @@ Thread[Rest[indices] > Most[indices]],
    $Failed,
    indices
  ]
];

fittingKNeighborhoodSpec[All] := <|"Mode" -> "All"|>;
fittingKNeighborhoodSpec[input_] := Module[
  {allowedKeys, center, periodic, hasRadius, hasRange, radius, range},
  If[!AssociationQ[input], Return[$Failed]];
  allowedKeys = {"Center", "Radius", "Range", "Periodic"};
  If[Complement[Keys[input], allowedKeys] =!= {}, Return[$Failed]];
  If[!KeyExistsQ[input, "Center"], Return[$Failed]];
  center = Lookup[input, "Center"];
  periodic = Lookup[input, "Periodic", True];
  hasRadius = KeyExistsQ[input, "Radius"];
  hasRange = KeyExistsQ[input, "Range"];
  If[
    Length[center] =!= 3 || !VectorQ[center, fittingRealNumericQ] ||
      !MemberQ[{True, False}, periodic] ||
      Boole[hasRadius] + Boole[hasRange] =!= 1,
    Return[$Failed]
  ];
  If[hasRadius,
    radius = Lookup[input, "Radius"];
    If[!fittingRealNumericQ[radius] || !TrueQ[radius >= 0],
      Return[$Failed]
    ];
    <|
      "Mode" -> "Radius",
      "Center" -> N[center],
      "Radius" -> N[radius],
      "Periodic" -> periodic
    |>,
    range = Lookup[input, "Range"];
    If[
      Length[range] =!= 3 ||
        !VectorQ[range, fittingRealNumericQ] ||
        !And @@ (TrueQ[# >= 0] & /@ range),
      Return[$Failed]
    ];
    <|
      "Mode" -> "Range",
      "Center" -> N[center],
      "Range" -> N[range],
      "Periodic" -> periodic
    |>
  ]
];

fittingKPointSelectedQ[_, <|"Mode" -> "All"|>] := True;
fittingKPointSelectedQ[point_List, specification_Association] := Module[
  {delta},
  delta = N[point - specification["Center"]];
  If[TrueQ[specification["Periodic"]], delta -= Round[delta]];
  Switch[
    specification["Mode"],
    "Radius", TrueQ[Norm[delta] <= specification["Radius"]],
    "Range", And @@ Thread[Abs[delta] <= specification["Range"]],
    _, False
  ]
];

(* The residual pairing is positional: sorted model eigenvalue j is paired
   with reference energy j at the same selected k-point. *)

fittingLinearHamiltonianModel[
    hamiltonian_List,
    params_List,
    momenta_List
  ] := Module[
  {
    coefficientRuleMatrix,
    coefficientAssociationMatrix,
    zeroPowers,
    constantTemplate,
    coefficientTemplates,
    constantMatrices,
    coefficientMatrices
  },
  coefficientRuleMatrix = Quiet@Check[
    Map[CoefficientRules[#, params] &, hamiltonian, {2}],
    $Failed
  ];
  If[
    coefficientRuleMatrix === $Failed ||
      !And @@ Flatten@Map[
        AllTrue[First /@ #, Total[#] <= 1 &] &,
        coefficientRuleMatrix,
        {2}
      ],
    Return[Failure["NonlinearParameterModel", <||>]]
  ];
  coefficientAssociationMatrix =
    Map[Association, coefficientRuleMatrix, {2}];
  zeroPowers = ConstantArray[0, Length[params]];
  constantTemplate = Map[
    Lookup[#, Key[zeroPowers], 0] &,
    coefficientAssociationMatrix,
    {2}
  ];
  coefficientTemplates = Table[
    Map[
      Lookup[#, Key[UnitVector[Length[params], index]], 0] &,
      coefficientAssociationMatrix,
      {2}
    ],
    {index, Length[params]}
  ];
  constantMatrices = N@Map[
    Function[k,
      constantTemplate /. Thread[{kx, ky, kz} -> k]
    ],
    momenta
  ];
  coefficientMatrices = N@Map[
    Function[template,
      Map[
        Function[k,
          template /. Thread[{kx, ky, kz} -> k]
        ],
        momenta
      ]
    ],
    coefficientTemplates
  ];
  If[
    !And @@ (MatrixQ[#, NumericQ] & /@ constantMatrices) ||
      !And @@ Flatten[
        Map[MatrixQ[#, NumericQ] &, coefficientMatrices, {2}]
      ],
    Return[Failure["NonNumericLinearModel", <||>]]
  ];
  <|
    "ConstantMatrices" -> constantMatrices,
    "CoefficientMatrices" -> coefficientMatrices
  |>
];

fittingLevenbergMarquardt[
    residualFunction_,
    initialValues_List,
    maxIterations_Integer,
    timeConstraint_,
    finiteDifferenceStep_?NumericQ,
    initialDamping_?NumericQ,
    tolerance_?NumericQ
  ] := Module[
  {
    values = N[initialValues],
    residual,
    loss,
    damping = N[initialDamping],
    startTime = AbsoluteTime[],
    timedOutQ,
    iteration,
    stepSizes,
    unitVectors,
    columns,
    columnFailure,
    jacobian,
    normalMatrix,
    gradient,
    dampingScale,
    delta,
    candidateValues,
    candidateResidual,
    candidateLoss,
    improvement,
    converged = False
  },
  timedOutQ[] := timeConstraint =!= Infinity &&
    AbsoluteTime[] - startTime >= timeConstraint;
  residual = Quiet@Check[N@residualFunction[values], $Failed];
  If[!VectorQ[residual, NumericQ],
    Return[Failure["InvalidResidual", <||>]]
  ];
  loss = residual.residual;
  unitVectors = IdentityMatrix[Length[values]];

  For[iteration = 1, iteration <= maxIterations, iteration++,
    If[timedOutQ[],
      Return[Failure[
        "TimeConstraint",
        <|"Values" -> values, "Loss" -> loss,
          "Iterations" -> iteration - 1|>
      ]]
    ];
    stepSizes = finiteDifferenceStep (Max[1., Abs[#]] & /@ values);
    columns = Table[
      If[timedOutQ[],
        Return[Failure[
          "TimeConstraint",
          <|"Values" -> values, "Loss" -> loss,
            "Iterations" -> iteration - 1|>
        ]],
        Module[{plusResidual, minusResidual},
          plusResidual = Quiet@Check[
            N@residualFunction[
              values + stepSizes[[index]] unitVectors[[index]]
            ],
            $Failed
          ];
          minusResidual = Quiet@Check[
            N@residualFunction[
              values - stepSizes[[index]] unitVectors[[index]]
            ],
            $Failed
          ];
          If[
            !VectorQ[plusResidual, NumericQ] ||
              !VectorQ[minusResidual, NumericQ] ||
              Length[plusResidual] =!= Length[residual] ||
              Length[minusResidual] =!= Length[residual],
            Return[Failure["InvalidResidual", <||>], Module]
          ];
          (plusResidual - minusResidual)/(2 stepSizes[[index]])
        ]
      ],
      {index, Length[values]}
    ];
    columnFailure = SelectFirst[
      columns,
      FailureQ,
      Missing["NoFailure"]
    ];
    If[!MissingQ[columnFailure], Return[columnFailure]];
    jacobian = Transpose[columns];
    normalMatrix = Transpose[jacobian].jacobian;
    gradient = Transpose[jacobian].residual;
    If[Norm[gradient, Infinity] <= tolerance,
      converged = True;
      Break[]
    ];
    dampingScale = Max[1., #] & /@ Diagonal[normalMatrix];
    delta = Quiet@Check[
      LinearSolve[
        normalMatrix + damping DiagonalMatrix[dampingScale],
        -gradient
      ],
      $Failed
    ];
    If[!VectorQ[delta, NumericQ],
      damping *= 10.;
      Continue[]
    ];
    candidateValues = values + delta;
    candidateResidual = Quiet@Check[
      N@residualFunction[candidateValues],
      $Failed
    ];
    If[!VectorQ[candidateResidual, NumericQ],
      damping *= 10.;
      Continue[]
    ];
    candidateLoss = candidateResidual.candidateResidual;
    If[candidateLoss < loss,
      improvement = loss - candidateLoss;
      values = candidateValues;
      residual = candidateResidual;
      loss = candidateLoss;
      damping = Max[damping/3., 10.^-15];
      If[
        Norm[delta] <= tolerance (1. + Norm[values]) ||
          improvement <= tolerance Max[1., loss],
        converged = True;
        Break[]
      ],
      damping = Min[damping*10., 10.^15]
    ]
  ];

  <|
    "Values" -> values,
    "Loss" -> loss,
    "Iterations" -> Min[iteration, maxIterations],
    "Converged" -> converged,
    "FinalDamping" -> damping
  |>
];

bandManipulateEig::matrix =
  "The Hamiltonian must be a nonempty square matrix after display wrappers are removed.";
bandManipulateEig::data =
  "Reference data must be a nonempty list of {{kx,ky,kz},{energy1,...}} numeric records.";
bandManipulateEig::parameter =
  "Subscript and Indexed expressions are not supported as fitting parameters: `1`.";
bandManipulateEig::function =
  "The Hamiltonian contains unsupported unresolved function heads: `1`.";
bandManipulateEig::nonnumeric =
  "The Hamiltonian remains nonnumeric after assigning all parameters and momentum: `1`.";

bandManipulateEig[h_, eigdata_] := Module[
  {
    hamiltonian,
    unsupportedParameters,
    unsupportedFunctions,
    params,
    controls,
    controlSpecifications,
    momenta,
    referenceBands,
    plotRange
  },
  hamiltonian = normalizeFittingHamiltonian[h];
  If[!validFittingHamiltonianQ[hamiltonian],
    Message[bandManipulateEig::matrix];
    Return[Failure[
      "InvalidHamiltonian",
      <|"Head" -> Head[h], "Dimensions" -> Dimensions[hamiltonian]|>
    ]]
  ];
  If[!validEigenvalueDataQ[eigdata],
    Message[bandManipulateEig::data];
    Return[Failure["InvalidReferenceBands", <||>]]
  ];
  unsupportedParameters =
    unsupportedBandParameterExpressions[hamiltonian];
  If[unsupportedParameters =!= {},
    Message[bandManipulateEig::parameter, unsupportedParameters];
    Return[Failure[
      "UnsupportedFittingParameter",
      <|"Expressions" -> unsupportedParameters|>
    ]]
  ];
  unsupportedFunctions = unsupportedBandFunctionHeads[hamiltonian];
  If[unsupportedFunctions =!= {},
    Message[bandManipulateEig::function, unsupportedFunctions];
    Return[Failure[
      "NonNumericHamiltonian",
      <|"UnsupportedFunctionHeads" -> unsupportedFunctions|>
    ]]
  ];

  params = fittingParameterSymbols[hamiltonian];
  controls = Table[Unique["fitParameter$"], {Length[params]}];
  controlSpecifications = MapThread[
    {{#1, 0, #2}, -1, 1} &,
    {controls, params}
  ];
  momenta = Chop@N[2 Pi eigdata[[All, 1]]];
  referenceBands = N[eigdata[[All, 2]]];
  plotRange = {Min[#] - 0.2, Max[#] + 0.2} &@Flatten[referenceBands];

  Print["Number of params:", Length[params]];
  Print["params:", params];

  With[
    {
      parameterSymbols = params,
      dynamicSymbols = controls,
      specifications = controlSpecifications,
      kList = momenta,
      targetBands = referenceBands,
      displayRange = plotRange,
      matrix = hamiltonian
    },
    Manipulate[
      Module[{matrices, invalid, modelBands},
        matrices = Map[
          Function[k,
            N[
              matrix /.
                Thread[parameterSymbols -> dynamicSymbols] /.
                Thread[{kx, ky, kz} -> k]
            ]
          ],
          kList
        ];
        invalid = SelectFirst[
          matrices,
          !MatrixQ[#, NumericQ] &,
          Missing[]
        ];
        If[!MissingQ[invalid],
          Message[bandManipulateEig::nonnumeric, invalid];
          Return@Style[
            "Hamiltonian is not numeric; inspect the remaining symbols.",
            Red
          ]
        ];
        modelBands = fittingNumericEigenvalues /@ matrices;
        Show[
          ListPlot[
            Transpose[modelBands],
            PlotRange -> All,
            PlotStyle -> Black
          ],
          ListPlot[
            Transpose[targetBands],
            PlotStyle -> Red
          ],
          PlotRange -> displayRange
        ]
      ],
      Evaluate[Sequence @@ specifications],
      Button[
        "ExportData",
        Print[Thread[parameterSymbols -> dynamicSymbols]]
      ]
    ]
  ]
];

fittingTB::matrix =
  "The Hamiltonian must be a nonempty square matrix after display wrappers are removed.";
fittingTB::data =
  "Reference data must contain numeric three-component k points and exactly `1` energies per selected point.";
fittingTB::range =
  "kRange must be a nonempty list of valid 1-based reference-data indices.";
fittingTB::energywindow =
  "EnergyWindow must be All or an inclusive real numeric interval {emin,emax} with emin <= emax.";
fittingTB::kneighborhood =
  "KPointNeighborhood must be All or an Association with Center, exactly one of Radius or Range, and optional Periodic -> True or False.";
fittingTB::bands =
  "BandSelection must be All, a valid 1-based band index, a strictly increasing list of unique indices, or a positive-step Span.";
fittingTB::weights =
  "ResidualWeights must be Automatic or a nonnegative real numeric matrix with the same record-by-band shape as eigenvalueData.";
fittingTB::selection =
  "The combined k-point, band, energy-window, and positive-weight selection contains no residuals.";
fittingTB::insufficient =
  "The selection contains only `1` residuals for `2` fitting parameters.";
fittingTB::init =
  "initialRules must assign a real numeric initial value to every fitting parameter: `1`.";
fittingTB::parameter =
  "The Hamiltonian contains unsupported parameter expressions or unresolved function heads.";
fittingTB::nonlinear =
  "The Hamiltonian must be affine-linear in every fitting parameter; no slower symbolic fallback was used.";
fittingTB::fit =
  "The Levenberg-Marquardt fit failed for the supplied Hamiltonian and initial values.";
fittingTB::option =
  "MaxIterations must be a positive integer; TimeConstraint, FiniteDifferenceStep, InitialDamping, and FitTolerance must be positive (TimeConstraint may also be Infinity).";
fittingTB::timeout =
  "The fit exceeded the requested TimeConstraint of `1` seconds.";

Options[fittingTB] = {
  MaxIterations -> 100,
  TimeConstraint -> 60,
  "EnergyWindow" -> All,
  "KPointNeighborhood" -> All,
  "BandSelection" -> All,
  "ResidualWeights" -> Automatic,
  "FiniteDifferenceStep" -> 10^-5,
  "InitialDamping" -> 10^-3,
  "FitTolerance" -> 10^-8
};

fittingTB[h_, eigdata_, krange_, initparms_, OptionsPattern[]] := Module[
  {
    hamiltonian,
    dimensions,
    bandCount,
    params,
    initialRules,
    initialValues,
    energyWindowOption,
    energyWindow,
    kNeighborhoodOption,
    kNeighborhood,
    bandSelectionOption,
    bandIndices,
    weightOption,
    allWeights,
    candidateKPointCount,
    selectedKIndices,
    selectedData,
    targetBands,
    selectedWeights,
    bandMask,
    selectionMask,
    usedKPointMask,
    selectionMaskFlat,
    residualWeights,
    residualCount,
    selectionRestrictedQ,
    selectedBandIndicesByKPoint,
    selectionSummary,
    objectiveName,
    momenta,
    linearModel,
    constantMatrices,
    coefficientMatrices,
    numericHamiltonians,
    residual,
    initialResidual,
    initialLoss,
    initialBands,
    maxIterations,
    timeConstraint,
    finiteDifferenceStep,
    initialDamping,
    fitTolerance,
    fitResult,
    fittedRules,
    fittedBands,
    comparisonPlot
  },
  hamiltonian = normalizeFittingHamiltonian[h];
  If[!validFittingHamiltonianQ[hamiltonian],
    Message[fittingTB::matrix];
    Return[$Failed]
  ];
  dimensions = Dimensions[hamiltonian];
  bandCount = dimensions[[1]];
  If[!validEigenvalueDataQ[eigdata],
    Message[fittingTB::data, bandCount];
    Return[$Failed]
  ];
  If[!And @@ (Length[#] === bandCount & /@ eigdata[[All, 2]]),
    Message[fittingTB::data, bandCount];
    Return[$Failed]
  ];
  If[
    !ListQ[krange] || krange === {} ||
      !VectorQ[krange, IntegerQ] ||
      !And @@ (Between[#, {1, Length[eigdata]}] & /@ krange),
    Message[fittingTB::range];
    Return[$Failed]
  ];

  If[
    unsupportedBandParameterExpressions[hamiltonian] =!= {} ||
      unsupportedBandFunctionHeads[hamiltonian] =!= {},
    Message[fittingTB::parameter];
    Return[$Failed]
  ];

  params = fittingParameterSymbols[hamiltonian];
  initialRules = fittingRulesList[initparms];
  If[initialRules === $Failed,
    Message[fittingTB::init, params];
    Return[$Failed]
  ];
  initialValues = params /. initialRules;
  If[
    !VectorQ[
      initialValues,
      NumericQ[#] && TrueQ[Chop[Im[N[#]]] === 0] &
    ],
    Message[fittingTB::init, params];
    Return[$Failed]
  ];
  maxIterations = OptionValue[MaxIterations];
  timeConstraint = OptionValue[TimeConstraint];
  finiteDifferenceStep = OptionValue["FiniteDifferenceStep"];
  initialDamping = OptionValue["InitialDamping"];
  fitTolerance = OptionValue["FitTolerance"];
  If[
    !IntegerQ[maxIterations] || maxIterations <= 0 ||
      !(timeConstraint === Infinity ||
        (NumericQ[timeConstraint] && TrueQ[timeConstraint > 0])) ||
      !NumericQ[finiteDifferenceStep] ||
      !TrueQ[finiteDifferenceStep > 0] ||
      !NumericQ[initialDamping] || !TrueQ[initialDamping > 0] ||
      !NumericQ[fitTolerance] || !TrueQ[fitTolerance > 0],
    Message[fittingTB::option];
    Return[$Failed]
  ];

  energyWindowOption = OptionValue["EnergyWindow"];
  energyWindow = Which[
    energyWindowOption === All,
      All,
    ListQ[energyWindowOption] && Length[energyWindowOption] === 2 &&
      VectorQ[energyWindowOption, fittingRealNumericQ] &&
      TrueQ[energyWindowOption[[1]] <= energyWindowOption[[2]]],
      N[energyWindowOption],
    True,
      $Failed
  ];
  If[energyWindow === $Failed,
    Message[fittingTB::energywindow];
    Return[$Failed]
  ];

  kNeighborhoodOption = OptionValue["KPointNeighborhood"];
  kNeighborhood = fittingKNeighborhoodSpec[kNeighborhoodOption];
  If[kNeighborhood === $Failed,
    Message[fittingTB::kneighborhood];
    Return[$Failed]
  ];

  bandSelectionOption = OptionValue["BandSelection"];
  bandIndices = fittingBandIndices[bandSelectionOption, bandCount];
  If[bandIndices === $Failed,
    Message[fittingTB::bands];
    Return[$Failed]
  ];

  weightOption = OptionValue["ResidualWeights"];
  allWeights = Which[
    weightOption === Automatic,
      ConstantArray[1., {Length[eigdata], bandCount}],
    Dimensions[weightOption] === {Length[eigdata], bandCount} &&
      MatrixQ[
        weightOption,
        fittingRealNumericQ[#] && TrueQ[# >= 0] &
      ],
      N[weightOption],
    True,
      $Failed
  ];
  If[allWeights === $Failed,
    Message[fittingTB::weights];
    Return[$Failed]
  ];

  selectedKIndices = Select[
    krange,
    fittingKPointSelectedQ[eigdata[[#, 1]], kNeighborhood] &
  ];
  If[selectedKIndices === {},
    Message[fittingTB::selection];
    Return[$Failed]
  ];
  candidateKPointCount = Length[selectedKIndices];
  selectedData = eigdata[[selectedKIndices]];
  targetBands = N[selectedData[[All, 2]]];
  selectedWeights = allWeights[[selectedKIndices]];
  bandMask = MemberQ[bandIndices, #] & /@ Range[bandCount];
  selectionMask = MapThread[
    Function[{energies, weights},
      MapThread[
        Function[{bandAllowed, energy, weight},
          TrueQ[
            bandAllowed && weight > 0 &&
              (energyWindow === All || Between[energy, energyWindow])
          ]
        ],
        {bandMask, energies, weights}
      ]
    ],
    {targetBands, selectedWeights}
  ];
  usedKPointMask = Or @@ # & /@ selectionMask;
  If[!MemberQ[usedKPointMask, True],
    Message[fittingTB::selection];
    Return[$Failed]
  ];
  selectedKIndices = Pick[selectedKIndices, usedKPointMask, True];
  selectedData = Pick[selectedData, usedKPointMask, True];
  targetBands = Pick[targetBands, usedKPointMask, True];
  selectedWeights = Pick[selectedWeights, usedKPointMask, True];
  selectionMask = Pick[selectionMask, usedKPointMask, True];
  selectionMaskFlat = Flatten[selectionMask];
  residualWeights = Pick[
    Flatten[selectedWeights],
    selectionMaskFlat,
    True
  ];
  residualCount = Length[residualWeights];
  If[residualCount === 0,
    Message[fittingTB::selection];
    Return[$Failed]
  ];
  selectionRestrictedQ = !(
    energyWindowOption === All &&
      kNeighborhoodOption === All &&
      bandSelectionOption === All &&
      weightOption === Automatic
  );
  If[selectionRestrictedQ && residualCount < Length[params],
    Message[fittingTB::insufficient, residualCount, Length[params]];
    Return[$Failed]
  ];
  selectedBandIndicesByKPoint =
    Pick[Range[bandCount], #, True] & /@ selectionMask;
  objectiveName = If[
    weightOption === Automatic,
    "SumSquaredBandResiduals",
    "WeightedSumSquaredBandResiduals"
  ];
  selectionSummary = <|
    "CandidateKPointCount" -> candidateKPointCount,
    "UsedKPointCount" -> Length[selectedKIndices],
    "UniqueKPointCount" -> Length[DeleteDuplicates[selectedKIndices]],
    "UsedKPointIndices" -> selectedKIndices,
    "ResidualCount" -> residualCount,
    "SelectionRestricted" -> selectionRestrictedQ,
    "SelectedBandIndicesByKPoint" -> selectedBandIndicesByKPoint,
    "BandSelection" -> bandIndices,
    "EnergyWindow" -> energyWindow,
    "KPointNeighborhood" -> kNeighborhood,
    "Weighting" -> If[
      weightOption === Automatic,
      "Uniform",
      "ExplicitMatrix"
    ],
    "BandPairing" ->
      "SortedModelEigenvalueToReferencePosition"
  |>;

  momenta = N[2 Pi selectedData[[All, 1]]];

  If[params === {},
    constantMatrices = N@Map[
      Function[k,
        hamiltonian /. Thread[{kx, ky, kz} -> k]
      ],
      momenta
    ];
    If[!And @@ (MatrixQ[#, NumericQ] & /@ constantMatrices),
      Message[fittingTB::parameter];
      Return[$Failed]
    ];
    initialBands = fittingNumericEigenvalues /@ constantMatrices;
    initialResidual = Sqrt[residualWeights] Pick[
      Flatten[initialBands - targetBands],
      selectionMaskFlat,
      True
    ];
    initialLoss = initialResidual.initialResidual;
    comparisonPlot = fittingComparisonPlot[
      selectedKIndices,
      targetBands,
      initialBands,
      initialBands,
      selectionMask
    ];
    Return[Join[<|
      "FittedParams" -> {},
      "Objective" -> objectiveName,
      "Optimizer" -> "None",
      "ParameterModel" -> "Constant",
      "Converged" -> True,
      "Iterations" -> 0,
      "LSQForInitParams" -> initialLoss,
      "LSQForFittedParams" -> initialLoss,
      "ComparisonPlot" -> comparisonPlot
    |>, selectionSummary]]
  ];

  linearModel = fittingLinearHamiltonianModel[
    hamiltonian,
    params,
    momenta
  ];
  If[MatchQ[linearModel, Failure["NonlinearParameterModel", _Association]],
    Message[fittingTB::nonlinear];
    Return[$Failed]
  ];
  If[!AssociationQ[linearModel],
    Message[fittingTB::parameter];
    Return[$Failed]
  ];
  constantMatrices = linearModel["ConstantMatrices"];
  coefficientMatrices = linearModel["CoefficientMatrices"];
  numericHamiltonians[values_List] /; VectorQ[values, NumericQ] :=
    constantMatrices + Total[
      MapThread[Times, {N[values], coefficientMatrices}]
    ];

  residual[values_List] /; VectorQ[values, NumericQ] := Module[
    {modelBands},
    modelBands = fittingNumericEigenvalues /@ numericHamiltonians[values];
    Sqrt[residualWeights] Pick[
      Flatten[modelBands - targetBands],
      selectionMaskFlat,
      True
    ]
  ];

  initialBands = fittingNumericEigenvalues /@
    numericHamiltonians[N[initialValues]];
  initialResidual = residual[N[initialValues]];
  initialLoss = initialResidual.initialResidual;
  fitResult = fittingLevenbergMarquardt[
    Function[values, residual[values]],
    N[initialValues],
    maxIterations,
    timeConstraint,
    finiteDifferenceStep,
    initialDamping,
    fitTolerance
  ];
  If[MatchQ[fitResult, Failure["TimeConstraint", _Association]],
    Message[fittingTB::timeout, timeConstraint];
    Return[$Failed]
  ];
  If[!AssociationQ[fitResult],
    Message[fittingTB::fit];
    Return[$Failed]
  ];
  fittedRules = Thread[
    params -> fitResult["Values"]
  ];
  fittedBands = fittingNumericEigenvalues /@
    numericHamiltonians[fitResult["Values"]];
  comparisonPlot = fittingComparisonPlot[
    selectedKIndices,
    targetBands,
    initialBands,
    fittedBands,
    selectionMask
  ];

  Join[<|
    "FittedParams" -> fittedRules,
    "Objective" -> objectiveName,
    "Optimizer" -> "LevenbergMarquardt",
    "ParameterModel" -> "CachedAffineCoefficientMatrices",
    "Converged" -> fitResult["Converged"],
    "Iterations" -> fitResult["Iterations"],
    "FinalDamping" -> fitResult["FinalDamping"],
    "LSQForInitParams" -> initialLoss,
    "LSQForFittedParams" -> fitResult["Loss"],
    "ComparisonPlot" -> comparisonPlot
  |>, selectionSummary]
];

End[]

EndPackage[]
