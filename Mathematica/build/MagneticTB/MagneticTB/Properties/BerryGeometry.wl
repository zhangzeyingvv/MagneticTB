(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  berryOptionNames,
  berryUnknownOptions,
  berryPhase,
  berryCurvature
];

berryOptionNames[function_] := First /@ Options[function];
berryUnknownOptions[function_, supplied_List] := Complement[
  DeleteDuplicates[
    First /@ Cases[supplied, _Rule | _RuleDelayed, {1}]
  ],
  berryOptionNames[function]
];

berryPhase::path =
  "The path must contain at least two finite real momentum vectors of equal positive dimension, including an endpoint that differs from the start by 2 Pi times an integer reciprocal vector.";
berryPhase::centers =
  "Wannier centers must contain one finite real coordinate vector per Hamiltonian row, with the same coordinate dimension as the path.";
berryPhase::occ =
  "The occupied-band count must be an integer between 1 and the Hamiltonian dimension; received `1`.";
berryPhase::ham =
  "The Hamiltonian at path point `1` must be a finite numeric square Hermitian matrix of constant dimension.";
berryPhase::gap =
  "The occupied subspace is not isolated at path point `1`; the direct gap `2` does not exceed tolerance `3`.";
berryPhase::covariance =
  "The endpoint Hamiltonians do not obey h(k+G)=V(G)^(-1).h(k).V(G); residual `1` exceeds tolerance `2`.";
berryPhase::overlap =
  "Occupied subspaces at adjacent or closure path positions have a singular overlap (minimum singular value `1` at link `2`).";
berryPhase::option = "Invalid berryPhase option value(s): `1`.";

Options[berryPhase] = {
  "HermitianTolerance" -> 10^-10,
  "GapTolerance" -> 10^-9,
  "CovarianceTolerance" -> 10^-8,
  "OverlapTolerance" -> 10^-10,
  "Output" -> "Phase"
};

berryPhase[
    hamiltonian_, centers_, occupied_, path_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, hermitianTolerance,
    gapTolerance, covarianceTolerance, overlapTolerance, output,
    coordinateDimension, closureVector, reciprocalCoordinates,
    matrices, dimension, hermitianResiduals, endpointSewing,
    covarianceResidual, spectraAndFrames, spectra, frames, gaps,
    minimumGap, overlaps, closureOverlap, linkData, unitaryLinks,
    minimumSingularValue, wilsonMatrix, unitarityResidual,
    wilsonDeterminant, phase, data
  },
  unknown = berryUnknownOptions[berryPhase, supplied];
  If[unknown =!= {},
    Message[berryPhase::option, unknown];
    Return[$Failed]
  ];
  hermitianTolerance = OptionValue["HermitianTolerance"];
  gapTolerance = OptionValue["GapTolerance"];
  covarianceTolerance = OptionValue["CovarianceTolerance"];
  overlapTolerance = OptionValue["OverlapTolerance"];
  output = OptionValue["Output"];
  If[
    !And @@ (
      wilsonLoopFiniteNumberQ[#] && TrueQ[Im[N[#]] == 0] &&
        TrueQ[Re[N[#]] >= 0] & /@ {
        hermitianTolerance, gapTolerance, covarianceTolerance,
        overlapTolerance
      }
    ) ||
      !MemberQ[
        {"Phase", "PhaseOverPi", "WilsonDeterminant", "Data"},
        output
      ],
    Message[berryPhase::option, {
      hermitianTolerance, gapTolerance, covarianceTolerance,
      overlapTolerance, output
    }];
    Return[$Failed]
  ];
  If[
    !ListQ[path] || Length[path] < 2 || !ListQ[First[path]] ||
      Length[First[path]] < 1,
    Message[berryPhase::path];
    Return[$Failed]
  ];
  coordinateDimension = Length[First[path]];
  If[
    !And @@ (
      ListQ[#] && Length[#] === coordinateDimension &&
        wilsonLoopRealVectorQ[#] & /@ path
    ),
    Message[berryPhase::path];
    Return[$Failed]
  ];
  closureVector = N[Last[path] - First[path]];
  reciprocalCoordinates = closureVector/(2 Pi);
  If[
    Max[Abs[
      reciprocalCoordinates - Round[reciprocalCoordinates]
    ]] > covarianceTolerance,
    Message[berryPhase::path];
    Return[$Failed]
  ];
  If[
    !ListQ[centers] || centers === {} ||
      !And @@ (
        ListQ[#] && Length[#] === coordinateDimension &&
          wilsonLoopRealVectorQ[#] & /@ centers
      ),
    Message[berryPhase::centers];
    Return[$Failed]
  ];
  matrices = wilsonLoopNumericMatrix[hamiltonian, #] & /@ path;
  If[MemberQ[matrices, $Failed],
    Message[
      berryPhase::ham,
      path[[First@FirstPosition[matrices, $Failed]]]
    ];
    Return[$Failed]
  ];
  dimension = Length[First[matrices]];
  If[Length[centers] =!= dimension,
    Message[berryPhase::centers];
    Return[$Failed]
  ];
  If[!IntegerQ[occupied] || !Between[occupied, {1, dimension}],
    Message[berryPhase::occ, occupied];
    Return[$Failed]
  ];
  If[!And @@ (Dimensions[#] === {dimension, dimension} & /@ matrices),
    Message[
      berryPhase::ham,
      path[[First@FirstPosition[
        Dimensions /@ matrices,
        dimensions_ /; dimensions =!= {dimension, dimension}
      ]]]
    ];
    Return[$Failed]
  ];
  hermitianResiduals =
    Max[Abs[Flatten[# - ConjugateTranspose[#]]]] & /@ matrices;
  If[Max[hermitianResiduals] > hermitianTolerance,
    Message[
      berryPhase::ham,
      path[[First@Ordering[hermitianResiduals, -1]]]
    ];
    Return[$Failed]
  ];
  endpointSewing = DiagonalMatrix[
    Exp[I closureVector.#] & /@ N[centers]
  ];
  covarianceResidual = Max@Abs@Flatten[
    Last[matrices] -
      ConjugateTranspose[endpointSewing] . First[matrices] .
        endpointSewing
  ];
  If[covarianceResidual > covarianceTolerance,
    Message[
      berryPhase::covariance,
      covarianceResidual,
      covarianceTolerance
    ];
    Return[$Failed]
  ];
  spectraAndFrames = wilsonLoopSortedFrame[#, occupied] & /@ matrices;
  If[MemberQ[spectraAndFrames, $Failed],
    Message[
      berryPhase::ham,
      path[[First@FirstPosition[spectraAndFrames, $Failed]]]
    ];
    Return[$Failed]
  ];
  spectra = spectraAndFrames[[All, 1]];
  frames = spectraAndFrames[[All, 2]];
  gaps = If[
    occupied < dimension,
    Re[#[[occupied + 1]] - #[[occupied]]] & /@ spectra,
    ConstantArray[Infinity, Length[path]]
  ];
  minimumGap = Min[gaps];
  If[minimumGap <= gapTolerance,
    Message[
      berryPhase::gap,
      path[[First@FirstPosition[gaps, minimumGap]]],
      minimumGap,
      gapTolerance
    ];
    Return[$Failed]
  ];
  overlaps = MapThread[
    Conjugate[#2] . Transpose[#1] &,
    {Most[frames], Rest[frames]}
  ];
  closureOverlap =
    Conjugate[First[frames]] . endpointSewing . Transpose[Last[frames]];
  linkData = Map[
    wilsonLoopUnitaryPart[#, overlapTolerance] &,
    Append[overlaps, closureOverlap]
  ];
  If[MemberQ[linkData, $Failed],
    Message[berryPhase::ham, "SVD"];
    Return[$Failed]
  ];
  If[MemberQ[linkData, {"Singular", _}],
    With[{position = First@FirstPosition[linkData, {"Singular", _}]},
      Message[
        berryPhase::overlap,
        linkData[[position, 2]],
        If[position <= Length[overlaps], position, "closure"]
      ]
    ];
    Return[$Failed]
  ];
  unitaryLinks = linkData[[All, 1]];
  minimumSingularValue = Min[linkData[[All, 2]]];
  wilsonMatrix = Last[unitaryLinks] .
    (Dot @@ Reverse[Most[unitaryLinks]]);
  wilsonMatrix = Chop[wilsonMatrix, hermitianTolerance];
  unitarityResidual = Max@Abs@Flatten[
    ConjugateTranspose[wilsonMatrix] . wilsonMatrix -
      IdentityMatrix[occupied]
  ];
  wilsonDeterminant = Det[wilsonMatrix];
  phase = Chop[Arg[wilsonDeterminant], hermitianTolerance];
  data = <|
    "Phase" -> phase,
    "PhaseOverPi" -> phase/Pi,
    "WilsonDeterminant" -> wilsonDeterminant,
    "WilsonMatrix" -> wilsonMatrix,
    "UnitarityResidual" -> unitarityResidual,
    "MinimumLinkSingularValue" -> minimumSingularValue,
    "MinimumDirectGap" -> minimumGap,
    "Path" -> path,
    "OccupiedBands" -> occupied,
    "ClosureVector" -> closureVector,
    "ReciprocalCoordinates" -> Round[reciprocalCoordinates],
    "EndpointCovarianceResidual" -> covarianceResidual,
    "Convention" ->
      "occupied-subspace links follow the forward path; phase=Arg Det[WilsonMatrix] in radians"
  |>;
  Switch[output,
    "Phase", data["Phase"],
    "PhaseOverPi", data["PhaseOverPi"],
    "WilsonDeterminant", data["WilsonDeterminant"],
    "Data", data
  ]
];

berryCurvature::point =
  "The point must be a finite real momentum vector of positive dimension; received `1`.";
berryCurvature::directions =
  "Directions must contain two distinct valid coordinate indices; received `1` for a point of dimension `2`.";
berryCurvature::step =
  "StepSize must be a positive finite real number or a pair of such numbers; received `1`.";
berryCurvature::option =
  "Invalid berryCurvature option value(s): `1`.";

Options[berryCurvature] = {
  "Directions" -> {1, 2},
  "StepSize" -> 10^-3,
  "HermitianTolerance" -> 10^-10,
  "GapTolerance" -> 10^-9,
  "CovarianceTolerance" -> 10^-8,
  "OverlapTolerance" -> 10^-10,
  "Output" -> "Curvature"
};

berryCurvature[
    hamiltonian_, centers_, occupied_, point_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, directions, stepSetting,
    steps, hermitianTolerance, gapTolerance, covarianceTolerance,
    overlapTolerance, output, firstDisplacement, secondDisplacement,
    path, phaseData, area, flux, curvature, data
  },
  unknown = berryUnknownOptions[berryCurvature, supplied];
  If[unknown =!= {},
    Message[berryCurvature::option, unknown];
    Return[$Failed]
  ];
  directions = OptionValue["Directions"];
  stepSetting = OptionValue["StepSize"];
  hermitianTolerance = OptionValue["HermitianTolerance"];
  gapTolerance = OptionValue["GapTolerance"];
  covarianceTolerance = OptionValue["CovarianceTolerance"];
  overlapTolerance = OptionValue["OverlapTolerance"];
  output = OptionValue["Output"];
  If[
    !ListQ[point] || Length[point] < 2 ||
      !wilsonLoopRealVectorQ[point],
    Message[berryCurvature::point, point];
    Return[$Failed]
  ];
  If[
    !ListQ[directions] || Length[directions] =!= 2 ||
      !VectorQ[directions, IntegerQ] ||
      !DuplicateFreeQ[directions] ||
      !And @@ (Between[#, {1, Length[point]}] & /@ directions),
    Message[berryCurvature::directions, directions, Length[point]];
    Return[$Failed]
  ];
  steps = Which[
    wilsonLoopFiniteNumberQ[stepSetting] &&
      TrueQ[Im[N[stepSetting]] == 0] && TrueQ[stepSetting > 0],
      ConstantArray[stepSetting, 2],
    ListQ[stepSetting] && Length[stepSetting] === 2 &&
      And @@ (
        wilsonLoopFiniteNumberQ[#] && TrueQ[Im[N[#]] == 0] &&
          TrueQ[# > 0] & /@ stepSetting
      ),
      stepSetting,
    True,
      Message[berryCurvature::step, stepSetting];
      Return[$Failed]
  ];
  If[
    !And @@ (
      wilsonLoopFiniteNumberQ[#] && TrueQ[Im[N[#]] == 0] &&
        TrueQ[Re[N[#]] >= 0] & /@ {
        hermitianTolerance, gapTolerance, covarianceTolerance,
        overlapTolerance
      }
    ) || !MemberQ[{"Curvature", "Flux", "Data"}, output],
    Message[berryCurvature::option, {
      hermitianTolerance, gapTolerance, covarianceTolerance,
      overlapTolerance, output
    }];
    Return[$Failed]
  ];
  firstDisplacement = ConstantArray[0, Length[point]];
  secondDisplacement = ConstantArray[0, Length[point]];
  firstDisplacement[[directions[[1]]]] = steps[[1]]/2;
  secondDisplacement[[directions[[2]]]] = steps[[2]]/2;
  path = {
    point - firstDisplacement - secondDisplacement,
    point + firstDisplacement - secondDisplacement,
    point + firstDisplacement + secondDisplacement,
    point - firstDisplacement + secondDisplacement,
    point - firstDisplacement - secondDisplacement
  };
  phaseData = berryPhase[
    hamiltonian,
    centers,
    occupied,
    path,
    "HermitianTolerance" -> hermitianTolerance,
    "GapTolerance" -> gapTolerance,
    "CovarianceTolerance" -> covarianceTolerance,
    "OverlapTolerance" -> overlapTolerance,
    "Output" -> "Data"
  ];
  If[phaseData === $Failed, Return[$Failed]];
  area = Times @@ steps;
  flux = phaseData["Phase"];
  curvature = flux/area;
  data = <|
    "Curvature" -> curvature,
    "Flux" -> flux,
    "Area" -> area,
    "Point" -> point,
    "Directions" -> directions,
    "StepSize" -> steps,
    "PlaquettePath" -> path,
    "PhaseData" -> phaseData,
    "Convention" ->
      "positive orientation follows Directions[[1]] then Directions[[2]]"
  |>;
  Switch[output,
    "Curvature", data["Curvature"],
    "Flux", data["Flux"],
    "Data", data
  ]
];

End[]

EndPackage[]
