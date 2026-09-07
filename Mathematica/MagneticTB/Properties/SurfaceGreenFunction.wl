(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  surfaceSliceAssociation,
  surfacePrincipalLayerBlocks,
  surfaceGreenFunction,
  surfaceSpectralFunction
];

surfaceSliceAssociation[data_Association, momentum_List] := Module[
  {translations, matrices, slices, translation, normalIndex, contribution},
  translations = data["Translations"];
  matrices = hoppingEffectiveMatrices[data];
  slices = <||>;
  Do[
    translation = translations[[index]];
    normalIndex = translation[[3]];
    contribution = matrices[[index]] Exp[
      2 Pi I momentum.translation[[;; 2]]
    ];
    AssociateTo[
      slices,
      ToString[normalIndex, InputForm] ->
        (Lookup[
          slices,
          ToString[normalIndex, InputForm],
          SparseArray[{}, Dimensions[contribution]]
        ] + contribution)
    ],
    {index, Length[translations]}
  ];
  slices
];

surfacePrincipalLayerBlocks[data_Association, momentum_List] := Module[
  {
    orbitalDimension, normalRange, layerThickness, slices, zero,
    slice, h00, h01, h10
  },
  orbitalDimension = data["NumWannier"];
  normalRange = Max[Abs[data["Translations"][[All, 3]]]];
  layerThickness = Max[1, normalRange];
  slices = surfaceSliceAssociation[data, momentum];
  zero = SparseArray[{}, {orbitalDimension, orbitalDimension}];
  slice[normalIndex_Integer] := Lookup[
    slices,
    ToString[normalIndex, InputForm],
    zero
  ];
  h00 = Total@Flatten[Table[
    hoppingEmbedBlock[
      slice[targetLayer - sourceLayer],
      sourceLayer + 1,
      targetLayer + 1,
      orbitalDimension,
      layerThickness
    ],
    {sourceLayer, 0, layerThickness - 1},
    {targetLayer, 0, layerThickness - 1}
  ], 1];
  h01 = Total@Flatten[Table[
    hoppingEmbedBlock[
      slice[layerThickness + targetLayer - sourceLayer],
      sourceLayer + 1,
      targetLayer + 1,
      orbitalDimension,
      layerThickness
    ],
    {sourceLayer, 0, layerThickness - 1},
    {targetLayer, 0, layerThickness - 1}
  ], 1];
  h10 = Total@Flatten[Table[
    hoppingEmbedBlock[
      slice[-layerThickness + targetLayer - sourceLayer],
      sourceLayer + 1,
      targetLayer + 1,
      orbitalDimension,
      layerThickness
    ],
    {sourceLayer, 0, layerThickness - 1},
    {targetLayer, 0, layerThickness - 1}
  ], 1];
  <|
    "H00" -> h00,
    "H01" -> h01,
    "H10" -> h10,
    "PrincipalLayerThickness" -> layerThickness,
    "OrbitalsPerCell" -> orbitalDimension,
    "PrincipalLayerDimension" -> layerThickness orbitalDimension
  |>
];

surfaceGreenFunction::data = transformHoppings::data;
surfaceGreenFunction::momentum =
  "The surface momentum must be a real numeric two-vector in reciprocal fractional coordinates; received `1`.";
surfaceGreenFunction::energy =
  "The energy must be a finite real number; received `1`.";
surfaceGreenFunction::cell =
  "CellMatrix must be Automatic or a nonsingular 3 by 3 integer matrix whose first two rows span the surface plane and whose third row is the stacking direction; received `1`.";
surfaceGreenFunction::hermitian = transformHoppings::hermitian;
surfaceGreenFunction::blocks =
  "The principal-layer hopping blocks are inconsistent with Hermiticity: H00 residual `1`, H10-H01^dagger residual `2`, tolerance `3`.";
surfaceGreenFunction::solve =
  "The surface Green-function linear solve failed at decimation iteration `1`.";
surfaceGreenFunction::convergence =
  "Surface Green-function decimation did not converge after `1` iterations; final coupling residual was `2` with tolerance `3`.";
surfaceGreenFunction::option =
  "Invalid surfaceGreenFunction option value(s): `1`.";

Options[surfaceGreenFunction] = {
  "CellMatrix" -> Automatic,
  "Broadening" -> 10^-3,
  "Tolerance" -> 10^-10,
  "MaxIterations" -> 200,
  "Surface" -> "Positive",
  "HermiticityTolerance" -> 10^-9,
  "Output" -> "GreenFunction"
};

surfaceGreenFunction[
    input_, momentum_, energy_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, cellMatrix, broadening,
    tolerance, maxIterations, surface, hermiticityTolerance, output,
    data, hoppingResidual, blocks, h00, h01, h10, h00Residual,
    couplingResidual, dimension, identity, zMatrix, epsilon,
    epsilonSurface, alpha, beta, green, forward, backward,
    alphaNew, betaNew, iteration = 0, residual = Infinity,
    converged = False, solveFailure = None, surfaceGreen,
    dysonResidual, spectralWeight
  },
  unknown = hoppingUnknownOptions[surfaceGreenFunction, supplied];
  If[unknown =!= {},
    Message[surfaceGreenFunction::option, unknown];
    Return[$Failed]
  ];
  cellMatrix = OptionValue["CellMatrix"];
  broadening = OptionValue["Broadening"];
  tolerance = OptionValue["Tolerance"];
  maxIterations = OptionValue["MaxIterations"];
  surface = OptionValue["Surface"];
  hermiticityTolerance = OptionValue["HermiticityTolerance"];
  output = OptionValue["Output"];
  If[!hoppingRealVectorQ[momentum, 2],
    Message[surfaceGreenFunction::momentum, momentum];
    Return[$Failed]
  ];
  If[
    !hoppingFiniteNumberQ[energy] || !TrueQ[Im[N[energy]] == 0],
    Message[surfaceGreenFunction::energy, energy];
    Return[$Failed]
  ];
  If[cellMatrix =!= Automatic && !hoppingCellMatrixQ[cellMatrix],
    Message[surfaceGreenFunction::cell, cellMatrix];
    Return[$Failed]
  ];
  If[
    !hoppingFiniteNumberQ[broadening] ||
      !TrueQ[Im[N[broadening]] == 0] || !TrueQ[broadening > 0] ||
      !hoppingFiniteNumberQ[tolerance] ||
      !TrueQ[Im[N[tolerance]] == 0] || !TrueQ[tolerance > 0] ||
      !IntegerQ[maxIterations] || maxIterations < 1 ||
      !MemberQ[{"Positive", "Negative"}, surface] ||
      !hoppingFiniteNumberQ[hermiticityTolerance] ||
      !TrueQ[Im[N[hermiticityTolerance]] == 0] ||
      !TrueQ[hermiticityTolerance >= 0] ||
      !MemberQ[{"GreenFunction", "SpectralWeight", "Data"}, output],
    Message[
      surfaceGreenFunction::option,
      {
        broadening, tolerance, maxIterations, surface,
        hermiticityTolerance, output
      }
    ];
    Return[$Failed]
  ];
  data = hoppingResolveData[input];
  If[data === $Failed || !hoppingDataQ[data],
    Message[surfaceGreenFunction::data];
    Return[$Failed]
  ];
  If[cellMatrix =!= Automatic,
    data = transformHoppings[
      data,
      cellMatrix,
      "HermiticityTolerance" -> hermiticityTolerance
    ];
    If[data === $Failed, Return[$Failed]]
  ];
  hoppingResidual = hoppingHermiticityResidual[data];
  If[TrueQ[hoppingResidual > hermiticityTolerance],
    Message[
      surfaceGreenFunction::hermitian,
      hoppingResidual,
      hermiticityTolerance
    ];
    Return[$Failed]
  ];
  blocks = surfacePrincipalLayerBlocks[data, momentum];
  h00 = N@Normal@blocks["H00"];
  h01 = N@Normal@blocks["H01"];
  h10 = N@Normal@blocks["H10"];
  h00Residual = Norm[h00 - ConjugateTranspose[h00], "Frobenius"];
  couplingResidual = Norm[
    h10 - ConjugateTranspose[h01],
    "Frobenius"
  ];
  If[
    TrueQ[h00Residual > hermiticityTolerance] ||
      TrueQ[couplingResidual > hermiticityTolerance],
    Message[
      surfaceGreenFunction::blocks,
      h00Residual,
      couplingResidual,
      hermiticityTolerance
    ];
    Return[$Failed]
  ];
  dimension = blocks["PrincipalLayerDimension"];
  identity = IdentityMatrix[dimension];
  zMatrix = (N[energy] + I N[broadening]) identity;
  epsilon = h00;
  epsilonSurface = h00;
  alpha = h01;
  beta = h10;
  Do[
    green = Quiet@Check[
      LinearSolve[zMatrix - epsilon, identity],
      $Failed
    ];
    If[green === $Failed || !MatrixQ[green, hoppingFiniteNumberQ],
      solveFailure = currentIteration;
      Break[]
    ];
    forward = alpha . green . beta;
    backward = beta . green . alpha;
    alphaNew = alpha . green . alpha;
    betaNew = beta . green . beta;
    epsilonSurface = epsilonSurface + If[
      surface === "Positive",
      forward,
      backward
    ];
    epsilon = epsilon + forward + backward;
    residual = Max[
      Norm[alphaNew, "Frobenius"],
      Norm[betaNew, "Frobenius"]
    ]/Max[1., Norm[epsilon, "Frobenius"]];
    alpha = alphaNew;
    beta = betaNew;
    iteration = currentIteration;
    If[TrueQ[residual <= tolerance],
      converged = True;
      Break[]
    ],
    {currentIteration, maxIterations}
  ];
  If[solveFailure =!= None,
    Message[surfaceGreenFunction::solve, solveFailure];
    Return[$Failed]
  ];
  If[!converged,
    Message[
      surfaceGreenFunction::convergence,
      maxIterations,
      residual,
      tolerance
    ];
    Return[$Failed]
  ];
  surfaceGreen = Quiet@Check[
    LinearSolve[zMatrix - epsilonSurface, identity],
    $Failed
  ];
  If[
    surfaceGreen === $Failed ||
      !MatrixQ[surfaceGreen, hoppingFiniteNumberQ],
    Message[surfaceGreenFunction::solve, "final"];
    Return[$Failed]
  ];
  dysonResidual = Norm[
    (zMatrix - epsilonSurface) . surfaceGreen - identity,
    "Frobenius"
  ]/Max[1., Norm[identity, "Frobenius"]];
  spectralWeight = Re@Chop[-Im[Tr[surfaceGreen]]/Pi];
  Switch[
    output,
    "GreenFunction",
      surfaceGreen,
    "SpectralWeight",
      spectralWeight,
    "Data",
      <|
        "Schema" -> "MagneticTBSurfaceGreenFunction",
        "SchemaVersion" -> 1,
        "GreenFunction" -> surfaceGreen,
        "SpectralWeight" -> spectralWeight,
        "Energy" -> energy,
        "Broadening" -> broadening,
        "SurfaceMomentum" -> momentum,
        "Surface" -> surface,
        "PrincipalLayerThickness" ->
          blocks["PrincipalLayerThickness"],
        "OrbitalsPerCell" -> blocks["OrbitalsPerCell"],
        "PrincipalLayerDimension" -> dimension,
        "H00" -> h00,
        "H01" -> h01,
        "H10" -> h10,
        "Iterations" -> iteration,
        "CouplingResidual" -> residual,
        "DysonResidual" -> dysonResidual,
        "HoppingHermitianResidual" -> hoppingResidual,
        "BlockHermitianResidual" -> Max[h00Residual, couplingResidual],
        "Convention" ->
          "surface plane=cell axes 1,2; stacking direction=cell axis 3"
      |>
  ]
];

Options[surfaceSpectralFunction] = DeleteCases[
  Options[surfaceGreenFunction],
  HoldPattern["Output" -> _]
];

surfaceSpectralFunction[
    input_, momentum_, energy_, suppliedOptions : OptionsPattern[]
  ] := surfaceGreenFunction[
  input,
  momentum,
  energy,
  suppliedOptions,
  "Output" -> "SpectralWeight"
];

End[]

EndPackage[]
