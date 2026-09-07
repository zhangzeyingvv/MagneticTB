(* ::Package:: *)

BeginPackage["RepresentationValidation`"]

HamiltonianSymmetryReport::usage =
  "HamiltonianSymmetryReport[H,representation,kSymbols,parameters] numerically certifies Gamma semilinear invariance and generic-k covariance.";
VerifyHamiltonianSymmetry::usage =
  "VerifyHamiltonianSymmetry[...] returns the Verified field of HamiltonianSymmetryReport.";

Begin["`Private`"]

validationMaximum[values_] :=
  If[Flatten[values] === {}, 0, Max[Flatten[values]]];

automaticParameterRuleSets[parameters_List] :=
  If[parameters === {}, {{}},
    (Thread[parameters -> #] &) /@ IdentityMatrix[Length[parameters]]];

Options[HamiltonianSymmetryReport] = {
  "SampleK" -> Automatic,
  "ParameterRuleSets" -> Automatic,
  "WorkingPrecision" -> 40,
  "Tolerance" -> 10^-24
};

HamiltonianSymmetryReport::hamiltonian = "hamiltonian must be a square matrix.";
HamiltonianSymmetryReport::representation =
  "representation must supply same-length square RepresentationMatrices, SpatialActions, and Boolean AntiunitaryFlags.";
HamiltonianSymmetryReport::k =
  "kSymbols and SampleK must have the coordinate dimension used by every rotation.";
HamiltonianSymmetryReport::rules =
  "ParameterRuleSets must be a list of replacement-rule lists.";

HamiltonianSymmetryReport[
    hamiltonian_, representation_Association, kSymbols_List,
    parameters_List, OptionsPattern[]
  ] := Module[
  {matrices, spatialActions, flags, operationIndices, unitaryIndices,
   antiunitaryIndices, sampleK, parameterRuleSets, precision, tolerance,
   gammaSemilinearResiduals, gammaUnitaryCommutatorResiduals,
   genericCovarianceResiduals, hGamma, hAtK, hAtImage, matrix, operation,
   imageK, gammaSemilinearByOperation, gammaCommutatorByOperation,
   genericByOperation, maximumResidual},
  If[!MatrixPredicates`ExactSquareMatrixQ[hamiltonian],
    Message[HamiltonianSymmetryReport::hamiltonian]; Return[$Failed]
  ];
  matrices = Lookup[representation, "RepresentationMatrices", $Failed];
  spatialActions = Lookup[representation, "SpatialActions", $Failed];
  flags = Lookup[representation, "AntiunitaryFlags", $Failed];
  If[!ListQ[matrices] || matrices === {} || !ListQ[spatialActions] ||
      !ListQ[flags] || !And @@ (BooleanQ /@ flags) ||
      Length[matrices] =!= Length[spatialActions] ||
      Length[matrices] =!= Length[flags] ||
      !And @@ (MatrixPredicates`ExactSquareMatrixQ /@ matrices) ||
      !And @@ (Dimensions[#] === Dimensions[hamiltonian] & /@ matrices),
    Message[HamiltonianSymmetryReport::representation]; Return[$Failed]
  ];
  sampleK = Replace[OptionValue["SampleK"],
    Automatic :> Take[{Pi/7, Pi/11, Pi/13}, Length[kSymbols]]];
  If[!VectorQ[sampleK] || Length[sampleK] =!= Length[kSymbols] ||
      !And @@ (SymmetryAlgebra`ValidSpatialOperationQ[#, Length[kSymbols]] & /@
        spatialActions),
    Message[HamiltonianSymmetryReport::k]; Return[$Failed]
  ];
  parameterRuleSets = Replace[OptionValue["ParameterRuleSets"],
    Automatic :> automaticParameterRuleSets[parameters]];
  If[!ListQ[parameterRuleSets] ||
      !And @@ (ListQ[#] && VectorQ[#, MatchQ[#, _Rule] &] & /@
        parameterRuleSets),
    Message[HamiltonianSymmetryReport::rules]; Return[$Failed]
  ];
  precision = OptionValue["WorkingPrecision"];
  tolerance = OptionValue["Tolerance"];
  operationIndices = Range[Length[spatialActions]];
  unitaryIndices = Pick[operationIndices, flags, False];
  antiunitaryIndices = Pick[operationIndices, flags, True];

  gammaSemilinearResiduals = Table[
    hGamma = N[hamiltonian /. rules /. Thread[kSymbols -> 0], precision];
    Table[
      matrix = matrices[[operationIndex]];
      validationMaximum[Abs@Flatten@N[
        matrix.If[flags[[operationIndex]], Conjugate[hGamma], hGamma].
          ConjugateTranspose[matrix] - hGamma, precision - 10]],
      {operationIndex, operationIndices}],
    {rules, parameterRuleSets}];
  gammaUnitaryCommutatorResiduals = Table[
    hGamma = N[hamiltonian /. rules /. Thread[kSymbols -> 0], precision];
    Table[validationMaximum[Abs@Flatten@N[
      matrices[[operationIndex]].hGamma - hGamma.matrices[[operationIndex]],
      precision - 10]], {operationIndex, unitaryIndices}],
    {rules, parameterRuleSets}];
  genericCovarianceResiduals = Table[
    hAtK = N[hamiltonian /. rules /. Thread[kSymbols -> sampleK], precision];
    Table[
      operation = spatialActions[[operationIndex]];
      matrix = matrices[[operationIndex]];
      imageK = If[flags[[operationIndex]], -1, 1] *
        Inverse[Transpose[operation["Rotation"]]].sampleK;
      hAtImage = N[hamiltonian /. rules /. Thread[kSymbols -> imageK], precision];
      validationMaximum[Abs@Flatten@N[
        matrix.If[flags[[operationIndex]], Conjugate[hAtK], hAtK].
          ConjugateTranspose[matrix] - hAtImage, precision - 10]],
      {operationIndex, operationIndices}],
    {rules, parameterRuleSets}];

  gammaSemilinearByOperation = If[gammaSemilinearResiduals === {}, {},
    Max /@ Transpose[gammaSemilinearResiduals]];
  gammaCommutatorByOperation = If[gammaUnitaryCommutatorResiduals === {} ||
      First[gammaUnitaryCommutatorResiduals] === {}, {},
    Max /@ Transpose[gammaUnitaryCommutatorResiduals]];
  genericByOperation = If[genericCovarianceResiduals === {}, {},
    Max /@ Transpose[genericCovarianceResiduals]];
  maximumResidual = validationMaximum[{
    gammaSemilinearByOperation, gammaCommutatorByOperation,
    genericByOperation}];
  <|
    "Convention" ->
      "U_g H(k)^(eta_g) U_g^dagger = H((-1)^eta_g p_g^(-T) k)",
    "SampleK" -> sampleK,
    "UnitaryOperationIndices" -> unitaryIndices,
    "AntiunitaryOperationIndices" -> antiunitaryIndices,
    "GammaSemilinearResidualsByOperation" -> gammaSemilinearByOperation,
    "GammaUnitaryCommutatorResidualsByOperation" ->
      gammaCommutatorByOperation,
    "GenericCovarianceResidualsByOperation" -> genericByOperation,
    "MaximumResidual" -> maximumResidual,
    "Tolerance" -> tolerance,
    "Verified" -> TrueQ[maximumResidual < tolerance]
  |>
];
HamiltonianSymmetryReport[___] :=
  (Message[HamiltonianSymmetryReport::hamiltonian]; $Failed);

Options[VerifyHamiltonianSymmetry] = Options[HamiltonianSymmetryReport];
VerifyHamiltonianSymmetry[
    hamiltonian_, representation_Association, kSymbols_List,
    parameters_List, options : OptionsPattern[]
  ] := Module[{report},
  report = HamiltonianSymmetryReport[
    hamiltonian, representation, kSymbols, parameters, options];
  AssociationQ[report] && TrueQ[report["Verified"]]
];

End[]

EndPackage[]
