(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  spinSpaceGroupPointOperationRecord,
  compileSpinSpaceGroupBasisAction
];

spinSpaceGroupPointOperationRecord[
    element_Association,
    evaluatedLattice_?MatrixQ,
    variables_List,
    spinorQ_
  ] := Module[{record, spinRotation},
  record = <|
    "CoordinateImages" -> nativeCoordinateImages[
      element["space"][[1]],
      evaluatedLattice,
      variables
    ],
    "Antiunitary" -> (element["spin"][[2]] === 1)
  |>;
  If[!TrueQ[spinorQ], Return[record]];
  spinRotation =
    Transpose[evaluatedLattice] . element["spin"][[1]] .
      Inverse[Transpose[evaluatedLattice]];
  Join[
    record,
    <|
      "ValueMatrix" -> ExpToTrig@nativeSpinMatrix[spinRotation],
      "AntiunitaryMatrix" -> I PauliMatrix[2]
    |>
  ]
];

compileSpinSpaceGroupBasisAction::input =
  "Expected a nonempty scalar or two-component spinor basis, a nonempty spin-space-group element list, three coordinate variables, and a nonsingular 3 by 3 evaluated lattice.";
compileSpinSpaceGroupBasisAction::compile =
  "The spin-space-group action on the supplied ordered basis could not be compiled into exact unitary local matrices.";

compileSpinSpaceGroupBasisAction[
    rawBasis_List,
    groupElements_List,
    evaluatedLattice_?MatrixQ,
    variables_List
  ] := Module[
  {basis, scalarQ, spinorQ, operationRecords, matrices,
   antiunitaryFlags, actionRecords},
  basis = resolveMagneticTBBasisEntry /@ rawBasis;
  scalarQ = basis =!= {} && And @@ (!ListQ[#] & /@ basis);
  spinorQ = basis =!= {} &&
    And @@ (ListQ[#] && Length[#] === 2 & /@ basis);
  If[
    (!scalarQ && !spinorQ) || groupElements === {} ||
      !And @@ (
        SymmetryAlgebra`ValidSpinSpaceGroupElementQ /@ groupElements
      ) ||
      Dimensions[evaluatedLattice] =!= {3, 3} ||
      TrueQ[PossibleZeroQ[Det[evaluatedLattice]]] ||
      Length[variables] =!= 3 || !DuplicateFreeQ[variables] ||
      !And @@ (MatchQ[#, _Symbol] & /@ variables),
    Message[compileSpinSpaceGroupBasisAction::input];
    Return[$Failed]
  ];
  operationRecords = spinSpaceGroupPointOperationRecord[
      #,
      evaluatedLattice,
      variables,
      spinorQ
    ] & /@ groupElements;
  matrices =
    FunctionBasisRepresentation`PointOperationRepresentationMatrices[
      basis,
      operationRecords,
      variables
    ];
  If[matrices === $Failed,
    Message[compileSpinSpaceGroupBasisAction::compile];
    Return[$Failed]
  ];
  antiunitaryFlags =
    (TrueQ[Lookup[#, "Antiunitary", False]] & /@ operationRecords);
  actionRecords = MapThread[
    <|
      "GroupElement" -> #1,
      "LocalMatrix" -> #2,
      "Antiunitary" -> #3
    |> &,
    {groupElements, matrices, antiunitaryFlags}
  ];
  <|
    "Basis" -> basis,
    "PointOperationRecords" -> operationRecords,
    "LocalMatrices" -> matrices,
    "AntiunitaryFlags" -> antiunitaryFlags,
    "ActionRecords" -> actionRecords
  |>
];
compileSpinSpaceGroupBasisAction[___] := (
  Message[compileSpinSpaceGroupBasisAction::input];
  $Failed
);

End[]

EndPackage[]
