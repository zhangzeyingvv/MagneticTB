(* ::Package:: *)

BeginPackage["MagneticTBLinearAlgebra`"]

RectangularTransportMatrix::usage =
  "RectangularTransportMatrix[left,right] represents T |-> left.T.right^dagger on row-major real coordinates. Set \"Antiunitary\"->True for T |-> left.Conjugate[T].right^dagger.";
RectangularGeneratorConstraintMatrix::usage =
  "RectangularGeneratorConstraintMatrix[left,right] represents left.T-T.right==0 on row-major real coordinates.";
CompileRectangularConstraintBlocks::usage =
  "CompileRectangularConstraintBlocks[{m,n},operations,continuous] compiles finite Same/Reverse operations and an optional continuous-generator pair into pure real constraint blocks.";
StabilizerConstraintMatrix::usage =
  "StabilizerConstraintMatrix[actions] stacks action-identity for a common fixed-space problem.";
InvariantBasis::usage =
  "InvariantBasis[actions] returns columns spanning the common fixed space.";
SolveBasisStabilizer::usage =
  "SolveBasisStabilizer[basis,transforms] solves fixed-space constraints in an arbitrary real operator basis.";

Begin["`Private`"]

Options[RectangularTransportMatrix] = {"Antiunitary" -> False};
RectangularTransportMatrix::left = "The left representation matrix must be square.";
RectangularTransportMatrix::right = "The right representation matrix must be square.";
RectangularTransportMatrix::antiunitary = "The Antiunitary option must be True or False; received `1`.";

RectangularTransportMatrix[left_, right_, OptionsPattern[]] := Module[
  {complexAction, realPart, imaginaryPart, antiunitary},
  If[!MatrixPredicates`ExactSquareMatrixQ[left],
    Message[RectangularTransportMatrix::left];
    Return[$Failed]
  ];
  If[!MatrixPredicates`ExactSquareMatrixQ[right],
    Message[RectangularTransportMatrix::right];
    Return[$Failed]
  ];
  antiunitary = OptionValue["Antiunitary"];
  If[!BooleanQ[antiunitary],
    Message[RectangularTransportMatrix::antiunitary, antiunitary];
    Return[$Failed]
  ];
  complexAction = KroneckerProduct[left, Conjugate[right]];
  realPart = Simplify@ComplexExpand@Re[complexAction];
  imaginaryPart = Simplify@ComplexExpand@Im[complexAction];
  If[antiunitary,
    ArrayFlatten[{{realPart, imaginaryPart}, {imaginaryPart, -realPart}}],
    ArrayFlatten[{{realPart, -imaginaryPart}, {imaginaryPart, realPart}}]
  ]
];

RectangularGeneratorConstraintMatrix::left =
  "The left continuous generator must be square.";
RectangularGeneratorConstraintMatrix::right =
  "The right continuous generator must be square.";

RectangularGeneratorConstraintMatrix[left_, right_] := Module[
  {leftDimension, rightDimension, complexConstraint, realPart, imaginaryPart},
  If[!MatrixPredicates`ExactSquareMatrixQ[left],
    Message[RectangularGeneratorConstraintMatrix::left];
    Return[$Failed]
  ];
  If[!MatrixPredicates`ExactSquareMatrixQ[right],
    Message[RectangularGeneratorConstraintMatrix::right];
    Return[$Failed]
  ];
  leftDimension = Length[left];
  rightDimension = Length[right];
  complexConstraint =
    KroneckerProduct[left, IdentityMatrix[rightDimension]] -
      KroneckerProduct[IdentityMatrix[leftDimension], Transpose[right]];
  realPart = Simplify@ComplexExpand@Re[complexConstraint];
  imaginaryPart = Simplify@ComplexExpand@Im[complexConstraint];
  ArrayFlatten[{{realPart, -imaginaryPart}, {imaginaryPart, realPart}}]
];

CompileRectangularConstraintBlocks::data =
  "Expected dimensions {m,n}, valid Left/Right/Antiunitary/Target operation records, and an optional matching Hermitian generator pair.";
CompileRectangularConstraintBlocks[
    dims : {m_Integer?Positive, n_Integer?Positive},
    operations_List, continuousPair_ : None
  ] := Module[{coordinateDimension = 2 m n, blocks, continuousBlock},
  blocks = Map[
    Function[operation,
      Module[{left, right, antiunitary, target, action},
        If[!AssociationQ[operation], Return[$Failed, Module]];
        left = Lookup[operation, "Left", $Failed];
        right = Lookup[operation, "Right", $Failed];
        antiunitary = Lookup[operation, "Antiunitary", False];
        target = Lookup[operation, "Target", "Same"];
        If[!BooleanQ[antiunitary] ||
            !MatrixPredicates`ExactSquareMatrixQ[left] ||
            !MatrixPredicates`ExactSquareMatrixQ[right] ||
            Dimensions[left] =!= {m, m} || Dimensions[right] =!= {n, n} ||
            !MemberQ[{"Same", "Reverse"}, target] ||
            (target === "Reverse" && m =!= n),
          Return[$Failed, Module]
        ];
        action = RectangularTransportMatrix[
          left, right, "Antiunitary" -> antiunitary];
        If[action === $Failed, Return[$Failed, Module]];
        Switch[target,
          "Same", action - IdentityMatrix[coordinateDimension],
          "Reverse", action - OperatorSpaceBasis`DaggerCoordinateMatrix[dims]
        ]
      ]
    ],
    operations
  ];
  If[MemberQ[blocks, $Failed],
    Message[CompileRectangularConstraintBlocks::data]; Return[$Failed]
  ];
  If[continuousPair =!= None,
    If[!MatchQ[continuousPair, {_?MatrixPredicates`ExactHermitianMatrixQ,
          _?MatrixPredicates`ExactHermitianMatrixQ}] ||
        Dimensions[continuousPair[[1]]] =!= {m, m} ||
        Dimensions[continuousPair[[2]]] =!= {n, n},
      Message[CompileRectangularConstraintBlocks::data]; Return[$Failed]
    ];
    continuousBlock = RectangularGeneratorConstraintMatrix @@ continuousPair;
    If[continuousBlock === $Failed, Return[$Failed]];
    blocks = Prepend[blocks, continuousBlock]
  ];
  blocks
];
CompileRectangularConstraintBlocks[___] :=
  (Message[CompileRectangularConstraintBlocks::data]; $Failed);

StabilizerConstraintMatrix::actions =
  "All actions must be nonempty square matrices of the same dimension.";
StabilizerConstraintMatrix[actions_List] /; actions =!= {} := Module[
  {dimensions, dimension},
  dimensions = Dimensions /@ actions;
  If[!And @@ (MatrixPredicates`ExactSquareMatrixQ /@ actions) ||
      !SameQ @@ dimensions,
    Message[StabilizerConstraintMatrix::actions]; Return[$Failed]
  ];
  dimension = First[First[dimensions]];
  Join @@ (Simplify[# - IdentityMatrix[dimension]] & /@ actions)
];
StabilizerConstraintMatrix[_] :=
  (Message[StabilizerConstraintMatrix::actions]; $Failed);

Options[InvariantBasis] = Options[LinearConstraintKernel`SolveConstraintKernel];
InvariantBasis::actions = StabilizerConstraintMatrix::actions;
InvariantBasis[actions_List, OptionsPattern[]] /; actions =!= {} := Module[
  {dimensions, coordinateDimension, blocks, kernel},
  dimensions = Dimensions /@ actions;
  If[!And @@ (MatrixPredicates`ExactSquareMatrixQ /@ actions) ||
      !SameQ @@ dimensions,
    Message[InvariantBasis::actions]; Return[$Failed]
  ];
  coordinateDimension = First[First[dimensions]];
  blocks = Simplify[# - IdentityMatrix[coordinateDimension]] & /@ actions;
  kernel = LinearConstraintKernel`SolveConstraintKernel[
    blocks, coordinateDimension,
    "NullSpaceFunction" -> OptionValue["NullSpaceFunction"],
    "Method" -> OptionValue["Method"],
    "ValidationLevel" -> OptionValue["ValidationLevel"]];
  If[kernel === $Failed, $Failed, kernel["BasisMatrix"]]
];
InvariantBasis[_, OptionsPattern[]] :=
  (Message[InvariantBasis::actions]; $Failed);

Options[SolveBasisStabilizer] =
  Options[LinearConstraintKernel`SolveConstraintKernel];
SolveBasisStabilizer::basis =
  "The operator basis must be nonempty and every transform must preserve its span.";
SolveBasisStabilizer[
    basis_List, transforms_List, OptionsPattern[]
  ] /; basis =!= {} := Module[{actions, blocks, kernel},
  actions = OperatorSpaceBasis`ActionMatrixInBasis[basis, #] & /@ transforms;
  If[MemberQ[actions, $Failed],
    Message[SolveBasisStabilizer::basis]; Return[$Failed]
  ];
  blocks = (# - IdentityMatrix[Length[basis]]) & /@ actions;
  kernel = LinearConstraintKernel`SolveConstraintKernel[
    blocks, Length[basis],
    "NullSpaceFunction" -> OptionValue["NullSpaceFunction"],
    "Method" -> OptionValue["Method"],
    "ValidationLevel" -> OptionValue["ValidationLevel"]];
  If[kernel === $Failed, Return[$Failed]];
  Join[<|"OperatorBasis" -> basis, "ActionMatrices" -> actions|>, kernel]
];
SolveBasisStabilizer[_, _, OptionsPattern[]] :=
  (Message[SolveBasisStabilizer::basis]; $Failed);

End[]

EndPackage[]
