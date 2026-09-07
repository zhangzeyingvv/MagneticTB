(* ::Package:: *)

BeginPackage["LinearConstraintKernel`"]

AssembleConstraintMatrix::usage =
  "AssembleConstraintMatrix[blocks,d] vertically stacks linear constraint blocks with d columns. An empty block list produces a 0 by d sparse matrix.";
SolveConstraintKernel::usage =
  "SolveConstraintKernel[blocks,d] returns the common null space of linear constraint blocks as a diagnostic Association. Method \"Iterative\" successively restricts the candidate subspace; Method \"Stacked\" solves the vertically stacked matrix once. The requested method fails directly and is never replaced by the other method. ValidationLevel \"None\", \"Basic\", or \"Full\" selects no mathematical certification, residual certification, or residual plus independence and rank-nullity certification. SolveConstraintKernel[matrix] solves one already assembled constraint matrix. The kernel is independent of MagneticTB geometry and symmetry objects.";

Begin["`Private`"]

ClearAll[
  constraintBlockQ,
  exactZeroArrayQ,
  runNullSpaceFunction,
  validNullSpaceRowsQ
];

constraintBlockQ[block_, coordinateDimension_Integer?Positive] :=
  MatrixQ[block] && Length[Dimensions[block]] == 2 &&
    Last[Dimensions[block]] == coordinateDimension;

exactZeroArrayQ[array_] := And @@ Map[
  TrueQ[# === 0] || TrueQ[PossibleZeroQ[#]] &,
  Flatten[Normal[array]]
];

validNullSpaceRowsQ[rows_, coordinateDimension_Integer?Positive] :=
  rows === {} ||
    (MatrixQ[rows] && Length[Dimensions[rows]] == 2 &&
      Last[Dimensions[rows]] == coordinateDimension);

(* A user-supplied backend must return normally.  Messages, an unmatched
   Throw, or a malformed return value are all solver failures; callers never
   reinterpret them as an empty constraint or switch algorithms. *)
runNullSpaceFunction[nullSpaceFunction_, matrix_] := Module[{outcome},
  (* The inner Catch handles an untagged Throw.  The outer Catch handles any
     explicitly tagged Throw.  A normal return is wrapped so that it cannot
     be confused with a value carried by Throw. *)
  outcome = Catch[
    Catch[
      {"Returned", Check[nullSpaceFunction[matrix], $Failed]}
    ],
    _,
    Function[{value, tag}, {"Thrown"}]
  ];
  If[MatchQ[outcome, {"Returned", _}], Last[outcome], $Failed]
];

AssembleConstraintMatrix::dimension =
  "The coordinate dimension must be a positive integer; received `1`.";
AssembleConstraintMatrix::block =
  "Every constraint block must be a matrix with exactly `1` columns.";

AssembleConstraintMatrix[
    blocks_List,
    coordinateDimension_Integer?Positive
  ] := Module[{},
  If[!And @@ (constraintBlockQ[#, coordinateDimension] & /@ blocks),
    Message[AssembleConstraintMatrix::block, coordinateDimension];
    Return[$Failed]
  ];
  If[
    blocks === {},
    SparseArray[{}, {0, coordinateDimension}],
    Join @@ blocks
  ]
];
AssembleConstraintMatrix[_, coordinateDimension_] :=
  (Message[AssembleConstraintMatrix::dimension, coordinateDimension]; $Failed);

Options[SolveConstraintKernel] = {
  "NullSpaceFunction" -> Automatic,
  "Method" -> "Iterative",
  "ValidationLevel" -> "Basic"
};

SolveConstraintKernel::dimension =
  "The coordinate dimension must be a positive integer; received `1`.";
SolveConstraintKernel::matrix =
  "The supplied constraint matrix must have a positive number of columns.";
SolveConstraintKernel::solver =
  "The configured null-space function failed for method `1` at step `2`, or did not return row vectors of length `3`; no alternative method was attempted.";
SolveConstraintKernel::verify =
  "The returned vectors failed the certification requested by ValidationLevel `1`.";
SolveConstraintKernel::method =
  "Method must be \"Iterative\", \"Stacked\", or \"Cyclotomic\"; received `1`.";
SolveConstraintKernel::validation =
  "ValidationLevel must be \"None\", \"Basic\", or \"Full\"; received `1`.";

SolveConstraintKernel[
    blocks_List,
    coordinateDimension_Integer?Positive,
    OptionsPattern[]
  ] := Module[
  {constraintMatrix, nullSpaceFunction, nullSpaceRows, basisMatrix,
   rank, nullity, residual, independentQ, completeQ, verifiedQ,
   method, iterationNullities, restrictedMatrix, restrictedRows,
   currentDimension, validationLevel, residualQ, blockIndex,
   iterationResult, iterationFailureTag},
  constraintMatrix = AssembleConstraintMatrix[blocks, coordinateDimension];
  If[constraintMatrix === $Failed, Return[$Failed]];

  method = OptionValue["Method"];
  If[!MemberQ[{"Iterative", "Stacked", "Cyclotomic"}, method],
    Message[SolveConstraintKernel::method, method];
    Return[$Failed]
  ];
  validationLevel = OptionValue["ValidationLevel"];
  If[!MemberQ[{"None", "Basic", "Full"}, validationLevel],
    Message[SolveConstraintKernel::validation, validationLevel];
    Return[$Failed]
  ];

  nullSpaceFunction = Replace[
    OptionValue["NullSpaceFunction"],
    Automatic -> NullSpace
  ];

  If[blocks === {},
    nullSpaceRows = IdentityMatrix[coordinateDimension];
    iterationNullities = {coordinateDimension},
    Switch[
      method,
      "Cyclotomic",
        nullSpaceRows = MagneticTB`CyclotomicCommonNullSpace[
          Normal /@ blocks,
          "CoordinateDimension" -> coordinateDimension
        ];
        If[
          FailureQ[nullSpaceRows] ||
            !validNullSpaceRowsQ[nullSpaceRows, coordinateDimension],
          Message[
            SolveConstraintKernel::solver,
            "Cyclotomic",
            "common kernel",
            coordinateDimension
          ];
          Return[$Failed]
        ];
        iterationNullities = {coordinateDimension, Length[nullSpaceRows]},
      "Stacked",
        nullSpaceRows = runNullSpaceFunction[
          nullSpaceFunction,
          constraintMatrix
        ];
        If[!validNullSpaceRowsQ[nullSpaceRows, coordinateDimension],
          Message[
            SolveConstraintKernel::solver,
            "Stacked",
            "assembled matrix",
            coordinateDimension
          ];
          Return[$Failed]
        ];
        iterationNullities = {coordinateDimension, Length[nullSpaceRows]},
      "Iterative",
        basisMatrix = IdentityMatrix[coordinateDimension, SparseArray];
        iterationNullities = {coordinateDimension};
        blockIndex = 0;
        iterationFailureTag = Unique["SolveConstraintKernelFailure"];
        iterationResult = Catch[
          Do[
            blockIndex++;
            currentDimension = Last[Dimensions[basisMatrix]];
            If[currentDimension == 0,
              AppendTo[iterationNullities, 0];
              Continue[]
            ];
            restrictedMatrix = Simplify[block.basisMatrix];
            restrictedRows = runNullSpaceFunction[
              nullSpaceFunction,
              restrictedMatrix
            ];
            If[!validNullSpaceRowsQ[restrictedRows, currentDimension],
              Message[
                SolveConstraintKernel::solver,
                "Iterative",
                blockIndex,
                currentDimension
              ];
              Throw[$Failed, iterationFailureTag]
            ];
            basisMatrix = If[
              restrictedRows === {},
              SparseArray[{}, {coordinateDimension, 0}],
              Simplify[basisMatrix.Transpose[restrictedRows]]
            ];
            AppendTo[iterationNullities, Last[Dimensions[basisMatrix]]],
            {block, blocks}
          ];
          basisMatrix,
          iterationFailureTag
        ];
        If[iterationResult === $Failed, Return[$Failed]];
        basisMatrix = iterationResult;
        nullSpaceRows = If[
          Last[Dimensions[basisMatrix]] == 0,
          {},
          Transpose[basisMatrix]
        ]
    ]
  ];

  nullity = Length[nullSpaceRows];
  basisMatrix = If[
    nullSpaceRows === {},
    SparseArray[{}, {coordinateDimension, 0}],
    Transpose[nullSpaceRows]
  ];
  rank = coordinateDimension - nullity;
  residual = Missing["NotComputed"];
  residualQ = Missing["NotChecked"];
  independentQ = Missing["NotChecked"];
  completeQ = Missing["NotChecked"];
  verifiedQ = Missing["NotChecked"];

  If[validationLevel =!= "None",
    residual = If[
      blocks === {},
      SparseArray[{}, {0, nullity}],
      Simplify[constraintMatrix.basisMatrix]
    ];
    residualQ = exactZeroArrayQ[residual];
    verifiedQ = residualQ
  ];

  If[validationLevel === "Full",
    rank = If[blocks === {}, 0, MatrixRank[constraintMatrix]];
    independentQ = nullSpaceRows === {} ||
      MatrixRank[nullSpaceRows] == nullity;
    completeQ = rank + nullity == coordinateDimension;
    verifiedQ = residualQ && independentQ && completeQ
  ];

  If[validationLevel =!= "None" && !TrueQ[verifiedQ],
    Message[SolveConstraintKernel::verify, validationLevel];
    Return[$Failed]
  ];

  <|
    "CoordinateDimension" -> coordinateDimension,
    "Method" -> method,
    "ValidationLevel" -> validationLevel,
    "IterationNullities" -> iterationNullities,
    "ConstraintBlocks" -> blocks,
    "ConstraintMatrix" -> constraintMatrix,
    "NullSpaceRows" -> nullSpaceRows,
    "BasisMatrix" -> basisMatrix,
    "Rank" -> rank,
    "RankSource" -> If[
      validationLevel === "Full",
      "IndependentMatrixRank",
      "NullityInference"
    ],
    "Nullity" -> nullity,
    "Residual" -> residual,
    "ResidualVerified" -> residualQ,
    "IndependentVerified" -> independentQ,
    "RankNullityVerified" -> completeQ,
    "Verified" -> verifiedQ
  |>
];
SolveConstraintKernel[_, coordinateDimension_Integer, OptionsPattern[]] /;
    coordinateDimension <= 0 :=
  (Message[SolveConstraintKernel::dimension, coordinateDimension]; $Failed);

SolveConstraintKernel[matrix_?MatrixQ, OptionsPattern[]] := Module[
  {dimensions = Dimensions[matrix], coordinateDimension},
  If[Length[dimensions] =!= 2 || Last[dimensions] < 1,
    Message[SolveConstraintKernel::matrix];
    Return[$Failed]
  ];
  coordinateDimension = Last[dimensions];
  SolveConstraintKernel[
    {matrix},
    coordinateDimension,
    "NullSpaceFunction" -> OptionValue["NullSpaceFunction"],
    "Method" -> OptionValue["Method"],
    "ValidationLevel" -> OptionValue["ValidationLevel"]
  ]
];
SolveConstraintKernel[_, OptionsPattern[]] :=
  (Message[SolveConstraintKernel::matrix]; $Failed);

End[]

EndPackage[]
