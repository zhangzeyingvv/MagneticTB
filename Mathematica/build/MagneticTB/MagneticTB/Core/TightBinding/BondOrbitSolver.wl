(* ::Package:: *)

BeginPackage["MagneticTBLinearAlgebra`"]

SolveRectangularStabilizer::usage =
  "SolveRectangularStabilizer[{m,n},operations] compiles and solves one representative hopping fixed-space problem.";
SolveBondConstraintProblem::usage =
  "SolveBondConstraintProblem[problem] solves the pure constraint blocks in a BondConstraintProblem.";

Begin["`Private`"]

ClearAll[rectangularDimensionsQ, kernelSummary];

rectangularDimensionsQ[{_Integer?Positive, _Integer?Positive}] := True;
rectangularDimensionsQ[_] := False;

kernelSummary[kernel_Association] := KeyTake[kernel, {
  "ConstraintBlocks", "ConstraintMatrix", "NullSpaceRows", "BasisMatrix",
  "Rank", "Nullity", "Method", "ValidationLevel", "ResidualVerified",
  "IndependentVerified", "RankNullityVerified", "IterationNullities",
  "Verified"
}];

Options[SolveRectangularStabilizer] = Join[
  Options[LinearConstraintKernel`SolveConstraintKernel],
  {"ContinuousGeneratorPair" -> None}
];
SolveRectangularStabilizer::dims =
  "Expected positive matrix dimensions {m,n}; received `1`.";
SolveRectangularStabilizer::operation =
  "The finite operation records or continuous generator pair are incompatible with dimensions `1`.";

SolveRectangularStabilizer[
    dims_?rectangularDimensionsQ, operations_List, OptionsPattern[]
  ] := Module[{blocks, coordinateDimension, kernel, actions},
  coordinateDimension = 2 Times @@ dims;
  blocks = CompileRectangularConstraintBlocks[
    dims, operations, OptionValue["ContinuousGeneratorPair"]];
  If[blocks === $Failed,
    Message[SolveRectangularStabilizer::operation, dims]; Return[$Failed]
  ];
  kernel = LinearConstraintKernel`SolveConstraintKernel[
    blocks, coordinateDimension,
    "NullSpaceFunction" -> OptionValue["NullSpaceFunction"],
    "Method" -> OptionValue["Method"],
    "ValidationLevel" -> OptionValue["ValidationLevel"]];
  If[kernel === $Failed, Return[$Failed]];
  actions = RectangularTransportMatrix[
      Lookup[#, "Left", $Failed], Lookup[#, "Right", $Failed],
      "Antiunitary" -> Lookup[#, "Antiunitary", False]] & /@ operations;
  Join[
    <|
      "Dimensions" -> dims,
      "CoordinateOrder" ->
        "row-major real parts, followed by row-major imaginary parts",
      "Operations" -> operations,
      "ActionMatrices" -> actions,
      "KernelMethod" -> kernel["Method"],
      "KernelVerified" -> kernel["Verified"]
    |>,
    kernelSummary[kernel]
  ]
];
SolveRectangularStabilizer[dims_, _, OptionsPattern[]] :=
  (Message[SolveRectangularStabilizer::dims, dims]; $Failed);

Options[SolveBondConstraintProblem] =
  Options[LinearConstraintKernel`SolveConstraintKernel];
SolveBondConstraintProblem::compiled =
  "Expected a BondConstraintProblem whose orbits contain valid Dimensions, CoordinateDimension, and pure ConstraintBlocks.";

SolveBondConstraintProblem[problem_Association, OptionsPattern[]] := Module[
  {orbits, solutions},
  orbits = Lookup[problem, "Orbits", $Failed];
  If[Lookup[problem, "Schema", None] =!= "BondConstraintProblem" ||
      !ListQ[orbits] ||
      !And @@ Map[
        Function[orbit,
          Module[{dims, coordinateDimension, blocks},
            If[!AssociationQ[orbit], Return[False, Module]];
            dims = Lookup[orbit, "Dimensions", $Failed];
            coordinateDimension = Lookup[orbit, "CoordinateDimension", $Failed];
            blocks = Lookup[orbit, "ConstraintBlocks", $Failed];
            rectangularDimensionsQ[dims] &&
              coordinateDimension === 2 Times @@ dims &&
              ListQ[blocks] &&
              And @@ (MatrixQ[#] && Last[Dimensions[#]] === coordinateDimension & /@ blocks)
          ]
        ],
        orbits
      ],
    Message[SolveBondConstraintProblem::compiled]; Return[$Failed]
  ];
  solutions = Map[
    Function[orbit,
      Module[{kernel},
        kernel = LinearConstraintKernel`SolveConstraintKernel[
          orbit["ConstraintBlocks"], orbit["CoordinateDimension"],
          "NullSpaceFunction" -> OptionValue["NullSpaceFunction"],
          "Method" -> OptionValue["Method"],
          "ValidationLevel" -> OptionValue["ValidationLevel"]];
        If[kernel === $Failed, $Failed,
          Join[
            <|
              "Dimensions" -> orbit["Dimensions"],
              "CoordinateOrder" ->
                "row-major real parts, followed by row-major imaginary parts",
              "KernelMethod" -> kernel["Method"],
              "KernelVerified" -> kernel["Verified"]
            |>,
            kernelSummary[kernel]
          ]
        ]
      ]
    ],
    orbits
  ];
  If[MemberQ[solutions, $Failed], Return[$Failed]];
  <|
    "Schema" -> "SolvedBondConstraintProblem",
    "SchemaVersion" -> 1,
    "StaticShellKey" -> Lookup[problem, "StaticShellKey",
      Lookup[problem, "Shell", Missing["NotAvailable"]]],
    "Shell" -> Lookup[problem, "Shell", Missing["NotAvailable"]],
    "Hermitian" -> Lookup[problem, "Hermitian", Missing["NotAvailable"]],
    "Solutions" -> solutions
  |>
];
SolveBondConstraintProblem[_, OptionsPattern[]] :=
  (Message[SolveBondConstraintProblem::compiled]; $Failed);

End[]

EndPackage[]
