(* ::Package:: *)

BeginPackage["MagneticTBLinearAlgebra`"]

ReconstructRectangularHopping::usage =
  "ReconstructRectangularHopping[basisMatrix,parameters,{m,n}] reconstructs an m by n complex hopping matrix.";
ReconstructFromBasis::usage =
  "ReconstructFromBasis[basisMatrix,parameters,operatorBasis] reconstructs an operator from a null-space basis.";
PropagateHopping::usage =
  "PropagateHopping[T,left,right] gives left.T.right^\[Dagger]; use \"Antiunitary\"->True to conjugate T first.";
ReconstructSolvedBondModel::usage =
  "ReconstructSolvedBondModel[model,p] reconstructs every directed hopping term in a solved bond model using real parameters p[orbit,index].";
AssembleBlochHamiltonian::usage =
  "AssembleBlochHamiltonian[blockDimensions,terms,k] assembles Sum Exp[I k.d] T_d directly from block terms, without parsing trigonometric expressions or globally simplifying the assembled matrix.";
BuildSolvedBondHamiltonian::usage =
  "BuildSolvedBondHamiltonian[model,k,p] reconstructs a compatibility solved model. BuildSolvedBondHamiltonian[static,problem,solved,k,p] assembles directly from separated static, constraint, and solution layers.";

Begin["`Private`"]

ClearAll[builderRectangularDimensionsQ];

builderRectangularDimensionsQ[{m_Integer?Positive, n_Integer?Positive}] := True;
builderRectangularDimensionsQ[_] := False;

ReconstructRectangularHopping::shape = "The basis matrix must have 2 m n rows and one column per parameter.";
ReconstructRectangularHopping[
    basisMatrix_?MatrixQ,
    parameters_List,
    dims_?builderRectangularDimensionsQ
  ] := Module[{m, n, coordinateDimension, coordinates},
  {m, n} = dims;
  coordinateDimension = 2 m n;
  If[
    Dimensions[basisMatrix] =!= {coordinateDimension, Length[parameters]},
    Message[ReconstructRectangularHopping::shape];
    Return[$Failed]
  ];
  coordinates = Simplify[basisMatrix.parameters];
  OperatorSpaceBasis`RealCoordinatesToMatrix[coordinates, dims]
];
ReconstructRectangularHopping[___] :=
  (Message[ReconstructRectangularHopping::shape]; $Failed);

ReconstructFromBasis::shape = "The basis matrix must have one row per operator-basis element and one column per parameter.";
ReconstructFromBasis[basisMatrix_?MatrixQ, parameters_List, operatorBasis_List] := Module[
  {coordinates},
  If[
    Dimensions[basisMatrix] =!= {Length[operatorBasis], Length[parameters]},
    Message[ReconstructFromBasis::shape];
    Return[$Failed]
  ];
  coordinates = Simplify[basisMatrix.parameters];
  Simplify@Total[MapThread[Times, {coordinates, operatorBasis}]]
];
ReconstructFromBasis[___] :=
  (Message[ReconstructFromBasis::shape]; $Failed);

Options[PropagateHopping] = {"Antiunitary" -> False};
PropagateHopping::data =
  "T, left, and right must have compatible matrix dimensions and Antiunitary must be True or False.";
PropagateHopping[T_?MatrixQ, left_?MatrixQ, right_?MatrixQ, OptionsPattern[]] := Module[
  {antiunitary = OptionValue["Antiunitary"]},
  If[!BooleanQ[antiunitary] ||
      !MatrixPredicates`ExactSquareMatrixQ[left] ||
      !MatrixPredicates`ExactSquareMatrixQ[right] ||
      Dimensions[T] =!= {Length[left], Length[right]},
    Message[PropagateHopping::data]; Return[$Failed]
  ];
  Simplify[left . If[antiunitary, ComplexExpand[Conjugate[T]], T] .
    ConjugateTranspose[right]]
];
PropagateHopping[___] := (Message[PropagateHopping::data]; $Failed);

ReconstructSolvedBondModel::solution =
  "Expected solved bond data containing Bonds and Orbits with Solution records.";
ReconstructSolvedBondModel[
    solvedModel_Association,
    parameterHead_Symbol
  ] := Module[
  {orbits, bonds, parametersByOrbit, representatives, terms},
  orbits = Lookup[solvedModel, "Orbits", $Failed];
  bonds = Lookup[solvedModel, "Bonds", $Failed];
  If[
    orbits === $Failed || bonds === $Failed ||
      !And @@ (
        AssociationQ[#] && KeyExistsQ[#, "Solution"] &&
          KeyExistsQ[#, "Dimensions"] &&
          KeyExistsQ[#, "MemberRecords"] & /@ orbits
      ),
    Message[ReconstructSolvedBondModel::solution];
    Return[$Failed]
  ];
  parametersByOrbit = MapIndexed[
    Function[{orbit, orbitIndex},
      Table[
        parameterHead[First[orbitIndex], parameterIndex],
        {parameterIndex, orbit["Solution"]["Nullity"]}
      ]
    ],
    orbits
  ];
  representatives = MapThread[
    ReconstructRectangularHopping[
      #1["Solution"]["BasisMatrix"],
      #2,
      #1["Dimensions"]
    ] &,
    {orbits, parametersByOrbit}
  ];
  If[MemberQ[representatives, $Failed], Return[$Failed]];
  terms = Flatten[
    MapIndexed[
      Function[{orbit, orbitPosition},
        Map[
          Function[member,
            Module[{operation, hopping, bond},
              operation = Lookup[member, "Operation", $Failed];
              If[!AssociationQ[operation] ||
                  !MemberQ[{"Direct", "Dagger"}, Lookup[member, "Mode", $Failed]] ||
                  !BooleanQ[Lookup[operation, "Antiunitary", $Failed]] ||
                  !IntegerQ[Lookup[member, "BondIndex", $Failed]] ||
                  !Between[member["BondIndex"], {1, Length[bonds]}],
                $Failed,
                hopping = PropagateHopping[
                  representatives[[First[orbitPosition]]],
                  operation["Left"],
                  operation["Right"],
                  "Antiunitary" -> operation["Antiunitary"]
                ];
                If[hopping === $Failed,
                  $Failed,
                  If[member["Mode"] === "Dagger",
                    hopping = Simplify@ComplexExpand@ConjugateTranspose[hopping]
                  ];
                  bond = bonds[[member["BondIndex"]]];
                  <|
                    "BondIndex" -> member["BondIndex"],
                    "OrbitIndex" -> First[orbitPosition],
                    "RowBlock" -> bond["RowSite"],
                    "ColumnBlock" -> bond["ColumnSite"],
                    "Displacement" -> bond["Displacement"],
                    "Matrix" -> hopping
                  |>
                ]
              ]
            ]
          ],
          orbit["MemberRecords"]
        ]
      ],
      orbits
    ],
    1
  ];
  If[MemberQ[terms, $Failed, Infinity], Return[$Failed]];
  <|
    "ParametersByOrbit" -> parametersByOrbit,
    "Parameters" -> Flatten[parametersByOrbit],
    "RepresentativeHoppings" -> representatives,
    "Terms" -> SortBy[terms, Lookup[#, "BondIndex"] &]
  |>
];
ReconstructSolvedBondModel[___] :=
  (Message[ReconstructSolvedBondModel::solution]; $Failed);

AssembleBlochHamiltonian::blocks = "Block dimensions must be a nonempty list of positive integers.";
AssembleBlochHamiltonian::term = "Each term needs valid \"RowBlock\", \"ColumnBlock\", \"Displacement\", and \"Matrix\" entries.";
validBlochTermQ[term_, blockDimensions_List, kDimension_Integer?NonNegative] :=
  Module[{row, column, displacement, hopping, phaseSign},
    If[!AssociationQ[term], Return[False]];
    row = Lookup[term, "RowBlock", $Failed];
    column = Lookup[term, "ColumnBlock", $Failed];
    displacement = Lookup[term, "Displacement", $Failed];
    hopping = Lookup[term, "Matrix", $Failed];
    phaseSign = Lookup[term, "PhaseSign", 1];
    IntegerQ[row] && IntegerQ[column] &&
      Between[row, {1, Length[blockDimensions]}] &&
      Between[column, {1, Length[blockDimensions]}] &&
      VectorQ[displacement] && Length[displacement] === kDimension &&
      MemberQ[{-1, 1}, phaseSign] && MatrixQ[hopping] &&
      Dimensions[hopping] === {
        blockDimensions[[row]], blockDimensions[[column]]}
  ];
validBlochTermQ[___] := False;

AssembleBlochHamiltonian[
    blockDimensions_List,
    terms_List,
    kVector_List
  ] := Module[{blocks, row, column, displacement, hopping, phaseSign},
  If[blockDimensions == {} || !And @@ (IntegerQ[#] && Positive[#] & /@ blockDimensions),
    Message[AssembleBlochHamiltonian::blocks];
    Return[$Failed]
  ];
  If[!VectorQ[kVector] ||
      !And @@ (validBlochTermQ[#, blockDimensions, Length[kVector]] & /@ terms),
    Message[AssembleBlochHamiltonian::term];
    Return[$Failed]
  ];
  blocks = Table[
    ConstantArray[0, {blockDimensions[[i]], blockDimensions[[j]]}],
    {i, Length[blockDimensions]}, {j, Length[blockDimensions]}
  ];
  Do[
    row = Lookup[term, "RowBlock", 0];
    column = Lookup[term, "ColumnBlock", 0];
    displacement = Lookup[term, "Displacement", $Failed];
    hopping = Lookup[term, "Matrix", $Failed];
    phaseSign = Lookup[term, "PhaseSign", 1];
    blocks[[row, column]] = blocks[[row, column]] +
      Exp[I phaseSign kVector.displacement] hopping,
    {term, terms}
  ];
  ArrayFlatten[blocks]
];
AssembleBlochHamiltonian[___] :=
  (Message[AssembleBlochHamiltonian::blocks]; $Failed);

BuildSolvedBondHamiltonian[
    solvedModel_Association,
    kVector_List,
    parameterHead_Symbol
  ] := Module[{reconstruction, hamiltonian, blockDimensions},
  reconstruction = ReconstructSolvedBondModel[solvedModel, parameterHead];
  If[reconstruction === $Failed, Return[$Failed]];
  blockDimensions = Lookup[solvedModel, "BlockDimensions", $Failed];
  If[blockDimensions === $Failed, Return[$Failed]];
  hamiltonian = AssembleBlochHamiltonian[
    blockDimensions,
    reconstruction["Terms"],
    kVector
  ];
  If[hamiltonian === $Failed, Return[$Failed]];
  Join[reconstruction, <|"Hamiltonian" -> hamiltonian|>]
];

BuildSolvedBondHamiltonian::layers =
  "Expected matching DirectedBondOrbitData, BondConstraintProblem, and SolvedBondConstraintProblem layers.";
BuildSolvedBondHamiltonian[
    staticData_Association,
    constraintProblem_Association,
    solvedProblem_Association,
    kVector_List,
    parameterHead_Symbol
  ] := Module[{orbits, solutions, transientSolvedModel},
  orbits = Lookup[constraintProblem, "Orbits", $Failed];
  solutions = Lookup[solvedProblem, "Solutions", $Failed];
  If[
    Lookup[staticData, "Schema", None] =!= "DirectedBondOrbitData" ||
      Lookup[constraintProblem, "Schema", None] =!=
        "BondConstraintProblem" ||
      Lookup[solvedProblem, "Schema", None] =!=
        "SolvedBondConstraintProblem" ||
      !ListQ[orbits] || !ListQ[solutions] ||
      Length[orbits] =!= Length[solutions] ||
      constraintProblem["StaticShellKey"] =!= staticData["Shell"] ||
      solvedProblem["StaticShellKey"] =!= staticData["Shell"],
    Message[BuildSolvedBondHamiltonian::layers];
    Return[$Failed]
  ];
  transientSolvedModel = <|
    "Bonds" -> staticData["Bonds"],
    "BlockDimensions" -> staticData["BlockDimensions"],
    "Orbits" -> MapThread[
      Join[#1, <|"Solution" -> #2|>] &,
      {orbits, solutions}
    ]
  |>;
  BuildSolvedBondHamiltonian[
    transientSolvedModel,
    kVector,
    parameterHead
  ]
];
BuildSolvedBondHamiltonian[___] :=
  (Message[BuildSolvedBondHamiltonian::layers]; $Failed);

End[]

EndPackage[]
