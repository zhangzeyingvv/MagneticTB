(* ::Package:: *)

(* Native implementations of the documented pre-refactor helper APIs.
   They use the in-package linear-algebra pipeline and never load Old or
   LegacyDataAdapter. *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  nativeCompatibilityBlockDiagonal,
  nativeCompatibilitySpinorQ,
  nativeUnconstrainedHopping
];

nativeCompatibilityBlockDiagonal[matrices_List] := Module[{dimensions},
  If[matrices === {}, Return[{}]];
  dimensions = Dimensions /@ matrices;
  If[
    !And @@ (Length[#] == 2 && SameQ @@ # & /@ dimensions),
    Return[$Failed]
  ];
  ArrayFlatten@Table[
    If[
      row == column,
      matrices[[row]],
      ConstantArray[0, {dimensions[[row, 1]], dimensions[[column, 2]]}]
    ],
    {row, Length[matrices]},
    {column, Length[matrices]}
  ]
];

nativeCompatibilitySpinorQ[basis_List] := Which[
  basis === {}, $Failed,
  And @@ (!ListQ[#] & /@ basis), False,
  And @@ (ListQ[#] && Length[#] == 2 & /@ basis), True,
  True, $Failed
];

(* Public compatibility names stay at the Interface boundary; the finite-group
   algorithms and their failure semantics remain owned by Core. *)
Options[MagneticTB`GenerateGroup] = Options[GroupEnumeration`GenerateGroup];
MagneticTB`GenerateGroup[
    generators_List,
    identityElement_,
    multiply_,
    OptionsPattern[]
  ] := GroupEnumeration`GenerateGroup[
  generators,
  identityElement,
  multiply,
  SameTest -> OptionValue[SameTest]
];
MagneticTB`GenerateGroup[___] := $Failed;

Options[MagneticTB`getGenerator] =
  Options[GroupEnumeration`FindConcreteGeneratorElements];
MagneticTB`getGenerator[
    elements_List,
    identityElement_,
    multiply_,
    OptionsPattern[]
  ] := GroupEnumeration`FindConcreteGeneratorElements[
  elements,
  identityElement,
  multiply,
  SameTest -> OptionValue[SameTest]
];
MagneticTB`getGenerator[___] := $Failed;

pointMatrix::input =
  "pointMatrix expects a nonempty symmetry-operation list, a nonempty scalar or two-component spinor basis, and a nonsingular 3 by 3 evaluated lattice.";
pointMatrix::compile =
  "The supplied basis could not be compiled into exact unitary point-operation matrices.";

pointMatrix[
    symmetryInformation_List,
    orbitBasis_List,
    evaluatedLattice_?MatrixQ
  ] := Module[
  {resolvedBasis, spinorQ, seitzOperations, spatialActions, spinActions,
   antiunitaryFlags, operationRecords, matrices},
  If[
    symmetryInformation === {} || orbitBasis === {} ||
      Dimensions[evaluatedLattice] =!= {3, 3} ||
      TrueQ[PossibleZeroQ[Det[evaluatedLattice]]],
    Message[pointMatrix::input];
    Return[$Failed]
  ];
  resolvedBasis = resolveMagneticTBBasisEntry /@ orbitBasis;
  spinorQ = nativeCompatibilitySpinorQ[resolvedBasis];
  If[spinorQ === $Failed,
    Message[pointMatrix::input];
    Return[$Failed]
  ];
  seitzOperations = nativeSeitzOperations[symmetryInformation];
  spatialActions = nativeSpatialActions[seitzOperations];
  spinActions = nativeMSGSpinActions[spatialActions, evaluatedLattice];
  antiunitaryFlags =
    (TrueQ[Lookup[#, "Antiunitary", False]] & /@ seitzOperations);
  operationRecords = nativePointOperationRecords[
    spatialActions,
    spinActions,
    antiunitaryFlags,
    evaluatedLattice,
    {x, y, z},
    spinorQ
  ];
  matrices =
    FunctionBasisRepresentation`PointOperationRepresentationMatrices[
      resolvedBasis,
      operationRecords,
      {x, y, z}
    ];
  If[matrices === $Failed,
    Message[pointMatrix::compile];
    Return[$Failed]
  ];
  matrices
];
pointMatrix[___] := (Message[pointMatrix::input]; $Failed);

nativeUnconstrainedHopping[
    bondIndex_Integer?Positive,
    rowDimension_Integer?Positive,
    columnDimension_Integer?Positive
  ] := Table[
  ToExpression[
    "Global`tr" <> ToString[bondIndex] <>
      ToString[row] <> ToString[column]
  ] + I ToExpression[
    "Global`ti" <> ToString[bondIndex] <>
      ToString[row] <> ToString[column]
  ],
  {row, rowDimension},
  {column, columnDimension}
];

unsymham::shell =
  "The shell index must be a positive integer; received `1`.";
unsymham::shellrange =
  "Bond shell `1` was not prepared by init; only `2` shells are available. Rerun init with InitialBondShells -> `1`.";
unsymham::build =
  "The unconstrained Hamiltonian could not be built for shell `1`.";

unsymham[shell_Integer?Positive] := Module[
  {compiled, blockDimensions, terms, hamiltonian, availableShells},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  availableShells = Length[
    $CurrentModelSession["ModelSpecification", "BondClasses"]
  ];
  If[bondShellAvailableQ[shell] =!= True,
    Message[unsymham::shellrange, shell, availableShells];
    Return[$Failed]
  ];
  compiled = $CurrentModelSession["CompiledBondShells"][[shell]];
  If[!AssociationQ[compiled],
    Message[unsymham::build, shell];
    Return[$Failed]
  ];
  blockDimensions = compiled["BlockDimensions"];
  terms = Map[
    Function[bond,
      <|
        "BondIndex" -> bond["BondIndex"],
        "RowBlock" -> bond["RowSite"],
        "ColumnBlock" -> bond["ColumnSite"],
        "Displacement" -> bond["Displacement"],
        "Matrix" -> nativeUnconstrainedHopping[
          bond["BondIndex"],
          blockDimensions[[bond["RowSite"]]],
          blockDimensions[[bond["ColumnSite"]]]
        ]
      |>
    ],
    compiled["Bonds"]
  ];
  hamiltonian = MagneticTBLinearAlgebra`AssembleBlochHamiltonian[
    blockDimensions,
    terms,
    {kx, ky, kz}
  ];
  If[hamiltonian === $Failed,
    Message[unsymham::build, shell];
    Return[$Failed]
  ];
  hamiltonian
];
unsymham[shell_] := (Message[unsymham::shell, shell]; $Failed);

Options[symmetrizationHRInit] = {
  "Lattice" -> {{a, 0, 0}, {-(a/2), (Sqrt[3] a)/2, 0}, {0, 0, c}},
  "LattPar" -> {a -> 1, c -> 3},
  "WyckoffPosition" -> {{{2/3, 1/3, 0}, {0, 0, 1/2}}},
  "Symminformation" -> {{
    "1",
    IdentityMatrix[3],
    {0, 0, 0},
    "F"
  }},
  "BasisFunctions" -> {"s"},
  "Software" -> "VASP"
};

symmetrizationHRInit::compile =
  "The HR symmetrization data could not be prepared from the supplied lattice, sites, symmetry operations, and basis functions.";
symmetrizationHRInit::software =
  "Only the documented software convention \"VASP\" is currently supported; received `1`.";
symmetrizationHRInit::spin =
  "The VASP spinor reordering requires an even number of Wannier centers.";

symmetrizationHRInit[OptionsPattern[]] := Module[
  {software, input, compiled, data, atomPositions, pointRepresentations,
   operationCount, orbitMatrices, representationMatrices, centers,
   basisSpecification, resolvedFirstBasis, spinorQ, ordering},
  software = OptionValue["Software"];
  If[software =!= "VASP",
    Message[symmetrizationHRInit::software, software];
    Return[$Failed]
  ];
  input = <|
    "Lattice" -> OptionValue["Lattice"],
    "LatticeParameters" -> OptionValue["LattPar"],
    "WyckoffPosition" -> OptionValue["WyckoffPosition"],
    "SymmetryInformation" -> OptionValue["Symminformation"],
    "BasisFunctions" -> OptionValue["BasisFunctions"],
    "Debug" -> False,
    "InitialBondShells" -> 1,
    "RepresentationMode" -> "DirectProduct",
    "SiteLocalData" -> Automatic
  |>;
  compiled = CompileMagneticTBInput[input];
  If[compiled === $Failed,
    Message[symmetrizationHRInit::compile];
    Return[$Failed]
  ];
  data = compileCompatibilityProjectionData[compiled];
  If[data === $Failed,
    Message[symmetrizationHRInit::compile];
    Return[$Failed]
  ];
  atomPositions = data["AtomPositions"];
  pointRepresentations = data["PointRepresentations"];
  operationCount = Length[data["SymmetryInformation"]];
  representationMatrices = Table[
    orbitMatrices = Table[
      KroneckerProduct[
        IdentityMatrix[Length[atomPositions[[orbit]]]],
        pointRepresentations[[orbit, operation]]
      ],
      {orbit, Length[atomPositions]}
    ];
    nativeCompatibilityBlockDiagonal[orbitMatrices],
    {operation, operationCount}
  ];
  If[MemberQ[representationMatrices, $Failed],
    Message[symmetrizationHRInit::compile];
    Return[$Failed]
  ];
  centers = data["WannierCenters"];
  basisSpecification = data["BasisSpecification"];
  resolvedFirstBasis = resolveMagneticTBBasisEntry[
    basisSpecification[[1, 1]]
  ];
  spinorQ = ListQ[resolvedFirstBasis];
  If[spinorQ,
    If[OddQ[Length[centers]],
      Message[symmetrizationHRInit::spin];
      Return[$Failed]
    ];
    ordering = Join[
      Range[1, Length[centers], 2],
      Range[2, Length[centers], 2]
    ];
    centers = centers[[ordering]];
    representationMatrices =
      N[#[[ordering, ordering]]] & /@ representationMatrices
  ];
  <|
    "wcc" -> centers,
    "DR" -> representationMatrices,
    "symmetry" -> data["SymmetryInformation"]
  |>
];

End[]

EndPackage[]
