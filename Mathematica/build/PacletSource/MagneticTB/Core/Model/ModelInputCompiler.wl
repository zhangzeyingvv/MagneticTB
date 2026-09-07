(* ::Package:: *)

BeginPackage["MagneticTB`"]

CompileMagneticTBInput::usage =
  "CompileMagneticTBInput[input] compiles the unchanged init option syntax into a canonical model, verified representation data, and static directed-bond orbits.";

Begin["`Private`"]

ClearAll[
  nativeRationalizeInput,
  nativeReciprocalLattice
];

nativeRationalizeInput[expression_] := expression /. {
  value_?InexactNumberQ :> Rationalize[Round[value, .001]]
};

nativeReciprocalLattice[lattice_?MatrixQ] := Module[{volume},
  If[Dimensions[lattice] =!= {3, 3}, Return[$Failed]];
  volume = Cross[lattice[[1]], lattice[[2]]].lattice[[3]];
  If[TrueQ[PossibleZeroQ[volume]], Return[$Failed]];
  2 Pi/volume {
    Cross[lattice[[2]], lattice[[3]]],
    Cross[lattice[[3]], lattice[[1]]],
    Cross[lattice[[1]], lattice[[2]]]
  }
];

CompileMagneticTBInput::input =
  "The init options are incomplete or mutually incompatible. symminformation must be either traditional {label,R,t,\"F\"|\"T\"} records or a discrete spin-space-group label-><|\"space\"->{R,t},\"spin\"->{S,0|1}|> Association.";
CompileMagneticTBInput::lattice =
  "The evaluated lattice must be a real, nonsingular numerical 3 by 3 matrix.";
CompileMagneticTBInput::basis =
  "basisFunctions must contain one nonempty scalar or two-component spinor basis list per Wyckoff seed.";
CompileMagneticTBInput::representation =
  "Explicit repinformation is incompatible with the complete group, Wyckoff-orbit count, representation mode, or orbitalLabels. DirectProduct requires one full ordered matrix list per Wyckoff orbit; Induced requires one valid SiteLocalData Association per orbit.";
CompileMagneticTBInput::continuous =
  "The initfromrep continuous input is invalid. Exactly one internal C_infty record is allowed; its space action must be the identity, and Continuous SiteMatrices must be exact one-parameter unitary representations matching every compiled site dimension and normalizing the complete finite representation.";
CompileMagneticTBInput::symmetry =
  "The supplied symminformation generators could not be expanded by GenerateGroup.";
CompileMagneticTBInput::fullgroup =
  "initfromrep requires symminformation to be a complete closed finite group in the same order as repinformation; generator-only input is not accepted.";
CompileMagneticTBInput::generate =
  "GenerateSymmetryGroup must be True or False; received `1`.";
CompileMagneticTBInput::mode =
  "RepresentationMode must be \"DirectProduct\" or \"Induced\"; received `1`.";
CompileMagneticTBInput::induced =
  "The induced representation could not be compiled. SiteLocalData must be Automatic or one valid site-symmetry Association per Wyckoff orbit; invalid local data fail explicitly.";
CompileMagneticTBInput::bonds =
  "The static directed-bond orbit data could not be compiled for all bond shells requested by InitialBondShells.";

CompileMagneticTBInput[input_Association] := Module[
  {
    lattice, latticeParameters, latticePlot, reciprocalLattice,
    wyckoffSeeds, rawSymmetryInformation, suppliedSymmetryInformation,
    symmetryInformation, explicitSymmetryParts, continuousSymmetryInput,
    spinSpaceGroupQ, spinSpaceGroupInformation, spinSpaceGroupElements,
    generateSymmetryGroup, basisSpecification, representationSource,
    explicitRepresentationQ, rawRepresentationInformation,
    explicitRepresentationParts, representationInformation,
    continuousRepresentationInformation, continuousSymmetry,
    continuousCompatibility,
    suppliedOrbitalLabels, explicitDimensions,
    initialBondShells, representationMode, suppliedSiteLocalData,
    atomPositions, siteOrbits, sitePermutationData,
    localBases, variables, spinorFlags,
    pointOperationRecords, compiledLocalRepresentations,
    spinSpaceGroupBasisCompilations,
    symmetryCompilation, spatialActions, spinActions,
    groupAlgebra, antiunitaryFlags,
    bondSearch, bondClasses,
    discoveredShellCount,
    model, prepared, actionCompilation, localData, representation,
    generatorIndices, pointRepresentations,
    metadata, compiledBondShells,
    basisOrderingData
  },
  lattice = nativeRationalizeInput@Lookup[input, "Lattice", $Failed];
  latticeParameters = nativeRationalizeInput@Lookup[
    input,
    "LatticeParameters",
    $Failed
  ];
  wyckoffSeeds = nativeRationalizeInput@Lookup[
    input,
    "WyckoffPosition",
    $Failed
  ];
  rawSymmetryInformation = nativeRationalizeInput@Lookup[
    input,
    "SymmetryInformation",
    $Failed
  ];
  representationSource = Lookup[
    input,
    "RepresentationSource",
    "BasisFunctions"
  ];
  explicitRepresentationQ = representationSource === "Matrices";
  If[
    !explicitRepresentationQ && AssociationQ[rawSymmetryInformation] &&
      AnyTrue[
        Values[rawSymmetryInformation],
        AssociationQ[#] && KeyExistsQ[#, "continuous"] &
      ],
    Message[CompileMagneticTBInput::continuous];
    Return[$Failed]
  ];
  basisSpecification = If[
    explicitRepresentationQ,
    Missing["CompiledFromRepresentation"],
    Lookup[input, "BasisFunctions", $Failed]
  ];
  rawRepresentationInformation = Lookup[
    input,
    "RepresentationInformation",
    $Failed
  ];
  If[explicitRepresentationQ,
    explicitSymmetryParts = nativeSplitExplicitSymmetryInformation[
      rawSymmetryInformation
    ];
    explicitRepresentationParts =
      nativeSplitExplicitRepresentationInformation[
        rawRepresentationInformation
      ];
    If[
      explicitSymmetryParts === $Failed ||
        explicitRepresentationParts === $Failed,
      Message[CompileMagneticTBInput::continuous];
      Return[$Failed]
    ];
    suppliedSymmetryInformation = explicitSymmetryParts["Finite"];
    continuousSymmetryInput = explicitSymmetryParts["Continuous"];
    representationInformation =
      explicitRepresentationParts["Discrete"];
    continuousRepresentationInformation =
      explicitRepresentationParts["Continuous"],
    suppliedSymmetryInformation = rawSymmetryInformation;
    continuousSymmetryInput = None;
    representationInformation = rawRepresentationInformation;
    continuousRepresentationInformation = None
  ];
  suppliedOrbitalLabels = Lookup[input, "OrbitalLabels", Automatic];
  initialBondShells = Lookup[input, "InitialBondShells", 10];
  representationMode = Lookup[input, "RepresentationMode", "DirectProduct"];
  suppliedSiteLocalData = Lookup[input, "SiteLocalData", Automatic];
  generateSymmetryGroup = If[
    explicitRepresentationQ,
    False,
    Lookup[input, "GenerateSymmetryGroup", False]
  ];
  spinSpaceGroupQ = AssociationQ[suppliedSymmetryInformation];

  If[
    !explicitRepresentationQ &&
      ListQ[wyckoffSeeds] && Length[wyckoffSeeds] == 1 &&
      ListQ[basisSpecification] && basisSpecification =!= {} &&
      !ListQ[First[basisSpecification]],
    basisSpecification = {basisSpecification}
  ];

  If[
    !MatrixQ[lattice] || Dimensions[lattice] =!= {3, 3} ||
      !ListQ[latticeParameters] || !ListQ[wyckoffSeeds] ||
      !(ListQ[suppliedSymmetryInformation] || spinSpaceGroupQ) ||
      Length[suppliedSymmetryInformation] === 0 ||
      !MemberQ[{"BasisFunctions", "Matrices"}, representationSource] ||
      If[
        explicitRepresentationQ,
        !ListQ[representationInformation] ||
          Length[wyckoffSeeds] =!= Length[representationInformation],
        !ListQ[basisSpecification] ||
          Length[wyckoffSeeds] =!= Length[basisSpecification]
      ] ||
      !IntegerQ[initialBondShells] || initialBondShells < 1,
    Message[CompileMagneticTBInput::input];
    Return[$Failed]
  ];
  If[!MemberQ[{"DirectProduct", "Induced"}, representationMode],
    Message[CompileMagneticTBInput::mode, representationMode];
    Return[$Failed]
  ];
  If[!MemberQ[{True, False}, generateSymmetryGroup],
    Message[CompileMagneticTBInput::generate, generateSymmetryGroup];
    Return[$Failed]
  ];

  latticePlot = FullSimplify[lattice /. latticeParameters];
  reciprocalLattice = nativeReciprocalLattice[lattice];
  If[
    reciprocalLattice === $Failed ||
      !MatrixQ[latticePlot, NumericQ] || Dimensions[latticePlot] =!= {3, 3} ||
      !TrueQ[Chop[Det[N[latticePlot]]] != 0] ||
      Max[Abs[Im[N[latticePlot]]]] != 0.,
    Message[CompileMagneticTBInput::lattice];
    Return[$Failed]
  ];

  If[
    !And @@ (
      MatchQ[#, {_List, _List}] && Length[#[[1]]] == 3 &&
        Length[#[[2]]] == 3 & /@ wyckoffSeeds
    ) ||
      If[
        spinSpaceGroupQ,
        !And @@ (
          SymmetryAlgebra`ValidSpinSpaceGroupElementQ /@
            Values[suppliedSymmetryInformation]
        ),
        !And @@ (
          MatchQ[#, {_, _?MatrixQ, _List, "F" | "T"}] &&
            Dimensions[#[[2]]] === {3, 3} && Length[#[[3]]] == 3 & /@
            suppliedSymmetryInformation
        )
      ] ||
      If[
        explicitRepresentationQ,
        False,
        !And @@ (ListQ[#] && # =!= {} & /@ basisSpecification)
      ],
    Message[CompileMagneticTBInput::input];
    Return[$Failed]
  ];

  symmetryCompilation = If[
    spinSpaceGroupQ,
    compileSpinSpaceGroupInput[
      suppliedSymmetryInformation,
      latticePlot,
      generateSymmetryGroup
    ],
    compileMagneticGroupInput[
      suppliedSymmetryInformation,
      latticePlot,
      generateSymmetryGroup
    ]
  ];
  If[symmetryCompilation === $Failed,
    Message[
      If[
        explicitRepresentationQ,
        CompileMagneticTBInput::fullgroup,
        CompileMagneticTBInput::symmetry
      ]
    ];
    Return[$Failed]
  ];
  symmetryInformation = symmetryCompilation["SymmetryInformation"];
  groupAlgebra = symmetryCompilation["GroupAlgebra"];
  spatialActions = symmetryCompilation["SpatialActions"];
  spinActions = symmetryCompilation["SpinActions"];
  antiunitaryFlags = symmetryCompilation["AntiunitaryFlags"];
  If[spinSpaceGroupQ,
    spinSpaceGroupInformation =
      symmetryCompilation["SpinSpaceGroupInformation"];
    spinSpaceGroupElements = symmetryCompilation["SpinSpaceGroupElements"];
    atomPositions = nativeSpinSpaceGroupAtomPositions[
      wyckoffSeeds,
      spinSpaceGroupElements
    ],
    atomPositions = nativeAtomPositions[wyckoffSeeds, symmetryInformation]
  ];
  If[atomPositions === $Failed || MemberQ[atomPositions, $Failed, Infinity],
    Message[CompileMagneticTBInput::input];
    Return[$Failed]
  ];
  siteOrbits = atomPositions[[All, All, 1]];
  sitePermutationData = SitePermutationCompiler`CompileSitePermutationData[
    siteOrbits, spatialActions];
  If[sitePermutationData === $Failed,
    Message[CompileMagneticTBInput::input]; Return[$Failed]
  ];
  variables = {x, y, z};
  If[explicitRepresentationQ,
    explicitDimensions = nativeExplicitRepresentationDimensions[
      representationInformation,
      representationMode,
      Length[siteOrbits],
      Length[spatialActions]
    ];
    If[explicitDimensions === $Failed,
      Message[CompileMagneticTBInput::representation];
      Return[$Failed]
    ];
    basisSpecification = nativeExplicitOrbitalLabels[
      suppliedOrbitalLabels,
      explicitDimensions
    ];
    If[basisSpecification === $Failed,
      Message[CompileMagneticTBInput::representation];
      Return[$Failed]
    ];
    localBases = basisSpecification;
    pointOperationRecords = nativeAbstractPointOperationRecords[
      Length[siteOrbits],
      symmetryInformation[[All, 1]],
      antiunitaryFlags
    ];
    compiledLocalRepresentations = If[
      representationMode === "DirectProduct",
      representationInformation,
      Missing["NotDirectProduct"]
    ],
    localBases = Map[resolveMagneticTBBasisEntry, basisSpecification, {2}];
    spinorFlags = Map[
      Function[orbitBasis,
        If[
          And @@ (ListQ[#] && Length[#] == 2 & /@ orbitBasis),
          True,
          If[And @@ (!ListQ[#] & /@ orbitBasis), False, $Failed]
        ]
      ],
      localBases
    ];
    If[MemberQ[spinorFlags, $Failed],
      Message[CompileMagneticTBInput::basis];
      Return[$Failed]
    ];
    If[spinSpaceGroupQ,
      spinSpaceGroupBasisCompilations = Map[
        compileSpinSpaceGroupBasisAction[
          #,
          spinSpaceGroupElements,
          latticePlot,
          variables
        ] &,
        basisSpecification
      ];
      If[
        MemberQ[spinSpaceGroupBasisCompilations, $Failed] ||
          !And @@ (AssociationQ /@ spinSpaceGroupBasisCompilations),
        Message[CompileMagneticTBInput::basis];
        Return[$Failed]
      ];
      pointOperationRecords = Lookup[
        spinSpaceGroupBasisCompilations,
        "PointOperationRecords"
      ];
      compiledLocalRepresentations = Lookup[
        spinSpaceGroupBasisCompilations,
        "LocalMatrices"
      ],
      pointOperationRecords = Map[
        nativePointOperationRecords[
          spatialActions,
          spinActions,
          antiunitaryFlags,
          latticePlot,
          variables,
          #
        ] &,
        spinorFlags
      ];
      compiledLocalRepresentations = Missing["NotSupplied"]
    ]
  ];
  continuousSymmetry = If[
    explicitRepresentationQ,
    nativeCompileExplicitContinuousRepresentation[
      continuousSymmetryInput,
      continuousRepresentationInformation,
      siteOrbits,
      explicitDimensions
    ],
    None
  ];
  If[continuousSymmetry === $Failed,
    Message[CompileMagneticTBInput::continuous];
    Return[$Failed]
  ];
  bondSearch = MagneticTBLinearAlgebra`FindPeriodicBondClasses[
    siteOrbits,
    latticePlot,
    initialBondShells,
    "ReturnMetadata" -> True
  ];
  If[!AssociationQ[bondSearch], Return[$Failed]];
  discoveredShellCount = Length[bondSearch["BondClasses"]];
  bondClasses = Take[bondSearch["BondClasses"], initialBondShells];
  bondSearch = Join[
    bondSearch,
    <|
      "DiscoveredShellCount" -> discoveredShellCount,
      "ShellCount" -> Length[bondClasses],
      "BondClasses" -> bondClasses
    |>
  ];
  metadata = Join[
    <|
      "Source" -> "MagneticTBInitSyntax",
      "SymmetryInputMode" -> If[
        spinSpaceGroupQ,
        "DiscreteSpinSpaceGroup",
        "MagneticSpaceGroup"
      ],
      "SymbolicLattice" -> lattice,
      "LatticeParameters" -> latticeParameters,
      "LatticePlot" -> latticePlot,
      "ReciprocalLattice" -> reciprocalLattice,
      "WyckoffPosition" -> wyckoffSeeds,
      "AtomPositions" -> atomPositions,
      "BasisSpecification" -> basisSpecification,
      "SymmetryInformation" -> symmetryInformation,
      "RepresentationMode" -> representationMode,
      "RepresentationSource" -> representationSource,
      "ContinuousSymmetryLabel" -> If[
        AssociationQ[continuousSymmetry],
        continuousSymmetry["Label"],
        None
      ]
    |>,
    If[
      spinSpaceGroupQ,
      <|"SpinSpaceGroupInformation" -> spinSpaceGroupInformation|>,
      <||>
    ]
  ];
  model = ModelSchema`CreateModelSpecification[
    <|
      "Lattice" -> latticePlot,
      "SiteOrbits" -> siteOrbits,
      "LocalBases" -> localBases,
      "Variables" -> variables,
      "PointOperationRecords" -> pointOperationRecords,
      "Symmetry" -> <|
        "GroupAlgebra" -> groupAlgebra,
        "SpatialActions" -> spatialActions,
        "SpinActions" -> spinActions,
        "AntiunitaryFlags" -> antiunitaryFlags
      |>,
      "ContinuousSymmetry" -> continuousSymmetry,
      "BondClasses" -> bondClasses,
      "Metadata" -> metadata
    |>
  ];
  If[model === $Failed, Return[$Failed]];

  Switch[representationMode,
    "DirectProduct",
      prepared = If[
        explicitRepresentationQ || spinSpaceGroupQ,
        ModelPreparation`PrepareDirectProductModel[
          model,
          compiledLocalRepresentations,
          sitePermutationData
        ],
        ModelPreparation`PrepareDirectProductModel[model, sitePermutationData]
      ],
    "Induced",
      localData = If[
        explicitRepresentationQ,
        representationInformation,
        If[suppliedSiteLocalData === Automatic,
        nativeInducedLocalData[
          siteOrbits,
          localBases,
          pointOperationRecords,
          sitePermutationData,
          variables,
          groupAlgebra,
          antiunitaryFlags
        ],
        suppliedSiteLocalData
        ]
      ];
      If[
        localData === $Failed || !ListQ[localData] ||
          Length[localData] =!= Length[siteOrbits] ||
          !And @@ (AssociationQ /@ localData),
        Message[CompileMagneticTBInput::induced];
        Return[$Failed]
      ];
      actionCompilation = PhysicalRepresentation`CompileInducedActionData[
        siteOrbits,
        localData,
        spatialActions,
        groupAlgebra,
        antiunitaryFlags,
        sitePermutationData
      ];
      If[actionCompilation === $Failed,
        Message[CompileMagneticTBInput::induced];
        Return[$Failed]
      ];
      prepared = ModelSchema`CreatePreparedModel[model, actionCompilation]
  ];
  If[prepared === $Failed, Return[$Failed]];
  continuousCompatibility =
    nativeValidateExplicitContinuousCompatibility[
      prepared,
      continuousSymmetry
    ];
  If[continuousCompatibility === $Failed,
    Message[CompileMagneticTBInput::continuous];
    Return[$Failed]
  ];
  If[AssociationQ[continuousSymmetry],
    continuousSymmetry = Join[
      continuousSymmetry,
      continuousCompatibility
    ];
    model = Join[
      model,
      <|"ContinuousSymmetry" -> continuousSymmetry|>
    ];
    prepared = Join[
      prepared,
      <|
        "ModelSpecification" -> model,
        "ContinuousSymmetry" -> continuousSymmetry
      |>
    ]
  ];
  basisOrderingData = If[
    representationMode === "Induced",
    If[
      explicitRepresentationQ,
      nativeAbstractInducedBasisOrderingData[
        actionCompilation,
        basisSpecification
      ],
      nativeInducedBasisOrderingData[
        actionCompilation,
        localBases,
        pointOperationRecords,
        variables
      ]
    ],
    <|"Mode" -> "DirectProduct"|>
  ];
  If[basisOrderingData === $Failed, Return[$Failed]];
  representation = ModelPreparation`AssemblePreparedRepresentation[prepared];
  If[representation === $Failed, Return[$Failed]];

  generatorIndices = GroupAlgebra`FindGeneratorIndices[groupAlgebra];
  If[generatorIndices === $Failed, Return[$Failed]];
  pointRepresentations = If[
    representationMode === "DirectProduct",
    prepared["LocalRepresentationMatrices"],
    Lookup[representation, "OrbitRepresentationMatrices", Missing["NotAvailable"]]
  ];
  (* The high-level compiler result must not retain PreparedModel or action
     compilation copies after static bond-orbit compilation has consumed them. *)
  representation = KeyTake[
    representation,
    {
      "Method",
      "SpatialActions",
      "SpinActions",
      "GroupAlgebra",
      "Convention",
      "AntiunitaryFlags",
      "RepresentationMatrices",
      "Dimension",
      "UnitaryVerified"
    }
  ];
  compiledBondShells =
    MagneticTBLinearAlgebra`CompileDirectedBondOrbits[
      prepared,
      Range[Length[bondClasses]]
    ];
  If[compiledBondShells === $Failed ||
      MemberQ[compiledBondShells, $Failed, Infinity],
    Message[CompileMagneticTBInput::bonds];
    Return[$Failed]
  ];
  <|
    "ModelSpecification" -> model,
    "FullRepresentation" -> representation,
    "RepresentationMode" -> representationMode,
    "RepresentationSource" -> representationSource,
    "BasisOrderingData" -> basisOrderingData,
    "GeneratorIndices" -> generatorIndices,
    (* Used only while installing traditional projected symbols or preparing
       the compatibility HR interface; it is not retained in the session. *)
    "PointRepresentations" -> pointRepresentations,
    (* BondClasses already live in ModelSpecification.  Keep only
       adaptive-search metadata here. *)
    "BondSearch" -> KeyDrop[bondSearch, "BondClasses"],
    "CompiledBondShells" -> compiledBondShells
  |>
];

CompileMagneticTBInput[_] :=
  (Message[CompileMagneticTBInput::input]; $Failed);

End[]

EndPackage[]
