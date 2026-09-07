(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  validHamiltonianBasisSiteQ,
  analyzeOrbitalBasisState,
  buildHamiltonianBasisOrderData,
  renderOrbitalTableGrid
];

orbitalTable::data =
  "The orbital table stored by init is missing or inconsistent.";

validHamiltonianBasisSiteQ[site_Association] :=
  And @@ (KeyExistsQ[site, #] & /@ {
    "SiteIndex",
    "OrbitIndex",
    "EquivalentIndex",
    "Position",
    "Dimension",
    "BlockRange"
  }) &&
    IntegerQ[site["SiteIndex"]] && site["SiteIndex"] > 0 &&
    IntegerQ[site["OrbitIndex"]] && site["OrbitIndex"] > 0 &&
    IntegerQ[site["EquivalentIndex"]] &&
      site["EquivalentIndex"] > 0 &&
    ListQ[site["Position"]] && Length[site["Position"]] == 3 &&
    IntegerQ[site["Dimension"]] && site["Dimension"] > 0 &&
    ListQ[site["BlockRange"]] &&
      Length[site["BlockRange"]] == site["Dimension"] &&
    And @@ (IntegerQ[#] && # > 0 & /@ site["BlockRange"]);

analyzeOrbitalBasisState[basisState_, variables_List] := Module[
  {
    nonzeroPosition, spatialOrbital, spinState,
    coordinatePattern, factorizedQ
  },
  If[
    !ListQ[basisState],
    Return[<|
      "SpatialOrbital" -> basisState,
      "SpinState" -> None,
      "SpinStructure" -> "Spinless"
    |>]
  ];
  If[Length[basisState] =!= 2,
    Return[<|
      "SpatialOrbital" -> Missing["EncodedInInternalState"],
      "SpinState" -> basisState,
      "SpinStructure" -> "GeneralInternalState"
    |>]
  ];
  nonzeroPosition = SelectFirst[
    Range[Length[basisState]],
    !MatrixPredicates`ExactZeroExpressionQ[
      basisState[[#]]
    ] &,
    Missing["NotFound"]
  ];
  If[MissingQ[nonzeroPosition],
    Return[<|
      "SpatialOrbital" -> Missing["ZeroSpinor"],
      "SpinState" -> basisState,
      "SpinStructure" -> "SpatialSpinor"
    |>]
  ];
  spatialOrbital = basisState[[nonzeroPosition]];
  spinState = FullSimplify[
    Cancel@Together[#/spatialOrbital] & /@ basisState
  ];
  coordinatePattern = Alternatives @@ variables;
  factorizedQ = FreeQ[spinState, coordinatePattern] &&
    And @@ (
      MatrixPredicates`ExactZeroExpressionQ /@
        FullSimplify[basisState - spatialOrbital spinState]
    );
  If[
    TrueQ[factorizedQ],
    <|
      "SpatialOrbital" -> spatialOrbital,
      "SpinState" -> spinState,
      "SpinStructure" -> "FactorizedSpinor"
    |>,
    <|
      "SpatialOrbital" -> Missing["EncodedInSpinor"],
      "SpinState" -> basisState,
      "SpinStructure" -> "SpatialSpinor"
    |>
  ]
];

buildHamiltonianBasisOrderData[session_Association] := Module[
  {
    compiledShells, staticShell, sites, metadata,
    basisSpecification, representationMode, basisOrderingData,
    representationSource,
    symmetryInformation, operationLabels, inducedOrbitData,
    latticePlot, modelSpecification, variables,
    records, failure, indices
  },
  compiledShells = Lookup[session, "CompiledBondShells", $Failed];
  modelSpecification = Lookup[session, "ModelSpecification", $Failed];
  metadata = If[
    AssociationQ[modelSpecification],
    Lookup[modelSpecification, "Metadata", $Failed],
    $Failed
  ];
  If[
    !ListQ[compiledShells] || compiledShells === {} ||
      !AssociationQ[modelSpecification] || !AssociationQ[metadata],
    Return[$Failed]
  ];
  staticShell = First[compiledShells];
  sites = Lookup[staticShell, "Sites", $Failed];
  basisSpecification = Lookup[
    metadata,
    "BasisSpecification",
    $Failed
  ];
  latticePlot = Lookup[metadata, "LatticePlot", $Failed];
  variables = Lookup[modelSpecification, "Variables", $Failed];
  representationMode = Lookup[session, "RepresentationMode", $Failed];
  representationSource = Lookup[
    session,
    "RepresentationSource",
    "BasisFunctions"
  ];
  basisOrderingData = Lookup[session, "BasisOrderingData", $Failed];
  symmetryInformation = Lookup[
    metadata,
    "SymmetryInformation",
    $Failed
  ];
  If[
    !ListQ[sites] || sites === {} ||
      !And @@ (validHamiltonianBasisSiteQ /@ sites) ||
      !ListQ[basisSpecification] || basisSpecification === {} ||
      !MatrixQ[latticePlot] || Dimensions[latticePlot] =!= {3, 3} ||
      !ListQ[variables] || variables === {} ||
      !MemberQ[{"DirectProduct", "Induced"}, representationMode] ||
      !MemberQ[{"BasisFunctions", "Matrices"}, representationSource] ||
      !AssociationQ[basisOrderingData] ||
      Lookup[basisOrderingData, "Mode", $Failed] =!=
        representationMode ||
      !ListQ[symmetryInformation] || symmetryInformation === {},
    Return[$Failed]
  ];
  operationLabels = symmetryInformation[[All, 1]];
  inducedOrbitData = If[
    representationMode === "Induced",
    Lookup[basisOrderingData, "OrbitData", $Failed],
    {}
  ];
  If[
    representationMode === "Induced" &&
      (
        !ListQ[inducedOrbitData] ||
        Length[inducedOrbitData] =!= Length[basisSpecification] ||
        !And @@ (AssociationQ /@ inducedOrbitData)
      ),
    Return[$Failed]
  ];
  records = Flatten@Map[
    Function[site,
      Module[
        {
          orbitIndex, localBasis, referenceSiteIndex,
          transportOperationIndices, transportOperationIndex,
          transportOperationLabel, basisConvention,
          transportedBases, actualLocalBasis
        },
        orbitIndex = site["OrbitIndex"];
        If[
          !Between[orbitIndex, {1, Length[basisSpecification]}],
          Return[$Failed, Module]
        ];
        localBasis = basisSpecification[[orbitIndex]];
        If[
          !ListQ[localBasis] ||
            Length[localBasis] =!= site["Dimension"],
          Return[$Failed, Module]
        ];
        If[representationMode === "Induced",
          referenceSiteIndex = Lookup[
            inducedOrbitData[[orbitIndex]],
            "ReferenceSiteIndex",
            $Failed
          ];
          transportOperationIndices = Lookup[
            inducedOrbitData[[orbitIndex]],
            "CosetRepresentativeIndices",
            $Failed
          ];
          If[
            !IntegerQ[referenceSiteIndex] ||
              !ListQ[transportOperationIndices] ||
              !Between[
                site["EquivalentIndex"],
                {1, Length[transportOperationIndices]}
              ],
            Return[$Failed, Module]
          ];
          transportOperationIndex =
            transportOperationIndices[[site["EquivalentIndex"]]];
          If[
            !IntegerQ[transportOperationIndex] ||
              !Between[
                transportOperationIndex,
                {1, Length[operationLabels]}
              ],
            Return[$Failed, Module]
          ];
          transportOperationLabel =
            operationLabels[[transportOperationIndex]];
          transportedBases = Lookup[
            inducedOrbitData[[orbitIndex]],
            "TransportedBases",
            $Failed
          ];
          If[
            !ListQ[transportedBases] ||
              !Between[
                site["EquivalentIndex"],
                {1, Length[transportedBases]}
              ] ||
              !ListQ[
                transportedBases[[site["EquivalentIndex"]]]
              ] ||
              Length[transportedBases[[site["EquivalentIndex"]]]] =!=
                site["Dimension"],
            Return[$Failed, Module]
          ];
          actualLocalBasis =
            transportedBases[[site["EquivalentIndex"]]];
          basisConvention = "Induced",
          referenceSiteIndex = None;
          transportOperationIndex = None;
          transportOperationLabel = None;
          actualLocalBasis = If[
            representationSource === "Matrices",
            localBasis,
            resolveMagneticTBBasisEntry /@ localBasis
          ];
          basisConvention = "DirectProduct"
        ];
        MapThread[
          Function[{
            hamiltonianIndex,
            localOrbitalIndex,
            basisFunction,
            transportedBasisFunction
          },
            Join[
              <|
                "OrbitalID" -> hamiltonianIndex,
                "HamiltonianIndex" -> hamiltonianIndex,
                "SiteIndex" -> site["SiteIndex"],
                "WyckoffOrbitIndex" -> orbitIndex,
                "EquivalentAtomIndex" -> site["EquivalentIndex"],
                "FractionalPosition" -> site["Position"],
                "CartesianPosition" ->
                  Simplify[site["Position"].latticePlot],
                "LocalBasisIndex" -> localOrbitalIndex,
                "LocalOrbitalIndex" -> localOrbitalIndex,
                "BasisFunctionInput" -> basisFunction,
                "BasisFunction" -> basisFunction,
                "BasisState" -> transportedBasisFunction,
                "TransportedBasisFunction" -> transportedBasisFunction,
                "BasisConvention" -> basisConvention,
                "ReferenceEquivalentAtomIndex" -> referenceSiteIndex,
                "TransportOperationIndex" -> transportOperationIndex,
                "TransportOperationLabel" -> transportOperationLabel
              |>,
              analyzeOrbitalBasisState[
                transportedBasisFunction,
                variables
              ]
            ]
          ],
          {
            site["BlockRange"],
            Range[site["Dimension"]],
            localBasis,
            actualLocalBasis
          }
        ]
      ]
    ],
    sites
  ];
  If[MemberQ[records, $Failed], Return[$Failed]];
  failure = SelectFirst[records, !AssociationQ[#] &, Missing["NotFound"]];
  If[!MissingQ[failure], Return[$Failed]];
  indices = Lookup[records, "OrbitalID", $Failed];
  If[
    indices =!= Range[Length[records]] ||
      Lookup[records, "HamiltonianIndex", $Failed] =!= indices ||
      Length[records] =!= Total[Lookup[sites, "Dimension"]],
    Return[$Failed]
  ];
  <|
    "Schema" -> "MagneticTBHamiltonianBasisOrder",
    "SchemaVersion" -> 4,
    "Dimension" -> Length[records],
    "RepresentationMode" -> representationMode,
    "Convention" -> "H[[i,j]] = <basis i|H|basis j>",
    "Orbitals" -> records
  |>
];

renderOrbitalTableGrid[data_Association] := Grid[
  Join[
    {{
      Row[{
        "Hamiltonian dimension = ", data["Dimension"],
        "; representation mode = ", data["RepresentationMode"]
      }],
      SpanFromLeft, SpanFromLeft, SpanFromLeft, SpanFromLeft,
      SpanFromLeft, SpanFromLeft, SpanFromLeft, SpanFromLeft,
      SpanFromLeft, SpanFromLeft, SpanFromLeft
    }},
    {{
      "H index",
      "Site",
      "Wyckoff orbit",
      "Atom",
      "Fractional position",
      "Cartesian position",
      "Local orbital",
      "Input basis",
      "Actual basis state",
      "Spin structure",
      "Spin state",
      "Transport operation"
    }},
    Map[
      {
        #["HamiltonianIndex"],
        #["SiteIndex"],
        #["WyckoffOrbitIndex"],
        #["EquivalentAtomIndex"],
        #["FractionalPosition"],
        #["CartesianPosition"],
        #["LocalOrbitalIndex"],
        #["BasisFunctionInput"],
        #["BasisState"],
        #["SpinStructure"],
        #["SpinState"],
        If[
          #["TransportOperationIndex"] === None,
          None,
          Row[{
            #["TransportOperationIndex"],
            ": ",
            #["TransportOperationLabel"]
          }]
        ]
      } &,
      data["Orbitals"]
    ]
  ],
  Frame -> All,
  Alignment -> Left
];

orbitalTable[] := Module[{data},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  data = buildHamiltonianBasisOrderData[$CurrentModelSession];
  If[data === $Failed,
    Message[orbitalTable::data];
    Return[$Failed]
  ];
  renderOrbitalTableGrid[data]
];

End[]

EndPackage[]
