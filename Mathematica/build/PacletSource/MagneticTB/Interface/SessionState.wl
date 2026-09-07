(* ::Package:: *)

BeginPackage["MagneticTB`"]

CurrentModelSession::usage =
  "CurrentModelSession[] returns the read-only Association prepared by init, including the model, compact verified representation data, static CompiledBondShells, solver caches, and real-space shell results.";

Begin["`Private`"]

If[
  !ValueQ[$CurrentModelSession],
  $CurrentModelSession = Missing["NotInitialized"]
];

CurrentModelSession::notinit =
  "MagneticTB has not been initialized. Call init[...] or initfromrep[...] first.";

CurrentModelSession[] := If[
  AssociationQ[$CurrentModelSession],
  $CurrentModelSession,
  Message[CurrentModelSession::notinit];
  $Failed
];

clearCurrentModelSession[] := Module[{},
  $CurrentModelSession = Missing["NotInitialized"];
  Clear[
    latt,
    latpar,
    lattplot,
    reclatt,
    wyckoff,
    symminfo,
    ops,
    atompos,
    basis,
    pointops,
    symmetryops,
    symmcompile,
    bondclassify,
    wcc,
    genindex,
    matrixRep,
    tran,
    timere
  ];
  Null
];

(* The session keeps only the full matrices and their exact unitarity result.
   Geometry and abstract group data have one owner: ModelSpecification. *)
compactFullRepresentationForSession[representation_Association] := KeyTake[
  representation,
  {
    "Method",
    "Convention",
    "RepresentationMatrices",
    "Dimension",
    "UnitaryVerified"
  }
];

(* ImageTable, SymmetryIndices and BondToDirectedOrbit are compilation audit
   intermediates.  Once DirectedOrbits and ReverseIndices have been verified,
   neither symham nor any read-only/report/IO consumer reads them. *)
compactBondShellForSession[shell_Association] := KeyDrop[
  shell,
  {"SymmetryIndices", "ImageTable", "BondToDirectedOrbit"}
];

compatibilityWannierCenters[
    siteOrbits_List,
    localDimensions_List
  ] := Flatten[
  MapThread[
    Function[{orbit, dimension},
      Flatten[ConstantArray[#, dimension] & /@ orbit, 1]
    ],
    {siteOrbits, localDimensions}
  ],
  1
];

compatibilitySymmetryCompile[
    symmetryInformation_List,
    representationMatrices_List,
    lattice_?MatrixQ,
    reciprocalLattice_?MatrixQ,
    latticeParameters_List
  ] := Table[
  {
    operationIndex,
    symmetryInformation[[operationIndex]],
    representationMatrices[[operationIndex]],
    (
      Inverse[Transpose[reciprocalLattice]] .
        ((Transpose[lattice].symmetryInformation[[operationIndex, 2]]) .
          Inverse[Transpose[lattice]]) .
        Transpose[reciprocalLattice]
    ) /. latticeParameters
  },
  {operationIndex, Length[symmetryInformation]}
];

(* CompatibilityData is built at the interface boundary from canonical
   compiler output.  It is consumed by legacy symbol projection or the HR
   compatibility helper and is never retained in the native session. *)
compileCompatibilityProjectionData[compiled_Association] := Module[
  {
    model, metadata, representation, symmetryInformation,
    representationMatrices, pointRepresentations, lattice,
    reciprocalLattice, latticeParameters, centers
  },
  model = Lookup[compiled, "ModelSpecification", $Failed];
  representation = Lookup[compiled, "FullRepresentation", $Failed];
  pointRepresentations = Lookup[compiled, "PointRepresentations", $Failed];
  If[!AssociationQ[model] || !AssociationQ[representation], Return[$Failed]];
  metadata = Lookup[model, "Metadata", $Failed];
  If[!AssociationQ[metadata], Return[$Failed]];
  symmetryInformation = Lookup[metadata, "SymmetryInformation", $Failed];
  representationMatrices = Lookup[
    representation,
    "RepresentationMatrices",
    $Failed
  ];
  lattice = Lookup[metadata, "SymbolicLattice", $Failed];
  reciprocalLattice = Lookup[metadata, "ReciprocalLattice", $Failed];
  latticeParameters = Lookup[metadata, "LatticeParameters", $Failed];
  If[
    !ListQ[symmetryInformation] || symmetryInformation === {} ||
      !ListQ[representationMatrices] ||
      Length[representationMatrices] =!= Length[symmetryInformation] ||
      !ListQ[pointRepresentations] ||
      !MatrixQ[lattice] || !MatrixQ[reciprocalLattice] ||
      !ListQ[latticeParameters],
    Return[$Failed]
  ];
  centers = compatibilityWannierCenters[
    model["SiteOrbits"],
    model["LocalDimensions"]
  ];
  Join[
    <|
      "Lattice" -> lattice,
      "LatticeParameters" -> latticeParameters,
      "LatticePlot" -> metadata["LatticePlot"],
      "ReciprocalLattice" -> reciprocalLattice,
      "WyckoffPosition" -> metadata["WyckoffPosition"],
      "SymmetryInformation" -> symmetryInformation,
      "OperationLabels" -> symmetryInformation[[All, 1]],
      "AtomPositions" -> metadata["AtomPositions"],
      "BasisSpecification" -> metadata["BasisSpecification"],
      "PointRepresentations" -> pointRepresentations,
      "RepresentationMatrices" -> representationMatrices,
      "SymmetryCompile" -> compatibilitySymmetryCompile[
        symmetryInformation,
        representationMatrices,
        lattice,
        reciprocalLattice,
        latticeParameters
      ],
      "BondClasses" -> model["BondClasses"],
      "WannierCenters" -> centers
    |>,
    KeyTake[metadata, "SpinSpaceGroupInformation"]
  ]
];
compileCompatibilityProjectionData[_] := $Failed;

installCompatibilityProjection[
    data_Association,
    generatorIndices_List
  ] := Module[{labels, symmetryInformation},
  symmetryInformation = data["SymmetryInformation"];
  labels = data["OperationLabels"];

  latt = data["Lattice"];
  latpar = data["LatticeParameters"];
  lattplot = data["LatticePlot"];
  reclatt = data["ReciprocalLattice"];
  wyckoff = data["WyckoffPosition"];
  symminfo = symmetryInformation;
  ops = labels;
  atompos = data["AtomPositions"];
  basis = data["BasisSpecification"];
  pointops = data["PointRepresentations"];
  symmetryops = data["RepresentationMatrices"];
  symmcompile = data["SymmetryCompile"];
  bondclassify = data["BondClasses"];
  wcc = data["WannierCenters"];
  genindex = generatorIndices;

  Clear[matrixRep, tran, timere];
  Do[
    matrixRep[labels[[index]]] = symmetryInformation[[index, 2]];
    tran[labels[[index]]] = symmetryInformation[[index, 3]];
    timere[labels[[index]]] = symmetryInformation[[index, 4]],
    {index, Length[labels]}
  ];
];

installCurrentModelSession[
    compiled_Association,
    initInput_Association
  ] := Module[{storedInitInput, compatibilityData},
  compatibilityData = compileCompatibilityProjectionData[compiled];
  If[compatibilityData === $Failed, Return[$Failed]];
  storedInitInput = If[
    Lookup[compiled, "RepresentationSource", "BasisFunctions"] ===
      "Matrices",
    KeyDrop[initInput, "RepresentationInformation"],
    initInput
  ];
  $CurrentModelSession = <|
    "Schema" -> "MagneticTBModelSession",
    "SchemaVersion" -> 7,
    "InitInput" -> storedInitInput,
    "RepresentationMode" -> compiled["RepresentationMode"],
    "RepresentationSource" -> compiled["RepresentationSource"],
    "BasisOrderingData" -> compiled["BasisOrderingData"],
    "ModelSpecification" -> compiled["ModelSpecification"],
    "FullRepresentation" -> compactFullRepresentationForSession[
      compiled["FullRepresentation"]
    ],
    "GeneratorIndices" -> compiled["GeneratorIndices"],
    "BondSearch" -> compiled["BondSearch"],
    "CompiledBondShells" ->
      (compactBondShellForSession /@ compiled["CompiledBondShells"]),
    "BondConstraintCache" -> <||>,
    "SolvedShellCache" -> <||>,
    "RealSpaceShellCache" -> <||>
  |>;
  installCompatibilityProjection[
    compatibilityData,
    compiled["GeneratorIndices"]
  ];
  Null
];

ensureCurrentModelSession[] := If[
  AssociationQ[$CurrentModelSession],
  True,
  Message[CurrentModelSession::notinit];
  False
];

bondShellAvailableQ[requestedShells_Integer?Positive] := Module[
  {currentCount},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  currentCount = Length[
    $CurrentModelSession["ModelSpecification", "BondClasses"]
  ];
  currentCount >= requestedShells
];

End[]

EndPackage[]
