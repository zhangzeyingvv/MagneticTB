(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

origin = {IdentityMatrix[3], {0, 0, 0}};

Options[init] = {
  lattice -> {{a, 0, 0}, {-(a/2), (Sqrt[3] a)/2, 0}, {0, 0, c}},
  lattpar -> {a -> 1, c -> 3},
  wyckoffposition -> {{{2/3, 1/3, 0}, {0, 0, 1/2}}},
  symminformation -> {{
    "1",
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
    {0, 0, 0},
    "F"
  }},
  basisFunctions -> {"s"},
  debugQ -> False,
  MagneticTB`InitialBondShells -> 10,
  MagneticTB`GenerateSymmetryGroup -> False,
  MagneticTB`RepresentationMode -> "DirectProduct"
};

init::compile =
  "The model could not be prepared from the supplied init options; the current model session remains uninitialized.";
init::option = "Unknown init option(s): `1`.";
init::debug = "The debugQ option must be True or False; received `1`.";

init[opts : OptionsPattern[]] := Module[
  {providedOptionNames, allowedOptionNames, unknownOptionNames, debugValue,
   input, compiled},
  clearCurrentModelSession[];
  providedOptionNames = First /@ Cases[{opts}, _Rule | _RuleDelayed, {1}];
  allowedOptionNames = First /@ Options[init];
  unknownOptionNames = DeleteDuplicates@Select[
    providedOptionNames, !MemberQ[allowedOptionNames, #] &];
  If[unknownOptionNames =!= {},
    Message[init::option, unknownOptionNames]; Return[$Failed]
  ];
  debugValue = OptionValue[debugQ];
  If[!BooleanQ[debugValue],
    Message[init::debug, debugValue]; Return[$Failed]
  ];
  input = <|
    "Lattice" -> OptionValue[lattice],
    "LatticeParameters" -> OptionValue[lattpar],
    "WyckoffPosition" -> OptionValue[wyckoffposition],
    "SymmetryInformation" -> OptionValue[symminformation],
    "BasisFunctions" -> OptionValue[basisFunctions],
    "Debug" -> debugValue,
    "InitialBondShells" -> OptionValue[MagneticTB`InitialBondShells],
    "GenerateSymmetryGroup" ->
      OptionValue[MagneticTB`GenerateSymmetryGroup],
    "RepresentationMode" -> OptionValue[MagneticTB`RepresentationMode]
  |>;
  compiled = CompileMagneticTBInput[input];
  If[compiled === $Failed,
    Message[init::compile];
    Return[$Failed]
  ];
  If[installCurrentModelSession[compiled, input] === $Failed,
    Message[init::compile];
    Return[$Failed]
  ];
  Print[
    "Generators:",
    compiled["ModelSpecification", "Metadata", "SymmetryInformation"][[
      compiled["GeneratorIndices"],
      {1, 4}
    ]]
  ];
  Null
];

Options[initfromrep] = {
  lattice -> {{a, 0, 0}, {-(a/2), (Sqrt[3] a)/2, 0}, {0, 0, c}},
  lattpar -> {a -> 1, c -> 3},
  wyckoffposition -> {{{2/3, 1/3, 0}, {0, 0, 1/2}}},
  symminformation -> {{
    "1",
    IdentityMatrix[3],
    {0, 0, 0},
    "F"
  }},
  repinformation -> {{{1}}},
  orbitalLabels -> Automatic,
  debugQ -> False,
  MagneticTB`InitialBondShells -> 10,
  MagneticTB`RepresentationMode -> "DirectProduct"
};

initfromrep::option =
  "Unknown initfromrep option(s): `1`. basisFunctions, SiteLocalData, and GenerateSymmetryGroup are intentionally not supported.";
initfromrep::compile =
  "The complete-group representation input could not be prepared; the current model session remains uninitialized.";
initfromrep::debug =
  "The debugQ option must be True or False; received `1`.";

initfromrep[opts : OptionsPattern[]] := Module[
  {providedOptionNames, allowedOptionNames, unknownOptionNames,
   debugValue, input, compiled},
  clearCurrentModelSession[];
  providedOptionNames = First /@ Cases[
    {opts},
    _Rule | _RuleDelayed,
    {1}
  ];
  allowedOptionNames = First /@ Options[initfromrep];
  unknownOptionNames = DeleteDuplicates@Select[
    providedOptionNames,
    !MemberQ[allowedOptionNames, #] &
  ];
  If[unknownOptionNames =!= {},
    Message[initfromrep::option, unknownOptionNames];
    Return[$Failed]
  ];
  debugValue = OptionValue[debugQ];
  If[!BooleanQ[debugValue],
    Message[initfromrep::debug, debugValue]; Return[$Failed]
  ];
  input = <|
    "Lattice" -> OptionValue[lattice],
    "LatticeParameters" -> OptionValue[lattpar],
    "WyckoffPosition" -> OptionValue[wyckoffposition],
    "SymmetryInformation" -> OptionValue[symminformation],
    "RepresentationInformation" -> OptionValue[repinformation],
    "OrbitalLabels" -> OptionValue[orbitalLabels],
    "RepresentationSource" -> "Matrices",
    "Debug" -> debugValue,
    "InitialBondShells" -> OptionValue[MagneticTB`InitialBondShells],
    "RepresentationMode" -> OptionValue[MagneticTB`RepresentationMode]
  |>;
  compiled = CompileMagneticTBInput[input];
  If[compiled === $Failed,
    Message[initfromrep::compile];
    Return[$Failed]
  ];
  If[installCurrentModelSession[compiled, input] === $Failed,
    Message[initfromrep::compile];
    Return[$Failed]
  ];
  Print[
    "Generators:",
    compiled["ModelSpecification", "Metadata", "SymmetryInformation"][[
      compiled["GeneratorIndices"],
      {1, 4}
    ]]
  ];
  Null
];

Options[symham] = {
  "CartesianCoordinates" -> False,
  "ValidationLevel" -> "Basic",
  "KernelMethod" -> "Iterative",
  "Hermitian" -> True
};

symham::shell = "The shell index must be a positive integer; received `1`.";
symham::shellrange =
  "Bond shell `1` was not prepared by init; only `2` shells are available. Rerun init with InitialBondShells -> `1`.";
symham::hermitian =
  "The Hermitian option must be True or False; received `1`.";
symham::cartesian =
  "The CartesianCoordinates option must be True or False; received `1`.";
symham::option = "Unknown symham option(s): `1`.";
symham::build = "The linear-algebra Hamiltonian construction failed for shell `1`.";

shellParameterSymbols[shell_Integer?Positive, count_Integer?NonNegative] :=
  Table[
    ToExpression@Switch[
      shell,
      1, "e" <> ToString[index],
      2, "t" <> ToString[index],
      3, "r" <> ToString[index],
      4, "s" <> ToString[index],
      _, "p" <> ToString[shell] <> "n" <> ToString[index]
    ],
    {index, count}
  ];

realSpaceShellCacheKey[
    shell_Integer?Positive,
    hermitianQ : (True | False),
    method_String,
    validationLevel_String
  ] := ToString[
  {
    shell,
    hermitianQ,
    method,
    validationLevel
  },
  InputForm
];

(* Cache only the part of a compiled constraint problem consumed by the
   solver, Hamiltonian builder and real-space provenance layer.  The complete
   low-level CompileBondConstraints result is still returned to callers of
   that function; only the session cache is compact. *)
compactBondConstraintProblemForCache[problem_Association] := Join[
  KeyTake[
    problem,
    {"Schema", "SchemaVersion", "StaticShellKey", "Shell", "Hermitian"}
  ],
  <|
    "Orbits" -> Map[
      KeyTake[
        #,
        {
          "RepresentativeBondIndex",
          "Dimensions",
          "CoordinateDimension",
          "ConstraintBlocks",
          "ConstraintSources",
          "HasContinuousConstraint",
          "MemberRecords",
          "SourceDirectedOrbitIndices"
        }
      ] &,
      problem["Orbits"]
    ]
  |>
];

(* Validation matrices are transient.  Once SolveConstraintKernel has
   certified a basis, rebuilding hopping matrices needs only Dimensions,
   BasisMatrix and Nullity; the small rank/method/verification summary is kept
   for diagnostics. *)
compactSolvedShellForCache[solved_Association] := Join[
  KeyTake[
    solved,
    {"Schema", "SchemaVersion", "StaticShellKey", "Shell", "Hermitian"}
  ],
  <|
    "Solutions" -> Map[
      KeyTake[
        #,
        {
          "Dimensions",
          "BasisMatrix",
          "Rank",
          "Nullity",
          "KernelMethod",
          "ValidationLevel",
          "ResidualVerified",
          "IndependentVerified",
          "RankNullityVerified",
          "KernelVerified"
        }
      ] &,
      solved["Solutions"]
    ]
  |>
];

cacheRealSpaceShellResult[
    shell_Integer?Positive,
    hermitianQ : (True | False),
    method_String,
    validationLevel_String,
    constraintRecord_Association,
    built_Association,
    parameters_List,
    parameterRules_List
  ] := Module[{cache, key, record},
  cache = $CurrentModelSession["RealSpaceShellCache"];
  key = realSpaceShellCacheKey[
    shell,
    hermitianQ,
    method,
    validationLevel
  ];
  record = <|
    "Schema" -> "MagneticTBRealSpaceShellResult",
    "SchemaVersion" -> 3,
    "Shell" -> shell,
    "Hermitian" -> hermitianQ,
    "KernelMethod" -> method,
    "ValidationLevel" -> validationLevel,
    "Parameters" -> parameters,
    "ParametersByOrbit" ->
      (built["ParametersByOrbit"] /. parameterRules),
    "RepresentativeBondIndices" -> Lookup[
      constraintRecord["Data", "Orbits"],
      "RepresentativeBondIndex"
    ],
    "Terms" -> (built["Terms"] /. parameterRules)
  |>;
  AssociateTo[cache, key -> record];
  $CurrentModelSession = Join[
    $CurrentModelSession,
    <|"RealSpaceShellCache" -> cache|>
  ];
  record
];

cachedBondConstraint[
    shell_Integer?Positive,
    hermitianQ : (True | False)
  ] := Module[{cache, key, staticShell, problem},
  cache = $CurrentModelSession["BondConstraintCache"];
  key = ToString[
    {
      shell,
      hermitianQ
    },
    InputForm
  ];
  problem = Lookup[cache, key, Missing["NotCached"]];
  If[MissingQ[problem],
    staticShell = $CurrentModelSession["CompiledBondShells"][[shell]];
    problem = MagneticTBLinearAlgebra`CompileBondConstraints[
      staticShell,
      "Hermitian" -> hermitianQ
    ];
    If[problem === $Failed, Return[$Failed]];
    problem = compactBondConstraintProblemForCache[problem];
    AssociateTo[cache, key -> problem];
    $CurrentModelSession = Join[
      $CurrentModelSession,
      <|"BondConstraintCache" -> cache|>
    ]
  ];
  <|"Key" -> key, "Data" -> problem|>
];

cachedSolvedShell[
    constraintRecord_Association,
    method_String,
    validationLevel_String
  ] := Module[{cache, key, solved},
  cache = $CurrentModelSession["SolvedShellCache"];
  key = ToString[
    {constraintRecord["Key"], method, validationLevel},
    InputForm
  ];
  solved = Lookup[cache, key, Missing["NotCached"]];
  If[MissingQ[solved],
    solved = MagneticTBLinearAlgebra`SolveBondConstraintProblem[
      constraintRecord["Data"],
      "Method" -> method,
      "ValidationLevel" -> validationLevel
    ];
    If[solved === $Failed, Return[$Failed]];
    solved = compactSolvedShellForCache[solved];
    AssociateTo[cache, key -> solved];
    $CurrentModelSession = Join[
      $CurrentModelSession,
      <|"SolvedShellCache" -> cache|>
    ]
  ];
  solved
];

symham[shell_Integer?Positive, opts : OptionsPattern[]] := Module[
  {providedOptionNames, unknownOptionNames, allowedOptionNames,
   cartesianQ, validationLevel, method, hermitianQ, constraintRecord,
   staticShell, solved, built, parameters, parameterRules, hamiltonian,
   reciprocalLattice, availableShells},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  providedOptionNames = First /@ Cases[
    {opts},
    _Rule | _RuleDelayed,
    {1}
  ];
  allowedOptionNames = First /@ Options[symham];
  unknownOptionNames = Select[
    providedOptionNames,
    !MemberQ[allowedOptionNames, #] &
  ];
  If[unknownOptionNames =!= {},
    Message[symham::option, DeleteDuplicates[unknownOptionNames]];
    Return[$Failed]
  ];
  availableShells = Length[
    $CurrentModelSession["ModelSpecification", "BondClasses"]
  ];
  If[bondShellAvailableQ[shell] =!= True,
    Message[symham::shellrange, shell, availableShells];
    Return[$Failed]
  ];
  cartesianQ = OptionValue["CartesianCoordinates"];
  If[!BooleanQ[cartesianQ],
    Message[symham::cartesian, cartesianQ]; Return[$Failed]
  ];
  validationLevel = OptionValue["ValidationLevel"];
  method = OptionValue["KernelMethod"];
  hermitianQ = OptionValue["Hermitian"];
  If[!MemberQ[{True, False}, hermitianQ],
    Message[symham::hermitian, hermitianQ];
    Return[$Failed]
  ];
  staticShell = $CurrentModelSession["CompiledBondShells"][[shell]];
  constraintRecord = cachedBondConstraint[
    shell,
    hermitianQ
  ];
  If[constraintRecord === $Failed,
    Message[symham::build, shell];
    Return[$Failed]
  ];
  solved = cachedSolvedShell[
    constraintRecord,
    method,
    validationLevel
  ];
  If[solved === $Failed,
    Message[symham::build, shell];
    Return[$Failed]
  ];
  built = MagneticTBLinearAlgebra`BuildSolvedBondHamiltonian[
    staticShell,
    constraintRecord["Data"],
    solved,
    {kx, ky, kz},
    magneticTBInternalParameter
  ];
  If[built === $Failed,
    Message[symham::build, shell];
    Return[$Failed]
  ];
  parameters = shellParameterSymbols[shell, Length[built["Parameters"]]];
  parameterRules = If[
    parameters === {},
    {},
    Thread[built["Parameters"] -> parameters]
  ];
  hamiltonian = built["Hamiltonian"] /. parameterRules;
  cacheRealSpaceShellResult[
    shell,
    hermitianQ,
    method,
    validationLevel,
    constraintRecord,
    built,
    parameters,
    parameterRules
  ];
  If[cartesianQ,
    reciprocalLattice =
      $CurrentModelSession[
        "ModelSpecification", "Metadata", "ReciprocalLattice"
      ];
    hamiltonian = hamiltonian /. Thread[
      {kx, ky, kz} ->
        ({kx, ky, kz} (2 Pi)).Inverse[reciprocalLattice]
    ]
  ];
  Print["params:", parameters];
  hamiltonian
];

symham[shell_, OptionsPattern[]] :=
  (Message[symham::shell, shell]; $Failed);

End[]

EndPackage[]
