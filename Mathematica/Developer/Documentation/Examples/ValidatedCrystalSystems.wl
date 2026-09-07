Needs["MagneticTB`"];

ClearAll[
  wyckoffLetterName,
  wyckoffDatabaseSeed,
  wyckoffLattice,
  wyckoffBandPath,
  wyckoffParameterRules,
  runWyckoffMode
];

(* Read the formal database shipped by the installed paclet. *)
magneticTBPackageRoot = DirectoryName@DirectoryName@FindFile["MagneticTB`"];
wyckoffDatabase = Import[
  FileNameJoin[{magneticTBPackageRoot, "Data", "wyckoffMSG.mx"}]
];

wyckoffLetterName[value_] := If[
  StringQ[value],
  value,
  SymbolName[Unevaluated[value]]
];

wyckoffDatabaseSeed[key_List, letter_String] := Module[{entry},
  entry = SelectFirst[
    wyckoffDatabase[key],
    wyckoffLetterName[#[[3]]] === letter &
  ];
  entry[[4, 1]] /. {
    x -> 1/7, y -> 2/11, z -> 3/13,
    mx -> 1/5, my -> 2/7, mz -> 3/11
  }
];

wyckoffLattice["Triclinic"] := {
  {a, 0, 0},
  {b Cos[gamma], b Sin[gamma], 0},
  {
    c Cos[beta],
    c (Cos[alpha] - Cos[beta] Cos[gamma]) Csc[gamma],
    c Sqrt[
      1 - Cos[alpha]^2 - Cos[beta]^2 +
        2 Cos[alpha] Cos[beta] Cos[gamma] - Cos[gamma]^2
    ] Csc[gamma]
  }
};
wyckoffLattice["Monoclinic"] := {
  {a, 0, 0},
  {0, 0, b},
  {c Cos[beta], c Sin[beta], 0}
};
wyckoffLattice["Orthorhombic"] := {
  {a, 0, 0}, {0, b, 0}, {0, 0, c}
};
wyckoffLattice["Tetragonal"] := {
  {a, 0, 0}, {0, a, 0}, {0, 0, c}
};
wyckoffLattice["Trigonal"] := {
  {a, 0, 0}, {-a/2, Sqrt[3] a/2, 0}, {0, 0, c}
};
wyckoffLattice["Hexagonal"] := wyckoffLattice["Trigonal"];
wyckoffLattice["Cubic"] := {
  {a, 0, 0}, {0, a, 0}, {0, 0, a}
};

wyckoffLatticeParameters["Triclinic"] := {
  a -> 1, b -> 6/5, c -> 7/5,
  alpha -> Pi/2, beta -> Pi/3, gamma -> Pi/4
};
wyckoffLatticeParameters["Monoclinic"] := {
  a -> 1, b -> 6/5, c -> 7/5, beta -> Pi/3
};
wyckoffLatticeParameters["Orthorhombic"] := {
  a -> 1, b -> 6/5, c -> 7/5
};
wyckoffLatticeParameters["Tetragonal" | "Trigonal" | "Hexagonal"] := {
  a -> 1, c -> 3/2
};
wyckoffLatticeParameters["Cubic"] := {a -> 1};

wyckoffBandPath = {
  {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
  {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
  {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
};

wyckoffParameterRules[hamiltonian_] := Module[{parameters},
  parameters = SortBy[
    DeleteDuplicates@Cases[
      hamiltonian,
      symbol_Symbol /; StringMatchQ[
        SymbolName[symbol], ("e" | "t") ~~ DigitCharacter ..
      ],
      Infinity
    ],
    SymbolName
  ];
  Thread[parameters -> N[Range[Length[parameters]]/(Length[parameters] + 2)]]
];

wyckoffCases = {
  <|"ID" -> "1.3-a", "CrystalSystem" -> "Triclinic",
    "Group" -> typeIV[{1, 3}], "GroupCode" -> "typeIV[{1, 3}]",
    "Key" -> {1, 3}, "Letter" -> "a"|>,
  <|"ID" -> "2.6-i", "CrystalSystem" -> "Triclinic",
    "Group" -> typeIII[{2, 6}], "GroupCode" -> "typeIII[{2, 6}]",
    "Key" -> {2, 6}, "Letter" -> "i"|>,
  <|"ID" -> "3.3-e", "CrystalSystem" -> "Monoclinic",
    "Group" -> typeIII[{3, 3}], "GroupCode" -> "typeIII[{3, 3}]",
    "Key" -> {3, 3}, "Letter" -> "e"|>,
  <|"ID" -> "3.4-a", "CrystalSystem" -> "Monoclinic",
    "Group" -> typeIV[{3, 4}], "GroupCode" -> "typeIV[{3, 4}]",
    "Key" -> {3, 4}, "Letter" -> "a"|>,
  <|"ID" -> "16.3-q", "CrystalSystem" -> "Orthorhombic",
    "Group" -> typeIII[{16, 3}], "GroupCode" -> "typeIII[{16, 3}]",
    "Key" -> {16, 3}, "Letter" -> "q"|>,
  <|"ID" -> "16.3-r", "CrystalSystem" -> "Orthorhombic",
    "Group" -> typeIII[{16, 3}], "GroupCode" -> "typeIII[{16, 3}]",
    "Key" -> {16, 3}, "Letter" -> "r"|>,
  <|"ID" -> "75.3-c", "CrystalSystem" -> "Tetragonal",
    "Group" -> typeIII[{75, 3}], "GroupCode" -> "typeIII[{75, 3}]",
    "Key" -> {75, 3}, "Letter" -> "c"|>,
  <|"ID" -> "143.3-a", "CrystalSystem" -> "Trigonal",
    "Group" -> typeIV[{143, 3}], "GroupCode" -> "typeIV[{143, 3}]",
    "Key" -> {143, 3}, "Letter" -> "a"|>,
  <|"ID" -> "168.111-b", "CrystalSystem" -> "Hexagonal",
    "Group" -> typeIII[{168, 111}], "GroupCode" -> "typeIII[{168, 111}]",
    "Key" -> {168, 111}, "Letter" -> "b"|>,
  <|"ID" -> "195.3-a", "CrystalSystem" -> "Cubic",
    "Group" -> typeIV[{195, 3}], "GroupCode" -> "typeIV[{195, 3}]",
    "Key" -> {195, 3}, "Letter" -> "a"|>
};

runWyckoffMode[case_Association, mode_String] := Module[
  {seed, operations, hamiltonian},
  seed = wyckoffDatabaseSeed[case["Key"], case["Letter"]];
  operations = Block[{Print = Function[Null]}, msgop[case["Group"]]];
  init[
    lattice -> wyckoffLattice[case["CrystalSystem"]],
    lattpar -> wyckoffLatticeParameters[case["CrystalSystem"]],
    wyckoffposition -> {seed},
    symminformation -> operations,
    basisFunctions -> {{"s"}},
    InitialBondShells -> 2,
    GenerateSymmetryGroup -> False,
    RepresentationMode -> mode
  ];
  hamiltonian = Sum[symham[i], {i, 2}];
  <|
    "DatabaseSeed" -> seed,
    "OrbitalTable" -> orbitalTable[],
    "Hamiltonian" -> hamiltonian,
    "Representation" -> showSymmetryRepresentations[],
    "BandManipulate" -> bandManipulate[
      wyckoffBandPath, 25, hamiltonian
    ],
    "BandPlot" -> bandplot[
      wyckoffBandPath,
      40,
      hamiltonian,
      wyckoffParameterRules[hamiltonian]
    ]
  |>
];
