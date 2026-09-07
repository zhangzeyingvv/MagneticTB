ClearAll["Global`*"];

version = Get["PacletInfo.wl"]["Version"];
archive = File[ExpandFileName[
  FileNameJoin[{"build", "MagneticTB-" <> version <> ".paclet"}]
]];
installed = PacletInstall[archive];
loadOK = Check[Needs["MagneticTB`"] ; True, False];

initResult = If[
  loadOK,
  Check[
    MagneticTB`init[
      MagneticTB`lattice -> IdentityMatrix[3],
      MagneticTB`lattpar -> {},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> {
        {"E", IdentityMatrix[3], {0, 0, 0}, "F"}
      },
      MagneticTB`basisFunctions -> {{"s"}},
      MagneticTB`InitialBondShells -> 2
    ],
    $Failed
  ],
  $Failed
];

hamiltonian = If[
  initResult =!= $Failed,
  Check[MagneticTB`symham[2], $Failed],
  $Failed
];

englishHelp = Documentation`ResolveLink[
  "paclet:MagneticTB/guide/MagneticTB",
  "English"
];
chineseHelp = Documentation`ResolveLink[
  "paclet:MagneticTB/guide/MagneticTB",
  "ChineseSimplified"
];
englishHelpExists = Check[FileExistsQ[englishHelp], False];
chineseHelpExists = Check[FileExistsQ[chineseHelp], False];
chineseNotebook = If[
  chineseHelpExists,
  Check[Get[chineseHelp], $Failed],
  $Failed
];
chineseContentPresent = If[
  Head[chineseNotebook] === Notebook,
  AnyTrue[
    Cases[chineseNotebook, text_String, Infinity],
    StringContainsQ[#, "\:78c1\:6027"] &
  ],
  False
];

report = <|
  "InstalledVersion" -> If[
    Head[installed] === PacletObject,
    installed["Version"],
    $Failed
  ],
  "LoadNew" -> loadOK,
  "OldContextLoaded" -> MemberQ[$Packages, "MagneticTBOld`"],
  "InitSucceeded" -> (initResult =!= $Failed),
  "HamiltonianDimensions" -> If[MatrixQ[hamiltonian], Dimensions[hamiltonian], $Failed],
  "EnglishHelpResolved" -> englishHelpExists,
  "ChineseHelpResolved" -> chineseHelpExists,
  "ChineseContentPresent" -> chineseContentPresent
|>;

Print[InputForm[report]];
Exit[If[
  TrueQ[
    loadOK &&
    initResult =!= $Failed &&
    MatrixQ[hamiltonian] &&
    !MemberQ[$Packages, "MagneticTBOld`"] &&
    englishHelpExists &&
    chineseHelpExists &&
    chineseContentPresent
  ],
  0,
  1
]];
