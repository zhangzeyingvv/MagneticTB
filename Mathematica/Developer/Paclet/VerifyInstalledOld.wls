ClearAll["Global`*"];

version = Get["PacletInfo.wl"]["Version"];
archive = File[ExpandFileName[
  FileNameJoin[{"build", "MagneticTB-" <> version <> ".paclet"}]
]];
installed = PacletInstall[archive];
loadOK = Check[Needs["MagneticTBOld`"] ; True, False];

initResult = If[
  loadOK,
  Check[
    MagneticTBOld`initold[
      MagneticTBOld`latticeold -> IdentityMatrix[3],
      MagneticTBOld`lattparold -> {},
      MagneticTBOld`wyckoffpositionold -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTBOld`symminformationold -> {
        {"1", IdentityMatrix[3], {0, 0, 0}, "F"}
      },
      MagneticTBOld`basisFunctionsold -> {{"s"}}
    ],
    $Failed
  ],
  $Failed
];

hamiltonian = If[
  initResult =!= $Failed,
  Check[MagneticTBOld`symhamold[2], $Failed],
  $Failed
];

report = <|
  "InstalledVersion" -> If[
    Head[installed] === PacletObject,
    installed["Version"],
    $Failed
  ],
  "LoadOld" -> loadOK,
  "NewContextLoaded" -> MemberQ[$Packages, "MagneticTB`"],
  "InitOldSucceeded" -> (initResult =!= $Failed),
  "HamiltonianDimensions" -> If[MatrixQ[hamiltonian], Dimensions[hamiltonian], $Failed]
|>;

Print[InputForm[report]];
Exit[If[
  TrueQ[
    loadOK &&
    initResult =!= $Failed &&
    MatrixQ[hamiltonian] &&
    !MemberQ[$Packages, "MagneticTB`"]
  ],
  0,
  1
]];
