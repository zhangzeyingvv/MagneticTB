ClearAll["Global`*"];

projectRoot = Directory[];
outputPath = FileNameJoin[{
  projectRoot,
  "Developer",
  "Documentation",
  ".InitFromRepDocumentationFixture.wl"
}];
Get[FileNameJoin[{
  projectRoot, "MagneticTB", "Kernel", "init.wl"
}]];

{fixture, fixtureMessages} = Block[{$MessageList = {}, Print = Function[Null]},
  value = Module[
    {
      session, orbitalTableData, orbitalTableGrid, onsite,
      representationMatrices, representationGrid, bondGrid,
      inducedOnsite, inducedRepresentationGrid,
      modeDirectOrbitalTable, modeInducedOrbitalTable,
      c4Operations, c4Matrices
    },
    c4Operations = {
      {"E", IdentityMatrix[3], {0, 0, 0}, "F"},
      {
        "C4z",
        {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
        {0, 0, 0},
        "F"
      },
      {"C2z", DiagonalMatrix[{-1, -1, 1}], {0, 0, 0}, "F"},
      {
        "C4z^-1",
        {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
        {0, 0, 0},
        "F"
      }
    };
    c4Matrices = {
      IdentityMatrix[2],
      {{0, -1}, {1, 0}},
      -IdentityMatrix[2],
      {{0, 1}, {-1, 0}}
    };

    MagneticTB`initfromrep[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> c4Operations,
      MagneticTB`repinformation -> {c4Matrices},
      MagneticTB`orbitalLabels -> {{"px", "py"}},
      MagneticTB`debugQ -> False,
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`RepresentationMode -> "DirectProduct"
    ];
    session = MagneticTB`CurrentModelSession[];
    orbitalTableData =
      MagneticTB`Private`buildHamiltonianBasisOrderData[session];
    orbitalTableGrid = MagneticTB`orbitalTable[];
    onsite = MagneticTB`symham[1];
    representationMatrices =
      session["FullRepresentation", "RepresentationMatrices"];
    representationGrid = MagneticTB`showSymmetryRepresentations[All];
    bondGrid = MagneticTB`showbonds[2];
    MagneticTB`initfromrep[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> c4Operations,
      MagneticTB`repinformation -> {<|
        "ReferenceSiteIndex" -> 1,
        "SiteSymmetryMatrices" -> c4Matrices
      |>},
      MagneticTB`orbitalLabels -> {{"px", "py"}},
      MagneticTB`debugQ -> False,
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`RepresentationMode -> "Induced"
    ];
    inducedOnsite = MagneticTB`symham[1];
    inducedRepresentationGrid =
      MagneticTB`showSymmetryRepresentations[All];

    (* Use one off-origin C4 orbit to make the difference between the two
       representation constructions visible in the reference page.  Direct
       product keeps px and py on both sites; induction starts from one px
       coordinate and transports it to the C4-related site. *)
    MagneticTB`initfromrep[
      MagneticTB`lattice -> DiagonalMatrix[{1, 1, 2}],
      MagneticTB`lattpar -> {},
      MagneticTB`wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> c4Operations,
      MagneticTB`repinformation -> {c4Matrices},
      MagneticTB`orbitalLabels -> {{"px", "py"}},
      MagneticTB`RepresentationMode -> "DirectProduct"
    ];
    modeDirectOrbitalTable = MagneticTB`orbitalTable[];

    MagneticTB`initfromrep[
      MagneticTB`lattice -> DiagonalMatrix[{1, 1, 2}],
      MagneticTB`lattpar -> {},
      MagneticTB`wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> c4Operations,
      MagneticTB`repinformation -> {<|
        "ReferenceSiteIndex" -> 1,
        "SiteSymmetryMatrices" -> {{{1}}, {{-1}}}
      |>},
      MagneticTB`orbitalLabels -> {{"px"}},
      MagneticTB`RepresentationMode -> "Induced"
    ];
    modeInducedOrbitalTable = MagneticTB`orbitalTable[];

    <|
      "OrbitalTableData" -> orbitalTableData,
      "OrbitalTableGrid" -> orbitalTableGrid,
      "RepresentationMatrices" -> representationMatrices,
      "RepresentationGrid" -> representationGrid,
      "BondGrid" -> bondGrid,
      "Onsite" -> onsite,
      "InducedOnsite" -> inducedOnsite,
      "InducedRepresentationGrid" -> inducedRepresentationGrid,
      "ModeDirectOrbitalTable" -> modeDirectOrbitalTable,
      "ModeInducedOrbitalTable" -> modeInducedOrbitalTable
    |>
  ];
  {value, $MessageList}
];

If[
  !AssociationQ[fixture] || fixtureMessages =!= {},
  Print[InputForm@<|
    "DocumentationFixture" -> "InitFromRep",
    "Value" -> fixture,
    "Messages" -> fixtureMessages
  |>];
  Exit[3]
];
Put[fixture, outputPath];
If[!FileExistsQ[outputPath], Exit[4]];
