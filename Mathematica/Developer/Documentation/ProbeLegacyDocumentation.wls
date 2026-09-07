ClearAll["Global`*"];

projectRoot = Directory[];
outputPath = FileNameJoin[{
  projectRoot,
  "Developer",
  "Documentation",
  ".LegacyDocumentationFixture.wl"
}];
Get[FileNameJoin[{
  projectRoot, "MagneticTB", "Old", "MagneticTBOld", "Kernel", "init.wl"
}]];

fixture = Block[{Print = Function[Null]},
operations = MagneticTBOld`msgopold[MagneticTBOld`grayold[191]];
MagneticTBOld`initold[
  MagneticTBOld`latticeold -> {
    {MagneticTBOld`aold, 0, 0},
    {-MagneticTBOld`aold/2, Sqrt[3] MagneticTBOld`aold/2, 0},
    {0, 0, MagneticTBOld`cold}
  },
  MagneticTBOld`lattparold -> {
    MagneticTBOld`aold -> 1,
    MagneticTBOld`cold -> 3
  },
  MagneticTBOld`wyckoffpositionold -> {
    {{1/3, 2/3, 0}, {0, 0, 0}}
  },
  MagneticTBOld`symminformationold -> operations,
  MagneticTBOld`basisFunctionsold -> {{"pz"}},
  MagneticTBOld`debugQold -> False
];
h1 = MagneticTBOld`symhamold[1];
h2 = MagneticTBOld`symhamold[2];
h2Cartesian = MagneticTBOld`symhamold[
  2,
  MagneticTBOld`symmetrysetold -> All,
  "CartesianCoordinates" -> True
];
path = {
  {{{0, 0, 0}, {0, 1/2, 0}}, {"G", "M"}},
  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
  {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}}
};
<|
  "Onsite" -> h1,
  "Nearest" -> h2,
  "CartesianNearest" -> h2Cartesian,
  "Hamiltonian" -> h1 + h2,
  "BandUI" -> MagneticTBOld`bandManipulateold[path, 30, h1 + h2],
  "BandPlot" -> MagneticTBOld`bandplotold[
    path,
    60,
    h1 + h2,
    {Global`e1old -> 0.05, Global`t1old -> 0.5}
  ]
|>
];

putResult = Put[fixture, outputPath];
If[!FileExistsQ[outputPath], Exit[3]];
