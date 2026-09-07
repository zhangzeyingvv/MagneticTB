(* Generate only the physical figure used by the concise Chinese quick-start
   PDF.  The example loads the formal MagneticTB entry point and uses the same
   verified graphene parameters as the Mathematica help notebooks. *)

ClearAll["Global`*"];

scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
projectRoot = DirectoryName[DirectoryName[DirectoryName[scriptDirectory]]];
figureDirectory = FileNameJoin[{projectRoot, "Documentation", "figures"}];

Get[FileNameJoin[{projectRoot, "MagneticTB", "Kernel", "init.wl"}]];

grapheneOperations = MagneticTB`msgop[MagneticTB`gray[191]];

MagneticTB`init[
  MagneticTB`lattice -> {
    {MagneticTB`a, 0, 0},
    {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
    {0, 0, MagneticTB`c}
  },
  MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 3},
  MagneticTB`wyckoffposition -> {
    {{1/3, 2/3, 0}, {0, 0, 0}}
  },
  MagneticTB`symminformation -> grapheneOperations,
  MagneticTB`basisFunctions -> {{"pz"}},
  MagneticTB`InitialBondShells -> 3,
  MagneticTB`GenerateSymmetryGroup -> False
];

grapheneHamiltonian = Total[MagneticTB`symham /@ {1, 2, 3}];
grapheneParameters = {
  Global`e1 -> 0.05,
  Global`r1 -> 0.02,
  Global`t1 -> 0.5
};
graphenePath = {
  {{{0, 0, 0}, {0, 1/2, 0}}, {"G", "M"}},
  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
  {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}}
};

grapheneBandPlot = MagneticTB`bandplot[
  graphenePath,
  80,
  grapheneHamiltonian,
  grapheneParameters
];

grapheneKPoint = {1/3, 1/3, 0};
grapheneGammaHamiltonian = Chop@N[
  grapheneHamiltonian /. grapheneParameters /.
    {MagneticTB`kx -> 0, MagneticTB`ky -> 0, MagneticTB`kz -> 0}
];
grapheneKEnergies = Sort@N@Eigenvalues[
  grapheneHamiltonian /. grapheneParameters /.
    Thread[
      {MagneticTB`kx, MagneticTB`ky, MagneticTB`kz} ->
        2 Pi grapheneKPoint
    ]
];
grapheneKGap = Chop[grapheneKEnergies[[2]] - grapheneKEnergies[[1]]];

figurePath = FileNameJoin[{figureDirectory, "quickstart-graphene-band.pdf"}];
report = <|
  "HamiltonianDimensions" -> Dimensions[grapheneHamiltonian],
  "GammaHamiltonian" -> grapheneGammaHamiltonian,
  "KPointEnergies" -> grapheneKEnergies,
  "KPointDirectGap" -> grapheneKGap,
  "BandFigure" -> figurePath
|>;

(* Exporting a typeset Graphics to PDF requires the real FrontEnd on macOS.
   Keep the numerical calculation in the clean kernel, and use the FrontEnd
   only for this final vector-graphics serialization step. *)
Needs["JLink`"];
UsingFrontEnd[Export[figurePath, grapheneBandPlot]];

Print[InputForm[report]];
Quit[0];
