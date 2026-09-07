(* Incremental documentation specifications for plotting functions and the
   finite-system data contract added on the Mathematica testing line.  This
   file is intentionally self-contained: the targeted documentation build can
   load it without evaluating unrelated model fixtures. *)

ClearAll[
  plotPropertyEvaluatedExample,
  plotPropertyFunctionSpec,
  plotPropertyFailureNote,
  plotPropertyCurrentPage,
  plotPropertyOptionKey,
  plotPropertyCallCode,
  plotPropertyGraphicsOptionExamples,
  testingApplyPlotPropertyTutorialUpdate
];

plotPropertyEvaluatedExample[text_String, code_String] := Module[
  {requestedPage, evaluateQ},
  requestedPage = Environment["MAGNETICTB_DOC_PAGE"];
  evaluateQ = !StringQ[requestedPage] || StringTrim[requestedPage] === "" ||
    FileBaseName[requestedPage] === plotPropertyCurrentPage;
  {
    text,
    code,
    If[evaluateQ,
      With[{
          label = "TestingPlotProperty-" <> IntegerString[Hash[code], 36],
          source = code
        },
        checkedDocumentationEvaluation[label, ToExpression[source, InputForm]]
      ],
      Missing["NotEvaluatedForIncrementalPage", plotPropertyCurrentPage]
    ]
  }
];

plotPropertyFunctionSpec[
    name_String,
    usage_List,
    notes_List,
    optionTable_List,
    examples_List,
    optionExamples_List,
    seeAlso_List,
    keywords_List
  ] := <|
  "Name" -> name,
  "Usage" -> usage,
  "Notes" -> notes,
  "OptionTable" -> optionTable,
  "Examples" -> examples,
  "OptionExamples" -> optionExamples,
  "SeeAlso" -> seeAlso,
  "Keywords" -> keywords,
  "Tutorials" -> {
    {"Real-space Hamiltonians, surfaces, and topology",
      "MagneticTB/tutorial/RealSpaceAndTopology"}
  },
  "Expected" -> (
    StringTrim[#, "\""] & /@ optionTable[[All, 1]]
  )
|>;

plotPropertyFailureNote[function_String, conditions_String] :=
  function <> " returns $Failed and issues a message when " <> conditions <>
    ". It does not repair the input, reduce the grid silently, or switch to an approximate fallback.";

plotPropertyOptionKey[rule_String] := StringTrim@First@StringSplit[rule, "->"];

plotPropertyCallCode[
    function_String,
    arguments_String,
    baseRules_List,
    rule_String : ""
  ] := Module[{rules, targetKey},
  rules = baseRules;
  If[rule =!= "",
    targetKey = plotPropertyOptionKey[rule];
    rules = Select[rules, plotPropertyOptionKey[#] =!= targetKey &];
    AppendTo[rules, rule]
  ];
  function <> "[\n  " <> arguments <>
    If[rules === {}, "", ",\n  " <> StringRiffle[rules, ",\n  "]] <>
    "\n]"
];

(* Plotting option examples must display the plot itself.  The Basic Example
   immediately above them constructs the physical model once; each option
   subsection then shows only the public plotting call that changes that
   option.  This follows the standard reference-page reading order and avoids
   burying one useful line under a repeated model definition.  For the
   explicit "Output" -> "Data" example, select the embedded Graphics rather
   than replacing the figure with an Association summary. *)
plotPropertyGraphicsOptionExamples[
    records_List,
    setup_String,
    function_String,
    arguments_String,
    baseRules_List,
    setupFirst_: False
  ] := Module[{graphicsRules},
  graphicsRules = Select[
    baseRules,
    plotPropertyOptionKey[#] =!= "\"Output\"" &
  ];
  MapIndexed[
    Function[{record, index},
      With[{
          call = plotPropertyCallCode[
            function,
            arguments,
            graphicsRules,
            record[[3]]
          ]
        },
        plotPropertyEvaluatedExample[
          record[[1]] <> " - " <> record[[2]],
          If[TrueQ[setupFirst] && First[index] === 1, setup <> "\n", ""] <>
            If[
              plotPropertyOptionKey[record[[3]]] === "\"Output\"",
              "(" <> call <> ")[\"Graphics\"]",
              call
            ]
        ]
      ]
    ],
    records
  ]
];

(* -------------------------------------------------------------------------
   Shared minimal inputs.

   These strings are copied into every example that needs them.  Keeping the
   setups small makes the option examples fast enough for fresh-kernel replay.
   ------------------------------------------------------------------------- *)

(* Every numerical hopping table below comes from the public symmetry workflow.
   Even the smallest option example starts from a nontrivial magnetic group;
   no hand-written Hamiltonian or documentation-only Association is used. *)
plotPropertyCubicDataCode = StringRiffle[{
  "cubicOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {{a,0,0},{0,a,0},{0,0,a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {{{0,0,0},{0,0,0}}},",
  "  symminformation -> cubicOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "cubicHamiltonian = symham[1] + symham[2];",
  "cubicData = hoppingData[{1,2},{e1 -> 0.,t1 -> 1.}];"
}, "\n"];

plotPropertyChainDataCode = StringRiffle[{
  "chainOperations = msgop[gray[75]];",
  "init[",
  "  lattice -> {{a,0,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 3,c -> 1},",
  "  wyckoffposition -> {{{0,0,0},{0,0,0}}},",
  "  symminformation -> chainOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "chainHamiltonian = symham[1] + symham[2];",
  "plotPropertyChain = hoppingData[",
  "  {1,2},",
  "  {e1 -> 0.,t1 -> 1.}",
  "];"
}, "\n"];

plotPropertySquareDataCode = plotPropertyCubicDataCode <> "\n" <>
  StringRiffle[{
    "plotPropertySquare = cubicData;",
    "plotPropertyLCells = {{0,0,0},{1,0,0},{0,1,0}};"
  }, "\n"];

plotPropertyLWavefunctionCode = StringRiffle[{
  plotPropertySquareDataCode,
  "plotPropertyLSystem = buildRealSpaceHamiltonian[",
  "  plotPropertySquare, plotPropertyLCells, \"Output\" -> \"Data\"",
  "];",
  "plotPropertyLState = {1,I,-1};"
}, "\n"];

(* A C4-symmetric breathing-square higher-order topological insulator.  Its
   finite Hamiltonian and corner state are both derived from gray group 75. *)
plotPropertyHOTIWavefunctionCode = StringRiffle[{
  "hotiOperations = msgop[gray[75]];",
  "init[",
  "  lattice -> {{1,0,0},{0,1,0},{0,0,10}},",
  "  lattpar -> {},",
  "  wyckoffposition -> {{{1/5,1/5,0},{0,0,0}}},",
  "  symminformation -> hotiOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 4,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "symham[2];",
  "symham[4];",
  "hotiData = hoppingData[{2,4},{t1 -> 1.,s1 -> 0.2}];",
  "hotiFinite = buildRealSpaceHamiltonian[",
  "  hotiData,{6,6,1},\"Output\" -> \"Data\"",
  "];",
  "{hotiEnergies,hotiStates} = Eigensystem[",
  "  N@Normal@hotiFinite[\"Hamiltonian\"]",
  "];",
  "hotiCornerSubspace = hotiStates[[Ordering[Abs[hotiEnergies],4]]];",
  "hotiPositions = Lookup[",
  "  hotiFinite[\"BasisRecords\"],\"CartesianPosition\"",
  "];",
  "hotiCornerSeed = UnitVector[",
  "  hotiFinite[\"Dimension\"],",
  "  First@Ordering[Total /@ hotiPositions,1]",
  "];",
  "hotiCornerState = Normalize@Total[",
  "  (Conjugate[#].hotiCornerSeed) # & /@ hotiCornerSubspace",
  "];"
}, "\n"];

plotPropertyTwoOrbitalDataCode = StringRiffle[{
  "twoOrbitalOperations = msgop[gray[75]];",
  "init[",
  "  lattice -> {{1,0,0},{0,1,0},{0,0,5}},",
  "  lattpar -> {},",
  "  wyckoffposition -> {{{0,0,0},{0,0,0}}},",
  "  symminformation -> twoOrbitalOperations,",
  "  basisFunctions -> {{\"px\",\"py\"}},",
  "  InitialBondShells -> 1,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "symham[1];",
  "plotPropertyTwoOrbital = hoppingData[",
  "  1,{e1 -> 0.,e2 -> 0.}",
  "];",
  "plotPropertyTwoOrbitalSystem = buildRealSpaceHamiltonian[",
  "  plotPropertyTwoOrbital, {{0,0,0}}, \"Output\" -> \"Data\"",
  "];",
  "plotPropertyState = {1,I};"
}, "\n"];

(* A symmetry-generated P4 Chern insulator used by the Wilson, Berry, and
   surface pages.  The numeric Hamiltonian is only a wrapper around symham. *)
plotPropertyChernSetupCode = StringRiffle[{
  "chernOperations = msgop[bnsdict[{75,1}]];",
  "init[",
  "  lattice -> {{a,0,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 1,c -> 5},",
  "  wyckoffposition -> {{{0,0,0},{0,0,1}}},",
  "  symminformation -> chernOperations,",
  "  basisFunctions -> {{\"s\",x + I y}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "chernMatrix = symham[1] + symham[2];",
  "chernRules = {",
  "  e1 -> -1.,e2 -> 1.,t1 -> 0.5,",
  "  t2 -> 0.5,t3 -> -0.5,t4 -> 0.",
  "};",
  "chernData = hoppingData[{1,2},chernRules];",
  "chernHamiltonian[{qx_?NumericQ,qy_?NumericQ}] := N[",
  "  chernMatrix /. chernRules /.",
  "    {kx -> qx,ky -> qy,kz -> 0}",
  "];",
  "chernCenters = chernData[\"WannierCenters\"][[All,1 ;; 2]];",
  "chernWilsonPath = {",
  "  {{{0,-Pi},{0,0}},{\"-Pi\",\"0\"}},",
  "  {{{0,0},{0,Pi}},{\"0\",\"Pi\"}}",
  "};"
}, "\n"];

plotPropertyChernWilsonSetupCode = plotPropertyChernSetupCode;

(* The four-band Z2 model is generated from gray group 75 and its current
   symham parameter coordinates. *)
plotPropertyZ2SetupCode = StringRiffle[{
  "z2Operations = msgop[gray[75]];",
  "init[",
  "  lattice -> {{a,0,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 1,c -> 1},",
  "  wyckoffposition -> {",
  "    {{0,0,1/4},{0,0,0}},",
  "    {{0,0,0},{0,0,0}}",
  "  },",
  "  symminformation -> z2Operations,",
  "  basisFunctions -> {{\"px\",\"py\"},{\"px\",\"py\"}},",
  "  InitialBondShells -> 7,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "z2Shells = Association@Table[i -> symham[i],{i,{2,3,4,5,7}}];",
  "z2Rules = {",
  "  t1->0,t2->2.5,r1->0,r2->2,",
  "  s1->0,s2->0,s3->-1,s4->0,s5->0,",
  "  s6->0,s7->0,s8->1,s9->0,s10->0,",
  "  p5n1->0.5,p5n2->0,p5n3->0,p5n4->0.5,",
  "  p7n1->0,p7n2->0,p7n3->0,p7n4->0,",
  "  p7n5->0.25,p7n6->-0.25,p7n7->0.25,",
  "  p7n8->0,p7n9->0,p7n10->0,p7n11->0,",
  "  p7n12->-0.25,p7n13->-0.25,p7n14->-0.25",
  "};",
  "z2Hamiltonian = Total[Values[z2Shells]] /. z2Rules;",
  "z2Data = hoppingData[Keys[z2Shells],z2Rules];",
  "z2Slice[{qx_?NumericQ,qy_?NumericQ}] := N[",
  "  z2Hamiltonian /. {kx -> qx,ky -> qy,kz -> 0}",
  "];",
  "z2Centers = z2Data[\"WannierCenters\"][[All,1 ;; 2]];",
  "z2WilsonPath = {",
  "  {{{0,0},{0,Pi}},{\"0\",\"Pi\"}}",
  "};",
  "z2SurfacePath = {",
  "  {{{-1/2,0},{0,0}},{\"-Pi\",\"0\"}},",
  "  {{{0,0},{1/2,0}},{\"0\",\"Pi\"}}",
  "};"
}, "\n"];

plotPropertyZ2WilsonSetupCode = plotPropertyZ2SetupCode;

plotPropertySurfacePathCode = StringRiffle[{
  "topologicalSurfacePath = {",
  "  {{{-1/2,0},{0,0}},{\"-Pi\",\"0\"}},",
  "  {{{0,0},{1/2,0}},{\"0\",\"Pi\"}}",
  "};"
}, "\n"];

plotPropertyChernSurfaceSetupCode = StringRiffle[{
  plotPropertyChernSetupCode,
  "chernSurfaceData = transformHoppings[",
  "  chernData,{{1,0,0},{0,0,1},{0,1,0}}",
  "];",
  plotPropertySurfacePathCode
}, "\n"];

plotPropertyZ2SurfaceSetupCode = plotPropertyZ2SetupCode;

(* BNS 143.3 supplies the magnetic Weyl Hamiltonian used by the 3D page. *)
plotPropertyWeylSetupCode = StringRiffle[{
  "weylOperations = msgop[bnsdict[{143,3}]];",
  "init[",
  "  lattice -> {{Sqrt[3] a/2,-a/2,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 1,c -> 2},",
  "  wyckoffposition -> {{{0,0,0},{0,0,1}}},",
  "  symminformation -> weylOperations,",
  "  basisFunctions -> {{{x + I y,0},{0,x - I y}}},",
  "  InitialBondShells -> 4,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "weylBlock = Sum[symham[i],{i,{2,3,4}}][[{1,4},{1,4}]] /. {",
  "  r11 -> -r2,r15 -> -r10,r3 -> -r6,r7 -> -r14,",
  "  s4 -> -s4,t4 -> t1,t8 -> t9,t11 -> -t7",
  "};",
  "weylParameters = {",
  "  r10 -> -0.09,r14 -> 0.25,r2 -> 0,r6 -> 0.5,",
  "  s4 -> -0.38,s8 -> 0,t1 -> 0,t15 -> 0.735,",
  "  t7 -> 0.05,t9 -> 0.015",
  "};",
  "weylHamiltonian[q_List] := N[",
  "  weylBlock /. weylParameters /. Thread[{kx,ky,kz} -> q]",
  "];",
  "weylCenters = {{0,0,0},{0,0,1/2}};"
}, "\n"];

(* -------------------------------------------------------------------------
   Reference-page specifications.

   Each spec owns its usage forms, notes, option table, primary examples, and
   one executable example per public option.
   ------------------------------------------------------------------------- *)

plotPropertyCurrentPage = "buildRealSpaceHamiltonian";
plotPropertyBuildRealSpaceSpec = plotPropertyFunctionSpec[
  "buildRealSpaceHamiltonian",
  {
    {"buildRealSpaceHamiltonian[data,{n1,n2,n3}]",
      "assembles a rectangular sparse finite Hamiltonian from numerical hopping blocks."},
    {"buildRealSpaceHamiltonian[data,cells]",
      "assembles an arbitrary open-boundary shape from an explicit ordered list of integer cell coordinates."}
  },
  {
    "For rectangular input, cell indices run from zero to ni-1 along each direct-lattice axis. BoundaryConditions may independently make those axes Open or Periodic.",
    StringJoin[
      "An explicit cell list must be nonempty, duplicate free, and contain ",
      "integer three-vectors. Its input order is the matrix cell order, and ",
      "only fully open boundaries are supported."
    ],
    StringJoin[
      "The matrix convention is row=(source cell,orbital), ",
      "column=(target cell,orbital). Output -> \"Data\" returns SchemaVersion 2, ",
      "Shape, Cells, CellBounds, and matrix-aligned BasisRecords."
    ],
    StringJoin[
      "When the hopping data contain a numerical 3 by 3 Lattice and one ",
      "WannierCenters vector per orbital, every BasisRecords entry includes ",
      "the fractional and Cartesian position used by plotRealSpaceWavefunction."
    ],
    plotPropertyFailureNote[
      "buildRealSpaceHamiltonian",
      StringJoin[
        "the hopping data are incomplete, the geometry is malformed or duplicated, ",
        "an explicit shape requests a periodic boundary, Hermiticity exceeds ",
        "the tolerance, or Output is unsupported"
      ]
    ]
  },
  {
    {"\"BoundaryConditions\"", "{\"Open\",\"Open\",\"Open\"}",
      "Open or Periodic independently along the three axes of a rectangular geometry; explicit cell lists require all Open."},
    {"\"HermiticityTolerance\"", "10^-9",
      "Nonnegative tolerance applied to the hopping table and assembled finite matrix."},
    {"\"Output\"", "\"Matrix\"",
      "Matrix returns SparseArray; Data returns the matrix, geometry, basis ordering, and available position metadata."}
  },
  {
    plotPropertyEvaluatedExample[
      "Build a three-site L shape and display its finite Hamiltonian:",
      plotPropertySquareDataCode <>
        "\nplotPropertyData = buildRealSpaceHamiltonian[plotPropertySquare, plotPropertyLCells, \"Output\" -> \"Data\"];\n" <>
        "MatrixPlot[Normal[plotPropertyData[\"Hamiltonian\"]],FrameTicks -> None]"
    ],
    plotPropertyEvaluatedExample[
      "Plot the matrix-aligned Cartesian positions stored in BasisRecords:",
      plotPropertySquareDataCode <>
        "\nplotPropertyData = buildRealSpaceHamiltonian[plotPropertySquare, plotPropertyLCells, \"Output\" -> \"Data\"];\n" <>
        "ListPlot[Lookup[plotPropertyData[\"BasisRecords\"],\"CartesianPosition\"][[All,1 ;; 2]],PlotStyle -> Directive[Blue,PointSize[0.035]],AspectRatio -> 1,AxesLabel -> {\"x\",\"y\"}]"
    ]
  },
  {
    plotPropertyEvaluatedExample[
      StringJoin[
        "\"BoundaryConditions\" -> {\"Open\",\"Open\",\"Periodic\"} wraps the chain ",
        "periodically. MatrixPlot displays the actual finite Hamiltonian, including the ",
        "corner couplings created by the periodic boundary:"
      ],
      plotPropertyChainDataCode <> StringRiffle[{
        "",
        "periodicChain = buildRealSpaceHamiltonian[",
        "  plotPropertyChain, {1,1,3},",
        "  \"BoundaryConditions\" -> {\"Open\",\"Open\",\"Periodic\"}",
        "];",
        "MatrixPlot[Normal[periodicChain], FrameTicks -> None]"
      }, "\n"]
    ],
    plotPropertyEvaluatedExample[
      StringJoin[
        "\"HermiticityTolerance\" -> 10^-12 applies a stricter residual threshold to the ",
        "same exact hopping data. MatrixPlot shows the accepted finite Hamiltonian:"
      ],
      plotPropertyChainDataCode <> StringRiffle[{
        "",
        "strictChain = buildRealSpaceHamiltonian[",
        "  plotPropertyChain, {1,1,3},",
        "  \"HermiticityTolerance\" -> 10^-12",
        "];",
        "MatrixPlot[Normal[strictChain], FrameTicks -> None]"
      }, "\n"]
    ],
    plotPropertyEvaluatedExample[
      StringJoin[
        "\"Output\" -> \"Data\" returns the Hamiltonian together with its geometry and basis ",
        "metadata. Extract the public Hamiltonian field and plot the matrix itself:"
      ],
      plotPropertyChainDataCode <> StringRiffle[{
        "",
        "chainDataResult = buildRealSpaceHamiltonian[",
        "  plotPropertyChain, {1,1,3},",
        "  \"Output\" -> \"Data\"",
        "];",
        "MatrixPlot[",
        "  Normal[chainDataResult[\"Hamiltonian\"]],",
        "  FrameTicks -> None",
        "]"
      }, "\n"]
    ]
  },
  {
    {"plotRealSpaceWavefunction", "MagneticTB/ref/plotRealSpaceWavefunction"},
    {"buildSlabHamiltonian", "MagneticTB/ref/buildSlabHamiltonian"},
    {"buildBlochHamiltonian", "MagneticTB/ref/buildBlochHamiltonian"}
  },
  {"finite Hamiltonian", "arbitrary shape", "basis records", "SparseArray"}
];

plotPropertyWilsonOptionTable = {
  {"\"ParameterSubdivisions\"", "60", "Positive number of intervals per labeled parameter-path segment."},
  {"\"LoopSubdivisions\"", "50", "Positive number of forward links in every translated Wilson loop."},
  {"\"HermitianTolerance\"", "10^-10", "Nonnegative Hamiltonian Hermiticity and Wilson-unitarity tolerance."},
  {"\"GapTolerance\"", "10^-9", "Nonnegative minimum direct gap for a partially occupied subspace."},
  {"\"CovarianceTolerance\"", "10^-8", "Nonnegative endpoint sewing-covariance tolerance."},
  {"\"OverlapTolerance\"", "10^-10", "Nonnegative lower bound for occupied-link singular values."},
  {"\"PhaseConvention\"", "\"PhaseOverPi\"", "PhaseOverPi gives phases on [-1,1]; WannierCenters gives centers modulo one."},
  {"\"Output\"", "\"Graphics\"", "Graphics returns the branch plot; Data also returns the sampled path, values, labels, and Graphics."},
  {"yTicks", "Automatic", "Automatic, None, numeric positions, or standard vertical tick specifications."},
  {"Joined", "False", "Whether sampled points in each Wilson branch are connected."},
  {"PlotMarkers", "Automatic", "Marker specification passed to ListPlot."},
  {"PlotStyle", "Black", "Style or list of styles passed to ListPlot."},
  {"PlotRange", "Automatic", "Displayed x/y range; Automatic uses the full path and the selected phase convention."},
  {"GridLines", "Automatic", "Grid-line specification; Automatic marks path boundaries and reference phases."},
  {"FrameLabel", "Automatic", "Two frame labels; Automatic names the parameter path and phase convention."},
  {"FontSize", "18", "Positive size for ticks, path labels, and frame labels."},
  {"FontFamily", "\"Times\"", "Font family for ticks, path labels, and frame labels."},
  {"ImageSize", "Large", "Final Graphics image size."}
};

plotPropertyWilsonBaseRules = {
  "\"ParameterSubdivisions\" -> 1",
  "\"LoopSubdivisions\" -> 4",
  "\"Output\" -> \"Data\""
};

plotPropertyCurrentPage = "plotWilsonLoop";
plotPropertyWilsonSpec = plotPropertyFunctionSpec[
  "plotWilsonLoop",
  {
    {"plotWilsonLoop[H,centers,occupied,start,end,parameterPath]",
      "plots Wilson-loop branches while translating one reciprocal-lattice-closed loop along a parameter path."}
  },
  {
    StringJoin[
      "start and end have the same momentum-coordinate dimension, and end-start ",
      "must be 2 Pi times an integer reciprocal vector. centers contains one ",
      "fractional orbital center per Hamiltonian basis state."
    ],
    StringJoin[
      "parameterPath may be an explicit ordered list of distinct momentum offsets ",
      "or connected labeled segments {{{p1,p2},{label1,label2}},...}. Labeled ",
      "segments use ParameterSubdivisions intervals per segment."
    ],
    StringJoin[
      "Every sample calls wilsonLoop with forward occupied-subspace links, ",
      "SVD unitarization, and orbital-center endpoint sewing. The plot introduces ",
      "no gauge, gap, or numerical fallback."
    ],
    StringJoin[
      "The symmetry-generated P4 Chern example has one occupied band whose Wilson phase winds once over a full ",
      "transverse Brillouin zone. The gray-P4 Z2 example has two occupied bands; their Wannier-center ",
      "branches exchange partners between the two time-reversal-invariant endpoints of half the zone. ",
      "plotWilsonLoop displays this flow but does not return a Z2 integer."
    ],
    "The default return is Graphics. Output -> \"Data\" returns ParameterPath, PathCoordinate, boundary labels, PhaseConvention, Values, BranchCount, and the same Graphics.",
    plotPropertyFailureNote[
      "plotWilsonLoop",
      "the loop or parameter path is malformed, a translated Wilson loop fails its Hermiticity/gap/covariance/overlap checks, or an option is invalid"
    ]
  },
  plotPropertyWilsonOptionTable,
  {
    plotPropertyEvaluatedExample[
      "The occupied band of the symmetry-generated P4 Chern insulator has one Wilson branch that winds once across the transverse Brillouin zone:",
      plotPropertyChernWilsonSetupCode <> "\n" <>
        "plotWilsonLoop[\n" <>
        "  chernHamiltonian, chernCenters, 1, {0,0}, {2 Pi,0},\n" <>
        "  chernWilsonPath,\n" <>
        "  \"ParameterSubdivisions\" -> 30,\n" <>
        "  \"LoopSubdivisions\" -> 60,\n" <>
        "  Joined -> False, PlotMarkers -> Automatic, PlotStyle -> Blue\n" <>
        "]"
    ],
    plotPropertyEvaluatedExample[
      "For the symmetry-generated gray-P4 Z2 insulator, the two occupied Wannier-center branches exchange partners between 0 and Pi:",
      plotPropertyZ2WilsonSetupCode <> "\n" <>
        "plotWilsonLoop[\n" <>
        "  z2Slice, z2Centers, 2, {0,0}, {2 Pi,0},\n" <>
        "  z2WilsonPath,\n" <>
        "  \"ParameterSubdivisions\" -> 40,\n" <>
        "  \"LoopSubdivisions\" -> 60,\n" <>
        "  \"PhaseConvention\" -> \"WannierCenters\",\n" <>
        "  Joined -> False, PlotMarkers -> Automatic,\n" <>
        "  PlotStyle -> {Blue,Red}\n" <>
        "]"
    ]
  },
  plotPropertyGraphicsOptionExamples[
    {
      {"\"ParameterSubdivisions\"", "uses two intervals for the labeled parameter segment.", "\"ParameterSubdivisions\" -> 2"},
      {"\"LoopSubdivisions\"", "uses six forward links in each translated loop.", "\"LoopSubdivisions\" -> 6"},
      {"\"HermitianTolerance\"", "uses a stricter Hermiticity and unitarity check.", "\"HermitianTolerance\" -> 10^-12"},
      {"\"GapTolerance\"", "sets the occupied direct-gap threshold.", "\"GapTolerance\" -> 10^-12"},
      {"\"CovarianceTolerance\"", "sets the endpoint covariance threshold.", "\"CovarianceTolerance\" -> 10^-10"},
      {"\"OverlapTolerance\"", "sets the minimum link singular value.", "\"OverlapTolerance\" -> 10^-12"},
      {"\"PhaseConvention\"", "plots Wannier centers modulo one.", "\"PhaseConvention\" -> \"WannierCenters\""},
      {"\"Output\"", "returns sampled values and Graphics together.", "\"Output\" -> \"Data\""},
      {"yTicks", "sets explicit vertical tick positions.", "yTicks -> {0,1/2,1}"},
      {"Joined", "connects successive samples in each branch.", "Joined -> True"},
      {"PlotMarkers", "suppresses point markers.", "PlotMarkers -> None"},
      {"PlotStyle", "draws the branch in red.", "PlotStyle -> Red"},
      {"PlotRange", "sets an explicit displayed range.", "PlotRange -> {{0,1},{-1,1}}"},
      {"GridLines", "suppresses reference grid lines.", "GridLines -> None"},
      {"FrameLabel", "sets both frame labels explicitly.", "FrameLabel -> {\"offset\",\"phase\"}"},
      {"FontSize", "sets a smaller label and tick size.", "FontSize -> 14"},
      {"FontFamily", "uses Helvetica for plot text.", "FontFamily -> \"Helvetica\""},
      {"ImageSize", "sets an explicit image width.", "ImageSize -> 320"}
    },
    plotPropertyChernWilsonSetupCode,
    "plotWilsonLoop",
    "chernHamiltonian, chernCenters, 1, {0,0}, {2 Pi,0}, chernWilsonPath",
    plotPropertyWilsonBaseRules
  ],
  {
    {"wilsonLoop", "MagneticTB/ref/wilsonLoop"},
    {"plotSurfaceSpectrum", "MagneticTB/ref/plotSurfaceSpectrum"},
    {"plotBerryCurvature2D", "MagneticTB/ref/plotBerryCurvature2D"}
  },
  {
    "Wilson loop plot", "Wannier centers", "topological branch flow",
    "Chern insulator", "Z2 insulator", "quantum spin Hall effect"
  }
];

AssociateTo[
  plotPropertyWilsonSpec,
  "Expected" -> DeleteDuplicates@Join[
    plotPropertyWilsonSpec["Expected"],
    {"symmetry-generated P4 Chern insulator", "symmetry-generated gray-P4 Z2 insulator"}
  ]
];

plotPropertySurfaceOptionTable = {
  {"\"MomentumSubdivisions\"", "60", "Positive number of intervals per labeled surface-momentum segment."},
  {"\"EnergyPoints\"", "201", "Number of uniformly spaced energies, at least two."},
  {"\"CellMatrix\"", "Automatic", "Automatic keeps the hopping cell; a nonsingular integer 3 by 3 matrix defines the surface cell first."},
  {"\"Broadening\"", "10^-3", "Positive imaginary retarded-energy broadening."},
  {"\"Tolerance\"", "10^-10", "Positive relative coupling tolerance for principal-layer decimation."},
  {"\"MaxIterations\"", "200", "Positive decimation iteration limit."},
  {"\"Surface\"", "\"Positive\"", "Positive or Negative termination along cell axis 3."},
  {"\"HermiticityTolerance\"", "10^-9", "Nonnegative tolerance for hopping and principal-layer Hermiticity."},
  {"\"Output\"", "\"Graphics\"", "Graphics returns the density plot; Data also returns the sampled path, energies, spectral weights, and Graphics."},
  {"ColorFunction", "\"SolarColors\"", "Color scheme or color function passed to ListDensityPlot."},
  {"ColorFunctionScaling", "True", "Whether ListDensityPlot rescales values before applying ColorFunction."},
  {"PlotLegends", "Automatic", "Legend specification passed to ListDensityPlot."},
  {"PlotRange", "All", "Displayed plot range passed to ListDensityPlot."},
  {"InterpolationOrder", "1", "Zero or one for the sampled spectral grid."},
  {"Mesh", "None", "Mesh specification passed to ListDensityPlot."},
  {"FrameLabel", "Automatic", "Two frame labels; Automatic uses surface k path and energy."},
  {"FontSize", "18", "Positive size for frame ticks and labels."},
  {"FontFamily", "\"Times\"", "Font family for frame ticks and labels."},
  {"AspectRatio", "1/GoldenRatio", "Height-to-width ratio of the density plot."},
  {"ImageSize", "Large", "Final Graphics image size."}
};

plotPropertySurfaceBaseRules = {
  "\"MomentumSubdivisions\" -> 1",
  "\"EnergyPoints\" -> 2",
  "\"Broadening\" -> 0.05",
  "PlotLegends -> None",
  "\"Output\" -> \"Data\""
};

plotPropertyCurrentPage = "plotSurfaceSpectrum";
plotPropertySurfaceSpec = plotPropertyFunctionSpec[
  "plotSurfaceSpectrum",
  {
    {"plotSurfaceSpectrum[data,path,{emin,emax}]",
      "plots the semi-infinite surface spectral weight along a two-dimensional surface-momentum path."}
  },
  {
    "data is the numerical hopping Association accepted by surfaceGreenFunction. Cell axes 1 and 2 span the surface, and axis 3 is the semi-infinite stacking direction.",
    "path may be explicit surface-momentum points or connected labeled two-dimensional segments. The plotted value is -Im Tr[Gsurface]/Pi on a uniform energy grid.",
    "Every grid point uses the same principal-layer decimation and convergence checks as surfaceSpectralFunction; nonconvergence is reported instead of hidden.",
    StringJoin[
      "For the P4 Chern example, transformHoppings moves the original second in-plane direction to cell ",
      "axis 3; axis 1 is then the conserved edge momentum and axis 3 is semi-infinite. The gray-P4 ",
      "four-band example uses its stored cell with axis 3 as the stacking direction. The resulting plots ",
      "show one chiral gap-crossing branch and a counterpropagating surface-state pair, respectively."
    ],
    "Output -> \"Data\" returns SurfaceMomentumPath, PathCoordinate, boundary labels, Energies, SpectralWeight, Broadening, Surface, and Graphics.",
    plotPropertyFailureNote[
      "plotSurfaceSpectrum",
      "the hopping data, path, ordered energy range, surface options, decimation controls, or display options are invalid, or a grid-point calculation does not converge"
    ]
  },
  plotPropertySurfaceOptionTable,
  {
    plotPropertyEvaluatedExample[
      "The symmetry-generated P4 Chern insulator has a single chiral edge branch that crosses the bulk gap:",
      plotPropertyChernSurfaceSetupCode <> "\n" <>
        "plotSurfaceSpectrum[\n" <>
        "  chernSurfaceData, topologicalSurfacePath, {-1.5,1.5},\n" <>
        "  \"MomentumSubdivisions\" -> 35,\n" <>
        "  \"EnergyPoints\" -> 121,\n" <>
        "  \"Broadening\" -> 0.025, PlotLegends -> None\n" <>
        "]"
    ],
    plotPropertyEvaluatedExample[
      "The symmetry-generated gray-P4 Z2 model has a pair of helical surface branches crossing at zero momentum:",
      plotPropertyZ2SurfaceSetupCode <> "\n" <>
        "plotSurfaceSpectrum[\n" <>
        "  z2Data, z2SurfacePath, {-3,3},\n" <>
        "  \"MomentumSubdivisions\" -> 35,\n" <>
        "  \"EnergyPoints\" -> 121,\n" <>
        "  \"Broadening\" -> 0.025, PlotLegends -> None\n" <>
        "]"
    ]
  },
  plotPropertyGraphicsOptionExamples[
    {
      {"\"MomentumSubdivisions\"", "uses two intervals on the labeled path.", "\"MomentumSubdivisions\" -> 2"},
      {"\"EnergyPoints\"", "samples three uniformly spaced energies.", "\"EnergyPoints\" -> 3"},
      {"\"CellMatrix\"", "keeps the axes through an explicit integer identity cell.", "\"CellMatrix\" -> IdentityMatrix[3]"},
      {"\"Broadening\"", "uses a larger positive retarded broadening.", "\"Broadening\" -> 0.1"},
      {"\"Tolerance\"", "sets the decimation convergence threshold.", "\"Tolerance\" -> 10^-8"},
      {"\"MaxIterations\"", "sets an explicit iteration bound.", "\"MaxIterations\" -> 100"},
      {"\"Surface\"", "plots the negative termination.", "\"Surface\" -> \"Negative\""},
      {"\"HermiticityTolerance\"", "uses a stricter hopping residual threshold.", "\"HermiticityTolerance\" -> 10^-12"},
      {"\"Output\"", "returns the numerical grid and Graphics together.", "\"Output\" -> \"Data\""},
      {"ColorFunction", "uses TemperatureMap for the density colors.", "ColorFunction -> \"TemperatureMap\""},
      {"ColorFunctionScaling", "passes unscaled spectral values to the color function.", "ColorFunctionScaling -> False"},
      {"PlotLegends", "suppresses the color legend.", "PlotLegends -> None"},
      {"PlotRange", "uses the complete sampled range explicitly.", "PlotRange -> All"},
      {"InterpolationOrder", "uses piecewise-constant cells.", "InterpolationOrder -> 0"},
      {"Mesh", "draws the sampling mesh.", "Mesh -> All"},
      {"FrameLabel", "sets explicit path and energy labels.", "FrameLabel -> {\"surface path\",\"energy\"}"},
      {"FontSize", "sets a smaller frame text size.", "FontSize -> 14"},
      {"FontFamily", "uses Helvetica for frame text.", "FontFamily -> \"Helvetica\""},
      {"AspectRatio", "uses a square density plot.", "AspectRatio -> 1"},
      {"ImageSize", "sets an explicit image width.", "ImageSize -> 320"}
    },
    plotPropertyChernSurfaceSetupCode,
    "plotSurfaceSpectrum",
    "chernSurfaceData, topologicalSurfacePath, {-1.5,1.5}",
    plotPropertySurfaceBaseRules
  ],
  {
    {"surfaceGreenFunction", "MagneticTB/ref/surfaceGreenFunction"},
    {"surfaceSpectralFunction", "MagneticTB/ref/surfaceSpectralFunction"},
    {"plotWilsonLoop", "MagneticTB/ref/plotWilsonLoop"}
  },
  {
    "surface spectrum", "surface Green function", "spectral weight",
    "Chern insulator", "Z2 insulator", "helical edge state"
  }
];

AssociateTo[
  plotPropertySurfaceSpec,
  "Expected" -> DeleteDuplicates@Join[
    plotPropertySurfaceSpec["Expected"],
    {"symmetry-generated P4 Chern insulator", "symmetry-generated gray-P4 Z2 model"}
  ]
];

plotPropertyBerry2DOptionTable = {
  {"\"Directions\"", "{1,2}", "Two distinct ordered momentum-coordinate indices; reversing them reverses curvature orientation."},
  {"\"FixedCoordinates\"", "Automatic", "Automatic uses zeros; an explicit vector fixes all unplotted momentum coordinates."},
  {"\"GridSize\"", "31", "One integer or two integers, each at least two, for the plotted grid."},
  {"\"StepSize\"", "Automatic", "Automatic derives plaquette steps from grid spacing; a positive scalar or full coordinate vector is also accepted."},
  {"\"HermitianTolerance\"", "10^-10", "Nonnegative Hamiltonian Hermiticity tolerance."},
  {"\"GapTolerance\"", "10^-9", "Nonnegative occupied direct-gap threshold."},
  {"\"CovarianceTolerance\"", "10^-8", "Nonnegative plaquette closure-covariance tolerance."},
  {"\"OverlapTolerance\"", "10^-10", "Nonnegative occupied-link singular-value threshold."},
  {"\"SymmetricColorRange\"", "True", "Whether Automatic PlotRange is symmetric about zero."},
  {"\"Output\"", "\"Graphics\"", "Graphics returns the density plot; Data also returns the grid, curvature, and Graphics."},
  {"ColorFunction", "Automatic", "Automatic uses TemperatureMap; a named ColorData scheme or color function may be supplied."},
  {"ColorFunctionScaling", "True", "Whether ListDensityPlot rescales curvature values before coloring."},
  {"PlotLegends", "Automatic", "Legend specification passed to ListDensityPlot."},
  {"PlotRange", "Automatic", "Displayed curvature range; Automatic follows SymmetricColorRange."},
  {"InterpolationOrder", "1", "Zero or one for the sampled curvature grid."},
  {"Mesh", "None", "Mesh specification passed to ListDensityPlot."},
  {"FrameLabel", "Automatic", "Two frame labels; Automatic uses the selected momentum-coordinate indices."},
  {"FontSize", "18", "Positive size for frame ticks and labels."},
  {"FontFamily", "\"Times\"", "Font family for frame ticks and labels."},
  {"AspectRatio", "1", "Height-to-width ratio of the density plot."},
  {"ImageSize", "Large", "Final Graphics image size."}
};

plotPropertyBerry2DBaseRules = {
  "\"GridSize\" -> 2",
  "\"StepSize\" -> 10^-3",
  "PlotLegends -> None",
  "\"Output\" -> \"Data\""
};

plotPropertyCurrentPage = "plotBerryCurvature2D";
plotPropertyBerry2DSpec = plotPropertyFunctionSpec[
  "plotBerryCurvature2D",
  {
    {"plotBerryCurvature2D[H,centers,occupied,{{u1,u2},{v1,v2}}]",
      "samples occupied-subspace Berry curvature on an oriented momentum slice and returns a signed density plot."}
  },
  {
    StringJoin[
      "centers determines the Hamiltonian momentum dimension and contains one ",
      "fractional center per basis state. Directions selects two distinct ordered ",
      "coordinates; FixedCoordinates sets the remaining slice coordinates."
    ],
    "Each sample calls berryCurvature on a forward Wilson plaquette with the displayed Hermiticity, gap, covariance, and overlap checks. Reversing Directions reverses the sign.",
    StringJoin[
      "The two plotting intervals and GridSize define the sampled coordinates. ",
      "StepSize controls the smaller Wilson plaquettes used at those coordinates ",
      "and need not equal the display-grid spacing."
    ],
    "Output -> \"Data\" returns Directions, FixedCoordinates, GridValues, GridSize, StepSize, Curvature, MaximumAbsoluteCurvature, and Graphics.",
    plotPropertyFailureNote[
      "plotBerryCurvature2D",
      "the centers, ranges, directions, fixed slice, grid, Wilson-plaquette controls, color settings, or display options are invalid, or a curvature sample fails"
    ]
  },
  plotPropertyBerry2DOptionTable,
  {
    plotPropertyEvaluatedExample[
      "Plot the occupied-band Berry curvature of the symmetry-generated P4 Chern model across the full Brillouin zone:",
      plotPropertyChernSetupCode <> "\n" <>
        "plotBerryCurvature2D[\n" <>
        "  chernHamiltonian, chernCenters, 1,\n" <>
        "  {{-Pi,Pi},{-Pi,Pi}},\n" <>
        "  \"GridSize\" -> 31, \"StepSize\" -> 10^-3,\n" <>
        "  PlotLegends -> None\n" <>
        "]"
    ]
  },
  plotPropertyGraphicsOptionExamples[
    {
      {"\"Directions\"", "reverses the oriented plane and curvature sign.", "\"Directions\" -> {2,1}"},
      {"\"FixedCoordinates\"", "sets the full momentum slice explicitly.", "\"FixedCoordinates\" -> {0.,0.}"},
      {"\"GridSize\"", "uses separate point counts along the two displayed axes.", "\"GridSize\" -> {2,3}"},
      {"\"StepSize\"", "uses a full coordinate-step vector for the Wilson plaquette.", "\"StepSize\" -> {0.001,0.002}"},
      {"\"HermitianTolerance\"", "uses a stricter Hamiltonian residual threshold.", "\"HermitianTolerance\" -> 10^-12"},
      {"\"GapTolerance\"", "sets the occupied direct-gap threshold.", "\"GapTolerance\" -> 10^-10"},
      {"\"CovarianceTolerance\"", "sets the plaquette closure threshold.", "\"CovarianceTolerance\" -> 10^-10"},
      {"\"OverlapTolerance\"", "sets the occupied-link singular-value threshold.", "\"OverlapTolerance\" -> 10^-12"},
      {"\"SymmetricColorRange\"", "allows an asymmetric Automatic color range.", "\"SymmetricColorRange\" -> False"},
      {"\"Output\"", "returns the sampled curvature and Graphics together.", "\"Output\" -> \"Data\""},
      {"ColorFunction", "uses SolarColors for the density plot.", "ColorFunction -> \"SolarColors\""},
      {"ColorFunctionScaling", "passes unscaled curvature values to the color function.", "ColorFunctionScaling -> False"},
      {"PlotLegends", "suppresses the color legend.", "PlotLegends -> None"},
      {"PlotRange", "uses the full sampled range explicitly.", "PlotRange -> All"},
      {"InterpolationOrder", "uses piecewise-constant cells.", "InterpolationOrder -> 0"},
      {"Mesh", "draws the sampling mesh.", "Mesh -> All"},
      {"FrameLabel", "sets explicit slice-coordinate labels.", "FrameLabel -> {\"u\",\"v\"}"},
      {"FontSize", "sets a smaller frame text size.", "FontSize -> 14"},
      {"FontFamily", "uses Helvetica for frame text.", "FontFamily -> \"Helvetica\""},
      {"AspectRatio", "uses a four-to-five height ratio.", "AspectRatio -> 4/5"},
      {"ImageSize", "sets an explicit image width.", "ImageSize -> 320"}
    },
    plotPropertyChernSetupCode,
    "plotBerryCurvature2D",
    "chernHamiltonian, chernCenters, 1, {{-0.01,0.01},{-0.01,0.01}}",
    plotPropertyBerry2DBaseRules
  ],
  {
    {"berryCurvature", "MagneticTB/ref/berryCurvature"},
    {"plotBerryCurvature3D", "MagneticTB/ref/plotBerryCurvature3D"},
    {"pointChernNumber", "MagneticTB/ref/pointChernNumber"}
  },
  {"Berry curvature plot", "Wilson plaquette", "Chern band"}
];

plotPropertyBerry3DOptionTable = {
  {"\"GridSize\"", "7", "One integer or three integers, each at least two, for the vector grid."},
  {"\"StepSize\"", "Automatic", "Automatic derives plaquette steps from grid spacing; a positive scalar or three-vector is also accepted."},
  {"\"HermitianTolerance\"", "10^-10", "Nonnegative Hamiltonian Hermiticity tolerance."},
  {"\"GapTolerance\"", "10^-9", "Nonnegative occupied direct-gap threshold."},
  {"\"CovarianceTolerance\"", "10^-8", "Nonnegative plaquette closure-covariance tolerance."},
  {"\"OverlapTolerance\"", "10^-10", "Nonnegative occupied-link singular-value threshold."},
  {"\"VectorScale\"", "Automatic", "Automatic uses 0.45 times the smallest grid spacing; a positive value sets maximum arrow length."},
  {"\"MagnitudeThreshold\"", "0", "Nonnegative threshold below which curvature vectors are hidden."},
  {"\"ArrowheadSize\"", "0.025", "Positive Arrowheads size."},
  {"\"ArrowThickness\"", "1.6", "Positive AbsoluteThickness for arrows."},
  {"\"Output\"", "\"Graphics\"", "Graphics returns the vector plot; Data also returns the complete vector grid, magnitudes, and Graphics."},
  {"ColorFunction", "\"SolarColors\"", "Named ColorData scheme or color function applied to relative vector magnitude."},
  {"Axes", "True", "Whether the three momentum axes are displayed."},
  {"AxesLabel", "Automatic", "Three axis labels; Automatic uses k1, k2, and k3."},
  {"FontSize", "14", "Positive base font size for the Graphics3D object."},
  {"ViewPoint", "Automatic", "Graphics3D viewpoint; the projection remains orthographic and Boxed remains False."},
  {"PlotRange", "All", "Three-dimensional displayed plot range."},
  {"ImageSize", "Large", "Final Graphics3D image size."}
};

plotPropertyBerry3DBaseRules = {
  "\"GridSize\" -> 2",
  "\"StepSize\" -> 10^-3",
  "\"Output\" -> \"Data\""
};

plotPropertyCurrentPage = "plotBerryCurvature3D";
plotPropertyBerry3DSpec = plotPropertyFunctionSpec[
  "plotBerryCurvature3D",
  {
    {"plotBerryCurvature3D[H,centers,occupied,{{k1min,k1max},{k2min,k2max},{k3min,k3max}}]",
      "plots the occupied-subspace Berry-curvature vector field {F23,F31,F12}."}
  },
  {
    "centers must contain one finite fractional three-vector per Hamiltonian basis state. The three ordered ranges and GridSize define the sample points.",
    "At every point, berryCurvature evaluates the oriented planes {2,3}, {3,1}, and {1,2}; these values form the Cartesian display vector {F23,F31,F12}.",
    "Arrow direction shows curvature orientation. Arrow length and color show magnitude; MagnitudeThreshold may hide small vectors without deleting them from Data output.",
    "The result uses Axes as requested, Boxed -> False, and an Orthographic parallel projection. Output -> \"Data\" returns the complete vector and magnitude grids plus Graphics.",
    plotPropertyFailureNote[
      "plotBerryCurvature3D",
      "the centers, ranges, grid, Wilson-plaquette controls, arrow settings, color function, or display options are invalid, or a curvature component fails"
    ]
  },
  plotPropertyBerry3DOptionTable,
  {
    plotPropertyEvaluatedExample[
      "Plot the radial {F23,F31,F12} Berry-curvature field around the BNS 143.3 symmetry-generated Weyl point:",
      plotPropertyWeylSetupCode <> "\n" <>
        "plotBerryCurvature3D[\n" <>
        "  weylHamiltonian, weylCenters, 1,\n" <>
        "  ConstantArray[{-0.2,0.2},3],\n" <>
        "  \"GridSize\" -> 6, \"StepSize\" -> 0.01,\n" <>
        "  ViewPoint -> {2,-2,1}\n" <>
        "]"
    ]
  },
  plotPropertyGraphicsOptionExamples[
    {
      {"\"GridSize\"", "uses separate point counts along the three axes.", "\"GridSize\" -> {2,2,3}"},
      {"\"StepSize\"", "uses one Wilson-plaquette step per coordinate.", "\"StepSize\" -> {0.001,0.001,0.002}"},
      {"\"HermitianTolerance\"", "uses a stricter Hamiltonian residual threshold.", "\"HermitianTolerance\" -> 10^-12"},
      {"\"GapTolerance\"", "sets the occupied direct-gap threshold.", "\"GapTolerance\" -> 10^-10"},
      {"\"CovarianceTolerance\"", "sets the plaquette closure threshold.", "\"CovarianceTolerance\" -> 10^-10"},
      {"\"OverlapTolerance\"", "sets the occupied-link singular-value threshold.", "\"OverlapTolerance\" -> 10^-12"},
      {"\"VectorScale\"", "sets the maximum arrow length explicitly.", "\"VectorScale\" -> 0.01"},
      {
        "\"MagnitudeThreshold\"",
        "hides the weaker of the two nonzero vectors on this small validation grid.",
        "\"MagnitudeThreshold\" -> 0.0001284768"
      },
      {"\"ArrowheadSize\"", "sets the Graphics3D arrowhead size.", "\"ArrowheadSize\" -> 0.02"},
      {"\"ArrowThickness\"", "sets absolute arrow thickness.", "\"ArrowThickness\" -> 1.2"},
      {"\"Output\"", "returns the complete vector grid and Graphics together.", "\"Output\" -> \"Data\""},
      {"ColorFunction", "uses TemperatureMap for relative vector magnitude.", "ColorFunction -> \"TemperatureMap\""},
      {"Axes", "hides the three momentum axes.", "Axes -> False"},
      {"AxesLabel", "sets three explicit momentum labels.", "AxesLabel -> {\"kx\",\"ky\",\"kz\"}"},
      {"FontSize", "sets a smaller base text size.", "FontSize -> 12"},
      {"ViewPoint", "sets an explicit parallel-view direction.", "ViewPoint -> {2,-2,1}"},
      {"PlotRange", "uses all sampled points explicitly.", "PlotRange -> All"},
      {"ImageSize", "sets an explicit image width.", "ImageSize -> 320"}
    },
    plotPropertyWeylSetupCode,
    "plotBerryCurvature3D",
    "weylHamiltonian, weylCenters, 1, ConstantArray[{-0.01,0.01},3]",
    plotPropertyBerry3DBaseRules
  ],
  {
    {"berryCurvature", "MagneticTB/ref/berryCurvature"},
    {"plotBerryCurvature2D", "MagneticTB/ref/plotBerryCurvature2D"},
    {"pointChernNumber", "MagneticTB/ref/pointChernNumber"}
  },
  {"Berry curvature vector field", "orthographic plot", "Weyl point"}
];

plotPropertyRealSpaceOptionTable = {
  {"\"Aggregation\"", "\"Atom\"", "Atom combines basis states at the same Cartesian position; Orbital keeps every matrix basis state separate."},
  {"\"Normalize\"", "True", "Whether the input state is normalized before probabilities are calculated."},
  {"\"WeightThreshold\"", "10^-8", "Nonnegative probability threshold for visible wavefunction spheres."},
  {"\"PositionTolerance\"", "10^-8", "Positive Cartesian-distance tolerance used only for Atom aggregation."},
  {"\"PhaseColoring\"", "False", "Whether visible spheres are colored by representative complex phase instead of probability."},
  {"\"ShowAllSites\"", "True", "Whether faint background spheres are drawn for every site record."},
  {"\"ShowSiteLabels\"", "False", "Whether site indices are drawn beside every site record."},
  {"\"AtomRadius\"", "Automatic", "Automatic uses 0.055 times the shortest lattice vector; a positive value sets background-site radius."},
  {"\"WavefunctionScale\"", "Automatic", "Automatic uses 0.22 times the shortest lattice vector; a positive value sets maximum wavefunction-sphere radius."},
  {"\"Output\"", "\"Graphics\"", "Graphics returns the real-space plot; Data also returns amplitudes, site records, weights, and Graphics."},
  {"ColorFunction", "Automatic", "Automatic uses SolarColors; a named ColorData scheme or color function colors probability when PhaseColoring is False."},
  {"FontSize", "12", "Positive site-label font size."},
  {"ViewPoint", "Automatic", "Graphics3D viewpoint; projection remains orthographic and axes/box remain hidden."},
  {"ImageSize", "Large", "Final Graphics3D image size."}
};

plotPropertyRealSpaceBaseRules = {"\"Output\" -> \"Data\""};

plotPropertyCurrentPage = "plotRealSpaceWavefunction";
plotPropertyRealSpaceSpec = plotPropertyFunctionSpec[
  "plotRealSpaceWavefunction",
  {
    {"plotRealSpaceWavefunction[data,state]",
      "maps a finite-system state onto the matrix-aligned BasisRecords returned by buildRealSpaceHamiltonian."}
  },
  {
    "data must be buildRealSpaceHamiltonian[...,\"Output\"->\"Data\"] with PositionDataAvailable -> True, complete BasisRecords, a numerical Lattice, and WannierCenters metadata.",
    StringJoin[
      "state is a finite numerical vector whose length equals data[\"Dimension\"]. ",
      "Sphere volume represents probability. Atom aggregation combines all basis ",
      "states at the same Cartesian position; Orbital aggregation preserves the ",
      "matrix basis order."
    ],
    "PhaseColoring uses the phase of the largest-amplitude basis state in each aggregate. Weight remains the sum of squared amplitudes in that aggregate.",
    StringJoin[
      "The Graphics3D result always uses Axes -> False, Boxed -> False, and an ",
      "Orthographic parallel projection. Output -> \"Data\" returns normalized ",
      "amplitudes, SiteRecords, visible count, maximum weight, and Graphics."
    ],
    StringJoin[
      "The second example constructs a two-dimensional C4 breathing-square ",
      "higher-order topological insulator directly with gray group 75 and symham. ",
      "With |s1/t1| = 0.2 < 1, it projects one corner ",
      "orbital onto the four lowest-absolute-energy finite-system states to ",
      "select a reproducible representative corner mode. The plot visualizes ",
      "localization; it does not by itself calculate the quadrupole invariant."
    ],
    plotPropertyFailureNote[
      "plotRealSpaceWavefunction",
      StringJoin[
        "the finite-system data lack position metadata, the state has the wrong length or zero norm, an ",
        "aggregation/display option is invalid, or the color function does not return ",
        "colors"
      ]
    ]
  },
  plotPropertyRealSpaceOptionTable,
  {
    plotPropertyEvaluatedExample[
      "Plot the probability and complex phase of a wavefunction on a three-site L-shaped finite system:",
      plotPropertyLWavefunctionCode <> "\n" <>
        "plotRealSpaceWavefunction[\n" <>
        "  plotPropertyLSystem, plotPropertyLState,\n" <>
        "  \"PhaseColoring\" -> True,\n" <>
        "  \"ShowSiteLabels\" -> True,\n" <>
        "  ViewPoint -> {2,-2,1}\n" <>
        "]"
    ],
    plotPropertyEvaluatedExample[
      "Plot a representative corner-localized state of the symmetry-generated two-dimensional C4 higher-order topological insulator:",
      plotPropertyHOTIWavefunctionCode <> "\n" <>
        "plotRealSpaceWavefunction[\n" <>
        "  hotiFinite, hotiCornerState,\n" <>
        "  \"WeightThreshold\" -> 10^-5,\n" <>
        "  ColorFunction -> \"TemperatureMap\",\n" <>
        "  ViewPoint -> {0,0,2}\n" <>
        "]"
    ]
  },
  plotPropertyGraphicsOptionExamples[
    {
      {"\"Aggregation\"", "keeps the two matrix orbitals as separate site records.", "\"Aggregation\" -> \"Orbital\""},
      {"\"Normalize\"", "preserves the input amplitudes instead of normalizing them.", "\"Normalize\" -> False"},
      {"\"WeightThreshold\"", "hides records whose probability is at most the threshold.", "\"WeightThreshold\" -> 0.6"},
      {"\"PositionTolerance\"", "uses an explicit tolerance when grouping coincident orbitals.", "\"PositionTolerance\" -> 10^-6"},
      {"\"PhaseColoring\"", "colors the aggregate by its representative complex phase.", "\"PhaseColoring\" -> True"},
      {"\"ShowAllSites\"", "suppresses faint background site spheres.", "\"ShowAllSites\" -> False"},
      {"\"ShowSiteLabels\"", "draws the site index beside each record.", "\"ShowSiteLabels\" -> True"},
      {"\"AtomRadius\"", "sets the background-site sphere radius explicitly.", "\"AtomRadius\" -> 0.04"},
      {"\"WavefunctionScale\"", "sets the maximum probability-sphere radius explicitly.", "\"WavefunctionScale\" -> 0.2"},
      {"\"Output\"", "returns site records and Graphics together.", "\"Output\" -> \"Data\""},
      {"ColorFunction", "uses TemperatureMap for probability colors.", "ColorFunction -> \"TemperatureMap\""},
      {"FontSize", "sets a smaller site-label size.", "FontSize -> 10"},
      {"ViewPoint", "sets an explicit parallel-view direction.", "ViewPoint -> {2,-2,1}"},
      {"ImageSize", "sets an explicit image width.", "ImageSize -> 320"}
    },
    plotPropertyTwoOrbitalDataCode,
    "plotRealSpaceWavefunction",
    "plotPropertyTwoOrbitalSystem, plotPropertyState",
    plotPropertyRealSpaceBaseRules,
    True
  ],
  {
    {"buildRealSpaceHamiltonian", "MagneticTB/ref/buildRealSpaceHamiltonian"},
    {"plotBerryCurvature3D", "MagneticTB/ref/plotBerryCurvature3D"},
    {"showCrystalStructure", "MagneticTB/ref/showCrystalStructure"}
  },
  {"real-space wavefunction", "finite system", "probability", "phase"}
];

AssociateTo[
  plotPropertyRealSpaceSpec,
  "Expected" -> DeleteDuplicates@Join[
    plotPropertyRealSpaceSpec["Expected"],
    {
      "C4 breathing-square higher-order topological insulator",
      "representative corner-localized state"
    }
  ]
];

testingPlotPropertySpecifications = {
  plotPropertyBuildRealSpaceSpec,
  plotPropertyWilsonSpec,
  plotPropertySurfaceSpec,
  plotPropertyBerry2DSpec,
  plotPropertyBerry3DSpec,
  plotPropertyRealSpaceSpec
};

(* -------------------------------------------------------------------------
   Idempotent tutorial insertion.

   Tagged cells are removed before the refreshed blocks are inserted.  Surface
   spectra are attached to the semi-infinite example, Wilson plots are attached
   to the Wilson/Berry example, and secondary plots remain before History.
   Full and incremental generation share one path.
   ------------------------------------------------------------------------- *)

testingApplyPlotPropertyTutorialUpdate[
    notebook : Notebook[cells_List, options___]
  ] := Module[{
    tag,
    cleanCells,
    wilsonSectionPosition,
    pointSectionPosition,
    historyPosition,
    surfaceFigureCells,
    wilsonFigureCells,
    additionalExampleCells
  },
  tag = "MagneticTBPlotPropertyTutorial";
  cleanCells = Select[cells, FreeQ[#, CellTags -> tag] &];
  wilsonSectionPosition = FirstPosition[
    cleanCells,
    Cell["Wilson and Berry geometry", "Section", ___],
    Missing["WilsonSection"]
  ];
  pointSectionPosition = FirstPosition[
    cleanCells,
    Cell["Point charge and gapless-point search", "Section", ___],
    Missing["PointSection"]
  ];
  historyPosition = FirstPosition[
    cleanCells,
    cell_ /; !FreeQ[cell, Cell["History", "HistorySection", ___]],
    Missing["History"]
  ];
  If[
    MissingQ[wilsonSectionPosition] ||
      MissingQ[pointSectionPosition] ||
      MissingQ[historyPosition],
    Return[$Failed]
  ];

  surfaceFigureCells = {
    Cell[
      StringJoin[
        "The following plotSurfaceSpectrum examples use numerical H(R) blocks produced by ",
        "hoppingData from symmetry-generated Hamiltonians. Their Graphics outputs are generated ",
        "by the current MagneticTB plotting API."
      ],
      "Text", CellTags -> tag
    ],
    Cell["Chern-insulator surface spectrum", "Subsection", CellTags -> tag],
    Cell[
      "The symmetry-generated P4 Chern model has one chiral branch crossing the bulk gap.",
      "Text", CellTags -> tag
    ],
    Cell[BoxData[codeBoxes[
      plotPropertySurfaceSpec["Examples"][[1, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertySurfaceSpec["Examples"][[1, 3]]]], "Output",
      CellTags -> tag],
    Cell["Z2-insulator surface spectrum", "Subsection", CellTags -> tag],
    Cell[
      "The symmetry-generated gray-P4 model has a counterpropagating pair of surface branches.",
      "Text", CellTags -> tag
    ],
    Cell[BoxData[codeBoxes[
      plotPropertySurfaceSpec["Examples"][[2, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertySurfaceSpec["Examples"][[2, 3]]]], "Output",
      CellTags -> tag]
  };

  wilsonFigureCells = {
    Cell[
      StringJoin[
        "The following plotWilsonLoop examples show the sampled Wilson spectra as unjoined points. ",
        "Each Graphics output is generated directly by the current MagneticTB plotting API."
      ],
      "Text", CellTags -> tag
    ],
    Cell["Chern-insulator Wilson flow", "Subsection", CellTags -> tag],
    Cell[
      "The occupied band of the symmetry-generated P4 Chern model winds once across the transverse Brillouin zone.",
      "Text", CellTags -> tag
    ],
    Cell[BoxData[codeBoxes[
      plotPropertyWilsonSpec["Examples"][[1, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyWilsonSpec["Examples"][[1, 3]]]], "Output",
      CellTags -> tag],
    Cell["Z2-insulator Wilson flow", "Subsection", CellTags -> tag],
    Cell[
      StringJoin[
        "The two occupied branches of the symmetry-generated gray-P4 model exchange partners between ",
        "the two time-reversal-invariant endpoints of half the Brillouin zone."
      ],
      "Text", CellTags -> tag
    ],
    Cell[BoxData[codeBoxes[
      plotPropertyWilsonSpec["Examples"][[2, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyWilsonSpec["Examples"][[2, 3]]]], "Output",
      CellTags -> tag]
  };

  additionalExampleCells = {
    Cell["Arbitrary finite shapes and wavefunctions", "Section",
      CellTags -> tag],
    Cell[
      StringJoin[
        "An explicit duplicate-free list of integer cells creates an arbitrary ",
        "open-boundary shape. Data output preserves that order in BasisRecords, ",
        "which is the position-aware input required by plotRealSpaceWavefunction."
      ],
      "Text", CellTags -> tag],
    Cell[BoxData[codeBoxes[
      plotPropertyBuildRealSpaceSpec["Examples"][[1, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyBuildRealSpaceSpec["Examples"][[1, 3]]]], "Output",
      CellTags -> tag],
    Cell[BoxData[codeBoxes[
      plotPropertyRealSpaceSpec["Examples"][[1, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyRealSpaceSpec["Examples"][[1, 3]]]], "Output",
      CellTags -> tag],
    Cell["C4 higher-order topological-insulator corner state", "Subsection",
      CellTags -> tag],
    Cell[
      StringJoin[
        "For |s1/t1| = 0.2 < 1 the open symmetry-generated C4 breathing-square model has four ",
        "near-zero-energy corner modes. Projecting a corner orbital onto that four-state ",
        "subspace selects one representative mode without relying on an ",
        "arbitrary eigenvector inside the nearly degenerate subspace."
      ],
      "Text", CellTags -> tag],
    Cell[BoxData[codeBoxes[
      plotPropertyRealSpaceSpec["Examples"][[2, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyRealSpaceSpec["Examples"][[2, 3]]]], "Output",
      CellTags -> tag],
    Cell["Berry-curvature density and vector plots", "Section",
      CellTags -> tag],
    Cell[
      StringJoin[
        "The two-dimensional plot shows signed curvature on an oriented slice. ",
        "The three-dimensional plot uses the component convention {F23,F31,F12}; ",
        "arrow direction shows orientation and both arrow length and color show ",
        "magnitude."
      ],
      "Text", CellTags -> tag],
    Cell[BoxData[codeBoxes[
      plotPropertyBerry2DSpec["Examples"][[1, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyBerry2DSpec["Examples"][[1, 3]]]], "Output",
      CellTags -> tag],
    Cell[BoxData[codeBoxes[
      plotPropertyBerry3DSpec["Examples"][[1, 2]]]], "Input",
      CellTags -> tag],
    Cell[BoxData[outputBoxes[
      plotPropertyBerry3DSpec["Examples"][[1, 3]]]], "Output",
      CellTags -> tag]
  };
  Notebook[
    Join[
      Take[cleanCells, First[wilsonSectionPosition] - 1],
      surfaceFigureCells,
      Take[
        cleanCells,
        {First[wilsonSectionPosition], First[pointSectionPosition] - 1}
      ],
      wilsonFigureCells,
      Take[
        cleanCells,
        {First[pointSectionPosition], First[historyPosition] - 1}
      ],
      additionalExampleCells,
      Drop[cleanCells, First[historyPosition] - 1]
    ],
    options
  ]
];

(* Keep model construction separate from the plotting workflow.  Each model
   is generated once with msgop -> init -> symham, while every visible example
   below is one direct MagneticTB property or plotting call. *)
ClearAll[
  tutorialCodeBefore,
  tutorialCodeFrom,
  tutorialInputContaining,
  tutorialOutputBoxesAfterInput
];

tutorialCodeBefore[code_String, marker_String] := Module[{position},
  position = StringPosition[code, marker, 1];
  If[position === {}, Return[code]];
  StringTrim@StringTake[code, First[position][[1]] - 1]
];

tutorialCodeFrom[code_String, marker_String] := Module[{position},
  position = StringPosition[code, marker, 1];
  If[position === {}, Return[code]];
  StringTrim@StringDrop[code, First[position][[1]] - 1]
];

tutorialInputContaining[
    notebook_Notebook,
    markers : (_String | {__String})
  ] := FirstCase[
  notebook,
  Cell[BoxData[code_String], "Input", ___] /;
      AnyTrue[Flatten@{markers}, StringContainsQ[code, #] &] :> code,
  $Failed,
  Infinity
];

tutorialOutputBoxesAfterInput[
    Notebook[cells_List, ___],
    marker_String
  ] := Module[{position},
  position = FirstPosition[
    cells,
    Cell[BoxData[code_String], "Input", ___] /;
      StringContainsQ[code, marker]
  ];
  If[MissingQ[position], Return[$Failed]];
  FirstCase[
    Drop[cells, First[position]],
    Cell[BoxData[boxes_], "Output", ___] :> boxes,
    $Failed
  ]
];

Clear[testingApplyPlotPropertyTutorialUpdate];
testingApplyPlotPropertyTutorialUpdate[seedNotebook_Notebook] := Module[
  {
    chernSetup,
    z2Setup,
    hotiSetup,
    weylSetup,
    chernWilsonCall,
    z2WilsonCall,
    chernSurfaceCall,
    z2SurfaceCall,
    hotiPlotCall,
    berry2DCall,
    weylPointCall,
    weylPointOutput,
    weylPlotCall,
    weylPlotOutput
  },
  chernSetup = tutorialCodeBefore[
    plotPropertySurfaceSpec["Examples"][[1, 2]],
    "plotSurfaceSpectrum["
  ];
  z2Setup = StringReplace[
    tutorialInputContaining[
      seedNotebook,
      {"c4TCIShells", "z2Shells"}
    ],
    "c4TCI" -> "z2"
  ];
  hotiSetup = tutorialCodeBefore[
    plotPropertyRealSpaceSpec["Examples"][[2, 2]],
    "plotRealSpaceWavefunction["
  ];
  weylSetup = tutorialInputContaining[seedNotebook, "weylOperations ="];

  chernWilsonCall = tutorialCodeFrom[
    plotPropertyWilsonSpec["Examples"][[1, 2]],
    "plotWilsonLoop["
  ];
  z2WilsonCall = StringReplace[
    tutorialCodeFrom[
      plotPropertyWilsonSpec["Examples"][[2, 2]],
      "plotWilsonLoop["
    ],
    "z2Centers" -> "z2SliceCenters"
  ];
  chernSurfaceCall = tutorialCodeFrom[
    plotPropertySurfaceSpec["Examples"][[1, 2]],
    "plotSurfaceSpectrum["
  ];
  z2SurfaceCall = tutorialCodeFrom[
    plotPropertySurfaceSpec["Examples"][[2, 2]],
    "plotSurfaceSpectrum["
  ];
  hotiPlotCall = tutorialCodeFrom[
    plotPropertyRealSpaceSpec["Examples"][[2, 2]],
    "plotRealSpaceWavefunction["
  ];
  berry2DCall = tutorialCodeFrom[
    plotPropertyBerry2DSpec["Examples"][[1, 2]],
    "plotBerryCurvature2D["
  ];
  weylPointCall = tutorialInputContaining[
    seedNotebook,
    "pointChernNumber[weylHamiltonian"
  ];
  weylPointOutput = tutorialOutputBoxesAfterInput[
    seedNotebook,
    "pointChernNumber[weylHamiltonian"
  ];
  weylPlotCall = tutorialInputContaining[
    seedNotebook,
    "plotBerryCurvature3D["
  ];
  weylPlotOutput = tutorialOutputBoxesAfterInput[
    seedNotebook,
    "plotBerryCurvature3D["
  ];

  Notebook[{
    Cell["Real-space Hamiltonians, surfaces, and topology", "Title"],
    Cell[
      StringJoin[
        "Construct each Hamiltonian once, then pass it directly to a MagneticTB plotting function. ",
        "The model-construction groups are collapsed so that the main workflow remains one or two ",
        "lines per figure; expand a group only when you want to inspect or rerun that Hamiltonian."
      ],
      "Text"
    ],
    Cell[BoxData[codeBoxes["Needs[\"MagneticTB`\"]"]], "Input"],

    Cell["Hamiltonians used below", "Section"],
    Cell[
      "Evaluate a model group once before running the figures that use it.",
      "Text"
    ],
    Cell[CellGroupData[{
      Cell["Chern-insulator Hamiltonian", "Subsection"],
      Cell[BoxData[codeBoxes[chernSetup]], "Input"]
    }, Closed]],
    Cell[CellGroupData[{
      Cell["Z2-insulator Hamiltonian", "Subsection"],
      Cell[BoxData[codeBoxes[z2Setup]], "Input"]
    }, Closed]],
    Cell[CellGroupData[{
      Cell["C4 higher-order topological-insulator Hamiltonian", "Subsection"],
      Cell[BoxData[codeBoxes[hotiSetup]], "Input"]
    }, Closed]],
    Cell[CellGroupData[{
      Cell["Magnetic Weyl Hamiltonian", "Subsection"],
      Cell[BoxData[codeBoxes[weylSetup]], "Input"]
    }, Closed]],

    Cell["Wilson-loop plots", "Section"],
    Cell[
      "One call gives the winding Chern branch; the second gives the Z2 partner exchange.",
      "Text"
    ],
    Cell["Chern insulator", "Subsection"],
    Cell[BoxData[codeBoxes[chernWilsonCall]], "Input"],
    Cell[BoxData[outputBoxes[
      plotPropertyWilsonSpec["Examples"][[1, 3]]]], "Output"],
    Cell["Z2 insulator", "Subsection"],
    Cell[BoxData[codeBoxes[z2WilsonCall]], "Input"],
    Cell[BoxData[outputBoxes[
      plotPropertyWilsonSpec["Examples"][[2, 3]]]], "Output"],

    Cell["Surface spectra", "Section"],
    Cell[
      "plotSurfaceSpectrum uses the H(R) data returned by hoppingData and draws the semi-infinite surface spectrum directly.",
      "Text"
    ],
    Cell["Chern insulator", "Subsection"],
    Cell[BoxData[codeBoxes[chernSurfaceCall]], "Input"],
    Cell[BoxData[outputBoxes[
      plotPropertySurfaceSpec["Examples"][[1, 3]]]], "Output"],
    Cell["Z2 insulator", "Subsection"],
    Cell[BoxData[codeBoxes[z2SurfaceCall]], "Input"],
    Cell[BoxData[outputBoxes[
      plotPropertySurfaceSpec["Examples"][[2, 3]]]], "Output"],

    Cell["C4 higher-order topological-insulator corner state", "Section"],
    Cell[
      "After the finite Hamiltonian and corner state are prepared, one call draws the real-space probability density.",
      "Text"
    ],
    Cell[BoxData[codeBoxes[hotiPlotCall]], "Input"],
    Cell[BoxData[outputBoxes[
      plotPropertyRealSpaceSpec["Examples"][[2, 3]]]], "Output"],

    Cell["Magnetic Weyl point", "Section"],
    Cell[
      "The first line evaluates the node charge; the second draws its Berry-curvature monopole field.",
      "Text"
    ],
    Cell[BoxData[codeBoxes[weylPointCall]], "Input"],
    Cell[BoxData[weylPointOutput], "Output"],
    Cell[BoxData[codeBoxes[weylPlotCall]], "Input"],
    Cell[BoxData[weylPlotOutput], "Output"],

    Cell["Berry-curvature density", "Section"],
    Cell[
      "The same Chern Hamiltonian produces the signed two-dimensional curvature map in one call.",
      "Text"
    ],
    Cell[BoxData[codeBoxes[berry2DCall]], "Input"],
    Cell[BoxData[outputBoxes[
      plotPropertyBerry2DSpec["Examples"][[1, 3]]]], "Output"],

    Cell["MagneticTB help home", "Section"],
    Cell[TextData[{
      ButtonBox[
        "MagneticTB",
        BaseStyle -> "Link",
        ButtonData -> "paclet:MagneticTB/guide/MagneticTB"
      ]
    }], "Text"],
    historyCells[],
    categorizationCells[
      "Tech Note",
      "MagneticTB`",
      "MagneticTB/tutorial/RealSpaceAndTopology"
    ],
    keywordCells[{
      "surface spectrum", "Chern insulator", "Z2 insulator",
      "higher-order topological insulator", "Weyl point", "Berry curvature"
    }]
  },
    TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
    WindowTitle -> "Real-space Hamiltonians, surfaces, and topology",
    StyleDefinitions -> FrontEnd`FileName[
      {"Wolfram"},
      "TechNotePageStylesExt.nb",
      CharacterEncoding -> "UTF-8"
    ]
  ]
];
