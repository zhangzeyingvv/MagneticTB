(* ::Package:: *)

(******************************************************************************
  Generate the user-facing fitting comparison notebook.

  The executable cells below are kept as readable strings so the generated
  notebook shows exactly the code evaluated by this generator.  Saved outputs
  are always produced by the current testing runtime; none are hand-written.
******************************************************************************)

ClearAll[inputCell, outputCell, evaluateExample];

scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
projectRoot = ExpandFileName@FileNameJoin[{scriptDirectory, "..", "..", ".."}];
outputDirectory = FileNameJoin[{projectRoot, "output"}];
sourcePath = FileNameJoin[{outputDirectory, "FittingOptimizationDemo.source.wl"}];
runtimeLoader = FileNameJoin[{projectRoot, "MagneticTB", "Kernel", "init.wl"}];

If[!FileExistsQ[runtimeLoader],
  Print["Missing MagneticTB runtime loader: ", runtimeLoader];
  Exit[2]
];

Get[runtimeLoader];

inputCell[code_String] := Cell[
  BoxData[code],
  "Input",
  Evaluatable -> True
];

outputCell[value_] := Cell[
  BoxData@ToBoxes[value, StandardForm],
  "Output"
];

(* These strings are trusted documentation source.  Evaluating the same string
   that is placed in the Input cell prevents the saved output from drifting
   away from the code visible to the reader. *)
evaluateExample[code_String] := Block[
  {Print = Function[Null]},
  ToExpression[code, InputForm]
];

loadCode = StringRiffle[{
  "projectRoot = ExpandFileName[FileNameJoin[{NotebookDirectory[], \"..\"}]];",
  "Get[FileNameJoin[{projectRoot, \"MagneticTB\", \"Kernel\", \"init.wl\"}]];",
  "eigenvalFile = FileNameJoin[{projectRoot, \"Examples\", \"EIGENVAL\"}];"
}, "\n"];

(* The generator already loaded the same runtime from projectRoot. *)
eigenvalFile = FileNameJoin[{projectRoot, "Examples", "EIGENVAL"}];

modelCode = StringRiffle[{
  "vaspBands = vaspEig[eigenvalFile, -0.0072, 1, 21, 23];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 10},",
  "  wyckoffposition -> {{{2/3, 1/3, 0}, {0, 0, 0}}},",
  "  symminformation -> Take[msgop[gray[156]], 6],",
  "  basisFunctions -> {{\"dz2\", \"dxy\", \"dx2-y2\"}},",
  "  InitialBondShells -> 3",
  "];",
  "fittingHamiltonian = symham[1] + symham[2];",
  "fittingPath = {",
  "  {{{0, 0, 0}, {1/2, 0, 0}}, {\"\\[CapitalGamma]\", \"M\"}},",
  "  {{{1/2, 0, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
  "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"\\[CapitalGamma]\"}}",
  "};",
  "initialRules = {",
  "  e1 -> 0.915, e2 -> 0.685,",
  "  t1 -> 0.1, t2 -> 0.09, t3 -> 0.1,",
  "  t4 -> -0.38, t5 -> 0.34, t6 -> 0.245,",
  "  t7 -> 0.565, t8 -> 0.305, t9 -> 0.365",
  "};"
}, "\n"];

fitCode = StringRiffle[{
  "allKPoints = Range[Length[vaspBands]];",
  "gammaNeighborhood = <|",
  "  \"Center\" -> {0, 0, 0},",
  "  \"Radius\" -> 0.35,",
  "  \"Periodic\" -> True",
  "|>;",
  "fitAll = fittingTB[",
  "  fittingHamiltonian, vaspBands, allKPoints, initialRules",
  "];",
  "fitFermi = fittingTB[",
  "  fittingHamiltonian, vaspBands, allKPoints, initialRules,",
  "  \"EnergyWindow\" -> {-0.5, 0.5}",
  "];",
  "fitGamma = fittingTB[",
  "  fittingHamiltonian, vaspBands, allKPoints, initialRules,",
  "  \"KPointNeighborhood\" -> gammaNeighborhood",
  "];",
  "fitLocal = fittingTB[",
  "  fittingHamiltonian, vaspBands, allKPoints, initialRules,",
  "  \"EnergyWindow\" -> {-0.5, 0.5},",
  "  \"KPointNeighborhood\" -> gammaNeighborhood",
  "];",
  "secondNeighborHamiltonian = fittingHamiltonian + symham[3];",
  "secondNeighborParameters = {r1, r2, r3, r4, r5, r6, r7, r8, r9, r10};",
  "fitSecondAll = fittingTB[",
  "  secondNeighborHamiltonian, vaspBands, allKPoints,",
  "  Join[",
  "    fitAll[\"FittedParams\"],",
  "    Thread[secondNeighborParameters -> 0.]",
  "  ]",
  "];",
  "fitSecondFermi = fittingTB[",
  "  secondNeighborHamiltonian, vaspBands, allKPoints,",
  "  Join[",
  "    fitFermi[\"FittedParams\"],",
  "    Thread[secondNeighborParameters -> 0.]",
  "  ],",
  "  \"EnergyWindow\" -> {-0.5, 0.5},",
  "  MaxIterations -> 300",
  "];",
  "fitCases = <|",
  "  \"All data\" -> fitAll,",
  "  \"All data with second neighbors\" -> fitSecondAll,",
  "  \"Fermi window\" -> fitFermi,",
  "  \"Fermi window with second neighbors\" -> fitSecondFermi,",
  "  \"Gamma neighborhood\" -> fitGamma,",
  "  \"Fermi window near Gamma\" -> fitLocal",
  "|>;"
}, "\n"];

comparisonHamiltonian[name_String] := If[
  MemberQ[{
    "All data with second neighbors",
    "Fermi window with second neighbors"
  }, name],
  "secondNeighborHamiltonian",
  "fittingHamiltonian"
];

comparisonCode[name_String] := StringRiffle[{
  "compareBand[",
  "  fittingPath, 30, " <> comparisonHamiltonian[name] <> ",",
  "  fitCases[\"" <> name <> "\"][\"FittedParams\"],",
  "  vaspBands,",
  "  plotRange -> {-1.5, 3}",
  "]"
}, "\n"];

evaluateExample[modelCode];
evaluateExample[fitCode];

comparisonNames = {
  "All data",
  "All data with second neighbors",
  "Fermi window",
  "Fermi window with second neighbors",
  "Gamma neighborhood",
  "Fermi window near Gamma"
};
comparisonOutputs = Association@Table[
  name -> evaluateExample[comparisonCode[name]],
  {name, comparisonNames}
];

notebook = Notebook[
  {
    Cell["MagneticTB 能带拟合：费米能窗与 k 点邻域", "Title"],
    Cell[
      "这个示例使用同一份 VASP EIGENVAL 和由 msgop → init → symham 生成的三轨道 Hamiltonian，直接比较全数据拟合、费米能窗、Γ 点邻域和两者组合的实际效果。",
      "Text"
    ],
    Cell["载入 testing 版", "Section"],
    inputCell[loadCode],
    Cell["构造模型并读取 VASP 能带", "Section"],
    Cell[
      "EIGENVAL 中的能量已由 vaspEig 减去费米能级，因此零能量就是费米能级。",
      "Text"
    ],
    inputCell[modelCode],
    Cell["四种拟合方式", "Section"],
    Cell[
      "EnergyWindow 只保留窗口内的参考能量；KPointNeighborhood 在分数倒坐标中选择Γ 点周围的周期邻域。所有选择器取交集。",
      "Text"
    ],
    inputCell[fitCode],
    Cell["费米能级附近的实际曲线", "Section"],
    Cell[
      "下面六张图使用完全相同的 -1.5 至 3 eV 纵轴范围，完整显示 EIGENVAL 中选取的三条 VASP 能带。紫色为紧束缚能带，蓝色为 VASP 能带。",
      "Text"
    ],
    Cell[
      "symham[1] 是 onsite 层，symham[2] 是最近邻 hopping，symham[3] 是第二近邻 hopping。下面先比较最近邻与加入第二近邻后的结果。",
      "Text"
    ],
    Cell["全数据拟合：最近邻", "Subsection"],
    inputCell[comparisonCode["All data"]],
    outputCell[comparisonOutputs["All data"]],
    Cell["全数据拟合：加入第二近邻", "Subsection"],
    inputCell[comparisonCode["All data with second neighbors"]],
    outputCell[comparisonOutputs["All data with second neighbors"]],
    Cell["费米能窗拟合：最近邻", "Subsection"],
    inputCell[comparisonCode["Fermi window"]],
    outputCell[comparisonOutputs["Fermi window"]],
    Cell["费米能窗拟合：加入第二近邻", "Subsection"],
    inputCell[comparisonCode["Fermi window with second neighbors"]],
    outputCell[comparisonOutputs["Fermi window with second neighbors"]],
    Cell["局部 k 点选择", "Section"],
    Cell["只拟合Γ 点半径 0.35 的邻域", "Subsection"],
    inputCell[comparisonCode["Gamma neighborhood"]],
    outputCell[comparisonOutputs["Gamma neighborhood"]],
    Cell["只拟合Γ 点附近的费米能带", "Subsection"],
    inputCell[comparisonCode["Fermi window near Gamma"]],
    outputCell[comparisonOutputs["Fermi window near Gamma"]],
    Cell[
      "组合选择能够把Γ 点附近的费米能带拟合得很紧，但它只约束一小部分数据。若要得到全路径可用的模型，应扩大选择范围、加入第二近邻或更多 hopping 层、增加轨道，或使用 ResidualWeights 保留窗外的约束。",
      "Text"
    ]
  },
  WindowTitle -> "MagneticTB fitting selections and second neighbors",
  StyleDefinitions -> "Default.nb"
];

If[!DirectoryQ[outputDirectory], CreateDirectory[outputDirectory]];
Export[
  sourcePath,
  ToString[notebook, InputForm, PageWidth -> 120],
  "Text",
  CharacterEncoding -> "UTF-8"
];
Print[sourcePath];
Exit[0];
