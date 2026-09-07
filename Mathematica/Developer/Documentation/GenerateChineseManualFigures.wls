(* Figures for Documentation/MagneticTB中文手册.tex.
   Run one group with the real kernel, for example:
   MAGNETICTB_MANUAL_FIGURES=kagome WolframKernel -noinit -script
     /absolute/path/GenerateChineseManualFigures.wls
   This does not build notebooks or run the documentation release suite. *)

projectRoot = DirectoryName[DirectoryName[DirectoryName[
  ExpandFileName[$InputFileName]
]]];
figureDirectory = FileNameJoin[{projectRoot, "Documentation", "figures"}];
group = Environment["MAGNETICTB_MANUAL_FIGURES"];

Needs["JLink`"];
saveFigure[name_, figure_] := Module[{path},
  path = FileNameJoin[{figureDirectory, "manual-" <> name <> ".pdf"}];
  Print["Exporting ", path];
  If[UsingFrontEnd[Export[path, figure, "PDF"]] === $Failed, Quit[1]];
];

If[group === "kagome",
  (* Reuse the actual saved outputs of the unchanged tutorial.
     No model or documentation fixture is replayed. *)
  notebook = Get[FileNameJoin[{
    projectRoot, "Documentation", "English", "Tutorials",
    "CrystalAndKPaths.nb"
  }]];
  outputBoxes = Cases[notebook,
    Cell[BoxData[boxes_], "Output", ___] :> boxes, Infinity
  ];
  Scan[
    Function[item,
      saveFigure[First[item],
        ToExpression[outputBoxes[[Last[item]]], StandardForm]]
    ],
    {
      "kagome-structure" -> 1,
      "kagome-top" -> 2,
      "kagome-bz" -> 5,
      "kagome-band" -> 7
    }
  ];
  Quit[0];
];

Get[FileNameJoin[{projectRoot, "MagneticTB", "Kernel", "init.wl"}]];

If[group === "graphene",
  sgop = msgop[gray[191]];
  init[
    lattice -> {{a, 0, 0}, {-a/2, Sqrt[3] a/2, 0}, {0, 0, c}},
    lattpar -> {a -> 1, c -> 3},
    wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
    symminformation -> sgop,
    basisFunctions -> {{"pz"}}
  ];
  path = Take[standardKPath[], 3];
  Print["Graphene path: ", InputForm[path]];
  saveFigure["graphene-structure", Show[
    showCrystalStructure[
      "CellRange" -> {{0, 2}, {0, 2}, {0, 0}}, ImageSize -> 440
    ], ViewPoint -> {0, 0, 3}, ViewVertical -> {0, 1, 0}
  ]];
  saveFigure["graphene-bz",
    showBrillouinZone["KPath" -> path, ImageSize -> 440]];
  ham = Sum[symham[i], {i, 3}];
  Print["Graphene Hamiltonian: ", InputForm[ham]];
  parameters = {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5};
  saveFigure["graphene-band",
    bandplot[path, 80, ham, parameters, ImageSize -> 440, FontSize -> 16]];
  Quit[0];
];

If[group === "graphene-soc",
  (* The SOC comparison uses the same lattice and positions, with two
     spin states per carbon.  Print the full generated matrix so that
     the manual can identify the SOC parameter from the real result. *)
  sgop = msgop[gray[191]];
  init[
    lattice -> {{a, 0, 0}, {-a/2, Sqrt[3] a/2, 0}, {0, 0, c}},
    lattpar -> {a -> 1, c -> 3},
    wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
    symminformation -> sgop,
    basisFunctions -> {{"pzup", "pzdn"}}
  ];
  hamsoc = Sum[symham[i], {i, 3}];
  Print["Spinful graphene Hamiltonian: ", InputForm[hamsoc]];
  path = Take[standardKPath[], 3];
  saveFigure["graphene-soc-off", bandplot[
    path, 80, hamsoc, {e1 -> 0.05, t1 -> 0.5, r1 -> 0, r2 -> 0.02},
    ImageSize -> 440, FontSize -> 16
  ]];
  saveFigure["graphene-soc-on", bandplot[
    path, 80, hamsoc, {e1 -> 0.05, t1 -> 0.5, r1 -> 0.05, r2 -> 0.02},
    ImageSize -> 440, FontSize -> 16
  ]];
  Quit[0];
];

If[group === "cubic",
  init[
    lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},
    lattpar -> {a -> 1},
    wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},
    symminformation -> msgop[gray[221]],
    basisFunctions -> {{"px"}},
    RepresentationMode -> "Induced"
  ];
  Print["Induced orbitals: ", InputForm[orbitalTable[]]];
  saveFigure["cubic-structure", showCrystalStructure[
    "CellRange" -> {{0, 1}, {0, 1}, {0, 1}}, ImageSize -> 440
  ]];
  hamCubic = Simplify[ComplexExpand[symham[1] + symham[2]]];
  Print["Cubic Hamiltonian: ", InputForm[hamCubic]];
  pathCubic = standardKPath[];
  Print["Cubic path: ", InputForm[pathCubic]];
  saveFigure["cubic-band", bandplot[
    pathCubic, 60, hamCubic, {e1 -> 0, t1 -> -1},
    ImageSize -> 440, FontSize -> 16
  ]];
  Quit[0];
];

Print["Choose one figure group: kagome, graphene, graphene-soc, or cubic."];
Quit[1];
