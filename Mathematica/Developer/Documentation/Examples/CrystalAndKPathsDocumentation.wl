(* ::Package:: *)

(* Human-maintained CrystalAndKPaths tutorial.
   One 120-degree noncollinear antiferromagnetic kagome model supplies every
   structure, path, and band.  The MSG transports both the reference magnetic
   moment and the induced local spinor; no arrows are added by hand.
   Evaluate the inputs in reading order; do not initialize a separate model
   for each picture.  The short targeted build uses this same source. *)

crystalKagomeInitCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 4},",
  "  wyckoffposition -> {{{1/2, 0, 0}, {1, 0, 0}}},",
  "  symminformation -> msgop[typeIII[{191, 239}]],",
  "  basisFunctions -> {{{1, -1}/Sqrt[2]}},",
  "  RepresentationMode -> \"Induced\"",
  "];"
}, "\n"];

crystalKagomeStructureCode = StringRiffle[{
  "showCrystalStructure[",
  "  \"CellRange\" -> {{0, 1}, {0, 1}, {0, 0}},",
  "  \"MomentScale\" -> 0.3,",
  "  ImageSize -> 400",
  "]"
}, "\n"];

crystalKagomeTopViewCode = StringRiffle[{
  "Show[",
  "  showCrystalStructure[",
  "    \"CellRange\" -> {{0, 2}, {0, 2}, {0, 0}},",
  "    \"MomentScale\" -> 0.3,",
  "    ImageSize -> 400",
  "  ],",
  "  ViewPoint -> {0, 0, 3}, ViewVertical -> {0, 1, 0}",
  "]"
}, "\n"];

crystalKagomePathCode = "kagomePath = standardKPath[]";
crystalKagomePlanePathCode = "planePath = Take[kagomePath, 3]";
crystalKagomeBZCode =
  "showBrillouinZone[\"KPath\" -> planePath, ImageSize -> 400]";

crystalKagomeHamiltonianCode = StringRiffle[{
  "kagomeHamiltonian = Simplify[ComplexExpand[symham[1] + symham[2]]];",
  "kagomeHamiltonian // MatrixForm"
}, "\n"];

crystalKagomeBandCode = StringRiffle[{
  "bandplot[",
  "  planePath, 60, kagomeHamiltonian, {e1 -> 0, t1 -> -1},",
  "  ImageSize -> 400, FontSize -> 16",
  "]"
}, "\n"];

(* Capture the actual formal-API results once, in the order shown below. *)
checkedDocumentationEvaluation[
  "CrystalAndKPaths initialization", ToExpression[crystalKagomeInitCode]
];
crystalKagomeOutputs = Association@Map[
  Function[item,
    With[{label = First[item], code = Last[item]},
      label -> checkedDocumentationEvaluation[label, ToExpression[code]]
    ]
  ],
  {
    "Structure" -> crystalKagomeStructureCode,
    "TopView" -> crystalKagomeTopViewCode,
    "Path" -> crystalKagomePathCode,
    "PlanePath" -> crystalKagomePlanePathCode,
    "BrillouinZone" -> crystalKagomeBZCode,
    "Hamiltonian" -> crystalKagomeHamiltonianCode,
    "Bands" -> crystalKagomeBandCode
  }
];

testingCrystalTutorialPage = Notebook[{
  Cell["Crystal structure, Brillouin zone, and k paths", "Title"],
  gettingStartedHelpHomeCell[],
  Cell[
    StringJoin[
      "Use a kagome antiferromagnet to draw the crystal structure, its magnetic moments, and ",
      "the first Brillouin zone. The three in-plane moments are 120 degrees apart. After ",
      "setting up this model, we obtain a k path and use it to plot the three bands."
    ],
    "Text"
  ],
  Cell["Initialize the model", "Section"],
  Cell[
    StringJoin[
      "Choose magnetic space group 191.239, P6'/m'mm', with a = 1 and c = 4. One ",
      "representative atom generates {1/2,0,0}, {0,1/2,0}, and {1/2,1/2,0}, in this order. ",
      "Their moments in the lattice-vector basis are {1,0,0}, {0,1,0}, and {-1,-1,0}. In ",
      "Cartesian coordinates, these are unit vectors 120 degrees apart and their sum is zero."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomeInitCode]], "Input"],
  Cell[
    StringJoin[
      "Here {1,-1}/Sqrt[2] is one spinor at the representative site, not two independent basis ",
      "states. Induced generates its counterparts at the other two sites. There is one state ",
      "per site, giving a three-band Hamiltonian; keeping both spin states on every site would ",
      "instead give six basis states."
    ],
    "Text"
  ],
  Cell["Crystal structure", "Section"],
  Cell[
    StringJoin[
      "Draw 2 by 2 primitive cells. The picture contains twelve atoms and twelve red ",
      "magnetic-moment arrows. MomentScale changes the arrow length. The gray lines mark cell ",
      "edges, not hopping bonds. The default view uses an orthographic projection."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomeStructureCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["Structure"]]], "Output"],
  Cell[
    StringJoin[
      "View 3 by 3 cells from above to see the repeating 120-degree magnetic pattern. There ",
      "are twenty-seven atoms in this picture. CellRange changes how many cells are drawn, not ",
      "the Hamiltonian:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomeTopViewCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["TopView"]]], "Output"],
  Cell["Conventional path and first Brillouin zone", "Section"],
  Cell[
    StringJoin[
      "Use standardKPath to obtain the conventional path for this hexagonal lattice. The ",
      "returned points are in reciprocal fractional coordinates. The function chooses a ",
      "built-in path from the lattice metric; it does not calculate little-group irreducible ",
      "representations."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomePathCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["Path"]]], "Output"],
  Cell[
    StringJoin[
      "For the in-plane bands, keep the first three segments, Gamma-M-K-Gamma at kz = 0. Use ",
      "the same path in the Brillouin-zone and band plots:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomePlanePathCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["PlanePath"]]], "Output"],
  Cell[BoxData[codeBoxes[crystalKagomeBZCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["BrillouinZone"]]], "Output"],
  Cell[
    StringJoin[
      "The first Brillouin zone is a hexagonal prism, and the colored path is in its central ",
      "plane. In this plot the fractional k points have been converted to Cartesian reciprocal ",
      "coordinates using reciprocal vectors that include 2 Pi."
    ],
    "Text"
  ],
  Cell["Hamiltonian and bands", "Section"],
  Cell[
    StringJoin[
      "Add symham[1] and symham[2] to keep onsite and nearest-neighbor terms. ComplexExpand ",
      "and Simplify rewrite the paired exponentials as cosines for real k. The three rows and ",
      "columns follow the site order given above:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomeHamiltonianCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["Hamiltonian"]]], "Output"],
  Cell[
    StringJoin[
      "The common onsite energy is e1, and t1 is the nearest-neighbor hopping parameter. The ",
      "relative minus sign in the (2,3) entry comes from the induced spinor basis used here. ",
      "These terms contain no interlayer hopping, so the bands are independent of kz. This is ",
      "an ideal kagome model, not a fit to a particular material."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[crystalKagomeBandCode]], "Input"],
  Cell[BoxData[outputBoxes[crystalKagomeOutputs["Bands"]]], "Output"],
  Cell[
    StringJoin[
      "Set e1 = 0 and t1 = -1 and plot the bands. The flat band has energy -2. The two ",
      "dispersive bands meet at K at energy 1, while the flat band touches the lower ",
      "dispersive band at Gamma. The energy unit is the magnitude of the hopping."
    ],
    "Text"
  ],
  historyCells[],
  categorizationCells["Tech Note", "MagneticTB`", "MagneticTB/tutorial/CrystalAndKPaths"],
  keywordCells[{"kagome", "noncollinear antiferromagnet", "induced spinor", "Brillouin zone", "standard k path"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Crystal structure, Brillouin zone, and k paths",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

crystalKagomeExpectedFragments = {
  "Crystal structure, Brillouin zone, and k paths",
  "Initialize the model",
  "Crystal structure",
  "Hamiltonian and bands",
  "kagome",
  "typeIII[{191, 239}]",
  "MagneticTB/tutorial/CrystalAndKPaths",
  "paclet:MagneticTB/guide/MagneticTB"
};
