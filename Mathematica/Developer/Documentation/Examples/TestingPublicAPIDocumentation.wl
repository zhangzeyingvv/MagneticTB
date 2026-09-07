(* Documentation specifications for public APIs added or changed on the
   Mathematica testing line after main@a57ad04c.  This file is read only by
   GenerateDocumentationSources.wls, after the current MagneticTB runtime and
   the common documentation helpers have been loaded. *)

ClearAll[
  evaluatedTestingExample,
  evaluatedTestingSetup,
  testingDocumentationTargetedQ,
  testingCurrentPage,
  testingFunctionSpec,
  testingRealSpaceTutorialLinks,
  testingOptionExample,
  testingOptionExamples,
  testingFailureNote
];

testingRealSpaceTutorialLinks[] := If[
  ValueQ[deferredDocumentationTutorialNames] &&
    MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
  {},
  {{
    "Real-space Hamiltonians, surfaces, and topology",
    "MagneticTB/tutorial/RealSpaceAndTopology"
  }}
];

testingIncrementalPageSetups = <||>;
testingEvaluatedPageSetups = <||>;

ensureTestingIncrementalSetup[] := Module[{setupCode},
  If[
    KeyExistsQ[testingIncrementalPageSetups, testingCurrentPage] &&
      !TrueQ[Lookup[testingEvaluatedPageSetups, testingCurrentPage, False]],
    setupCode = testingIncrementalPageSetups[testingCurrentPage];
    With[{
        label = "TestingPublicAPISetup-" <> testingCurrentPage,
        source = setupCode
      },
      checkedDocumentationEvaluation[
        label,
        ToExpression[source, InputForm]
      ]
    ];
    AssociateTo[testingEvaluatedPageSetups, testingCurrentPage -> True]
  ]
];

evaluatedTestingExample[text_String, code_String] := Module[
  {requestedPage, evaluateQ},
  requestedPage = Environment["MAGNETICTB_DOC_PAGE"];
  evaluateQ = !StringQ[requestedPage] || StringTrim[requestedPage] === "" ||
    FileBaseName[requestedPage] === testingCurrentPage;
  {
    text,
    code,
    If[evaluateQ,
      ensureTestingIncrementalSetup[];
      With[{
          label = "TestingPublicAPI-" <> IntegerString[Hash[code], 36],
          source = code
        },
        checkedDocumentationEvaluation[label, ToExpression[source, InputForm]]
      ],
      Missing["NotEvaluatedForIncrementalPage", testingCurrentPage]
    ]
  }
];

(* Basic Examples often need a short model definition before the public call.
   Evaluate that setup while generating the expected outputs, but do not add a
   meaningless Null output cell to the user-facing notebook. *)
evaluatedTestingSetup[text_String, code_String] := Module[
  {requestedPage, evaluateQ, label, source},
  requestedPage = Environment["MAGNETICTB_DOC_PAGE"];
  evaluateQ = !StringQ[requestedPage] || StringTrim[requestedPage] === "" ||
    FileBaseName[requestedPage] === testingCurrentPage;
  If[evaluateQ,
    If[
      StringQ[requestedPage] && StringTrim[requestedPage] =!= "",
      AssociateTo[testingIncrementalPageSetups, testingCurrentPage -> code];
      ensureTestingIncrementalSetup[],
      label = "TestingPublicAPISetup-" <> IntegerString[Hash[code], 36];
      source = code;
      checkedDocumentationEvaluation[label, ToExpression[source, InputForm]]
    ]
  ];
  {text, code}
];

testingOptionExample[
    option_String,
    explanation_String,
    code_String
  ] := evaluatedTestingExample[
  option <> " - " <> explanation,
  code
];

testingOptionExamples[records_List] := testingOptionExample @@@ records;

testingFunctionSpec[
    name_String,
    usage_List,
    notes_List,
    optionTable_List,
    examples_List,
    optionExamples_List,
    seeAlso_List,
    keywords_List,
    tutorials_: Automatic
  ] := Join[
  <|
    "Name" -> name,
    "Usage" -> usage,
    "Notes" -> notes,
    "OptionTable" -> optionTable,
    "Examples" -> examples,
    "OptionExamples" -> optionExamples,
    "SeeAlso" -> seeAlso,
    "Keywords" -> keywords,
    "Expected" -> (
      StringTrim[#, "\""] & /@ optionTable[[All, 1]]
    )
  |>,
  If[tutorials === Automatic, <||>, <|"Tutorials" -> tutorials|>]
];

testingFailureNote[function_String, conditions_String] :=
  function <> " returns $Failed and issues a message when " <> conditions <>
    ". It does not repair the input or switch to an approximate fallback.";

testingDocumentationTargetedQ[names_List] := Module[{requestedPage},
  requestedPage = Environment["MAGNETICTB_DOC_PAGE"];
  !StringQ[requestedPage] ||
    StringTrim[requestedPage] === "" ||
    MemberQ[names, FileBaseName[requestedPage]]
];

(* -------------------------------------------------------------------------
   Plotting and reciprocal-space paths.

   The setup strings below are copied into each standalone example so every
   page can be replayed in its own fresh kernel.
   ------------------------------------------------------------------------- *)

testingPlotInitCode = StringRiffle[{
  "plotOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}},",
  "    {{1/2, 1/2, 1/2}, {0, 0, 0}}",
  "  },",
  "  symminformation -> plotOperations,",
  "  basisFunctions -> {{\"s\"}, {\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];"
}, "\n"];

(* One-orbital simple-cubic state shared by band-path and Brillouin-zone
   examples.  The Hamiltonian is generated from two prepared bond shells, so
   the initialized model must provide exactly one basis state and two shells. *)
testingCubicInitCode = StringRiffle[{
  "cubicOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> IdentityMatrix[3],",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> cubicOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];"
}, "\n"];

testingBandCode = StringRiffle[{
  testingCubicInitCode,
  "plotHamiltonian = symham[2];",
  "plotRules = {t1 -> 1};"
}, "\n"];

testingExplicitPathCode = StringRiffle[{
  "plotPath = {",
  "  {{{0, 0, 0}, {0, 1/2, 0}}, {\"G\", \"X\"}},",
  "  {{{0, 1/2, 0}, {1/2, 1/2, 0}}, {\"X\", \"M\"}}",
  "};"
}, "\n"];

testingBandOptionRecords[function_String] := {
  {
    "plotRange",
    "restricts only the displayed energy interval; the sampled eigenvalues are unchanged.",
    testingBandCode <> "\n" <> testingExplicitPathCode <> "\n" <>
      function <> "[plotPath, 4, plotHamiltonian, plotRules, plotRange -> {-1, 1}]"
  },
  {
    "yTicks",
    "accepts explicit positions and labels for the vertical frame axis.",
    testingBandCode <> "\n" <> testingExplicitPathCode <> "\n" <>
      function <> "[plotPath, 4, plotHamiltonian, plotRules, " <>
        "yTicks -> {{-1, \"low\"}, {0, \"zero\"}, {1, \"high\"}}]"
  },
  {
    "FontSize",
    "sets the path labels and tick-label size.",
    testingBandCode <> "\n" <> testingExplicitPathCode <> "\n" <>
      function <> "[plotPath, 4, plotHamiltonian, plotRules, FontSize -> 16]"
  },
  {
    "FontFamily",
    "sets one font family for path labels and frame ticks.",
    testingBandCode <> "\n" <> testingExplicitPathCode <> "\n" <>
      function <> "[plotPath, 4, plotHamiltonian, plotRules, FontFamily -> \"Helvetica\"]"
  },
  {
    "ImageSize",
    "sets the final Graphics size without changing the band data.",
    testingBandCode <> "\n" <> testingExplicitPathCode <> "\n" <>
      function <> "[plotPath, 4, plotHamiltonian, plotRules, ImageSize -> 360]"
  }
};

(* The testing implementation's four-argument showband[path,npoint,H,rules]
   form is currently captured by the automatic-path overload.  Until the
   program window resolves that Usage/implementation discrepancy, document
   and replay only the verified three-argument automatic-path form. *)
testingBandOptionRecords["showband"] := {
  {
    "plotRange",
    "restricts only the displayed energy interval; the sampled eigenvalues are unchanged.",
    "showband[4, plotHamiltonian, plotRules, plotRange -> {-1, 1}]"
  },
  {
    "yTicks",
    "accepts explicit positions and labels for the vertical frame axis.",
    "showband[4, plotHamiltonian, plotRules, " <>
      "yTicks -> {{-1, \"low\"}, {0, \"zero\"}, {1, \"high\"}}]"
  },
  {
    "FontSize",
    "sets the path labels and tick-label size.",
    "showband[4, plotHamiltonian, plotRules, FontSize -> 16]"
  },
  {
    "FontFamily",
    "sets one font family for path labels and frame ticks.",
    "showband[4, plotHamiltonian, plotRules, FontFamily -> \"Helvetica\"]"
  },
  {
    "ImageSize",
    "sets the final Graphics size without changing the band data.",
    "showband[4, plotHamiltonian, plotRules, ImageSize -> 360]"
  }
};

testingBandOptionTable = {
  {
    "plotRange",
    "All",
    "Vertical energy range passed to the band Graphics."
  },
  {
    "yTicks",
    "Automatic",
    "Automatic, None, numeric positions, or standard tick specifications for the vertical frame axis."
  },
  {
    "FontSize",
    "24",
    "Positive numeric size used by path and tick labels."
  },
  {
    "FontFamily",
    "\"Times\"",
    "Font family used by path and tick labels."
  },
  {"ImageSize", "Automatic", "Final Graphics image size."}
};

testingCurrentPage = "bandplot";
testingBandplotSpec = testingFunctionSpec[
  "bandplot",
  {
    {"bandplot[path,npoint,H,rules]", "draws numerical eigenvalues of H along an explicit reciprocal-fractional path."},
    {"bandplot[Automatic,npoint,H,rules]", "uses standardKPath[] for the current initialized lattice."}
  },
  {
    StringJoin[
      "Each path segment is {{k1,k2},{label1,label2}} in reciprocal fractional coordinates. MagneticTB ",
      "samples npoint intervals per segment and evaluates the Bloch variables at 2 Pi k."
    ],
    "The result is a framed Graphics object. Automatic path selection requires a current model prepared by init or initfromrep.",
    testingFailureNote["bandplot", "the path, sample count, parameter rules, vertical ticks, or font settings are invalid"]
  },
  testingBandOptionTable,
  {
    evaluatedTestingSetup[
      "Initialize a one-orbital simple-cubic model and define a two-segment reciprocal-space path:",
      testingBandCode <> "\n" <> testingExplicitPathCode
    ],
    evaluatedTestingExample[
      "Draw the one-band dispersion along the explicit path:",
      "bandplot[plotPath, 6, plotHamiltonian, plotRules]"
    ],
    evaluatedTestingExample[
      "After initialization, Automatic selects the built-in conventional path of the recognized cubic lattice:",
      "bandplot[Automatic, 4, plotHamiltonian, plotRules]"
    ]
  },
  testingOptionExamples[testingBandOptionRecords["bandplot"]],
  {{"showband", "MagneticTB/ref/showband"}, {"standardKPath", "MagneticTB/ref/standardKPath"}},
  {
    "band structure",
    "static plot",
    "high-symmetry path"
  },
  {{"Crystal structure, Brillouin zone, and k paths", "MagneticTB/tutorial/CrystalAndKPaths"}}
];

AssociateTo[
  testingIncrementalPageSetups,
  "showband" -> testingBandCode
];
testingCurrentPage = "showband";
testingShowbandSpec = testingFunctionSpec[
  "showband",
  {{"showband[npoint,H,rules]", "draws a static band plot along standardKPath[] for the current initialized lattice."}},
  {
    "The verified three-argument form returns the same Graphics form and shares all display options with bandplot.",
    "A current initialized model is required because the path is selected from the lattice metric. Use bandplot[path,npoint,H,rules] when an explicit path is required.",
    StringJoin[
      "Usage.wl also advertises showband[path,npoint,H,rules], but the current testing implementation ",
      "misroutes that exact four-argument form through the automatic-path overload. This page does not ",
      "present the unresolved form as working behavior."
    ],
    testingFailureNote["showband", "automatic path recognition fails or any shared bandplot input is invalid"]
  },
  testingBandOptionTable,
  {
    evaluatedTestingSetup[
      "Initialize the one-orbital simple-cubic model used by the automatic path:",
      testingBandCode
    ],
    evaluatedTestingExample[
      "Let showband choose the conventional simple-cubic path and draw the dispersion:",
      "showband[4, plotHamiltonian, plotRules]"
    ]
  },
  testingOptionExamples[testingBandOptionRecords["showband"]],
  {{"bandplot", "MagneticTB/ref/bandplot"}, {"standardKPath", "MagneticTB/ref/standardKPath"}},
  {
    "band structure",
    "automatic path",
    "static plot"
  },
  {{"Crystal structure, Brillouin zone, and k paths", "MagneticTB/tutorial/CrystalAndKPaths"}}
];

AssociateTo[
  testingShowbandSpec,
  "Expected" -> DeleteDuplicates@Join[
    testingShowbandSpec["Expected"],
    {"msgop[gray[221]]", "symham[2]"}
  ]
];

AssociateTo[
  testingIncrementalPageSetups,
  "standardKPath" -> testingCubicInitCode
];
testingCurrentPage = "standardKPath";
testingStandardKPathSpec = testingFunctionSpec[
  "standardKPath",
  {{"standardKPath[]", "returns a conventional high-symmetry path for the current lattice in reciprocal fractional coordinates."}},
  {
    "The algorithm recognizes the lattice metric and then reads a built-in path table covering 25 Brillouin-zone variations of the 14 three-dimensional Bravais lattice types.",
    StringJoin[
      "This is not a space-group little-group or irreducible-representation calculation. Nonstandard ",
      "primitive bases and centred-lattice conventions can prevent metric recognition; use an explicitly ",
      "standardized basis or the matching BravaisType value."
    ],
    "The returned path uses Bradley-Cracknell display labels where the correspondence is unambiguous, with Setyawan-Curtarolo labels retained as a fallback.",
    testingFailureNote[
      "standardKPath",
      StringJoin[
        "no current model exists, the evaluated lattice is singular or nonnumeric, or the requested ",
        "BravaisType conflicts with the recognized metric"
      ]
    ]
  },
  {
    {
      "\"BravaisType\"",
      "Automatic",
      "Automatic metric recognition, or the matching recognized family/BZ string such as \"CubicPrimitive\" or \"CUB\"."
    },
    {"\"Tolerance\"", "10^-6", "Positive numerical tolerance used only in metric classification."}
  },
  {
    evaluatedTestingSetup[
      "Initialize a one-orbital simple-cubic model before requesting its conventional path:",
      testingCubicInitCode
    ],
    evaluatedTestingExample[
      "For a simple-cubic lattice, the path contains the conventional Gamma-X-M-Gamma-R-X and M-R segments:",
      "standardKPath[]"
    ]
  },
  testingOptionExamples[{
    {
      "\"BravaisType\"",
      "makes the recognized convention explicit without changing the cubic path.",
      "standardKPath[\"BravaisType\" -> \"CubicPrimitive\"]"
    },
    {
      "\"Tolerance\"",
      "controls only metric classification and leaves this exact cubic lattice unchanged.",
      "standardKPath[\"Tolerance\" -> 10^-7]"
    }
  }],
  {{"showband", "MagneticTB/ref/showband"}, {"showBrillouinZone", "MagneticTB/ref/showBrillouinZone"}},
  {
    "Bravais lattice",
    "high-symmetry path",
    "reciprocal coordinates"
  },
  {{"Crystal structure, Brillouin zone, and k paths", "MagneticTB/tutorial/CrystalAndKPaths"}}
];

AssociateTo[
  testingStandardKPathSpec,
  "Expected" -> DeleteDuplicates@Join[
    testingStandardKPathSpec["Expected"],
    {"msgop[gray[221]]"}
  ]
];

testingCrystalMomentInitCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, c}},",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {{{1/4, 1/4, 1/2}, {0, 0, 1}}},",
  "  symminformation -> msgop[bnsdict[{75, 1}]],",
  "  basisFunctions -> {{\"s\"}}",
  "];"
}, "\n"];

testingCrystalSymbolicMomentInitCode = StringReplace[
  testingCrystalMomentInitCode,
  "{0, 0, 1}" -> "{0, 0, m}"
];

(* Show the magnetic ordering across four cells rather than only a few
   isolated sites.  These code strings belong only to this reference page. *)
testingFerromagneticCrystalCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 1}},",
  "    {{1/2, 1/2, 1/2}, {0, 0, 1}}",
  "  },",
  "  symminformation -> msgop[bnsdict[{75, 1}]],",
  "  basisFunctions -> {{\"s\"}, {\"s\"}}",
  "];",
  "showCrystalStructure[",
  "  \"CellRange\" -> {{0, 1}, {0, 1}, {0, 0}},",
  "  ImageSize -> 400",
  "]"
}, "\n"];

testingAntiferromagneticCrystalCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 1}},",
  "    {{1/2, 1/2, 1/2}, {0, 0, -1}}",
  "  },",
  "  symminformation -> msgop[bnsdict[{75, 1}]],",
  "  basisFunctions -> {{\"s\"}, {\"s\"}}",
  "];",
  "showCrystalStructure[",
  "  \"CellRange\" -> {{0, 1}, {0, 1}, {0, 0}},",
  "  ImageSize -> 400",
  "]"
}, "\n"];

testingFerrimagneticCrystalCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 1}},",
  "    {{1/2, 1/2, 1/3}, {0, 0, 1}},",
  "    {{1/2, 1/2, 2/3}, {0, 0, -1}}",
  "  },",
  "  symminformation -> msgop[bnsdict[{75, 1}]],",
  "  basisFunctions -> {{\"s\"}, {\"s\"}, {\"s\"}}",
  "];",
  "showCrystalStructure[",
  "  \"CellRange\" -> {{0, 1}, {0, 1}, {0, 0}},",
  "  ImageSize -> 400",
  "]"
}, "\n"];

testingNoncollinearCrystalCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{1/4, 1/8, 0}, {1, 0, 0}}",
  "  },",
  "  symminformation -> msgop[bnsdict[{75, 1}]],",
  "  basisFunctions -> {{\"s\"}}",
  "];",
  "showCrystalStructure[",
  "  \"CellRange\" -> {{0, 1}, {0, 1}, {0, 0}},",
  "  ImageSize -> 400",
  "]"
}, "\n"];

AssociateTo[
  testingIncrementalPageSetups,
  "showCrystalStructure" -> testingCrystalMomentInitCode
];
testingCurrentPage = "showCrystalStructure";
testingShowCrystalStructureSpec = testingFunctionSpec[
  "showCrystalStructure",
  {{"showCrystalStructure[]", "draws the current primitive cell, atoms, and nonzero magnetic moments."}},
  {
    "Atomic positions and moments come from the model prepared by init. A moment is transformed with the direct lattice and displayed as a red arrow.",
    "Symmetry expands each Wyckoff seed into its equivalent atoms. CellRange repeats these atoms in translated cells. Sphere colors distinguish Wyckoff orbits; arrows show moment directions, not their magnitudes.",
    StringJoin[
      "The result is Graphics3D with Axes -> False, Boxed -> False, and ViewProjection -> \"Orthographic\".",
      " These three presentation choices are fixed current behavior rather than public options."
    ],
    testingFailureNote[
      "showCrystalStructure",
      StringJoin[
        "the current model lacks valid lattice/atom data, CellRange is malformed, a symbolic moment ",
        "remains unresolved, or a positive size option is invalid"
      ]
    ]
  },
  {
    {
      "\"CellRange\"",
      "{{0,0},{0,0},{0,0}}",
      "Integer cell bounds, or n for translations from -n through n on every axis."
    },
    {
      "\"MomentRules\"",
      "{}",
      "Rules that make symbolic moment components numerical."
    },
    {
      "\"MomentScale\"",
      "Automatic",
      "Automatic or a positive numerical arrow length scale."
    },
    {
      "\"AtomRadius\"",
      "Automatic",
      "Automatic or a positive numerical sphere radius."
    },
    {
      "\"ShowAtomLabels\"",
      "False",
      "Whether to label every atom by orbit and equivalent-atom index."
    },
    {
      "FontSize",
      "14",
      "Font size used when atom labels are shown."
    },
    {"ImageSize", "Large", "Final Graphics3D image size."}
  },
  {
    evaluatedTestingExample[
      "Ferromagnetic order: draw 2 by 2 cells in the ab plane. All eight atomic moments point along +z:",
      testingFerromagneticCrystalCode
    ],
    evaluatedTestingExample[
      "Antiferromagnetic order: the same eight atoms form two sublattices with equal and opposite moments. The red arrows alternate between +z and -z:",
      testingAntiferromagneticCrystalCode
    ],
    evaluatedTestingExample[
      "Ferrimagnetic order: each cell has two +z moments and one -z moment. The four-cell picture contains twelve atoms and has a nonzero net moment:",
      testingFerrimagneticCrystalCode
    ],
    evaluatedTestingExample[
      "Noncollinear order: P4 generates four in-plane moments from one seed. Repeating the cell shows sixteen atoms with moments along +x, -x, +y, and -y:",
      testingNoncollinearCrystalCode
    ]
  },
  testingOptionExamples[{
    {
      "\"CellRange\"",
      "show two neighbouring cells of the noncollinear model above. Each cell contains four atoms:",
      "showCrystalStructure[\"CellRange\" -> {{0, 1}, {0, 0}, {0, 0}}, ImageSize -> 400]"
    },
    {
      "\"MomentRules\"",
      "P4 expands the seed into four atoms in one cell. Set m = -1 to point all four magnetic moments along -z:",
      testingCrystalSymbolicMomentInitCode <>
        "\nshowCrystalStructure[\"MomentRules\" -> {m -> -1}, ImageSize -> 400]"
    },
    {
      "\"MomentScale\"",
      "use the same four-atom arrangement with +z moments and shorten the arrows. This changes only their displayed length:",
      testingCrystalMomentInitCode <> "\nshowCrystalStructure[\"MomentScale\" -> 0.35, ImageSize -> 400]"
    },
    {
      "\"AtomRadius\"",
      "reduce the sphere radius for all four atoms:",
      "showCrystalStructure[\"AtomRadius\" -> 0.06, ImageSize -> 400]"
    },
    {
      "\"ShowAtomLabels\"",
      "label the four equivalent atoms 1.1 through 1.4. The first number identifies the Wyckoff orbit:",
      "showCrystalStructure[\"ShowAtomLabels\" -> True, ImageSize -> 400]"
    },
    {
      "FontSize",
      "enlarge the four atom labels:",
      "showCrystalStructure[\"ShowAtomLabels\" -> True, FontSize -> 18, ImageSize -> 400]"
    },
    {"ImageSize", "show the four-atom cell at a smaller size:", "showCrystalStructure[ImageSize -> 360]"}
  }],
  {{"showBrillouinZone", "MagneticTB/ref/showBrillouinZone"}, {"orbitalTable", "MagneticTB/ref/orbitalTable"}},
  {
    "crystal structure",
    "magnetic moment",
    "orthographic view"
  },
  {{"Crystal structure, Brillouin zone, and k paths", "MagneticTB/tutorial/CrystalAndKPaths"}}
];

AssociateTo[
  testingShowCrystalStructureSpec,
  "Expected" -> DeleteDuplicates@Join[
    testingShowCrystalStructureSpec["Expected"],
    {
      "msgop[bnsdict[{75, 1}]]",
      "Ferromagnetic order",
      "Antiferromagnetic order",
      "Ferrimagnetic order",
      "Noncollinear order"
    }
  ]
];

(* The hexagonal model is shared by Basic Examples and Options.  Other
   crystal systems below have their own complete, visible initialization. *)
testingHexagonalBZInitCode = StringRiffle[{
  "Needs[\"MagneticTB`\"];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 3/2},",
  "  wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},",
  "  symminformation -> msgop[gray[191]],",
  "  basisFunctions -> {{\"pz\"}}",
  "];"
}, "\n"];

AssociateTo[
  testingIncrementalPageSetups,
  "showBrillouinZone" -> testingHexagonalBZInitCode
];
testingCurrentPage = "showBrillouinZone";
testingShowBrillouinZoneSpec = testingFunctionSpec[
  "showBrillouinZone",
  {{"showBrillouinZone[]", "draws the first Brillouin zone and the current conventional high-symmetry path."}},
  {
    StringJoin[
      "Run init with the desired lattice before plotting. The first Brillouin zone depends on the ",
      "lattice, so no Hamiltonian or hopping parameters are needed."
    ],
    "The lattice vectors are the rows of lattice. The plot uses Cartesian reciprocal-space coordinates, with reciprocal vectors including the factor 2 Pi.",
    "\"KPath\" -> Automatic uses standardKPath[]. An explicit path is a list of {{kStart,kEnd},{labelStart,labelEnd}} segments in reciprocal fractional coordinates; equivalent points are folded into the first zone.",
    "The result is Graphics3D with Axes -> False, Boxed -> False, and ViewProjection -> \"Orthographic\". These are fixed current presentation properties.",
    testingFailureNote["showBrillouinZone", "the reciprocal lattice cannot define a bounded first zone, the supplied path is malformed, or TranslationRange/Tolerance is invalid"]
  },
  {
    {
      "\"KPath\"",
      "Automatic",
      "Automatic conventional path, None, or an explicit MagneticTB band path."
    },
    {
      "\"ShowKPath\"",
      "True",
      "Whether the selected path and point labels are drawn."
    },
    {
      "\"TranslationRange\"",
      "2",
      "Positive integer reciprocal-lattice search range used to construct and fold the zone."
    },
    {
      "\"Tolerance\"",
      "10^-8",
      "Positive geometric tolerance for zone vertices and path folding."
    },
    {
      "FontSize",
      "14",
      "Font size used for high-symmetry point labels."
    },
    {"ImageSize", "Large", "Final Graphics3D image size."}
  },
  {
    evaluatedTestingSetup[
      "Hexagonal lattice: initialize two symmetry-related pz sites using gray group 191. The in-plane primitive vectors enclose 120 degrees:",
      testingHexagonalBZInitCode
    ],
    evaluatedTestingExample[
      "The first Brillouin zone is a hexagonal prism. Red lines mark the automatically selected high-symmetry path; Gamma is at the zone center:",
      "showBrillouinZone[]"
    ]
  },
  testingOptionExamples[{
    {
      "\"KPath\"",
      "draw only the G-M-K-G path in the kz = 0 plane of the hexagonal zone. Here G denotes Gamma, and coordinates are fractions of the reciprocal vectors:",
      StringRiffle[{
        "showBrillouinZone[",
        "  \"KPath\" -> {",
        "    {{{0, 0, 0}, {1/2, 0, 0}}, {\"G\", \"M\"}},",
        "    {{{1/2, 0, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
        "    {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"G\"}}",
        "  }",
        "]"
      }, "\n"]
    },
    {
      "\"ShowKPath\"",
      "hide the path and labels to see the hexagonal prism by itself:",
      "showBrillouinZone[\"ShowKPath\" -> False]"
    },
    {
      "\"TranslationRange\"",
      "include more reciprocal-lattice translations. This hexagonal zone already has the same shape at the default value 2:",
      "showBrillouinZone[\"TranslationRange\" -> 3]"
    },
    {
      "\"Tolerance\"",
      "change the geometric comparison tolerance. This well-conditioned hexagonal example keeps the same visible shape:",
      "showBrillouinZone[\"Tolerance\" -> 10^-7]"
    },
    {
      "FontSize",
      "enlarge the high-symmetry point labels on the hexagonal zone:",
      "showBrillouinZone[FontSize -> 18]"
    },
    {"ImageSize", "draw a smaller hexagonal-zone image:", "showBrillouinZone[ImageSize -> 360]"}
  }],
  {{"standardKPath", "MagneticTB/ref/standardKPath"}, {"showCrystalStructure", "MagneticTB/ref/showCrystalStructure"}},
  {
    "Brillouin zone",
    "high-symmetry path",
    "orthographic view"
  },
  {{"Crystal structure, Brillouin zone, and k paths", "MagneticTB/tutorial/CrystalAndKPaths"}}
];

AssociateTo[
  testingShowBrillouinZoneSpec,
  "Applications" -> {
    evaluatedTestingExample[
      "Primitive tetragonal lattice (gray group 123): a = b and c = 3 a/2 give a square prism. Its height in reciprocal space is 2 Pi/c:",
      StringRiffle[{
        "Needs[\"MagneticTB`\"];",
        "init[",
        "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, c}},",
        "  lattpar -> {a -> 1, c -> 3/2},",
        "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
        "  symminformation -> msgop[gray[123]],",
        "  basisFunctions -> {{\"s\"}}",
        "];",
        "showBrillouinZone[ImageSize -> 400]"
      }, "\n"]
    ],
    evaluatedTestingExample[
      "Primitive orthorhombic lattice (gray group 47): three unequal lattice constants give a rectangular box with edge lengths 2 Pi/a, 2 Pi/b, and 2 Pi/c:",
      StringRiffle[{
        "Needs[\"MagneticTB`\"];",
        "init[",
        "  lattice -> {{a, 0, 0}, {0, b, 0}, {0, 0, c}},",
        "  lattpar -> {a -> 1, b -> 3/2, c -> 2},",
        "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
        "  symminformation -> msgop[gray[47]],",
        "  basisFunctions -> {{\"s\"}}",
        "];",
        "showBrillouinZone[ImageSize -> 400]"
      }, "\n"]
    ],
    evaluatedTestingExample[
      "Monoclinic lattice (gray group 10): the first and third primitive vectors enclose 70 degrees. Hide the path to inspect the zone faces; the Wigner-Seitz zone is not simply the skew reciprocal primitive cell:",
      StringRiffle[{
        "Needs[\"MagneticTB`\"];",
        "init[",
        "  lattice -> {",
        "    {a, 0, 0},",
        "    {0, b, 0},",
        "    {c Cos[beta], 0, c Sin[beta]}",
        "  },",
        "  lattpar -> {a -> 1, b -> 5/4, c -> 3/2, beta -> 70 Degree},",
        "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
        "  symminformation -> msgop[gray[10]],",
        "  basisFunctions -> {{\"s\"}}",
        "];",
        "showBrillouinZone[\"KPath\" -> None, ImageSize -> 400]"
      }, "\n"]
    ]
  }
];

AssociateTo[
  testingShowBrillouinZoneSpec,
  "Expected" -> DeleteDuplicates@Join[
    testingShowBrillouinZoneSpec["Expected"],
    {
      "msgop[gray[191]]", "msgop[gray[123]]",
      "msgop[gray[47]]", "msgop[gray[10]]"
    }
  ]
];

(* -------------------------------------------------------------------------
   Hopping, finite-system, surface, and topology specifications.

   Every Hamiltonian below is produced by the public MagneticTB workflow.
   No documentation-only Pauli model, constant matrix, hand-written hopping
   Association, or identity-only placeholder group is used.
   ------------------------------------------------------------------------- *)

testingCubicDataCode = StringRiffle[{
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
  "cubicData = hoppingData[{1,2}, {e1 -> 0., t1 -> 1.}];"
}, "\n"];

testingChainDataCode = StringRiffle[{
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
  "chainData = hoppingData[",
  "  {1,2},",
  "  {e1 -> 0., t1 -> 1.}",
  "];"
}, "\n"];

testingWeylModelCode = StringRiffle[{
  "weylOperations = msgop[bnsdict[{143,3}]];",
  "init[",
  "  lattice -> {{Sqrt[3] a/2,-a/2,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {{{0,0,0},{0,0,1}}},",
  "  symminformation -> weylOperations,",
  "  basisFunctions -> {{{x + I y,0},{0,x - I y}}},",
  "  InitialBondShells -> 4,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "weylBlock = Sum[symham[i], {i,{2,3,4}}][[{1,4},{1,4}]] /. {",
  "  r11 -> -r2, r15 -> -r10, r3 -> -r6, r7 -> -r14,",
  "  s4 -> -s4, t4 -> t1, t8 -> t9, t11 -> -t7",
  "};",
  "weylParameters = {",
  "  r10 -> -0.09, r14 -> 0.25, r2 -> 0, r6 -> 0.5,",
  "  s4 -> -0.38, s8 -> 0, t1 -> 0, t15 -> 0.735,",
  "  t7 -> 0.05, t9 -> 0.015",
  "};",
  "weylHamiltonian[q_List] := N[",
  "  weylBlock /. weylParameters /.",
  "    Thread[{kx,ky,kz} -> q]",
  "];",
  "weylCenters = {{0,0,0},{0,0,1/2}};"
}, "\n"];

testingChernSliceCode = testingWeylModelCode <> StringRiffle[{
  "",
  "chernSlice[{qx_?NumericQ,qy_?NumericQ}] :=",
  "  weylHamiltonian[{qx,qy,0.25}];",
  "chernCenters = weylCenters[[All,1 ;; 2]];",
  "chernParameterPath = {",
  "  {{{0,-Pi},{0,0}},{\"-Pi\",\"0\"}},",
  "  {{{0,0},{0,Pi}},{\"0\",\"Pi\"}}",
  "};"
}, "\n"];

testingC4TopologicalCode = StringRiffle[{
  "z2Operations = msgop[gray[75]];",
  "init[",
  "  lattice -> {{a,0,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 1, c -> 1},",
  "  wyckoffposition -> {",
  "    {{0,0,1/4},{0,0,0}},",
  "    {{0,0,0},{0,0,0}}",
  "  },",
  "  symminformation -> z2Operations,",
  "  basisFunctions -> {{\"px\",\"py\"},{\"px\",\"py\"}},",
  "  InitialBondShells -> 7,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "z2Shells = Association@Table[i -> symham[i], {i,{2,3,4,5,7}}];",
  "z2Rules = {",
  "  t1->0, t2->2.5, r1->0, r2->2,",
  "  s1->0, s2->0, s3->-1, s4->0, s5->0,",
  "  s6->0, s7->0, s8->1, s9->0, s10->0,",
  "  p5n1->0.5, p5n2->0, p5n3->0, p5n4->0.5,",
  "  p7n1->0, p7n2->0, p7n3->0, p7n4->0,",
  "  p7n5->0.25, p7n6->-0.25, p7n7->0.25,",
  "  p7n8->0, p7n9->0, p7n10->0, p7n11->0,",
  "  p7n12->-0.25, p7n13->-0.25, p7n14->-0.25",
  "};",
  "z2Hamiltonian = Total[Values[z2Shells]] /. z2Rules;",
  "z2Data = hoppingData[Keys[z2Shells], z2Rules];",
  "z2Centers = z2Data[\"WannierCenters\"];",
  "z2Slice[{qx_?NumericQ,qy_?NumericQ}] := N[",
  "  z2Hamiltonian /. {kx->qx,ky->qy,kz->0}",
  "];",
  "z2SliceCenters = z2Centers[[All,1 ;; 2]];",
  "z2SurfacePath = {",
  "  {{{-1/2,0},{0,0}},{\"-Pi\",\"0\"}},",
  "  {{{0,0},{1/2,0}},{\"0\",\"Pi\"}}",
  "};"
}, "\n"];

(* The surface Green-function pages use the same symmetry-generated four-band
   model as a C4 topological crystalline insulator.  Give that presentation
   its own descriptive symbols so readers do not mistake it for a generic or
   purely Z2 example; the Wilson-loop pages retain their established names. *)
testingC4TCISurfaceCode = StringReplace[
  testingC4TopologicalCode,
  "z2" -> "c4TCI"
];

testingHoppingDataCode[optionRules_String] := StringRiffle[{
  "cubicOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> cubicOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "symham[2" <> If[optionRules === "", "", ", " <> optionRules] <> "];",
  "nearestHoppings = hoppingData[",
  "  {2},",
  "  {t1 -> 1}" <>
    If[optionRules === "", "", ",\n  " <> optionRules],
  "];",
  "Grid[",
  "  {",
  "    {\"Translations\", nearestHoppings[\"Translations\"]},",
  "    {\"HoppingMatrices\", nearestHoppings[\"HoppingMatrices\"]}",
  "  },",
  "  Frame -> All, Alignment -> Left",
  "]"
}, "\n"];

testingCurrentPage = "hoppingData";
testingHoppingDataSpec = Join[
  testingFunctionSpec[
    "hoppingData",
  {
    {"hoppingData[n,rules]", "returns numerical real-space hopping blocks for cached shells 1 through n."},
    {"hoppingData[{n1,n2,...},rules]", "returns numerical blocks for the explicitly selected cached shells."}
  },
  {
    "Call init and symham first. hoppingData reads the solved real-space shell cache directly; it never reconstructs hopping amplitudes by Fourier-inverting H(k).",
    "The returned Association uses H(R)[a,b]=<0,a|H|R,b> and H(k)=Sum_R H(R) Exp[2 Pi I k.R]. Lattice and WannierCenters are included when the current model supplies them.",
    testingFailureNote[
      "hoppingData",
      StringJoin[
        "the shell selection is invalid or unprepared, parameter values remain symbolic, the requested ",
        "symham option combination was not cached, or cached real-space data are ",
        "inconsistent"
      ]
    ]
  },
  {
    {
      "\"Hermitian\"",
      "True",
      "Must match the Hermitian setting used when the selected shells were solved by symham."
    },
    {
      "\"KernelMethod\"",
      "\"Iterative\"",
      "Must match the exact common-kernel method used by symham."
    },
    {"\"ValidationLevel\"", "\"Basic\"", "Must match the validation level used by symham."}
  },
  {
    evaluatedTestingExample[
      "Solve the nearest-neighbor shell of a one-orbital cubic model and read its numerical real-space blocks from the cache:",
      testingHoppingDataCode[""]
    ]
  },
  {},
  {{"symham", "MagneticTB/ref/symham"}, {"transformHoppings", "MagneticTB/ref/transformHoppings"}},
  {
    "real-space hopping",
    "shell cache",
    "Wannier90 convention"
  },
    testingRealSpaceTutorialLinks[]
  ],
  <|"AllowPartialOptionExamples" -> True|>
];

testingCurrentPage = "transformHoppings";
testingTransformHoppingsSpec = testingFunctionSpec[
  "transformHoppings",
  {{"transformHoppings[data,T]", "rewrites complete numerical hopping blocks in the exact integer cell whose direct-lattice rows are T times the old rows."}},
  {
    "T must be a nonsingular 3 by 3 integer matrix. Exact Smith decomposition enumerates all cell representatives; no real-space cutoff is introduced.",
    "The returned Association has Abs[Det[T]] times the original orbital dimension and, when available, transformed lattice vectors and Wannier centers.",
    testingFailureNote[
      "transformHoppings",
      StringJoin[
        "the hopping Association/path is incomplete, T is not a nonsingular integer matrix, ",
        "center/lattice metadata are malformed, or the Hermiticity residual exceeds the requested ",
        "tolerance"
      ]
    ]
  },
  {
    {
      "\"WannierCenters\"",
      "Automatic",
      "Automatic uses centers stored in data; None omits them; an explicit list supplies one fractional three-vector per input orbital."
    },
    {
      "\"Lattice\"",
      "Automatic",
      "Automatic uses the stored row-vector lattice; None omits it; an explicit real numerical 3 by 3 matrix replaces it."
    },
    {"\"HermiticityTolerance\"", "10^-9", "Nonnegative numerical tolerance for H(R)=H(-R)^dagger before the exact cell transform."}
  },
  {
    evaluatedTestingSetup[
      "Build the tetragonal gray-group 75 hopping model used below:",
      testingChainDataCode
    ],
    evaluatedTestingExample[
      "Double the chain cell along axis 3. The transformed model contains two orbitals per new cell:",
      "transformed = transformHoppings[chainData, DiagonalMatrix[{1, 1, 2}]];\n" <>
        "KeyTake[transformed, {\"NumWannier\", \"CellVolumeFactor\", \"Lattice\", \"WannierCenters\"}]"
    ]
  },
  testingOptionExamples[{
    {
      "\"WannierCenters\"",
      "uses explicit orbital centers instead of metadata stored in the input.",
      "KeyTake[transformHoppings[chainData, DiagonalMatrix[{1,1,2}], \"WannierCenters\" -> {{1/4,0,0}}], {\"NumWannier\", \"WannierCenters\"}]"
    },
    {
      "\"Lattice\"",
      "uses an explicit old direct lattice and returns the transformed row vectors.",
      "KeyTake[transformHoppings[chainData, DiagonalMatrix[{1,1,2}], \"Lattice\" -> DiagonalMatrix[{2,2,3}]], {\"Lattice\"}]"
    },
    {
      "\"HermiticityTolerance\"",
      "checks the input hopping table against the stated numerical tolerance.",
      "KeyTake[transformHoppings[chainData, IdentityMatrix[3], \"HermiticityTolerance\" -> 10^-12], {\"NumWannier\", \"HermitianResidual\"}]"
    }
  }],
  {{"hoppingData", "MagneticTB/ref/hoppingData"}, {"buildBlochHamiltonian", "MagneticTB/ref/buildBlochHamiltonian"}},
  {
    "integer supercell",
    "Smith decomposition",
    "real-space hopping"
  },
  testingRealSpaceTutorialLinks[]
];

AssociateTo[
  testingTransformHoppingsSpec,
  "Expected" -> DeleteDuplicates@Join[
    testingTransformHoppingsSpec["Expected"],
    {"msgop[gray[75]]"}
  ]
];

testingCurrentPage = "buildRealSpaceHamiltonian";
testingBuildRealSpaceHamiltonianSpec = testingFunctionSpec[
  "buildRealSpaceHamiltonian",
  {{"buildRealSpaceHamiltonian[data,{n1,n2,n3}]", "assembles a sparse finite Hamiltonian from numerical hopping blocks."}},
  {
    "Cell indices run from zero to ni-1 along each direct-lattice axis. The default is open on all three axes; periodic directions wrap the target cell index.",
    "The matrix convention is row=(source cell,orbital), column=(target cell,orbital). Output -> \"Data\" also returns the cell ordering and boundary conditions.",
    testingFailureNote[
      "buildRealSpaceHamiltonian",
      StringJoin[
        "data are incomplete, the size or boundary list is invalid, the hopping table is non-Hermitian ",
        "at the stated tolerance, or Output is unsupported"
      ]
    ]
  },
  {
    {
      "\"BoundaryConditions\"",
      "{\"Open\",\"Open\",\"Open\"}",
      "Open or Periodic independently along the three cell axes."
    },
    {
      "\"HermiticityTolerance\"",
      "10^-9",
      "Nonnegative tolerance applied to the input hopping table and assembled matrix."
    },
    {"\"Output\"", "\"Matrix\"", "Matrix returns SparseArray; Data returns an Association containing the matrix and cell ordering."}
  },
  {
    evaluatedTestingExample[
      "Build a finite symmetry-generated cubic sample and plot its sparse Hamiltonian:",
      testingCubicDataCode <> StringRiffle[{
        "",
        "finiteHamiltonian = buildRealSpaceHamiltonian[cubicData,{4,4,1}];",
        "MatrixPlot[Normal[finiteHamiltonian], FrameTicks -> None]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"BoundaryConditions\"",
      "closes the one-cell third direction periodically, so the two opposite hoppings add on the onsite block.",
      testingChainDataCode <> "\nNormal@buildRealSpaceHamiltonian[chainData, {1,1,1}, \"BoundaryConditions\" -> {\"Open\",\"Open\",\"Periodic\"}]"
    },
    {
      "\"HermiticityTolerance\"",
      "uses a stricter residual threshold for an exactly Hermitian table.",
      testingChainDataCode <> "\nNormal@buildRealSpaceHamiltonian[chainData, {1,1,3}, \"HermiticityTolerance\" -> 10^-12]"
    },
    {
      "\"Output\"",
      "returns the cell order and matrix metadata together with the SparseArray.",
      testingChainDataCode <> StringJoin[
        "\nKeyTake[buildRealSpaceHamiltonian[chainData, {1,1,3}, \"Output\" -> \"Data\"], {\"Cells\",",
        "\"Dimension\",\"BoundaryConditions\",\"HermitianResidual\"}]"
      ]
    }
  }],
  {{"buildSlabHamiltonian", "MagneticTB/ref/buildSlabHamiltonian"}, {"buildBlochHamiltonian", "MagneticTB/ref/buildBlochHamiltonian"}},
  {
    "finite Hamiltonian",
    "open boundary",
    "SparseArray"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "buildSlabHamiltonian";
testingBuildSlabHamiltonianSpec = testingFunctionSpec[
  "buildSlabHamiltonian",
  {{"buildSlabHamiltonian[data,{n1,n2,n3},k]", "builds a hybrid-space Hamiltonian with continuous momentum along selected periodic axes and finite cells along the others."}},
  {
    "With PeriodicDirections -> {1,2}, sizes n1 and n2 must be 1 and k={k1,k2} is in the reciprocal fractional basis; axis 3 remains a finite open stack.",
    "One periodic direction produces a ribbon or rod. CellMatrix may first reorient the exact integer calculation cell.",
    testingFailureNote[
      "buildSlabHamiltonian",
      StringJoin[
        "periodic axes are repeated or out of range, a periodic size is not one, momentum dimension is ",
        "wrong, cell data are invalid, or Hermiticity/Output checks ",
        "fail"
      ]
    ]
  },
  {
    {
      "\"PeriodicDirections\"",
      "{1,2}",
      "One or two distinct direct-lattice axis indices carrying Bloch momentum."
    },
    {
      "\"CellMatrix\"",
      "Automatic",
      "Automatic keeps the input cell; a nonsingular integer 3 by 3 matrix reorients or enlarges it first."
    },
    {
      "\"HermiticityTolerance\"",
      "10^-9",
      "Nonnegative tolerance for the hopping table and assembled hybrid matrix."
    },
    {"\"Output\"", "\"Matrix\"", "Matrix returns SparseArray; Data also reports open/periodic axes, cells, and embedded momentum."}
  },
  {
    evaluatedTestingExample[
      "Build a four-layer symmetry-generated cubic slab and plot its hybrid-space Hamiltonian:",
      testingCubicDataCode <> StringRiffle[{
        "",
        "slabHamiltonian = buildSlabHamiltonian[cubicData,{1,1,4},{1/7,2/11}];",
        "MatrixPlot[Normal[slabHamiltonian], FrameTicks -> None]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"PeriodicDirections\"",
      "keeps only axis 1 periodic and makes a finite 3 by 2 cross-section.",
      testingCubicDataCode <> StringJoin[
        "\nKeyTake[buildSlabHamiltonian[cubicData, {1,3,2}, {1/7}, \"PeriodicDirections\" -> {1}, ",
        "\"Output\" -> \"Data\"], {\"Dimension\",\"PeriodicDirections\",",
        "\"OpenDirections\"}]"
      ]
    },
    {
      "\"CellMatrix\"",
      "reorients the calculation axes with an exact integer permutation.",
      testingCubicDataCode <> StringJoin[
        "\nKeyTake[buildSlabHamiltonian[cubicData, {1,1,3}, {0,0}, \"CellMatrix\" -> {{0,1,0},{0,0,",
        "1},{1,0,0}}, \"Output\" -> \"Data\"], {\"Dimension\",",
        "\"CellMatrix\"}]"
      ]
    },
    {
      "\"HermiticityTolerance\"",
      "uses a stricter residual threshold for exact input data.",
      testingCubicDataCode <> "\nNormal@buildSlabHamiltonian[cubicData, {1,1,3}, {0,0}, \"HermiticityTolerance\" -> 10^-12]"
    },
    {
      "\"Output\"",
      "returns the hybrid-space cell ordering and embedded three-component momentum.",
      testingCubicDataCode <> StringJoin[
        "\nKeyTake[buildSlabHamiltonian[cubicData, {1,1,3}, {1/7,2/11}, \"Output\" -> \"Data\"], ",
        "{\"Cells\",\"CrystalMomentum\",\"EmbeddedCrystalMomentum\",",
        "\"Dimension\"}]"
      ]
    }
  }],
  {{"buildRealSpaceHamiltonian", "MagneticTB/ref/buildRealSpaceHamiltonian"}, {"surfaceGreenFunction", "MagneticTB/ref/surfaceGreenFunction"}},
  {
    "slab Hamiltonian",
    "ribbon",
    "hybrid space"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "buildBlochHamiltonian";
testingBuildBlochHamiltonianSpec = testingFunctionSpec[
  "buildBlochHamiltonian",
  {{"buildBlochHamiltonian[data,{k1,k2,k3}]", "sums numerical hopping blocks into the fully periodic Bloch Hamiltonian at reciprocal-fractional momentum k."}},
  {
    "The convention is H(k)=Sum_R H(R) Exp[2 Pi I k.R]. k is expressed in the reciprocal basis of the active cell.",
    "CellMatrix first constructs the exact integer supercell. The matrix dimension is then multiplied by Abs[Det[CellMatrix]].",
    testingFailureNote["buildBlochHamiltonian", "the hopping table, momentum, cell transform, Hermiticity tolerance, or Output value is invalid"]
  },
  {
    {
      "\"CellMatrix\"",
      "Automatic",
      "Automatic keeps the input cell; a nonsingular integer matrix constructs an exact supercell first."
    },
    {
      "\"HermiticityTolerance\"",
      "10^-9",
      "Nonnegative tolerance for input and final matrix Hermiticity."
    },
    {"\"Output\"", "\"Matrix\"", "Matrix returns the Bloch matrix; Data also reports momentum, cell representatives, and residual."}
  },
  {
    evaluatedTestingExample[
      "Plot the symmetry-generated cubic band along Gamma-X:",
      testingCubicDataCode <> StringRiffle[{
        "",
        "blochBands = Table[",
        "  First@Eigenvalues@buildBlochHamiltonian[cubicData,{q,0,0}],",
        "  {q,0,1/2,1/80}",
        "];",
        "ListLinePlot[blochBands, Frame -> True,",
        "  FrameLabel -> {\"Gamma-X\",\"Energy\"}]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"CellMatrix\"",
      "doubles the chain cell and folds two primitive-cell branches into a 2 by 2 Bloch matrix.",
      testingChainDataCode <> "\nbuildBlochHamiltonian[chainData, {1/7,1/11,1/13}, \"CellMatrix\" -> DiagonalMatrix[{1,1,2}]]"
    },
    {
      "\"HermiticityTolerance\"",
      "uses a stricter numerical residual bound.",
      testingCubicDataCode <> "\nbuildBlochHamiltonian[cubicData, {1/7,1/11,1/13}, \"HermiticityTolerance\" -> 10^-12]"
    },
    {
      "\"Output\"",
      "returns the cell and momentum metadata together with the matrix.",
      testingCubicDataCode <> StringJoin[
        "\nKeyTake[buildBlochHamiltonian[cubicData, {1/7,1/11,1/13}, \"Output\" -> \"Data\"], ",
        "{\"Hamiltonian\",\"Dimension\",\"CrystalMomentum\",\"CellVolumeFactor\",",
        "\"HermitianResidual\"}]"
      ]
    }
  }],
  {{"transformHoppings", "MagneticTB/ref/transformHoppings"}, {"buildSlabHamiltonian", "MagneticTB/ref/buildSlabHamiltonian"}},
  {
    "Bloch Hamiltonian",
    "supercell",
    "reciprocal fractional momentum"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "surfaceGreenFunction";
testingSurfaceGreenFunctionSpec = testingFunctionSpec[
  "surfaceGreenFunction",
  {{"surfaceGreenFunction[data,{k1,k2},energy]", "computes the retarded Green function of a semi-infinite surface by principal-layer decimation."}},
  {
    "Cell axes 1 and 2 span the surface; axis 3 is the stacking direction. Longer-range hopping along axis 3 is grouped into an exact principal layer.",
    "The default Positive surface is terminated toward increasing layer index. Output -> \"Data\" includes the H00/H01/H10 blocks, iterations, residuals, and spectral weight.",
    testingFailureNote[
      "surfaceGreenFunction",
      StringJoin[
        "data, surface momentum, real energy, cell transform, decimation controls, surface side, ",
        "Hermiticity, or output form is invalid; nonconvergence is reported rather than ",
        "hidden"
      ]
    ]
  },
  {
    {
      "\"CellMatrix\"",
      "Automatic",
      "Automatic keeps the input axes; a nonsingular integer matrix defines the surface plane and stacking direction first."
    },
    {
      "\"Broadening\"",
      "10^-3",
      "Positive imaginary retarded-energy broadening eta."
    },
    {
      "\"Tolerance\"",
      "10^-10",
      "Positive relative coupling tolerance for principal-layer decimation."
    },
    {
      "\"MaxIterations\"",
      "200",
      "Positive integer decimation iteration limit."
    },
    {
      "\"Surface\"",
      "\"Positive\"",
      "Positive or Negative termination along the stacking axis."
    },
    {
      "\"HermiticityTolerance\"",
      "10^-9",
      "Nonnegative residual tolerance for hopping and principal-layer blocks."
    },
    {"\"Output\"", "\"GreenFunction\"", "GreenFunction, SpectralWeight, or a diagnostic Data Association."}
  },
  {
    evaluatedTestingSetup[
      "Build the symmetry-generated four-band C4 topological crystalline insulator (TCI):",
      testingC4TCISurfaceCode
    ],
    evaluatedTestingExample[
      "Compute the retarded surface Green matrix at the surface Gamma point:",
      "surfaceGreenFunction[c4TCIData,{0,0},0]"
    ],
    evaluatedTestingExample[
      "Plot the surface spectral function of the same C4 TCI:",
      StringRiffle[{
        "plotSurfaceSpectrum[",
        "  c4TCIData, c4TCISurfacePath, {-3,3},",
        "  \"MomentumSubdivisions\" -> 30,",
        "  \"EnergyPoints\" -> 101,",
        "  \"Broadening\" -> 0.05, PlotLegends -> None",
        "]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"CellMatrix\"",
      "uses an explicit exact calculation cell before constructing the principal layers.",
      StringJoin[
        "KeyTake[surfaceGreenFunction[c4TCIData, {0,0}, 0, \"CellMatrix\" -> IdentityMatrix[3], ",
        "\"Output\" -> \"Data\"], {\"PrincipalLayerDimension\",",
        "\"SurfaceMomentum\"}]"
      ]
    },
    {
      "\"Broadening\"",
      "changes the positive retarded imaginary part of energy.",
      "surfaceGreenFunction[c4TCIData,{0,0},0,\"Broadening\"->10^-2]"
    },
    {
      "\"Tolerance\"",
      "sets the decimation coupling convergence threshold.",
      "KeyTake[surfaceGreenFunction[c4TCIData,{0,0},0,\"Tolerance\"->10^-8,\"Output\"->\"Data\"],{\"Iterations\",\"CouplingResidual\"}]"
    },
    {
      "\"MaxIterations\"",
      "sets an explicit upper bound while leaving sufficient iterations for this chain.",
      "KeyTake[surfaceGreenFunction[c4TCIData,{0,0},0,\"MaxIterations\"->300,\"Output\"->\"Data\"],{\"Iterations\",\"CouplingResidual\"}]"
    },
    {
      "\"Surface\"",
      "selects the negative termination of the topological model.",
      "surfaceGreenFunction[c4TCIData,{0,0},0,\"Surface\"->\"Negative\"]"
    },
    {
      "\"HermiticityTolerance\"",
      "uses a stricter residual threshold for exact hopping blocks.",
      StringJoin[
        "KeyTake[surfaceGreenFunction[c4TCIData, {0,0}, 0, \"HermiticityTolerance\" -> 10^-12, ",
        "\"Output\" -> \"Data\"], {\"HoppingHermitianResidual\",",
        "\"BlockHermitianResidual\"}]"
      ]
    },
    {"\"Output\"", "returns the scalar surface spectral weight directly.", "surfaceGreenFunction[c4TCIData,{0,0},0,\"Output\"->\"SpectralWeight\"]"}
  }],
  {{"surfaceSpectralFunction", "MagneticTB/ref/surfaceSpectralFunction"}, {"buildSlabHamiltonian", "MagneticTB/ref/buildSlabHamiltonian"}},
  {
    "surface Green function",
    "principal layer",
    "decimation"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "surfaceSpectralFunction";
testingSurfaceSpectralFunctionSpec = testingFunctionSpec[
  "surfaceSpectralFunction",
  {{"surfaceSpectralFunction[data,{k1,k2},energy]", "returns -Im Tr[Gsurface]/Pi for the semi-infinite surface."}},
  {
    "The geometry, principal-layer construction, and convergence rules are identical to surfaceGreenFunction.",
    "The return value is a real spectral weight, not the Green matrix. Use surfaceGreenFunction with Output -> \"Data\" for block and convergence diagnostics.",
    testingFailureNote["surfaceSpectralFunction", "any shared surface Green-function input or convergence check fails"]
  },
  {
    {
      "\"CellMatrix\"",
      "Automatic",
      "Optional exact integer cell defining surface and stacking axes."
    },
    {
      "\"Broadening\"",
      "10^-3",
      "Positive retarded broadening eta."
    },
    {
      "\"Tolerance\"",
      "10^-10",
      "Positive decimation convergence tolerance."
    },
    {
      "\"MaxIterations\"",
      "200",
      "Positive decimation iteration limit."
    },
    {
      "\"Surface\"",
      "\"Positive\"",
      "Positive or Negative termination along cell axis 3."
    },
    {"\"HermiticityTolerance\"", "10^-9", "Nonnegative residual tolerance for hopping and principal-layer blocks."}
  },
  {
    evaluatedTestingSetup[
      "Build the symmetry-generated four-band C4 topological crystalline insulator (TCI):",
      testingC4TCISurfaceCode
    ],
    evaluatedTestingExample[
      "Plot the semi-infinite surface spectral weight of the C4 TCI:",
      StringRiffle[{
        "plotSurfaceSpectrum[",
        "  c4TCIData, c4TCISurfacePath, {-3,3},",
        "  \"MomentumSubdivisions\" -> 30,",
        "  \"EnergyPoints\" -> 101,",
        "  \"Broadening\" -> 0.05, PlotLegends -> None",
        "]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"CellMatrix\"",
      "uses an explicit exact calculation cell.",
      "surfaceSpectralFunction[c4TCIData,{0,0},0,\"CellMatrix\"->IdentityMatrix[3]]"
    },
    {
      "\"Broadening\"",
      "broadens the surface pole by eta=0.01.",
      "surfaceSpectralFunction[c4TCIData,{0,0},0,\"Broadening\"->10^-2]"
    },
    {
      "\"Tolerance\"",
      "sets the principal-layer coupling convergence threshold.",
      "surfaceSpectralFunction[c4TCIData,{0,0},0,\"Tolerance\"->10^-8]"
    },
    {
      "\"MaxIterations\"",
      "sets an explicit safe iteration bound for this chain.",
      "surfaceSpectralFunction[c4TCIData,{0,0},0,\"MaxIterations\"->300]"
    },
    {
      "\"Surface\"",
      "selects the negative termination.",
      "surfaceSpectralFunction[c4TCIData,{0,0},0,\"Surface\"->\"Negative\"]"
    },
    {
      "\"HermiticityTolerance\"",
      "uses a stricter residual threshold for exact hopping blocks.",
      "surfaceSpectralFunction[c4TCIData,{0,0},0,\"HermiticityTolerance\"->10^-12]"
    }
  }],
  {{"surfaceGreenFunction", "MagneticTB/ref/surfaceGreenFunction"}, {"buildSlabHamiltonian", "MagneticTB/ref/buildSlabHamiltonian"}},
  {
    "surface spectral function",
    "surface density of states",
    "decimation"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingChernModelCode = testingChernSliceCode;
testingWeylCode = testingWeylModelCode;

testingCurrentPage = "wilsonLoop";
testingWilsonLoopSpec = testingFunctionSpec[
  "wilsonLoop",
  {{"wilsonLoop[h,wcc,occ,start,end]", "computes Wilson-loop eigenphases for an isolated occupied subspace along a forward reciprocal-lattice loop."}},
  {
    StringJoin[
      "start and end have the same coordinate dimension and end-start must be 2 Pi times an integer ",
      "reciprocal vector. wcc contains one fractional orbital center per Hamiltonian basis state."
    ],
    StringJoin[
      "Occupied-subspace overlaps are unitarized by SVD and the endpoint is closed by the orbital-center ",
      "sewing matrix. The production calculation does not use an eigenvector-gauge fallback."
    ],
    "The default return is a sorted list of phases divided by Pi. Alternate outputs expose Wannier centers modulo one, Wilson eigenvalues, or a diagnostic Association.",
    testingFailureNote[
      "wilsonLoop",
      StringJoin[
        "the loop is not reciprocal-lattice closed, centers/occupation are inconsistent, H is ",
        "non-Hermitian, the direct gap closes, endpoint covariance fails, an overlap is singular, or an ",
        "option is invalid"
      ]
    ]
  },
  {
    {
      "\"PathSubdivisions\"",
      "50",
      "Positive number of equal forward-path intervals."
    },
    {
      "\"HermitianTolerance\"",
      "10^-10",
      "Nonnegative Hermiticity and final unitarity tolerance."
    },
    {
      "\"GapTolerance\"",
      "10^-9",
      "Nonnegative minimum direct-gap threshold for an occupied subspace smaller than the full Hilbert space."
    },
    {
      "\"CovarianceTolerance\"",
      "10^-8",
      "Nonnegative endpoint sewing-covariance tolerance."
    },
    {
      "\"OverlapTolerance\"",
      "10^-10",
      "Nonnegative lower bound for link singular values."
    },
    {"\"Output\"", "\"PhaseOverPi\"", "PhaseOverPi, WannierCenters, Eigenvalues, or Data."}
  },
  {
    evaluatedTestingSetup[
      "Build a gapped Chern slice of the BNS 143.3 symmetry-generated Weyl model:",
      testingChernSliceCode
    ],
    evaluatedTestingExample[
      "Compute one occupied-band Wilson loop on the Chern slice:",
      "wilsonLoop[chernSlice, chernCenters, 1, {0,0}, {2 Pi,0}]"
    ],
    evaluatedTestingExample[
      "Plot the Wilson phase while the loop is translated across the transverse Brillouin zone:",
      StringRiffle[{
        "plotWilsonLoop[",
        "  chernSlice, chernCenters, 1,",
        "  {0,0}, {2 Pi,0}, chernParameterPath,",
        "  \"ParameterSubdivisions\" -> 24,",
        "  \"LoopSubdivisions\" -> 40",
        "]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"PathSubdivisions\"",
      "uses twelve equal forward links on the Chern slice.",
      "wilsonLoop[chernSlice,chernCenters,1,{0,0},{2 Pi,0},\"PathSubdivisions\"->12]"
    },
    {
      "\"HermitianTolerance\"",
      "uses a stricter Hermiticity and final unitarity threshold.",
      "wilsonLoop[chernSlice,chernCenters,1,{0,0},{2 Pi,0},\"HermitianTolerance\"->10^-12]"
    },
    {
      "\"GapTolerance\"",
      "sets the occupied-band direct-gap threshold.",
      "wilsonLoop[chernSlice,chernCenters,1,{0,0},{2 Pi,0},\"GapTolerance\"->10^-10]"
    },
    {
      "\"CovarianceTolerance\"",
      "sets the endpoint sewing-covariance threshold.",
      "wilsonLoop[chernSlice,chernCenters,1,{0,0},{2 Pi,0},\"CovarianceTolerance\"->10^-10]"
    },
    {
      "\"OverlapTolerance\"",
      "sets the minimum acceptable occupied-link singular value.",
      "wilsonLoop[chernSlice,chernCenters,1,{0,0},{2 Pi,0},\"OverlapTolerance\"->10^-12]"
    },
    {
      "\"Output\"",
      "returns the diagnostic Wilson matrix, residuals, gap, path, and closure vector.",
      StringJoin[
        "KeyTake[wilsonLoop[chernSlice, chernCenters, 1, {0,0}, {2 Pi,0}, ",
        "\"PathSubdivisions\" -> 12, \"Output\" -> \"Data\"], {\"PhasesOverPi\",\"WannierCenters\",",
        "\"WilsonMatrix\",\"ClosureVector\",\"EndpointCovarianceResidual\"}]"
      ]
    }
  }],
  {{"berryPhase", "MagneticTB/ref/berryPhase"}, {"berryCurvature", "MagneticTB/ref/berryCurvature"}},
  {
    "Wilson loop",
    "Wannier center",
    "orbital sewing"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

AssociateTo[
  testingWilsonLoopSpec,
  "Expected" -> DeleteDuplicates@Join[
    testingWilsonLoopSpec["Expected"],
    {"msgop[bnsdict[{143,3}]]", "plotWilsonLoop["}
  ]
];

testingCurrentPage = "berryPhase";
testingBerryPhaseSpec = testingFunctionSpec[
  "berryPhase",
  {{"berryPhase[h,wcc,occ,path]", "computes the total gauge-invariant occupied-subspace Berry phase along an explicitly sampled closed path."}},
  {
    StringJoin[
      "path is an ordered list of momentum points. Its endpoint may equal the start or differ by 2 Pi ",
      "times an integer reciprocal vector; orbital centers provide the closure sewing matrix."
    ],
    "The default phase is in radians on the principal Arg branch. The algorithm uses the same SVD-unitarized forward links and strict covariance checks as wilsonLoop.",
    testingFailureNote[
      "berryPhase",
      StringJoin[
        "the sampled path is open, centers or occupied count are invalid, H is non-Hermitian or gapless, ",
        "endpoint covariance fails, a link is singular, or an option is ",
        "invalid"
      ]
    ]
  },
  {
    {
      "\"HermitianTolerance\"",
      "10^-10",
      "Nonnegative Hamiltonian Hermiticity and Wilson unitarity tolerance."
    },
    {
      "\"GapTolerance\"",
      "10^-9",
      "Nonnegative minimum direct-gap threshold."
    },
    {
      "\"CovarianceTolerance\"",
      "10^-8",
      "Nonnegative endpoint sewing-covariance tolerance."
    },
    {
      "\"OverlapTolerance\"",
      "10^-10",
      "Nonnegative lower bound for occupied-link singular values."
    },
    {"\"Output\"", "\"Phase\"", "Phase, PhaseOverPi, WilsonDeterminant, or Data."}
  },
  {
    evaluatedTestingSetup[
      "Use the same BNS 143.3 Chern slice for the Berry-phase examples:",
      testingChernSliceCode
    ],
    evaluatedTestingExample[
      "Plot the occupied-band Berry phase for reciprocal loops translated across the slice:",
      StringRiffle[{
        "berryPhaseFlow = Table[",
        "  {qy/Pi, berryPhase[",
        "    chernSlice, chernCenters, 1,",
        "    Subdivide[{0,qy},{2 Pi,qy},24]",
        "  ]/Pi},",
        "  {qy,-Pi,Pi,Pi/12}",
        "];",
        "ListPlot[berryPhaseFlow, Frame -> True,",
        "  FrameLabel -> {\"ky/Pi\",\"Berry phase/Pi\"}]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"HermitianTolerance\"",
      "uses a stricter Hamiltonian and final Wilson-matrix tolerance.",
      "berryPhase[chernSlice,chernCenters,1,Subdivide[{0,0},{2 Pi,0},12],\"HermitianTolerance\"->10^-12]"
    },
    {
      "\"GapTolerance\"",
      "sets the minimum direct gap for a partially occupied space.",
      "berryPhase[chernSlice,chernCenters,1,Subdivide[{0,0},{2 Pi,0},12],\"GapTolerance\"->10^-10]"
    },
    {
      "\"CovarianceTolerance\"",
      "sets the endpoint sewing-covariance threshold.",
      "berryPhase[chernSlice,chernCenters,1,Subdivide[{0,0},{2 Pi,0},12],\"CovarianceTolerance\"->10^-10]"
    },
    {
      "\"OverlapTolerance\"",
      "sets the occupied-link singular-value threshold.",
      "berryPhase[chernSlice,chernCenters,1,Subdivide[{0,0},{2 Pi,0},12],\"OverlapTolerance\"->10^-12]"
    },
    {
      "\"Output\"",
      "returns phase, Wilson determinant, residuals, minimum gap, and closure data together.",
      StringJoin[
        "KeyTake[berryPhase[chernSlice, chernCenters, 1, Subdivide[{0,0},{2 Pi,0},12], ",
        "\"Output\" -> \"Data\"], {\"Phase\",\"PhaseOverPi\",",
        "\"WilsonDeterminant\",\"ClosureVector\",\"EndpointCovarianceResidual\"}]"
      ]
    }
  }],
  {{"wilsonLoop", "MagneticTB/ref/wilsonLoop"}, {"berryCurvature", "MagneticTB/ref/berryCurvature"}},
  {
    "Berry phase",
    "closed momentum path",
    "occupied subspace"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "berryCurvature";
testingBerryCurvatureSpec = testingFunctionSpec[
  "berryCurvature",
  {{"berryCurvature[h,wcc,occ,k]", "computes occupied-subspace Berry curvature from the oriented flux through a small Wilson plaquette at k."}},
  {
    "Directions selects the ordered coordinate plane. Reversing the two indices reverses the sign. StepSize may be one positive length or two side lengths.",
    "Output -> \"Flux\" returns the plaquette Berry phase; Curvature divides it by the oriented plaquette area; Data includes the sampled path and Berry-phase diagnostics.",
    testingFailureNote["berryCurvature", "k, Directions, StepSize, shared Berry-phase tolerances, or Output is invalid"]
  },
  {
    {
      "\"Directions\"",
      "{1,2}",
      "Two distinct ordered momentum-coordinate indices."
    },
    {
      "\"StepSize\"",
      "10^-3",
      "One positive plaquette side length or a pair of positive side lengths."
    },
    {
      "\"HermitianTolerance\"",
      "10^-10",
      "Nonnegative Hamiltonian Hermiticity tolerance."
    },
    {
      "\"GapTolerance\"",
      "10^-9",
      "Nonnegative occupied direct-gap threshold."
    },
    {
      "\"CovarianceTolerance\"",
      "10^-8",
      "Nonnegative closure covariance tolerance."
    },
    {
      "\"OverlapTolerance\"",
      "10^-10",
      "Nonnegative occupied-link singular-value threshold."
    },
    {"\"Output\"", "\"Curvature\"", "Curvature, Flux, or Data."}
  },
  {
    evaluatedTestingSetup[
      "Build a gapped Chern slice of the BNS 143.3 symmetry-generated Weyl model:",
      testingChernModelCode
    ],
    evaluatedTestingExample[
      StringJoin[
        "Sample berryCurvature on a momentum grid and display the resulting field. ",
        "The plot shows a density over the Brillouin zone rather than one isolated number:"
      ],
      StringRiffle[{
        "curvatureSamples = Flatten[Table[",
        "  {qx,qy,berryCurvature[",
        "    chernSlice, chernCenters, 1, {qx,qy},",
        "    \"StepSize\" -> 0.02",
        "  ]},",
        "  {qx,-Pi,Pi,Pi/12},{qy,-Pi,Pi,Pi/12}",
        "] ,1];",
        "ListDensityPlot[curvatureSamples,",
        "  FrameLabel -> {\"kx\",\"ky\"}, PlotLegends -> Automatic]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"Directions\"",
      "reverses the plaquette orientation and hence the curvature sign.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"Directions\"->{2,1}]"
    },
    {
      "\"StepSize\"",
      "uses different side lengths along the two ordered directions.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"StepSize\"->{0.002,0.003}]"
    },
    {
      "\"HermitianTolerance\"",
      "uses a stricter Hamiltonian residual tolerance.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"HermitianTolerance\"->10^-12]"
    },
    {
      "\"GapTolerance\"",
      "sets the occupied direct-gap threshold on the plaquette.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"GapTolerance\"->10^-10]"
    },
    {
      "\"CovarianceTolerance\"",
      "sets the plaquette closure covariance threshold.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"CovarianceTolerance\"->10^-10]"
    },
    {
      "\"OverlapTolerance\"",
      "sets the occupied-link singular-value threshold.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"OverlapTolerance\"->10^-12]"
    },
    {
      "\"Output\"",
      "returns the oriented plaquette Berry phase before division by area.",
      "berryCurvature[chernSlice,chernCenters,1,{0.,0.},\"Output\"->\"Flux\"]"
    }
  }],
  {{"berryPhase", "MagneticTB/ref/berryPhase"}, {"pointChernNumber", "MagneticTB/ref/pointChernNumber"}},
  {
    "Berry curvature",
    "Wilson plaquette",
    "Chern band"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

AssociateTo[
  testingIncrementalPageSetups,
  "pointChernNumber" -> testingWeylCode
];
testingCurrentPage = "pointChernNumber";
testingPointChernNumberSpec = testingFunctionSpec[
  "pointChernNumber",
  {{"pointChernNumber[H,occupied,k0,r]", "computes the integer Chern charge enclosed by an oriented cube around a three-dimensional momentum point."}},
  {
    "H may be a numerical function H[k] or a matrix containing the configured momentum symbols. The occupied subspace must remain gapped on the cube surface.",
    "By default the center must be gapless within CenterGapTolerance. The raw surface flux is accepted only when it lies within IntegerTolerance of an integer.",
    testingFailureNote[
      "pointChernNumber",
      StringJoin[
        "the point/radius/occupation/Hamiltonian is invalid, the required center is gapped, the ",
        "enclosing surface gap closes, an overlap is singular, or flux quantization/options ",
        "fail"
      ]
    ]
  },
  {
    {
      "\"SurfaceSubdivisions\"",
      "8",
      "Positive number of subdivisions on each cube edge."
    },
    {
      "\"HermitianTolerance\"",
      "10^-10",
      "Nonnegative Hamiltonian Hermiticity tolerance."
    },
    {
      "\"SurfaceGapTolerance\"",
      "10^-8",
      "Positive minimum direct gap required on the enclosing cube."
    },
    {
      "\"CenterGapTolerance\"",
      "10^-6",
      "Nonnegative maximum center gap when RequireGaplessCenter is True."
    },
    {
      "\"OverlapTolerance\"",
      "10^-10",
      "Nonnegative occupied-overlap singular-value threshold."
    },
    {
      "\"IntegerTolerance\"",
      "0.005",
      "Nonnegative allowed difference between raw flux/2 Pi and the nearest integer."
    },
    {
      "\"RequireGaplessCenter\"",
      "True",
      "Whether k0 must be gapless within CenterGapTolerance."
    },
    {
      "\"MomentumSymbols\"",
      "{kx,ky,kz}",
      "Three distinct symbols used when H is supplied as a matrix."
    },
    {"\"Output\"", "\"ChernNumber\"", "ChernNumber or a diagnostic Data Association."}
  },
  {
    evaluatedTestingSetup[
      "Build the BNS 143.3 symmetry-generated magnetic Weyl model:",
      testingWeylCode
    ],
    evaluatedTestingExample[
      "Compute the charge enclosed around the symmetry-generated Weyl point:",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->4]"
    ],
    evaluatedTestingExample[
      "Plot the Berry-curvature vector field around the same Weyl point:",
      StringRiffle[{
        "plotBerryCurvature3D[",
        "  weylHamiltonian, weylCenters, 1,",
        "  ConstantArray[{-0.2,0.2},3],",
        "  \"GridSize\" -> 6, \"StepSize\" -> 0.01,",
        "  ViewPoint -> {2,-2,1}",
        "]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"SurfaceSubdivisions\"",
      "uses a three-interval mesh on every cube edge.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3]"
    },
    {
      "\"HermitianTolerance\"",
      "uses a stricter Hermiticity threshold.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3,\"HermitianTolerance\"->10^-12]"
    },
    {
      "\"SurfaceGapTolerance\"",
      "requires a finite direct gap on every sampled surface point.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3,\"SurfaceGapTolerance\"->10^-6]"
    },
    {
      "\"CenterGapTolerance\"",
      "sets the maximum accepted gap at the enclosed center.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3,\"CenterGapTolerance\"->10^-8]"
    },
    {
      "\"OverlapTolerance\"",
      "sets the minimum occupied-overlap singular value on a surface triangle.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3,\"OverlapTolerance\"->10^-12]"
    },
    {
      "\"IntegerTolerance\"",
      "sets the allowed flux quantization error.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3,\"IntegerTolerance\"->0.01]"
    },
    {
      "\"RequireGaplessCenter\"",
      "turns off the independent center-gap requirement while retaining all surface checks.",
      "pointChernNumber[weylHamiltonian,1,{0,0,0},0.2,\"SurfaceSubdivisions\"->3,\"RequireGaplessCenter\"->False]"
    },
    {
      "\"MomentumSymbols\"",
      "maps an explicit symbolic matrix to three chosen momentum variables.",
      StringJoin[
        "pointChernNumber[weylBlock /. weylParameters,1,{0,0,0},0.2,",
        "\"SurfaceSubdivisions\"->3,\"MomentumSymbols\"->{kx,ky,kz}]"
      ]
    },
    {
      "\"Output\"",
      "returns raw charge, quantization error, surface gap, mesh size, and overlap diagnostics.",
      StringJoin[
        "KeyTake[pointChernNumber[weylHamiltonian, 1, {0,0,0}, 0.2, \"SurfaceSubdivisions\" -> 3, ",
        "\"Output\" -> \"Data\"], {\"ChernNumber\",\"RawChernNumber\",\"QuantizationError\",",
        "\"MinimumSurfaceGap\",\"TriangleCount\"}]"
      ]
    }
  }],
  {{"findGaplessPoints", "MagneticTB/ref/findGaplessPoints"}, {"berryCurvature", "MagneticTB/ref/berryCurvature"}},
  {
    "Weyl point",
    "Chern charge",
    "oriented surface flux"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "findGaplessPoints";
testingFindGaplessPointsSpec = testingFunctionSpec[
  "findGaplessPoints",
  {{"findGaplessPoints[H,occupied]", "searches a periodic Brillouin zone for zeros of the direct gap between bands occupied and occupied+1."}},
  {
    StringJoin[
      "The function first samples an eigenvalue-only periodic grid, keeps a bounded number of low-gap ",
      "seeds, and refines gap squared locally. CandidateCount is a cap, so the result does not claim ",
      "mathematical completeness."
    ],
    "H may be a numerical momentum function or a symbolic matrix. Output -> \"Data\" reports seeds, failed refinements, minimum gap, points, and CompletenessGuaranteed -> False.",
    testingFailureNote["findGaplessPoints", "the Brillouin zone, grid, occupation, Hamiltonian, refinement controls, momentum symbols, or Output form is invalid"]
  },
  {
    {
      "\"BrillouinZone\"",
      "{{-Pi,Pi},{-Pi,Pi},{-Pi,Pi}}",
      "One to three ordered finite intervals."
    },
    {
      "\"GridSize\"",
      "15",
      "Integer at least two, or one such integer per momentum dimension."
    },
    {
      "\"CandidateCount\"",
      "32",
      "Positive cap on grid seeds sent to local refinement."
    },
    {
      "\"GapTolerance\"",
      "10^-7",
      "Nonnegative final direct-gap threshold for accepting a point."
    },
    {
      "\"MergeTolerance\"",
      "10^-4",
      "Nonnegative periodic distance used to merge duplicate refined points."
    },
    {
      "\"HermitianTolerance\"",
      "10^-10",
      "Nonnegative Hamiltonian Hermiticity tolerance."
    },
    {
      "\"MaxIterations\"",
      "500",
      "Positive local-refinement iteration limit."
    },
    {
      "\"RefinementMethod\"",
      "\"PrincipalAxis\"",
      "PrincipalAxis or QuasiNewton for local minimization of gap squared."
    },
    {
      "\"MomentumSymbols\"",
      "{kx,ky,kz}",
      "Distinct variables used when H is a symbolic matrix."
    },
    {"\"Output\"", "\"Points\"", "Points or a diagnostic Data Association."}
  },
  {
    evaluatedTestingSetup[
      "Build the BNS 143.3 symmetry-generated Weyl Hamiltonian used in the search:",
      testingWeylCode
    ],
    evaluatedTestingExample[
      "Find the Weyl node and display the returned momentum points:",
      StringRiffle[{
        "gaplessPoints = findGaplessPoints[",
        "  weylHamiltonian, 1,",
        "  \"BrillouinZone\" -> ConstantArray[{-0.2,0.2},3],",
        "  \"GridSize\" -> 5, \"CandidateCount\" -> 16,",
        "  \"GapTolerance\" -> 10^-4,",
        "  \"MergeTolerance\" -> 0.01,",
        "  \"MaxIterations\" -> 200",
        "];",
        "ListPointPlot3D[gaplessPoints,",
        "  PlotRange -> ConstantArray[{-0.2,0.2},3],",
        "  AxesLabel -> {\"kx\",\"ky\",\"kz\"},",
        "  PlotStyle -> Directive[Red,PointSize[0.03]]]"
      }, "\n"]
    ]
  },
  testingOptionExamples[{
    {
      "\"BrillouinZone\"",
      "sets the three-dimensional search box around the Weyl point.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.15,0.15},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4]"
    },
    {
      "\"GridSize\"",
      "uses five samples along every momentum direction.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4]"
    },
    {
      "\"CandidateCount\"",
      "caps the refinement seed list at sixteen records.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4]"
    },
    {
      "\"GapTolerance\"",
      "sets the final gap threshold for accepting refined points.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4]"
    },
    {
      "\"MergeTolerance\"",
      "sets the periodic distance used to merge duplicate accepted points.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4,\"MergeTolerance\"->0.01]"
    },
    {
      "\"HermitianTolerance\"",
      "uses a stricter Hamiltonian residual threshold.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4,\"HermitianTolerance\"->10^-12]"
    },
    {
      "\"MaxIterations\"",
      "sets an explicit local-refinement limit for the Weyl search.",
      "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4,\"MaxIterations\"->200]"
    },
    {
      "\"RefinementMethod\"",
      "uses QuasiNewton instead of PrincipalAxis for gap-squared refinement.",
      StringJoin[
        "findGaplessPoints[weylHamiltonian,1,\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],",
        "\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4,",
        "\"RefinementMethod\"->\"QuasiNewton\"]"
      ]
    },
    {
      "\"MomentumSymbols\"",
      "maps the symmetry-generated symbolic matrix to its three momentum variables.",
      StringJoin[
        "findGaplessPoints[weylBlock /. weylParameters,1,",
        "\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],",
        "\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4,",
        "\"MomentumSymbols\"->{kx,ky,kz}]"
      ]
    },
    {
      "\"Output\"",
      "returns the bounded-search diagnostics and explicitly records that completeness is not guaranteed.",
      StringJoin[
        "KeyTake[findGaplessPoints[weylHamiltonian,1,",
        "\"BrillouinZone\"->ConstantArray[{-0.2,0.2},3],",
        "\"GridSize\"->5,\"CandidateCount\"->16,\"GapTolerance\"->10^-4,",
        "\"Output\"->\"Data\"], {\"Points\",\"GridSize\",",
        "\"SeedCount\",\"MinimumGap\",\"RefinementMethod\",\"CompletenessGuaranteed\"}]"
      ]
    }
  }],
  {{"pointChernNumber", "MagneticTB/ref/pointChernNumber"}, {"berryCurvature", "MagneticTB/ref/berryCurvature"}},
  {
    "gapless points",
    "Brillouin-zone search",
    "local refinement"
  },
  {{"Real-space Hamiltonians, surfaces, and topology", "MagneticTB/tutorial/RealSpaceAndTopology"}}
];

testingCurrentPage = "orbitalTable";
testingOrbitalTableSpec = testingFunctionSpec[
  "orbitalTable",
  {{"orbitalTable[]", "displays the Hamiltonian row and column basis prepared by init or initfromrep."}},
  {
    StringJoin[
      "The public return value is a presentation Grid. Its rows show the Hamiltonian index, site and ",
      "Wyckoff-orbit indices, equivalent atom, fractional and Cartesian positions, local orbital, actual ",
      "transported basis state, spin structure, and transport operation."
    ],
    StringJoin[
      "The Grid is only the display layer. MagneticTB still keeps the model data internally as Association ",
      "and List records; those private records are not a second public return format and should not be ",
      "indexed through the Grid."
    ],
    "orbitalTable never solves a bond shell and never changes the current model.",
    testingFailureNote["orbitalTable", "no compatible current model exists or the stored Hamiltonian basis order is inconsistent"]
  },
  {},
  {
    evaluatedTestingExample[
      "Initialize a two-site cubic model and display the exact row/column ordering shared by all Hamiltonians and representation matrices:",
      testingPlotInitCode <> "\norbitalTable[]"
    ]
  },
  {},
  {{"showHamiltonianBasis", "MagneticTB/ref/showHamiltonianBasis"}, {"showCrystalStructure", "MagneticTB/ref/showCrystalStructure"}},
  {
    "orbital order",
    "Hamiltonian basis",
    "Grid"
  },
  {{"Crystal structure, Brillouin zone, and k paths", "MagneticTB/tutorial/CrystalAndKPaths"}}
];

(* -------------------------------------------------------------------------
   Optional MSGCorep interface.

   These examples keep the dependency boundary explicit and use a nontrivial
   magnetic group for both single- and double-valued routes.
   ------------------------------------------------------------------------- *)

testingCorepOperationsCode = StringRiffle[{
  "Needs[\"MSGCorep`\"];",
  "corepOperations = getMSGElemFromMSGCorep[{164, 87}];"
}, "\n"];

testingCorepDoubleCode = StringRiffle[{
  testingCorepOperationsCode,
  "init[",
  "  lattice -> {",
  "    {0, -a, 0},",
  "    {Sqrt[3] a/2, a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> corepOperations,",
  "  basisFunctions -> {{\"sup\", \"sdn\"}},",
  "  InitialBondShells -> 1,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "corepHamiltonian = symham[1];",
  "getTBBandCorep[",
  "  {164, 87}, corepHamiltonian, {e1 -> 0},",
  "  {{0, 0, 0}, {1/3, 1/3, 0}, {0, 0, 1/2}}",
  "]"
}, "\n"];

testingCorepSingleCode = StringReplace[
  testingCorepDoubleCode,
  "basisFunctions -> {{\"sup\", \"sdn\"}}" ->
    "basisFunctions -> {{\"s\"}}"
];

testingCurrentPage = "getMSGElemFromMSGCorep";
testingGetMSGElemFromMSGCorepSpec = testingFunctionSpec[
  "getMSGElemFromMSGCorep",
  {{
    "getMSGElemFromMSGCorep[{number,index}]",
    StringJoin[
      "converts the magnetic-space-group operations of an explicitly loaded MSGCorep application ",
      "to MagneticTB's ordered operation records."
    ]
  }},
  {
    "This interface is defined only for operations returned by MSGCorep. It is not a general converter for arbitrary lists that merely resemble MSGCorep data.",
    StringJoin[
      "Each returned record contains the MSGCorep label, its rotation matrix, fractional translation, and ",
      "MagneticTB's \"F\"/\"T\" unitary flag. The function also prints the primitive lattice vectors ",
      "supplied by SpaceGroupIrep."
    ],
    "MSGCorep and SpaceGroupIrep are optional external Mathematica applications. The user must load MSGCorep explicitly with Needs[\"MSGCorep`\"] before calling this function; MagneticTB does not load it automatically.",
    testingFailureNote["getMSGElemFromMSGCorep", "MSGCorep or its SpaceGroupIrep dependency has not already been loaded"]
  },
  {},
  {
    evaluatedTestingExample[
      "Load the nontrivial magnetic group {164,87}:",
      "Needs[\"MSGCorep`\"];\ngetMSGElemFromMSGCorep[{164, 87}]"
    ]
  },
  {},
  {{"getTBBandCorep", "MagneticTB/ref/getTBBandCorep"}, {"init", "MagneticTB/ref/init"}},
  {
    "magnetic corepresentation",
    "MSGCorep",
    "optional dependency"
  },
  {{"Band corepresentations with MSGCorep", "MagneticTB/tutorial/Corepresentations"}}
];

testingCurrentPage = "getTBBandCorep";
testingGetTBBandCorepSpec = Join[
  testingFunctionSpec[
    "getTBBandCorep",
    {{"getTBBandCorep[MSG,H,rules,kset]", "prints the MSGCorep band-corepresentation table at every supplied reciprocal-fractional k point."}},
    {
      "Call getMSGElemFromMSGCorep[MSG], then use exactly that ordered operation list in init. getTBBandCorep is supported only for this MSGCorep-generated route.",
      "The basis selected by init fixes single- or double-valued treatment: a scalar basis such as s is single valued, while the explicit {sup,sdn} spinor basis is double valued.",
      StringJoin[
        "kset is a nonempty list of three-component reciprocal-fractional points. The function returns Null ",
        "and prints one k-name line and one Grid per point; the Grid, rather than an Association, is the ",
        "public presentation result."
      ],
      "The optional MSGCorep and SpaceGroupIrep applications must be installed, and the user must execute Needs[\"MSGCorep`\"] before either MagneticTB wrapper is called. MagneticTB neither loads MSGCorep automatically nor replaces its calculation with an internal approximation.",
      testingFailureNote["getTBBandCorep", "the dependency is unavailable, the current init session is incompatible, H/rules/kset are malformed, or MSGCorep rejects the model"]
    },
    {},
    {
      evaluatedTestingSetup[
        "For MSG {164,87}, initialize a double-valued spinor model and print corepresentations at Gamma, K, and A:",
        testingCorepDoubleCode
      ]
    },
    {},
    {{"getMSGElemFromMSGCorep", "MagneticTB/ref/getMSGElemFromMSGCorep"}, {"orbitalTable", "MagneticTB/ref/orbitalTable"}},
    {
      "band corepresentation",
      "double group",
      "high-symmetry k point"
    },
    {{"Band corepresentations with MSGCorep", "MagneticTB/tutorial/Corepresentations"}}
  ],
  <|
    "Applications" -> {
      evaluatedTestingSetup[
        "The same nontrivial group with one scalar s orbital follows the single-valued MSGCorep route:",
        testingCorepSingleCode
      ]
    }
  |>
];

(* graylayer and grayrod are lookup tables; their useful public workflow is
   to pass the returned identifier directly to the matching operation reader. *)

testingCurrentPage = "graylayer";
testingGrayLayerSpec = Join[
  testingFunctionSpec[
    "graylayer",
    {{
      "graylayer[n]",
      "returns the database identifier of the gray magnetic layer group associated with layer-group number n."
    }},
    {
      "graylayer is an Association whose integer keys run from 1 through 80.",
      "Use mlgop[graylayer[n]] to read the complete ordered operation list; graylayer[n] alone returns only the three-integer database identifier."
    },
    {},
    {
      evaluatedTestingExample[
        "Look up gray layer group 25 and immediately read its complete eight-operation list:",
        "mlgop[graylayer[25]]"
      ]
    },
    {},
    {
      {"mlgop", "MagneticTB/ref/mlgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"grayrod", "MagneticTB/ref/grayrod"}
    },
    {"gray magnetic layer group", "database identifier", "layer-group operations"}
  ],
  <|
    "Issues" -> {{
      "A number outside 1 through 80 is not a key of the association:",
      "graylayer[81]",
      Missing["KeyAbsent", 81]
    }}
  |>
];

testingCurrentPage = "grayrod";
testingGrayRodSpec = Join[
  testingFunctionSpec[
    "grayrod",
    {{
      "grayrod[n]",
      "returns the database identifier of the gray magnetic rod group associated with rod-group number n."
    }},
    {
      "grayrod is an Association whose integer keys run from 1 through 75.",
      "Use mrgop[grayrod[n]] to read the complete ordered operation list; grayrod[n] alone returns only the three-integer database identifier."
    },
    {},
    {
      evaluatedTestingExample[
        "Look up gray rod group 25 and immediately read its complete eight-operation list:",
        "mrgop[grayrod[25]]"
      ]
    },
    {},
    {
      {"mrgop", "MagneticTB/ref/mrgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"graylayer", "MagneticTB/ref/graylayer"}
    },
    {"gray magnetic rod group", "database identifier", "rod-group operations"}
  ],
  <|
    "Issues" -> {{
      "A number outside 1 through 75 is not a key of the association:",
      "grayrod[76]",
      Missing["KeyAbsent", 76]
    }}
  |>
];

(* hop is intentionally documented as a short user workflow: first construct
   the model Hamiltonian, then either export that matrix or select solved
   neighbour contributions.  The option examples retain only the two choices
   needed to understand file output and orbital centers. *)

testingCurrentPage = "hop";
testingHopSpec = Join[
  testingFunctionSpec[
  "hop",
  {
    {
      "hop[n,rules]",
      "writes the real-space hoppings obtained from solved shells 1 through n in Wannier90 HR format."
    },
    {
      "hop[{n1,n2,...},rules]",
      "writes only the selected solved shells in Wannier90 HR format."
    },
    {
      "hop[H,rules]",
      "Fourier-decomposes a supplied finite exponential Bloch Hamiltonian and writes the resulting real-space matrices in Wannier90 HR format."
    }
  },
  {
    StringJoin[
      "There are two ways to use hop. In hop[n,rules] and hop[{n1,n2,...},rules], first calculate the ",
      "required shells with symham. The integer form exports shells 1 through n, while the list form ",
      "exports only the listed shells. The real-space hopping matrices already obtained by symham are ",
      "written directly."
    ],
    "For hop[H,rules], H must be written as a finite sum of exponential Bloch phases so that integer lattice translations can be identified.",
    StringJoin[
      "Set \"hrExport\" to a path to write wannier90_hr.dat; with None, the text is shown in the notebook. ",
      "Use \"wcc\" -> Automatic when the orbital centers should come from init."
    ]
  },
  {
    {
      "\"Hermitian\"",
      "True",
      "For shell export, select the solution previously obtained with the same Hermitian setting. This option is not used for a supplied matrix H(k)."
    },
    {
      "\"KernelMethod\"",
      "\"Iterative\"",
      "For shell export, select the solution previously obtained with the same \"Iterative\" or \"Stacked\" method."
    },
    {
      "\"ValidationLevel\"",
      "\"Basic\"",
      "For shell export, select the solution previously obtained with the same \"None\", \"Basic\", or \"Full\" validation level."
    },
    {
      "\"hrExport\"",
      "None",
      "None returns and prints Wannier90 text; a file or directory path writes wannier90_hr.dat."
    },
    {
      "\"RealDigits\"",
      "12",
      "Integer at least 6 giving the number of decimal digits written for real and imaginary hopping components."
    },
    {
      "\"wcc\"",
      "None",
      StringJoin[
        "For a supplied H(k), None uses zero centers, Automatic uses the current initialized orbital centers ",
        "when available, and an explicit list supplies one fractional center per orbital."
      ]
    },
    {
      "\"TranslationTolerance\"",
      "10^-9",
      "Positive numerical tolerance used only when extracting integer real-space translations from exponential phases in a supplied H(k)."
    }
  },
  {
    evaluatedTestingExample[
      "Initialize graphene and construct its Hamiltonian from the first three neighbour contributions:",
      simpleInitCode <> "\nham = Sum[symham[i], {i, {1, 2, 3}}]"
    ],
    evaluatedTestingExample[
      "Convert the Hamiltonian directly to a Wannier90 HR file:",
      StringRiffle[{
        "hrFile = FileNameJoin[{$TemporaryDirectory, \"graphene_hr.dat\"}];",
        "hop[",
        "  ham, {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5},",
        "  \"hrExport\" -> hrFile, \"wcc\" -> Automatic",
        "];",
        "Take[Import[hrFile, \"Lines\"], 12]"
      }, "\n"]
    ],
    evaluatedTestingExample[
      "Alternatively, export only a selected neighbour contribution:",
      StringRiffle[{
        "hop[",
        "  {3}, {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5},",
        "  \"hrExport\" -> hrFile",
        "];",
        "Take[Import[hrFile, \"Lines\"], 12]"
      }, "\n"]
    ]
  },
  {
    evaluatedTestingExample[
      "\"hrExport\" writes the HR data to the specified file. With None, hop returns the text directly in the notebook:",
      StringRiffle[{
        "FileNameTake[hop[",
        "  {1}, {e1 -> 0.05},",
        "  \"hrExport\" -> $TemporaryDirectory",
        "]]"
      }, "\n"]
    ],
    evaluatedTestingExample[
      StringJoin[
        "When converting a supplied Hamiltonian, set \"wcc\" to the fractional center of each orbital. ",
        "This example uses the graphene centers; for a periodic-gauge Hamiltonian, replace them with one zero vector per orbital:"
      ],
      StringRiffle[{
        "grapheneNearest = symham[2];",
        "readHR[",
        "  hop[",
        "    grapheneNearest, {t1 -> 0.5},",
        "    \"hrExport\" -> FileNameJoin[",
        "      {$TemporaryDirectory, \"hop-explicit-wcc_hr.dat\"}",
        "    ],",
        "    \"wcc\" -> {",
        "      {1/3, 2/3, 0}, {2/3, 1/3, 0}",
        "    }",
        "  ]",
        "][\"Translations\"]"
      }, "\n"]
    ]
  },
  {
    {"readHR", "MagneticTB/ref/readHR"},
    {"showHoppingParameters", "MagneticTB/ref/showHoppingParameters"},
    {"symham", "MagneticTB/ref/symham"}
  },
    {"Wannier90", "hr.dat", "real-space hopping", "export"}
  ],
  <|"AllowPartialOptionExamples" -> True|>
];

(* -------------------------------------------------------------------------
   VASP import, interactive tuning, numerical fitting, and final comparison.

   These four pages deliberately form one workflow.  They share the real
   Examples/EIGENVAL sample and a three-orbital Hamiltonian generated through
   msgop -> init -> symham.  No hand-written Hamiltonian is used on the VASP
   pages, and every saved output comes from the current testing runtime.
   ------------------------------------------------------------------------- *)

testingVaspBandsCode = StringRiffle[{
  "eigenvalFile = FileNameJoin[{NotebookDirectory[], \"EIGENVAL\"}];",
  "vaspBands = vaspEig[eigenvalFile, -0.0072, 1, 21, 23];"
}, "\n"];

testingFittingModelCode = StringRiffle[{
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
  "fittingHamiltonian = symham[1] + symham[2];"
}, "\n"];

testingFittingPathCode = StringRiffle[{
  "fittingPath = {",
  "  {{{0, 0, 0}, {1/2, 0, 0}}, {\"\\[CapitalGamma]\", \"M\"}},",
  "  {{{1/2, 0, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
  "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"\\[CapitalGamma]\"}}",
  "};"
}, "\n"];

testingFittingInitialRulesCode = StringRiffle[{
  "initialRules = {",
  "  e1 -> 0.915, e2 -> 0.685,",
  "  t1 -> 0.1, t2 -> 0.09, t3 -> 0.1,",
  "  t4 -> -0.38, t5 -> 0.34, t6 -> 0.245,",
  "  t7 -> 0.565, t8 -> 0.305, t9 -> 0.365",
  "};"
}, "\n"];

testingVaspEigFixture = If[
  testingDocumentationTargetedQ[{
    "vaspEig", "bandManipulateEig", "compareBand", "fittingTB",
    "BandFitting"
  }],
  checkedDocumentationEvaluation[
    "TestingVaspEig",
    Module[{bands},
      bands = MagneticTB`vaspEig[
        FileNameJoin[{projectRoot, "Examples", "EIGENVAL"}],
        -0.0072,
        1,
        21,
        23
      ];
      <|
        "Bands" -> bands,
        "First" -> First[bands],
        "Dimensions" -> Dimensions[bands],
        "Plot" -> ListLinePlot[
          Transpose[bands[[All, 2]]],
          Frame -> True,
          FrameLabel -> {"k-point index", "Energy relative to efermi"},
          PlotRange -> All
        ]
      |>
    ]
  ],
  <|
    "Bands" -> Missing["NotEvaluatedForIncrementalPage"],
    "First" -> Missing["NotEvaluatedForIncrementalPage"],
    "Dimensions" -> Missing["NotEvaluatedForIncrementalPage"],
    "Plot" -> Missing["NotEvaluatedForIncrementalPage"]
  |>
];

testingFittingModelFixture = If[
  testingDocumentationTargetedQ[{
    "bandManipulateEig", "compareBand", "fittingTB", "BandFitting"
  }],
  checkedDocumentationEvaluation[
    "TestingFittingModel",
    Module[{operations, hamiltonian, path, initialRules},
      operations = Take[
        MagneticTB`msgop[MagneticTB`gray[156]],
        6
      ];
      Block[{Print = Function[Null]},
        MagneticTB`init[
          MagneticTB`lattice -> {
            {MagneticTB`a, 0, 0},
            {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
            {0, 0, MagneticTB`c}
          },
          MagneticTB`lattpar -> {
            MagneticTB`a -> 1,
            MagneticTB`c -> 10
          },
          MagneticTB`wyckoffposition -> {
            {{2/3, 1/3, 0}, {0, 0, 0}}
          },
          MagneticTB`symminformation -> operations,
          MagneticTB`basisFunctions -> {{"dz2", "dxy", "dx2-y2"}},
          MagneticTB`InitialBondShells -> 3
        ];
        hamiltonian = MagneticTB`symham[1] + MagneticTB`symham[2]
      ];
      path = {
        {{{0, 0, 0}, {1/2, 0, 0}}, {"\[CapitalGamma]", "M"}},
        {{{1/2, 0, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
        {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "\[CapitalGamma]"}}
      };
      initialRules = {
        Global`e1 -> 0.915,
        Global`e2 -> 0.685,
        Global`t1 -> 0.1,
        Global`t2 -> 0.09,
        Global`t3 -> 0.1,
        Global`t4 -> -0.38,
        Global`t5 -> 0.34,
        Global`t6 -> 0.245,
        Global`t7 -> 0.565,
        Global`t8 -> 0.305,
        Global`t9 -> 0.365
      };
      <|
        "Hamiltonian" -> hamiltonian,
        "Path" -> path,
        "InitialRules" -> initialRules
      |>
    ]
  ],
  <|
    "Hamiltonian" -> Missing["NotEvaluatedForIncrementalPage"],
    "Path" -> Missing["NotEvaluatedForIncrementalPage"],
    "InitialRules" -> Missing["NotEvaluatedForIncrementalPage"]
  |>
];

testingBandManipulateEigFixture = If[
  testingDocumentationTargetedQ[{"bandManipulateEig"}],
  checkedDocumentationEvaluation[
    "TestingBandManipulateEig",
    MagneticTB`bandManipulateEig[
      testingFittingModelFixture["Hamiltonian"],
      testingVaspEigFixture["Bands"]
    ]
  ],
  Missing["NotEvaluatedForIncrementalPage", "bandManipulateEig"]
];

testingVaspFittingFixture = If[
  testingDocumentationTargetedQ[{
    "compareBand", "fittingTB", "BandFitting"
  }],
  checkedDocumentationEvaluation[
    "TestingVaspFitting",
    Module[
      {
        nearestHamiltonian, nearestInitialRules, nearestFit,
        secondHamiltonian, secondInitialRules, secondFit
      },
      nearestHamiltonian = testingFittingModelFixture["Hamiltonian"];
      nearestInitialRules = testingFittingModelFixture["InitialRules"];
      nearestFit = MagneticTB`fittingTB[
        nearestHamiltonian,
        testingVaspEigFixture["Bands"],
        Range[Length[testingVaspEigFixture["Bands"]]],
        nearestInitialRules
      ];
      secondHamiltonian = nearestHamiltonian + MagneticTB`symham[3];
      secondInitialRules = Join[
        nearestFit["FittedParams"],
        Thread[{
          Global`r1, Global`r2, Global`r3, Global`r4, Global`r5,
          Global`r6, Global`r7, Global`r8, Global`r9, Global`r10
        } -> 0.]
      ];
      secondFit = MagneticTB`fittingTB[
        secondHamiltonian,
        testingVaspEigFixture["Bands"],
        Range[Length[testingVaspEigFixture["Bands"]]],
        secondInitialRules
      ];
      <|
        "NearestHamiltonian" -> nearestHamiltonian,
        "NearestInitialRules" -> nearestInitialRules,
        "NearestFit" -> nearestFit,
        "NearestPlot" -> MagneticTB`compareBand[
          testingFittingModelFixture["Path"], 30,
          nearestHamiltonian, nearestFit["FittedParams"],
          testingVaspEigFixture["Bands"],
          MagneticTB`plotRange -> {-1.5, 3}
        ],
        "Hamiltonian" -> secondHamiltonian,
        "InitialRules" -> secondInitialRules,
        "Fit" -> secondFit,
        "Parameters" -> secondFit["FittedParams"],
        "SecondPlot" -> MagneticTB`compareBand[
          testingFittingModelFixture["Path"], 30,
          secondHamiltonian, secondFit["FittedParams"],
          testingVaspEigFixture["Bands"],
          MagneticTB`plotRange -> {-1.5, 3}
        ]
      |>
    ]
  ],
  <|
    "Hamiltonian" -> Missing["NotEvaluatedForIncrementalPage"],
    "InitialRules" -> Missing["NotEvaluatedForIncrementalPage"],
    "Fit" -> Missing["NotEvaluatedForIncrementalPage"],
    "Parameters" -> Missing["NotEvaluatedForIncrementalPage"],
    "NearestHamiltonian" -> Missing["NotEvaluatedForIncrementalPage"],
    "NearestInitialRules" -> Missing["NotEvaluatedForIncrementalPage"],
    "NearestFit" -> Missing["NotEvaluatedForIncrementalPage"],
    "NearestPlot" -> Missing["NotEvaluatedForIncrementalPage"],
    "SecondPlot" -> Missing["NotEvaluatedForIncrementalPage"]
  |>
];

(* The tutorial keeps two focused fits separate from the function-page option
   fixture.  This avoids evaluating every fittingTB option merely to obtain
   the two plots that explain EnergyWindow and KPointNeighborhood. *)
testingBandFittingTutorialFixture = If[
  testingDocumentationTargetedQ[{"BandFitting"}],
  checkedDocumentationEvaluation[
    "TestingBandFittingTutorial",
    Module[
      {
        hamiltonian, referenceBands, indices, initialRules,
        energyWindowFit, kFit
      },
      hamiltonian = testingVaspFittingFixture["NearestHamiltonian"];
      referenceBands = testingVaspEigFixture["Bands"];
      indices = Range[Length[referenceBands]];
      initialRules = testingVaspFittingFixture["NearestInitialRules"];
      energyWindowFit = MagneticTB`fittingTB[
        hamiltonian,
        referenceBands,
        indices,
        initialRules,
        "EnergyWindow" -> {-0.5, 0.5}
      ];
      kFit = MagneticTB`fittingTB[
        hamiltonian,
        referenceBands,
        indices,
        initialRules,
        "KPointNeighborhood" -> <|
          "Center" -> {1/3, 1/3, 0},
          "Radius" -> 0.15,
          "Periodic" -> True
        |>
      ];
      <|
        "EnergyWindowPlot" -> MagneticTB`compareBand[
          testingFittingModelFixture["Path"],
          30,
          hamiltonian,
          energyWindowFit["FittedParams"],
          referenceBands,
          MagneticTB`plotRange -> {-1.5, 3}
        ],
        "KPointPlot" -> kFit["ComparisonPlot"]
      |>
    ]
  ],
  <|
    "EnergyWindowPlot" -> Missing["NotEvaluatedForIncrementalPage"],
    "KPointPlot" -> Missing["NotEvaluatedForIncrementalPage"]
  |>
];

testingFittingTBFixture = If[
  testingDocumentationTargetedQ[{"fittingTB"}],
  checkedDocumentationEvaluation[
    "TestingFittingTB",
    Module[
      {hamiltonian, referenceBands, initialRules, indices, fit, weights},
      hamiltonian = testingVaspFittingFixture["NearestHamiltonian"];
      referenceBands = testingVaspEigFixture["Bands"];
      initialRules = testingVaspFittingFixture["NearestInitialRules"];
      indices = Range[Length[referenceBands]];
      fit = testingVaspFittingFixture["NearestFit"];
      weights = ConstantArray[1., Dimensions[referenceBands[[All, 2]]]];
      weights[[1 ;; 30, All]] = 5.;
      <|
        "MaxIterations" -> KeyTake[
          MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            MaxIterations -> 1
          ],
          {"Converged", "Iterations"}
        ],
        "TimeConstraint" -> MagneticTB`fittingTB[
          hamiltonian, referenceBands, indices, initialRules,
          TimeConstraint -> Infinity
        ]["Converged"],
        "FiniteDifferenceStep" -> KeyTake[
          MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "FiniteDifferenceStep" -> 10^-4
          ],
          {"Iterations", "LSQForFittedParams"}
        ],
        "InitialDamping" -> KeyTake[
          MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "InitialDamping" -> 10^-2
          ],
          {"Iterations", "FinalDamping"}
        ],
        "FitTolerance" -> KeyTake[
          MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "FitTolerance" -> 10^-6
          ],
          {"Converged", "Iterations"}
        ],
        "EnergyWindow" -> MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "EnergyWindow" -> {-0.5, 0.5}
          ]["ComparisonPlot"],
        "KPointNeighborhood" -> MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "KPointNeighborhood" -> <|
              "Center" -> {0, 0, 0},
              "Radius" -> 0.35,
              "Periodic" -> True
            |>
          ]["ComparisonPlot"],
        "BandSelection" -> MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "BandSelection" -> 1
          ]["ComparisonPlot"],
        "ResidualWeights" -> MagneticTB`fittingTB[
            hamiltonian, referenceBands, indices, initialRules,
            "ResidualWeights" -> weights
          ]["ComparisonPlot"]
      |>
    ]
  ],
  <|
    "MaxIterations" -> Missing["NotEvaluatedForIncrementalPage"],
    "TimeConstraint" -> Missing["NotEvaluatedForIncrementalPage"],
    "FiniteDifferenceStep" -> Missing["NotEvaluatedForIncrementalPage"],
    "InitialDamping" -> Missing["NotEvaluatedForIncrementalPage"],
    "FitTolerance" -> Missing["NotEvaluatedForIncrementalPage"],
    "EnergyWindow" -> Missing["NotEvaluatedForIncrementalPage"],
    "KPointNeighborhood" -> Missing["NotEvaluatedForIncrementalPage"],
    "BandSelection" -> Missing["NotEvaluatedForIncrementalPage"],
    "ResidualWeights" -> Missing["NotEvaluatedForIncrementalPage"]
  |>
];

testingCompareBandFixture = If[
  testingDocumentationTargetedQ[{"compareBand"}],
  checkedDocumentationEvaluation[
    "TestingCompareBand",
    <|
      "Default" -> MagneticTB`compareBand[
        testingFittingModelFixture["Path"],
        30,
        testingVaspFittingFixture["Hamiltonian"],
        testingVaspFittingFixture["Parameters"],
        testingVaspEigFixture["Bands"],
        MagneticTB`plotRange -> {-1.5, 3}
      ],
      "PlotRange" -> MagneticTB`compareBand[
        testingFittingModelFixture["Path"],
        30,
        testingVaspFittingFixture["Hamiltonian"],
        testingVaspFittingFixture["Parameters"],
        testingVaspEigFixture["Bands"],
        MagneticTB`plotRange -> {-0.5, 0.5}
      ]
    |>
  ],
  <|
    "Default" -> Missing["NotEvaluatedForIncrementalPage"],
    "PlotRange" -> Missing["NotEvaluatedForIncrementalPage"]
  |>
];

testingCurrentPage = "vaspEig";
testingVaspEigSpec = Join[
  testingFunctionSpec[
    "vaspEig",
    {{
      "vaspEig[file,efermi,spin,startBand,endBand]",
      "reads a selected spin channel and inclusive band interval from a complete numeric VASP EIGENVAL file."
    }},
    {
      "The result contains one {{kx,ky,kz},{energy1,energy2,...}} record per k point. The reciprocal coordinates remain fractional and every energy is shifted by subtracting efermi.",
      "Band numbers are 1 based. startBand and endBand are both included. spin -> 1 selects the first energy column; spin -> 2 is available only in a two-spin EIGENVAL file.",
      "The parser accepts only the complete numeric VASP layout. It does not evaluate file contents or silently keep a partial record.",
      testingFailureNote[
        "vaspEig",
        "the file is unreadable or incomplete, a numeric record is malformed, spin is unavailable, or the inclusive band interval is invalid"
      ]
    },
    {},
    {
      {
        "Read bands 21 through 23 from the included EIGENVAL file and show the first k-point record:",
        testingVaspBandsCode <> "\nFirst[vaspBands]",
        testingVaspEigFixture["First"]
      },
      {
        "Plot the three imported bands before fitting:",
        StringRiffle[{
          "ListLinePlot[",
          "  Transpose[vaspBands[[All, 2]]],",
          "  Frame -> True,",
          "  FrameLabel -> {\"k-point index\", \"Energy relative to efermi\"},",
          "  PlotRange -> All",
          "]"
        }, "\n"],
        testingVaspEigFixture["Plot"]
      }
    },
    {},
    {
      {"bandManipulateEig", "MagneticTB/ref/bandManipulateEig"},
      {"fittingTB", "MagneticTB/ref/fittingTB"},
      {"compareBand", "MagneticTB/ref/compareBand"}
    },
    {"VASP", "EIGENVAL", "band import", "Fermi energy"},
    {{
      "Fitting MagneticTB bands to VASP",
      "MagneticTB/tutorial/BandFitting"
    }}
  ],
  <|
    "Properties" -> {{
      "The included sample contains 90 k points and three selected bands:",
      "Dimensions[vaspBands]",
      testingVaspEigFixture["Dimensions"]
    }}
  |>
];

testingCurrentPage = "bandManipulateEig";
testingBandManipulateEigSpec = testingFunctionSpec[
  "bandManipulateEig",
  {{
    "bandManipulateEig[H,referenceBands]",
    "opens a Manipulate panel for tuning the symbolic parameters of H against numeric reference bands."
  }},
  {
    "referenceBands uses the {{kx,ky,kz},{energy1,energy2,...}} record format returned by vaspEig. H is evaluated at 2 Pi times each reciprocal-fractional k point.",
    "Every ordinary symbolic parameter in H receives a slider from -1 to 1 with initial value 0. The ExportData button prints the current parameter rules; it does not write a file.",
    "Tight-binding bands are black and reference bands are red. The displayed energy range is the reference range expanded by 0.2 at each end.",
    "bandManipulateEig has no public options and returns a Manipulate expression.",
    "A malformed Hamiltonian, malformed reference data, Subscript or Indexed parameters, or unresolved function heads return a descriptive Failure; the function does not substitute a numerical fallback."
  },
  {},
  {
    {
      "Read the included VASP bands and construct the symmetry-generated three-orbital hexagonal Hamiltonian:",
      testingVaspBandsCode <> "\n\n" <> testingFittingModelCode
    },
    {
      "Move the sliders until the black tight-binding bands follow the red VASP bands:",
      "bandManipulateEig[fittingHamiltonian, vaspBands]",
      portableBandManipulateOutput[testingBandManipulateEigFixture]
    }
  },
  {},
  {
    {"vaspEig", "MagneticTB/ref/vaspEig"},
    {"fittingTB", "MagneticTB/ref/fittingTB"},
    {"compareBand", "MagneticTB/ref/compareBand"}
  },
  {"interactive bands", "VASP", "EIGENVAL", "parameter fitting"},
  {{
    "Fitting MagneticTB bands to VASP",
    "MagneticTB/tutorial/BandFitting"
  }}
];

testingCurrentPage = "fittingTB";
testingFittingTBSpec = testingFunctionSpec[
  "fittingTB",
  {{
    "fittingTB[H,referenceBands,kRange,initialRules]",
    "fits every affine-linear symbolic parameter of H to selected numeric reference bands."
  }},
  {
    "referenceBands uses the {{kx,ky,kz},{energy1,energy2,...}} format returned by vaspEig. kRange is a nonempty list of 1-based record indices, and every selected record must contain exactly Dimensions[H][[1]] energies.",
    "The EnergyWindow, KPointNeighborhood, BandSelection, and ResidualWeights selectors are intersected with kRange. EnergyWindow masks individual reference energies, not complete k-point records.",
    "KPointNeighborhood uses fractional reciprocal coordinates. With Periodic -> True, the displacement is reduced by the componentwise minimum-image convention before applying a spherical Radius or componentwise Range.",
    "BandSelection uses 1-based reference positions. At every retained k point, sorted model eigenvalue j is paired with reference position j; the reference bands are not dynamically reordered.",
    "ResidualWeights is Automatic or a nonnegative record-by-band matrix with the same shape as referenceBands[[All,2]]. A zero weight removes that residual; positive weights minimize the weighted sum of squared band differences.",
    "initialRules can be a rule list or Association and must assign a real numeric starting value to every fitting parameter. H is evaluated at 2 Pi times each reciprocal-fractional k point.",
    "The Hamiltonian must be affine-linear in the fitting parameters. MagneticTB caches its constant and parameter-coefficient matrices, then uses a finite-difference Levenberg-Marquardt optimizer; it does not switch to a slower nonlinear symbolic fallback.",
    "The objective minimizes sorted energy residuals only. It does not preserve a band gap, band character, or a chosen band ordering; always inspect ComparisonPlot and reject a lower-loss fit when its band structure is physically wrong.",
    StringJoin[
      "The result reports FittedParams, convergence diagnostics, ComparisonPlot, CandidateKPointCount, ",
      "UsedKPointCount, UsedKPointIndices, ResidualCount, the effective selectors, weighting, and band pairing."
    ],
    testingFailureNote[
      "fittingTB",
      "the Hamiltonian or reference records are invalid, kRange or a selector is invalid, the combined selection is empty or underdetermined, initial values are incomplete, the model is nonlinear in its parameters, the time limit is exceeded, or the optimizer fails"
    ]
  },
  {
    {"MaxIterations", "100", "Positive integer upper bound on Levenberg-Marquardt iterations."},
    {"TimeConstraint", "60", "Positive fitting time limit in seconds, or Infinity."},
    {"\"EnergyWindow\"", "All", "All keeps every reference energy; {emin,emax} keeps energies in the inclusive interval."},
    {"\"KPointNeighborhood\"", "All", "All keeps every k point; an Association selects a periodic or nonperiodic Radius or componentwise Range around Center."},
    {"\"BandSelection\"", "All", "All keeps every reference position; an integer, increasing integer list, or positive-step Span selects 1-based band positions."},
    {"\"ResidualWeights\"", "Automatic", "Automatic assigns unit weight; a nonnegative record-by-band matrix weights individual residuals, with zero removing a residual."},
    {"\"FiniteDifferenceStep\"", "10^-5", "Positive relative step used for central finite-difference Jacobian columns."},
    {"\"InitialDamping\"", "10^-3", "Positive initial Levenberg-Marquardt damping factor."},
    {"\"FitTolerance\"", "10^-8", "Positive stopping tolerance for the gradient, step size, and objective improvement."}
  },
  {
    {
      "Read the VASP three-band data and construct the symmetry-generated hexagonal model:",
      testingVaspBandsCode <> "\n\n" <> testingFittingModelCode <>
        "\n" <> testingFittingPathCode <>
        "\n" <> testingFittingInitialRulesCode
    },
    {
      "Fit all three VASP bands with onsite and nearest-neighbour hopping, then draw the full-range comparison:",
      StringRiffle[{
        "nearestFit = fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules",
        "];",
        "compareBand[",
        "  fittingPath, 30, fittingHamiltonian,",
        "  nearestFit[\"FittedParams\"], vaspBands,",
        "  plotRange -> {-1.5, 3}",
        "]"
      }, "\n"],
      testingVaspFittingFixture["NearestPlot"]
    },
    {
      "Add second-neighbour hopping, refit the same three VASP bands, and draw the full-range comparison again:",
      StringRiffle[{
        "secondNeighborHamiltonian = fittingHamiltonian + symham[3];",
        "secondNeighborInitialRules = Join[",
        "  nearestFit[\"FittedParams\"],",
        "  Thread[{r1, r2, r3, r4, r5, r6, r7, r8, r9, r10} -> 0.]",
        "];",
        "secondNeighborFit = fittingTB[",
        "  secondNeighborHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], secondNeighborInitialRules",
        "];",
        "compareBand[",
        "  fittingPath, 30, secondNeighborHamiltonian,",
        "  secondNeighborFit[\"FittedParams\"], vaspBands,",
        "  plotRange -> {-1.5, 3}",
        "]"
      }, "\n"],
      testingVaspFittingFixture["SecondPlot"]
    }
  },
  {
    {
      "\"EnergyWindow\" - fit only reference energies from -0.5 through 0.5:",
      StringRiffle[{
        "fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules,",
        "  \"EnergyWindow\" -> {-0.5, 0.5}",
        "][\"ComparisonPlot\"]"
      }, "\n"],
      testingFittingTBFixture["EnergyWindow"]
    },
    {
      "\"KPointNeighborhood\" - fit the periodic radius 0.35 around Gamma in fractional reciprocal coordinates:",
      StringRiffle[{
        "fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules,",
        "  \"KPointNeighborhood\" -> <|",
        "    \"Center\" -> {0, 0, 0},",
        "    \"Radius\" -> 0.35,",
        "    \"Periodic\" -> True",
        "  |>",
        "][\"ComparisonPlot\"]"
      }, "\n"],
      testingFittingTBFixture["KPointNeighborhood"]
    },
    {
      "\"BandSelection\" - select the first sorted reference position:",
      StringRiffle[{
        "fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules,",
        "  \"BandSelection\" -> 1",
        "][\"ComparisonPlot\"]"
      }, "\n"],
      testingFittingTBFixture["BandSelection"]
    },
    {
      "\"ResidualWeights\" - give the first path segment five times the weight of the remaining points:",
      StringRiffle[{
        "fitWeights = ConstantArray[1., Dimensions[vaspBands[[All, 2]]]];",
        "fitWeights[[1 ;; 30, All]] = 5.;",
        "fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules,",
        "  \"ResidualWeights\" -> fitWeights",
        "][\"ComparisonPlot\"]"
      }, "\n"],
      testingFittingTBFixture["ResidualWeights"]
    },
    {
      "MaxIterations - stop this example after one optimizer iteration:",
      StringRiffle[{
        "KeyTake[",
        "  fittingTB[",
        "    fittingHamiltonian, vaspBands,",
        "    Range[Length[vaspBands]], initialRules,",
        "    MaxIterations -> 1",
        "  ],",
        "  {\"Converged\", \"Iterations\"}",
        "]"
      }, "\n"],
      testingFittingTBFixture["MaxIterations"]
    },
    {
      "TimeConstraint - remove the elapsed-time limit while retaining the iteration limit:",
      StringRiffle[{
        "fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules,",
        "  TimeConstraint -> Infinity",
        "][\"Converged\"]"
      }, "\n"],
      testingFittingTBFixture["TimeConstraint"]
    },
    {
      "\"FiniteDifferenceStep\" - use a larger relative step for the central-difference Jacobian:",
      StringRiffle[{
        "KeyTake[",
        "  fittingTB[",
        "    fittingHamiltonian, vaspBands,",
        "    Range[Length[vaspBands]], initialRules,",
        "    \"FiniteDifferenceStep\" -> 10^-4",
        "  ],",
        "  {\"Iterations\", \"LSQForFittedParams\"}",
        "]"
      }, "\n"],
      testingFittingTBFixture["FiniteDifferenceStep"]
    },
    {
      "\"InitialDamping\" - start from a stronger damping factor:",
      StringRiffle[{
        "KeyTake[",
        "  fittingTB[",
        "    fittingHamiltonian, vaspBands,",
        "    Range[Length[vaspBands]], initialRules,",
        "    \"InitialDamping\" -> 10^-2",
        "  ],",
        "  {\"Iterations\", \"FinalDamping\"}",
        "]"
      }, "\n"],
      testingFittingTBFixture["InitialDamping"]
    },
    {
      "\"FitTolerance\" - use a looser stopping tolerance:",
      StringRiffle[{
        "KeyTake[",
        "  fittingTB[",
        "    fittingHamiltonian, vaspBands,",
        "    Range[Length[vaspBands]], initialRules,",
        "    \"FitTolerance\" -> 10^-6",
        "  ],",
        "  {\"Converged\", \"Iterations\"}",
        "]"
      }, "\n"],
      testingFittingTBFixture["FitTolerance"]
    }
  },
  {
    {"vaspEig", "MagneticTB/ref/vaspEig"},
    {"bandManipulateEig", "MagneticTB/ref/bandManipulateEig"},
    {"compareBand", "MagneticTB/ref/compareBand"}
  },
  {"band fitting", "Levenberg-Marquardt", "VASP", "least squares"},
  {{
    "Fitting MagneticTB bands to VASP",
    "MagneticTB/tutorial/BandFitting"
  }}
];

testingCurrentPage = "compareBand";
testingCompareBandSpec = testingFunctionSpec[
  "compareBand",
  {{
    "compareBand[path,npoint,H,rules,referenceBands]",
    "overlays numerical tight-binding bands with supplied reference bands along a labeled path."
  }},
  {
    "path is a nonempty list of numeric three-component reciprocal-fractional segment endpoints and endpoint labels. npoint is a positive integer.",
    "referenceBands uses the vaspEig record format and must contain exactly Length[path] npoint records in the same segment order. The tight-binding curve includes both endpoints of every segment.",
    "rules can be a rule list or Association. After applying rules and 2 Pi times each path momentum, H must be a numerical square matrix.",
    "The result is a Labeled plot: tight-binding bands are purple and reference bands are blue. plotRange changes only the visible energy interval.",
    testingFailureNote[
      "compareBand",
      "H, path, npoint, rules, or referenceBands are malformed, their record counts disagree, or H remains nonnumeric"
    ]
  },
  {{
    "plotRange",
    "All",
    "All displays the complete overlay; {emin,emax} restricts only the visible energy interval."
  }},
  {
    {
      "Read the VASP bands and construct the symmetry-generated model:",
      testingVaspBandsCode <> "\n\n" <> testingFittingModelCode <>
        "\n" <> testingFittingPathCode <>
        "\n" <> testingFittingInitialRulesCode
    },
    {
      "Fit the nearest-neighbour model, add second-neighbour hopping, and refit all 90 records:",
      StringRiffle[{
        "nearestFit = fittingTB[",
        "  fittingHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], initialRules",
        "];",
        "secondNeighborHamiltonian = fittingHamiltonian + symham[3];",
        "secondNeighborInitialRules = Join[",
        "  nearestFit[\"FittedParams\"],",
        "  Thread[{r1, r2, r3, r4, r5, r6, r7, r8, r9, r10} -> 0.]",
        "];",
        "fitResult = fittingTB[",
        "  secondNeighborHamiltonian, vaspBands,",
        "  Range[Length[vaspBands]], secondNeighborInitialRules",
        "];"
      }, "\n"]
    },
    {
      "Display all three selected VASP bands together with the fitted second-neighbour model from -1.5 to 3 eV:",
      StringRiffle[{
        "compareBand[",
        "  fittingPath, 30, secondNeighborHamiltonian,",
        "  fitResult[\"FittedParams\"], vaspBands,",
        "  plotRange -> {-1.5, 3}",
        "]"
      }, "\n"],
      testingCompareBandFixture["Default"]
    }
  },
  {{
    "plotRange - zoom into the interval from -0.5 to 0.5 after inspecting all three bands:",
    StringRiffle[{
      "compareBand[",
      "  fittingPath, 30, secondNeighborHamiltonian,",
      "  fitResult[\"FittedParams\"], vaspBands,",
      "  plotRange -> {-0.5, 0.5}",
      "]"
    }, "\n"],
    testingCompareBandFixture["PlotRange"]
  }},
  {
    {"vaspEig", "MagneticTB/ref/vaspEig"},
    {"bandManipulateEig", "MagneticTB/ref/bandManipulateEig"},
    {"fittingTB", "MagneticTB/ref/fittingTB"},
    {"bandplot", "MagneticTB/ref/bandplot"}
  },
  {"band comparison", "VASP", "EIGENVAL", "tight binding"},
  {{
    "Fitting MagneticTB bands to VASP",
    "MagneticTB/tutorial/BandFitting"
  }}
];

(* -------------------------------------------------------------------------
   A complete VASP fitting workflow.

   Keep the reference pages compact: the tutorial is the one place where the
   nearest-neighbor model, its second-neighbor extension, and two focused-fit
   selectors are compared through the actual three-band plots.
   ------------------------------------------------------------------------- *)

testingBandFittingTutorialPage = Notebook[{
  Cell["Fitting MagneticTB bands to VASP", "Title"],
  Cell[
    StringJoin[
      "To fit a tight-binding model to VASP bands, first read EIGENVAL and construct the ",
      "Hamiltonian from symmetry. Here we fit three bands with a three-orbital model, compare ",
      "nearest-neighbor and second-neighbor results, and then fit a selected energy or k-point ",
      "range. Plot the fitted and VASP bands together to see both the agreement and any ",
      "remaining differences."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["Needs[\"MagneticTB`\"]"]], "Input"],

  Cell["Read and inspect the three VASP bands", "Section"],
  Cell[
    StringJoin[
      "Use vaspEig to read bands 21-23 and subtract the Fermi energy. First plot all three ",
      "bands to see the reference band structure:"
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "eigenvalFile = FileNameJoin[{",
      "  NotebookDirectory[],",
      "  \"..\", \"ReferencePages\", \"Symbols\", \"EIGENVAL\"",
      "}];",
      "vaspBands = vaspEig[eigenvalFile, -0.0072, 1, 21, 23];",
      "ListLinePlot[",
      "  Transpose[vaspBands[[All, 2]]],",
      "  Frame -> True,",
      "  FrameLabel -> {\"k-point index\", \"Energy relative to efermi\"},",
      "  PlotRange -> All",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[BoxData[outputBoxes[testingVaspEigFixture["Plot"]]], "Output"],

  Cell["Build the symmetry-generated model", "Section"],
  Cell[
    StringJoin[
      "For this three-orbital model, take the first six operations from msgop[gray[156]], ",
      "which form its complete unitary subgroup. The onsite and nearest-neighbor terms contain ",
      "e1, e2, and t1-t9. These are the parameters to be fitted:"
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      testingFittingModelCode,
      testingFittingPathCode,
      testingFittingInitialRulesCode
    }, "\n"]]],
    "Input"
  ],

  Cell["Fit the onsite and nearest-neighbor model", "Section"],
  Cell[
    "Fit the three bands at all 90 k points. Then draw them together with the VASP bands:",
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "nearestFit = fittingTB[",
      "  fittingHamiltonian, vaspBands,",
      "  Range[Length[vaspBands]], initialRules",
      "];",
      "compareBand[",
      "  fittingPath, 30, fittingHamiltonian,",
      "  nearestFit[\"FittedParams\"], vaspBands,",
      "  plotRange -> {-1.5, 3}",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[BoxData[outputBoxes[testingVaspFittingFixture["NearestPlot"]]], "Output"],

  Cell["Add second-neighbor hopping", "Section"],
  Cell[
    StringJoin[
      "To include second-neighbor hopping, keep the fitted nearest-neighbor parameters as the ",
      "initial values and set the new parameters r1-r10 to zero. Fit again and compare the two ",
      "band plots to see where the extra hoppings improve the agreement."
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "secondNeighborHamiltonian = fittingHamiltonian + symham[3];",
      "secondNeighborInitialRules = Join[",
      "  nearestFit[\"FittedParams\"],",
      "  Thread[{r1, r2, r3, r4, r5, r6, r7, r8, r9, r10} -> 0.]",
      "];",
      "secondNeighborFit = fittingTB[",
      "  secondNeighborHamiltonian, vaspBands,",
      "  Range[Length[vaspBands]], secondNeighborInitialRules",
      "];",
      "compareBand[",
      "  fittingPath, 30, secondNeighborHamiltonian,",
      "  secondNeighborFit[\"FittedParams\"], vaspBands,",
      "  plotRange -> {-1.5, 3}",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[BoxData[outputBoxes[testingVaspFittingFixture["SecondPlot"]]], "Output"],

  Cell["Fit the energy range of interest", "Section"],
  Cell[
    StringJoin[
      "If only bands near the Fermi energy are of interest, set EnergyWindow to the desired ",
      "energy interval. Only reference energies in this interval enter the fit. The plot still ",
      "shows the full bands, so changes outside the interval can also be seen."
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "fermiFit = fittingTB[",
      "  fittingHamiltonian, vaspBands,",
      "  Range[Length[vaspBands]], initialRules,",
      "  \"EnergyWindow\" -> {-0.5, 0.5}",
      "];",
      "compareBand[",
      "  fittingPath, 30, fittingHamiltonian,",
      "  fermiFit[\"FittedParams\"], vaspBands,",
      "  plotRange -> {-1.5, 3}",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[
      testingBandFittingTutorialFixture["EnergyWindowPlot"]
    ]],
    "Output"
  ],

  Cell["Fit near the K point", "Section"],
  Cell[
    StringJoin[
      "To fit near K, set KPointNeighborhood to K and radius 0.15 in reciprocal fractional ",
      "coordinates. Distances include periodic equivalents of k points. ComparisonPlot ",
      "displays only this neighborhood: blue points are VASP, gray dashed curves are the ",
      "initial model, and red curves are the fitted model."
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "kFit = fittingTB[",
      "  fittingHamiltonian, vaspBands,",
      "  Range[Length[vaspBands]], initialRules,",
      "  \"KPointNeighborhood\" -> <|",
      "    \"Center\" -> {1/3, 1/3, 0},",
      "    \"Radius\" -> 0.15,",
      "    \"Periodic\" -> True",
      "  |>",
      "];",
      "kFit[\"ComparisonPlot\"]"
    }, "\n"]]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[testingBandFittingTutorialFixture["KPointPlot"]]],
    "Output"
  ],

  Cell["Judge the plotted bands, not only the loss", "Section"],
  Cell[
    StringJoin[
      "fittingTB compares sorted model eigenvalues with the corresponding reference energies. ",
      "It does not impose a gap or track orbital character and band connectivity. A small ",
      "energy error alone is therefore not enough: also compare gaps and crossings in the ",
      "plot. Use physically reasonable initial parameters, and retain the fitted values when ",
      "adding more neighbor hoppings."
    ],
    "Text"
  ],
  Cell["Related functions", "Section"],
  Cell[
    TextData[{
      documentationLink["vaspEig", "MagneticTB/ref/vaspEig"],
      "   ",
      documentationLink["bandManipulateEig", "MagneticTB/ref/bandManipulateEig"],
      "   ",
      documentationLink["fittingTB", "MagneticTB/ref/fittingTB"],
      "   ",
      documentationLink["compareBand", "MagneticTB/ref/compareBand"]
    }],
    "Text"
  ],
  historyCells[],
  categorizationCells[
    "Tech Note",
    "MagneticTB`",
    "MagneticTB/tutorial/BandFitting"
  ],
  keywordCells[{
    "VASP", "EIGENVAL", "band fitting", "nearest neighbor",
    "second neighbor", "EnergyWindow", "KPointNeighborhood"
  }]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Fitting MagneticTB bands to VASP",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"},
    "TechNotePageStylesExt.nb",
    CharacterEncoding -> "UTF-8"
  ]
];

(* -------------------------------------------------------------------------
   Exported specification list and shared tutorials.

   GenerateDocumentationSources.wls validates this list by name before it
   merges the pages with the stable documentation surface.
   ------------------------------------------------------------------------- *)

testingFunctionSpecifications = {
  testingBandplotSpec,
  testingShowbandSpec,
  testingStandardKPathSpec,
  testingShowCrystalStructureSpec,
  testingShowBrillouinZoneSpec,
  testingOrbitalTableSpec,
  testingHoppingDataSpec,
  testingTransformHoppingsSpec,
  testingBuildBlochHamiltonianSpec,
  testingBuildRealSpaceHamiltonianSpec,
  testingBuildSlabHamiltonianSpec,
  testingSurfaceGreenFunctionSpec,
  testingSurfaceSpectralFunctionSpec,
  testingWilsonLoopSpec,
  testingBerryPhaseSpec,
  testingBerryCurvatureSpec,
  testingPointChernNumberSpec,
  testingFindGaplessPointsSpec,
  testingGetMSGElemFromMSGCorepSpec,
  testingGetTBBandCorepSpec,
  testingGrayLayerSpec,
  testingGrayRodSpec,
  testingHopSpec,
  testingVaspEigSpec,
  testingBandManipulateEigSpec,
  testingFittingTBSpec,
  testingCompareBandSpec
};

(* A targeted reference-page build has everything it needs at this point.
   Avoid constructing or evaluating unrelated tutorial fixtures. *)
If[
  !AnyTrue[
    {
      Environment["MAGNETICTB_DOC_TESTING_PUBLIC_ONLY"],
      Environment["MAGNETICTB_DOC_BAND_FITTING_TUTORIAL_ONLY"]
    },
    MatchQ[#, 1 | "1" | True] &
  ],

(* Share the same small tutorial source with the targeted build. *)
Get[FileNameJoin[{scriptDirectory, "Examples", "CrystalAndKPathsDocumentation.wl"}]];

testingRealSpaceTutorialPage = Notebook[{
  Cell["Real-space Hamiltonians, surfaces, and topology", "Title"],
  Cell[
    StringJoin[
      "MagneticTB uses one explicit hopping convention from finite samples through Bloch, slab, and ",
      "semi-infinite calculations: H(R)[a,b]=<0,a|H|R,b> and H(k)=Sum_R H(R) Exp[2 Pi I k.R]."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["Needs[\"MagneticTB`\"]"]], "Input"],
  Cell["Symmetry-generated real-space Hamiltonian", "Section"],
  Cell[BoxData[codeBoxes[
    testingBuildRealSpaceHamiltonianSpec["Examples"][[1, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[testingBuildRealSpaceHamiltonianSpec["Examples"][[1, 3]]]], "Output"],
  Cell["Bloch bands from the same hopping data", "Section"],
  Cell[BoxData[codeBoxes[
    testingBuildBlochHamiltonianSpec["Examples"][[1, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    testingBuildBlochHamiltonianSpec["Examples"][[1, 3]]
  ]], "Output"],
  Cell["Topological-insulator surface spectrum", "Section"],
  Cell[BoxData[codeBoxes[
    testingSurfaceSpectralFunctionSpec["Examples"][[2, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    testingSurfaceSpectralFunctionSpec["Examples"][[2, 3]]
  ]], "Output"],
  Cell["Wilson and Berry geometry", "Section"],
  Cell[
    StringJoin[
      "Wilson loops and Berry phases use forward occupied-subspace links, SVD unitarization, and ",
      "orbital-center endpoint sewing. A closed path must be closed exactly modulo a reciprocal ",
      "vector."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    testingWilsonLoopSpec["Examples"][[3, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    testingWilsonLoopSpec["Examples"][[3, 3]]
  ]], "Output"],
  Cell["Point charge and gapless-point search", "Section"],
  Cell[
    StringJoin[
      "pointChernNumber integrates occupied-subspace flux on a closed cube. findGaplessPoints is a ",
      "bounded grid-and-refinement search and reports that completeness is not guaranteed."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    testingPointChernNumberSpec["Examples"][[2, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    testingPointChernNumberSpec["Examples"][[2, 3]]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    testingPointChernNumberSpec["Examples"][[3, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    testingPointChernNumberSpec["Examples"][[3, 3]]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    testingFindGaplessPointsSpec["Examples"][[2, 2]]
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    testingFindGaplessPointsSpec["Examples"][[2, 3]]
  ]], "Output"],
  historyCells[],
  categorizationCells["Tech Note", "MagneticTB`", "MagneticTB/tutorial/RealSpaceAndTopology"],
  keywordCells[{"real-space Hamiltonian", "surface Green function", "Wilson loop", "Chern number"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Real-space Hamiltonians, surfaces, and topology",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"},
    "TechNotePageStylesExt.nb",
    CharacterEncoding -> "UTF-8"
  ]
];

Get[FileNameJoin[{
  scriptDirectory, "Examples", "CorepresentationsDocumentation.wl"
}]];

];
