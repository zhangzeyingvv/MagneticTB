(* Generate editable authoring-notebook sources for MagneticTB 2.0.
   The generated .source.wl files are canonicalized with notebook_guard.py
   before being placed in Documentation/English. *)

ClearAll["Global`*"];

(* This tutorial explicitly loads MSGCorep before MagneticTB, matching its
   visible example and the program's existing corepresentation tests. *)
If[Environment["MAGNETICTB_DOC_PAGE"] === "Tutorials/Corepresentations",
  Needs["MSGCorep`"]
];

(* -------------------------------------------------------------------------
   Phase 0: resolve the authoring, expected-output, and installed-page trees.
   ------------------------------------------------------------------------- *)

(* Keep the generation path relative to the repository working directory.
   WolframScript 1.x can misdecode a non-ASCII absolute $InputFileName. *)
scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
projectRoot = DirectoryName[DirectoryName[scriptDirectory]];
documentationRoot = FileNameJoin[{projectRoot, "Documentation", "English"}];
sourceRoot = FileNameJoin[{scriptDirectory, "Sources"}];
expectedRoot = FileNameJoin[{scriptDirectory, "Expected"}];

(* These pages are intentionally outside the active help system while their
   examples are rewritten one at a time.  Their existing bilingual sources,
   expected text, and notebooks are preserved under DeferredPages/. *)
deferredDocumentationPageNames = {
  "buildSlabHamiltonian",
  "buildRealSpaceHamiltonian",
  "buildBlochHamiltonian",
  "plotBerryCurvature3D",
  "plotBerryCurvature2D",
  "plotRealSpaceWavefunction",
  "plotSurfaceSpectrum",
  "findGaplessPoints",
  "pointChernNumber",
  "plotWilsonLoop",
  "symmetrizationHRInit",
  "berryPhase",
  "berryCurvature",
  "surfaceGreenFunction",
  "surfaceSpectralFunction",
  "wilsonLoop"
};

deferredDocumentationTutorialNames = {
  "ValidatedModels",
  "WyckoffCrystalSystems",
  "RealSpaceAndTopology"
};

activeDocumentationSpecificationQ[spec_Association] :=
  !MemberQ[
    deferredDocumentationPageNames,
    Lookup[spec, "Name", ""]
  ];

Scan[
  If[!DirectoryQ[#], CreateDirectory[#, CreateIntermediateDirectories -> True]] &,
  {
    FileNameJoin[{documentationRoot, "ReferencePages", "Symbols"}],
    FileNameJoin[{documentationRoot, "Guides"}],
    FileNameJoin[{documentationRoot, "Tutorials"}],
    FileNameJoin[{sourceRoot, "ReferencePages", "Symbols"}],
    FileNameJoin[{sourceRoot, "Guides"}],
    FileNameJoin[{sourceRoot, "Tutorials"}],
    FileNameJoin[{expectedRoot, "ReferencePages", "Symbols"}],
    FileNameJoin[{expectedRoot, "Guides"}],
    FileNameJoin[{expectedRoot, "Tutorials"}]
  }
];

(* The VASP import, interactive tuning, fitting, and band-comparison pages use one
   real EIGENVAL file.  Keep the sample beside the generated reference pages
   so NotebookDirectory[] works both here and in an installed paclet. *)
compareBandEigenvalSource = FileNameJoin[{projectRoot, "Examples", "EIGENVAL"}];
compareBandEigenvalTarget = FileNameJoin[{
  documentationRoot,
  "ReferencePages",
  "Symbols",
  "EIGENVAL"
}];
If[
  FileExistsQ[compareBandEigenvalSource] &&
    (!FileExistsQ[compareBandEigenvalTarget] ||
      FileHash[compareBandEigenvalSource, "SHA256"] =!=
        FileHash[compareBandEigenvalTarget, "SHA256"]),
  CopyFile[
    compareBandEigenvalSource,
    compareBandEigenvalTarget,
    OverwriteTarget -> True
  ]
];

writeUTF8TextIfChanged[path_String, text_String] := Module[{existing},
  existing = If[
    FileExistsQ[path],
    FromCharacterCode[BinaryReadList[path], "UTF8"],
    Missing["NotFound"]
  ];
  If[existing =!= text,
    Export[path, text, "Text", CharacterEncoding -> "UTF-8"]
  ];
];

writeGuideSourceOnly[notebook_Notebook, expected_List] := Module[
  {sourcePath, expectedPath, sourceText, expectedText},
  sourcePath = FileNameJoin[{sourceRoot, "Guides", "MagneticTB.source.wl"}];
  expectedPath = FileNameJoin[{expectedRoot, "Guides", "MagneticTB.expected.txt"}];
  sourceText = ToString[notebook, InputForm, PageWidth -> Infinity];
  expectedText = StringRiffle[expected, "\n"] <> "\n";
  writeUTF8TextIfChanged[sourcePath, sourceText];
  writeUTF8TextIfChanged[expectedPath, expectedText];
  Print[sourcePath]
];

(* Keep executable source as a string inside BoxData.  Reconstructing boxes
   with ToExpression/ToBoxes discards the author's semantic line breaks and
   turns readable multi-line examples into one very wide row. *)
codeBoxes[code_String] := code;

outputBoxes[value_String] := codeBoxes[value];
outputBoxes[value_] := ToBoxes[value, StandardForm];

(* A saved Manipulate DynamicModule is tied to the kernel that created it and
   is not a portable documentation result.  Keep every executable Input cell,
   but omit only the immediately following saved dynamic Output cell. *)
ClearAll[
  bandManipulateInputCellQ,
  savedDynamicOutputCellQ,
  stripSavedBandManipulateOutputCells,
  stripSavedBandManipulateOutputs
];

bandManipulateInputCellQ[
    Cell[BoxData[code_String], "Input", ___]
  ] := AnyTrue[
  {"bandManipulate[", "bandManipulateEig[", "bandManipulateold["},
  StringContainsQ[code, #] &
];
bandManipulateInputCellQ[_] := False;

savedDynamicOutputCellQ[
    Cell[BoxData[boxes_], "Output", ___]
  ] := !FreeQ[Unevaluated[boxes], _DynamicModuleBox, Infinity];
savedDynamicOutputCellQ[_] := False;

stripSavedBandManipulateOutputCells[cells_List] := Module[
  {result = {}, index = 1, cell},
  While[index <= Length[cells],
    cell = cells[[index]] /. HoldPattern[
      Cell[CellGroupData[groupCells_List, state_], options___]
    ] :> Cell[
      CellGroupData[
        stripSavedBandManipulateOutputCells[groupCells],
        state
      ],
      options
    ];
    AppendTo[result, cell];
    If[
      bandManipulateInputCellQ[cell] &&
        index < Length[cells] &&
        savedDynamicOutputCellQ[cells[[index + 1]]],
      index += 2,
      index++
    ]
  ];
  result
];

stripSavedBandManipulateOutputs[
    Notebook[cells_List, options___]
  ] := Notebook[
  stripSavedBandManipulateOutputCells[cells],
  options
];

(* A few exact multi-orbital matrices are wider than the printable help-page
   body.  Keep the full physical result and reduce only its displayed size. *)
compactMatrixOutput[value_] := Style[
  value,
  Magnification -> 0.58
];

(* Store the exact boxes returned by Manipulate behind a short compressed
   payload.  The outer DynamicBox reconstructs them once when the help page is
   displayed, avoiding both a hand-drawn imitation and a giant raw box tree in
   the textual Notebook source. *)
portableBandManipulateOutput[panel_] := Module[
  {payload = Compress[ToBoxes[panel, StandardForm]]},
  RawBoxes@With[
    {encoded = payload},
    DynamicBox[
      Uncompress[encoded],
      SingleEvaluation -> True,
      TrackedSymbols :> {}
    ]
  ]
];

(* The expanded MoS2 and cubic 195.3-a matrices are substantially wider than
   the other tutorial matrices.  A separate print-only scale keeps every term
   inside an A4 help-page body without shrinking the ordinary outputs. *)
veryWideMatrixOutput[value_] := Style[
  value,
  Magnification -> 0.24
];

documentationLink[label_String, uri_String] := ButtonBox[
  label,
  BaseStyle -> "Link",
  ButtonData -> "paclet:" <> uri
];

(* A guide entry is always one paclet link followed by an em-dash summary.
   Keeping this shape in one helper makes the guide definition below read as
   content rather than repeated low-level Cell/TextData boilerplate. *)
guideEntryCell[
    label_String,
    uri_String,
    description_String
  ] := Cell[
  TextData[{
    documentationLink[label, uri],
    " \[LongDash] " <> description
  }],
  "GuideText"
];

(* Use the official guide tutorial styles for a separate, always-open
   tutorial list above the function index.  Keep each title on its own line. *)
guideTutorialCells[label_String, uri_String, description_String] := Sequence[
  Cell[TextData[{documentationLink[label, uri]}], "GuideTutorial"],
  Cell[description, "GuideText"]
];

(* Mathematica's guide stylesheet intentionally hides GuideFunctionsSection:
   it is the outer container, not a visible category heading.  Each visible
   category must therefore be a grouped GuideFunctionsSubsection. *)
guideSubsectionCell[title_String, entries_List] := Cell[CellGroupData[
  Prepend[entries, Cell[title, "GuideFunctionsSubsection"]],
  Open
]];

guideFunctionBlock[sections_List] := Cell[CellGroupData[
  Prepend[
    Riffle[sections, Cell["\t", "GuideDelimiter"]],
    Cell["", "GuideFunctionsSection"]
  ],
  Open
]];

historyCells[] := Cell[CellGroupData[{
  Cell["History", "HistorySection"],
  Cell[
    TextData[{
      "New in: ", Cell["2.0.0", "HistoryData", CellTags -> "New"],
      " | Modified in: ", Cell[" ", "HistoryData", CellTags -> "Modified"],
      " | Obsolete in: ", Cell[" ", "HistoryData", CellTags -> "Obsolete"]
    }],
    "History"
  ]
}, Closed]];

categorizationCells[entityType_String, context_String, uri_String] :=
  Cell[CellGroupData[{
    Cell["Categorization", "CategorizationSection"],
    Cell[entityType, "Categorization", CellLabel -> "Entity Type"],
    Cell["MagneticTB", "Categorization", CellLabel -> "Paclet Name"],
    Cell[context, "Categorization", CellLabel -> "Context"],
    Cell[uri, "Categorization", CellLabel -> "URI"]
  }, Closed]];

keywordCells[keywords_List] := Cell[CellGroupData[
  Prepend[Cell[#, "Keywords"] & /@ keywords, Cell["Keywords", "KeywordsSection"]],
  Closed
]];

guideExpectedFragments = {
  "MagneticTB constructs analytical tight-binding Hamiltonians",
  "Tutorials: start here",
  "Model setup and inspection",
  "Crystal, symmetry, and group theory",
  "Hamiltonians and band calculations",
  "Band fitting and VASP comparison",
  "Fitting MagneticTB bands to VASP",
  "MagneticTB/tutorial/BandFitting",
  "Plotting and reciprocal-space paths",
  If[
    MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
    Nothing,
    "Physical properties and topology"
  ],
  "The MagneticTB 1.0 functions are included under MagneticTBOld`",
  "Getting Started with MagneticTB",
  "General MagneticTB Examples",
  "Crystal structure, Brillouin zone, and k paths",
  If[
    MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
    Nothing,
    "Real-space Hamiltonians, surfaces, and topology"
  ],
  "Band corepresentations with MSGCorep",
  "High-symmetry tetragonal antiferromagnet: spin symmetry and SOC",
  "MagneticTB/tutorial/ContinuousSymmetry",
  "MagneticTB/tutorial/LegacyBackend",
  If[
    MemberQ[deferredDocumentationPageNames, "plotSurfaceSpectrum"],
    Nothing,
    "plotSurfaceSpectrum"
  ],
  "gray",
  "bnsdict",
  "typeI",
  "typeIII",
  "typeIV",
  "graylayer",
  "grayrod",
  "fittingTB",
  "MagneticTB/guide/MagneticTB"
};

writeGuidePage[page_Notebook] := writeGuideSourceOnly[
  page,
  guideExpectedFragments
];

usageCell[syntax_String, description_String] := Cell[
  TextData[{
    Cell["   ", "ModInfo"],
    Cell[BoxData[codeBoxes[syntax]], "InlineFormula"],
    " \[LineSeparator]",
    description
  }],
  "Usage"
];

optionTableCell[rows_List] := Cell[
  BoxData@GridBox[
    Map[
      {
        Cell["      ", "TableRowIcon"],
        codeBoxes[#[[1]]],
        codeBoxes[#[[2]]],
        Cell[#[[3]], "TableText"]
      } &,
      rows
    ],
    GridBoxAlignment -> {"Columns" -> {{Left}}},
    GridBoxDividers -> {"Rows" -> {{False}}},
    GridBoxSpacings -> {
      "Columns" -> {
        Offset[0.27999999999999997],
        {Offset[0.7]},
        Offset[0.27999999999999997]
      },
      "Rows" -> {
        Offset[0.2],
        {Offset[0.8]},
        Offset[0.2]
      }
    }
  ],
  "3ColumnTableMod"
];

descriptionTableCell[rows_List] := Cell[
  BoxData@GridBox[
    Map[
      {
        Cell["      ", "TableRowIcon"],
        codeBoxes[#[[1]]],
        Cell[#[[2]], "TableText"],
        Cell[#[[3]], "TableText"]
      } &,
      rows
    ],
    GridBoxAlignment -> {"Columns" -> {{Left}}},
    GridBoxDividers -> {"Rows" -> {{False}}},
    GridBoxSpacings -> {
      "Columns" -> {
        Offset[0.27999999999999997],
        {Offset[0.7]},
        Offset[0.27999999999999997]
      },
      "Rows" -> {
        Offset[0.2],
        {Offset[0.4]},
        Offset[0.2]
      }
    }
  ],
  "3ColumnTableMod"
];

exampleMetadata[example_List] := If[
  Length[example] >= 4 && AssociationQ[example[[4]]],
  example[[4]],
  <||>
];

exampleInputCell[
    example_List,
    initialization_String
  ] := Module[{metadata, messages, tags, initializationLine, code},
  metadata = exampleMetadata[example];
  messages = Lookup[metadata, "ExpectedMessages", {}];
  tags = "ExpectedMessage:" <> # & /@ messages;
  initializationLine = StringTrim[initialization];
  code = If[
    initializationLine === "",
    example[[2]],
    If[!StringEndsQ[initializationLine, ";"],
      initializationLine = initializationLine <> ";"
    ];
    initializationLine <> "\n" <> example[[2]]
  ];
  If[
    tags === {},
    Cell[BoxData[codeBoxes[code]], "Input"],
    Cell[BoxData[codeBoxes[code]], "Input", CellTags -> tags]
  ]
];

(* Tutorial failure examples do not need a repeated initialization string.
   Preserve their expected-message CellTags instead of leaving an unevaluated
   helper expression inside the Notebook source. *)
exampleInputCell[example_List] := exampleInputCell[example, ""];

exampleMessageCells[example_List] := Module[{metadata, messageText},
  metadata = exampleMetadata[example];
  messageText = Lookup[metadata, "MessageText", None];
  If[
    StringQ[messageText],
    {Cell[BoxData[codeBoxes[messageText]], "Message"]},
    {}
  ]
];

exampleCells[
    examples_List,
    initialization_String
  ] := Flatten@Map[
  Function[example,
    Join[
      {
        Cell[example[[1]], "ExampleText"],
        exampleInputCell[example, initialization]
      },
      exampleMessageCells[example],
      If[
        Length[example] >= 3,
        {Cell[BoxData[outputBoxes[example[[3]]]], "Output"]},
        {}
      ]
    ]
  ],
  examples
];

(* Match the Wolfram reference-page template: load the documented context
   once before any option or basic example is evaluated.  Keeping this in a
   dedicated ExampleInitialization cell makes a page runnable in a fresh
   kernel without repeating Needs in every individual example. *)
examplesInitializationSection[spec_Association] := Module[
  {initialization},
  initialization = Lookup[
    spec,
    "Initialization",
    "Needs[\"MagneticTB`\"]"
  ];
  {Cell[CellGroupData[{
    Cell["Examples Initialization", "ExamplesInitializationSection"],
    Cell[BoxData[codeBoxes[initialization]], "ExampleInitialization"]
  }, Open]]}
];

(* Follow the Wolfram function-page template exactly.  Additional examples
   live in one ExtendedExamplesSection immediately after Basic Examples.
   Options, applications, properties, and issues are child ExampleSection
   groups; the documentation builder then produces the standard More Examples
   link and keeps these sections ahead of See Also and Related Guides. *)
optionalExtendedExampleSection[
    spec_Association,
    key_String,
    title_String,
    initialization_String
  ] :=
  If[
    KeyExistsQ[spec, key] && Lookup[spec, key, {}] =!= {},
    {
      Cell[CellGroupData[
        Join[
          {Cell[title, "ExampleSection"]},
          exampleCells[spec[key], initialization]
        ],
        Open
      ]]
    },
    {}
  ];

(* Wolfram reference pages give each documented option its own subsection.
   Keep that structure here as well: the option table remains in Details,
   while the executable examples are grouped by option name.  Matching is
   deliberately strict so a newly added option cannot be hidden in an
   unrelated paragraph or silently shipped without an example. *)
ClearAll[
  optionExampleName,
  optionExampleSection,
  extendedExamplesSection
];

optionExampleName[example_List, optionNames_List] := Module[
  {text, starts, contains},
  text = StringTrim[example[[1]]];
  starts = Select[
    optionNames,
    Function[name, With[{spellings = {name, "\"" <> name <> "\""}},
      AnyTrue[
        spellings,
        Function[spelling,
          StringStartsQ[text, {
            spelling <> " ->",
            spelling <> " -",
            spelling <> ":"
          }]
        ]
      ]
    ]]
  ];
  If[Length[starts] === 1, Return[First[starts]]];
  contains = Select[
    optionNames,
    Function[name, With[{spellings = {name, "\"" <> name <> "\""}},
      AnyTrue[spellings, Function[spelling, StringContainsQ[text, spelling]]]
    ]]
  ];
  If[Length[contains] === 1, First[contains], Missing["NoOptionName"]]
];

optionExampleSection[
    spec_Association,
    initialization_String
  ] := Module[
  {
    optionNames, examples, assignments, unmatched, missing, groups,
    allowPartial
  },
  If[
    !KeyExistsQ[spec, "OptionTable"] || spec["OptionTable"] === {},
    Return[{}]
  ];
  optionNames = StringTrim[#, "\""] & /@ spec["OptionTable"][[All, 1]];
  examples = Lookup[spec, "OptionExamples", {}];
  assignments = optionExampleName[#, optionNames] & /@ examples;
  unmatched = Pick[examples, MissingQ /@ assignments];
  missing = Select[optionNames, FreeQ[assignments, #] &];
  allowPartial = TrueQ[Lookup[spec, "AllowPartialOptionExamples", False]];
  If[allowPartial && examples === {}, Return[{}]];
  If[unmatched =!= {} || (missing =!= {} && !allowPartial),
    Print[InputForm@<|
      "DocumentationPage" -> spec["Name"],
      "UnmatchedOptionExamples" -> unmatched,
      "OptionsWithoutExamples" -> missing
    |>];
    Exit[15]
  ];
  If[allowPartial,
    optionNames = Select[optionNames, MemberQ[assignments, #] &]
  ];
  groups = Map[
    Function[name,
      Cell[CellGroupData[
        Join[
          {Cell[name, "ExampleSubsection"]},
          exampleCells[
            Pick[examples, assignments, name],
            initialization
          ]
        ],
        Open
      ]]
    ],
    optionNames
  ];
  {
    Cell[CellGroupData[
      Prepend[groups, Cell["Options", "ExampleSection"]],
      Open
    ]]
  }
];

extendedExamplesSection[
    spec_Association,
    initialization_String
  ] := Module[{sections},
  sections = Join[
    optionExampleSection[spec, initialization],
    optionalExtendedExampleSection[
      spec,
      "Applications",
      "Applications",
      initialization
    ],
    optionalExtendedExampleSection[
      spec,
      "Properties",
      "Properties & Relations",
      initialization
    ],
    optionalExtendedExampleSection[
      spec,
      "Issues",
      "Possible Issues",
      initialization
    ]
  ];
  If[
    sections === {},
    {},
    {
      Cell[CellGroupData[
        Prepend[
          sections,
          Cell[
            "More Examples",
            "ExtendedExamplesSection",
            CellTags -> "ExtendedExamples"
          ]
        ],
        Open
      ]]
    }
  ]
];

linkListingCell[style_String, links_List] := Cell[
  TextData@Riffle[
    documentationLink[#[[1]], #[[2]]] & /@ links,
    " \[FilledVerySmallSquare] "
  ],
  style
];

(* Keep the installation workflow connected to the main guide at the point
   where readers need it.  This is a real paclet link, not example code that
   users must evaluate. *)
gettingStartedHelpHomeCell[] := Cell[
  TextData[{
    documentationLink[
      "MagneticTB help home",
      "MagneticTB/guide/MagneticTB"
    ]
  }],
  "Text"
];

ensureTutorialHelpHomeLinks[] := Scan[
  Function[path,
    sourceTextValue = FromCharacterCode[BinaryReadList[path], "UTF8"];
    If[!StringContainsQ[
        sourceTextValue,
        "paclet:MagneticTB/guide/MagneticTB"
      ],
      historyMarker =
        "Cell[CellGroupData[{Cell[\"History\", \"HistorySection\"]";
      If[!StringContainsQ[sourceTextValue, historyMarker],
        Print["Tutorial source has no History insertion marker: ", path];
        Quit[8]
      ];
      helpHomeFragment = StringJoin[
        "Cell[\"MagneticTB help home\", \"Section\"], ",
        "Cell[TextData[{ButtonBox[\"MagneticTB\", ",
        "BaseStyle -> \"Link\", ButtonData -> ",
        "\"paclet:MagneticTB/guide/MagneticTB\"]}], \"Text\"], "
      ];
      sourceTextValue = StringReplace[
        sourceTextValue,
        historyMarker -> helpHomeFragment <> historyMarker,
        1
      ];
      writeUTF8TextIfChanged[path, sourceTextValue]
    ]
  ],
  Sort@FileNames[
    "*.source.wl",
    FileNameJoin[{sourceRoot, "Tutorials"}],
    Infinity
  ]
];

If[MatchQ[
    Environment["MAGNETICTB_DOC_NAVIGATION_ONLY"],
    1 | "1" | True
  ],
  ensureTutorialHelpHomeLinks[];
  Quit[0]
];

If[MatchQ[
    Environment["MAGNETICTB_DOC_GUIDE_ONLY"],
    1 | "1" | True
  ],
  generatorHeld = ToExpression[
    FromCharacterCode[
      BinaryReadList[ExpandFileName[$InputFileName]],
      "UTF8"
    ],
    InputForm,
    HoldComplete
  ];
  guidePageOnly = FirstCase[
    generatorHeld,
    HoldPattern[Set[guidePage, page_Notebook]] :> page,
    Missing["GuidePageNotFound"],
    Infinity
  ];
  If[MissingQ[guidePageOnly],
    Print["Could not locate the held MagneticTB guide page."];
    Quit[9]
  ];
  If[Head[guidePageOnly] =!= Notebook,
    Print["The held MagneticTB guide page did not evaluate to Notebook."];
    Quit[10]
  ];
  writeGuidePage[guidePageOnly];
  Quit[0]
];

functionPage[spec_Association] := Notebook[
  Join[{
    Cell[CellGroupData[
      Join[
        {Cell[spec["Name"], "ObjectName"]},
        usageCell @@@ spec["Usage"],
        Cell[#, "Notes"] & /@ spec["Notes"],
        If[
          KeyExistsQ[spec, "OptionTable"] && spec["OptionTable"] =!= {},
          {optionTableCell[spec["OptionTable"]]},
          {}
        ]
      ],
      Open
    ]]
    },
    examplesInitializationSection[spec],
    {
      Cell[CellGroupData[
        Join[
          {Cell["Basic Examples", "PrimaryExamplesSection"]},
          exampleCells[
            spec["Examples"],
            ""
          ]
        ],
        Open
      ]]
    },
    extendedExamplesSection[spec, ""],
    {
    Cell[CellGroupData[{
      Cell["Related Functions", "SeeAlsoSection"],
      linkListingCell["SeeAlso", spec["SeeAlso"]]
    }, Open]],
    Cell[CellGroupData[{
      Cell["Related Guides", "MoreAboutSection"],
      linkListingCell["MoreAbout", {{"MagneticTB", "MagneticTB/guide/MagneticTB"}}]
    }, Open]],
    Cell[CellGroupData[{
      Cell["Related Tutorials", "TutorialsSection"],
      linkListingCell[
        "Tutorials",
        Lookup[
          spec,
          "Tutorials",
          {{"Getting Started with MagneticTB", "MagneticTB/tutorial/GettingStarted"}}
        ]
      ]
    }, Open]],
    historyCells[],
    categorizationCells[
      "Symbol",
      Lookup[spec, "Context", "MagneticTB`"],
      "MagneticTB/ref/" <> spec["Name"]
    ],
    keywordCells[spec["Keywords"]]
    }
  ],
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> spec["Name"],
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"},
    "FunctionPageStylesExt.nb",
    CharacterEncoding -> "UTF-8"
  ]
];

(* -------------------------------------------------------------------------
   Phase 1: stable reference-page specifications and reusable model snippets.

   Each specification below is data.  functionPage turns that data into the
   low-level Notebook expression, so page authors should normally edit the
   specification rather than hand-editing generated .source.wl or .nb files.

   Keep visible examples result-oriented.  Prefer the shortest complete call
   to the public API, remove test-only helpers and temporary bookkeeping, and
   show the physical result directly.  Do not hide model initialization when
   a copied example needs it in a fresh kernel.
   ------------------------------------------------------------------------- *)

grapheneInitCode = StringRiffle[{
  "grapheneOperations = msgop[gray[191]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 3},",
  "  wyckoffposition -> {",
  "    {{1/3, 2/3, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> grapheneOperations,",
  "  basisFunctions -> {{\"pz\"}},",
  "  InitialBondShells -> 3,",
  "  GenerateSymmetryGroup -> False",
  "];"
}, "\n"];

grapheneHamiltonianCode = StringRiffle[{
  "grapheneOnsite = symham[1];",
  "grapheneNearest = symham[2];",
  "grapheneNextNearest = symham[3];",
  "grapheneHamiltonian =",
  "  grapheneOnsite + grapheneNearest + grapheneNextNearest;",
  "MatrixForm[",
  "  FullSimplify[",
  "    ExpToTrig[grapheneHamiltonian],",
  "    Element[{kx, ky, kz}, Reals]",
  "  ]",
  "]"
}, "\n"];

graphenePathCode = StringRiffle[{
  "graphenePath = {",
  "  {{{0, 0, 0}, {0, 1/2, 0}}, {\"G\", \"M\"}},",
  "  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
  "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"G\"}}",
  "};"
}, "\n"];

simpleInitCode = StringReplace[
  grapheneInitCode,
  "  basisFunctions -> {{\"pz\"}},\n  InitialBondShells -> 3,\n  GenerateSymmetryGroup -> False" ->
    "  basisFunctions -> {{\"pz\"}}"
];

grapheneBondShellInitCode = StringReplace[
  grapheneInitCode,
  "  InitialBondShells -> 3,\n  GenerateSymmetryGroup -> False" ->
    "  InitialBondShells -> 3"
];

chainInitCode = StringRiffle[{
  "chainOperations = msgop[gray[75]];",
  "init[",
  "  lattice -> {{a,0,0},{0,a,0},{0,0,c}},",
  "  lattpar -> {a -> 3, c -> 1},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> chainOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];"
}, "\n"];

generatedC4InitCode = StringRiffle[{
  "c4z = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> {",
  "    {\"C4z\", c4z, {0, 0, 0}, \"F\"}",
  "  },",
  "  basisFunctions -> {{\"px\", \"py\"}},",
  "  GenerateSymmetryGroup -> True",
  "];"
}, "\n"];

cubicPOptionInitCode[mode_String] := StringRiffle[{
  "cubicOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, a}",
  "  },",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {",
  "    {{1/2, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> cubicOperations,",
  If[
    mode === "DirectProduct",
    "  basisFunctions -> {{\"px\", \"py\", \"pz\"}},",
    "  basisFunctions -> {{\"px\"}},"
  ],
  "  RepresentationMode -> \"" <> mode <> "\"",
  "];"
}, "\n"];

(* Three-band kagome QNP model from Feng et al., Phys. Rev. B 109,
   045423 (2024).  The paper labels an in-plane spinor on one kagome site;
   the present database representative is related to it by a global in-plane
   spin rotation.  Induction then transports this one local state to the two
   other sites.  In MagneticTB's ordered directed-bond gauge, the paper's
   (epsilon0,t1,t2)=(-0.15,0.155,-0.23) corresponds to
   (e1,t1,r1)=(-0.15,-0.155/2,0.23/2). *)
kagomeQNPInducedInitCode = StringRiffle[{
  "kagomeOperations = msgop[typeIII[{191, 239}]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 8},",
  "  wyckoffposition -> {",
  "    {{1/2, 0, 0}, {1, 0, 0}}",
  "  },",
  "  symminformation -> kagomeOperations,",
  "  basisFunctions -> {{{1, -1}/Sqrt[2]}},",
  "  RepresentationMode -> \"Induced\",",
  "  InitialBondShells -> 3",
  "];",
  "kagomeHamiltonian = Total[symham /@ Range[3]];",
  "kagomePath = {",
  "  {{{0, 1/2, 0}, {0, 0, 0}}, {\"M\", \"G\"}},",
  "  {{{0, 0, 0}, {1/3, 1/3, 0}}, {\"G\", \"K\"}},",
  "  {{{1/3, 1/3, 0}, {0, 1/2, 0}}, {\"K\", \"M\"}}",
  "};",
  "bandplot[",
  "  kagomePath,",
  "  80,",
  "  kagomeHamiltonian,",
  "  {e1 -> -0.15, t1 -> -0.155/2, r1 -> 0.23/2}",
  "]"
}, "\n"];

p3HermitianInitCode = StringRiffle[{
  "p3Operations = msgop[bnsdict[{143, 1}]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> p3Operations,",
  "  basisFunctions -> {{\"s\"}}",
  "];"
}, "\n"];

p3HermitianShellReference = {{
  (-I t1 + t2)/Exp[I MagneticTB`kx] +
    (-I t1 + t2)/Exp[I MagneticTB`ky] +
    Exp[I (MagneticTB`kx + MagneticTB`ky)] (-I t1 + t2) +
    Exp[I MagneticTB`kx] (I t1 + t2) +
    Exp[-I (MagneticTB`kx + MagneticTB`ky)] (I t1 + t2) +
    Exp[I MagneticTB`ky] (I t1 + t2)
}};

p3NonHermitianShellReference = {{
  Exp[I MagneticTB`kx] (I t1 + t2) +
    Exp[-I (MagneticTB`kx + MagneticTB`ky)] (I t1 + t2) +
    Exp[I MagneticTB`ky] (I t1 + t2) +
    (I t3 + t4)/Exp[I MagneticTB`kx] +
    (I t3 + t4)/Exp[I MagneticTB`ky] +
    Exp[I (MagneticTB`kx + MagneticTB`ky)] (I t3 + t4)
}};

p3HermitianComparisonReference = Column[{
  Row[{
    "Hermitian -> True: ",
    MatrixForm[p3HermitianShellReference]
  }],
  Row[{
    "Hermitian -> False: ",
    MatrixForm[p3NonHermitianShellReference]
  }]
}];

(* orbitalTable is now a presentation Grid.  Keep the full public display in
   examples; the underlying Association/List records remain private model
   data and are not exposed as a second public return form. *)
compactOrbitalTable[table_Grid] := table;
compactOrbitalTableCode = "orbitalTable[]";

(* Association values remain the internal data contract.  When a tutorial
   needs a compact human-readable summary, display those key/value pairs as a
   plain Grid instead of embedding an interactive Dataset/DynamicModule. *)
ClearAll[documentationKeyValueGrid];
documentationKeyValueGrid[association_Association] := Grid[
  List @@@ Normal[association],
  Frame -> All,
  Alignment -> Left
];

spinfulGrapheneInitCode = StringRiffle[{
  "spinfulGrapheneOperations = msgop[gray[191]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 3},",
  "  wyckoffposition -> {",
  "    {{1/3, 2/3, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> spinfulGrapheneOperations,",
  "  basisFunctions -> {{\"pzup\", \"pzdn\"}},",
  "  InitialBondShells -> 3,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "spinfulGrapheneHamiltonian = Sum[symham[i], {i, 3}];",
  "spinfulGraphenePath = {",
  "  {{{0, 0, 0}, {0, 1/2, 0}}, {\"G\", \"M\"}},",
  "  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
  "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"G\"}}",
  "};",
  "spinfulGrapheneParameters = {",
  "  e1 -> 0, r1 -> 0.02, r2 -> 0, t1 -> 0.5",
  "};"
}, "\n"];

explicitInitCode = StringRiffle[{
  "initfromrep[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> {",
  "    {\"E\", IdentityMatrix[3], {0, 0, 0}, \"F\"},",
  "    {",
  "      \"C4z\",",
  "      {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},",
  "      {0, 0, 0},",
  "      \"F\"",
  "    },",
  "    {\"C2z\", DiagonalMatrix[{-1, -1, 1}], {0, 0, 0}, \"F\"},",
  "    {",
  "      \"C4z^-1\",",
  "      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},",
  "      {0, 0, 0},",
  "      \"F\"",
  "    }",
  "  },",
  "  repinformation -> {{",
  "    IdentityMatrix[2],",
  "    {{0, -1}, {1, 0}},",
  "    -IdentityMatrix[2],",
  "    {{0, 1}, {-1, 0}}",
  "  }},",
  "  orbitalLabels -> {{\"px\", \"py\"}}",
  "];"
}, "\n"];

explicitBondShellInitCode = StringReplace[
  explicitInitCode,
  "  orbitalLabels -> {{\"px\", \"py\"}}" ->
    "  orbitalLabels -> {{\"px\", \"py\"}},\n" <>
      "  InitialBondShells -> 2"
];

explicitModeDirectInitCode = StringRiffle[{
  "initfromrep[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> {",
  "    {\"E\", IdentityMatrix[3], {0, 0, 0}, \"F\"},",
  "    {",
  "      \"C4z\",",
  "      {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},",
  "      {0, 0, 0},",
  "      \"F\"",
  "    },",
  "    {\"C2z\", DiagonalMatrix[{-1, -1, 1}], {0, 0, 0}, \"F\"},",
  "    {",
  "      \"C4z^-1\",",
  "      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},",
  "      {0, 0, 0},",
  "      \"F\"",
  "    }",
  "  },",
  "  repinformation -> {{",
  "    IdentityMatrix[2],",
  "    {{0, -1}, {1, 0}},",
  "    -IdentityMatrix[2],",
  "    {{0, 1}, {-1, 0}}",
  "  }},",
  "  orbitalLabels -> {{\"px\", \"py\"}},",
  "  RepresentationMode -> \"DirectProduct\"",
  "];"
}, "\n"];

explicitModeInducedInitCode = StringRiffle[{
  "initfromrep[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> {",
  "    {\"E\", IdentityMatrix[3], {0, 0, 0}, \"F\"},",
  "    {",
  "      \"C4z\",",
  "      {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},",
  "      {0, 0, 0},",
  "      \"F\"",
  "    },",
  "    {\"C2z\", DiagonalMatrix[{-1, -1, 1}], {0, 0, 0}, \"F\"},",
  "    {",
  "      \"C4z^-1\",",
  "      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},",
  "      {0, 0, 0},",
  "      \"F\"",
  "    }",
  "  },",
  "  repinformation -> {",
  "    <|",
  "      \"ReferenceSiteIndex\" -> 1,",
  "      \"SiteSymmetryMatrices\" -> {{{1}}, {{-1}}}",
  "    |>",
  "  },",
  "  orbitalLabels -> {{\"px\"}},",
  "  RepresentationMode -> \"Induced\"",
  "];"
}, "\n"];

symmetrizationOptionCode = StringRiffle[{
  "c4SymmetrizationOperations = msgop[gray[75]];",
  "c4SymmetrizationLattice = DiagonalMatrix[{a, a, c}];",
  "c4SymmetrizationParameters = {a -> 1, c -> 2};",
  "c4SymmetrizationSites = {{{0, 0, 0}, {0, 0, 0}}};",
  "c4SymmetrizationBasis = {\"px\", \"py\"};",
  "c4SymmetrizationData = symmetrizationHRInit[",
  "  \"Lattice\" -> c4SymmetrizationLattice,",
  "  \"LattPar\" -> c4SymmetrizationParameters,",
  "  \"WyckoffPosition\" -> c4SymmetrizationSites,",
  "  \"Symminformation\" -> c4SymmetrizationOperations,",
  "  \"BasisFunctions\" -> c4SymmetrizationBasis,",
  "  \"Software\" -> \"VASP\"",
  "];"
}, "\n"];

oldInitCode = StringRiffle[{
  "grapheneOperationsOld = msgopold[grayold[191]];",
  "initold[",
  "  latticeold -> {",
  "    {aold, 0, 0},",
  "    {-aold/2, Sqrt[3] aold/2, 0},",
  "    {0, 0, cold}",
  "  },",
  "  lattparold -> {aold -> 1, cold -> 3},",
  "  wyckoffpositionold -> {",
  "    {{1/3, 2/3, 0}, {0, 0, 0}}",
  "  },",
  "  symminformationold -> grapheneOperationsOld,",
  "  basisFunctionsold -> {{\"pz\"}},",
  "  debugQold -> False",
  "];"
}, "\n"];

oldDebugInitCode = StringReplace[
  oldInitCode,
  "  debugQold -> False" -> "  debugQold -> True"
];

(* Generate the visible outputs from the same source package that is being
   documented.  This prevents hand-written Hamiltonians, tables, grids, and
   plots from drifting away from the executable examples. *)
runtimeLoader = FileNameJoin[{projectRoot, "MagneticTB", "Kernel", "init.wl"}];
If[!FileExistsQ[runtimeLoader],
  Print["Missing MagneticTB runtime loader: ", runtimeLoader];
  Exit[4]
];

{runtimeLoadResult, runtimeLoadMessages} = Block[{$MessageList = {}},
  {Get[runtimeLoader], $MessageList}
];
If[runtimeLoadResult === $Failed || runtimeLoadMessages =!= {},
  Print["Unable to load the documented runtime: ", runtimeLoadMessages];
  Exit[5]
];

ClearAll[checkedDocumentationEvaluation];
SetAttributes[checkedDocumentationEvaluation, HoldAllComplete];
checkedDocumentationEvaluation[label_String, expression_] := Module[
  {},
  Block[{Print = Function[Null]}, expression]
];

(* Generate only the crystal tutorial when requested; do not evaluate any
   unrelated reference-page or tutorial examples. *)
If[Environment["MAGNETICTB_DOC_PAGE"] === "Tutorials/CrystalAndKPaths",
  Get[FileNameJoin[{
    scriptDirectory, "Examples", "CrystalAndKPathsDocumentation.wl"
  }]];
  If[Head[testingCrystalTutorialPage] =!= Notebook, Exit[10]];
  selectedCrystalSourcePath = FileNameJoin[{
    sourceRoot, "Tutorials", "CrystalAndKPaths.source.wl"
  }];
  writeUTF8TextIfChanged[
    selectedCrystalSourcePath,
    ToString[testingCrystalTutorialPage, InputForm, PageWidth -> Infinity]
  ];
  writeUTF8TextIfChanged[
    FileNameJoin[{expectedRoot, "Tutorials", "CrystalAndKPaths.expected.txt"}],
    StringRiffle[crystalKagomeExpectedFragments, "\n"] <> "\n"
  ];
  Print[selectedCrystalSourcePath];
  Exit[0]
];

(* Corepresentation tables are Print results, captured by this page's own
   authoring source.  Do not obtain them from an unevaluated reference spec. *)
If[Environment["MAGNETICTB_DOC_PAGE"] === "Tutorials/Corepresentations",
  Get[FileNameJoin[{
    scriptDirectory, "Examples", "CorepresentationsDocumentation.wl"
  }]];
  If[Head[testingCorepTutorialPage] =!= Notebook, Exit[10]];
  selectedCorepSourcePath = FileNameJoin[{
    sourceRoot, "Tutorials", "Corepresentations.source.wl"
  }];
  writeUTF8TextIfChanged[
    selectedCorepSourcePath,
    ToString[testingCorepTutorialPage, InputForm, PageWidth -> Infinity]
  ];
  writeUTF8TextIfChanged[
    FileNameJoin[{expectedRoot, "Tutorials", "Corepresentations.expected.txt"}],
    StringRiffle[corepTutorialExpectedFragments, "\n"] <> "\n"
  ];
  Print[selectedCorepSourcePath];
  Exit[0]
];

(* Incremental build for the ContinuousSymmetry tutorial.  The dedicated
   human-maintained source evaluates only the type-IV antiferromagnetic
   comparison used by this page. *)
If[MatchQ[
    Environment["MAGNETICTB_DOC_CONTINUOUS_TUTORIAL_ONLY"],
    1 | "1" | True
  ],
  continuousDocumentationSource = FileNameJoin[{
    scriptDirectory,
    "Examples",
    "ContinuousSymmetryDocumentation.wl"
  }];
  If[!FileExistsQ[continuousDocumentationSource],
    Print[
      "Missing ContinuousSymmetry documentation source: ",
      continuousDocumentationSource
    ];
    Exit[9]
  ];
  Get[continuousDocumentationSource];
  If[Head[continuousSymmetryDocumentationPage] =!= Notebook,
    Print[
      "The ContinuousSymmetry tutorial did not evaluate to a Notebook."
    ];
    Exit[10]
  ];
  selectedContinuousSourcePath = FileNameJoin[{
    sourceRoot,
    "Tutorials",
    "ContinuousSymmetry.source.wl"
  }];
  selectedContinuousExpectedPath = FileNameJoin[{
    expectedRoot,
    "Tutorials",
    "ContinuousSymmetry.expected.txt"
  }];
  Scan[
    If[!DirectoryQ[#],
      CreateDirectory[#, CreateIntermediateDirectories -> True]
    ] &,
    {
      DirectoryName[selectedContinuousSourcePath],
      DirectoryName[selectedContinuousExpectedPath]
    }
  ];
  writeUTF8TextIfChanged[
    selectedContinuousSourcePath,
    ToString[
      continuousSymmetryDocumentationPage,
      InputForm,
      PageWidth -> 120
    ]
  ];
  writeUTF8TextIfChanged[
    selectedContinuousExpectedPath,
    StringRiffle[continuousSymmetryExpectedFragments, "\n"] <> "\n"
  ];
  ensureTutorialHelpHomeLinks[];
  Print[selectedContinuousSourcePath];
  Exit[0]
];

(* Incremental build for the VASP fitting tutorial.  Its source specification
   evaluates only the fitting fixtures needed by this page and skips all
   unrelated tutorial fixtures. *)
If[MatchQ[
    Environment["MAGNETICTB_DOC_BAND_FITTING_TUTORIAL_ONLY"],
    1 | "1" | True
  ],
  testingDocumentationSource = FileNameJoin[{
    scriptDirectory,
    "Examples",
    "TestingPublicAPIDocumentation.wl"
  }];
  If[!FileExistsQ[testingDocumentationSource],
    Print["Missing testing documentation source: ",
      testingDocumentationSource];
    Exit[9]
  ];
  Get[testingDocumentationSource];
  If[Head[testingBandFittingTutorialPage] =!= Notebook,
    Print["The VASP fitting tutorial did not evaluate to a Notebook."];
    Exit[10]
  ];
  selectedTestingSourcePath = FileNameJoin[{
    sourceRoot,
    "Tutorials",
    "BandFitting.source.wl"
  }];
  selectedTestingExpectedPath = FileNameJoin[{
    expectedRoot,
    "Tutorials",
    "BandFitting.expected.txt"
  }];
  Scan[
    If[!DirectoryQ[#],
      CreateDirectory[#, CreateIntermediateDirectories -> True]
    ] &,
    {
      DirectoryName[selectedTestingSourcePath],
      DirectoryName[selectedTestingExpectedPath]
    }
  ];
  writeUTF8TextIfChanged[
    selectedTestingSourcePath,
    ToString[
      testingBandFittingTutorialPage,
      InputForm,
      PageWidth -> 120
    ]
  ];
  writeUTF8TextIfChanged[
    selectedTestingExpectedPath,
    StringRiffle[{
      "Fitting MagneticTB bands to VASP",
      "Read and inspect the three VASP bands",
      "Build the symmetry-generated model",
      "Fit the onsite and nearest-neighbor model",
      "Add second-neighbor hopping",
      "Fit the energy range of interest",
      "Fit near the K point",
      "Judge the plotted bands, not only the loss",
      "MagneticTB/tutorial/BandFitting"
    }, "\n"] <> "\n"
  ];
  ensureTutorialHelpHomeLinks[];
  Print[selectedTestingSourcePath];
  Exit[0]
];

(* Targeted build for one testing-line public API.  The specification file
   uses MAGNETICTB_DOC_PAGE to evaluate only that page and skips its unrelated
   tutorial fixtures in this mode. *)
If[MatchQ[
    Environment["MAGNETICTB_DOC_TESTING_PUBLIC_ONLY"],
    1 | "1" | True
  ],
  testingDocumentationSource = FileNameJoin[{
    scriptDirectory,
    "Examples",
    "TestingPublicAPIDocumentation.wl"
  }];
  If[!FileExistsQ[testingDocumentationSource],
    Print["Missing testing public-API documentation source: ",
      testingDocumentationSource];
    Exit[9]
  ];
  Get[testingDocumentationSource];
  requestedTestingPage = Environment["MAGNETICTB_DOC_PAGE"];
  If[
    !StringQ[requestedTestingPage] ||
      StringTrim[requestedTestingPage] === "",
    Print["Set MAGNETICTB_DOC_PAGE for the targeted public-API build."];
    Exit[10]
  ];
  requestedTestingName = FileBaseName[requestedTestingPage];
  selectedTestingSpecifications = Select[
    testingFunctionSpecifications,
    Lookup[#, "Name", ""] === requestedTestingName &&
      activeDocumentationSpecificationQ[#] &
  ];
  If[Length[selectedTestingSpecifications] =!= 1,
    Print["The requested testing public API is unavailable: ",
      requestedTestingPage];
    Exit[10]
  ];
  selectedTestingSpec = First[selectedTestingSpecifications];
  selectedTestingPage = functionPage[selectedTestingSpec];
  selectedTestingExpected = DeleteDuplicates@Join[
    {
      selectedTestingSpec["Name"],
      First[selectedTestingSpec["Usage"]][[2]],
      "Examples",
      "Related Functions",
      "MagneticTB/ref/" <> selectedTestingSpec["Name"]
    },
    If[
      Lookup[selectedTestingSpec, "OptionTable", {}] === {},
      {},
      Join[
        If[
          Lookup[selectedTestingSpec, "OptionExamples", {}] === {},
          {},
          {"Options"}
        ],
        StringTrim[#, "\""] & /@
          Lookup[selectedTestingSpec, "OptionTable"][[All, 1]]
      ]
    ],
    Lookup[selectedTestingSpec, "Expected", {}],
    If[Lookup[selectedTestingSpec, "Applications", {}] === {},
      {}, {"Applications"}],
    If[Lookup[selectedTestingSpec, "Properties", {}] === {},
      {}, {"Properties & Relations"}],
    If[Lookup[selectedTestingSpec, "Issues", {}] === {},
      {}, {"Possible Issues"}]
  ];
  selectedTestingSourcePath = FileNameJoin[{
    sourceRoot,
    "ReferencePages",
    "Symbols",
    selectedTestingSpec["Name"] <> ".source.wl"
  }];
  selectedTestingExpectedPath = FileNameJoin[{
    expectedRoot,
    "ReferencePages",
    "Symbols",
    selectedTestingSpec["Name"] <> ".expected.txt"
  }];
  writeUTF8TextIfChanged[
    selectedTestingSourcePath,
    ToString[selectedTestingPage, InputForm, PageWidth -> Infinity]
  ];
  writeUTF8TextIfChanged[
    selectedTestingExpectedPath,
    StringRiffle[selectedTestingExpected, "\n"] <> "\n"
  ];
  Print[selectedTestingSourcePath];
  Exit[0]
];

(* The plotting/property documentation is deliberately buildable without the
   unrelated cyclotomic, crystal-system, and MSGCorep fixtures below.  This is
   the supported incremental path for the six affected reference pages and
   their shared real-space/topology tutorial. *)
If[MatchQ[
    Environment["MAGNETICTB_DOC_PLOT_PROPERTY_ONLY"],
    1 | "1" | True
  ],
  plotPropertyDocumentationSource = FileNameJoin[{
    scriptDirectory,
    "Examples",
    "TestingPlotPropertyDocumentation.wl"
  }];
  If[!FileExistsQ[plotPropertyDocumentationSource],
    Print["Missing plotting/property documentation source: ",
      plotPropertyDocumentationSource];
    Exit[13]
  ];
  Get[plotPropertyDocumentationSource];
  If[
    !ListQ[testingPlotPropertySpecifications] ||
      Length[testingPlotPropertySpecifications] =!= 6 ||
      Sort[Lookup[testingPlotPropertySpecifications, "Name", ""]] =!=
        Sort[{
          "buildRealSpaceHamiltonian", "plotWilsonLoop",
          "plotSurfaceSpectrum", "plotBerryCurvature2D",
          "plotBerryCurvature3D", "plotRealSpaceWavefunction"
        }],
    Print["The incremental plotting/property specification is incomplete."];
    Exit[14]
  ];
  requestedPlotPropertyPage = Environment["MAGNETICTB_DOC_PAGE"];
  If[
    StringQ[requestedPlotPropertyPage] &&
      FileBaseName[requestedPlotPropertyPage] === "RealSpaceAndTopology",
    If[
      MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
      Print[
        "The real-space/topology tutorial is currently preserved under ",
        "DeferredPages and is not part of the active help system."
      ];
      Exit[14]
    ];
    plotPropertyTutorialPath = FileNameJoin[{
      sourceRoot,
      "Tutorials",
      "RealSpaceAndTopology.source.wl"
    }];
    If[!FileExistsQ[plotPropertyTutorialPath],
      Print[
        "Missing existing real-space/topology tutorial source: ",
        plotPropertyTutorialPath
      ];
      Exit[15]
    ];
    plotPropertyTutorialSeedPath = FileNameJoin[{
      documentationRoot,
      "Tutorials",
      "RealSpaceAndTopology.nb"
    }];
    If[!FileExistsQ[plotPropertyTutorialSeedPath],
      Print[
        "Missing installed real-space/topology tutorial seed: ",
        plotPropertyTutorialSeedPath
      ];
      Exit[15]
    ];
    plotPropertyTutorialPage = testingApplyPlotPropertyTutorialUpdate[
      Get[plotPropertyTutorialSeedPath]
    ];
    If[Head[plotPropertyTutorialPage] =!= Notebook,
      Print["Could not update the real-space/topology tutorial."];
      Exit[16]
    ];
    writeUTF8TextIfChanged[
      plotPropertyTutorialPath,
      ToString[plotPropertyTutorialPage, InputForm, PageWidth -> Infinity]
    ];
    writeUTF8TextIfChanged[
      FileNameJoin[{
        expectedRoot,
        "Tutorials",
        "RealSpaceAndTopology.expected.txt"
      }],
      StringRiffle[{
        "Real-space Hamiltonians, surfaces, and topology",
        "Hamiltonians used below",
        "Wilson-loop plots",
        "Surface spectra",
        "C4 higher-order topological-insulator corner state",
        "Magnetic Weyl point",
        "Berry-curvature density",
        "MagneticTB/tutorial/RealSpaceAndTopology"
      }, "\n"] <> "\n"
    ];
    Print[plotPropertyTutorialPath];
    Exit[0]
  ];
  selectedPlotPropertySpecifications = Select[
    testingPlotPropertySpecifications,
    activeDocumentationSpecificationQ
  ];
  If[
    StringQ[requestedPlotPropertyPage] &&
      StringTrim[requestedPlotPropertyPage] =!= "",
    requestedPlotPropertyName = FileBaseName[requestedPlotPropertyPage];
    selectedPlotPropertySpecifications = Select[
      testingPlotPropertySpecifications,
      Lookup[#, "Name", ""] === requestedPlotPropertyName &&
        activeDocumentationSpecificationQ[#] &
    ];
    If[selectedPlotPropertySpecifications === {},
      Print[
        "The requested plotting/property page is not in the incremental set: ",
        requestedPlotPropertyPage
      ];
      Exit[14]
    ]
  ];
  writePlotPropertySource[
      relativePath_String,
      notebook_Notebook,
      expected_List
    ] := Module[{sourcePath, expectedPath},
    sourcePath = FileNameJoin[{sourceRoot, relativePath <> ".source.wl"}];
    expectedPath = FileNameJoin[{expectedRoot, relativePath <> ".expected.txt"}];
    writeUTF8TextIfChanged[
      sourcePath,
      ToString[notebook, InputForm, PageWidth -> Infinity]
    ];
    writeUTF8TextIfChanged[
      expectedPath,
      StringRiffle[expected, "\n"] <> "\n"
    ];
    Print[sourcePath]
  ];
  Print[
    "Writing ",
    Length[selectedPlotPropertySpecifications],
    " plotting/property reference pages."
  ];
  Scan[
    Function[spec, Module[
      {page, baseExpected, optionExpected, additionalExpected, expected},
      Print["Constructing reference page: ", spec["Name"]];
      page = functionPage[spec];
      If[Head[page] =!= Notebook,
        Print[
          "Could not construct plotting/property reference page: ",
          spec["Name"],
          "; head was ",
          Head[page]
        ];
        Exit[17]
      ];
      baseExpected = {
        spec["Name"],
        First[spec["Usage"]][[2]],
        "Examples",
        "Options",
        "Related Functions",
        "MagneticTB/ref/" <> spec["Name"]
      };
      optionExpected =
        StringTrim[#, "\""] & /@ spec["OptionTable"][[All, 1]];
      additionalExpected = Lookup[spec, "Expected", {}];
      expected = DeleteDuplicates@Join[
        baseExpected,
        optionExpected,
        additionalExpected
      ];
      If[!ListQ[expected],
        Print[
          "Could not construct expected-output markers for ",
          spec["Name"],
          "; component heads were ",
          Head /@ {baseExpected, optionExpected, additionalExpected},
          "; result head was ", Head[expected]
        ];
        Exit[18]
      ];
      writePlotPropertySource[
        FileNameJoin[{"ReferencePages", "Symbols", spec["Name"]}],
        page,
        expected
      ]
    ]],
    selectedPlotPropertySpecifications
  ];
  If[
    StringQ[requestedPlotPropertyPage] &&
      StringTrim[requestedPlotPropertyPage] =!= "",
    Exit[0]
  ];
  If[
    MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
    Exit[0]
  ];
  plotPropertyTutorialPath =
    FileNameJoin[{sourceRoot, "Tutorials", "RealSpaceAndTopology.source.wl"}];
  If[!FileExistsQ[plotPropertyTutorialPath],
    Print["Missing existing real-space/topology tutorial source: ",
      plotPropertyTutorialPath];
    Exit[15]
  ];
  plotPropertyTutorialPage = testingApplyPlotPropertyTutorialUpdate[
    Get[plotPropertyTutorialPath]
  ];
  If[Head[plotPropertyTutorialPage] =!= Notebook,
    Print["Could not update the real-space/topology tutorial."];
    Exit[16]
  ];
  writePlotPropertySource[
    FileNameJoin[{"Tutorials", "RealSpaceAndTopology"}],
    plotPropertyTutorialPage,
    {
      "Real-space Hamiltonians, surfaces, and topology",
      "C4 TCI surface Green function",
      "Wilson-loop plots",
      "Surface spectra",
      "C4 higher-order topological-insulator corner state",
      "Magnetic Weyl point",
      "Berry-curvature density",
      "MagneticTB/tutorial/RealSpaceAndTopology"
    }
  ];
  Exit[0]
];

ClearAll[documentationRestoredNested];
documentationRestoredNested[matrix_Association] := Module[
  {restored, rows, columns},
  restored = MagneticTB`RestoreCyclotomicExpression[matrix];
  If[FailureQ[restored], Return[restored]];
  rows = restored["rows"];
  columns = restored["columns"];
  Which[
    rows === 0, {},
    columns === 0, ConstantArray[{}, rows],
    True, Partition[restored["entries"], columns]
  ]
];

cyclotomicDocumentationFixture = checkedDocumentationEvaluation[
  "CyclotomicNullSpace",
  Module[
    {
      iRoot, structuredCompiled, rationalCompiled, rationalMatrix, rref,
      kernel, commonCompiled, commonResult, emptyCompiled, emptyResult,
      zeroCompiled, zeroResult, fullRankCompiled, fullRankResult,
      restoredKernel, lowDegreeFailure, nonExactFailure,
      unsupportedFailure, dimensionFailure, missingDimensionFailure,
      malformedFailure
    },
    iRoot = <|
      "kind" -> "root_sum",
      "constant" -> 0,
      "terms" -> {
        <|"coefficient" -> 1, "order" -> 4, "power" -> 1|>
      }
    |>;
    structuredCompiled = MagneticTB`CompileCyclotomicMatrices[<|
      "coordinate_dimension" -> 3,
      "constraints" -> {
        {{1, iRoot, 0}},
        {{0, 1, 0}}
      },
      "max_cyclotomic_degree" -> 128
    |>];
    rationalCompiled = MagneticTB`CompileCyclotomicMatrices[<|
      "coordinate_dimension" -> 3,
      "constraints" -> {{{1, 2, 3}, {2, 4, 6}}}
    |>];
    rationalMatrix = rationalCompiled["constraints"][[1]];
    rref = MagneticTB`CyclotomicRREF[rationalMatrix];
    kernel = MagneticTB`CyclotomicNullSpace[rationalMatrix];
    commonCompiled = MagneticTB`CompileCyclotomicMatrices[<|
      "coordinate_dimension" -> 3,
      "constraints" -> {{{1, 0, 0}}, {{0, 1, 0}}}
    |>];
    commonResult = MagneticTB`CyclotomicCommonKernel[commonCompiled];
    emptyCompiled = MagneticTB`CompileCyclotomicMatrices[<|
      "coordinate_dimension" -> 3,
      "constraints" -> {}
    |>];
    emptyResult = MagneticTB`CyclotomicCommonKernel[emptyCompiled];
    zeroCompiled = MagneticTB`CompileCyclotomicMatrices[<|
      "coordinate_dimension" -> 0,
      "constraints" -> {}
    |>];
    zeroResult = MagneticTB`CyclotomicCommonKernel[zeroCompiled];
    fullRankCompiled = MagneticTB`CompileCyclotomicMatrices[<|
      "coordinate_dimension" -> 2,
      "constraints" -> {IdentityMatrix[2]}
    |>];
    fullRankResult = MagneticTB`CyclotomicNullSpace[
      fullRankCompiled["constraints"][[1]]
    ];
    restoredKernel = MagneticTB`RestoreCyclotomicExpression[kernel];
    lowDegreeFailure = MagneticTB`CyclotomicCommonNullSpace[
      {{{I, 1}}},
      "MaxCyclotomicDegree" -> 1
    ];
    nonExactFailure = MagneticTB`CyclotomicCommonNullSpace[{{{0.5}}}];
    unsupportedFailure = MagneticTB`CyclotomicCommonNullSpace[
      {{{2^(1/3)}}}
    ];
    dimensionFailure = MagneticTB`CyclotomicCommonNullSpace[{
      {{1, 0}},
      {{1, 0, 0}}
    }];
    missingDimensionFailure = MagneticTB`CyclotomicCommonNullSpace[{}];
    malformedFailure = MagneticTB`CyclotomicNullSpace[{{1, 2, 3}}];
    <|
      "IRoot" -> iRoot,
      "StructuredCompiled" -> structuredCompiled,
      "StructuredSummary" -> KeyTake[
        structuredCompiled,
        {"type", "conductor", "coordinate_dimension", "constraint_count"}
      ],
      "StructuredMatrices" -> (
        documentationRestoredNested /@ structuredCompiled["constraints"]
      ),
      "RationalCompiled" -> rationalCompiled,
      "RREFSummary" -> {
        MatrixForm@documentationRestoredNested[rref["reduced_matrix"]],
        KeyTake[rref, {"pivot_columns", "free_columns", "rank", "nullity"}]
      },
      "CyclotomicRREFSummary" -> With[{
        result = MagneticTB`CyclotomicRREF[
          structuredCompiled["constraints"][[1]]
        ]
      },
        {
          MatrixForm@documentationRestoredNested[result["reduced_matrix"]],
          KeyTake[result, {"pivot_columns", "free_columns", "rank", "nullity"}]
        }
      ],
      "KernelSummary" -> {
        MatrixForm@documentationRestoredNested[kernel["nullspace_rows"]],
        MatrixForm@documentationRestoredNested[kernel["basis_matrix"]],
        KeyTake[
          kernel,
          {
            "pivot_columns", "free_columns", "rank", "nullity",
            "exact_residual_verified"
          }
        ]
      },
      "FullRankSummary" -> {
        KeyTake[
          fullRankResult,
          {"rank", "nullity", "exact_residual_verified"}
        ],
        KeyTake[
          MagneticTB`RestoreCyclotomicExpression[
            fullRankResult["basis_matrix"]
          ],
          {"rows", "columns", "entries"}
        ]
      },
      "CommonSummary" -> {
        MatrixForm@documentationRestoredNested[commonResult["nullspace_rows"]],
        MatrixForm@documentationRestoredNested[commonResult["basis_matrix"]],
        KeyTake[
          commonResult,
          {
            "iteration_nullities", "rank", "nullity",
            "exact_residual_verified"
          }
        ]
      },
      "EmptyCommonSummary" -> {
        MatrixForm@documentationRestoredNested[emptyResult["nullspace_rows"]],
        KeyTake[
          emptyResult,
          {"iteration_nullities", "rank", "nullity"}
        ]
      },
      "ZeroCommonSummary" -> {
        KeyTake[
          zeroResult,
          {"iteration_nullities", "rank", "nullity"}
        ],
        KeyTake[
          MagneticTB`RestoreCyclotomicExpression[zeroResult["basis_matrix"]],
          {"rows", "columns", "entries"}
        ]
      },
      "DirectRows" -> MagneticTB`CyclotomicCommonNullSpace[{
        {
          {1, 0, 0, -Sqrt[2], 0},
          {0, 1, 0, -Sqrt[3], 0}
        },
        {{0, 0, 1, -Sqrt[5], 0}}
      }],
      "DirectResiduals" -> With[{
        blocks = {
          {
            {1, 0, 0, -Sqrt[2], 0},
            {0, 1, 0, -Sqrt[3], 0}
          },
          {{0, 0, 1, -Sqrt[5], 0}}
        },
        rows = MagneticTB`CyclotomicCommonNullSpace[{
          {
            {1, 0, 0, -Sqrt[2], 0},
            {0, 1, 0, -Sqrt[3], 0}
          },
          {{0, 0, 1, -Sqrt[5], 0}}
        }]
      },
        MatrixForm /@ (# . Transpose[rows] & /@ blocks)
      ],
      "DirectEmpty" -> MagneticTB`CyclotomicCommonNullSpace[
        {}, "CoordinateDimension" -> 3
      ],
      "DirectZero" -> MagneticTB`CyclotomicCommonNullSpace[
        {}, "CoordinateDimension" -> 0
      ],
      "DegreeSuccess" -> MagneticTB`CyclotomicCommonNullSpace[
        {{{I, 1}}}, "MaxCyclotomicDegree" -> 2
      ],
      "LowDegreeFailure" -> lowDegreeFailure[[1]],
      "NonExactFailure" -> nonExactFailure[[1]],
      "UnsupportedFailure" -> unsupportedFailure[[1]],
      "DimensionFailure" -> dimensionFailure[[1]],
      "MissingDimensionFailure" -> missingDimensionFailure[[1]],
      "MalformedFailure" -> malformedFailure[[1]],
      "RestoredStructuredMatrix" -> MatrixForm@
        documentationRestoredNested[structuredCompiled["constraints"][[1]]],
      "RestoredKernel" -> {
        MatrixForm@Partition[
          restoredKernel["nullspace_rows", "entries"],
          restoredKernel["nullspace_rows", "columns"]
        ],
        MatrixForm@Partition[
          restoredKernel["basis_matrix", "entries"],
          restoredKernel["basis_matrix", "columns"]
        ]
      },
      "RestoredIElement" -> MagneticTB`RestoreCyclotomicExpression[
        structuredCompiled["context"],
        structuredCompiled["constraints"][[1, "entries", 2]]
      ]
    |>
  ]
];

grapheneFixture = checkedDocumentationEvaluation[
  "Graphene",
  Module[
    {
      operations, session, shells, hamiltonian, path, standardPath, parameters,
      hrPath, hrOptionPath, hrData, bandDataPath, brokenRules,
      covarianceParameters,
      covarianceRepresentation, covarianceReport
    },
    operations = MagneticTB`msgop[MagneticTB`gray[191]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 3},
      MagneticTB`wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"pz"}},
      MagneticTB`InitialBondShells -> 3,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    session = MagneticTB`CurrentModelSession[];
    shells = MagneticTB`symham /@ Range[3];
    hamiltonian = Total[shells];
    covarianceParameters = Complement[
      Variables[hamiltonian],
      {MagneticTB`kx, MagneticTB`ky, MagneticTB`kz}
    ];
    covarianceRepresentation = Join[
      session["ModelSpecification", "Symmetry"],
      session["FullRepresentation"]
    ];
    covarianceReport =
      RepresentationValidation`HamiltonianSymmetryReport[
        hamiltonian,
        covarianceRepresentation,
        {MagneticTB`kx, MagneticTB`ky, MagneticTB`kz},
        covarianceParameters,
        "SampleK" -> {Pi/7, Pi/11, Pi/13}
      ];
    path = {
      {{{0, 0, 0}, {0, 1/2, 0}}, {"G", "M"}},
      {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
      {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}}
    };
    standardPath = MagneticTB`standardKPath[];
    parameters = {Global`e1 -> 0.05, Global`r1 -> 0.02, Global`t1 -> 0.5};
    hrPath = FileNameJoin[{
      $TemporaryDirectory,
      "magnetictb-documentation-graphene_hr.dat"
    }];
    If[FileExistsQ[hrPath], DeleteFile[hrPath]];
    MagneticTB`hop[
      {1, 2, 3},
      parameters,
      "hrExport" -> hrPath,
      "wcc" -> Automatic
    ];
    hrOptionPath = FileNameJoin[{
      $TemporaryDirectory,
      "magnetictb-documentation-graphene-options_hr.dat"
    }];
    If[FileExistsQ[hrOptionPath], DeleteFile[hrOptionPath]];
    MagneticTB`hop[
      {1, 2, 3},
      parameters,
      "Hermitian" -> True,
      "KernelMethod" -> "Iterative",
      "ValidationLevel" -> "Basic",
      "hrExport" -> hrOptionPath,
      "RealDigits" -> 10,
      "wcc" -> Automatic,
      "TranslationTolerance" -> 10^-9
    ];
    hrData = MagneticTB`readHR[hrPath];
    bandDataPath = FileNameJoin[{
      $TemporaryDirectory,
      "magnetictb-documentation-band.csv"
    }];
    If[FileExistsQ[bandDataPath], DeleteFile[bandDataPath]];
    MagneticTB`banddata[
      path,
      10,
      hamiltonian,
      parameters,
      bandDataPath
    ];
    brokenRules = MagneticTB`brokenSymmetryInitRules[{2}];
    <|
      "SelectedOperations" -> operations[[{1, 2, 25}]],
      "AtomPositions" -> session["ModelSpecification", "Metadata", "AtomPositions"],
      "GeneratorLabels" ->
        session["ModelSpecification", "Metadata", "SymmetryInformation"][[
          session["GeneratorIndices"], {1, 4}
        ]],
      "SessionSummary" -> documentationKeyValueGrid@<|
        "RepresentationSource" -> session["RepresentationSource"],
        "RepresentationMode" -> session["RepresentationMode"],
        "SiteOrbits" -> session["ModelSpecification", "SiteOrbits"],
        "LocalBases" -> session["ModelSpecification", "LocalBases"],
        "PreparedBondShells" -> Lookup[session["CompiledBondShells"], "Shell"]
      |>,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "BondGrid" -> MagneticTB`showbonds[2],
      "NextBondGrid" -> MagneticTB`showbonds[3],
      "Onsite" -> shells[[1]],
      "Nearest" -> shells[[2]],
      "NextNearest" -> shells[[3]],
      "Hamiltonian" -> hamiltonian,
      "ReadableHamiltonian" -> FullSimplify[
        ExpToTrig[hamiltonian],
        Element[{MagneticTB`kx, MagneticTB`ky, MagneticTB`kz}, Reals]
      ],
      "ConventionII" -> MagneticTB`symhamII[hamiltonian],
      "NonHermitianNearest" ->
        MagneticTB`symham[2, "Hermitian" -> False],
      "StackedNearest" ->
        MagneticTB`symham[2, "KernelMethod" -> "Stacked"],
      "FullValidationNearest" ->
        MagneticTB`symham[2, "ValidationLevel" -> "Full"],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, hamiltonian],
      "StandardKPath" -> standardPath,
      "AutomaticBandUI" -> MagneticTB`bandManipulate[
        standardPath,
        4,
        hamiltonian
      ],
      "BandPlot" -> MagneticTB`bandplot[
        path, 60, hamiltonian, parameters
      ],
      "BandPlotRange" -> MagneticTB`bandplot[
        path, 60, hamiltonian, parameters,
        MagneticTB`plotRange -> {-0.8, 0.8}
      ],
      "SymmetryRepresentations" ->
        MagneticTB`showSymmetryRepresentations[],
      "SelectedSymmetryRepresentations" ->
        MagneticTB`showSymmetryRepresentations[{1, 3, 25}],
      "TimeReversalRepresentation" ->
        MagneticTB`showSymmetryRepresentations[25],
      "GenericCovarianceReport" -> documentationKeyValueGrid@KeyTake[
        covarianceReport,
        {
          "Convention",
          "SampleK",
          "UnitaryOperationIndices",
          "AntiunitaryOperationIndices",
          "GenericCovarianceResidualsByOperation",
          "MaximumResidual",
          "Verified"
        }
      ],
      "HamiltonianBasis" -> MagneticTB`showHamiltonianBasis[],
      "HamiltonianElementBasis" -> MagneticTB`showHamiltonianBasis[1, 2],
      "HoppingParameters" -> MagneticTB`showHoppingParameters[2],
      "HoppingOccurrences" -> MagneticTB`showHoppingParameters[2, Global`t1],
      "NonHermitianHoppingParameters" ->
        MagneticTB`showHoppingParameters[2, "Hermitian" -> False],
      "BandDataRows" -> Take[Import[bandDataPath, "CSV"], 5],
      "BandDataPlot" -> ListPlot[
        Transpose[Import[bandDataPath, "CSV"]]
      ],
      "HRLines" -> Take[Import[hrPath, "Lines"], 12],
      "HROptionLines" -> Take[Import[hrOptionPath, "Lines"], 12],
      "HRBlochHamiltonian" ->
        MagneticTB`buildBlochHamiltonian[hrData, {0, 0, 0}],
      "HRDataOnsite" -> With[{
        onsiteData = MagneticTB`readHR[hrPath, "ncell" -> 0]
      },
        documentationKeyValueGrid@<|
          "Translations" -> onsiteData["Translations"],
          "HoppingMatrices" -> onsiteData["HoppingMatrices"]
        |>
      ],
      "BrokenSymmetryRules" -> brokenRules
    |>
  ]
];

initOptionFixture = checkedDocumentationEvaluation[
  "InitOptions",
  Module[{c4z},
    c4z = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {
        MagneticTB`a -> 1,
        MagneticTB`c -> 2
      },
      MagneticTB`wyckoffposition -> {
        {{0, 0, 0}, {0, 0, 0}}
      },
      MagneticTB`symminformation -> {
        {"C4z", c4z, {0, 0, 0}, "F"}
      },
      MagneticTB`basisFunctions -> {{"px", "py"}},
      MagneticTB`InitialBondShells -> 1,
      MagneticTB`GenerateSymmetryGroup -> True
    ];
    <|
      "GeneratedC4Representations" ->
        MagneticTB`showSymmetryRepresentations[All]
    |>
  ]
];

kagomeQNPDocumentationFixture = <|
  "BandPlot" -> checkedDocumentationEvaluation[
    "KagomeQNPInduced",
    ToExpression[kagomeQNPInducedInitCode, InputForm]
  ]
|>;

Print["Evaluating documentation fixture: ExplicitRepresentation (isolated)"];
explicitFixturePath = FileNameJoin[{
  scriptDirectory,
  ".InitFromRepDocumentationFixture.wl"
}];
If[FileExistsQ[explicitFixturePath], DeleteFile[explicitFixturePath]];
explicitFixtureProcess = RunProcess[
  {
    FileNameJoin[{
      $InstallationDirectory,
      "MacOS",
      "WolframKernel"
    }],
    "-script",
    FileNameJoin[{
      scriptDirectory,
      "ProbeInitFromRepDocumentation.wls"
    }]
  },
  All,
  ProcessDirectory -> projectRoot,
  ProcessEnvironment -> <|
    "OMP_NUM_THREADS" -> "1",
    "MKL_NUM_THREADS" -> "1"
  |>
];
If[
  explicitFixtureProcess["ExitCode"] =!= 0 ||
    explicitFixtureProcess["StandardError"] =!= "" ||
    !FileExistsQ[explicitFixturePath],
  Print[InputForm@<|
    "DocumentationFixture" -> "ExplicitRepresentation",
    "Process" -> explicitFixtureProcess
  |>];
  Exit[7]
];
explicitRepresentationFixture = Get[explicitFixturePath];
DeleteFile[explicitFixturePath];
If[
    !AssociationQ[explicitRepresentationFixture] ||
    !And @@ (KeyExistsQ[explicitRepresentationFixture, #] & /@ {
      "OrbitalTableData", "OrbitalTableGrid", "RepresentationMatrices",
      "RepresentationGrid", "BondGrid", "Onsite", "InducedOnsite",
      "InducedRepresentationGrid", "ModeDirectOrbitalTable",
      "ModeInducedOrbitalTable"
    }),
  Print["The isolated initfromrep documentation fixture is incomplete."];
  Exit[8]
];

ClearAll[documentationParameterRules];
documentationParameterRules[hamiltonian_, specified_List] := Module[
  {parameters, specifiedSymbols, defaults},
  parameters = Variables[hamiltonian];
  specifiedSymbols = First /@ specified;
  defaults = Thread[Complement[parameters, specifiedSymbols] -> 0];
  Join[specified, defaults]
];

advancedModelFixture = checkedDocumentationEvaluation[
  "GeneralExamples",
  Module[
    {operations, moS2Operations, shells, hamiltonian, path, rules, spinfulGraphene,
     spinfulGrapheneHRPath, moS2,
     moS2SpinOrbitA, moS2SpinOrbitB,
     moS2SpinRulesA, moS2SpinRulesB,
     c3, c3Expansion, c3Numeric,
     c3ChernNumber, nodal, c4, c4Raw},

    operations = MagneticTB`msgop[MagneticTB`gray[191]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 3},
      MagneticTB`wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"pzup", "pzdn"}},
      MagneticTB`InitialBondShells -> 3,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    hamiltonian = Total[MagneticTB`symham /@ Range[3]];
    path = {
      {{{0, 0, 0}, {0, 1/2, 0}}, {"G", "M"}},
      {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
      {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}}
    };
    rules = documentationParameterRules[hamiltonian, {
      Global`e1 -> 0, Global`r1 -> 0.02,
      Global`r2 -> 0, Global`t1 -> 0.5
    }];
    spinfulGrapheneHRPath = FileNameJoin[{
      $TemporaryDirectory,
      "magnetictb-documentation-spinful-graphene_hr.dat"
    }];
    If[FileExistsQ[spinfulGrapheneHRPath],
      DeleteFile[spinfulGrapheneHRPath]
    ];
    MagneticTB`hop[
      {1, 2, 3},
      rules,
      "hrExport" -> spinfulGrapheneHRPath
    ];
    spinfulGraphene = <|
      "Hamiltonian" -> hamiltonian,
      "BandPlot" -> MagneticTB`bandplot[path, 60, hamiltonian, rules],
      "HRLines" -> Take[Import[spinfulGrapheneHRPath, "Lines"], 12]
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[187]];
    moS2Operations = MapAt[
      Function[operation,
        FullSimplify[
          {{1, -1, 0}, {0, 1, 0}, {0, 0, 1}} . operation .
            Inverse[{{1, -1, 0}, {0, 1, 0}, {0, 0, 1}}]
        ]
      ],
      operations,
      {All, 2}
    ];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {1, 0, 0},
        {1/2, Sqrt[3]/2, 0},
        {0, 0, 100}
      },
      MagneticTB`lattpar -> {},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> moS2Operations,
      MagneticTB`basisFunctions -> {{"dz2", "dxy", "dx2-y2"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    moS2 = Total[MagneticTB`symham /@ {1, 2}] /. {
      Global`e1 -> Global`e2, Global`e2 -> Global`e1,
      Global`t1 -> Global`t6, Global`t2 -> Global`t5,
      Global`t3 -> Global`t4, Global`t4 -> Global`t3,
      Global`t5 -> Global`t2, Global`t6 -> Global`t1
    };
    path = {
      {{{0, 0, 0}, {2/3, 1/3, 0}}, {"G", "K"}},
      {{{2/3, 1/3, 0}, {1/2, 1/2, 0}}, {"K", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    rules = documentationParameterRules[moS2, {
      Global`e1 -> 1.046, Global`e2 -> 2.104,
      Global`t1 -> -0.184, Global`t2 -> 0.401,
      Global`t3 -> 0.218, Global`t4 -> 0.507,
      Global`t5 -> 0.338, Global`t6 -> 0.057
    }];
    moS2 = <|
      "BandPlot" -> MagneticTB`bandplot[path, 60, moS2, rules]
    |>;

    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 4},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{
        "dz2up", "dxyup", "dx2-y2up",
        "dz2dn", "dxydn", "dx2-y2dn"
      }},
      MagneticTB`InitialBondShells -> 1,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    moS2SpinRulesA = MagneticTB`brokenSymmetryInitRules[{9, 11, 13}];
    MagneticTB`init @@ moS2SpinRulesA;
    moS2SpinOrbitA = <|
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Hamiltonian" -> MagneticTB`symham[1]
    |>;

    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 4},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{
        "dz2up", "dz2dn", "dxyup",
        "dxydn", "dx2-y2up", "dx2-y2dn"
      }},
      MagneticTB`InitialBondShells -> 1,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    moS2SpinRulesB = MagneticTB`brokenSymmetryInitRules[{9, 11, 13}];
    MagneticTB`init @@ moS2SpinRulesB;
    moS2SpinOrbitB = <|
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Hamiltonian" -> MagneticTB`symham[1]
    |>;

    operations = MagneticTB`msgop[MagneticTB`bnsdict[{143, 3}]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {Sqrt[3] MagneticTB`a/2, -MagneticTB`a/2, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 1}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{{
        MagneticTB`x + I MagneticTB`y, 0
      }, {
        0, MagneticTB`x - I MagneticTB`y
      }}},
      MagneticTB`InitialBondShells -> 4,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    c3 = Total[MagneticTB`symham /@ {2, 3, 4}];
    c3 = c3[[{1, 4}, {1, 4}]];
    c3 = c3 /. {
      Global`r11 -> -Global`r2,
      Global`r15 -> -Global`r10,
      Global`r3 -> -Global`r6,
      Global`r7 -> -Global`r14,
      Global`s4 -> -Global`s4,
      Global`t4 -> Global`t1,
      Global`t8 -> Global`t9,
      Global`t11 -> -Global`t7
    };
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "M"}},
      {{{1/2, 0, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
      {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}},
      {{{0, 0, 0}, {0, 0, 1/2}}, {"G", "A"}},
      {{{0, 0, 1/2}, {1/2, 0, 1/2}}, {"A", "L"}},
      {{{1/2, 0, 1/2}, {1/3, 1/3, 1/2}}, {"L", "H"}},
      {{{1/3, 1/3, 1/2}, {0, 0, 1/2}}, {"H", "A"}},
      {{{1/2, 0, 1/2}, {1/2, 0, 0}}, {"L", "M"}},
      {{{1/3, 1/3, 0}, {1/3, 1/3, 1/2}}, {"K", "H"}}
    };
    rules = documentationParameterRules[c3, {
      Global`r10 -> -0.09, Global`r14 -> 0.25,
      Global`r2 -> 0, Global`r6 -> 0.5,
      Global`s4 -> -0.38, Global`s8 -> 0,
      Global`t1 -> 0,
      Global`t15 -> 0.735, Global`t7 -> 0.05,
      Global`t9 -> 0.015
    }];
    c3Expansion = Expand@Normal@Chop@Series[
      c3,
      {MagneticTB`ky, 0, 3},
      {MagneticTB`kx, 0, 3},
      {MagneticTB`kz, 0, 3}
    ];
    c3Expansion = FromCoefficientRules[
      Map[
        Select[#, Total[First[#]] <= 3 &] &,
        CoefficientRules[
          c3Expansion,
          {MagneticTB`kx, MagneticTB`ky, MagneticTB`kz}
        ],
        {2}
      ],
      {MagneticTB`kx, MagneticTB`ky, MagneticTB`kz}
    ];
    c3Numeric = Function[{momentum},
      c3Expansion /. rules /. {
        MagneticTB`kx -> momentum[[1]],
        MagneticTB`ky -> momentum[[2]],
        MagneticTB`kz -> momentum[[3]]
      }
    ];
    c3ChernNumber = MagneticTB`pointChernNumber[
      c3Numeric,
      1,
      {0, 0, 0},
      2,
      "SurfaceSubdivisions" -> 8
    ];
    c3 = <|
      "Hamiltonian" -> c3,
      "GammaAxisHamiltonian" -> FullSimplify[
        c3 /. {MagneticTB`kx -> 0, MagneticTB`ky -> 0},
        Element[MagneticTB`kz, Reals]
      ],
      "BandPlot" -> MagneticTB`bandplot[path, 200, c3, rules],
      "ChernNumber" -> c3ChernNumber
    |>;

    operations = MagneticTB`msgop[MagneticTB`bnsdict[{184, 196}]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {Sqrt[3] MagneticTB`a/2, -MagneticTB`a/2, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 1}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{{
        MagneticTB`x + I MagneticTB`y, 0
      }, {
        0, MagneticTB`x - I MagneticTB`y
      }}},
      MagneticTB`InitialBondShells -> 5,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    nodal = Total[MagneticTB`symham /@ Range[5]];
    nodal = nodal[[{1, 4}, {1, 4}]] /. {
      Global`p5n1 -> Global`p5n4,
      Global`p5n8 -> Global`p5n2,
      Global`s1 -> Global`s3,
      Global`t1 -> Global`t2,
      Global`p5n6 -> -Global`p5n7,
      Global`r2 -> -Global`r1
    };
    path = {
      {{{0, 1/10, 1/4}, {0, 0, 1/4}}, {"Q", "P"}},
      {{{0, 0, 1/4}, {-1/10, 1/15, 1/4}}, {"P", "Q'"}}
    };
    nodal = <|
      "Hamiltonian" -> nodal
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[75]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 1},
      MagneticTB`wyckoffposition -> {
        {{0, 0, 1/4}, {0, 0, 0}},
        {{0, 0, 0}, {0, 0, 0}}
      },
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"px", "py"}, {"px", "py"}},
      MagneticTB`InitialBondShells -> 7,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    c4Raw = Association@Table[
      shell -> MagneticTB`symham[shell],
      {shell, {2, 3, 4, 5, 7}}
    ];
    c4 = FullSimplify[
      (c4Raw[2] /. {Global`t1 -> 0, Global`t2 -> Global`t1p}) +
      (c4Raw[3] /. {Global`r1 -> 0, Global`r2 -> Global`tzp}) +
      (c4Raw[4] /. {
        Global`s10 -> 0, Global`s7 -> 0, Global`s9 -> 0,
        Global`s2 -> 0, Global`s4 -> 0, Global`s1 -> 0,
        Global`s6 -> 0, Global`s5 -> 0,
        Global`s8 -> Global`t1a, Global`s3 -> Global`t1b
      }) +
      ((c4Raw[5] /. {
        Global`p5n2 -> 0, Global`p5n3 -> 0,
        Global`p5n1 -> Global`p5n4
      }) /. Global`p5n4 -> Global`t2p) +
      ((c4Raw[7] /. {
        Global`p7n7 -> Global`p7n5, Global`p7n1 -> 0,
        Global`p7n4 -> 0, Global`p7n2 -> 0,
        Global`p7n3 -> 0, Global`p7n14 -> Global`p7n12,
        Global`p7n11 -> 0, Global`p7n8 -> 0,
        Global`p7n10 -> 0, Global`p7n9 -> 0,
        Global`p7n6 -> -Global`p7n5,
        Global`p7n13 -> Global`p7n12
      }) /. {
        Global`p7n5 -> Global`t2a/2,
        Global`p7n12 -> Global`t2b/2
      })
    ];
    path = {
      {{{0, 0, 0}, {1/2, 1/2, 0}}, {"G", "M"}},
      {{{1/2, 1/2, 0}, {1/2, 1/2, 1/2}}, {"M", "A"}},
      {{{1/2, 1/2, 1/2}, {0, 0, 1/2}}, {"A", "Z"}},
      {{{0, 0, 1/2}, {0, 0, 0}}, {"Z", "G"}}
    };
    rules = {
      Global`t1a -> 1, Global`t1b -> -1,
      Global`t2a -> 0.5, Global`t2b -> -0.5,
      Global`t1p -> 2.5, Global`t2p -> 0.5,
      Global`tzp -> 2
    };
    c4 = <|
      "Hamiltonian" -> c4,
      "ConventionII" -> MagneticTB`symhamII[c4],
      "BandPlot" -> MagneticTB`bandplot[path, 50, c4, rules]
    |>;

    <|
      "SpinfulGraphene" -> spinfulGraphene,
      "MoS2" -> moS2,
      "MoS2SpinOrbitA" -> moS2SpinOrbitA,
      "MoS2SpinOrbitB" -> moS2SpinOrbitB,
      "C3Weyl" -> c3,
      "NodalLine" -> nodal,
      "C4TI" -> c4,
      "WyckoffGrid" -> MagneticTB`showMSGWyckoff[MagneticTB`gray[193]],
      "WyckoffGridBNS187211" -> MagneticTB`showMSGWyckoff[
        MagneticTB`bnsdict[{187, 211}]
      ],
      "LayerOperations" -> Take[
        MagneticTB`mlgop[MagneticTB`graylayer[25]], 4
      ],
      "LayerOperationsNumbered" -> Take[
        MagneticTB`mlgop[{25, 3, 145}], 4
      ],
      "RodOperations" -> Take[
        MagneticTB`mrgop[MagneticTB`grayrod[25]], 4
      ],
      "RodOperationsNumbered" -> Take[
        MagneticTB`mrgop[{25, 3, 130}], 4
      ]
    |>
  ]
];

(* The archived backend is deliberately impossible to load together with the
   new backend.  Generate its visible results in a separate child kernel and
   import only the inert expressions used by the documentation pages. *)
legacyFixturePath = FileNameJoin[{
  scriptDirectory,
  ".LegacyDocumentationFixture.wl"
}];
If[FileExistsQ[legacyFixturePath], DeleteFile[legacyFixturePath]];
legacyFixtureProcess = RunProcess[
  {
    FileNameJoin[{
      $InstallationDirectory,
      "MacOS",
      "WolframKernel"
    }],
    "-script",
    FileNameJoin[{scriptDirectory, "ProbeLegacyDocumentation.wls"}]
  },
  All,
  ProcessDirectory -> projectRoot,
  ProcessEnvironment -> <|
    "OMP_NUM_THREADS" -> "1",
    "MKL_NUM_THREADS" -> "1"
  |>
];
If[
  legacyFixtureProcess["ExitCode"] =!= 0 ||
    legacyFixtureProcess["StandardError"] =!= "" ||
    !FileExistsQ[legacyFixturePath],
  Print[InputForm@<|
    "DocumentationFixture" -> "LegacyGraphene",
    "Process" -> legacyFixtureProcess
  |>];
  Exit[7]
];
legacyGrapheneFixture = Get[legacyFixturePath];
DeleteFile[legacyFixturePath];
If[
  !AssociationQ[legacyGrapheneFixture] ||
    !And @@ (KeyExistsQ[legacyGrapheneFixture, #] & /@ {
      "Onsite", "Nearest", "Hamiltonian", "BandUI", "BandPlot"
    }),
  Print["The isolated legacy documentation fixture is incomplete."];
  Exit[8]
];

validationModelFixture = checkedDocumentationEvaluation[
  "ValidatedApplications",
  Module[
    {
      operations, h, path, rules, session, simpleCubic, cscl, fese,
      moS2ElevenBand, spinSquare,
      cubicPDirect, cubicPInduced, collinear, noncollinear, pt, c3zt,
      identity3, c4z, c3xyz, z3, z32, omega,
      collinearSSG, noncollinearSSG, c3ztOperations
    },
    identity3 = IdentityMatrix[3];
    c4z = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
    c3xyz = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};

    operations = MagneticTB`msgop[MagneticTB`gray[221]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"s"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    h = Total[MagneticTB`symham /@ {1, 2}];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    simpleCubic = <|
      "Hamiltonian" -> h,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h],
      "BandPlot" -> MagneticTB`bandplot[
        path, 40, h, {Global`e1 -> 0.1, Global`t1 -> -0.4}
      ]
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[221]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {
        {{0, 0, 0}, {0, 0, 0}},
        {{1/2, 1/2, 1/2}, {0, 0, 0}}
      },
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"s"}, {"px", "py", "pz"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    h = Total[MagneticTB`symham /@ {1, 2}];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    rules = documentationParameterRules[h, {
      Global`e1 -> -0.3, Global`e2 -> 0.2, Global`t1 -> 0.5
    }];
    cscl = <|
      "Hamiltonian" -> h,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "CrossBlock" -> h[[1 ;; 1, 2 ;; 4]],
      "HoppingSources" -> MagneticTB`showHoppingParameters[2],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h],
      "BandPlot" -> MagneticTB`bandplot[path, 40, h, rules]
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[129]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      MagneticTB`wyckoffposition -> {
        {{3/4, 1/4, 0}, {0, 0, 0}},
        {{1/4, 1/4, 267/1000}, {0, 0, 0}}
      },
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {
        {"dxz", "dyz"},
        {"px", "py", "pz"}
      },
      MagneticTB`InitialBondShells -> 3,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    h = MagneticTB`symham[1] + MagneticTB`symham[3];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    rules = documentationParameterRules[h, {
      Global`e1 -> -0.4, Global`e2 -> 0.1,
      Global`e3 -> 0.5, Global`t1 -> 0.3
    }];
    fese = <|
      "Hamiltonian" -> h,
      "Block" -> h[[1 ;; 3, 1 ;; 3]],
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h],
      "BandPlot" -> MagneticTB`bandplot[path, 40, h, rules]
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[187]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 4},
      MagneticTB`wyckoffposition -> {
        {{0, 0, 0}, {0, 0, 0}},
        {{1/3, 2/3, 1/4}, {0, 0, 0}}
      },
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {
        {
          "dz2",
          "dx2-y2",
          "dxy", "dxz", "dyz"
        },
        {"px", "py", "pz"}
      },
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    h = Total[MagneticTB`symham /@ {1, 2}];
    path = {
      {{{0, 0, 0}, {2/3, 1/3, 0}}, {"G", "K"}},
      {{{2/3, 1/3, 0}, {1/2, 1/2, 0}}, {"K", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    moS2ElevenBand = <|
      "Hamiltonian" -> h,
      "OnsiteDBlock" -> h[[1 ;; 5, 1 ;; 5]],
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h]
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[123]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`c}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 4},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{
        "pxup", "pxdn", "pyup", "pydn"
      }},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    h = Total[MagneticTB`symham /@ {1, 2}];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    spinSquare = <|
      "Hamiltonian" -> h,
      "Onsite" -> MagneticTB`symham[1],
      "Representation" -> MagneticTB`showSymmetryRepresentations[],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h]
    |>;

    operations = MagneticTB`msgop[MagneticTB`gray[221]];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"px", "py", "pz"}},
      MagneticTB`RepresentationMode -> "DirectProduct"
    ];
    h = MagneticTB`symham[1];
    cubicPDirect = <|
      "Hamiltonian" -> h,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[]
    |>;
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{1/2, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> operations,
      MagneticTB`basisFunctions -> {{"px"}},
      MagneticTB`RepresentationMode -> "Induced"
    ];
    h = MagneticTB`symham[1];
    cubicPInduced = <|
      "Hamiltonian" -> h,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[]
    |>;

    collinearSSG = <|
      "C4spin" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {c4z, 0}
      |>,
      "tauHalfT" -> <|
        "space" -> {identity3, {1/2, 0, 0}},
        "spin" -> {identity3, 1}
      |>
    |>;
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 1}}},
      MagneticTB`symminformation -> collinearSSG,
      MagneticTB`basisFunctions -> {{"sup", "sdn"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> True
    ];
    h = Total[MagneticTB`symham /@ {1, 2}];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    collinear = <|
      "Hamiltonian" -> h,
      "AtomPositions" -> MagneticTB`atompos,
      "CrystalStructure" -> MagneticTB`showCrystalStructure[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h]
    |>;

    noncollinearSSG = <|
      "C4zspin" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {c4z, 0}
      |>,
      "C3xyzspin" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {c3xyz, 0}
      |>
    |>;
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> noncollinearSSG,
      MagneticTB`basisFunctions -> {{"sup", "sdn"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> True
    ];
    h = Total[MagneticTB`symham /@ {1, 2}];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    noncollinear = <|
      "Hamiltonian" -> h,
      "Representation" -> MagneticTB`showSymmetryRepresentations[],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h]
    |>;

    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{1/4, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> {
        {"E", identity3, {0, 0, 0}, "F"},
        {"PT", -identity3, {0, 0, 0}, "T"}
      },
      MagneticTB`basisFunctions -> {{"s"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False,
      MagneticTB`RepresentationMode -> "Induced"
    ];
    h = MagneticTB`symham[2];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    pt = <|
      "Hamiltonian" -> h,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[All],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h]
    |>;

    z3 = {
      {-1/2, -Sqrt[3]/2, 0},
      {Sqrt[3]/2, -1/2, 0},
      {0, 0, 1}
    };
    z32 = MatrixPower[z3, 2];
    c3ztOperations = <|
      "E" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {identity3, 0}
      |>,
      "C3" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {z3, 0}
      |>,
      "C3^2" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {z32, 0}
      |>,
      "T" -> <|
        "space" -> {-identity3, {0, 0, 0}},
        "spin" -> {identity3, 1}
      |>,
      "C3T" -> <|
        "space" -> {-identity3, {0, 0, 0}},
        "spin" -> {z3, 1}
      |>,
      "C3^2T" -> <|
        "space" -> {-identity3, {0, 0, 0}},
        "spin" -> {z32, 1}
      |>
    |>;
    omega = Exp[2 Pi I/3];
    MagneticTB`initfromrep[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{1/4, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> c3ztOperations,
      MagneticTB`repinformation -> {
        <|
          "ReferenceSiteIndex" -> 1,
          "SiteSymmetryOperationIndices" -> {1, 2, 3},
          "SiteSymmetryMatrices" -> {
            {{1}}, {{omega}}, {{omega^2}}
          }
        |>
      },
      MagneticTB`orbitalLabels -> {{"omega"}},
      MagneticTB`RepresentationMode -> "Induced",
      MagneticTB`InitialBondShells -> 5
    ];
    h = Sum[MagneticTB`symham[i], {i, 5}];
    path = {
      {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
      {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}},
      {{{1/2, 1/2, 0}, {0, 0, 0}}, {"M", "G"}}
    };
    c3zt = <|
      "Hamiltonian" -> h,
      "Representation" -> MagneticTB`showSymmetryRepresentations[All],
      "BandUI" -> MagneticTB`bandManipulate[path, 30, h]
    |>;

    <|
      "SimpleCubicS" -> simpleCubic,
      "CsCl" -> cscl,
      "FeSe" -> fese,
      "MoS2ElevenBand" -> moS2ElevenBand,
      "SpinfulSquare" -> spinSquare,
      "CubicPDirect" -> cubicPDirect,
      "CubicPInduced" -> cubicPInduced,
      "CollinearSSG" -> collinear,
      "NoncollinearSSG" -> noncollinear,
      "EPT" -> pt,
      "C3Z2T" -> c3zt
    |>
  ]
];

spinLayerDocumentationFixture = Get[FileNameJoin[{
  scriptDirectory,
  "Examples",
  "SpinLayer14123L1Fixture.wl"
}]];
If[!AssociationQ[spinLayerDocumentationFixture],
  Print["The bundled 14.1.2.3.L.1 documentation fixture is invalid."];
  Exit[9]
];

continuousModelFixture = checkedDocumentationEvaluation[
  "ContinuousCInfinityModels",
  Module[
    {
      identity3, theta, continuousRotation, continuousInformation,
      continuousMatrix, h, simple, slgData, finiteRecords,
      modelLattice, seedPosition, finiteInformation, generatedSites,
      stabilizerIndices, corepResult, corepTable, finiteSpinorMatrices,
      siteSymmetryMatrices, fullInformation, direct, induced,
      directContinuousMatrices, inducedContinuousMatrices,
      directHamiltonian, inducedHamiltonian, path
    },
    identity3 = IdentityMatrix[3];
    theta = Global`theta;
    continuousRotation = {
      {Cos[theta], -Sin[theta], 0},
      {Sin[theta], Cos[theta], 0},
      {0, 0, 1}
    };
    continuousInformation = <|
      "E" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {identity3, 0}
      |>,
      "C_infty" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {continuousRotation, 0},
        "continuous" -> theta
      |>
    |>;
    continuousMatrix = DiagonalMatrix[{
      Exp[-I theta/2], Exp[I theta/2]
    }];
    MagneticTB`initfromrep[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> continuousInformation,
      MagneticTB`repinformation -> <|
        "Discrete" -> {{identity3[[1 ;; 2, 1 ;; 2]]}},
        "Continuous" -> <|
          "SiteMatrices" -> {{continuousMatrix}}
        |>
      |>,
      MagneticTB`orbitalLabels -> {{"m=+1/2", "m=-1/2"}},
      MagneticTB`InitialBondShells -> 2
    ];
    h = MagneticTB`symham[1];
    simple = <|
      "Hamiltonian" -> h,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[]
    |>;

    slgData = spinLayerDocumentationFixture;
    finiteRecords = slgData["FiniteRepresentatives"];
    modelLattice = DiagonalMatrix[{1, 6/5, 4}];
    seedPosition = {0, 0, 0};
    finiteInformation = AssociationThread[
      "g" <> ToString[#] & /@ Range[Length[finiteRecords]],
      Map[
        Function[record,
          <|
            "space" -> {
              record["SpaceMatrix"], record["Translation"]
            },
            "spin" -> {
              record["SpinMatrix"],
              Boole[TrueQ[record["AntiunitaryQ"]]]
            }
          |>
        ],
        finiteRecords
      ]
    ];
    generatedSites = DeleteDuplicates[
      Mod[
        #["SpaceMatrix"].seedPosition + #["Translation"],
        1
      ] & /@ finiteRecords
    ];
    stabilizerIndices = Select[
      Range[Length[finiteRecords]],
      Mod[
        finiteRecords[[#, "SpaceMatrix"]].seedPosition +
          finiteRecords[[#, "Translation"]],
        1
      ] === seedPosition &
    ];
    finiteSpinorMatrices = slgData["FiniteSpinorMatrices"];
    siteSymmetryMatrices = (
      {{#[[1, 1]]}} & /@
        finiteSpinorMatrices[[stabilizerIndices]]
    );
    fullInformation = Append[
      finiteInformation,
      "C_infty" -> <|
        "space" -> {identity3, {0, 0, 0}},
        "spin" -> {continuousRotation, 0},
        "continuous" -> theta
      |>
    ];
    path = {
      {{{0, 0, 0}, {1/2, 1/2, 0}}, {"G", "S"}},
      {{{1/2, 1/2, 0}, {1/2, 0, 0}}, {"S", "X"}},
      {{{1/2, 0, 0}, {0, 0, 0}}, {"X", "G"}}
    };

    directContinuousMatrices = ConstantArray[
      continuousMatrix, Length[generatedSites]
    ];
    MagneticTB`initfromrep[
      MagneticTB`lattice -> modelLattice,
      MagneticTB`lattpar -> {},
      MagneticTB`wyckoffposition -> {{seedPosition, {0, 0, 1}}},
      MagneticTB`symminformation -> fullInformation,
      MagneticTB`repinformation -> <|
        "Discrete" -> {finiteSpinorMatrices},
        "Continuous" -> <|
          "SiteMatrices" -> {directContinuousMatrices}
        |>
      |>,
      MagneticTB`orbitalLabels -> {{"m=+1/2", "m=-1/2"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`RepresentationMode -> "DirectProduct"
    ];
    directHamiltonian = Sum[MagneticTB`symham[i], {i, 2}];
    direct = <|
      "Hamiltonian" -> directHamiltonian,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[],
      "BandUI" -> MagneticTB`bandManipulate[path, 40, directHamiltonian]
    |>;

    inducedContinuousMatrices = {
      {{Exp[-I theta/2]}},
      {{Exp[I theta/2]}}
    };
    MagneticTB`initfromrep[
      MagneticTB`lattice -> modelLattice,
      MagneticTB`lattpar -> {},
      MagneticTB`wyckoffposition -> {{seedPosition, {0, 0, 1}}},
      MagneticTB`symminformation -> fullInformation,
      MagneticTB`repinformation -> <|
        "Discrete" -> {
          <|
            "ReferenceSiteIndex" -> 1,
            "SiteSymmetryOperationIndices" -> stabilizerIndices,
            "SiteSymmetryMatrices" -> siteSymmetryMatrices
          |>
        },
        "Continuous" -> <|
          "SiteMatrices" -> {inducedContinuousMatrices}
        |>
      |>,
      MagneticTB`orbitalLabels -> {{"m=+1/2"}},
      MagneticTB`InitialBondShells -> 5,
      MagneticTB`RepresentationMode -> "Induced"
    ];
    inducedHamiltonian = Sum[MagneticTB`symham[i], {i, 5}];
    induced = <|
      "Hamiltonian" -> inducedHamiltonian,
      "OrbitalTable" -> MagneticTB`orbitalTable[],
      "Representation" -> MagneticTB`showSymmetryRepresentations[],
      "HoppingSources" -> MagneticTB`showHoppingParameters[5],
      "BandUI" -> MagneticTB`bandManipulate[path, 40, inducedHamiltonian],
      "BandPlot" -> Plot[
        Evaluate@Eigenvalues[
          inducedHamiltonian /. {
            Global`e1 -> 0, Global`r1 -> 0, Global`s1 -> 0,
            Global`p5n1 -> 1, Global`p5n2 -> 0,
            MagneticTB`kx -> Pi Global`q,
            MagneticTB`ky -> Pi Global`q,
            MagneticTB`kz -> 0
          }
        ],
        {Global`q, 0, 1},
        Frame -> True,
        FrameLabel -> {"G-S", "Energy"},
        PlotRange -> All
      ]
    |>;

    <|
      "Simple" -> simple,
      "SLGSummary" -> <|
        "Name" -> "14.1.2.3.L.1",
        "MagneticType" -> slgData["MagneticType"],
        "FiniteRepresentativeCount" -> Length[finiteRecords],
        "ContinuousSpinPart" -> slgData["ContinuousSpinPart"],
        "GeneratedSites" -> generatedSites,
        "SiteStabilizerIndices" -> stabilizerIndices
      |>,
      "DirectProduct" -> direct,
      "Induced" -> induced
    |>
  ]
];

compilerInputCode = StringRiffle[{
  "compilerInput = <|",
  "  \"Lattice\" -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, a}",
  "  },",
  "  \"LatticeParameters\" -> {a -> 1},",
  "  \"WyckoffPosition\" -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  \"SymmetryInformation\" -> {",
  "    {\"E\", IdentityMatrix[3], {0, 0, 0}, \"F\"}",
  "  },",
  "  \"BasisFunctions\" -> {{\"s\"}},",
  "  \"Debug\" -> False,",
  "  \"InitialBondShells\" -> 2,",
  "  \"GenerateSymmetryGroup\" -> False,",
  "  \"RepresentationMode\" -> \"DirectProduct\",",
  "  \"SiteLocalData\" -> Automatic",
  "|>;"
}, "\n"];

advancedUtilityFixture = checkedDocumentationEvaluation[
  "AdvancedUtilities",
  Module[{compiled, unconstrained, constrained, pointOperations, pointMatrices,
    symmetrizationOperations, symmetrizationData,
    shiftedSymmetrizationData, unitarySymmetrizationData,
    scalarSymmetrizationData},
    compiled = MagneticTB`CompileMagneticTBInput[<|
      "Lattice" -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      "LatticeParameters" -> {MagneticTB`a -> 1},
      "WyckoffPosition" -> {{{0, 0, 0}, {0, 0, 0}}},
      "SymmetryInformation" -> MagneticTB`msgop[MagneticTB`gray[221]],
      "BasisFunctions" -> {{"s"}},
      "Debug" -> False,
      "InitialBondShells" -> 2,
      "GenerateSymmetryGroup" -> False,
      "RepresentationMode" -> "DirectProduct",
      "SiteLocalData" -> Automatic
    |>];
    MagneticTB`init[
      MagneticTB`lattice -> {
        {MagneticTB`a, 0, 0},
        {0, MagneticTB`a, 0},
        {0, 0, MagneticTB`a}
      },
      MagneticTB`lattpar -> {MagneticTB`a -> 1},
      MagneticTB`wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},
      MagneticTB`symminformation -> MagneticTB`msgop[MagneticTB`gray[221]],
      MagneticTB`basisFunctions -> {{"s"}},
      MagneticTB`InitialBondShells -> 2,
      MagneticTB`GenerateSymmetryGroup -> False
    ];
    unconstrained = MagneticTB`unsymham[2];
    constrained = MagneticTB`symham[2];
    pointOperations = Select[
      MagneticTB`msgop[MagneticTB`gray[75]],
      Last[#] === "F" &
    ];
    pointMatrices = MagneticTB`pointMatrix[
      pointOperations, {"px", "py"}, IdentityMatrix[3]
    ];
    symmetrizationOperations = MagneticTB`msgop[MagneticTB`gray[75]];
    symmetrizationData = MagneticTB`symmetrizationHRInit[
      "Lattice" -> DiagonalMatrix[{MagneticTB`a, MagneticTB`a, MagneticTB`c}],
      "LattPar" -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      "WyckoffPosition" -> {{{0, 0, 0}, {0, 0, 0}}},
      "Symminformation" -> symmetrizationOperations,
      "BasisFunctions" -> {"px", "py"},
      "Software" -> "VASP"
    ];
    shiftedSymmetrizationData = MagneticTB`symmetrizationHRInit[
      "Lattice" -> DiagonalMatrix[{MagneticTB`a, MagneticTB`a, MagneticTB`c}],
      "LattPar" -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      "WyckoffPosition" -> {{{1/2, 1/2, 0}, {0, 0, 0}}},
      "Symminformation" -> symmetrizationOperations,
      "BasisFunctions" -> {"px", "py"},
      "Software" -> "VASP"
    ];
    unitarySymmetrizationData = MagneticTB`symmetrizationHRInit[
      "Lattice" -> DiagonalMatrix[{MagneticTB`a, MagneticTB`a, MagneticTB`c}],
      "LattPar" -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      "WyckoffPosition" -> {{{0, 0, 0}, {0, 0, 0}}},
      "Symminformation" -> Select[
        symmetrizationOperations,
        Last[#] === "F" &
      ],
      "BasisFunctions" -> {"px", "py"},
      "Software" -> "VASP"
    ];
    scalarSymmetrizationData = MagneticTB`symmetrizationHRInit[
      "Lattice" -> DiagonalMatrix[{MagneticTB`a, MagneticTB`a, MagneticTB`c}],
      "LattPar" -> {MagneticTB`a -> 1, MagneticTB`c -> 2},
      "WyckoffPosition" -> {{{0, 0, 0}, {0, 0, 0}}},
      "Symminformation" -> symmetrizationOperations,
      "BasisFunctions" -> {"s", "pz"},
      "Software" -> "VASP"
    ];
    <|
      "CompilerSummary" -> documentationKeyValueGrid@<|
        "Lattice" -> compiled["ModelSpecification", "Lattice"],
        "SiteOrbits" -> compiled["ModelSpecification", "SiteOrbits"],
        "RepresentationMatrices" ->
          compiled["FullRepresentation", "RepresentationMatrices"],
        "FirstBondShell" ->
          compiled["ModelSpecification", "BondClasses"][[1]]
      |>,
      "Unconstrained" -> unconstrained,
      "Constrained" -> constrained,
      "PointMatrices" -> pointMatrices,
      "SymmetrizationData" -> symmetrizationData,
      "ShiftedSymmetrizationData" -> shiftedSymmetrizationData,
      "UnitarySymmetrizationData" -> unitarySymmetrizationData,
      "ScalarSymmetrizationData" -> scalarSymmetrizationData,
      "C4Group" -> MagneticTB`GenerateGroup[
        {1}, 0, Mod[#1 + #2, 4] &
      ],
      "C4Generators" -> MagneticTB`getGenerator[
        Range[0, 3], 0, Mod[#1 + #2, 4] &
      ]
    |>
  ]
];

simpleCubicUtilityInitCode = StringRiffle[{
  "cubicUtilityOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, a}",
  "  },",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> cubicUtilityOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];"
}, "\n"];

cyclotomicIRootCode = StringRiffle[{
  "iRoot = <|",
  "  \"kind\" -> \"root_sum\",",
  "  \"constant\" -> 0,",
  "  \"terms\" -> {",
  "    <|",
  "      \"coefficient\" -> 1,",
  "      \"order\" -> 4,",
  "      \"power\" -> 1",
  "    |>",
  "  }",
  "|>;"
}, "\n"];

cyclotomicStructuredCompileCode = StringRiffle[{
  cyclotomicIRootCode,
  "compiledCyclotomic = CompileCyclotomicMatrices[<|",
  "  \"coordinate_dimension\" -> 3,",
  "  \"constraints\" -> {",
  "    {{1, iRoot, 0}},",
  "    {{0, 1, 0}}",
  "  },",
  "  \"max_cyclotomic_degree\" -> 128",
  "|>];"
}, "\n"];

cyclotomicRationalCompileCode = StringRiffle[{
  "compiledRational = CompileCyclotomicMatrices[<|",
  "  \"coordinate_dimension\" -> 3,",
  "  \"constraints\" -> {",
  "    {",
  "      {1, 2, 3},",
  "      {2, 4, 6}",
  "    }",
  "  }",
  "|>];",
  "exactRationalMatrix = compiledRational[\"constraints\"][[1]];"
}, "\n"];

cyclotomicCommonCompileCode = StringRiffle[{
  "compiledCommon = CompileCyclotomicMatrices[<|",
  "  \"coordinate_dimension\" -> 3,",
  "  \"constraints\" -> {",
  "    {{1, 0, 0}},",
  "    {{0, 1, 0}}",
  "  }",
  "|>];"
}, "\n"];

functionSpecifications = {
  <|
    "Name" -> "CompileCyclotomicMatrices",
    "Usage" -> {{
      "CompileCyclotomicMatrices[input]",
      "puts a collection of exact constraint matrices into one common cyclotomic field for exact row reduction and null-space calculations."
    }},
    "Notes" -> {
      StringJoin[
        "Use this advanced entry when a program already stores roots of unity in the documented root_sum ",
        "form. For ordinary exact Wolfram Language matrices containing rationals, I, radicals, or ",
        "rational-angle roots of unity, CyclotomicCommonNullSpace is the simpler entry."
      ],
      StringJoin[
        "input is an Association with \"coordinate_dimension\", \"constraints\", and optional ",
        "\"max_cyclotomic_degree\" keys. Every constraint has the same number of columns. A constraint may ",
        "be a rectangular nested list or an Association with \"rows\", \"columns\", and flat \"entries\"."
      ],
      StringJoin[
        "A root_sum scalar has an exact rational constant and a list of terms with exact rational ",
        "\"coefficient\", positive integer \"order\", and integer \"power\". The compiler takes the least ",
        "common multiple of all root orders and stores every matrix in that one exact field."
      ],
      StringJoin[
        "The result contains exact internal matrix objects for CyclotomicRREF, CyclotomicNullSpace, and ",
        "CyclotomicCommonKernel. Use RestoreCyclotomicExpression before presenting their entries as ordinary ",
        "Wolfram Language expressions."
      ]
    },
    "Examples" -> {{
      "Compile two three-coordinate constraints containing the fourth root of unity i. The summary shows that both blocks share conductor 4:",
      cyclotomicStructuredCompileCode <>
        "\nKeyTake[\n" <>
        "  compiledCyclotomic,\n" <>
        "  {\"type\", \"conductor\", \"coordinate_dimension\", \"constraint_count\"}\n" <>
        "]",
      cyclotomicDocumentationFixture["StructuredSummary"]
    }},
    "Applications" -> {{
      "Restore the two compiled blocks to see the exact matrices represented by the shared field:",
      cyclotomicStructuredCompileCode <>
        "\nrestoredBlocks = RestoreCyclotomicExpression /@\n" <>
        "  compiledCyclotomic[\"constraints\"];\n" <>
        "MatrixForm[\n" <>
        "  Partition[#[\"entries\"], #[\"columns\"]]\n" <>
        "] & /@ restoredBlocks",
      MatrixForm /@ cyclotomicDocumentationFixture["StructuredMatrices"]
    }},
    "Issues" -> {
      {
        "Every matrix must have the declared coordinate dimension. A row with only two entries in a three-coordinate problem fails explicitly:",
        "dimensionFailure = CompileCyclotomicMatrices[<|\n" <>
          "  \"coordinate_dimension\" -> 3,\n" <>
          "  \"constraints\" -> {{{1, 2}}}\n" <>
          "|>];\n" <>
          "<|\"FailureTag\" -> dimensionFailure[[1]]|>",
        <|"FailureTag" -> cyclotomicDocumentationFixture["DimensionFailure"]|>
      },
      {
        "The field degree may not exceed max_cyclotomic_degree. This boundary stops the calculation instead of changing to approximate arithmetic:",
        "z257 = <|\n" <>
          "  \"kind\" -> \"root_sum\",\n" <>
          "  \"constant\" -> 0,\n" <>
          "  \"terms\" -> {\n" <>
          "    <|\"coefficient\" -> 1, \"order\" -> 257, \"power\" -> 1|>\n" <>
          "  }\n" <>
          "|>;\n" <>
          "degreeFailure = CompileCyclotomicMatrices[<|\n" <>
          "  \"coordinate_dimension\" -> 1,\n" <>
          "  \"constraints\" -> {{{z257}}},\n" <>
          "  \"max_cyclotomic_degree\" -> 128\n" <>
          "|>];\n" <>
          "<|\"FailureTag\" -> degreeFailure[[1]]|>",
        <|"FailureTag" -> cyclotomicDocumentationFixture["LowDegreeFailure"]|>
      }
    },
    "SeeAlso" -> {
      {"CyclotomicRREF", "MagneticTB/ref/CyclotomicRREF"},
      {"CyclotomicNullSpace", "MagneticTB/ref/CyclotomicNullSpace"},
      {"CyclotomicCommonKernel", "MagneticTB/ref/CyclotomicCommonKernel"},
      {"CyclotomicCommonNullSpace", "MagneticTB/ref/CyclotomicCommonNullSpace"},
      {"RestoreCyclotomicExpression", "MagneticTB/ref/RestoreCyclotomicExpression"}
    },
    "Keywords" -> {"cyclotomic field", "exact constraints", "root of unity"}
  |>,
  <|
    "Name" -> "CyclotomicRREF",
    "Usage" -> {{
      "CyclotomicRREF[matrix]",
      "computes a deterministic exact reduced row echelon form for one matrix returned by CompileCyclotomicMatrices."
    }},
    "Notes" -> {
      StringJoin[
        "matrix must be an exact matrix object from CompileCyclotomicMatrices. The function returns ",
        "\"reduced_matrix\", \"pivot_columns\", \"free_columns\", \"rank\", and \"nullity\"."
      ],
      StringJoin[
        "pivot_columns and free_columns use 0-based indices so that the result has the same deterministic ",
        "convention as the cross-language fixture format. Add one when using an index with Wolfram Language ",
        "Part."
      ],
      StringJoin[
        "The production calculation uses MagneticTB's exact cyclotomic row-reduction algorithm. It does not ",
        "call the built-in RowReduce and does not switch to floating point, SVD, or another fallback."
      ]
    },
    "Examples" -> {{
      "Compile a rank-one rational matrix and display its exact RREF together with its 0-based pivot and free columns:",
      cyclotomicRationalCompileCode <>
        "\nrrefResult = CyclotomicRREF[exactRationalMatrix];\n" <>
        "restoredRREF = RestoreCyclotomicExpression[\n" <>
        "  rrefResult[\"reduced_matrix\"]\n" <>
        "];\n" <>
        "{\n" <>
        "  MatrixForm@Partition[\n" <>
        "    restoredRREF[\"entries\"], restoredRREF[\"columns\"]\n" <>
        "  ],\n" <>
        "  KeyTake[\n" <>
        "    rrefResult,\n" <>
        "    {\"pivot_columns\", \"free_columns\", \"rank\", \"nullity\"}\n" <>
        "  ]\n" <>
        "}",
      cyclotomicDocumentationFixture["RREFSummary"]
    }},
    "Applications" -> {{
      "The same call works in a nontrivial cyclotomic field. Here the matrix entry i is kept exact throughout the reduction:",
      cyclotomicStructuredCompileCode <>
        "\ncyclotomicRREF = CyclotomicRREF[\n" <>
        "  compiledCyclotomic[\"constraints\"][[1]]\n" <>
        "];\n" <>
        "restoredCyclotomicRREF = RestoreCyclotomicExpression[\n" <>
        "  cyclotomicRREF[\"reduced_matrix\"]\n" <>
        "];\n" <>
        "{\n" <>
        "  MatrixForm@Partition[\n" <>
        "    restoredCyclotomicRREF[\"entries\"],\n" <>
        "    restoredCyclotomicRREF[\"columns\"]\n" <>
        "  ],\n" <>
        "  KeyTake[\n" <>
        "    cyclotomicRREF,\n" <>
        "    {\"pivot_columns\", \"free_columns\", \"rank\", \"nullity\"}\n" <>
        "  ]\n" <>
        "}",
      cyclotomicDocumentationFixture["CyclotomicRREFSummary"]
    }},
    "Issues" -> {{
      "An ordinary nested list is not a compiled exact matrix object. Use CompileCyclotomicMatrices first, or use CyclotomicCommonNullSpace for ordinary exact matrices:",
      "rrefFailure = CyclotomicRREF[{{1, 2, 3}}];\n" <>
        "<|\"FailureTag\" -> rrefFailure[[1]]|>",
      <|"FailureTag" -> cyclotomicDocumentationFixture["MalformedFailure"]|>
    }},
    "SeeAlso" -> {
      {"CompileCyclotomicMatrices", "MagneticTB/ref/CompileCyclotomicMatrices"},
      {"CyclotomicNullSpace", "MagneticTB/ref/CyclotomicNullSpace"},
      {"RestoreCyclotomicExpression", "MagneticTB/ref/RestoreCyclotomicExpression"}
    },
    "Keywords" -> {"exact RREF", "pivot columns", "cyclotomic linear algebra"}
  |>,
  <|
    "Name" -> "CyclotomicNullSpace",
    "Usage" -> {{
      "CyclotomicNullSpace[matrix]",
      "finds a deterministic exact null-space basis for one matrix returned by CompileCyclotomicMatrices."
    }},
    "Notes" -> {
      StringJoin[
        "The result contains both \"nullspace_rows\" and \"basis_matrix\". nullspace_rows has one basis ",
        "vector per row; basis_matrix is its transpose, so matrix . basis_matrix is exactly zero."
      ],
      StringJoin[
        "The result also gives the exact reduced matrix, 0-based pivot_columns and free_columns, rank, ",
        "nullity, and exact_residual_verified. The row basis is ordered by increasing free-column index."
      ],
      StringJoin[
        "The production algorithm does not call the built-in NullSpace or RowReduce. A failure is returned ",
        "directly; the function never changes to floating point, SVD, or another null-space algorithm."
      ]
    },
    "Examples" -> {{
      "For this rank-one matrix, the two row-basis vectors are {-2,1,0} and {-3,0,1}. Their transpose is the column-basis matrix annihilated by the input:",
      cyclotomicRationalCompileCode <>
        "\nkernelResult = CyclotomicNullSpace[exactRationalMatrix];\n" <>
        "restoredRows = RestoreCyclotomicExpression[\n" <>
        "  kernelResult[\"nullspace_rows\"]\n" <>
        "];\n" <>
        "restoredBasis = RestoreCyclotomicExpression[\n" <>
        "  kernelResult[\"basis_matrix\"]\n" <>
        "];\n" <>
        "{\n" <>
        "  MatrixForm@Partition[\n" <>
        "    restoredRows[\"entries\"], restoredRows[\"columns\"]\n" <>
        "  ],\n" <>
        "  MatrixForm@Partition[\n" <>
        "    restoredBasis[\"entries\"], restoredBasis[\"columns\"]\n" <>
        "  ],\n" <>
        "  KeyTake[\n" <>
        "    kernelResult,\n" <>
        "    {\"pivot_columns\", \"free_columns\", \"rank\", \"nullity\",\n" <>
        "      \"exact_residual_verified\"}\n" <>
        "  ]\n" <>
        "}",
      cyclotomicDocumentationFixture["KernelSummary"]
    }},
    "Properties" -> {{
      "Full column rank gives a zero-dimensional null space. The basis matrix retains the exact 2 by 0 shape even though it has no entries:",
      "fullRankCompiled = CompileCyclotomicMatrices[<|\n" <>
        "  \"coordinate_dimension\" -> 2,\n" <>
        "  \"constraints\" -> {IdentityMatrix[2]}\n" <>
        "|>];\n" <>
        "fullRankKernel = CyclotomicNullSpace[\n" <>
        "  fullRankCompiled[\"constraints\"][[1]]\n" <>
        "];\n" <>
        "{\n" <>
        "  KeyTake[\n" <>
        "    fullRankKernel,\n" <>
        "    {\"rank\", \"nullity\", \"exact_residual_verified\"}\n" <>
        "  ],\n" <>
        "  KeyTake[\n" <>
        "    RestoreCyclotomicExpression[fullRankKernel[\"basis_matrix\"]],\n" <>
        "    {\"rows\", \"columns\", \"entries\"}\n" <>
        "  ]\n" <>
        "}",
      cyclotomicDocumentationFixture["FullRankSummary"]
    }},
    "Issues" -> {{
      "CyclotomicNullSpace accepts a compiled exact matrix, not a raw nested list. The direct ordinary-matrix entry is CyclotomicCommonNullSpace:",
      "kernelFailure = CyclotomicNullSpace[{{1, 2, 3}}];\n" <>
        "<|\"FailureTag\" -> kernelFailure[[1]]|>",
      <|"FailureTag" -> cyclotomicDocumentationFixture["MalformedFailure"]|>
    }},
    "SeeAlso" -> {
      {"CompileCyclotomicMatrices", "MagneticTB/ref/CompileCyclotomicMatrices"},
      {"CyclotomicRREF", "MagneticTB/ref/CyclotomicRREF"},
      {"CyclotomicCommonKernel", "MagneticTB/ref/CyclotomicCommonKernel"},
      {"CyclotomicCommonNullSpace", "MagneticTB/ref/CyclotomicCommonNullSpace"},
      {"RestoreCyclotomicExpression", "MagneticTB/ref/RestoreCyclotomicExpression"}
    },
    "Keywords" -> {"exact null space", "row basis", "column basis"}
  |>,
  <|
    "Name" -> "CyclotomicCommonKernel",
    "Usage" -> {{
      "CyclotomicCommonKernel[compiled]",
      "finds the vectors annihilated by every exact constraint block in a compiled cyclotomic problem."
    }},
    "Notes" -> {
      "compiled must be the complete Association returned by CompileCyclotomicMatrices. All blocks have the same coordinate dimension and exact field context.",
      StringJoin[
        "The result stores the common basis by columns in \"basis_matrix\" and by rows in \"nullspace_rows\".",
        " iteration_nullities begins with the coordinate dimension and records the remaining nullity after ",
        "each input block."
      ],
      StringJoin[
        "The method is the production iterative common-kernel algorithm. Every block is checked against the ",
        "final exact basis. If any step or residual check fails, the function returns Failure; it never ",
        "stacks all blocks as an automatic fallback and never uses floating point or SVD."
      ]
    },
    "Examples" -> {{
      "The first constraint removes the first coordinate and the second removes the second. Their common kernel is therefore spanned by {0,0,1}:",
      cyclotomicCommonCompileCode <>
        "\ncommonResult = CyclotomicCommonKernel[compiledCommon];\n" <>
        "commonRows = RestoreCyclotomicExpression[\n" <>
        "  commonResult[\"nullspace_rows\"]\n" <>
        "];\n" <>
        "commonBasis = RestoreCyclotomicExpression[\n" <>
        "  commonResult[\"basis_matrix\"]\n" <>
        "];\n" <>
        "{\n" <>
        "  MatrixForm@Partition[\n" <>
        "    commonRows[\"entries\"], commonRows[\"columns\"]\n" <>
        "  ],\n" <>
        "  MatrixForm@Partition[\n" <>
        "    commonBasis[\"entries\"], commonBasis[\"columns\"]\n" <>
        "  ],\n" <>
        "  KeyTake[\n" <>
        "    commonResult,\n" <>
        "    {\"iteration_nullities\", \"rank\", \"nullity\",\n" <>
        "      \"exact_residual_verified\"}\n" <>
        "  ]\n" <>
        "}",
      cyclotomicDocumentationFixture["CommonSummary"]
    }},
    "Properties" -> {
      {
        "With no constraint blocks, every coordinate is free. The row basis is the identity and iteration_nullities contains only the initial dimension:",
        "emptyCompiled = CompileCyclotomicMatrices[<|\n" <>
          "  \"coordinate_dimension\" -> 3,\n" <>
          "  \"constraints\" -> {}\n" <>
          "|>];\n" <>
          "emptyResult = CyclotomicCommonKernel[emptyCompiled];\n" <>
          "emptyRows = RestoreCyclotomicExpression[\n" <>
          "  emptyResult[\"nullspace_rows\"]\n" <>
          "];\n" <>
          "{\n" <>
          "  MatrixForm@Partition[\n" <>
          "    emptyRows[\"entries\"], emptyRows[\"columns\"]\n" <>
          "  ],\n" <>
          "  KeyTake[emptyResult, {\"iteration_nullities\", \"rank\", \"nullity\"}]\n" <>
          "}",
        cyclotomicDocumentationFixture["EmptyCommonSummary"]
      },
      {
        "A problem with coordinate dimension zero has rank and nullity zero and retains a 0 by 0 basis matrix:",
        "zeroCompiled = CompileCyclotomicMatrices[<|\n" <>
          "  \"coordinate_dimension\" -> 0,\n" <>
          "  \"constraints\" -> {}\n" <>
          "|>];\n" <>
          "zeroResult = CyclotomicCommonKernel[zeroCompiled];\n" <>
          "{\n" <>
          "  KeyTake[zeroResult, {\"iteration_nullities\", \"rank\", \"nullity\"}],\n" <>
          "  KeyTake[\n" <>
          "    RestoreCyclotomicExpression[zeroResult[\"basis_matrix\"]],\n" <>
          "    {\"rows\", \"columns\", \"entries\"}\n" <>
          "  ]\n" <>
          "}",
        cyclotomicDocumentationFixture["ZeroCommonSummary"]
      }
    },
    "Issues" -> {{
      "The function accepts only a complete compiled problem. It does not infer missing context or dimensions from an arbitrary Association:",
      "commonFailure = CyclotomicCommonKernel[<||>];\n" <>
        "<|\"FailureTag\" -> commonFailure[[1]]|>",
      <|"FailureTag" -> cyclotomicDocumentationFixture["MalformedFailure"]|>
    }},
    "SeeAlso" -> {
      {"CompileCyclotomicMatrices", "MagneticTB/ref/CompileCyclotomicMatrices"},
      {"CyclotomicNullSpace", "MagneticTB/ref/CyclotomicNullSpace"},
      {"CyclotomicCommonNullSpace", "MagneticTB/ref/CyclotomicCommonNullSpace"},
      {"RestoreCyclotomicExpression", "MagneticTB/ref/RestoreCyclotomicExpression"}
    },
    "Keywords" -> {"common kernel", "exact constraints", "iterative null space"}
  |>,
  <|
    "Name" -> "CyclotomicCommonNullSpace",
    "Usage" -> {{
      "CyclotomicCommonNullSpace[{matrix1,matrix2,...}]",
      "accepts ordinary exact Wolfram Language matrices and returns a row basis for their common null space."
    }},
    "Notes" -> {
      StringJoin[
        "This is the main entry for ordinary exact matrices. If rows is the returned list, Transpose[rows] ",
        "is the column-basis matrix and every input block Ai satisfies Ai . Transpose[rows] == 0 exactly."
      ],
      StringJoin[
        "Supported scalar inputs include exact integers and rationals, I, square roots of exact rationals, ",
        "roots of unity written with rational-angle Exp, Sin, or Cos, Conjugate, and exact sums, products, ",
        "and supported powers of these expressions. Machine numbers, symbolic parameters, and unsupported ",
        "algebraic numbers such as general cube roots return Failure."
      ],
      "All matrices must have the same column count. Empty constraints require CoordinateDimension because no matrix is available from which to infer that count.",
      StringJoin[
        "The function compiles one exact field and calls the production iterative common-kernel algorithm. ",
        "It does not use the built-in NullSpace or RowReduce and never changes to floating point, SVD, ",
        "stacked constraints, or another fallback."
      ]
    },
    "OptionTable" -> {
      {
        "\"CoordinateDimension\"", "Automatic",
        "Number of vector coordinates. Automatic infers the common column count from nonempty input; an explicit nonnegative integer is required when the constraint list is empty."
      },
      {
        "\"MaxCyclotomicDegree\"", "128",
        "Positive integer upper bound on the exact cyclotomic field degree. Exceeding the bound returns Failure instead of using approximate arithmetic."
      }
    },
    "Examples" -> {{
      "Find the common null space of two ordinary exact constraint blocks containing three independent square roots. The two displayed residual matrices are exactly zero:",
      "constraintBlocks = {\n" <>
        "  {\n" <>
        "    {1, 0, 0, -Sqrt[2], 0},\n" <>
        "    {0, 1, 0, -Sqrt[3], 0}\n" <>
        "  },\n" <>
        "  {{0, 0, 1, -Sqrt[5], 0}}\n" <>
        "};\n" <>
        "commonRows = CyclotomicCommonNullSpace[constraintBlocks];\n" <>
        "{\n" <>
        "  MatrixForm[commonRows],\n" <>
        "  MatrixForm /@ (# . Transpose[commonRows] & /@ constraintBlocks)\n" <>
        "}",
      {
        MatrixForm[cyclotomicDocumentationFixture["DirectRows"]],
        cyclotomicDocumentationFixture["DirectResiduals"]
      }
    }},
    "OptionExamples" -> {
      {
        "CoordinateDimension fixes the vector space when the constraint list is empty. With three coordinates, the complete row basis is the identity:",
        "CyclotomicCommonNullSpace[\n" <>
          "  {},\n" <>
          "  \"CoordinateDimension\" -> 3\n" <>
          "]",
        cyclotomicDocumentationFixture["DirectEmpty"]
      },
      {
        "MaxCyclotomicDegree bounds the exact field degree, not a numerical tolerance. The fourth-root field has degree two, so limit two succeeds and limit one fails:",
        "{\n" <>
          "  CyclotomicCommonNullSpace[\n" <>
          "    {{{I, 1}}},\n" <>
          "    \"MaxCyclotomicDegree\" -> 2\n" <>
          "  ],\n" <>
          "  CyclotomicCommonNullSpace[\n" <>
          "    {{{I, 1}}},\n" <>
          "    \"MaxCyclotomicDegree\" -> 1\n" <>
          "  ][[1]]\n" <>
          "}",
        {
          cyclotomicDocumentationFixture["DegreeSuccess"],
          cyclotomicDocumentationFixture["LowDegreeFailure"]
        }
      }
    },
    "Properties" -> {{
      "For a zero-dimensional vector space, empty constraints return the unique empty row basis:",
      "CyclotomicCommonNullSpace[\n" <>
        "  {},\n" <>
        "  \"CoordinateDimension\" -> 0\n" <>
        "]",
      cyclotomicDocumentationFixture["DirectZero"]
    }},
    "Issues" -> {
      {
        "Without a nonempty matrix or an explicit CoordinateDimension, the vector-space dimension is undefined:",
        "missingDimension = CyclotomicCommonNullSpace[{}];\n" <>
          "<|\"FailureTag\" -> missingDimension[[1]]|>",
        <|"FailureTag" -> cyclotomicDocumentationFixture["MissingDimensionFailure"]|>
      },
      {
        "All blocks must have the same number of columns:",
        "dimensionFailure = CyclotomicCommonNullSpace[{\n" <>
          "  {{1, 0}},\n" <>
          "  {{1, 0, 0}}\n" <>
          "}];\n" <>
          "<|\"FailureTag\" -> dimensionFailure[[1]]|>",
        <|"FailureTag" -> cyclotomicDocumentationFixture["DimensionFailure"]|>
      },
      {
        "Machine numbers and unsupported exact algebraic expressions fail for different reasons; neither is approximated or sent to a numerical fallback:",
        "{\n" <>
          "  CyclotomicCommonNullSpace[{{{0.5}}}][[1]],\n" <>
          "  CyclotomicCommonNullSpace[{{{2^(1/3)}}}][[1]]\n" <>
          "}",
        {
          cyclotomicDocumentationFixture["NonExactFailure"],
          cyclotomicDocumentationFixture["UnsupportedFailure"]
        }
      }
    },
    "SeeAlso" -> {
      {"CompileCyclotomicMatrices", "MagneticTB/ref/CompileCyclotomicMatrices"},
      {"CyclotomicCommonKernel", "MagneticTB/ref/CyclotomicCommonKernel"},
      {"CyclotomicNullSpace", "MagneticTB/ref/CyclotomicNullSpace"},
      {"RestoreCyclotomicExpression", "MagneticTB/ref/RestoreCyclotomicExpression"}
    },
    "Keywords" -> {"exact common null space", "ordinary matrices", "cyclotomic field"}
  |>,
  <|
    "Name" -> "RestoreCyclotomicExpression",
    "Usage" -> {
      {
        "RestoreCyclotomicExpression[object]",
        "turns an exact matrix, null-space result, or common-kernel result into ordinary exact Wolfram Language expressions."
      },
      {
        "RestoreCyclotomicExpression[context,element]",
        "restores one exact cyclotomic element using its matching field context."
      }
    },
    "Notes" -> {
      StringJoin[
        "The one-argument form accepts an exact_matrix, kernel_result, or common_kernel_result produced by ",
        "the public cyclotomic functions. A restored matrix records \"rows\", \"columns\", and row-major ",
        "flat \"entries\"."
      ],
      "For an RREF result, restore its \"reduced_matrix\" field. The complete rref_result Association is metadata and is not itself a restorable object.",
      "The two-argument form requires the exact context to which element belongs. A context mismatch or malformed object returns Failure; no numerical approximation is used."
    },
    "Examples" -> {{
      "Compile a matrix containing the fourth root of unity and restore it as the ordinary exact expression I:",
      cyclotomicStructuredCompileCode <>
        "\nrestoredMatrix = RestoreCyclotomicExpression[\n" <>
        "  compiledCyclotomic[\"constraints\"][[1]]\n" <>
        "];\n" <>
        "MatrixForm@Partition[\n" <>
        "  restoredMatrix[\"entries\"],\n" <>
        "  restoredMatrix[\"columns\"]\n" <>
        "]",
      cyclotomicDocumentationFixture["RestoredStructuredMatrix"]
    }},
    "Applications" -> {{
      "Restore a complete kernel result, then display its row-basis and column-basis matrices:",
      cyclotomicRationalCompileCode <>
        "\nkernelResult = CyclotomicNullSpace[exactRationalMatrix];\n" <>
        "restoredKernel = RestoreCyclotomicExpression[kernelResult];\n" <>
        "{\n" <>
        "  MatrixForm@Partition[\n" <>
        "    restoredKernel[\"nullspace_rows\", \"entries\"],\n" <>
        "    restoredKernel[\"nullspace_rows\", \"columns\"]\n" <>
        "  ],\n" <>
        "  MatrixForm@Partition[\n" <>
        "    restoredKernel[\"basis_matrix\", \"entries\"],\n" <>
        "    restoredKernel[\"basis_matrix\", \"columns\"]\n" <>
        "  ]\n" <>
        "}",
      cyclotomicDocumentationFixture["RestoredKernel"]
    }},
    "Properties" -> {{
      "The two-argument form restores one element from its matching context. The second entry of this matrix is exactly I:",
      cyclotomicStructuredCompileCode <>
        "\nRestoreCyclotomicExpression[\n" <>
        "  compiledCyclotomic[\"context\"],\n" <>
        "  compiledCyclotomic[\"constraints\"][[1, \"entries\", 2]]\n" <>
        "]",
      cyclotomicDocumentationFixture["RestoredIElement"]
    }},
    "Issues" -> {{
      "An arbitrary Association is not an exact cyclotomic object and fails explicitly:",
      "restoreFailure = RestoreCyclotomicExpression[<||>];\n" <>
        "<|\"FailureTag\" -> restoreFailure[[1]]|>",
      <|"FailureTag" -> cyclotomicDocumentationFixture["MalformedFailure"]|>
    }},
    "SeeAlso" -> {
      {"CompileCyclotomicMatrices", "MagneticTB/ref/CompileCyclotomicMatrices"},
      {"CyclotomicRREF", "MagneticTB/ref/CyclotomicRREF"},
      {"CyclotomicNullSpace", "MagneticTB/ref/CyclotomicNullSpace"},
      {"CyclotomicCommonKernel", "MagneticTB/ref/CyclotomicCommonKernel"}
    },
    "Keywords" -> {"restore exact expression", "cyclotomic element", "exact matrix"}
  |>,
  <|
    "Name" -> "unsymham",
    "Usage" -> {{
      "unsymham[n]",
      "returns the translation-periodic Hamiltonian of shell n before point-group, magnetic, or Hermitian constraints are imposed."
    }},
    "Notes" -> {
      "Every directed hopping matrix element is independent, although the bond geometry and Bloch phases still come from the lattice prepared by init.",
      "This is mainly useful for understanding how symmetry reduces the parameter space. Use symham to construct the physical symmetry-constrained model."
    },
    "Examples" -> {{
      "Initialize the complete cubic gray-group 221 s-orbital model, then display the independent complex amplitude carried by every directed nearest-neighbour bond:",
      simpleCubicUtilityInitCode <>
        "\nMatrixForm[unsymham[2]]",
      MatrixForm[advancedUtilityFixture["Unconstrained"]]
    }},
    "Properties" -> {{
      "For the same prepared shell, symham reduces these independent entries to the symmetry-allowed Hermitian subspace:",
      "{MatrixForm[unsymham[2]], MatrixForm[symham[2]]}",
      {
        MatrixForm[advancedUtilityFixture["Unconstrained"]],
        MatrixForm[advancedUtilityFixture["Constrained"]]
      }
    }},
    "Issues" -> {{
      "The shell must already have been prepared by init; the function never extends the bond search:",
      "unsymham[3]",
      $Failed,
      <|
        "ExpectedMessages" -> {"unsymham::shellrange"},
        "MessageText" ->
          "unsymham::shellrange: Bond shell 3 was not prepared by init; only 2 shells are available. Rerun init with InitialBondShells -> 3."
      |>
    }},
    "SeeAlso" -> {{"symham", "MagneticTB/ref/symham"}, {"init", "MagneticTB/ref/init"}},
    "Expected" -> {"msgop[gray[221]]", "complete cubic gray-group 221"},
    "Keywords" -> {"unconstrained Hamiltonian", "directed hopping", "diagnostic"}
  |>,
  <|
    "Name" -> "CompileMagneticTBInput",
    "Usage" -> {{
      "CompileMagneticTBInput[input]",
      "checks and prepares program-generated MagneticTB input without making it the current model."
    }},
    "Notes" -> {
      StringJoin[
        "Use this function when another program has already assembled the lattice, Wyckoff positions, ",
        "symmetry operations, and local orbitals of a model. For an ordinary notebook calculation, init or ",
        "initfromrep is simpler."
      ],
      StringJoin[
        "The string keys of input describe the same physical quantities used by init. The result lets the ",
        "calling program inspect the generated sites, representation matrices, and bond geometry. It does ",
        "not calculate any hopping parameter or replace the model currently used by symham."
      ]
    },
    "Examples" -> {{
      "Here is a complete simple-cubic model in this input form. The output shows the lattice, generated sites, exact representation matrices, and first bond class:",
      compilerInputCode <>
        "\ncompiledModel = CompileMagneticTBInput[compilerInput];\n" <>
        "Grid[List @@@ Normal@<|\n" <>
        "  \"Lattice\" -> compiledModel[\"ModelSpecification\", \"Lattice\"],\n" <>
        "  \"SiteOrbits\" -> compiledModel[\"ModelSpecification\", \"SiteOrbits\"],\n" <>
        "  \"RepresentationMatrices\" -> compiledModel[\"FullRepresentation\", \"RepresentationMatrices\"],\n" <>
        "  \"FirstBondShell\" -> compiledModel[\"ModelSpecification\", \"BondClasses\"][[1]]\n" <>
        "|>, Frame -> All, Alignment -> Left]",
      advancedUtilityFixture["CompilerSummary"]
    }},
    "Applications" -> {{
      StringJoin[
        "When another program writes the model input, it can first inspect the exact symmetry matrices ",
        "returned here. The same physical input can then be passed to init or initfromrep for an interactive ",
        "calculation:"
      ],
      compilerInputCode <>
        "\nCompileMagneticTBInput[compilerInput][\"FullRepresentation\", \"RepresentationMatrices\"]",
      {{{1}}}
    }},
    "Issues" -> {{
      "All physical inputs must be present. For example, a lattice by itself does not define a model:",
      "CompileMagneticTBInput[<|\"Lattice\" -> IdentityMatrix[3]|>]",
      $Failed,
      <|
        "ExpectedMessages" -> {"CompileMagneticTBInput::input"},
        "MessageText" ->
          StringJoin[
            "CompileMagneticTBInput::input: The init options are incomplete or mutually incompatible. ",
            "symminformation must be either traditional {label,R,t,F|T} records or a discrete spin-space-group ",
            "label-><|space->{R,t},spin->{S,0|1}|> Association."
          ]
      |>
    }},
    "SeeAlso" -> {{"init", "MagneticTB/ref/init"}, {"initfromrep", "MagneticTB/ref/initfromrep"}},
    "Keywords" -> {"programmatic model input", "model preparation", "representation matrix"}
  |>,
  <|
    "Name" -> "pointMatrix",
    "Usage" -> {{
      "pointMatrix[operations, basis, lattice]",
      "gives the exact matrix of each point operation on a chosen local orbital or spinor basis."
    }},
    "Notes" -> {
      StringJoin[
        "Each operation is written as {label,R,t,flag}. Only its point action enters the local orbital ",
        "matrix; the translation is retained so the same operation data can be used by init."
      ],
      StringJoin[
        "basis may contain scalar orbitals such as px and py or two-component spinors. The matrices are ",
        "returned in the same operation order and act on the basis in the order supplied."
      ],
      "For an ordered basis {f1,...,fn}, column j contains the coefficients of the transformed state g fj in that same basis: g fj = Sum_i fi D_ij(g)."
    },
    "Examples" -> {{
      "The full unitary tetragonal C4 point group from gray group 75 acts on the ordered {px,py} basis as follows:",
      StringRiffle[{
        "pointOperations = Select[",
        "  msgop[gray[75]],",
        "  Last[#] === \"F\" &",
        "];",
        "MatrixForm /@ pointMatrix[",
        "  pointOperations,",
        "  {\"px\", \"py\"},",
        "  IdentityMatrix[3]",
        "]"
      }, "\n"],
      MatrixForm /@ advancedUtilityFixture["PointMatrices"]
    }},
    "Issues" -> {{
      "A singular lattice cannot define the Cartesian point action:",
      "pointMatrix[{{\"E\", IdentityMatrix[3], {0, 0, 0}, \"F\"}}, {\"s\"}, ConstantArray[0, {3, 3}]]",
      $Failed,
      <|
        "ExpectedMessages" -> {"pointMatrix::input"},
        "MessageText" ->
          "pointMatrix::input: pointMatrix expects a nonempty symmetry-operation list, a nonempty scalar or two-component spinor basis, and a nonsingular 3 by 3 evaluated lattice."
      |>
    }},
    "Expected" -> {"full unitary tetragonal C4 point group"},
    "SeeAlso" -> {{"showSymmetryRepresentations", "MagneticTB/ref/showSymmetryRepresentations"}, {"init", "MagneticTB/ref/init"}},
    "Keywords" -> {"point representation", "orbital basis", "unitary matrix"}
  |>,
  <|
    "Name" -> "GenerateGroup",
    "Usage" -> {{
      "GenerateGroup[generators, identity, multiply]",
      "generates a complete finite group from its generators and an explicitly supplied multiplication law."
    }},
    "Notes" -> {
      StringJoin[
        "For example, starting from one generator of C4, repeated multiplication produces its four elements. ",
        "The returned list always begins with the identity; the trivial group is therefore {E}, not an empty ",
        "list."
      ],
      "The caller supplies both multiplication and equality because magnetic-space-group and spin-space-group elements need not use the same concrete product."
    },
    "OptionTable" -> {{
      "SameTest", "SameQ",
      "Binary predicate used to decide whether two generated concrete elements are equal. It must return True or False for every pair."
    }},
    "Examples" -> {{
      "Generate the additive cyclic group C4 from the element 1:",
      "GenerateGroup[{1}, 0, Mod[#1 + #2, 4] &]",
      advancedUtilityFixture["C4Group"]
    }},
    "OptionExamples" -> {{
      "SameTest can provide a domain-specific exact equality predicate; here it reproduces the same C4 ordering:",
      "GenerateGroup[{1}, 0, Mod[#1 + #2, 4] &, SameTest -> (#1 === #2 &)]",
      advancedUtilityFixture["C4Group"]
    }},
    "SeeAlso" -> {{"getGenerator", "MagneticTB/ref/getGenerator"}, {"init", "MagneticTB/ref/init"}},
    "Keywords" -> {"finite group", "generators", "closure"}
  |>,
  <|
    "Name" -> "getGenerator",
    "Usage" -> {{
      "getGenerator[group, identity, multiply]",
      "finds a small generating subset of a complete finite group."
    }},
    "Notes" -> {
      "When the complete finite group is already known, getGenerator looks for a smaller set whose products recover every group element.",
      StringJoin[
        "Supply the complete group, its identity, and its multiplication law. The function does not infer a ",
        "product from the appearance of the elements. The name is retained so that existing notebooks ",
        "continue to run."
      ]
    },
    "OptionTable" -> {{
      "SameTest", "SameQ",
      "Binary predicate used for exact equality of concrete group elements."
    }},
    "Examples" -> {{
      "Recover one generator of the complete additive C4 group:",
      "getGenerator[Range[0, 3], 0, Mod[#1 + #2, 4] &]",
      advancedUtilityFixture["C4Generators"]
    }},
    "OptionExamples" -> {{
      "An explicit exact SameTest leaves the mathematical result unchanged:",
      "getGenerator[Range[0, 3], 0, Mod[#1 + #2, 4] &, SameTest -> (#1 === #2 &)]",
      advancedUtilityFixture["C4Generators"]
    }},
    "SeeAlso" -> {{"GenerateGroup", "MagneticTB/ref/GenerateGroup"}, {"init", "MagneticTB/ref/init"}},
    "Keywords" -> {"finite group", "generator search", "deprecated compatibility"}
  |>,
  <|
    "Name" -> "symmetrizationHRInit",
    "Usage" -> {{
      "symmetrizationHRInit[opts]",
      "prepares the orbital centers and symmetry matrices needed to symmetrize a Wannier90 HR model."
    }},
    "Notes" -> {
      StringJoin[
        "To symmetrize a Wannier90 HR model, one needs the center of every Wannier orbital and the matrix of ",
        "every symmetry operation in the same orbital order. This function constructs those data from the ",
        "lattice, Wyckoff positions, symmetry operations, and local orbitals."
      ],
      "The orbital order must agree with the external electronic-structure calculation. At present the function supports the documented VASP ordering.",
      "The result contains wcc for the fractional Wannier centers, DR for the full symmetry matrices, and symmetry for the corresponding ordered operations."
    },
    "OptionTable" -> {
      {
        "\"Lattice\"",
        "{{a,0,0},{-a/2,Sqrt[3] a/2,0},{0,0,c}}",
        "Symbolic direct-lattice vectors stored by rows."
      },
      {
        "\"LattPar\"",
        "{a->1,c->3}",
        "Numerical replacement rules for the symbolic lattice."
      },
      {
        "\"WyckoffPosition\"",
        "{{{2/3,1/3,0},{0,0,1/2}}}",
        "Wyckoff seeds and magnetic moments in the same nesting used by init."
      },
      {
        "\"Symminformation\"",
        "{{\"1\",IdentityMatrix[3],{0,0,0},\"F\"}}",
        "Complete ordered symmetry records."
      },
      {
        "\"BasisFunctions\"",
        "{\"s\"}",
        "Local basis used to build Wannier centers and point matrices."
      },
      {"\"Software\"", "\"VASP\"", "Orbital ordering convention; VASP is the only supported value."}
    },
    "Examples" -> {
      {
        "Prepare a full tetragonal C4 model with the complete gray-group 75 operation list and a px/py basis, then show a compact summary instead of a large raw Association:",
        symmetrizationOptionCode <> "\n" <> StringRiffle[{
          "Grid[{",
          "  {\"Wannier centers\", c4SymmetrizationData[\"wcc\"]},",
          "  {\"Operation count\", Length[c4SymmetrizationData[\"symmetry\"]]},",
          "  {\"Matrix dimensions\", Dimensions /@ c4SymmetrizationData[\"DR\"]}",
          "}, Frame -> All]"
        }, "\n"],
        documentationKeyValueGrid@<|
          "Wannier centers" -> advancedUtilityFixture["SymmetrizationData"]["wcc"],
          "Operation count" -> Length[advancedUtilityFixture["SymmetrizationData"]["symmetry"]],
          "Matrix dimensions" -> (Dimensions /@ advancedUtilityFixture["SymmetrizationData"]["DR"])
        |>
      }
    },
    "OptionExamples" -> {
      {
        "\"Lattice\" -> DiagonalMatrix[{2 a,2 a,3 c}] changes the tetragonal direct-lattice rows while preserving C4 symmetry. All eight operation matrices remain 2 by 2:",
        StringRiffle[{
          "Dimensions /@ symmetrizationHRInit[",
          "  \"Lattice\" -> DiagonalMatrix[{2 a, 2 a, 3 c}],",
          "  \"LattPar\" -> c4SymmetrizationParameters,",
          "  \"WyckoffPosition\" -> c4SymmetrizationSites,",
          "  \"Symminformation\" -> c4SymmetrizationOperations,",
          "  \"BasisFunctions\" -> c4SymmetrizationBasis,",
          "  \"Software\" -> \"VASP\"",
          "][\"DR\"]"
        }, "\n"],
        Dimensions /@ advancedUtilityFixture["SymmetrizationData"]["DR"]
      },
      {
        "\"LattPar\" -> {a->2,c->3} evaluates every symbolic tetragonal lattice entry before the symmetry matrices are constructed:",
        StringRiffle[{
          "MatrixForm@First@symmetrizationHRInit[",
          "  \"Lattice\" -> c4SymmetrizationLattice,",
          "  \"LattPar\" -> {a -> 2, c -> 3},",
          "  \"WyckoffPosition\" -> c4SymmetrizationSites,",
          "  \"Symminformation\" -> c4SymmetrizationOperations,",
          "  \"BasisFunctions\" -> c4SymmetrizationBasis,",
          "  \"Software\" -> \"VASP\"",
          "][\"DR\"]"
        }, "\n"],
        MatrixForm@First@advancedUtilityFixture["SymmetrizationData"]["DR"]
      },
      {
        "\"WyckoffPosition\" -> {{{1/2,1/2,0},{0,0,0}}} places both C4-related local orbitals on the invariant fractional center {1/2,1/2,0}:",
        StringRiffle[{
          "symmetrizationHRInit[",
          "  \"Lattice\" -> c4SymmetrizationLattice,",
          "  \"LattPar\" -> c4SymmetrizationParameters,",
          "  \"WyckoffPosition\" -> {{{1/2, 1/2, 0}, {0, 0, 0}}},",
          "  \"Symminformation\" -> c4SymmetrizationOperations,",
          "  \"BasisFunctions\" -> c4SymmetrizationBasis,",
          "  \"Software\" -> \"VASP\"",
          "][\"wcc\"]"
        }, "\n"],
        advancedUtilityFixture["ShiftedSymmetrizationData"]["wcc"]
      },
      {
        "\"Symminformation\" may select the complete unitary C4 subgroup from gray group 75. The returned records preserve all four ordered operations:",
        StringRiffle[{
          "symmetrizationHRInit[",
          "  \"Lattice\" -> c4SymmetrizationLattice,",
          "  \"LattPar\" -> c4SymmetrizationParameters,",
          "  \"WyckoffPosition\" -> c4SymmetrizationSites,",
          "  \"Symminformation\" -> Select[c4SymmetrizationOperations, Last[#] == \"F\" &],",
          "  \"BasisFunctions\" -> c4SymmetrizationBasis,",
          "  \"Software\" -> \"VASP\"",
          "][\"symmetry\"][[All, {1, 4}]]"
        }, "\n"],
        advancedUtilityFixture["UnitarySymmetrizationData"]["symmetry"][[All, {1, 4}]]
      },
      {
        "\"BasisFunctions\" -> {\"s\",\"pz\"} changes the VASP-ordered local basis while retaining the same complete C4 symmetry list:",
        StringRiffle[{
          "MatrixForm@symmetrizationHRInit[",
          "  \"Lattice\" -> c4SymmetrizationLattice,",
          "  \"LattPar\" -> c4SymmetrizationParameters,",
          "  \"WyckoffPosition\" -> c4SymmetrizationSites,",
          "  \"Symminformation\" -> c4SymmetrizationOperations,",
          "  \"BasisFunctions\" -> {\"s\", \"pz\"},",
          "  \"Software\" -> \"VASP\"",
          "][\"DR\"][[4]]"
        }, "\n"],
        MatrixForm[advancedUtilityFixture["ScalarSymmetrizationData"]["DR"][[4]]]
      },
      {
        "\"Software\" -> \"VASP\" selects the currently supported external orbital-order convention. The result has exactly the documented three keys:",
        StringRiffle[{
          "Keys@symmetrizationHRInit[",
          "  \"Lattice\" -> c4SymmetrizationLattice,",
          "  \"LattPar\" -> c4SymmetrizationParameters,",
          "  \"WyckoffPosition\" -> c4SymmetrizationSites,",
          "  \"Symminformation\" -> c4SymmetrizationOperations,",
          "  \"BasisFunctions\" -> c4SymmetrizationBasis,",
          "  \"Software\" -> \"VASP\"",
          "]"
        }, "\n"],
        Keys[advancedUtilityFixture["SymmetrizationData"]]
      }
    },
    "Issues" -> {{
      "Only the documented VASP orbital ordering is currently available:",
      "symmetrizationHRInit[\"Software\" -> \"Other\"]",
      $Failed,
      <|
        "ExpectedMessages" -> {"symmetrizationHRInit::software"},
        "MessageText" ->
          "symmetrizationHRInit::software: Only the documented software convention VASP is currently supported; received Other."
      |>
    }},
    "SeeAlso" -> {{"hop", "MagneticTB/ref/hop"}, {"readHR", "MagneticTB/ref/readHR"}},
    "Expected" -> {"msgop[gray[75]]", "complete gray-group 75"},
    "Keywords" -> {"Wannier90", "HR symmetrization", "Wannier centers"}
  |>,
  <|
    "Name" -> "gray",
    "Usage" -> {{
      "gray[n]",
      "returns the database identifier of the gray magnetic space group associated with ordinary space-group number n."
    }},
    "Notes" -> {
      "gray is an Association whose integer keys run from 1 through 230.",
      "Pass the returned integer to msgop or showMSGWyckoff. It is a database identifier, not an operation list."
    },
    "Examples" -> {
      {
        "Look up the gray magnetic-space-group identifier associated with ordinary space group 191:",
        "gray[191]",
        1464
      },
      {
        StringJoin[
          "Use the identifier returned by gray[191] directly with msgop. The selected records ",
          "show two unitary operations and pure time reversal:"
        ],
        "msgop[gray[191]][[{1, 2, 25}]]",
        grapheneFixture["SelectedOperations"]
      }
    },
    "Issues" -> {{
      "A number outside 1 through 230 is not a key of the association:",
      "gray[0]",
      Missing["KeyAbsent", 0]
    }},
    "SeeAlso" -> {
      {"msgop", "MagneticTB/ref/msgop"},
      {"showMSGWyckoff", "MagneticTB/ref/showMSGWyckoff"},
      {"bnsdict", "MagneticTB/ref/bnsdict"},
      {"typeI", "MagneticTB/ref/typeI"},
      {"typeIII", "MagneticTB/ref/typeIII"},
      {"typeIV", "MagneticTB/ref/typeIV"}
    },
    "Keywords" -> {"gray magnetic space group", "database identifier", "BNS"}
  |>,
  <|
    "Name" -> "bnsdict",
    "Usage" -> {{
      "bnsdict[{n,m}]",
      "returns the bundled database identifier of the magnetic space group with BNS number n.m."
    }},
    "Notes" -> {
      "bnsdict is an Association containing all 1651 BNS-number pairs in the bundled magnetic-space-group database.",
      "Pass the returned integer to msgop or showMSGWyckoff. It is a database identifier, not an operation list.",
      StringJoin[
        "The values of typeI, gray, typeIII, and typeIV partition the same 1651 identifiers. ",
        "Use bnsdict when the complete BNS pair is already known."
      ]
    },
    "Examples" -> {
      {
        "Look up the database identifier of the nontrivial type-III magnetic space group with BNS number 164.87:",
        "bnsdict[{164, 87}]",
        1317
      },
      {
        "The BNS pair 191.234 and gray[191] select the same gray magnetic space group:",
        "{bnsdict[{191, 234}], gray[191]}",
        {1464, 1464}
      },
      {
        StringJoin[
          "Use the identifier returned by bnsdict[{164, 87}] directly with msgop. The selected ",
          "records show unitary and antiunitary operations of this type-III group:"
        ],
        "msgop[bnsdict[{164, 87}]][[{1, 2, 7}]]",
        {
          {
            "1",
            IdentityMatrix[3],
            {0, 0, 0},
            "F"
          },
          {
            "3z",
            {{0, -1, 0}, {1, -1, 0}, {0, 0, 1}},
            {0, 0, 0},
            "F"
          },
          {
            "2x",
            {{1, -1, 0}, {0, -1, 0}, {0, 0, -1}},
            {0, 0, 0},
            "T"
          }
        }
      }
    },
    "Issues" -> {{
      "A BNS pair absent from the bundled association returns Missing:",
      "bnsdict[{75, 999}]",
      Missing["KeyAbsent", {75, 999}]
    }},
    "SeeAlso" -> {
      {"msgop", "MagneticTB/ref/msgop"},
      {"showMSGWyckoff", "MagneticTB/ref/showMSGWyckoff"},
      {"gray", "MagneticTB/ref/gray"},
      {"typeI", "MagneticTB/ref/typeI"},
      {"typeIII", "MagneticTB/ref/typeIII"},
      {"typeIV", "MagneticTB/ref/typeIV"}
    },
    "Keywords" -> {
      "BNS number", "magnetic space group", "database identifier"
    }
  |>,
  <|
    "Name" -> "typeI",
    "Usage" -> {{
      "typeI[{n,m}]",
      "returns the database identifier of the type-I magnetic space group with BNS number n.m."
    }},
    "Notes" -> {
      "typeI is an Association with 230 BNS-number pairs as keys.",
      "Pass the returned integer to msgop or showMSGWyckoff. The key must be the pair {n,m}, not n alone."
    },
    "Examples" -> {{
      "Look up the type-I magnetic-space-group identifier with BNS number 191.233:",
      "typeI[{191, 233}]",
      1463
    }},
    "Issues" -> {{
      "A single integer is not a typeI key; use the complete BNS pair:",
      "typeI[191]",
      Missing["KeyAbsent", 191]
    }},
    "SeeAlso" -> {
      {"msgop", "MagneticTB/ref/msgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"bnsdict", "MagneticTB/ref/bnsdict"},
      {"typeIII", "MagneticTB/ref/typeIII"},
      {"typeIV", "MagneticTB/ref/typeIV"}
    },
    "Keywords" -> {"type-I magnetic space group", "database identifier", "BNS"}
  |>,
  <|
    "Name" -> "typeIII",
    "Usage" -> {{
      "typeIII[{n,m}]",
      "returns the database identifier of the type-III magnetic space group with BNS number n.m."
    }},
    "Notes" -> {
      "typeIII is an Association with 674 BNS-number pairs as keys.",
      "Pass the returned integer to msgop or showMSGWyckoff. The key must be the pair {n,m}, not n alone."
    },
    "Examples" -> {{
      "Look up the type-III magnetic-space-group identifier with BNS number 164.87:",
      "typeIII[{164, 87}]",
      1317
    }},
    "Issues" -> {{
      "A BNS pair absent from the bundled association returns Missing:",
      "typeIII[{164, 999}]",
      Missing["KeyAbsent", {164, 999}]
    }},
    "SeeAlso" -> {
      {"msgop", "MagneticTB/ref/msgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"bnsdict", "MagneticTB/ref/bnsdict"},
      {"typeI", "MagneticTB/ref/typeI"},
      {"typeIV", "MagneticTB/ref/typeIV"}
    },
    "Keywords" -> {"type-III magnetic space group", "database identifier", "BNS"}
  |>,
  <|
    "Name" -> "typeIV",
    "Usage" -> {{
      "typeIV[{n,m}]",
      "returns the database identifier of the type-IV magnetic space group with BNS number n.m."
    }},
    "Notes" -> {
      "typeIV is an Association with 517 BNS-number pairs as keys.",
      "Pass the returned integer to msgop or showMSGWyckoff. The key must be the pair {n,m}, not n alone."
    },
    "Examples" -> {{
      "Look up the type-IV magnetic-space-group identifier with BNS number 75.4:",
      "typeIV[{75, 4}]",
      664
    }},
    "Issues" -> {{
      "A BNS pair absent from the bundled association returns Missing:",
      "typeIV[{75, 999}]",
      Missing["KeyAbsent", {75, 999}]
    }},
    "SeeAlso" -> {
      {"msgop", "MagneticTB/ref/msgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"bnsdict", "MagneticTB/ref/bnsdict"},
      {"typeI", "MagneticTB/ref/typeI"},
      {"typeIII", "MagneticTB/ref/typeIII"}
    },
    "Keywords" -> {"type-IV magnetic space group", "database identifier", "BNS"}
  |>,
  <|
    "Name" -> "graylayer",
    "Usage" -> {{
      "graylayer[n]",
      "returns the database identifier of the gray magnetic layer group associated with layer-group number n."
    }},
    "Notes" -> {
      "graylayer is an Association whose integer keys run from 1 through 80.",
      "Pass the returned three-integer identifier directly to mlgop."
    },
    "Examples" -> {{
      "Look up the gray magnetic-layer-group identifier associated with layer group 25:",
      "graylayer[25]",
      {25, 2, 144}
    }},
    "Issues" -> {{
      "A number outside 1 through 80 is not a key of the association:",
      "graylayer[81]",
      Missing["KeyAbsent", 81]
    }},
    "SeeAlso" -> {
      {"mlgop", "MagneticTB/ref/mlgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"grayrod", "MagneticTB/ref/grayrod"}
    },
    "Keywords" -> {"gray magnetic layer group", "database identifier"}
  |>,
  <|
    "Name" -> "grayrod",
    "Usage" -> {{
      "grayrod[n]",
      "returns the database identifier of the gray magnetic rod group associated with rod-group number n."
    }},
    "Notes" -> {
      "grayrod is an Association whose integer keys run from 1 through 75.",
      "Pass the returned three-integer identifier directly to mrgop."
    },
    "Examples" -> {{
      "Look up the gray magnetic-rod-group identifier associated with rod group 25:",
      "grayrod[25]",
      {25, 2, 129}
    }},
    "Issues" -> {{
      "A number outside 1 through 75 is not a key of the association:",
      "grayrod[76]",
      Missing["KeyAbsent", 76]
    }},
    "SeeAlso" -> {
      {"mrgop", "MagneticTB/ref/mrgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"graylayer", "MagneticTB/ref/graylayer"}
    },
    "Keywords" -> {"gray magnetic rod group", "database identifier"}
  |>,
  <|
    "Name" -> "msgop",
    "Usage" -> {
      {
        "msgop[group]",
        "reads the symmetry operations of a magnetic space group from the bundled database."
      }
    },
    "Notes" -> {
      "Use gray[n] for a gray group. Type-I, type-III, and type-IV groups use complete BNS pairs: typeI[{n,m}], typeIII[{n,m}], and typeIV[{n,m}].",
      StringJoin[
        "Each returned item is {label,R,t,flag}: R and t are the point and translation parts in fractional ",
        "coordinates, and flag states whether the operation is unitary or antiunitary."
      ],
      "The list may be passed directly to symminformation in init. Its order is also the order used when displaying the corresponding representation matrices."
    },
    "Examples" -> {
      {
        "Read gray magnetic space group 191 and inspect representative unitary and antiunitary records. Each record is {label, point matrix, translation, antiunitary flag}:",
        StringRiffle[{
          "operations191 = msgop[gray[191]];",
          "operations191[[{1, 2, 25}]]"
        }, "\n"],
        grapheneFixture["SelectedOperations"]
      }
    },
    "SeeAlso" -> {
      {"gray", "MagneticTB/ref/gray"},
      {"bnsdict", "MagneticTB/ref/bnsdict"},
      {"typeI", "MagneticTB/ref/typeI"},
      {"typeIII", "MagneticTB/ref/typeIII"},
      {"typeIV", "MagneticTB/ref/typeIV"},
      {"init", "MagneticTB/ref/init"}
    },
    "Keywords" -> {"magnetic space group", "symmetry operations", "BNS", "gray group"}
  |>,
  <|
    "Name" -> "init",
    "Usage" -> {
      {
        "init[opts]",
        "sets up a tight-binding model from its lattice, Wyckoff positions, symmetry, and local orbitals."
      }
    },
    "Notes" -> {
      StringJoin[
        "The five main inputs have a direct physical meaning: lattice and lattpar define the Bravais lattice,",
        " wyckoffposition specifies the inequivalent sites and their moments, symminformation gives the ",
        "symmetry operations, and basisFunctions specifies the local orbital space on each Wyckoff orbit."
      ],
      StringJoin[
        "After init finishes, first inspect the generated atoms, the orbital order, and the neighbour shells.",
        " The hopping parameters are not needed at this stage; symham determines them later, one shell at a ",
        "time."
      ],
      StringJoin[
        "By default init finds the first 10 complete bond shells. If a model needs a farther shell, choose a ",
        "larger InitialBondShells value and run init again; symham itself does not search for new bonds."
      ],
      StringJoin[
        "RepresentationMode selects how the full orbital representation is built. Use \"DirectProduct\" when ",
        "the action on the complete local basis is known, and \"Induced\" when the input starts from a ",
        "representation of a site-symmetry group."
      ],
      StringJoin[
        "Every call to init starts a new model. If the input is inconsistent, the calculation stops and asks ",
        "the user to correct it; MagneticTB does not silently continue with the preceding model."
      ]
    },
    "OptionTable" -> {
      {
        "lattice",
        "{{a,0,0},{-a/2,Sqrt[3] a/2,0},{0,0,c}}",
        "Primitive direct-lattice vectors stored by rows. After applying lattpar, the result must be a real, nonsingular numerical 3 by 3 matrix."
      },
      {
        "lattpar",
        "{a->1,c->3}",
        "Replacement rules that assign numerical values to every symbolic lattice parameter needed by lattice."
      },
      {
        "wyckoffposition",
        "{{{2/3,1/3,0},{0,0,1/2}}}",
        "One {fractional seed position, magnetic moment} pair for each Wyckoff orbit. The outer-list order also fixes the order of basisFunctions."
      },
      {
        "symminformation",
        "{{\"1\",IdentityMatrix[3],{0,0,0},\"F\"}}",
        "Either ordered magnetic-group records {label,R,t,\"F\"|\"T\"} or a discrete spin-space-group Association. The antiunitary flag belongs to each operation record."
      },
      {
        "basisFunctions",
        "{\"s\"}",
        "A nonempty local scalar-orbital or two-component spinor basis for every Wyckoff seed. Named orbitals and exact symbolic functions are accepted."
      },
      {
        "InitialBondShells",
        "10",
        "Number of complete bond shells found during initialization. Shell 1 contains the zero-distance onsite terms; shell 2 is the first nonzero-distance neighbour shell."
      },
      {
        "GenerateSymmetryGroup",
        "False",
        "False treats the supplied ordered records as the complete finite group. True closes supplied generators using the operation-specific multiplication law."
      },
      {
        "RepresentationMode",
        "\"DirectProduct\"",
        StringJoin[
          "\"DirectProduct\" uses the local action on the full basis together with the permutation of ",
          "equivalent sites. \"Induced\" derives the reference-site representation from basisFunctions and ",
          "transports it to the other sites. Use initfromrep when exact site-symmetry matrices are already known."
        ]
      }
    },
    "Examples" -> {
      {
        "Initialize graphene with the complete ordered list returned by msgop:",
        simpleInitCode
      },
      {
        "Inspect the Hamiltonian ordering. The two rows are the symmetry-related carbon sites generated from the single Wyckoff seed:",
        compactOrbitalTableCode,
        compactOrbitalTable[grapheneFixture["OrbitalTable"]]
      }
    },
    "OptionExamples" -> {
      {
        StringJoin[
          "lattice -> {{a,0,0},{-a/2,Sqrt[3]a/2,0},{0,0,c}} gives the three primitive ",
          "vectors of graphene by rows. The other four model inputs are supplied in the same call. ",
          "The second bond shell then displays the three nearest-neighbour carbon bonds fixed by ",
          "this lattice and Wyckoff position:"
        ],
        simpleInitCode <> "\nshowbonds[2]",
        grapheneFixture["BondGrid"]
      },
      {
        StringJoin[
          "lattpar -> {a->1,c->3} assigns numerical lengths to every symbol in lattice before ",
          "bonds are found. In these units showbonds[2] reports the graphene nearest-neighbour ",
          "distance and the three translated target positions:"
        ],
        simpleInitCode <> "\nshowbonds[2]",
        grapheneFixture["BondGrid"]
      },
      {
        StringJoin[
          "wyckoffposition -> {{{1/3,2/3,0},{0,0,0}}} supplies one nonmagnetic carbon seed in ",
          "fractional coordinates. The gray group generates the second carbon site; orbitalTable[] ",
          "shows both sites in the Hamiltonian order:"
        ],
        simpleInitCode <> "\norbitalTable[]",
        compactOrbitalTable[grapheneFixture["OrbitalTable"]]
      },
      {
        StringJoin[
          "symminformation -> msgop[gray[191]] supplies the complete ordered magnetic-space-group ",
          "records used for graphene. showSymmetryRepresentations[] displays the actual orbital ",
          "matrices of the selected generators, including their unitary or antiunitary type:"
        ],
        simpleInitCode <> "\nshowSymmetryRepresentations[]",
        grapheneFixture["SymmetryRepresentations"]
      },
      {
        StringJoin[
          "basisFunctions -> {{\"pz\"}} assigns one local pz orbital to the inequivalent carbon ",
          "seed. Symmetry transports it to the second site, so orbitalTable[] contains two pz basis ",
          "states and fixes their Hamiltonian row order:"
        ],
        simpleInitCode <> "\norbitalTable[]",
        compactOrbitalTable[grapheneFixture["OrbitalTable"]]
      },
      {
        StringJoin[
          "InitialBondShells -> 3 prepares onsite, nearest-neighbour, and next-nearest-neighbour ",
          "shells during the complete graphene initialization. showbonds[3] displays the six ",
          "next-nearest neighbours of each carbon site:"
        ],
        grapheneBondShellInitCode <> "\nshowbonds[3]",
        grapheneFixture["NextBondGrid"]
      },
      {
        StringJoin[
          "GenerateSymmetryGroup -> True closes the supplied C4z generator to the complete fourfold ",
          "group before constructing a px/py square-lattice model. The displayed matrices show the ",
          "four generated operations acting on the two local orbitals:"
        ],
        generatedC4InitCode <> "\nshowSymmetryRepresentations[All]",
        initOptionFixture["GeneratedC4Representations"]
      },
      {
        StringJoin[
          "In cubic symmetry, RepresentationMode -> \"DirectProduct\" requires the complete local ",
          "px, py, pz basis because rotations mix all three orbitals. The orbital table therefore ",
          "contains all three orbitals on every symmetry-equivalent site:"
        ],
        cubicPOptionInitCode["DirectProduct"] <>
          "\norbitalTable[]",
        validationModelFixture["CubicPDirect", "OrbitalTable"]
      },
      {
        StringJoin[
          "For the same cubic Wyckoff orbit, RepresentationMode -> \"Induced\" needs only the ",
          "reference px orbital at {1/2, 0, 0}. Cubic symmetry transports it to the corresponding ",
          "py and pz states at the other equivalent sites, as the orbital table shows:"
        ],
        cubicPOptionInitCode["Induced"] <>
          "\norbitalTable[]",
        validationModelFixture["CubicPInduced", "OrbitalTable"]
      },
      {
        StringJoin[
          "A spinful induced example uses the complete magnetic group P6'/m'mm' (BNS 191.239). ",
          "One exact in-plane spinor on the reference kagome site generates the other two local ",
          "states.  The first three neighbour-hopping layers reproduce the three-band quadratic-",
          "nodal-point model of Feng et al., Phys. Rev. B 109, 045423 (2024).  In MagneticTB's ",
          "basis and directed-bond convention, the paper parameters map to the three rules shown ",
          "in bandplot:"
        ],
        kagomeQNPInducedInitCode,
        kagomeQNPDocumentationFixture["BandPlot"]
      }
    },
    "SeeAlso" -> {
      {"initfromrep", "MagneticTB/ref/initfromrep"},
      {"symham", "MagneticTB/ref/symham"},
      {"CurrentModelSession", "MagneticTB/ref/CurrentModelSession"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"},
      {"showMSGWyckoff", "MagneticTB/ref/showMSGWyckoff"},
      {"msgop", "MagneticTB/ref/msgop"},
      {"gray", "MagneticTB/ref/gray"},
      {"bnsdict", "MagneticTB/ref/bnsdict"}
    },
    "Keywords" -> {"initialize model", "Wyckoff position", "tight binding", "bond shells"}
  |>,
  <|
    "Name" -> "initfromrep",
    "Usage" -> {
      {
        "initfromrep[opts]",
        "sets up a tight-binding model when the exact symmetry matrices are supplied directly instead of being derived from basis functions."
      }
    },
    "Notes" -> {
      StringJoin[
        "Use this interface when another calculation has already produced the symmetry matrices. It is ",
        "deliberately separate from init: basisFunctions, SiteLocalData, and GenerateSymmetryGroup are not ",
        "accepted."
      ],
      StringJoin[
        "symminformation must contain the complete finite group in a fixed order, and repinformation must ",
        "contain one exact unitary matrix for every operation in exactly that order."
      ],
      StringJoin[
        "For the supported continuous case, one pure-internal C_infinity generator may be supplied together ",
        "with the finite matrices through the documented Discrete and Continuous form. The continuous ",
        "operation constrains the local states but does not move atoms or generate bond orbits."
      ],
      StringJoin[
        "By default initfromrep finds the first 10 complete bond shells. If a model needs a farther ",
        "neighbour hopping, choose a larger InitialBondShells value and run initfromrep again; symham ",
        "does not search for new bonds."
      ],
      StringJoin[
        "Every call to initfromrep starts a new model. If the ordered group and representation data are ",
        "inconsistent, initialization fails and MagneticTB does not continue with the preceding model."
      ]
    },
    "OptionTable" -> {
      {
        "lattice",
        "{{a,0,0},{-a/2,Sqrt[3] a/2,0},{0,0,c}}",
        "Primitive direct-lattice vectors stored by rows. The numerical lattice obtained after lattpar must be real and nonsingular."
      },
      {
        "lattpar",
        "{a->1,c->3}",
        "Replacement rules assigning every symbolic lattice parameter used by lattice."
      },
      {
        "wyckoffposition",
        "{{{2/3,1/3,0},{0,0,1/2}}}",
        "One fractional seed and magnetic moment per Wyckoff orbit; the orbit order must match repinformation and orbitalLabels."
      },
      {
        "symminformation",
        "{{\"1\",IdentityMatrix[3],{0,0,0},\"F\"}}",
        "A complete, already closed, ordered finite magnetic or spin-space group. One additional pure-internal C_infty Association is allowed."
      },
      {
        "repinformation",
        "{{{1}}}",
        "Exact representation matrices in exactly the finite-group order. Use the documented Discrete/Continuous Association for C_infty."
      },
      {
        "orbitalLabels",
        "Automatic",
        "Names used to display the coordinates of the supplied representation. Automatic uses orb1, orb2, and so on; these labels do not change the matrices."
      },
      {
        "InitialBondShells",
        "10",
        "Number of complete bond shells prepared during initialization. Shell 1 contains onsite terms; shell 2 is the first nonzero-distance neighbour shell."
      },
      {
        "RepresentationMode",
        "\"DirectProduct\"",
        "Use \"DirectProduct\" when repinformation already gives the full local action, or \"Induced\" when it gives the representation of a reference site's symmetry group."
      }
    },
    "Examples" -> {
      {
        StringJoin[
          "Initialize the two-dimensional px/py representation of the tetragonal ",
          "{E,C4z,C2z,C4z^-1} group. The outer list ",
          "enumerates Wyckoff orbits and the inner list enumerates the complete ordered finite ",
          "group:"
        ],
        explicitInitCode
      },
      {
        "The user-supplied labels determine the displayed orbital order; no basis-function transformation is inferred:",
        "orbitalTable[]",
        MagneticTB`Private`renderOrbitalTableGrid[
          explicitRepresentationFixture["OrbitalTableData"]
        ]
      }
    },
    "OptionExamples" -> {
      {
        StringJoin[
          "lattice supplies the three tetragonal primitive vectors by rows. The same complete ",
          "initfromrep call also supplies its lattice parameters, Wyckoff seed, ordered C4 group, ",
          "exact matrices, and orbital labels. showbonds[2] displays the four in-plane nearest neighbours:"
        ],
        explicitInitCode <> "\nshowbonds[2]",
        explicitRepresentationFixture["BondGrid"]
      },
      {
        StringJoin[
          "lattpar -> {a->1,c->2} assigns the numerical tetragonal lattice constants before the exact ",
          "representation is installed. showbonds[2] then reports unit in-plane nearest-neighbour distance:"
        ],
        explicitInitCode <> "\nshowbonds[2]",
        explicitRepresentationFixture["BondGrid"]
      },
      {
        StringJoin[
          "wyckoffposition supplies one nonmagnetic seed at the origin. orbitalTable[] shows that the ",
          "two supplied representation coordinates px and py occupy this same physical site:"
        ],
        explicitInitCode <> "\norbitalTable[]",
        explicitRepresentationFixture["OrbitalTableGrid"]
      },
      {
        StringJoin[
          "symminformation supplies the complete ordered {E,C4z,C2z,C4z^-1} group. The representation ",
          "table shows all four operations together with their supplied two-dimensional orbital matrices:"
        ],
        explicitInitCode <> "\nshowSymmetryRepresentations[All]",
        explicitRepresentationFixture["RepresentationGrid"]
      },
      {
        StringJoin[
          "repinformation supplies the exact px/py matrices for E, C4z, C2z, and C4z^-1. ",
          "showSymmetryRepresentations[All] displays all four matrices in the completed model:"
        ],
        explicitInitCode <> "\nshowSymmetryRepresentations[All]",
        explicitRepresentationFixture["RepresentationGrid"]
      },
      {
        StringJoin[
          "orbitalLabels -> {{\"px\",\"py\"}} names the two coordinates of the supplied ",
          "representation. orbitalTable[] displays those labels on the model orbitals:"
        ],
        explicitInitCode <> "\norbitalTable[]",
        explicitRepresentationFixture["OrbitalTableGrid"]
      },
      {
        StringJoin[
          "InitialBondShells -> 2 prepares the onsite and nearest-neighbour shells for this complete ",
          "tetragonal model. showbonds[2] displays the four in-plane nearest neighbours:"
        ],
        explicitBondShellInitCode <> "\nshowbonds[2]",
        explicitRepresentationFixture["BondGrid"]
      },
      {
        StringJoin[
          "RepresentationMode -> \"DirectProduct\" keeps the supplied px and py coordinates on both ",
          "C4-related sites. orbitalTable[] therefore contains four Hamiltonian basis states:"
        ],
        explicitModeDirectInitCode <> "\norbitalTable[]",
        explicitRepresentationFixture["ModeDirectOrbitalTable"]
      },
      {
        StringJoin[
          "RepresentationMode -> \"Induced\" needs only the px representation of the reference-site ",
          "stabilizer {E,C2z}. C4 transports that one coordinate to the second site, so orbitalTable[] ",
          "contains only the two symmetry-required basis states and records C4z as the transport operation:"
        ],
        explicitModeInducedInitCode <> "\norbitalTable[]",
        explicitRepresentationFixture["ModeInducedOrbitalTable"]
      }
    },
    "Issues" -> {
      {
        "initfromrep deliberately rejects basisFunctions because exact ordered representation matrices, rather than transformed functions, define this input route:",
        "initfromrep[basisFunctions -> {{\"s\"}}]",
        $Failed,
        <|
          "ExpectedMessages" -> {"initfromrep::option"},
          "MessageText" ->
            "initfromrep::option: basisFunctions is not a supported initfromrep option."
        |>
      }
    },
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"symham", "MagneticTB/ref/symham"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"}
    },
    "Expected" -> {
      "InitialBondShells -> 2",
      "RepresentationMode -> \\\"Induced\\\""
    },
    "Keywords" -> {"explicit representation", "representation matrix", "continuous symmetry"}
  |>,
  <|
    "Name" -> "symham",
    "Usage" -> {
      {
        "symham[n]",
        "returns the symmetry-allowed contribution of bond shell n to the Bloch Hamiltonian."
      },
      {
        "symham[n, \"Hermitian\" -> False]",
        "constructs the same bond shell without imposing the Hermitian relation between a hopping and its reverse."
      }
    },
    "Notes" -> {
      StringJoin[
        "Call init first. For bond shell n, symham starts from the real-space hoppings in that shell, ",
        "relates symmetry-equivalent matrix elements, and keeps the independent amplitudes. Each hopping is ",
        "then multiplied by its Bloch phase and placed in the Hamiltonian."
      ],
      StringJoin[
        "Shell 1 contains onsite terms. Shell 2 is the first nonzero-distance neighbour shell, shell 3 is ",
        "the next one, and so forth. A complete model is obtained by adding the required shell contributions."
      ],
      StringJoin[
        "The result is Hermitian by default. The KernelMethod and ValidationLevel options control how the ",
        "common zero space is found and checked; they do not change the symmetry or orbitals of the model. ",
        "CartesianCoordinates -> False keeps reciprocal-fractional momentum components."
      ],
      "symham uses the atoms and bonds already found by init; it does not initialize the model again. If shell n was not included, rerun init with InitialBondShells -> n."
    },
    "OptionTable" -> {
      {
        "\"Hermitian\"",
        "True",
        "True imposes H_ij(R)=Conjugate[H_ji(-R)] and gives a Hermitian shell contribution. False treats the two directions as independent amplitudes."
      },
      {
        "\"KernelMethod\"",
        "\"Iterative\"",
        StringJoin[
          "Exact method used to find the common null space of the symmetry constraints. \"Iterative\" applies ",
          "the constraints successively; \"Stacked\" solves their vertically joined matrix. A failed method is ",
          "reported rather than replaced silently."
        ]
      },
      {
        "\"ValidationLevel\"",
        "\"Basic\"",
        StringJoin[
          "Amount of consistency checking after the null space is found. \"None\" skips these checks, ",
          "\"Basic\" checks the constraint residuals, and \"Full\" also checks independence and rank-nullity."
        ]
      },
      {
        "\"CartesianCoordinates\"",
        "False",
        "False uses reciprocal-fractional momentum components kx, ky, kz; True rewrites the result in Cartesian momentum coordinates."
      }
    },
    "Examples" -> {
      {
        "Initialize the complete graphene model used below. symham never performs initialization implicitly:",
        simpleInitCode
      },
      {
        "Shell 1 is the onsite contribution. The same onsite parameter appears on the two equivalent carbon sites:",
        "grapheneOnsite = symham[1];\nMatrixForm[grapheneOnsite]",
        MatrixForm[grapheneFixture["Onsite"]]
      },
      {
        "Shell 2 is the nearest-neighbour hopping. Its off-diagonal entries contain the three symmetry-related Bloch phases:",
        "grapheneNearest = symham[2];\nMatrixForm[grapheneNearest]",
        MatrixForm[grapheneFixture["Nearest"]]
      },
      {
        "Add onsite, nearest-neighbour, and next-nearest-neighbour shells to obtain the Bloch Hamiltonian used below:",
        grapheneHamiltonianCode,
        MatrixForm[grapheneFixture["ReadableHamiltonian"]]
      }
    },
    "OptionExamples" -> {
      {
        "CartesianCoordinates -> True rewrites the same shell in Cartesian momentum coordinates defined by the reciprocal lattice:",
        "grapheneCartesian = symham[2, \"CartesianCoordinates\" -> True];\nMatrixForm[grapheneCartesian]",
        MatrixForm@checkedDocumentationEvaluation[
          "GrapheneCartesianReference",
          MagneticTB`init[
            MagneticTB`lattice -> {
              {MagneticTB`a, 0, 0},
              {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
              {0, 0, MagneticTB`c}
            },
            MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 3},
            MagneticTB`wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
            MagneticTB`symminformation -> MagneticTB`msgop[MagneticTB`gray[191]],
            MagneticTB`basisFunctions -> {{"pz"}},
            MagneticTB`InitialBondShells -> 3,
            MagneticTB`GenerateSymmetryGroup -> False
          ];
          MagneticTB`symham[2, "CartesianCoordinates" -> True]
        ]
      },
      {
        StringJoin[
          "A chiral P3 model makes the Hermitian option visible. With Hermitian -> True the reverse ",
          "hoppings are fixed by t1 and t2; Hermitian -> False introduces independent reverse ",
          "amplitudes t3 and t4:"
        ],
        p3HermitianInitCode <> "\n" <> StringRiffle[{
          "p3Hermitian = symham[2];",
          "p3NonHermitian = symham[2, \"Hermitian\" -> False];",
          "Column[{",
          "  Row[{\"Hermitian -> True: \", MatrixForm[p3Hermitian]}],",
          "  Row[{\"Hermitian -> False: \", MatrixForm[p3NonHermitian]}]",
          "}]"
        }, "\n"],
        p3HermitianComparisonReference
      },
      {
        StringJoin[
          "KernelMethod -> \"Stacked\" solves the same exact symmetry constraints as one joined matrix. It ",
          "gives the same nearest-neighbour Hamiltonian as the default iterative ",
          "method:"
        ],
        "MatrixForm[symham[2, \"KernelMethod\" -> \"Stacked\"]]",
        MatrixForm[grapheneFixture["StackedNearest"]]
      },
      {
        "ValidationLevel -> \"Full\" additionally checks independence and rank-nullity after solving. These checks leave the exact nearest-neighbour Hamiltonian unchanged:",
        "MatrixForm[symham[2, \"ValidationLevel\" -> \"Full\"]]",
        MatrixForm[grapheneFixture["FullValidationNearest"]]
      }
    },
    "Issues" -> {
      {
        "A shell outside the range prepared by init fails; symham never extends the bond search implicitly:",
        "symham[4]",
        $Failed,
        <|
          "ExpectedMessages" -> {"symham::shellrange"},
          "MessageText" ->
            "symham::shellrange: Bond shell 4 was not prepared by init; only 3 shells are available. Rerun init with InitialBondShells -> 4."
        |>
      },
      {
        "The Hermitian option is strictly Boolean:",
        "symham[2, \"Hermitian\" -> Automatic]",
        $Failed,
        <|
          "ExpectedMessages" -> {"symham::hermitian"},
          "MessageText" ->
            "symham::hermitian: The Hermitian option must be True or False; received Automatic."
        |>
      }
    },
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"bandManipulate", "MagneticTB/ref/bandManipulate"},
      {"showbonds", "MagneticTB/ref/showbonds"}
    },
    "Expected" -> {
      "Hermitian -> False",
      "KernelMethod -> \\\"Stacked\\\"",
      "ValidationLevel -> \\\"Full\\\"",
      "CartesianCoordinates -> True"
    },
    "Keywords" -> {"Hamiltonian", "null space", "hopping", "Hermitian", "symmetry constraints"}
  |>,
  <|
    "Name" -> "bandManipulate",
    "Usage" -> {
      {
        "bandManipulate[path, npoint, H]",
        "draws the bands of H along path and provides sliders for the tight-binding parameters."
      }
    },
    "Notes" -> {
      "path is a list of reciprocal-fractional line segments and endpoint labels. npoint controls the sampling of each segment.",
      StringJoin[
        "Every symbol in H other than kx, ky, and kz becomes an interactive parameter. The sliders start at ",
        "zero, so the initial plot may show coincident flat bands when all Hamiltonian parameters vanish."
      ],
      StringJoin[
        "Use standardKPath[] as the first argument when MagneticTB should choose the built-in conventional ",
        "path from the current initialized lattice. bandManipulate receives the explicit path returned by ",
        "standardKPath[]."
      ],
      StringJoin[
        "No energy shift is added. For a genuinely non-Hermitian Hamiltonian, choose explicitly whether to ",
        "plot the real part, imaginary part, or complex spectrum before calling this function."
      ]
    },
    "Examples" -> {
      {
        "Initialize graphene and construct the Hamiltonian before plotting:",
        "Needs[\"MagneticTB`\"];\n" <>
          simpleInitCode <> "\n" <> grapheneHamiltonianCode,
        MatrixForm[grapheneFixture["ReadableHamiltonian"]]
      },
      {
        "Define the G-M-K-G path in reciprocal fractional coordinates:",
        graphenePathCode
      },
      {
        "Open the live interactive panel along the explicit G-M-K-G path. Its three sliders are the onsite, nearest-neighbour, and next-nearest-neighbour amplitudes; no hidden energy shift is added:",
        "bandManipulate[graphenePath, 30, grapheneHamiltonian]"
      },
      {
        "Display the conventional path selected automatically from the initialized graphene lattice:",
        "grapheneAutomaticPath = standardKPath[]",
        grapheneFixture["StandardKPath"]
      },
      {
        "Use the automatically selected path to open the live Manipulate panel:",
        "bandManipulate[grapheneAutomaticPath, 4, grapheneHamiltonian]",
        portableBandManipulateOutput[
          grapheneFixture["AutomaticBandUI"]
        ]
      },
      {
        "A fixed parameter choice gives the following static band plot, which is useful for printed documentation and direct comparison:",
        StringRiffle[{
          "bandplot[",
          "  graphenePath,",
          "  60,",
          "  grapheneHamiltonian,",
          "  {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}",
          "]"
        }, "\n"],
        grapheneFixture["BandPlot"]
      }
    },
    "Issues" -> {
      {
        "The Hamiltonian must be a nonempty square matrix. A rectangular array is rejected before the interactive panel is created:",
        "bandManipulate[graphenePath, 30, {{1, 2}}]",
        Failure[
          "InvalidHamiltonian",
          <|
            "InputHead" -> List,
            "NormalizedHead" -> List,
            "Dimensions" -> {1, 2}
          |>
        ],
        <|
          "ExpectedMessages" -> {"bandManipulate::matrix"},
          "MessageText" ->
            "bandManipulate::matrix: The Hamiltonian must be a nonempty square matrix after display wrappers are removed; received head List with dimensions {1, 2}."
        |>
      },
      {
        "Use ordinary symbols for adjustable parameters. Subscript and Indexed expressions are display constructs and are rejected explicitly:",
        "bandManipulate[graphenePath, 30, {{Subscript[p, 1]}}]",
        Failure[
          "UnsupportedBandParameter",
          <|"UnsupportedParameterExpressions" -> {Subscript[p, 1]}|>
        ],
        <|
          "ExpectedMessages" -> {"bandManipulate::parameter"},
          "MessageText" ->
            "bandManipulate::parameter: Subscript and Indexed expressions are not supported as band parameters. Use ordinary symbols instead: {Subscript[p, 1]}."
        |>
      },
      {
        "An unresolved user-defined function is not treated as a slider parameter. Replace it by an explicit matrix expression before calling bandManipulate:",
        "bandManipulate[graphenePath, 30, {{f[kx]}}]",
        Failure[
          "NonNumericHamiltonian",
          <|"UnsupportedFunctionHeads" -> {f}|>
        ],
        <|
          "ExpectedMessages" -> {"bandManipulate::function"},
          "MessageText" ->
            "bandManipulate::function: The Hamiltonian contains unsupported unresolved function heads: {f}."
        |>
      }
    },
    "SeeAlso" -> {
      {"symham", "MagneticTB/ref/symham"},
      {"init", "MagneticTB/ref/init"},
      {"standardKPath", "MagneticTB/ref/standardKPath"}
    },
    "Keywords" -> {"band structure", "interactive plot", "tight-binding parameters"}
  |>,
  <|
    "Name" -> "CurrentModelSession",
    "Usage" -> {
      {
        "CurrentModelSession[]",
        "returns the read-only data of the model most recently prepared by init or initfromrep."
      }
    },
    "Notes" -> {
      StringJoin[
        "Most users can inspect a model with orbitalTable, showbonds, and showSymmetryRepresentations. ",
        "CurrentModelSession is useful when a program needs all of that information together: the crystal ",
        "model, the full symmetry matrices, the bond shells found by init, and the shell results already ",
        "calculated by symham."
      ],
      "Calling CurrentModelSession before a successful initialization returns $Failed.",
      "Reading this value does not change the model, initialize another model, or solve a missing shell."
    },
    "Examples" -> {
      {
        "Initialize a complete model before inspecting it:",
        simpleInitCode
      },
      {
        StringJoin[
          "Read the representation mode, generated site orbits, local bases, and available shell numbers of ",
          "the graphene model. Merely viewing them does not calculate another ",
          "Hamiltonian:"
        ],
        StringRiffle[{
          "session = CurrentModelSession[];",
          "Grid[List @@@ Normal@<|",
          "  \"RepresentationSource\" -> session[\"RepresentationSource\"],",
          "  \"RepresentationMode\" -> session[\"RepresentationMode\"],",
          "  \"SiteOrbits\" -> session[\"ModelSpecification\", \"SiteOrbits\"],",
          "  \"LocalBases\" -> session[\"ModelSpecification\", \"LocalBases\"],",
          "  \"PreparedBondShells\" ->",
          "    Lookup[session[\"CompiledBondShells\"], \"Shell\"]",
          "|>, Frame -> All, Alignment -> Left]"
        }, "\n"],
        grapheneFixture["SessionSummary"]
      }
    },
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"initfromrep", "MagneticTB/ref/initfromrep"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"}
    },
    "Keywords" -> {"current model", "prepared bond shells", "representation data"}
  |>,
  <|
    "Name" -> "orbitalTable",
    "Usage" -> {
      {
        "orbitalTable[]",
        "lists the orbitals in the row and column order used by the Hamiltonian and symmetry matrices."
      }
    },
    "Notes" -> {
      "Each row identifies the Wyckoff orbit, the equivalent atom, its fractional and Cartesian coordinates, the local basis state, and any spin label.",
      "OrbitalID is the matrix index. The same order is used by symham, showHamiltonianBasis, and every full representation matrix.",
      "orbitalTable only displays the model prepared by init; it does not call symham."
    },
    "Examples" -> {
      {
        "Initialize a complete model before requesting its orbital table:",
        simpleInitCode
      },
      {
        "Display the complete Grid. OrbitalID is the row and column index used by both the Hamiltonian and every full representation matrix:",
        compactOrbitalTableCode,
        compactOrbitalTable[grapheneFixture["OrbitalTable"]]
      }
    },
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"CurrentModelSession", "MagneticTB/ref/CurrentModelSession"},
      {"symham", "MagneticTB/ref/symham"}
    },
    "Keywords" -> {"orbital order", "basis", "Wyckoff orbit", "Hamiltonian indices"}
  |>,
  <|
    "Name" -> "showbonds",
    "Usage" -> {
      {
        "showbonds[n]",
        "shows the bond length and all directed bonds in prepared shell n."
      }
    },
    "Notes" -> {
      "Shell 1 is the zero-distance onsite shell; shell 2 is the shortest nonzero-distance neighbour shell.",
      "The output identifies the starting site, the translated ending site, and the displacement of each directed bond.",
      "The shell must have been prepared by init. showbonds does not search farther and does not evaluate symham."
    },
    "Examples" -> {
      {
        "Initialize a complete model before displaying a prepared shell:",
        simpleInitCode
      },
      {
        "Display the first nonzero-distance shell. The table gives the source site, translated target site, displacement, and operation that produces every directed bond:",
        "showbonds[2]",
        grapheneFixture["BondGrid"]
      }
    },
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"symham", "MagneticTB/ref/symham"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"}
    },
    "Keywords" -> {"bond shell", "directed bond", "neighbour", "hopping geometry"}
  |>,
  <|
    "Name" -> "showSymmetryRepresentations",
    "Usage" -> {{
      "showSymmetryRepresentations[selection]",
      "shows how selected symmetry operations act on the complete Hamiltonian basis."
    }},
    "Notes" -> {
      "Automatic displays a generating set, while All displays every operation. An integer or list of integers selects operations by their order in symminformation.",
      StringJoin[
        "The table places the spatial operation, spin action, and unitary or antiunitary type beside its ",
        "matrix. Since antiunitarity is already part of the operation data, the matrix is shown without a ",
        "separate formal K."
      ],
      StringJoin[
        "Every matrix acts in the orbital order displayed by orbitalTable[]. Calling this function only ",
        "displays the matrices already determined by init or initfromrep; it does not solve a Hamiltonian."
      ]
    },
    "Examples" -> {
      {
        "Initialize graphene and display the matrices of four operations that generate its complete symmetry group:",
        simpleInitCode <> "\nshowSymmetryRepresentations[]",
        grapheneFixture["SymmetryRepresentations"]
      },
      {
        "Select a short ordered list containing unitary operations and pure time reversal:",
        "showSymmetryRepresentations[{1, 3, 25}]",
        grapheneFixture["SelectedSymmetryRepresentations"]
      },
      {
        "Select one operation by its position in symminformation. Operation 25 is pure time reversal in this ordered gray group:",
        "showSymmetryRepresentations[25]",
        grapheneFixture["TimeReversalRepresentation"]
      }
    },
    "Expected" -> {"showSymmetryRepresentations[{1, 3, 25}]"},
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"},
      {"showHamiltonianBasis", "MagneticTB/ref/showHamiltonianBasis"}
    },
    "Keywords" -> {"representation matrix", "antiunitary", "generator", "corepresentation"}
  |>,
  <|
    "Name" -> "showHamiltonianBasis",
    "Usage" -> {
      {
        "showHamiltonianBasis[]",
        "shows the basis states in the exact row and column order of the Hamiltonian."
      },
      {
        "showHamiltonianBasis[i,j]",
        "identifies the final bra state and initial ket state of matrix element H[[i,j]]."
      }
    },
    "Notes" -> {
      "The displayed order is identical to orbitalTable and to the full representation matrices.",
      "For H[[i,j]], column j is the initial ket and row i is the final bra. In Induced mode the table also shows how a reference-site state was transported to an equivalent site."
    },
    "Examples" -> {
      {
        "Initialize graphene and display the complete two-state Hamiltonian basis:",
        simpleInitCode <> "\nshowHamiltonianBasis[]",
        grapheneFixture["HamiltonianBasis"]
      },
      {
        "Identify the bra and ket states of the nearest-neighbour matrix element H[[1,2]]:",
        "showHamiltonianBasis[1, 2]",
        grapheneFixture["HamiltonianElementBasis"]
      }
    },
    "SeeAlso" -> {
      {"orbitalTable", "MagneticTB/ref/orbitalTable"},
      {"showHoppingParameters", "MagneticTB/ref/showHoppingParameters"},
      {"symham", "MagneticTB/ref/symham"}
    },
    "Keywords" -> {"Hamiltonian basis", "orbital order", "bra", "ket"}
  |>,
  <|
    "Name" -> "showHoppingParameters",
    "Usage" -> {
      {
        "showHoppingParameters[n]",
        "shows the representative real-space hopping associated with each independent parameter in a solved shell."
      },
      {
        "showHoppingParameters[n,p]",
        "shows every real-space matrix element generated from parameter p by symmetry and Hermiticity."
      }
    },
    "Notes" -> {
      "Evaluate symham[n] first. The table then connects a symbol such as t1 to a definite initial orbital, final orbital, and lattice translation.",
      "The two-argument form follows the same parameter through all symmetry-related bonds, making it possible to see where every term of the Bloch Hamiltonian comes from.",
      "This is an inspection function; it does not solve a missing shell."
    },
    "OptionTable" -> {
      {
        "\"Hermitian\"",
        "True",
        "Choose the earlier Hermitian or non-Hermitian shell calculation whose hopping origins are to be displayed."
      },
      {
        "\"KernelMethod\"",
        "\"Iterative\"",
        "Choose the earlier shell calculation made with the \"Iterative\" or \"Stacked\" null-space method. A missing calculation is not run automatically."
      },
      {
        "\"ValidationLevel\"",
        "\"Basic\"",
        "Choose the earlier shell calculation made with the corresponding \"None\", \"Basic\", or \"Full\" consistency check."
      }
    },
    "Examples" -> {
      {
        "Initialize graphene, solve its nearest-neighbour shell, and display which real-space hopping is represented by t1:",
        simpleInitCode <>
          "\ngrapheneNearest = symham[2];\nshowHoppingParameters[2]",
        grapheneFixture["HoppingParameters"]
      },
      {
        "Display every real-space occurrence generated from the representative t1 hopping:",
        "showHoppingParameters[2, t1]",
        grapheneFixture["HoppingOccurrences"]
      }
    },
    "OptionExamples" -> {
      {
        "\"Hermitian\" -> False selects a non-Hermitian shell with independent forward and reverse hoppings instead of identifying them by Hermiticity:",
        StringRiffle[{
          "symham[2, \"Hermitian\" -> False];",
          "showHoppingParameters[2, \"Hermitian\" -> False]"
        }, "\n"],
        grapheneFixture["NonHermitianHoppingParameters"]
      },
      {
        "\"KernelMethod\" -> \"Stacked\" selects the shell solved by the stacked exact null-space method. The same option is used when displaying its hopping origins:",
        StringRiffle[{
          "symham[2, \"KernelMethod\" -> \"Stacked\"];",
          "showHoppingParameters[2, \"KernelMethod\" -> \"Stacked\"]"
        }, "\n"],
        checkedDocumentationEvaluation[
          "StackedHoppingProvenance",
          MagneticTB`init[
            MagneticTB`lattice -> {
              {MagneticTB`a, 0, 0},
              {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
              {0, 0, MagneticTB`c}
            },
            MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 3},
            MagneticTB`wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
            MagneticTB`symminformation -> MagneticTB`msgop[MagneticTB`gray[191]],
            MagneticTB`basisFunctions -> {{"pz"}},
            MagneticTB`InitialBondShells -> 3,
            MagneticTB`GenerateSymmetryGroup -> False
          ];
          MagneticTB`symham[2, "KernelMethod" -> "Stacked"];
          MagneticTB`showHoppingParameters[
            2, "KernelMethod" -> "Stacked"
          ]
        ]
      },
      {
        "\"ValidationLevel\" -> \"Full\" selects the shell solved with the full residual, independence, and rank-nullity checks. The shell is solved explicitly before its origins are displayed:",
        StringRiffle[{
          "symham[2, \"ValidationLevel\" -> \"Full\"];",
          "showHoppingParameters[2, \"ValidationLevel\" -> \"Full\"]"
        }, "\n"],
        checkedDocumentationEvaluation[
          "FullHoppingProvenance",
          MagneticTB`init[
            MagneticTB`lattice -> {
              {MagneticTB`a, 0, 0},
              {-MagneticTB`a/2, Sqrt[3] MagneticTB`a/2, 0},
              {0, 0, MagneticTB`c}
            },
            MagneticTB`lattpar -> {MagneticTB`a -> 1, MagneticTB`c -> 3},
            MagneticTB`wyckoffposition -> {{{1/3, 2/3, 0}, {0, 0, 0}}},
            MagneticTB`symminformation -> MagneticTB`msgop[MagneticTB`gray[191]],
            MagneticTB`basisFunctions -> {{"pz"}},
            MagneticTB`InitialBondShells -> 3,
            MagneticTB`GenerateSymmetryGroup -> False
          ];
          MagneticTB`symham[
            2, "KernelMethod" -> "Iterative", "ValidationLevel" -> "Full"
          ];
          MagneticTB`showHoppingParameters[
            2, "KernelMethod" -> "Iterative", "ValidationLevel" -> "Full"
          ]
        ]
      }
    },
    "SeeAlso" -> {
      {"symham", "MagneticTB/ref/symham"},
      {"showHamiltonianBasis", "MagneticTB/ref/showHamiltonianBasis"},
      {"hop", "MagneticTB/ref/hop"}
    },
    "Keywords" -> {"hopping parameter", "real space", "representative bond"}
  |>,
  <|
    "Name" -> "bandplot",
    "Usage" -> {{
      "bandplot[path,npoint,H,rules]",
      "plots the numerical bands of H along a reciprocal-fractional momentum path after the model parameters are fixed."
    }},
    "Notes" -> {
      "Each path segment is given by two reciprocal-fractional endpoints and their labels. MagneticTB converts these coordinates to the phase convention of H.",
      "rules must assign numerical values to every tight-binding parameter. plotRange changes only the displayed energy window."
    },
    "OptionTable" -> {
      {
        "plotRange",
        "All",
        "All displays the complete numerical spectrum; {emin,emax} restricts only the visible energy interval."
      }
    },
    "Examples" -> {
      {
        "Construct the complete graphene Hamiltonian and plot a fixed numerical parameter choice:",
        simpleInitCode <> "\n" <> grapheneHamiltonianCode <> ";\n" <>
          graphenePathCode <> "\nbandplot[graphenePath, 60, grapheneHamiltonian, {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}]",
        grapheneFixture["BandPlot"]
      }
    },
    "OptionExamples" -> {{
      "plotRange changes only the displayed vertical energy interval; the Hamiltonian and sampled eigenvalues are unchanged:",
      simpleInitCode <> "\n" <> grapheneHamiltonianCode <> ";\n" <>
        graphenePathCode <>
        "\nbandplot[graphenePath, 60, grapheneHamiltonian, {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}, plotRange -> {-0.8, 0.8}]",
      grapheneFixture["BandPlotRange"]
    }},
    "SeeAlso" -> {
      {"bandManipulate", "MagneticTB/ref/bandManipulate"},
      {"banddata", "MagneticTB/ref/banddata"},
      {"symham", "MagneticTB/ref/symham"}
    },
    "Keywords" -> {"band structure", "static plot", "momentum path"}
  |>,
  <|
    "Name" -> "banddata",
    "Usage" -> {{
      "banddata[path,npoint,H,rules,file]",
      "writes the numerical band energies of H along path to a data file."
    }},
    "Notes" -> {
      StringJoin[
        "After applying rules, the function samples the same reciprocal-fractional path convention used by ",
        "bandplot. Each sampled momentum occupies one row of the file, and each band occupies one column."
      ],
      "Use FileNameJoin rather than a machine-specific literal path in notebooks that will be shared."
    },
    "Examples" -> {
      {
        "Initialize graphene, construct its Hamiltonian and reciprocal-space path, and choose a portable temporary CSV file:",
        simpleInitCode <>
          "\ngrapheneHamiltonian = Sum[symham[i], {i, 3}];\n" <>
          graphenePathCode <>
          "\ngrapheneBandFile = FileNameJoin[{$TemporaryDirectory, \"graphene-band.csv\"}];"
      },
      {
        "Export the sampled band energies. Each row of the CSV file is one momentum point and each column is one band:",
        StringJoin[
          "banddata[graphenePath, 10, grapheneHamiltonian, ",
          "{e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}, grapheneBandFile];"
        ]
      },
      {
        "Import the CSV file and plot the two graphene bands:",
        StringRiffle[{
          "grapheneBandEnergies = Import[grapheneBandFile, \"CSV\"];",
          "ListPlot[Transpose[grapheneBandEnergies]]"
        }, "\n"],
        grapheneFixture["BandDataPlot"]
      }
    },
    "SeeAlso" -> {
      {"bandplot", "MagneticTB/ref/bandplot"},
      {"bandManipulate", "MagneticTB/ref/bandManipulate"}
    },
    "Keywords" -> {"band data", "CSV", "export"}
  |>,
  <|
    "Name" -> "hop",
    "Usage" -> {
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
    "Notes" -> {
      StringJoin[
        "There are two ways to use hop. In hop[n,rules] and hop[{n1,n2,...},rules], first calculate the ",
        "required shells with symham. The integer form exports shells 1 through n, while the list form ",
        "exports only the listed shells. The real-space hopping matrices already obtained by symham are ",
        "written directly."
      ],
      "For hop[H,rules], H must be written as a finite sum of exponential Bloch phases so that integer lattice translations can be identified.",
      StringJoin[
        "The matrix form hop[H,rules] extracts the coefficients of the finite exponential phases in H. ",
        "The examples below obtain H from symham, so the model construction remains visible and reproducible. ",
        "Set \"hrExport\" to a path to write wannier90_hr.dat; with None, the text is shown in the notebook. ",
        "Use \"wcc\" -> Automatic when the orbital centers should come from init."
      ]
    },
    "OptionTable" -> {
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
    "Examples" -> {
      {
        "Initialize graphene, solve three shells, write a Wannier90 HR file, and inspect its header and first hopping entries:",
        simpleInitCode <>
          StringJoin[
            "\nScan[symham, {1, 2, 3}];\nhrFile = FileNameJoin[{$TemporaryDirectory, \"graphene_hr.dat\"}];",
            "\nhop[{1, 2, 3}, {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}, \"hrExport\" -> hrFile, \"wcc\" -> Automatic];",
            "\nTake[Import[hrFile, \"Lines\"], 12]"
          ],
        grapheneFixture["HRLines"]
      }
    },
    "OptionExamples" -> {
      {
        "\"Hermitian\" -> False exports the shell previously solved with independent directed hoppings. The same setting is supplied to symham and hop:",
        simpleInitCode <> "\n" <> StringRiffle[{
          "symham[2, \"Hermitian\" -> False];",
          "{FileNameTake[hop[",
          "  {2}, {t1 -> 0.5},",
          "  \"Hermitian\" -> False,",
          "  \"hrExport\" -> $TemporaryDirectory",
          "]]}"
        }, "\n"],
        {"wannier90_hr.dat"}
      },
      {
        "\"KernelMethod\" -> \"Stacked\" exports the shell solved by the stacked exact null-space method. The selector must match in symham and hop:",
        simpleInitCode <> "\n" <> StringRiffle[{
          "symham[2, \"KernelMethod\" -> \"Stacked\"];",
          "{FileNameTake[hop[",
          "  {2}, {t1 -> 0.5},",
          "  \"KernelMethod\" -> \"Stacked\",",
          "  \"hrExport\" -> $TemporaryDirectory",
          "]]}"
        }, "\n"],
        {"wannier90_hr.dat"}
      },
      {
        "\"ValidationLevel\" -> \"Full\" exports the shell solved with full exact validation. The selector must match in symham and hop:",
        simpleInitCode <> "\n" <> StringRiffle[{
          "symham[2, \"ValidationLevel\" -> \"Full\"];",
          "{FileNameTake[hop[",
          "  {2}, {t1 -> 0.5},",
          "  \"ValidationLevel\" -> \"Full\",",
          "  \"hrExport\" -> $TemporaryDirectory",
          "]]}"
        }, "\n"],
        {"wannier90_hr.dat"}
      },
      {
        "\"hrExport\" -> directory writes wannier90_hr.dat there instead of printing the text in the notebook:",
        simpleInitCode <> "\n" <> StringRiffle[{
          "symham[2];",
          "{FileNameTake[hop[",
          "  {2}, {t1 -> 0.5},",
          "  \"hrExport\" -> $TemporaryDirectory",
          "]]}"
        }, "\n"],
        {"wannier90_hr.dat"}
      },
      {
        "\"RealDigits\" -> 8 writes eight digits after the decimal point for both real and imaginary hopping components:",
        simpleInitCode <> "\n" <> StringRiffle[{
          "symham[2];",
          "realDigitsFile = hop[",
          "  {2}, {t1 -> 1/3},",
          "  \"hrExport\" -> FileNameJoin[",
          "    {$TemporaryDirectory, \"hop-real-digits_hr.dat\"}",
          "  ],",
          "  \"RealDigits\" -> 8",
          "];",
          "SelectFirst[",
          "  Import[realDigitsFile, \"Lines\"],",
          "  StringContainsQ[#, \"0.33333333\"] &",
          "]"
        }, "\n"],
        "-1 0 0 1 2 0.33333333 0.00000000"
      },
      {
        "\"wcc\" supplies one explicit fractional center per orbital when converting a Bloch Hamiltonian. This graphene Hamiltonian and both centers come from the model initialized above:",
        simpleInitCode <> "\n" <> StringRiffle[{
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
        }, "\n"],
        {{-1, 0, 0}, {0, -1, 0}, {0, 0, 0}, {0, 1, 0}, {1, 0, 0}}
      },
      {
        "\"TranslationTolerance\" -> 10^-6 accepts a 10^-8 numerical offset in the symmetry-generated tetragonal C4 chain phases and rounds them to the two nearest-neighbour cells:",
        chainInitCode <> "\n" <> StringRiffle[{
          "chainHamiltonian = symham[2] /. {",
          "  t1 -> 1,",
          "  kz -> (1 + 10^-8) kz",
          "};",
          "readHR[",
          "  hop[",
          "    chainHamiltonian, {},",
          "    \"hrExport\" -> FileNameJoin[",
          "      {$TemporaryDirectory, \"hop-translation-tolerance_hr.dat\"}",
          "    ],",
          "    \"TranslationTolerance\" -> 10^-6",
          "  ]",
          "][\"Translations\"]"
        }, "\n"],
        {{0, 0, -1}, {0, 0, 1}}
      }
    },
    "SeeAlso" -> {
      {"readHR", "MagneticTB/ref/readHR"},
      {"showHoppingParameters", "MagneticTB/ref/showHoppingParameters"},
      {"symham", "MagneticTB/ref/symham"}
    },
    "Keywords" -> {"Wannier90", "hr.dat", "real-space hopping", "export"}
  |>,
  <|
    "Name" -> "readHR",
    "Usage" -> {{
      "readHR[file]",
      "reads a Wannier90 HR file and returns its lattice translations and real-space hopping matrices."
    }},
    "Notes" -> {
      StringJoin[
        "Each lattice translation R in a Wannier90 HR file is paired with a real-space hopping matrix H(R). ",
        "Its row and column order is the Wannier-orbital order used when the file was written."
      ],
      StringJoin[
        "Before returning these matrices, readHR checks the file header, degeneracies, translation vectors, ",
        "orbital indices, and complex matrix elements. A malformed or incomplete file is rejected."
      ],
      StringJoin[
        "Pass the returned Association to buildBlochHamiltonian to recover H(k). For repeated evaluation, ",
        "define a function such as h[k_List] := buildBlochHamiltonian[hrData,k]."
      ],
      "The ncell option may then keep only translations inside a chosen real-space box."
    },
    "OptionTable" -> {
      {
        "\"prec\"",
        "10^-10",
        "Nonnegative numerical threshold used to Chop small real and imaginary matrix entries after parsing."
      },
      {
        "\"ncell\"",
        "All",
        "All keeps every translation; a nonnegative integer n keeps cells whose largest absolute translation component is at most n."
      }
    },
    "Examples" -> {
      {
        "Create a graphene HR file with hop, read it back, reconstruct H(k), and display the Hamiltonian at the Gamma point:",
        simpleInitCode <>
          StringJoin[
            "\nScan[symham, {1, 2, 3}];\nhrFile = FileNameJoin[{$TemporaryDirectory, \"graphene_hr.dat\"}];",
            "\nhop[{1, 2, 3}, {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}, \"hrExport\" -> hrFile, \"wcc\" -> Automatic];",
            "\nhrData = readHR[hrFile];",
            "\nhrHamiltonian[k_List] := buildBlochHamiltonian[hrData, k];",
            "\nMatrixForm@hrHamiltonian[{0, 0, 0}]"
          ],
        MatrixForm[grapheneFixture["HRBlochHamiltonian"]]
      }
    },
    "OptionExamples" -> {
      {
        "\"prec\" -> 10^-12 chops symmetry-generated hopping entries below that threshold. The 10^-13 non-onsite hopping becomes an exact zero:",
        StringRiffle[{
          "tinyHRFile = FileNameJoin[{$TemporaryDirectory, \"graphene_tiny_hr.dat\"}];",
          "hop[",
          "  {1,2,3}, {e1 -> 1., r1 -> 10^-13, t1 -> 10^-13},",
          "  \"hrExport\" -> tinyHRFile, \"wcc\" -> Automatic",
          "];",
          "tinyHRData = readHR[tinyHRFile, \"prec\" -> 10^-12];",
          "MatrixForm@First@Pick[",
          "  tinyHRData[\"HoppingMatrices\"],",
          "  tinyHRData[\"Translations\"],",
          "  Except[{0,0,0}]",
          "]"
        }, "\n"],
        MatrixForm[{{0}}]
      },
      {
        "\"ncell\" -> 0 keeps only the onsite translation after the complete HR file has been parsed and validated:",
        StringRiffle[{
          "readHR[hrFile, \"ncell\" -> 0][\"Translations\"]"
        }, "\n"],
        {{0, 0, 0}}
      }
    },
    "Expected" -> {
      "buildBlochHamiltonian",
      "symmetry-generated hopping entries"
    },
    "SeeAlso" -> {
      {"hop", "MagneticTB/ref/hop"},
      {"symham", "MagneticTB/ref/symham"}
    },
    "Keywords" -> {"Wannier90", "hr.dat", "import", "hopping matrices"}
  |>,
  <|
    "Name" -> "brokenSymmetryInitRules",
    "Usage" -> {{
      "brokenSymmetryInitRules[generators]",
      "builds a new set of init rules for the subgroup generated by selected symmetry operations."
    }},
    "Notes" -> {
      StringJoin[
        "Use this function when a term in the model preserves only part of the original symmetry. The ",
        "selected operations generate the symmetry subgroup that is to remain, while all atoms are kept."
      ],
      "If the subgroup no longer relates all atoms of a Wyckoff orbit, that orbit is split into separate entries of the new wyckoffposition input.",
      StringJoin[
        "The function only writes a new set of init rules; it does not change the model by itself. Check the ",
        "new Wyckoff positions and symmetry operations, then evaluate init@@rules. Use {} when only the ",
        "identity is to remain."
      ]
    },
    "Examples" -> {
      {
        "Initialize graphene and return the initialization rules for the subgroup generated by operation 2:",
        simpleInitCode <>
          "\nbrokenRules = brokenSymmetryInitRules[{2}];\nbrokenRules",
        grapheneFixture["BrokenSymmetryRules"]
      }
    },
    "SeeAlso" -> {
      {"init", "MagneticTB/ref/init"},
      {"msgop", "MagneticTB/ref/msgop"}
    },
    "Keywords" -> {"subgroup", "broken symmetry", "Wyckoff splitting"}
  |>,
  <|
    "Name" -> "showMSGWyckoff",
    "Usage" -> {{
      "showMSGWyckoff[group]",
      "shows the Wyckoff positions and symmetry-allowed magnetic moments of a magnetic space group."
    }},
    "Notes" -> {
      "group may be a database identifier such as gray[n] or the corresponding magnetic-group number.",
      StringJoin[
        "Each row gives a representative position. Symbols x, y, and z are free Wyckoff coordinates; the ",
        "moment entry gives the directions allowed by the site symmetry. The multiplicity is the number of ",
        "equivalent positions in the conventional cell."
      ]
    },
    "Examples" -> {
      {
        "Display the Wyckoff table of gray magnetic space group 193:",
        "showMSGWyckoff[gray[193]]",
        advancedModelFixture["WyckoffGrid"]
      }
    },
    "SeeAlso" -> {
      {"msgop", "MagneticTB/ref/msgop"},
      {"init", "MagneticTB/ref/init"}
    },
    "Keywords" -> {"Wyckoff position", "magnetic moment", "database"}
  |>,
  <|
    "Name" -> "mlgop",
    "Usage" -> {{
      "mlgop[group]",
      "reads the symmetry operations of a magnetic layer group from the bundled database."
    }},
    "Notes" -> {
      "Use graylayer[n] to select a gray magnetic layer group.",
      "As for msgop, each item contains the operation label, point matrix, fractional translation, and unitary or antiunitary flag."
    },
    "Examples" -> {
      {
        "Inspect the first four operations of gray layer group 25:",
        "Take[mlgop[graylayer[25]], 4]",
        advancedModelFixture["LayerOperations"]
      }
    },
    "SeeAlso" -> {
      {"mrgop", "MagneticTB/ref/mrgop"},
      {"msgop", "MagneticTB/ref/msgop"}
    },
    "Keywords" -> {"magnetic layer group", "symmetry database"}
  |>,
  <|
    "Name" -> "mrgop",
    "Usage" -> {{
      "mrgop[group]",
      "reads the symmetry operations of a magnetic rod group from the bundled database."
    }},
    "Notes" -> {
      "Use grayrod[n] to select a gray magnetic rod group.",
      "As for msgop, each item contains the operation label, point matrix, fractional translation, and unitary or antiunitary flag."
    },
    "Examples" -> {
      {
        "Inspect the first four operations of gray rod group 25:",
        "Take[mrgop[grayrod[25]], 4]",
        advancedModelFixture["RodOperations"]
      }
    },
    "SeeAlso" -> {
      {"mlgop", "MagneticTB/ref/mlgop"},
      {"msgop", "MagneticTB/ref/msgop"}
    },
    "Keywords" -> {"magnetic rod group", "symmetry database"}
  |>,
  <|
    "Name" -> "symhamII",
    "Usage" -> {{
      "symhamII[H]",
      "changes the Bloch basis convention of H by including the orbital-center phases."
    }},
    "Notes" -> {
      "In convention I, the Bloch phase is attached to the full displacement between orbitals. Convention II moves the intracell orbital-center phases into the basis states.",
      StringJoin[
        "Initialize the model first so that the fractional orbital centers are known, or supply them with ",
        "\"wcc\". An explicit list must follow the row and column order of H, and lets symhamII convert a ",
        "user-supplied Hamiltonian without relying on the current init session. The transformation changes ",
        "the matrix phases but not the energy spectrum."
      ]
    },
    "OptionTable" -> {{
      "\"wcc\"",
      "None",
      "None uses the orbital centers prepared by init. An explicit list supplies one fractional center per Hamiltonian orbital in the same order as the rows and columns of H."
    }},
    "Examples" -> {
      {
        "Construct graphene and display the same Hamiltonian in convention II:",
        simpleInitCode <> "\n" <> grapheneHamiltonianCode <>
          ";\nMatrixForm[symhamII[grapheneHamiltonian]]",
        MatrixForm[grapheneFixture["ConventionII"]]
      }
    },
    "OptionExamples" -> {{
      StringJoin[
        "\"wcc\" -> {{1/3,2/3,0},{2/3,1/3,0}} explicitly supplies the two graphene orbital ",
        "centers in the row and column order of grapheneHamiltonian. Use this form to convert a ",
        "user-supplied Hamiltonian that was not built by the current init session:"
      ],
      StringRiffle[{
        "grapheneCenters = {",
        "  {1/3, 2/3, 0},",
        "  {2/3, 1/3, 0}",
        "};",
        "MatrixForm[symhamII[",
        "  grapheneHamiltonian,",
        "  \"wcc\" -> grapheneCenters",
        "]]"
      }, "\n"],
      MatrixForm[grapheneFixture["ConventionII"]]
    }},
    "SeeAlso" -> {
      {"symham", "MagneticTB/ref/symham"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"}
    },
    "Keywords" -> {"Bloch convention", "Wannier center", "Hamiltonian phase"}
  |>,
  <|
    "Name" -> "initold",
    "Context" -> "MagneticTBOld`",
    "Initialization" -> "Needs[\"MagneticTBOld`\"]",
    "Tutorials" -> {{"Using MagneticTB 1.0", "MagneticTB/tutorial/LegacyBackend"}},
    "Usage" -> {
      {
        "initold[opts]",
        "sets up a model with the original MagneticTB implementation."
      }
    },
    "Notes" -> {
      "Use initold to run notebooks written for the original MagneticTB, or to compare the original and current constructions.",
      StringJoin[
        "Its physical inputs are the original lattice, Wyckoff, symmetry, and basis-function options, ",
        "renamed with the suffix old so that they cannot be confused with the new interface."
      ],
      "MagneticTBOld` and MagneticTB` do not load or call one another. Evaluate them in separate fresh Wolfram kernels when comparing results."
    },
    "OptionTable" -> {
      {
        "latticeold",
        "{{aold,0,0},{-aold/2,Sqrt[3] aold/2,0},{0,0,cold}}",
        "Primitive direct-lattice vectors used by the original version, stored by rows."
      },
      {
        "lattparold",
        "{aold->1,cold->3}",
        "Numerical substitutions for the symbolic lattice parameters in latticeold."
      },
      {
        "wyckoffpositionold",
        "{{{2/3,1/3,0},{0,0,1/2}}}",
        "Fractional Wyckoff seeds and magnetic moments used by the original version."
      },
      {
        "symminformationold",
        "{{\"1\",IdentityMatrix[3],{0,0,0},\"F\"}}",
        "Complete ordered magnetic-group operation records used by the original version."
      },
      {
        "basisFunctionsold",
        "{\"s\"}",
        "Local function basis, ordered exactly as wyckoffpositionold."
      },
      {
        "debugQold",
        "False",
        "Original Boolean diagnostic option, retained so that earlier notebooks still run."
      }
    },
    "Examples" -> {
      {
        "Initialize graphene with the original version and display its onsite Hamiltonian. Run this example in a kernel that loads only MagneticTBOld`:",
        oldInitCode <>
          "\nhOldOnsite = symhamold[1];\nMatrixForm[hOldOnsite]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      }
    },
    "OptionExamples" -> {
      {
        "latticeold supplies the original implementation's direct-lattice row vectors. This complete call uses the graphene lattice and displays the resulting onsite Hamiltonian:",
        oldInitCode <> "\nMatrixForm[symhamold[1]]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      },
      {
        "lattparold -> {aold->1,cold->3} evaluates the symbolic graphene lattice parameters before the original model is initialized:",
        oldInitCode <> "\nMatrixForm[symhamold[1]]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      },
      {
        "wyckoffpositionold supplies the original implementation's fractional seed and magnetic moment. The graphene seed generates the two-site onsite Hamiltonian:",
        oldInitCode <> "\nMatrixForm[symhamold[1]]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      },
      {
        "symminformationold supplies the complete ordered original-backend symmetry list returned by msgopold:",
        oldInitCode <> "\nMatrixForm[symhamold[1]]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      },
      {
        "basisFunctionsold -> {{\"pz\"}} assigns one original-backend pz orbital to the graphene Wyckoff seed:",
        oldInitCode <> "\nMatrixForm[symhamold[1]]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      },
      {
        "debugQold -> True uses the original Boolean diagnostic setting while leaving the onsite Hamiltonian unchanged:",
        oldDebugInitCode <> "\nMatrixForm[symhamold[1]]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      }
    },
    "SeeAlso" -> {
      {"symhamold", "MagneticTB/ref/symhamold"},
      {"bandManipulateold", "MagneticTB/ref/bandManipulateold"}
    },
    "Expected" -> {"debugQold -> False"},
    "Keywords" -> {"original version", "compatibility", "initialize model"}
  |>,
  <|
    "Name" -> "symhamold",
    "Context" -> "MagneticTBOld`",
    "Initialization" -> "Needs[\"MagneticTBOld`\"]",
    "Tutorials" -> {{"Using MagneticTB 1.0", "MagneticTB/tutorial/LegacyBackend"}},
    "Usage" -> {
      {
        "symhamold[n]",
        "returns bond shell n of the symmetry-allowed Hamiltonian obtained with the original MagneticTB algorithm."
      }
    },
    "Notes" -> {
      "Call initold first. The shell numbering and Bloch Hamiltonian follow the original package, while every symbol carries the old suffix where needed for isolation.",
      "symhamold is kept for existing notebooks and for independent comparisons with the current version. Its calculation is entirely the one used by the original MagneticTB."
    },
    "OptionTable" -> {
      {
        "symmetrysetold",
        "genindex",
        "Subset of symmetry-operation indices used by the original algorithm. This option exists only in MagneticTBOld`."
      },
      {
        "\"CartesianCoordinates\"",
        "False",
        "False uses the reciprocal-fractional momenta of the original version; True rewrites the result with its reciprocal lattice."
      }
    },
    "Examples" -> {
      {
        "Initialize the complete original model used in this example. symhamold does not perform initialization:",
        oldInitCode <>
          "\nhOldOnsite = symhamold[1];\nMatrixForm[hOldOnsite]",
        MatrixForm[legacyGrapheneFixture["Onsite"]]
      },
      {
        "Add the onsite and nearest-neighbour terms to obtain the Bloch Hamiltonian of the original version:",
        "hOldNearest = symhamold[2];\nhOld = hOldOnsite + hOldNearest;\nMatrixForm[hOld]",
        MatrixForm[legacyGrapheneFixture["Hamiltonian"]]
      }
    },
    "OptionExamples" -> {
      {
        "symmetrysetold -> All uses every original-backend symmetry operation supplied to initold when constructing the nearest-neighbour shell:",
        "MatrixForm[symhamold[2, symmetrysetold -> All]]",
        MatrixForm[legacyGrapheneFixture["Nearest"]]
      },
      {
        "\"CartesianCoordinates\" -> True writes the original nearest-neighbour term in Cartesian momentum coordinates:",
        "MatrixForm[symhamold[2, \"CartesianCoordinates\" -> True]]",
        MatrixForm[legacyGrapheneFixture["CartesianNearest"]]
      }
    },
    "SeeAlso" -> {
      {"initold", "MagneticTB/ref/initold"},
      {"bandManipulateold", "MagneticTB/ref/bandManipulateold"}
    },
    "Keywords" -> {"original version", "Hamiltonian", "compatibility"}
  |>,
  <|
    "Name" -> "bandManipulateold",
    "Context" -> "MagneticTBOld`",
    "Initialization" -> "Needs[\"MagneticTBOld`\"]",
    "Tutorials" -> {{"Using MagneticTB 1.0", "MagneticTB/tutorial/LegacyBackend"}},
    "Usage" -> {
      {
        "bandManipulateold[path, npoint, H]",
        "draws the bands of a Hamiltonian obtained with the original version and provides sliders for its parameters."
      }
    },
    "Notes" -> {
      "This function belongs to MagneticTBOld` and uses kxold, kyold, and kzold. Construct H with symhamold in the same fresh kernel before plotting it."
    },
    "Examples" -> {
      {
        "Initialize the original graphene model and construct its Hamiltonian before plotting:",
        oldInitCode <>
          "\nhOld = Total[symhamold /@ {1, 2}];\nMatrixForm[hOld]",
        MatrixForm[legacyGrapheneFixture["Hamiltonian"]]
      },
      {
        "Define the graphene path and open the interactive band panel of the original version:",
        StringRiffle[{
          "oldPath = {",
          "  {{{0, 0, 0}, {0, 1/2, 0}}, {\"G\", \"M\"}},",
          "  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
          "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"G\"}}",
          "};",
          "bandManipulateold[oldPath, 30, hOld]"
        }, "\n"],
        portableBandManipulateOutput[
          legacyGrapheneFixture["BandUI"]
        ]
      },
      {
        "Fix the original parameters to obtain a static band plot:",
        StringRiffle[{
          "bandplotold[",
          "  oldPath,",
          "  60,",
          "  hOld,",
          "  {e1old -> 0.05, t1old -> 0.5}",
          "]"
        }, "\n"],
        legacyGrapheneFixture["BandPlot"]
      }
    },
    "SeeAlso" -> {
      {"initold", "MagneticTB/ref/initold"},
      {"symhamold", "MagneticTB/ref/symhamold"}
    },
    "Keywords" -> {"original version", "band structure", "interactive plot"}
  |>
};

(* -------------------------------------------------------------------------
   Phase 2: merge testing-line API specifications.

   TestingPublicAPIDocumentation.wl owns the general testing API surface.
   TestingPlotPropertyDocumentation.wl is separately loadable so the plotting
   and property pages can use the small incremental gate above.
   ------------------------------------------------------------------------- *)

testingDocumentationSource = FileNameJoin[{
  scriptDirectory,
  "Examples",
  "TestingPublicAPIDocumentation.wl"
}];
If[!FileExistsQ[testingDocumentationSource],
  Print["Missing testing public-API documentation source: ",
    testingDocumentationSource];
  Exit[9]
];
Get[testingDocumentationSource];
testingBaseFunctionSpecifications = testingFunctionSpecifications;
plotPropertyDocumentationSource = FileNameJoin[{
  scriptDirectory,
  "Examples",
  "TestingPlotPropertyDocumentation.wl"
}];
If[!FileExistsQ[plotPropertyDocumentationSource],
  Print["Missing plotting/property documentation source: ",
    plotPropertyDocumentationSource];
  Exit[13]
];
Get[plotPropertyDocumentationSource];
testingFunctionSpecifications = Join[
  Select[
    testingBaseFunctionSpecifications,
    #["Name"] =!= "buildRealSpaceHamiltonian" &
  ],
  testingPlotPropertySpecifications
];
testingRealSpaceTutorialPage = testingApplyPlotPropertyTutorialUpdate[
  testingRealSpaceTutorialPage
];
If[Head[testingRealSpaceTutorialPage] =!= Notebook,
  Print["Could not update the generated real-space/topology tutorial."];
  Exit[16]
];
If[
  !ListQ[testingFunctionSpecifications] ||
    Length[testingFunctionSpecifications] =!= 32,
  Print["The testing public-API documentation specification is incomplete."];
  Exit[10]
];
testingDocumentationNames = Lookup[
  testingFunctionSpecifications,
  "Name",
  Missing["Name"]
];
If[
  !And @@ (AssociationQ /@ testingFunctionSpecifications) ||
    Sort[testingDocumentationNames] =!= Sort[{
      "bandplot", "showband", "standardKPath", "showCrystalStructure",
      "showBrillouinZone", "orbitalTable", "hoppingData",
      "transformHoppings", "buildBlochHamiltonian",
      "buildRealSpaceHamiltonian", "buildSlabHamiltonian",
      "surfaceGreenFunction", "surfaceSpectralFunction", "wilsonLoop",
      "berryPhase", "berryCurvature", "pointChernNumber",
      "findGaplessPoints", "getMSGElemFromMSGCorep", "getTBBandCorep",
      "graylayer", "grayrod", "hop", "vaspEig", "bandManipulateEig",
      "fittingTB", "compareBand",
      "plotWilsonLoop", "plotSurfaceSpectrum", "plotBerryCurvature2D",
      "plotBerryCurvature3D", "plotRealSpaceWavefunction"
    }],
  Print[InputForm@<|
    "TestingDocumentationNames" -> testingDocumentationNames,
    "Heads" -> (Head /@ testingFunctionSpecifications)
  |>];
  Exit[11]
];

(* Keep the complete authoring specifications above, but omit deferred pages
   from active source generation until their names are removed from the list. *)
testingFunctionSpecifications = Select[
  testingFunctionSpecifications,
  activeDocumentationSpecificationQ
];

(* A few stable pages have updated testing-line documentation contracts.
   Replace those pages and append every newly documented testing-line function
   without creating duplicate URIs. *)
functionSpecifications = Join[
  Select[
    functionSpecifications,
    !MemberQ[
      {
        "bandplot", "orbitalTable", "graylayer", "grayrod", "hop",
        "compareBand"
      },
      #["Name"]
    ] &
  ],
  testingFunctionSpecifications
];
functionSpecifications = Select[
  functionSpecifications,
  activeDocumentationSpecificationQ
];
testingDocumentationPageHeads = Head /@ (
  functionPage /@ testingFunctionSpecifications
);
If[!And @@ (# === Notebook & /@ testingDocumentationPageHeads),
  Print[InputForm@<|
    "TestingDocumentationNames" -> testingDocumentationNames,
    "PageHeads" -> testingDocumentationPageHeads
  |>];
  Exit[12]
];

functionSpecifications = Map[
  Function[spec,
    Switch[spec["Name"],
      "CompileMagneticTBInput",
        Join[spec, <|
          "Notes" -> Append[
            spec["Notes"],
            StringJoin[
              "When the program-generated input omits InitialBondShells, the testing-line default is 10 complete ",
              "shells, matching init and initfromrep. An explicit value in the input Association always takes ",
              "precedence."
            ]
          ]
        |>],
      _, spec
    ]
  ],
  functionSpecifications
];

(* Reused complete physical inputs for the main tutorial and the validated
   applications page.  Keep one authoring string per model so the two pages
   cannot silently drift apart. *)
csclApplicationCode = StringRiffle[{
  "csclOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, a}",
  "  },",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}},",
  "    {{1/2, 1/2, 1/2}, {0, 0, 0}}",
  "  },",
  "  symminformation -> csclOperations,",
  "  basisFunctions -> {{\"s\"}, {\"px\", \"py\", \"pz\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "csclHamiltonian = Sum[symham[i], {i, 2}];"
}, "\n"];

collinearApplicationCode = StringRiffle[{
  "identity3 = IdentityMatrix[3];",
  "c4z = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};",
  "collinearSSG = <|",
  "  \"C4spin\" -> <|",
  "    \"space\" -> {identity3, {0, 0, 0}},",
  "    \"spin\" -> {c4z, 0}",
  "  |>,",
  "  \"tauHalfT\" -> <|",
  "    \"space\" -> {identity3, {1/2, 0, 0}},",
  "    \"spin\" -> {identity3, 1}",
  "  |>",
  "|>;",
  "init[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 1}}},",
  "  symminformation -> collinearSSG,",
  "  basisFunctions -> {{\"sup\", \"sdn\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> True",
  "];",
  "collinearHamiltonian = Sum[symham[i], {i, 2}];"
}, "\n"];

validationBandPathCode = StringRiffle[{
  "validationBandPath = {",
  "  {{{0, 0, 0}, {1/2, 0, 0}}, {\"G\", \"X\"}},",
  "  {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {\"X\", \"M\"}},",
  "  {{{1/2, 1/2, 0}, {0, 0, 0}}, {\"M\", \"G\"}}",
  "};"
}, "\n"];

(* -------------------------------------------------------------------------
   Phase 3: guide and tutorial content.

   guideEntryCell keeps the guide list readable.  Tutorial Input/Output pairs
   still come from evaluated fixtures, never from hand-written output text.
   ------------------------------------------------------------------------- *)

guidePage = Notebook[{
  Cell["MagneticTB", "GuideTitle"],
  Cell[
    StringJoin[
      "MagneticTB constructs analytical tight-binding Hamiltonians from crystal symmetry and local orbital ",
      "information. It supports ordinary space groups, magnetic space groups, and finite discrete ",
      "spin-space groups."
    ],
    "GuideAbstract"
  ],
  Cell[CellGroupData[{
    Cell["Tutorials: start here", "GuideTutorialsSection"],
    guideTutorialCells[
      "Getting Started with MagneticTB",
      "MagneticTB/tutorial/GettingStarted",
      "Install, update, and load the package, then build a graphene Hamiltonian and plot its bands."
    ],
    guideTutorialCells[
      "General MagneticTB Examples",
      "MagneticTB/tutorial/GeneralExamples",
      "Follow complete model examples from symmetry operations to Hamiltonians and band structures."
    ],
    guideTutorialCells[
      "Crystal structure, Brillouin zone, and k paths",
      "MagneticTB/tutorial/CrystalAndKPaths",
      "Draw a noncollinear magnetic Kagome structure, its Brillouin zone, and its bands."
    ],
    guideTutorialCells[
      "Fitting MagneticTB bands to VASP",
      "MagneticTB/tutorial/BandFitting",
      "Fit VASP bands with nearest- and second-neighbor hopping, or select an energy window or k-point region."
    ],
    guideTutorialCells[
      "High-symmetry tetragonal antiferromagnet: spin symmetry and SOC",
      "MagneticTB/tutorial/ContinuousSymmetry",
      "Compare an antiferromagnet with SOC against the continuous spin symmetry required without SOC."
    ],
    guideTutorialCells[
      "Band corepresentations with MSGCorep",
      "MagneticTB/tutorial/Corepresentations",
      "Load MSGCorep explicitly and calculate single- and double-valued band corepresentations."
    ],
    guideTutorialCells[
      "Using MagneticTB 1.0",
      "MagneticTB/tutorial/LegacyBackend",
      "Load version 1.0 in a separate kernel to run existing models."
    ]
  }, Open]],
  guideFunctionBlock[{
    guideSubsectionCell[
      "Model setup and inspection",
      {
  guideEntryCell[
    "CurrentModelSession",
    "MagneticTB/ref/CurrentModelSession",
    "read all information already available for the current model."
  ],
  guideEntryCell[
    "orbitalTable",
    "MagneticTB/ref/orbitalTable",
    "inspect the fixed Hamiltonian orbital order."
  ],
  guideEntryCell[
    "CompileMagneticTBInput",
    "MagneticTB/ref/CompileMagneticTBInput",
    "check program-generated model input without making it the current model."
  ]
      }
    ],
    guideSubsectionCell[
      "Crystal, symmetry, and group theory",
      {
  guideEntryCell[
    "msgop",
    "MagneticTB/ref/msgop",
    "obtain magnetic-space-group operations from the bundled database."
  ],
  guideEntryCell[
    "gray",
    "MagneticTB/ref/gray",
    "look up gray magnetic-space-group database identifiers."
  ],
  guideEntryCell[
    "bnsdict",
    "MagneticTB/ref/bnsdict",
    "look up magnetic-space-group database identifiers by complete BNS pair."
  ],
  guideEntryCell[
    "typeI",
    "MagneticTB/ref/typeI",
    "look up type-I magnetic-space-group identifiers by BNS pair."
  ],
  guideEntryCell[
    "typeIII",
    "MagneticTB/ref/typeIII",
    "look up type-III magnetic-space-group identifiers by BNS pair."
  ],
  guideEntryCell[
    "typeIV",
    "MagneticTB/ref/typeIV",
    "look up type-IV magnetic-space-group identifiers by BNS pair."
  ],
  guideEntryCell[
    "graylayer",
    "MagneticTB/ref/graylayer",
    "look up gray magnetic-layer-group database identifiers."
  ],
  guideEntryCell[
    "grayrod",
    "MagneticTB/ref/grayrod",
    "look up gray magnetic-rod-group database identifiers."
  ],
  guideEntryCell[
    "init",
    "MagneticTB/ref/init",
    "initialize a model from basis functions."
  ],
  guideEntryCell[
    "initfromrep",
    "MagneticTB/ref/initfromrep",
    "initialize a model from exact representation matrices."
  ],
  guideEntryCell[
    "showSymmetryRepresentations",
    "MagneticTB/ref/showSymmetryRepresentations",
    "display full matrices for selected symmetry operations."
  ],
  guideEntryCell[
    "showHamiltonianBasis",
    "MagneticTB/ref/showHamiltonianBasis",
    "display the Hamiltonian bra-ket order."
  ],
  guideEntryCell[
    "brokenSymmetryInitRules",
    "MagneticTB/ref/brokenSymmetryInitRules",
    "prepare explicit subgroup initialization rules."
  ],
  guideEntryCell[
    "pointMatrix",
    "MagneticTB/ref/pointMatrix",
    "construct exact local point-operation matrices."
  ],
  guideEntryCell[
    "GenerateGroup",
    "MagneticTB/ref/GenerateGroup",
    "generate a finite group from an explicit product."
  ],
  guideEntryCell[
    "getGenerator",
    "MagneticTB/ref/getGenerator",
    "find a generating subset of a complete finite group (2.x compatibility function)."
  ],
  guideEntryCell[
    "showMSGWyckoff",
    "MagneticTB/ref/showMSGWyckoff",
    "display magnetic Wyckoff data."
  ],
  guideEntryCell[
    "mlgop",
    "MagneticTB/ref/mlgop",
    "read magnetic layer-group operations."
  ],
  guideEntryCell[
    "mrgop",
    "MagneticTB/ref/mrgop",
    "read magnetic rod-group operations."
  ],
  guideEntryCell[
    "getMSGElemFromMSGCorep",
    "MagneticTB/ref/getMSGElemFromMSGCorep",
    "convert the ordered operations of an installed MSGCorep package."
  ],
  guideEntryCell[
    "getTBBandCorep",
    "MagneticTB/ref/getTBBandCorep",
    "print single- or double-valued band corepresentations at supplied k points."
  ]
      }
    ],
    guideSubsectionCell[
      "Hamiltonians and band calculations",
      {
  guideEntryCell[
    "symham",
    "MagneticTB/ref/symham",
    "obtain the symmetry-allowed Hamiltonian of one bond shell."
  ],
  guideEntryCell[
    "unsymham",
    "MagneticTB/ref/unsymham",
    "inspect an unconstrained directed-bond Hamiltonian."
  ],
  guideEntryCell[
    "symhamII",
    "MagneticTB/ref/symhamII",
    "change the Bloch phase convention."
  ],
  guideEntryCell[
    "hoppingData",
    "MagneticTB/ref/hoppingData",
    "read numerical real-space hopping blocks from solved shells."
  ],
  guideEntryCell[
    "transformHoppings",
    "MagneticTB/ref/transformHoppings",
    "transform hopping blocks to an exact integer cell."
  ],
  guideEntryCell[
    "banddata",
    "MagneticTB/ref/banddata",
    "export numerical band eigenvalues."
  ],
  guideEntryCell[
    "hop",
    "MagneticTB/ref/hop",
    "export solved shells or an exponential H(k) to wannier90_hr.dat."
  ],
  guideEntryCell[
    "readHR",
    "MagneticTB/ref/readHR",
    "read and validate a Wannier90 HR file."
  ],
  Nothing
      }
    ],
    guideSubsectionCell[
      "Band fitting and VASP comparison",
      {
  guideEntryCell[
    "vaspEig",
    "MagneticTB/ref/vaspEig",
    "read selected bands from a complete numeric VASP EIGENVAL file."
  ],
  guideEntryCell[
    "bandManipulateEig",
    "MagneticTB/ref/bandManipulateEig",
    "tune a symbolic tight-binding Hamiltonian interactively against reference bands."
  ],
  guideEntryCell[
    "fittingTB",
    "MagneticTB/ref/fittingTB",
    "fit an affine-linear tight-binding Hamiltonian to reference bands."
  ],
  guideEntryCell[
    "compareBand",
    "MagneticTB/ref/compareBand",
    "overlay fitted tight-binding bands with supplied reference bands."
  ]
      }
    ],
    guideSubsectionCell[
      "Plotting and reciprocal-space paths",
      {
  guideEntryCell[
    "showbonds",
    "MagneticTB/ref/showbonds",
    "display the length and geometry of a bond shell."
  ],
  guideEntryCell[
    "showHoppingParameters",
    "MagneticTB/ref/showHoppingParameters",
    "trace each parameter to real-space hoppings."
  ],
  guideEntryCell[
    "bandManipulate",
    "MagneticTB/ref/bandManipulate",
    "explore the band structure interactively."
  ],
  guideEntryCell[
    "bandplot",
    "MagneticTB/ref/bandplot",
    "make a static numerical band plot."
  ],
  guideEntryCell[
    "showband",
    "MagneticTB/ref/showband",
    "draw a static band plot with an automatic conventional path."
  ],
  guideEntryCell[
    "standardKPath",
    "MagneticTB/ref/standardKPath",
    "select a built-in conventional path from the current lattice metric."
  ],
  guideEntryCell[
    "showCrystalStructure",
    "MagneticTB/ref/showCrystalStructure",
    "draw atoms, primitive cells, and magnetic-moment arrows."
  ],
  guideEntryCell[
    "showBrillouinZone",
    "MagneticTB/ref/showBrillouinZone",
    "draw the first Brillouin zone and selected k path."
  ],
  If[
    MemberQ[deferredDocumentationPageNames, "plotSurfaceSpectrum"],
    Nothing,
    guideEntryCell[
      "plotSurfaceSpectrum",
      "MagneticTB/ref/plotSurfaceSpectrum",
      "plot semi-infinite surface spectral weight on a momentum-energy grid."
    ]
  ]
      }
    ],
    If[
      MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
      Nothing,
      guideSubsectionCell[
        "Physical properties and topology",
        {
    guideEntryCell[
      "Real-space Hamiltonians, surfaces, and topology",
      "MagneticTB/tutorial/RealSpaceAndTopology",
      "build real-space, surface, and topological workflows."
    ],
    Nothing
        }
      ]
    ]
  }],
  guideFunctionBlock[{
    guideSubsectionCell["MagneticTB 1.0 functions", {
      Cell[
        StringJoin[
          "The MagneticTB 1.0 functions are included under MagneticTBOld`, with old appended to their ",
          "names. They can run earlier notebooks and provide a direct comparison with version 2.0. Load the ",
          "1.0 and 2.0 versions in separate kernels."
        ],
        "GuideText"
      ],
      guideEntryCell[
        "initold", "MagneticTB/ref/initold",
        "initialize a model with version 1.0."
      ],
      guideEntryCell[
        "symhamold", "MagneticTB/ref/symhamold",
        "construct a Hamiltonian with version 1.0."
      ],
      guideEntryCell[
        "bandManipulateold", "MagneticTB/ref/bandManipulateold",
        "explore bands obtained from version 1.0."
      ]
    }]
  }],
  historyCells[],
  categorizationCells["Guide", "MagneticTB`", "MagneticTB/guide/MagneticTB"],
  keywordCells[{"tight binding", "magnetic symmetry", "Hamiltonian", "band structure"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "MagneticTB",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "GuidePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

tutorialPage = Notebook[{
  Cell["Getting Started with MagneticTB", "Title"],
  Cell[
    StringJoin[
      "To construct a tight-binding model, first specify the lattice, atomic positions, ",
      "symmetry operations, and local orbitals. With these inputs, MagneticTB finds the ",
      "equivalent atoms and the symmetry-allowed hopping terms. We use graphene to show how to ",
      "obtain the Hamiltonian and plot its bands."
    ],
    "Text"
  ],
  Cell["Install, update, uninstall, and load", "Section"],
  Cell[
    StringJoin[
      "Obtain a local MagneticTB .paclet file, either supplied to you or built from the ",
      "source. The installation commands below use this file; they do not download a release ",
      "from the Internet."
    ],
    "Text"
  ],
  Cell["Install", "Subsection"],
  gettingStartedHelpHomeCell[],
  Cell[
    StringJoin[
      "Install the .paclet file with the following command. Replace the path and filename with ",
      "those of your file:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "PacletInstall[\"/absolute/path/MagneticTB-2.x.y.paclet\"]"
  ]], "Input"],
  Cell[
    StringJoin[
      "After installation, quit Mathematica completely and open it again before loading ",
      "MagneticTB or opening its help. This also removes any definitions left by the ",
      "previously loaded version."
    ],
    "Text"
  ],
  Cell["Update", "Subsection"],
  gettingStartedHelpHomeCell[],
  Cell[
    StringJoin[
      "To update, install the newer .paclet file in the same way, then quit and reopen ",
      "Mathematica. If PacletInstall reports samevers, that version is already installed. To ",
      "replace it with another file of the same version, uninstall that version first and then ",
      "install the new file."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "PacletInstall[\"/absolute/path/MagneticTB-2.x.y.paclet\"]"
  ]], "Input"],
  Cell["Uninstall or reinstall the same version", "Subsection"],
  gettingStartedHelpHomeCell[],
  Cell[
    "To uninstall MagneticTB, run the following command, then quit and reopen Mathematica:",
    "Text"
  ],
  Cell[BoxData[codeBoxes["PacletUninstall[\"MagneticTB\"]"]], "Input"],
  Cell["Load and open the documentation", "Subsection"],
  gettingStartedHelpHomeCell[],
  Cell[
    StringJoin[
      "After reopening Mathematica, load the package with Needs. To open the help, search for ",
      "MagneticTB in the Documentation Center, or run the command below to open its home page. ",
      "The functions from version 1.0 are loaded with MagneticTBOld`; use a separate kernel ",
      "for that version."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["Needs[\"MagneticTB`\"]"]], "Input"],
  Cell[BoxData[codeBoxes[
    "SystemOpen[\"paclet:MagneticTB/guide/MagneticTB\"]"
  ]], "Input"],

  Cell["The five physical inputs", "Section"],
  Cell[
    StringJoin[
      "init sets up the crystal and its local orbitals. Here lattice gives the three lattice ",
      "vectors as rows, and lattpar supplies their parameter values. wyckoffposition gives one ",
      "representative position and its magnetic moment for each inequivalent Wyckoff orbit. ",
      "symminformation gives the symmetry operations, and basisFunctions lists the orbitals in ",
      "the same order as the positions."
    ],
    "Text"
  ],
  Cell[
    StringJoin[
      "The two carbon atoms in graphene are symmetry equivalent. Enter only {1/3,2/3,0}; the ",
      "other atom at {2/3,1/3,0} is generated by symmetry. Keep one pz orbital on each carbon. ",
      "The lattice parameter c separates the periodically repeated graphene layers."
    ],
    "Text"
  ],
  descriptionTableCell[{
    {
      "lattice",
      "length",
      "Three primitive vectors stored as rows. Coordinates in wyckoffposition are fractional coordinates in this basis."
    },
    {
      "lattpar",
      "length",
      "Numerical assignments for every symbolic lattice constant. MagneticTB does not impose an angstrom unit; all lengths must use one consistent unit."
    },
    {
      "wyckoffposition",
      "fractional position and magnetic moment",
      "The outer list enumerates inequivalent Wyckoff orbits. Each entry is {seed position, magnetic moment}; the same outer order is used by basisFunctions."
    },
    {
      "symminformation",
      "symmetry operations",
      "Each magnetic-group record is {label,R,t,\"F\"|\"T\"}. The last field records whether the operation is unitary or antiunitary."
    },
    {
      "basisFunctions",
      "local orbital basis",
      "One nonempty local basis list per Wyckoff orbit. Different entries may have different dimensions, which produces rectangular inter-orbit hopping blocks."
    }
  }],
  Cell[BoxData[codeBoxes[simpleInitCode]], "Input"],
  Cell["View atomic positions and magnetization directions", "Section"],
  Cell[
    StringJoin[
      "After init, use atompos to view the generated atomic positions and magnetization ",
      "directions. The atoms are grouped by Wyckoff orbit. Each entry is {fractional position, ",
      "magnetization direction}, with both vectors expressed in the lattice-vector basis:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["atompos"]], "Input"],
  Cell[BoxData[outputBoxes[grapheneFixture["AtomPositions"]]], "Output"],
  Cell[
    StringJoin[
      "The two entries are the carbon atoms at {1/3,2/3,0} and {2/3,1/3,0}. Their ",
      "magnetization directions are both {0,0,0}, since this model is nonmagnetic."
    ],
    "Text"
  ],
  Cell["Check the neighbour shells", "Section"],
  Cell[
    StringJoin[
      "Hoppings are numbered by distance. Number 1 is the onsite term, number 2 is the nearest ",
      "neighbor, and number 3 is the second neighbor. For graphene, showbonds[2] displays the ",
      "three nearest-neighbor carbon-carbon bonds and their reverse directions. No hopping ",
      "parameters need to be specified to draw these bonds."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["showbonds[2]"]], "Input"],
  Cell[BoxData[outputBoxes[grapheneFixture["BondGrid"]]], "Output"],

  Cell["Solve the Hamiltonian shell by shell", "Section"],
  Cell[
    StringJoin[
      "First calculate the onsite term with symham[1]. The two carbon atoms are equivalent, so ",
      "both diagonal entries contain the same real parameter e1:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "grapheneOnsite = symham[1];\nMatrixForm[grapheneOnsite]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[grapheneFixture["Onsite"]]]], "Output"],
  Cell[
    StringJoin[
      "Next calculate the nearest-neighbor hopping with symham[2]. The three phase factors in ",
      "the upper-right entry come from the three nearest-neighbor bonds. For real t1, the ",
      "lower-left entry is its complex conjugate:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "grapheneNearest = symham[2];\nMatrixForm[grapheneNearest]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[grapheneFixture["Nearest"]]]], "Output"],
  Cell[
    StringJoin[
      "symham[3] gives the second-neighbor hopping. These bonds connect carbon atoms on the ",
      "same sublattice, so their six phase factors appear in the diagonal entries:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "grapheneNextNearest = symham[3];\nMatrixForm[grapheneNextNearest]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[grapheneFixture["NextNearest"]]]], "Output"],
  Cell[
    StringJoin[
      "Add the contributions to obtain the Hamiltonian. Here we keep the onsite, ",
      "nearest-neighbor, and second-neighbor terms:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[grapheneHamiltonianCode]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[grapheneFixture["ReadableHamiltonian"]]]], "Output"],

  Cell["Explore the bands", "Section"],
  Cell[
    StringJoin[
      "The path is given in reciprocal fractional coordinates. Each entry contains a segment ",
      "and its endpoint labels. bandManipulate draws the bands along this path and provides a ",
      "slider for each parameter. The sliders initially equal zero, and no extra energy shift ",
      "is added."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[graphenePathCode]], "Input"],
  Cell[BoxData[codeBoxes[
    "bandManipulate[graphenePath, 30, grapheneHamiltonian]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[grapheneFixture["BandUI"]]], "Output"],
  Cell[
    StringJoin[
      "To draw a band plot with fixed parameters, assign their values as below. The ",
      "second-neighbor term r1 breaks particle-hole symmetry, but the Dirac crossing at K ",
      "remains:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "bandplot[",
    "  graphenePath,",
    "  60,",
    "  grapheneHamiltonian,",
    "  {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}",
    "]"
  }, "\n"]]], "Input"],
  Cell[BoxData[outputBoxes[grapheneFixture["BandPlot"]]], "Output"],
  Cell[
    StringJoin[
      "The two bands meet at K. plotRange only changes the displayed energy range; it does not ",
      "change the energies or the position of the crossing."
    ],
    "Text"
  ],

  Cell["Multiple orbitals and multiple Wyckoff positions", "Section"],
  Cell[
    StringJoin[
      "For several inequivalent Wyckoff positions, list one representative of each in ",
      "wyckoffposition and place its orbitals in the corresponding entry of basisFunctions. In ",
      "the CsCl example below, the corner site carries s, and the body-center site carries px, ",
      "py, and pz. The full Hamiltonian is shown in the order {s,px,py,pz}, including the ",
      "onsite terms and hopping between the two sites."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    StringReplace[csclApplicationCode, {
      ",\n  InitialBondShells -> 2" -> "",
      ",\n  GenerateSymmetryGroup -> False" -> ""
    }] <> "\n" <> validationBandPathCode
  ]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[validationModelFixture["CsCl", "OrbitalTable"]]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    "csclHamiltonian = FullSimplify[ExpToTrig[csclHamiltonian]];\nMatrixForm[csclHamiltonian]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[
    FullSimplify[ExpToTrig[validationModelFixture["CsCl", "Hamiltonian"]]]
  ]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, csclHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CsCl", "BandUI"]]], "Output"],

  Cell["A spinful antiunitary model", "Section"],
  Cell[
    StringJoin[
      "A spin-space group allows the spin rotation to be specified separately from the spatial ",
      "operation. The following example has a pure spin C4 rotation and a half translation ",
      "combined with time reversal. In each dictionary entry, spin -> {S,0} denotes a unitary ",
      "operation and spin -> {S,1} an antiunitary one; S is the spin-rotation matrix."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    StringReplace[collinearApplicationCode,
      "  InitialBondShells -> 2,\n" -> ""
    ] <> "\n" <> validationBandPathCode
  ]], "Input"],
  Cell[
    "Use atompos again to view the two atoms and their magnetization directions:",
    "Text"
  ],
  Cell[BoxData[codeBoxes["atompos"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CollinearSSG", "AtomPositions"]]], "Output"],
  Cell[
    StringJoin[
      "The moments at {0,0,0} and {1/2,0,0} point along +z and -z. Use showCrystalStructure to ",
      "draw the atoms and their opposite moments:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["showCrystalStructure[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CollinearSSG", "CrystalStructure"]]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[collinearHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[compactMatrixOutput[
    MatrixForm[validationModelFixture["CollinearSSG", "Hamiltonian"]]
  ]]], "Output"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CollinearSSG", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, collinearHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CollinearSSG", "BandUI"]]], "Output"],
  Cell[
    StringJoin[
      "DirectProduct uses the symmetry action on the orbitals supplied at each site. Induced ",
      "starts from local states at one representative site and generates the symmetry-related ",
      "states. For example, in the cubic p-orbital model, DirectProduct requires px, py, and ",
      "pz explicitly, whereas Induced can start from px at the representative site. The ",
      "operation relating two sites may also include time reversal."
    ],
    "Text"
  ],

  Cell["Related functions and tutorials", "Section"],
  linkListingCell[
    "Text",
    {
      {"init", "MagneticTB/ref/init"},
      {"symham", "MagneticTB/ref/symham"},
      {"orbitalTable", "MagneticTB/ref/orbitalTable"},
      {"showbonds", "MagneticTB/ref/showbonds"},
      {"bandManipulate", "MagneticTB/ref/bandManipulate"},
      {"General Examples", "MagneticTB/tutorial/GeneralExamples"}
    }
  ],
  historyCells[],
  categorizationCells["Tech Note", "MagneticTB`", "MagneticTB/tutorial/GettingStarted"],
  keywordCells[{"getting started", "graphene", "tight binding", "Hamiltonian", "band structure"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Getting Started with MagneticTB",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

moS2InitCode = StringRiffle[{
  "moS2Operations = msgop[gray[187]];",
  "moS2BasisTransform = {{1, -1, 0}, {0, 1, 0}, {0, 0, 1}};",
  "moS2Operations = MapAt[",
  "  Function[operation,",
  "    FullSimplify[",
  "      moS2BasisTransform.operation.",
  "        Inverse[moS2BasisTransform]",
  "    ]",
  "  ],",
  "  moS2Operations,",
  "  {All, 2}",
  "];",
  "init[",
  "  lattice -> {",
  "    {1, 0, 0},",
  "    {1/2, Sqrt[3]/2, 0},",
  "    {0, 0, 100}",
  "  },",
  "  lattpar -> {},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> moS2Operations,",
  "  basisFunctions -> {{\"dz2\", \"dxy\", \"dx2-y2\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "moS2Hamiltonian = (symham[1] + symham[2]) /. {",
  "  e1 -> e2, e2 -> e1,",
  "  t1 -> t6, t2 -> t5, t3 -> t4,",
  "  t4 -> t3, t5 -> t2, t6 -> t1",
  "};"
}, "\n"];

ClearAll[moS2SpinOrbitCode];
moS2SpinOrbitCode[labels_List, prefix_String] := Module[{labelCode},
  labelCode = "{" <>
    StringRiffle[("\"" <> # <> "\"") & /@ labels, ", "] <>
    "}";
  StringRiffle[{
    prefix <> "Operations = msgop[gray[187]];",
    "init[",
    "  lattice -> {",
    "    {a, 0, 0},",
    "    {-a/2, Sqrt[3] a/2, 0},",
    "    {0, 0, c}",
    "  },",
    "  lattpar -> {a -> 1, c -> 4},",
    "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
    "  symminformation -> " <> prefix <> "Operations,",
    "  basisFunctions -> {" <> labelCode <> "},",
    "  InitialBondShells -> 1,",
    "  GenerateSymmetryGroup -> False",
    "];",
    prefix <> "Rules = brokenSymmetryInitRules[{9, 11, 13}];",
    "init @@ " <> prefix <> "Rules;",
    prefix <> "Hamiltonian = symham[1];"
  }, "\n"]
];

moS2PathCode = StringRiffle[{
  "moS2Path = {",
  "  {{{0, 0, 0}, {2/3, 1/3, 0}}, {\"G\", \"K\"}},",
  "  {{{2/3, 1/3, 0}, {1/2, 1/2, 0}}, {\"K\", \"M\"}},",
  "  {{{1/2, 1/2, 0}, {0, 0, 0}}, {\"M\", \"G\"}}",
  "};",
  "moS2Parameters = {",
  "  e1 -> 1.046, e2 -> 2.104,",
  "  t1 -> -0.184, t2 -> 0.401, t3 -> 0.218,",
  "  t4 -> 0.507, t5 -> 0.338, t6 -> 0.057",
  "};"
}, "\n"];

c3WeylInitCode = StringRiffle[{
  "c3Operations = msgop[bnsdict[{143, 3}]];",
  "init[",
  "  lattice -> {",
  "    {Sqrt[3] a/2, -a/2, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 1}}",
  "  },",
  "  symminformation -> c3Operations,",
  "  basisFunctions -> {{{x + I y, 0}, {0, x - I y}}},",
  "  InitialBondShells -> 4,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "c3Full = Sum[symham[i], {i, {2, 3, 4}}];",
  "c3OriginalParameterCoordinates = {",
  "  r11 -> -r2, r15 -> -r10,",
  "  r3 -> -r6, r7 -> -r14,",
  "  s4 -> -s4,",
  "  t4 -> t1, t8 -> t9, t11 -> -t7",
  "};",
  "c3Hamiltonian =",
  "  c3Full[[{1, 4}, {1, 4}]] /.",
  "    c3OriginalParameterCoordinates;"
}, "\n"];

c3PathCode = StringRiffle[{
  "c3Path = {",
  "  {{{0, 0, 0}, {1/2, 0, 0}}, {\"G\", \"M\"}},",
  "  {{{1/2, 0, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
  "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"G\"}},",
  "  {{{0, 0, 0}, {0, 0, 1/2}}, {\"G\", \"A\"}},",
  "  {{{0, 0, 1/2}, {1/2, 0, 1/2}}, {\"A\", \"L\"}},",
  "  {{{1/2, 0, 1/2}, {1/3, 1/3, 1/2}}, {\"L\", \"H\"}},",
  "  {{{1/3, 1/3, 1/2}, {0, 0, 1/2}}, {\"H\", \"A\"}},",
  "  {{{1/2, 0, 1/2}, {1/2, 0, 0}}, {\"L\", \"M\"}},",
  "  {{{1/3, 1/3, 0}, {1/3, 1/3, 1/2}}, {\"K\", \"H\"}}",
  "};",
  "c3Parameters = {",
  "  r10 -> -0.09, r14 -> 0.25, r2 -> 0, r6 -> 0.5,",
  "  s4 -> -0.38, s8 -> 0, t1 -> 0,",
  "  t15 -> 0.735, t7 -> 0.05, t9 -> 0.015",
  "};"
}, "\n"];

nodalLineInitCode = StringRiffle[{
  "nodalOperations = msgop[bnsdict[{184, 196}]];",
  "init[",
  "  lattice -> {",
  "    {Sqrt[3] a/2, -a/2, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 1}}",
  "  },",
  "  symminformation -> nodalOperations,",
  "  basisFunctions -> {{{x + I y, 0}, {0, x - I y}}},",
  "  InitialBondShells -> 5,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "nodalFull = Sum[symham[i], {i, 5}];",
  "nodalOriginalParameterCoordinates = {",
  "  p5n1 -> p5n4, p5n8 -> p5n2,",
  "  s1 -> s3, t1 -> t2,",
  "  p5n6 -> -p5n7, r2 -> -r1",
  "};",
  "nodalHamiltonian =",
  "  nodalFull[[{1, 4}, {1, 4}]] /.",
  "    nodalOriginalParameterCoordinates;"
}, "\n"];

c4InitCode = StringRiffle[{
  "c4Operations = msgop[gray[75]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 1},",
  "  wyckoffposition -> {",
  "    {{0, 0, 1/4}, {0, 0, 0}},",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> c4Operations,",
  "  basisFunctions -> {{\"px\", \"py\"}, {\"px\", \"py\"}},",
  "  InitialBondShells -> 7,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "c4Raw = Association@Table[",
  "  shell -> symham[shell],",
  "  {shell, {2, 3, 4, 5, 7}}",
  "];"
}, "\n"];

generalExamplesPage = Notebook[{
  Cell["General MagneticTB Examples", "Title"],
  Cell[
    StringJoin[
      "The following examples show how to construct models for graphene, MoS2, magnetic Weyl ",
      "points, nodal lines, and a C4 topological insulator. For each model, first specify the ",
      "lattice, atomic positions, symmetry, and orbitals. Then calculate the Hamiltonian, ",
      "assign its parameters, and plot the bands. Some examples also calculate a Chern number ",
      "or export Wannier90 hopping data."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["Needs[\"MagneticTB`\"]"]], "Input"],
  Cell["Graphene", "Section"],
  Cell[
    StringJoin[
      "Start with graphene, retaining one pz orbital on each carbon atom. The complete input ",
      "below gives the two-band Hamiltonian with onsite, nearest-neighbor, and second-neighbor ",
      "terms:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    grapheneInitCode <> "\n" <> grapheneHamiltonianCode
  ]], "Input"],
  Cell[
    BoxData[outputBoxes[MatrixForm[grapheneFixture["ReadableHamiltonian"]]]],
    "Output"
  ],
  Cell[BoxData[codeBoxes[graphenePathCode]], "Input"],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "bandplot[",
    "  graphenePath, 60, grapheneHamiltonian,",
    "  {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}",
    "]"
  }, "\n"]]], "Input"],
  Cell[BoxData[outputBoxes[grapheneFixture["BandPlot"]]], "Output"],

  Cell["Spinful graphene, interactive SOC bands, and Wannier90 export", "Section"],
  Cell[
    StringJoin[
      "To include spin-orbit coupling, replace pz by {pzup,pzdn} on each carbon atom. There ",
      "are now four basis states. The gray group includes time reversal, and symham gives the ",
      "spin-dependent terms allowed by this symmetry. We use this Hamiltonian below to plot ",
      "bands and export Wannier90 data."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[spinfulGrapheneInitCode]], "Input"],
  Cell[BoxData[codeBoxes[
    "MatrixForm[spinfulGrapheneHamiltonian]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[
    advancedModelFixture["SpinfulGraphene", "Hamiltonian"]
  ]]], "Output"],
  Cell[
    StringJoin[
      "Run the following command to display the parameter sliders. Move them to see how the ",
      "bands change:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "bandManipulate[spinfulGraphenePath, 30, spinfulGrapheneHamiltonian]"
  ]], "Input"],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "bandplot[",
    "  spinfulGraphenePath,",
    "  60,",
    "  spinfulGrapheneHamiltonian,",
    "  spinfulGrapheneParameters",
    "]"
  }, "\n"]]], "Input"],
  Cell[BoxData[outputBoxes[
    advancedModelFixture["SpinfulGraphene", "BandPlot"]
  ]], "Output"],
  Cell[
    StringJoin[
      "After symham has obtained the hopping terms, use hop to write them to a Wannier90 HR ",
      "file. The beginning of the file is shown below:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "spinfulGrapheneHR = FileNameJoin[{",
    "  $TemporaryDirectory,",
    "  \"spinful-graphene_hr.dat\"",
    "}];",
    "hop[",
    "  {1, 2, 3},",
    "  spinfulGrapheneParameters,",
    "  \"hrExport\" -> spinfulGrapheneHR",
    "];",
    "Take[Import[spinfulGrapheneHR, \"Lines\"], 12]"
  }, "\n"]]], "Input"],
  Cell[BoxData[outputBoxes[
    advancedModelFixture["SpinfulGraphene", "HRLines"]
  ]], "Output"],
  Cell[
    StringJoin[
      "The header gives the number of Wannier orbitals and translation vectors. Each hopping ",
      "entry then gives a translation vector, two orbital indices, and the real and imaginary ",
      "parts of the matrix element."
    ],
    "Text"
  ],

  Cell["Three-band MoS2 model", "Section"],
  Cell[
    StringJoin[
      "For monolayer MoS2, keep the three local orbitals dz2, dxy, and dx2-y2. The ",
      "substitution after symham converts the parameters to the convention used in the ",
      "MagneticTB 1.0 example. Assign moS2Parameters only after this conversion, so each value ",
      "corresponds to the intended matrix element."
    ],
    "Text"
  ],
  Cell[
    StringJoin[
      "Reference: G.-B. Liu, W.-Y. Shan, Y. Yao, W. Yao, and D. Xiao, \"Three-band ",
      "tight-binding model for monolayers of group-VIB transition metal dichalcogenides,\" ",
      "Physical Review B 88, 085433 (2013), https://doi.org/10.1103/PhysRevB.88.085433. The ",
      "example below uses the three-band model and parameter values from the MagneticTB 1.0 ",
      "example; the parameter conversion does not change its numerical Hamiltonian."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[moS2InitCode]], "Input"],
  Cell[BoxData[codeBoxes[moS2PathCode]], "Input"],
  Cell[BoxData[codeBoxes[
    "bandplot[moS2Path, 60, moS2Hamiltonian, moS2Parameters]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[advancedModelFixture["MoS2", "BandPlot"]]], "Output"],
  Cell["Spinful MoS2 and two explicit basis orders", "Subsection"],
  Cell[
    StringJoin[
      "For the spinful onsite Hamiltonian, first place all spin-up orbitals before the ",
      "spin-down orbitals. Specify the reduced symmetry in init before calculating the ",
      "Hamiltonian:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[moS2SpinOrbitCode[
    {
      "dz2up", "dxyup", "dx2-y2up",
      "dz2dn", "dxydn", "dx2-y2dn"
    },
    "moS2SpinA"
  ]]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[
      advancedModelFixture["MoS2SpinOrbitA", "OrbitalTable"]
    ]
  ]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[moS2SpinAHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[
    advancedModelFixture["MoS2SpinOrbitA", "Hamiltonian"]
  ]]], "Output"],
  Cell[
    StringJoin[
      "Alternatively, place the spin-up and spin-down states next to each other for each ",
      "orbital. The order shown by orbitalTable is also the row and column order of the ",
      "Hamiltonian:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[moS2SpinOrbitCode[
    {
      "dz2up", "dz2dn", "dxyup",
      "dxydn", "dx2-y2up", "dx2-y2dn"
    },
    "moS2SpinB"
  ]]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[
      advancedModelFixture["MoS2SpinOrbitB", "OrbitalTable"]
    ]
  ]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[moS2SpinBHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[
    advancedModelFixture["MoS2SpinOrbitB", "Hamiltonian"]
  ]]], "Output"],

  Cell["Magnetic C3 Weyl point and Chern diagnostic", "Section"],
  Cell[
    StringJoin[
      "For magnetic space group 143.3, the local spinor pair gives a four-band Hamiltonian. As ",
      "in the MagneticTB 1.0 example, select the two-band block with row and column indices ",
      "{1,4}. The replacement c3OriginalParameterCoordinates puts its parameters in the 1.0 ",
      "convention before numerical values are assigned."
    ],
    "Text"
  ],
  Cell[
    StringJoin[
      "The parameters in this Weyl model and the following nodal-line and C4 models illustrate ",
      "symmetry-allowed band structures. They are the values used in the MagneticTB 1.0 ",
      "examples, not fitted parameters for a particular material."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[c3WeylInitCode]], "Input"],
  Cell[
    StringJoin[
      "The full two-band c3Hamiltonian is used for the band plot and Chern-number calculation ",
      "below. To see its matrix elements more clearly, first display it on Gamma-A by setting ",
      "kx=ky=0:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "c3GammaAxisHamiltonian = FullSimplify[",
    "  c3Hamiltonian /. {kx -> 0, ky -> 0},",
    "  Element[kz, Reals]",
    "];",
    "MatrixForm[c3GammaAxisHamiltonian]"
  }, "\n"]]], "Input"],
  Cell[
    BoxData[outputBoxes[MatrixForm[
      advancedModelFixture["C3Weyl", "GammaAxisHamiltonian"]
    ]]],
    "Output"
  ],
  Cell[BoxData[codeBoxes[c3PathCode]], "Input"],
  Cell[
    "Use the following command to adjust the parameters and view the bands:",
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    "bandManipulate[c3Path, 20, c3Hamiltonian]"
  ]], "Input"],
  Cell[BoxData[codeBoxes[
    "bandplot[c3Path, 200, c3Hamiltonian, c3Parameters]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[advancedModelFixture["C3Weyl", "BandPlot"]]], "Output"],
  Cell[
    StringJoin[
      "Expand the Hamiltonian about Gamma through total order three, then pass it to ",
      "pointChernNumber. The function calculates the occupied-band Berry flux through a cube ",
      "surrounding the point. The occupied band must remain separated from the other band on ",
      "the cube surface; the enclosed Chern number is given below."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "c3Taylor = Expand@Normal@Chop@Series[",
    "  c3Hamiltonian, {ky, 0, 3}, {kx, 0, 3}, {kz, 0, 3}",
    "];",
    "c3Taylor = FromCoefficientRules[",
    "  Map[",
    "    Select[#, Total[First[#]] <= 3 &] & ,",
    "    CoefficientRules[c3Taylor, {kx, ky, kz}],",
    "    {2}",
    "  ],",
    "  {kx, ky, kz}",
    "];",
    "c3Numeric[momentum_List] := N[",
    "  c3Taylor /. c3Parameters /.",
    "    Thread[{kx, ky, kz} -> momentum]",
    "];"
  }, "\n"]]], "Input"],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "pointChernNumber[",
    "  c3Numeric, 1, {0, 0, 0}, 2,",
    "  \"SurfaceSubdivisions\" -> 8",
    "]"
  }, "\n"]]], "Input"],
  Cell[
    BoxData[outputBoxes[advancedModelFixture["C3Weyl", "ChernNumber"]]],
    "Output"
  ],

  Cell["Magnetic cubic nodal line", "Section"],
  Cell[
    StringJoin[
      "Next use magnetic space group 184.196 with the same local spinor basis. The operations ",
      "numbered {2,7,13} in the MagneticTB 1.0 example generate this 24-operation group; init ",
      "below receives the full group. Again select the {1,4} two-band block, and use ",
      "nodalOriginalParameterCoordinates to put its parameters in the 1.0 convention."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[nodalLineInitCode]], "Input"],
  Cell[BoxData[codeBoxes[
    "MatrixForm[nodalHamiltonian]"
  ]], "Input"],
  Cell[
    BoxData[outputBoxes[MatrixForm[
      advancedModelFixture["NodalLine", "Hamiltonian"]
    ]]],
    "Output"
  ],
  Cell[
    "Run the last line to display the sliders and plot the bands:",
    "Text"
  ],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "nodalPath = {",
    "  {{{0, 1/10, 1/4}, {0, 0, 1/4}}, {\"Q\", \"P\"}},",
    "  {{{0, 0, 1/4}, {-1/10, 1/15, 1/4}}, {\"P\", \"Q'\"}}",
    "};",
    "bandManipulate[nodalPath, 100, nodalHamiltonian]"
  }, "\n"]]], "Input"],

  Cell["C4 topological-insulator model", "Section"],
  Cell[
    StringJoin[
      "This model has two inequivalent Wyckoff positions and two local orbitals. Keep the ",
      "indicated neighbor hoppings and make the parameter substitutions shown below, as in the ",
      "MagneticTB 1.0 example:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[c4InitCode]], "Input"],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "c4Hamiltonian = FullSimplify[",
    "  (c4Raw[2] /. {t1 -> 0, t2 -> t1p}) +",
    "  (c4Raw[3] /. {r1 -> 0, r2 -> tzp}) +",
    "  (c4Raw[4] /. {",
    "    s10 -> 0, s7 -> 0, s9 -> 0, s2 -> 0,",
    "    s4 -> 0, s1 -> 0, s6 -> 0, s5 -> 0,",
    "    s8 -> t1a, s3 -> t1b",
    "  }) +",
    "  ((c4Raw[5] /. {p5n2 -> 0, p5n3 -> 0, p5n1 -> p5n4}) /.",
    "    p5n4 -> t2p) +",
    "  ((c4Raw[7] /. {",
    "    p7n7 -> p7n5, p7n1 -> 0, p7n4 -> 0,",
    "    p7n2 -> 0, p7n3 -> 0, p7n14 -> p7n12,",
    "    p7n11 -> 0, p7n8 -> 0, p7n10 -> 0,",
    "    p7n9 -> 0, p7n6 -> -p7n5, p7n13 -> p7n12",
    "  }) /. {p7n5 -> t2a/2, p7n12 -> t2b/2})",
    "];",
    "MatrixForm[symhamII[c4Hamiltonian]]"
  }, "\n"]]], "Input"],
  Cell[
    BoxData[outputBoxes[
      MatrixForm[advancedModelFixture["C4TI", "ConventionII"]]
    ]],
    "Output"
  ],
  Cell[BoxData[codeBoxes[StringRiffle[{
    "c4Path = {",
    "  {{{0, 0, 0}, {1/2, 1/2, 0}}, {\"G\", \"M\"}},",
    "  {{{1/2, 1/2, 0}, {1/2, 1/2, 1/2}}, {\"M\", \"A\"}},",
    "  {{{1/2, 1/2, 1/2}, {0, 0, 1/2}}, {\"A\", \"Z\"}},",
    "  {{{0, 0, 1/2}, {0, 0, 0}}, {\"Z\", \"G\"}}",
    "};",
    "bandplot[c4Path, 50, c4Hamiltonian, {",
    "  t1a -> 1, t1b -> -1, t2a -> 0.5, t2b -> -0.5,",
    "  t1p -> 2.5, t2p -> 0.5, tzp -> 2",
    "}]"
  }, "\n"]]], "Input"],
  Cell[BoxData[outputBoxes[advancedModelFixture["C4TI", "BandPlot"]]], "Output"],

  Cell["Wyckoff, layer-group, and rod-group databases", "Section"],
  Cell[
    StringJoin[
      "The symmetry databases can also be used on their own. showMSGWyckoff lists the ",
      "representative positions and allowed moments of a magnetic space group. mlgop and mrgop ",
      "give the operations of magnetic layer groups and magnetic rod groups, respectively."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["showMSGWyckoff[gray[193]]"]], "Input"],
  Cell[BoxData[outputBoxes[advancedModelFixture["WyckoffGrid"]]], "Output"],
  Cell[BoxData[codeBoxes[
    "showMSGWyckoff[bnsdict[{187, 211}]]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    advancedModelFixture["WyckoffGridBNS187211"]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    "Take[mlgop[graylayer[25]], 4]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[advancedModelFixture["LayerOperations"]]], "Output"],
  Cell[BoxData[codeBoxes[
    "Take[mlgop[{25, 3, 145}], 4]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    advancedModelFixture["LayerOperationsNumbered"]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    "Take[mrgop[grayrod[25]], 4]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[advancedModelFixture["RodOperations"]]], "Output"],
  Cell[BoxData[codeBoxes[
    "Take[mrgop[{25, 3, 130}], 4]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    advancedModelFixture["RodOperationsNumbered"]
  ]], "Output"],

  Cell["Changing the symmetry in version 2.0", "Section"],
  Cell[
    StringJoin[
      "To lower the symmetry, change the model inputs and run init again before calculating a ",
      "new Hamiltonian. brokenSymmetryInitRules gives the retained subgroup and the resulting ",
      "Wyckoff positions. Inspect these rules and pass them to init@@rules. The symmetry is ",
      "set by init, not by a symmetryset option in symham."
    ],
    "Text"
  ],
  historyCells[],
  categorizationCells[
    "Tech Note",
    "MagneticTB`",
    "MagneticTB/tutorial/GeneralExamples"
  ],
  keywordCells[{
    "Graphene", "MoS2", "Weyl point", "Chern number",
    "nodal line", "topological insulator", "magnetic space group"
  }]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "General MagneticTB Examples",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

legacyTutorialPage = Notebook[{
  Cell["Using MagneticTB 1.0", "Title"],
  Cell[
    StringJoin[
      "MagneticTB 2.0 includes the functions from version 1.0. To continue using a 1.0 ",
      "notebook, load MagneticTBOld` and use the function names ending in old. The graphene ",
      "example below shows how to initialize a model, calculate its Hamiltonian, and plot the ",
      "bands with these functions."
    ],
    "Text"
  ],
  Cell["Load MagneticTB 1.0", "Section"],
  Cell[
    StringJoin[
      "Start a fresh kernel and load MagneticTBOld` with the following command. Use separate ",
      "kernels for versions 1.0 and 2.0:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes["Needs[\"MagneticTBOld`\"]"]], "Input"],
  Cell["Initialize the model with version 1.0", "Section"],
  Cell[
    StringJoin[
      "For graphene, give the hexagonal lattice, one representative carbon position, and its ",
      "pz orbital. The symmetry operations generate the second carbon atom. In version 1.0, ",
      "the function and option names below all end in old:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[oldInitCode]], "Input"],
  Cell["Construct and display the Hamiltonian with version 1.0", "Section"],
  Cell[BoxData[codeBoxes[
    "hOldOnsite = symhamold[1];\nhOldNearest = symhamold[2];\nhOld = hOldOnsite + hOldNearest;\nMatrixForm[hOld]"
  ]], "Input"],
  Cell[
    BoxData[outputBoxes[MatrixForm[legacyGrapheneFixture["Hamiltonian"]]]],
    "Output"
  ],
  Cell["Plot the bands with version 1.0", "Section"],
  Cell[
    BoxData[codeBoxes[
      StringRiffle[{
        "oldPath = {",
        "  {{{0, 0, 0}, {0, 1/2, 0}}, {\"G\", \"M\"}},",
        "  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {\"M\", \"K\"}},",
        "  {{{1/3, 1/3, 0}, {0, 0, 0}}, {\"K\", \"G\"}}",
        "};",
        "bandManipulateold[oldPath, 30, hOld]"
      }, "\n"]
    ]],
    "Input"
  ],
  Cell[BoxData[outputBoxes[legacyGrapheneFixture["BandUI"]]], "Output"],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "bandplotold[",
      "  oldPath,",
      "  60,",
      "  hOld,",
      "  {e1old -> 0.05, t1old -> 0.5}",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[BoxData[outputBoxes[legacyGrapheneFixture["BandPlot"]]], "Output"],
  Cell["Compare results from versions 1.0 and 2.0", "Section"],
  Cell[
    StringJoin[
      "To compare versions 1.0 and 2.0, enter the same lattice, positions, symmetry, and ",
      "orbitals in two separate kernels. Then compare the Hamiltonians and bands. Each version ",
      "needs its own initialization; neither uses the model set up in the other kernel."
    ],
    "Text"
  ],
  Cell["Related functions", "Section"],
  linkListingCell[
    "Text",
    {
      {"initold", "MagneticTB/ref/initold"},
      {"symhamold", "MagneticTB/ref/symhamold"},
      {"bandManipulateold", "MagneticTB/ref/bandManipulateold"},
      {"MagneticTB", "MagneticTB/guide/MagneticTB"}
    }
  ],
  historyCells[],
  categorizationCells["Tech Note", "MagneticTBOld`", "MagneticTB/tutorial/LegacyBackend"],
  keywordCells[{"MagneticTB 1.0", "compatibility", "version comparison"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Using MagneticTB 1.0",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

continuousModelsExampleCode = Import[
  FileNameJoin[{scriptDirectory, "Examples", "ContinuousModels.wl"}],
  "Text",
  CharacterEncoding -> "UTF-8"
];
If[!StringQ[continuousModelsExampleCode],
  Print["Unable to read the continuous-symmetry documentation example."];
  Exit[3]
];
continuousSLGStart = First@First@StringPosition[
  continuousModelsExampleCode,
  "slgName ="
];
continuousDirectStart = First@First@StringPosition[
  continuousModelsExampleCode,
  "ClearAll[initializeSLGDirectProduct]"
];
continuousInducedStart = First@First@StringPosition[
  continuousModelsExampleCode,
  "ClearAll[initializeSLGInduced]"
];
continuousMinimalExampleCode = StringTrim@StringTake[
  continuousModelsExampleCode,
  continuousSLGStart - 1
];
continuousMinimalSetupCode = StringTrim@StringReplace[
  continuousMinimalExampleCode,
  "Needs[\"MagneticTB`\"];" -> ""
];
continuousSLGFixtureCode =
  "documentationSpinLayer14123L1Fixture = " <>
  ToString[
    spinLayerDocumentationFixture,
    InputForm,
    PageWidth -> 78
  ] <> ";";
continuousSLGPreparationCode = StringTrim@StringTake[
  continuousModelsExampleCode,
  {continuousSLGStart, continuousDirectStart - 1}
];
continuousSLGDirectCode = StringTrim@StringTake[
  continuousModelsExampleCode,
  {continuousDirectStart, continuousInducedStart - 1}
];
continuousSLGInducedCode = StringTrim@StringDrop[
  continuousModelsExampleCode,
  continuousInducedStart - 1
];

crystalSystemsExampleCode = Import[
  FileNameJoin[{scriptDirectory, "Examples", "ValidatedCrystalSystems.wl"}],
  "Text",
  CharacterEncoding -> "UTF-8"
];
If[!StringQ[crystalSystemsExampleCode],
  Print["Unable to read the seven-crystal-system documentation example."];
  Exit[3]
];

Get[FileNameJoin[{scriptDirectory, "Examples", "ValidatedCrystalSystems.wl"}]];
magneticTBPackageRoot = FileNameJoin[{projectRoot, "MagneticTB"}];
wyckoffDatabase = Import[
  FileNameJoin[{magneticTBPackageRoot, "Data", "wyckoffMSG.mx"}]
];
wyckoffDocumentationFixture = checkedDocumentationEvaluation[
  "TenWyckoffDatabaseCases",
  Table[
    <|
      "Case" -> case,
      "DirectProduct" -> runWyckoffMode[case, "DirectProduct"],
      "Induced" -> runWyckoffMode[case, "Induced"]
    |>,
    {case, wyckoffCases}
  ]
];

simpleCubicApplicationCode = StringRiffle[{
  "simpleCubicOperations = msgop[gray[221]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, a}",
  "  },",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}}",
  "  },",
  "  symminformation -> simpleCubicOperations,",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "simpleCubicHamiltonian = Sum[symham[i], {i, 2}];"
}, "\n"];

moS2ElevenBandApplicationCode = StringRiffle[{
  "moS2ElevenOperations = msgop[gray[187]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {-a/2, Sqrt[3] a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 4},",
  "  wyckoffposition -> {",
  "    {{0, 0, 0}, {0, 0, 0}},",
  "    {{1/3, 2/3, 1/4}, {0, 0, 0}}",
  "  },",
  "  symminformation -> moS2ElevenOperations,",
  "  basisFunctions -> {",
  "    {\"dz2\", \"dx2-y2\", \"dxy\", \"dxz\", \"dyz\"},",
  "    {\"px\", \"py\", \"pz\"}",
  "  },",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "moS2ElevenHamiltonian = Sum[symham[i], {i, 2}];"
}, "\n"];

feseApplicationCode = StringRiffle[{
  "feseOperations = msgop[gray[129]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {",
  "    {{3/4, 1/4, 0}, {0, 0, 0}},",
  "    {{1/4, 1/4, 267/1000}, {0, 0, 0}}",
  "  },",
  "  symminformation -> feseOperations,",
  "  basisFunctions -> {",
  "    {\"dxz\", \"dyz\"},",
  "    {\"px\", \"py\", \"pz\"}",
  "  },",
  "  InitialBondShells -> 3,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "feseHamiltonian = symham[1] + symham[3];"
}, "\n"];

spinSquareApplicationCode = StringRiffle[{
  "spinSquareOperations = msgop[gray[123]];",
  "init[",
  "  lattice -> {",
  "    {a, 0, 0},",
  "    {0, a, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 4},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> spinSquareOperations,",
  "  basisFunctions -> {{\"pxup\", \"pxdn\", \"pyup\", \"pydn\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False",
  "];",
  "spinSquareHamiltonian = Sum[symham[i], {i, 2}];"
}, "\n"];

cubicPApplicationCode[mode_String, name_String] :=
  cubicPOptionInitCode[mode] <>
    "\n" <> name <> " = symham[1];";

noncollinearApplicationCode = StringRiffle[{
  "identity3 = IdentityMatrix[3];",
  "c4z = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};",
  "c3xyz = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};",
  "noncollinearSSG = <|",
  "  \"C4zspin\" -> <|",
  "    \"space\" -> {identity3, {0, 0, 0}},",
  "    \"spin\" -> {c4z, 0}",
  "  |>,",
  "  \"C3xyzspin\" -> <|",
  "    \"space\" -> {identity3, {0, 0, 0}},",
  "    \"spin\" -> {c3xyz, 0}",
  "  |>",
  "|>;",
  "init[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> noncollinearSSG,",
  "  basisFunctions -> {{\"sup\", \"sdn\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> True",
  "];",
  "noncollinearHamiltonian = Sum[symham[i], {i, 2}];"
}, "\n"];

eptApplicationCode = StringRiffle[{
  "identity3 = IdentityMatrix[3];",
  "init[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {{{1/4, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> {",
  "    {\"E\", identity3, {0, 0, 0}, \"F\"},",
  "    {\"PT\", -identity3, {0, 0, 0}, \"T\"}",
  "  },",
  "  basisFunctions -> {{\"s\"}},",
  "  InitialBondShells -> 2,",
  "  GenerateSymmetryGroup -> False,",
  "  RepresentationMode -> \"Induced\"",
  "];",
  "eptHamiltonian = symham[2];"
}, "\n"];

c3z2tApplicationCode = StringRiffle[{
  "identity3 = IdentityMatrix[3];",
  "z3 = {",
  "  {-1/2, -Sqrt[3]/2, 0},",
  "  {Sqrt[3]/2, -1/2, 0},",
  "  {0, 0, 1}",
  "};",
  "z32 = MatrixPower[z3, 2];",
  "c3z2tOperations = <|",
  "  \"E\" -> <|",
  "    \"space\" -> {identity3, {0, 0, 0}},",
  "    \"spin\" -> {identity3, 0}",
  "  |>,",
  "  \"C3\" -> <|",
  "    \"space\" -> {identity3, {0, 0, 0}},",
  "    \"spin\" -> {z3, 0}",
  "  |>,",
  "  \"C3^2\" -> <|",
  "    \"space\" -> {identity3, {0, 0, 0}},",
  "    \"spin\" -> {z32, 0}",
  "  |>,",
  "  \"T\" -> <|",
  "    \"space\" -> {-identity3, {0, 0, 0}},",
  "    \"spin\" -> {identity3, 1}",
  "  |>,",
  "  \"C3T\" -> <|",
  "    \"space\" -> {-identity3, {0, 0, 0}},",
  "    \"spin\" -> {z3, 1}",
  "  |>,",
  "  \"C3^2T\" -> <|",
  "    \"space\" -> {-identity3, {0, 0, 0}},",
  "    \"spin\" -> {z32, 1}",
  "  |>",
  "|>;",
  "omega = Exp[2 Pi I/3];",
  "initfromrep[",
  "  lattice -> {{a, 0, 0}, {0, a, 0}, {0, 0, a}},",
  "  lattpar -> {a -> 1},",
  "  wyckoffposition -> {{{1/4, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> c3z2tOperations,",
  "  repinformation -> {",
  "    <|",
  "      \"ReferenceSiteIndex\" -> 1,",
  "      \"SiteSymmetryOperationIndices\" -> {1, 2, 3},",
  "      \"SiteSymmetryMatrices\" -> {",
  "        {{1}}, {{omega}}, {{omega^2}}",
  "      }",
  "    |>",
  "  },",
  "  orbitalLabels -> {{\"omega\"}},",
  "  RepresentationMode -> \"Induced\",",
  "  InitialBondShells -> 5",
  "];",
  "c3z2tHamiltonian = Sum[symham[i], {i, 5}];"
}, "\n"];

ClearAll[wyckoffExpressionCode, wyckoffIndentCode, wyckoffInitializationCode];
wyckoffExpressionCode[expression_] := StringReplace[
  ToString[expression, InputForm, PageWidth -> 58],
  {"MagneticTB`" -> "", "Global`" -> ""}
];
wyckoffIndentCode[code_String, count_Integer?NonNegative] :=
  StringReplace[code, "\n" -> "\n" <> StringRepeat[" ", count]];

wyckoffInitializationCode[
    case_Association,
    index_Integer,
    mode_String
  ] := Module[
  {prefix, latticeCode, parameterCode, keyCode, letter, groupCode},
  prefix = "wyckoff" <>
    If[mode === "DirectProduct", "Direct", "Induced"] <>
    ToString[index];
  latticeCode = wyckoffExpressionCode[
    wyckoffLattice[case["CrystalSystem"]]
  ];
  parameterCode = wyckoffExpressionCode[
    wyckoffLatticeParameters[case["CrystalSystem"]]
  ];
  keyCode = ToString[case["Key"], InputForm];
  letter = case["Letter"];
  groupCode = case["GroupCode"];
  StringRiffle[{
    "Needs[\"MagneticTB`\"];",
    prefix <> "PackageRoot = DirectoryName@DirectoryName@" <>
      "FindFile[\"MagneticTB`\"];",
    prefix <> "Database = Import[",
    "  FileNameJoin[{" <> prefix <> "PackageRoot, \"Data\", " <>
      "\"wyckoffMSG.mx\"}]",
    "];",
    prefix <> "LetterName[value_] := If[",
    "  StringQ[value], value, SymbolName[Unevaluated[value]]",
    "];",
    prefix <> "Entry = SelectFirst[",
    "  " <> prefix <> "Database[" <> keyCode <> "],",
    "  " <> prefix <> "LetterName[#[[3]]] === \"" <> letter <> "\" &",
    "];",
    prefix <> "Seed = " <> prefix <> "Entry[[4, 1]] /. {",
    "  x -> 1/7, y -> 2/11, z -> 3/13,",
    "  mx -> 1/5, my -> 2/7, mz -> 3/11",
    "};",
    prefix <> "Operations = msgop[" <> groupCode <> "];",
    "init[",
    "  lattice -> " <> wyckoffIndentCode[latticeCode, 2] <> ",",
    "  lattpar -> " <> wyckoffIndentCode[parameterCode, 2] <> ",",
    "  wyckoffposition -> {" <> prefix <> "Seed},",
    "  symminformation -> " <> prefix <> "Operations,",
    "  basisFunctions -> {{\"s\"}},",
    "  InitialBondShells -> 2,",
    "  GenerateSymmetryGroup -> False,",
    "  RepresentationMode -> \"" <> mode <> "\"",
    "];",
    prefix <> "Path = {",
    "  {{{0, 0, 0}, {1/2, 0, 0}}, {\"G\", \"X\"}},",
    "  {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {\"X\", \"M\"}},",
    "  {{{1/2, 1/2, 0}, {0, 0, 0}}, {\"M\", \"G\"}}",
    "};"
  }, "\n"]
];

(* The cubic 195.3-a phase sums are substantially wider than the other
   two-band Wyckoff examples.  Preserve the exact matrix, but scale only its
   displayed output so both DirectProduct and Induced results print in full. *)
wyckoffHamiltonianOutput[value_, case_Association] := If[
  case["ID"] === "195.3-a",
  veryWideMatrixOutput[MatrixForm[value]],
  MatrixForm[value]
];

wyckoffCaseCells[fixture_Association, index_Integer] := Module[
  {
    case, direct, induced, directName, inducedName,
    directCode, inducedCode
  },
  case = fixture["Case"];
  direct = fixture["DirectProduct"];
  induced = fixture["Induced"];
  directName = "wyckoffDirect" <> ToString[index];
  inducedName = "wyckoffInduced" <> ToString[index];
  directCode = wyckoffInitializationCode[
    case, index, "DirectProduct"
  ];
  inducedCode = wyckoffInitializationCode[
    case, index, "Induced"
  ];
  {
    Cell[
      case["CrystalSystem"] <> " \[LongDash] BNS " <> case["ID"],
      "Section"
    ],
    Cell[
      StringJoin[
        "For this Wyckoff position, the unitary subgroup leaves the representative site in its own orbit, ",
        "while an antiunitary operation generates the second site. The position shown below is read from the ",
        "magnetic Wyckoff database. We use it first in DirectProduct mode and then in Induced mode."
      ],
      "Text"
    ],
    Cell["DirectProduct representation", "Subsection"],
    Cell[BoxData[codeBoxes[directCode]], "Input"],
    Cell[BoxData[codeBoxes[directName <> "Seed"]], "Input"],
    Cell[BoxData[outputBoxes[direct["DatabaseSeed"]]], "Output"],
    Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
    Cell[BoxData[outputBoxes[
      compactOrbitalTable[direct["OrbitalTable"]]
    ]], "Output"],
    Cell[
      StringJoin[
        "The two rows of orbitalTable are the two sites generated from the database seed. They are also the ",
        "two basis states used by the Hamiltonian and by every displayed symmetry matrix."
      ],
      "Text"
    ],
    Cell[BoxData[codeBoxes[
      directName <> "Hamiltonian = Sum[symham[i], {i, 2}];"
    ]], "Input"],
    Cell[BoxData[codeBoxes[
      "MatrixForm[" <> directName <> "Hamiltonian]"
    ]], "Input"],
    Cell[BoxData[outputBoxes[
      wyckoffHamiltonianOutput[direct["Hamiltonian"], case]
    ]], "Output"],
    Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
    Cell[BoxData[outputBoxes[direct["Representation"]]], "Output"],
    Cell[BoxData[codeBoxes[
      "bandManipulate[" <> directName <> "Path, 25, " <>
        directName <> "Hamiltonian]"
    ]], "Input"],
    Cell[BoxData[outputBoxes[direct["BandManipulate"]]], "Output"],
    Cell["Induced representation", "Subsection"],
    Cell[BoxData[codeBoxes[inducedCode]], "Input"],
    Cell[BoxData[codeBoxes[inducedName <> "Seed"]], "Input"],
    Cell[BoxData[outputBoxes[induced["DatabaseSeed"]]], "Output"],
    Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
    Cell[BoxData[outputBoxes[
      compactOrbitalTable[induced["OrbitalTable"]]
    ]], "Output"],
    Cell[BoxData[codeBoxes[
      inducedName <> "Hamiltonian = Sum[symham[i], {i, 2}];"
    ]], "Input"],
    Cell[BoxData[codeBoxes[
      "MatrixForm[" <> inducedName <> "Hamiltonian]"
    ]], "Input"],
    Cell[BoxData[outputBoxes[
      wyckoffHamiltonianOutput[induced["Hamiltonian"], case]
    ]], "Output"],
    Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
    Cell[BoxData[outputBoxes[induced["Representation"]]], "Output"],
    Cell[
      StringJoin[
        "Induced starts from the site-symmetry representation of the first row and constructs the partner ",
        "state with a coset representative. Although its intermediate construction differs from ",
        "DirectProduct, the two Hamiltonians span the same allowed matrix space and their band degeneracies ",
        "agree."
      ],
      "Text"
    ],
    Cell[BoxData[codeBoxes[
      "bandManipulate[" <> inducedName <> "Path, 25, " <>
        inducedName <> "Hamiltonian]"
    ]], "Input"],
    Cell[BoxData[outputBoxes[induced["BandManipulate"]]], "Output"]
  }
];

wyckoffCrystalSystemsPage = Notebook[
  Join[
    {
      Cell["Wyckoff Database Models in Seven Crystal Systems", "Title"],
      Cell[
        StringJoin[
          "The following ten examples cover the seven crystal systems. Each representative position is read ",
          "directly from the magnetic Wyckoff database supplied with MagneticTB. The same orbit is then ",
          "described in two ways: DirectProduct uses the action on a complete local basis, whereas Induced ",
          "starts from the site-symmetry representation of one reference point. For both descriptions the ",
          "complete init input, Hamiltonian, symmetry matrices, and bands are shown."
        ],
        "Text"
      ],
      Cell[
        StringJoin[
          "Magnetic Wyckoff conventions follow the BNS magnetic-group tables described by S. V. Gallego et al.,",
          " \"Magnetic symmetry in the Bilbao Crystallographic Server,\" Journal of Applied Crystallography 45,",
          " 1236-1247 (2012), https://doi.org/10.1107/S0021889812042185. The bundled database is the ",
          "executable source used here; the citation documents the corresponding magnetic Wyckoff data ",
          "convention."
        ],
        "Text"
      ]
    },
    Flatten@MapIndexed[wyckoffCaseCells[#1, First[#2]] &, wyckoffDocumentationFixture],
    {
      Cell["How the two constructions are compared", "Section"],
      Cell[
        StringJoin[
          "The parameter names in the two descriptions need not be the same. What matters is that the two sets ",
          "of matrices describe the same allowed Hamiltonians and give the same bands. The comparison ",
          "therefore checks both directions between the two matrix spaces. It also applies every symmetry at ",
          "general momenta and at non-Gamma little-group points. When two momenta differ by a reciprocal ",
          "vector, the phase from the orbital positions is included. A commutator evaluated only at Gamma ",
          "would miss this information."
        ],
        "Text"
      ],
      historyCells[],
      categorizationCells[
        "Tech Note",
        "MagneticTB`",
        "MagneticTB/tutorial/WyckoffCrystalSystems"
      ],
      keywordCells[{
        "Wyckoff database", "crystal systems", "direct product",
        "induced representation", "antiunitary symmetry"
      }]
    }
  ],
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Wyckoff Database Models in Seven Crystal Systems",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

continuousDocumentationSource = FileNameJoin[{
  scriptDirectory,
  "Examples",
  "ContinuousSymmetryDocumentation.wl"
}];
If[!FileExistsQ[continuousDocumentationSource],
  Print[
    "Missing ContinuousSymmetry documentation source: ",
    continuousDocumentationSource
  ];
  Exit[18]
];
Get[continuousDocumentationSource];
continuousSymmetryPage = continuousSymmetryDocumentationPage;

validatedModelsPage = Notebook[{
  Cell["Representative MagneticTB Models", "Title"],
  Cell[
    StringJoin[
      "This page collects models that require more than the two-band graphene example. They include ",
      "unequal orbital dimensions on two atoms, several Wyckoff orbits, spinful magnetic symmetry, induced ",
      "representations, and an antiunitary operation that carries a state from one site to another. Each ",
      "section gives the complete input and then shows the resulting Hamiltonian, orbital order or hopping,",
      " symmetry matrices, and bands."
    ],
    "Text"
  ],
  Cell["Load MagneticTB", "Section"],
  Cell[BoxData[codeBoxes["Needs[\"MagneticTB`\"]"]], "Input"],

  Cell["SimpleCubicS: one s orbital on a simple-cubic site", "Section"],
  Cell[
    StringJoin[
      "Place one s orbital at the origin of a simple-cubic lattice. The gray Pm-3m group keeps this site ",
      "fixed, and the six bonds in shell 2 connect it to the nearest cells along plus or minus x, y, and z.",
      " The resulting one-band Hamiltonian contains one onsite energy and the usual sum of three cosine ",
      "dispersions."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    simpleCubicApplicationCode <> "\n" <> validationBandPathCode
  ]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[
      validationModelFixture["SimpleCubicS", "OrbitalTable"]
    ]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    "MatrixForm[simpleCubicHamiltonian]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[
    validationModelFixture["SimpleCubicS", "Hamiltonian"]
  ]]], "Output"],
  Cell[BoxData[codeBoxes[
    "bandManipulate[validationBandPath, 30, simpleCubicHamiltonian]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    validationModelFixture["SimpleCubicS", "BandUI"]
  ]], "Output"],

  Cell["CsCl: a 1 by 3 rectangular s-p hopping block", "Section"],
  Cell[
    StringJoin[
      "The corner site carries one s orbital and the body-centred site carries px, py, pz. The ",
      "off-diagonal row below is therefore the actual rectangular hopping block; it is not required to be ",
      "a square or Hermitian matrix by itself."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[csclApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes["MatrixForm[csclHamiltonian[[1 ;; 1, 2 ;; 4]]]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[validationModelFixture["CsCl", "CrossBlock"]]]], "Output"],
  Cell[BoxData[codeBoxes["showHoppingParameters[2]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CsCl", "HoppingSources"]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, csclHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CsCl", "BandUI"]]], "Output"],

  Cell["FeSe: two Wyckoff orbits with d and p orbitals", "Section"],
  Cell[
    StringJoin[
      "The first Wyckoff orbit carries dxz and dyz, while the second carries px, py, and pz. Because both ",
      "seeds generate equivalent partners, the complete basis has ten states. The orbital table fixes ",
      "their order; the displayed 3 by 3 corner is a readable part of the same ten-band Hamiltonian used ",
      "for the band plot."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[feseApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[validationModelFixture["FeSe", "OrbitalTable"]]
  ]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[feseHamiltonian[[1 ;; 3, 1 ;; 3]]]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[validationModelFixture["FeSe", "Block"]]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, feseHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["FeSe", "BandUI"]]], "Output"],

  Cell["MoS2: an eleven-band d-p orbital model", "Section"],
  Cell[
    StringJoin[
      "The first Wyckoff orbit carries the five metal d orbitals. The second seed expands to two chalcogen ",
      "sites, each with px, py, pz, giving eleven orbitals in total. The parameters below are illustrative ",
      "symmetry-allowed amplitudes; this example demonstrates the complete orbital hierarchy and does not ",
      "claim a fitted reproduction of a particular material data set."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[
    moS2ElevenBandApplicationCode <> "\n" <> validationBandPathCode
  ]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[
      validationModelFixture["MoS2ElevenBand", "OrbitalTable"]
    ]
  ]], "Output"],
  Cell[BoxData[codeBoxes[
    "MatrixForm[moS2ElevenHamiltonian[[1 ;; 5, 1 ;; 5]]]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[
    validationModelFixture["MoS2ElevenBand", "OnsiteDBlock"]
  ]]], "Output"],
  Cell[BoxData[codeBoxes[
    "bandManipulate[validationBandPath, 30, moS2ElevenHamiltonian]"
  ]], "Input"],
  Cell[BoxData[outputBoxes[
    validationModelFixture["MoS2ElevenBand", "BandUI"]
  ]], "Output"],

  Cell["Spinful square px/py model", "Section"],
  Cell[
    StringJoin[
      "The local ordering is px-up, px-down, py-up, py-down. The onsite Hamiltonian and full generator ",
      "matrices show how orbital and spin transformations are combined before the bond constraints are ",
      "solved."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[spinSquareApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes["MatrixForm[symham[1]]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[validationModelFixture["SpinfulSquare", "Onsite"]]]], "Output"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["SpinfulSquare", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, spinSquareHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["SpinfulSquare", "BandUI"]]], "Output"],

  Cell["Cubic p orbitals: DirectProduct and Induced representations", "Section"],
  Cell[
    StringJoin[
      "The cubic rotations of Pm-3m mix px, py, and pz. DirectProduct must therefore receive the complete ",
      "three-orbital basis on the representative site. Induced starts with only px at {1/2, 0, 0}; ",
      "the site-symmetry representation and coset transport generate the corresponding py and pz states ",
      "on the other two symmetry-equivalent sites. The two orbital tables make this reduction explicit."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[cubicPApplicationCode["DirectProduct", "cubicPDirectHamiltonian"]]], "Input"],
  Cell[BoxData[codeBoxes["orbitalTable[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CubicPDirect", "OrbitalTable"]]], "Output"],
  Cell[BoxData[codeBoxes[cubicPApplicationCode["Induced", "cubicPInducedHamiltonian"]]], "Input"],
  Cell[BoxData[codeBoxes["orbitalTable[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CubicPInduced", "OrbitalTable"]]], "Output"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CubicPInduced", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[cubicPInducedHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[
    MatrixForm[validationModelFixture["CubicPInduced", "Hamiltonian"]]
  ]], "Output"],

  Cell["Collinear discrete spin-space group", "Section"],
  Cell[
    StringJoin[
      "C4spin rotates spin by 90 degrees without rotating real space. tauHalfT combines a half translation ",
      "with an antiunitary spin action. With GenerateSymmetryGroup -> True, these two operations generate ",
      "the complete finite spin-space group used to constrain the Hamiltonian."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[collinearApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes["MatrixForm[collinearHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[compactMatrixOutput[
    MatrixForm[validationModelFixture["CollinearSSG", "Hamiltonian"]]
  ]]], "Output"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CollinearSSG", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, collinearHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["CollinearSSG", "BandUI"]]], "Output"],

  Cell["Noncollinear octahedral spin-space group", "Section"],
  Cell[
    StringJoin[
      "Two noncommuting pure-spin rotations generate the finite noncollinear group. The representation ",
      "output displays the spinor matrices actually used to constrain the onsite and hopping terms."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[noncollinearApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes["MatrixForm[noncollinearHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[validationModelFixture["NoncollinearSSG", "Hamiltonian"]]]], "Output"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["NoncollinearSSG", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, noncollinearHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["NoncollinearSSG", "BandUI"]]], "Output"],

  Cell["E plus PT: antiunitary induction between two sites", "Section"],
  Cell[
    StringJoin[
      "The identity E cannot move the site at x=1/4 to its partner at x=3/4. In this model the required ",
      "coset representative is PT, so the second site is reached by an antiunitary operation. ",
      "symminformation already records this fact; the representation table therefore shows only its matrix ",
      "U. Its action on the Hamiltonian is U Conjugate[H(k)] U-dagger."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[eptApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes[compactOrbitalTableCode]], "Input"],
  Cell[BoxData[outputBoxes[
    compactOrbitalTable[validationModelFixture["EPT", "OrbitalTable"]]
  ]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[eptHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[validationModelFixture["EPT", "Hamiltonian"]]]], "Output"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[All]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["EPT", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, eptHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["EPT", "BandUI"]]], "Output"],

  Cell["C3 times Z2^T over C3: a complex local character", "Section"],
  Cell[
    StringJoin[
      "The reference site has C3 symmetry and carries the one-dimensional character {1,omega,omega^2}. An ",
      "antiunitary coset maps it to the partner site. Antiunitarity complex-conjugates the local character ",
      "during this transport, so C3 is represented by omega on one site and omega^2 on the other."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[c3z2tApplicationCode <> "\n" <> validationBandPathCode]], "Input"],
  Cell[BoxData[codeBoxes["showSymmetryRepresentations[All]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["C3Z2T", "Representation"]]], "Output"],
  Cell[BoxData[codeBoxes["MatrixForm[c3z2tHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[MatrixForm[validationModelFixture["C3Z2T", "Hamiltonian"]]]], "Output"],
  Cell[BoxData[codeBoxes["bandManipulate[validationBandPath, 30, c3z2tHamiltonian]"]], "Input"],
  Cell[BoxData[outputBoxes[validationModelFixture["C3Z2T", "BandUI"]]], "Output"],

  Cell["Related applications", "Section"],
  Cell[
    StringJoin[
      "The General Examples tutorial develops graphene, the three-band MoS2 model, magnetic Weyl and ",
      "nodal-line models, and a C4 topological-insulator model. ",
      "The Wyckoff tutorial gives ten DirectProduct/Induced pairs drawn from all seven crystal systems."
    ],
    "Text"
  ],
  linkListingCell[
    "Text",
    {
      {"General MagneticTB Examples", "MagneticTB/tutorial/GeneralExamples"},
      {"Wyckoff Database Models", "MagneticTB/tutorial/WyckoffCrystalSystems"}
    }
  ],
  Cell["Original version", "Section"],
  Cell[
    StringJoin[
      "Run the original MagneticTB functions in a separate fresh kernel. The tutorial Using the Original ",
      "Version gives complete initold, symhamold, and bandManipulateold examples with their outputs."
    ],
    "Text"
  ],
  Cell[TextData[{documentationLink["Using MagneticTB 1.0", "MagneticTB/tutorial/LegacyBackend"]}], "Text"],
  historyCells[],
  categorizationCells["Tech Note", "MagneticTB`", "MagneticTB/tutorial/ValidatedModels"],
  keywordCells[{"validated models", "CsCl", "FeSe", "spin-space group", "induced representation", "antiunitary induction"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Representative MagneticTB Models",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

(* -------------------------------------------------------------------------
   Phase 4: write authoring sources and expected-text fixtures.

   This phase deliberately writes UTF-8 Notebook expressions, not final .nb
   files.  notebook_guard.py performs the canonical ASCII-Unicode conversion.
   ------------------------------------------------------------------------- *)

writeNotebookSource[relativePath_String, notebook_, expected_] := Module[
  {
    sourcePath,
    expectedPath,
    sourceText,
    expectedText,
    safeNotebook,
    sourcePageWidth
  },
  If[Head[notebook] =!= Notebook || !ListQ[expected],
    Print[InputForm@<|
      "DocumentationPage" -> relativePath,
      "NotebookHead" -> Head[notebook],
      "ExpectedHead" -> Head[expected],
      "ExpectedComponentHeads" -> If[
        Head[expected] === Join,
        Head /@ (List @@ expected),
        {}
      ]
    |>];
    Exit[17]
  ];
  sourcePath = FileNameJoin[{sourceRoot, relativePath <> ".source.wl"}];
  expectedPath = FileNameJoin[{expectedRoot, relativePath <> ".expected.txt"}];
  safeNotebook = stripSavedBandManipulateOutputs[notebook];
  (* Keep the interactive-band reference pages line-wrapped so future
     source edits remain readable without changing their Notebook structure. *)
  sourcePageWidth = If[
    MemberQ[
      FileNameJoin[{
        "ReferencePages",
        "Symbols",
        #
      }] & /@ {
        "bandManipulate",
        "bandManipulateEig",
        "bandManipulateold"
      },
      relativePath
    ],
    120,
    Infinity
  ];
  sourceText = ToString[
    safeNotebook,
    InputForm,
    PageWidth -> sourcePageWidth
  ];
  expectedText = StringRiffle[expected, "\n"] <> "\n";
  writeUTF8TextIfChanged[sourcePath, sourceText];
  writeUTF8TextIfChanged[expectedPath, expectedText];
  Print[sourcePath]
];

(* Render-QA occasionally exposes a printable-width defect in one of the two
   large model tutorials.  Keep their repair path genuinely incremental: both
   page definitions are evaluated above, but this mode writes only these two
   English authoring sources and then exits before the general page scan. *)
If[MatchQ[
    Environment["MAGNETICTB_DOC_RENDER_QA_TUTORIALS_ONLY"],
    1 | "1" | True
  ],
  writeNotebookSource[
    FileNameJoin[{"Tutorials", "GeneralExamples"}],
    generalExamplesPage,
    {
      "General MagneticTB Examples",
      "Three-band MoS2 model",
      "Magnetic C3 Weyl point and Chern diagnostic",
      "pointChernNumber",
      "Magnetic cubic nodal line",
      "C4 topological-insulator model",
      "Wyckoff, layer-group, and rod-group databases",
      "MagneticTB/tutorial/GeneralExamples"
    }
  ];
  If[
    !MemberQ[deferredDocumentationTutorialNames, "WyckoffCrystalSystems"],
    writeNotebookSource[
      FileNameJoin[{"Tutorials", "WyckoffCrystalSystems"}],
      wyckoffCrystalSystemsPage,
      Join[
        {
          "Wyckoff Database Models in Seven Crystal Systems",
          "DirectProduct representation",
          "Induced representation",
          "How the two constructions are compared",
          "MagneticTB/tutorial/WyckoffCrystalSystems"
        },
        Lookup[wyckoffCases, "ID"]
      ]
    ]
  ];
  ensureTutorialHelpHomeLinks[];
  Quit[0]
];

Scan[
  Function[spec,
    writeNotebookSource[
      FileNameJoin[{"ReferencePages", "Symbols", spec["Name"]}],
      functionPage[spec],
      DeleteDuplicates@Join[
        {
          spec["Name"],
          First[spec["Usage"]][[2]],
          "Examples",
          "Related Functions",
          "MagneticTB/ref/" <> spec["Name"]
        },
        If[
          Lookup[spec, "OptionTable", {}] === {},
          {},
          Join[
            {"Options"},
            StringTrim[#, "\""] & /@ Lookup[spec, "OptionTable"][[All, 1]]
          ]
        ],
        Lookup[spec, "Expected", {}],
        If[Lookup[spec, "Applications", {}] === {}, {}, {"Applications"}],
        If[Lookup[spec, "Properties", {}] === {}, {}, {"Properties & Relations"}],
        If[Lookup[spec, "Issues", {}] === {}, {}, {"Possible Issues"}]
      ]
    ]
  ],
  functionSpecifications
];

writeGuidePage[guidePage];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "GettingStarted"}],
  tutorialPage,
  {
    "Getting Started with MagneticTB",
    "Install, update, uninstall, and load",
    "MagneticTB help home",
    "PacletUninstall[\\\"MagneticTB\\\"]",
    "The five physical inputs",
    "Solve the Hamiltonian shell by shell",
    "Multiple orbitals and multiple Wyckoff positions",
    "MatrixForm[csclHamiltonian]",
    "MagneticTB/tutorial/GettingStarted"
  }
];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "LegacyBackend"}],
  legacyTutorialPage,
  {
    "Using MagneticTB 1.0",
    "Load MagneticTB 1.0",
    "Each version needs its own initialization",
    "two separate kernels",
    "MagneticTB/tutorial/LegacyBackend"
  }
];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "GeneralExamples"}],
  generalExamplesPage,
  {
    "General MagneticTB Examples",
    "Three-band MoS2 model",
    "Magnetic C3 Weyl point and Chern diagnostic",
    "pointChernNumber",
    "Magnetic cubic nodal line",
    "C4 topological-insulator model",
    "Wyckoff, layer-group, and rod-group databases",
    "MagneticTB/tutorial/GeneralExamples"
  }
];

If[
  !MemberQ[deferredDocumentationTutorialNames, "ValidatedModels"],
  writeNotebookSource[
    FileNameJoin[{"Tutorials", "ValidatedModels"}],
    validatedModelsPage,
    {
      "Representative MagneticTB Models",
      "CsCl: a 1 by 3 rectangular s-p hopping block",
      "FeSe: two Wyckoff orbits with d and p orbitals",
      "Cubic p orbitals: DirectProduct and Induced representations",
      "E plus PT: antiunitary induction between two sites",
      "MagneticTB/tutorial/ValidatedModels"
    }
  ]
];

If[
  !MemberQ[deferredDocumentationTutorialNames, "WyckoffCrystalSystems"],
  writeNotebookSource[
    FileNameJoin[{"Tutorials", "WyckoffCrystalSystems"}],
    wyckoffCrystalSystemsPage,
    Join[
      {
        "Wyckoff Database Models in Seven Crystal Systems",
        "DirectProduct representation",
        "Induced representation",
        "How the two constructions are compared",
        "MagneticTB/tutorial/WyckoffCrystalSystems"
      },
      Lookup[wyckoffCases, "ID"]
    ]
  ]
];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "ContinuousSymmetry"}],
  continuousSymmetryPage,
  continuousSymmetryExpectedFragments
];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "CrystalAndKPaths"}],
  testingCrystalTutorialPage,
  crystalKagomeExpectedFragments
];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "BandFitting"}],
  testingBandFittingTutorialPage,
  {
    "Fitting MagneticTB bands to VASP",
    "Read and inspect the three VASP bands",
    "Build the symmetry-generated model",
    "Fit the onsite and nearest-neighbor model",
    "Add second-neighbor hopping",
    "Fit the energy range of interest",
    "Fit near the K point",
    "Judge the plotted bands, not only the loss",
    "MagneticTB/tutorial/BandFitting"
  }
];

If[
  !MemberQ[deferredDocumentationTutorialNames, "RealSpaceAndTopology"],
  writeNotebookSource[
    FileNameJoin[{"Tutorials", "RealSpaceAndTopology"}],
    testingRealSpaceTutorialPage,
    {
      "Real-space Hamiltonians, surfaces, and topology",
      "Hamiltonians used below",
      "Wilson-loop plots",
      "Surface spectra",
      "C4 higher-order topological-insulator corner state",
      "Magnetic Weyl point",
      "Berry-curvature density",
      "MagneticTB/tutorial/RealSpaceAndTopology"
    }
  ]
];

writeNotebookSource[
  FileNameJoin[{"Tutorials", "Corepresentations"}],
  testingCorepTutorialPage,
  corepTutorialExpectedFragments
];

ensureTutorialHelpHomeLinks[];
