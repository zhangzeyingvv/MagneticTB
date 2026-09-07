(* ::Package:: *)

(* Corepresentations prints its tables and returns Null.  Capture those Print
   expressions in the authoring code and save them as real Print cells.
   The visible inputs remain ordinary public MagneticTB calls. *)

corepTutorialLoadCode = "Needs[\"MSGCorep`\"];\nNeeds[\"MagneticTB`\"];";
corepTutorialOperationsCode =
  "corepOperations = getMSGElemFromMSGCorep[{164, 87}];";

corepTutorialDoubleInitCode = StringRiffle[{
  "init[",
  "  lattice -> {",
  "    {0, -a, 0},",
  "    {Sqrt[3] a/2, a/2, 0},",
  "    {0, 0, c}",
  "  },",
  "  lattpar -> {a -> 1, c -> 2},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 0}}},",
  "  symminformation -> corepOperations,",
  "  basisFunctions -> {{\"sup\", \"sdn\"}}",
  "];",
  "corepHamiltonian = symham[1];",
  "corepHamiltonian // MatrixForm"
}, "\n"];

corepTutorialSingleInitCode = StringReplace[
  corepTutorialDoubleInitCode,
  "{{\"sup\", \"sdn\"}}" -> "{{\"s\"}}"
];

corepTutorialCalculationCode = StringRiffle[{
  "getTBBandCorep[",
  "  {164, 87}, corepHamiltonian, {e1 -> 0},",
  "  {{0, 0, 0}, {1/3, 1/3, 0}, {0, 0, 1/2}}",
  "]"
}, "\n"];

(* Preserve Print's evaluated arguments and Null return semantics.  Do not
   pass this through checkedDocumentationEvaluation, which suppresses Print.
   A failed calculation stops generation instead of writing a placeholder. *)
corepTutorialExampleCells[code_String] := Module[
  {printed = {}, value, messages},
  {value, messages} = Block[
    {$MessageList = {}, Print = Function[Null,
      AppendTo[printed, If[Length[{##}] === 1, First[{##}], SequenceForm[##]]];
      Null
    ]},
    {ToExpression[code, InputForm], $MessageList}
  ];
  If[value === $Failed || messages =!= {},
    Print["Corepresentations example failed: ", code, "\n", messages];
    Exit[10]
  ];
  Join[
    {Cell[BoxData[codeBoxes[code]], "Input"]},
    Cell[BoxData[outputBoxes[#]], "Print"] & /@ printed,
    If[value === Null, {}, {Cell[BoxData[outputBoxes[value]], "Output"]}]
  ]
];

(* Evaluate once in tutorial order.  The same operation list is retained for
   both initializations and both single-/double-valued calculations. *)
corepTutorialLoadCells = corepTutorialExampleCells[corepTutorialLoadCode];
corepTutorialOperationCells = corepTutorialExampleCells[corepTutorialOperationsCode];
corepTutorialDoubleInitCells = corepTutorialExampleCells[corepTutorialDoubleInitCode];
corepTutorialDoubleCells = corepTutorialExampleCells[corepTutorialCalculationCode];
corepTutorialSingleInitCells = corepTutorialExampleCells[corepTutorialSingleInitCode];
corepTutorialSingleCells = corepTutorialExampleCells[corepTutorialCalculationCode];

testingCorepTutorialPage = Notebook[Join[
  {
    Cell["Band corepresentations with MSGCorep", "Title"],
    gettingStartedHelpHomeCell[],
    Cell[
      StringJoin[
      "To identify the symmetry of a band at a high-symmetry k point, use its small ",
      "corepresentation. For magnetic space group 164.87, the examples below compare spinor ",
      "and scalar s bases at Gamma, K, and A."
    ],
      "Text"
    ],
    Cell[
      StringJoin[
      "This calculation uses the external MSGCorep package and its SpaceGroupIrep dependency. ",
      "Install both first, then load MSGCorep and MagneticTB in a fresh kernel as shown below. ",
      "MSGCorep must be loaded by the user:"
    ],
      "Text"
    ]
  },
  corepTutorialLoadCells,
  {
    Cell["Use one ordered MSGCorep operation list", "Section"],
    Cell[
      StringJoin[
      "Use getMSGElemFromMSGCorep to obtain the symmetry operations, and pass them to init ",
      "without changing their order. The printed primitive vectors give the lattice for this ",
      "example. The corepresentation calculation requires this operation list; do not replace ",
      "it with a list from msgop."
    ],
      "Text"
    ]
  },
  corepTutorialOperationCells,
  {
    Cell["Double-valued high-symmetry calculation", "Section"],
    Cell[
      StringJoin[
      "Place {sup,sdn} at the origin to calculate double-valued corepresentations. Here we ",
      "keep only the onsite term, with the same energy e1 for both states, to show how the ",
      "result depends on the chosen basis."
    ],
      "Text"
    ]
  },
  corepTutorialDoubleInitCells,
  {
    Cell[
      StringJoin[
      "Specify Gamma, K, and A in reciprocal fractional coordinates and call getTBBandCorep. ",
      "The function prints the k-point names and result tables. It returns Null, so no ",
      "additional Out cell is expected:"
    ],
      "Text"
    ]
  },
  corepTutorialDoubleCells,
  {
    Cell[
      StringJoin[
      "The columns list band indices, energy, degeneracy, and the small corepresentation. With ",
      "e1 = 0, the two bands have zero energy and are doubly degenerate at all three points. ",
      "Their labels are Gamma_4(2), K_4 K_5(2), and A_4(2), respectively; each label belongs ",
      "to its own k point."
    ],
      "Text"
    ],
    Cell["Single-valued basis", "Section"],
    Cell[
      StringJoin[
      "Now replace {sup,sdn} by a scalar s orbital. Use the same ordered symmetry operations, ",
      "run init again, and recalculate the Hamiltonian before finding its corepresentations:"
    ],
      "Text"
    ]
  },
  corepTutorialSingleInitCells,
  corepTutorialSingleCells,
  {
    Cell[
      StringJoin[
      "There is now one nondegenerate band, with labels Gamma_1(1), K_1(1), and A_1(1). Both ",
      "examples use zero onsite energy, but the spinor and scalar bases give double-valued and ",
      "single-valued corepresentations, respectively."
    ],
      "Text"
    ],
    historyCells[],
    categorizationCells["Tech Note", "MagneticTB`", "MagneticTB/tutorial/Corepresentations"],
    keywordCells[{"MSGCorep", "band corepresentation", "double group", "high-symmetry point"}]
  }
],
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle -> "Band corepresentations with MSGCorep",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"
  ]
];

corepTutorialExpectedFragments = {
  "Band corepresentations with MSGCorep",
  "Use one ordered MSGCorep operation list",
  "Double-valued high-symmetry calculation",
  "Single-valued basis",
  "getTBBandCorep[",
  "small corep",
  "MagneticTB/tutorial/Corepresentations",
  "paclet:MagneticTB/guide/MagneticTB"
};
