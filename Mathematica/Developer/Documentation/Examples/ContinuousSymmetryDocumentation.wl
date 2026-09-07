(* ::Package:: *)

(* Human-maintained source for the ContinuousSymmetry tutorial.

   The example uses the high-symmetry tetragonal type-IV magnetic space
   group 124.360 (P_c4/mcc).  Its 32 operations generate only two magnetic
   sites from the selected Wyckoff seed.  The complete p_x/p_y/p_z spinor
   basis is the smallest ordinary p shell in which onsite SOC can mix
   different orbital and spin states.

   The finite representation used by initfromrep is extracted from the
   immediately preceding, program-generated SOC model.  The spatial action
   is then written as a spin-space-group action and enlarged by the required
   pure-internal C_infinity rotation.  No hand-written representation matrix
   or external prototype is used.
*)

continuousAFMCommonCode = "Needs[\"MagneticTB`\"]";

continuousAFMSOCCode = StringRiffle[{
  "init[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 3/2},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 1}}},",
  "  symminformation -> msgop[typeIV[{124, 360}]],",
  "  basisFunctions -> {{",
  "    \"pxup\", \"pxdn\", \"pyup\", \"pydn\", \"pzup\", \"pzdn\"",
  "  }},",
  "  InitialBondShells -> 2",
  "];",
  "socHamiltonian = symham[1] + symham[2];",
  "socRepresentation = CurrentModelSession[][",
  "  \"FullRepresentation\", \"RepresentationMatrices\"",
  "];"
}, "\n"];

continuousAFMContinuousCode = StringRiffle[{
  "initfromrep[",
  "  lattice -> DiagonalMatrix[{a, a, c}],",
  "  lattpar -> {a -> 1, c -> 3/2},",
  "  wyckoffposition -> {{{0, 0, 0}, {0, 0, 1}}},",
  "  symminformation -> Append[",
  "    Association @ MapIndexed[",
  "      \"g\" <> ToString[#2[[1]]] -> <|",
  "        \"space\" -> #1[[{2, 3}]],",
  "        \"spin\" -> {",
  "          Det[#1[[2]]] #1[[2]], Boole[#1[[4]] == \"T\"]",
  "        }",
  "      |> &,",
  "      msgop[typeIV[{124, 360}]]",
  "    ],",
  "    \"C_infty\" -> <|",
  "      \"space\" -> {IdentityMatrix[3], {0, 0, 0}},",
  "      \"spin\" -> {{",
  "        {Cos[theta], -Sin[theta], 0},",
  "        {Sin[theta], Cos[theta], 0},",
  "        {0, 0, 1}",
  "      }, 0},",
  "      \"continuous\" -> theta",
  "    |>",
  "  ],",
  "  repinformation -> <|",
  "    \"Discrete\" -> {MapThread[",
  "      Function[{matrix, operation},",
  "        With[{site = First@FirstPosition[",
  "            First@CurrentModelSession[][",
  "              \"ModelSpecification\", \"SiteOrbits\"",
  "            ],",
  "            Mod[operation[[3]], 1]",
  "          ]},",
  "          matrix[[6 (site - 1) + Range[6], Range[6]]]",
  "        ]",
  "      ],",
  "      {socRepresentation, msgop[typeIV[{124, 360}]]}",
  "    ]},",
  "    \"Continuous\" -> <|",
  "      \"SiteMatrices\" -> {ConstantArray[",
  "        KroneckerProduct[",
  "          IdentityMatrix[3],",
  "          DiagonalMatrix[{Exp[-I theta/2], Exp[I theta/2]}]",
  "        ],",
  "        2",
  "      ]}",
  "    |>",
  "  |>,",
  "  orbitalLabels -> {{",
  "    \"pxup\", \"pxdn\", \"pyup\", \"pydn\", \"pzup\", \"pzdn\"",
  "  }},",
  "  InitialBondShells -> 2",
  "];",
  "noSOCHamiltonian = symham[1] + symham[2];"
}, "\n"];

continuousSymmetryDocumentationFixture =
  checkedDocumentationEvaluation[
    "ContinuousSymmetry high-symmetry tetragonal antiferromagnet",
    Module[
      {
        afmLatticeValue,
        afmSeedValue,
        afmBasisValue,
        operations,
        socOnsite,
        socHamiltonian,
        socSession,
        socFullMatrices,
        socStructure,
        identityValue,
        thetaValue,
        continuousSpinRotationValue,
        continuousLocalMatrixValue,
        finiteSpinSpaceOperationsValue,
        afmSitesValue,
        targetSiteIndicesValue,
        localDimensionValue,
        finiteLocalMatricesValue,
        continuousOperationsValue,
        noSOCOnsite,
        noSOCHamiltonian,
        bandPathValue,
        socRulesValue,
        noSOCRulesValue,
        socBandHamiltonian,
        noSOCBandHamiltonian,
        socBandPlot,
        noSOCBandPlot,
        socGammaMultiplicities,
        noSOCGammaMultiplicities
      },
      afmLatticeValue = DiagonalMatrix[{a, a, c}];
      afmSeedValue = {{{0, 0, 0}, {0, 0, 1}}};
      afmBasisValue = {
        "pxup", "pxdn", "pyup", "pydn", "pzup", "pzdn"
      };
      operations = msgop[typeIV[{124, 360}]];

      init[
        lattice -> afmLatticeValue,
        lattpar -> {a -> 1, c -> 3/2},
        wyckoffposition -> afmSeedValue,
        symminformation -> operations,
        basisFunctions -> {afmBasisValue},
        InitialBondShells -> 2
      ];
      socOnsite = symham[1];
      socHamiltonian = socOnsite + symham[2];
      socSession = CurrentModelSession[];
      socFullMatrices = socSession[
        "FullRepresentation", "RepresentationMatrices"
      ];
      socStructure = showCrystalStructure[];

      identityValue = IdentityMatrix[3];
      thetaValue = Global`theta;
      continuousSpinRotationValue = {
        {Cos[thetaValue], -Sin[thetaValue], 0},
        {Sin[thetaValue], Cos[thetaValue], 0},
        {0, 0, 1}
      };
      continuousLocalMatrixValue = KroneckerProduct[
        IdentityMatrix[3],
        DiagonalMatrix[{
          Exp[-I thetaValue/2], Exp[I thetaValue/2]
        }]
      ];
      finiteSpinSpaceOperationsValue = Association @ MapIndexed[
        Function[{record, index},
          "g" <> ToString[index[[1]]] -> <|
            "space" -> record[[{2, 3}]],
            "spin" -> {
              Det[record[[2]]] record[[2]],
              Boole[record[[4]] == "T"]
            }
          |>
        ],
        operations
      ];
      afmSitesValue = First@socSession[
        "ModelSpecification", "SiteOrbits"
      ];
      targetSiteIndicesValue = Map[
        Function[record,
          First@FirstPosition[
            afmSitesValue,
            Mod[record[[3]], 1]
          ]
        ],
        operations
      ];
      localDimensionValue = Length[afmBasisValue];
      finiteLocalMatricesValue = MapThread[
        Function[{matrix, targetSite},
          matrix[[
            localDimensionValue (targetSite - 1) +
              Range[localDimensionValue],
            Range[localDimensionValue]
          ]]
        ],
        {socFullMatrices, targetSiteIndicesValue}
      ];
      continuousOperationsValue = Append[
        finiteSpinSpaceOperationsValue,
        "C_infty" -> <|
          "space" -> {identityValue, {0, 0, 0}},
          "spin" -> {continuousSpinRotationValue, 0},
          "continuous" -> thetaValue
        |>
      ];

      initfromrep[
        lattice -> afmLatticeValue,
        lattpar -> {a -> 1, c -> 3/2},
        wyckoffposition -> afmSeedValue,
        symminformation -> continuousOperationsValue,
        repinformation -> <|
          "Discrete" -> {finiteLocalMatricesValue},
          "Continuous" -> <|
            "SiteMatrices" -> {{
              continuousLocalMatrixValue,
              continuousLocalMatrixValue
            }}
          |>
        |>,
        orbitalLabels -> {afmBasisValue},
        InitialBondShells -> 2
      ];
      noSOCOnsite = symham[1];
      noSOCHamiltonian = noSOCOnsite + symham[2];

      bandPathValue = {
        {{{0, 0, 0}, {0, 0, 1/2}}, {"G", "Z"}}
      };
      socRulesValue = Thread[
        {
          e1, e2, e3, e4, e5, e6, e7, e8,
          t1, t2, t3, t4, t5
        } -> {
          0.25, 0.25, 0, 0, -0.8, -0.4, -0.2, 0.2,
          -0.18, 0, 0, 0, 0
        }
      ];
      noSOCRulesValue = Thread[
        {e1, e2, e3, e4, e5, e6, t1, t2, t3} ->
          {-0.2, 0.2, 0.4, 0.8, 0, 0, -0.18, 0, 0}
      ];
      socBandHamiltonian = socHamiltonian /. socRulesValue;
      noSOCBandHamiltonian = noSOCHamiltonian /. noSOCRulesValue;
      socBandPlot = bandplot[
        bandPathValue,
        20,
        socBandHamiltonian,
        {},
        plotRange -> {-1.3, 1.3},
        FontSize -> 18,
        ImageSize -> 360
      ];
      noSOCBandPlot = bandplot[
        bandPathValue,
        20,
        noSOCBandHamiltonian,
        {},
        plotRange -> {-1.3, 1.3},
        FontSize -> 18,
        ImageSize -> 360
      ];
      socGammaMultiplicities = Last /@ Tally[
        Round[
          Sort@N@Eigenvalues[
            socBandHamiltonian /. {kx -> 0, ky -> 0, kz -> 0}
          ],
          10^-8
        ]
      ];
      noSOCGammaMultiplicities = Last /@ Tally[
        Round[
          Sort@N@Eigenvalues[
            noSOCBandHamiltonian /. {kx -> 0, ky -> 0, kz -> 0}
          ],
          10^-8
        ]
      ];

      <|
        "SOCStructure" -> socStructure,
        "SOCOnsiteBlock" -> socOnsite[[1 ;; 6, 1 ;; 6]],
        "NoSOCOnsiteBlock" -> noSOCOnsite[[1 ;; 6, 1 ;; 6]],
        "ParameterCounts" -> {
          Length[Variables[socHamiltonian]],
          Length[Variables[noSOCHamiltonian]]
        },
        "SOCBandPlot" -> socBandPlot,
        "NoSOCBandPlot" -> noSOCBandPlot,
        "GammaMultiplicities" -> <|
          "SOC" -> socGammaMultiplicities,
          "NoSOC" -> noSOCGammaMultiplicities
        |>
      |>
    ]
  ];

If[
  !AssociationQ[continuousSymmetryDocumentationFixture] ||
    MemberQ[continuousSymmetryDocumentationFixture, $Failed, Infinity],
  Print[
    "The ContinuousSymmetry high-symmetry antiferromagnetic example failed."
  ];
  Exit[18]
];

continuousSymmetryDocumentationPage = Notebook[{
  Cell[
    "High-symmetry tetragonal antiferromagnet: spin symmetry and SOC",
    "Title"
  ],
  Cell[
    StringJoin[
      "Consider the tetragonal antiferromagnet with magnetic space group 124.360, P_c4/mcc. ",
      "Its 32 operations generate two sites, at z=0 and z=1/2, with opposite moments. We ",
      "compare the Hamiltonians and bands with and without spin-orbit coupling (SOC). Keep all ",
      "px, py, and pz spin states so that onsite SOC can mix the different orbitals and spins."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[continuousAFMCommonCode]], "Input"],

  Cell["High-symmetry MSG 124.360 with SOC", "Section"],
  Cell[
    StringJoin[
      "With SOC, spin rotates together with the spatial operation. Pass the complete magnetic ",
      "space group to init, then calculate the onsite term and the nearest-neighbor hopping ",
      "between the two layers:"
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[continuousAFMSOCCode]], "Input"],
  Cell[BoxData[codeBoxes["showCrystalStructure[]"]], "Input"],
  Cell[
    BoxData[outputBoxes[
      continuousSymmetryDocumentationFixture["SOCStructure"]
    ]],
    "Output"
  ],
  Cell[
    StringJoin[
      "The two layers have opposite moments. The 6 by 6 matrix below is the onsite block of ",
      "the z=0 layer, with rows and columns ordered as {pxup,pxdn,pyup,pydn,pzup,pzdn}:"
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes["MatrixForm[symham[1][[1 ;; 6, 1 ;; 6]]]"]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[MatrixForm[
      continuousSymmetryDocumentationFixture["SOCOnsiteBlock"]
    ]]],
    "Output"
  ],

  Cell["Without SOC: add continuous C-infinity spin symmetry", "Section"],
  Cell[
    StringJoin[
      "Without SOC, spin can also rotate independently about the magnetic axis. Write the ",
      "finite operations as a spin-space group and add this continuous C_infinity symmetry. ",
      "The finite representation matrices below come from the preceding model; pass them to ",
      "initfromrep together with the continuous rotation."
    ],
    "Text"
  ],
  Cell[BoxData[codeBoxes[continuousAFMContinuousCode]], "Input"],
  Cell[
    BoxData[codeBoxes[
      "MatrixForm[symham[1][[1 ;; 6, 1 ;; 6]]]"
    ]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[MatrixForm[
      continuousSymmetryDocumentationFixture["NoSOCOnsiteBlock"]
    ]]],
    "Output"
  ],
  Cell[
    StringJoin[
      "Compare the two onsite matrices. C_infinity sets the px/py-pz spin-mixing terms to ",
      "zero. With onsite and nearest-neighbor hopping retained, the numbers of real parameters ",
      "are 13 with SOC and 9 without SOC:"
    ],
    "Text"
  ],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "{",
      "  Length[Variables[socHamiltonian]],",
      "  Length[Variables[noSOCHamiltonian]]",
      "}"
    }, "\n"]]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[
      continuousSymmetryDocumentationFixture["ParameterCounts"]
    ]],
    "Output"
  ],

  Cell["Band splitting along G-Z", "Section"],
  Cell[
    StringJoin[
      "The nearest neighbors lie along c, so plot the bands on G-Z. For clarity, keep one ",
      "nonzero hopping parameter in each model. In the SOC model, also turn on the two onsite ",
      "spin-orbital mixing terms shown in the matrix above."
    ],
    "Text"
  ],
  Cell["The bands with SOC are:", "Text"],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "bandplot[",
      "  {{{{0, 0, 0}, {0, 0, 1/2}}, {\"G\", \"Z\"}}},",
      "  20,",
      "  socHamiltonian /. Thread[",
      "    {e1, e2, e3, e4, e5, e6, e7, e8, t1, t2, t3, t4, t5} ->",
      "    {0.25, 0.25, 0, 0, -0.8, -0.4, -0.2, 0.2, -0.18, 0, 0, 0, 0}",
      "  ],",
      "  {},",
      "  plotRange -> {-1.3, 1.3}, FontSize -> 18, ImageSize -> 360",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[
      continuousSymmetryDocumentationFixture["SOCBandPlot"]
    ]],
    "Output"
  ],
  Cell["With C_infinity symmetry and no SOC, the bands are:", "Text"],
  Cell[
    BoxData[codeBoxes[StringRiffle[{
      "bandplot[",
      "  {{{{0, 0, 0}, {0, 0, 1/2}}, {\"G\", \"Z\"}}},",
      "  20,",
      "  noSOCHamiltonian /. Thread[",
      "    {e1, e2, e3, e4, e5, e6, t1, t2, t3} ->",
      "    {-0.2, 0.2, 0.4, 0.8, 0, 0, -0.18, 0, 0}",
      "  ],",
      "  {},",
      "  plotRange -> {-1.3, 1.3}, FontSize -> 18, ImageSize -> 360",
      "]"
    }, "\n"]]],
    "Input"
  ],
  Cell[
    BoxData[outputBoxes[
      continuousSymmetryDocumentationFixture["NoSOCBandPlot"]
    ]],
    "Output"
  ],
  Cell[
    StringJoin[
      "For these parameter values, the SOC model has six doubly degenerate levels at G. ",
      "Without SOC, the degeneracies at G are {4,2,4,2}. The fourfold degeneracies occur at ",
      "this high-symmetry point; they do not extend along the whole path."
    ],
    "Text"
  ],

  Cell["Scope of the continuous interface", "Section"],
  Cell[
    StringJoin[
      "The continuous operation used by initfromrep is one pure spin C_infinity rotation with ",
      "no spatial action. Supply the complete finite group and its matrices in matching order. ",
      "This rotation restricts the local states; it does not generate more atoms, complete the ",
      "finite group, or add neighbor hoppings."
    ],
    "Text"
  ],
  linkListingCell[
    "Text",
    {
      {"typeIV", "MagneticTB/ref/typeIV"},
      {"init", "MagneticTB/ref/init"},
      {"initfromrep", "MagneticTB/ref/initfromrep"},
      {"symham", "MagneticTB/ref/symham"},
      {"bandplot", "MagneticTB/ref/bandplot"},
      {"showCrystalStructure", "MagneticTB/ref/showCrystalStructure"},
      {"MagneticTB", "MagneticTB/guide/MagneticTB"}
    }
  ],
  historyCells[],
  categorizationCells[
    "Tech Note",
    "MagneticTB`",
    "MagneticTB/tutorial/ContinuousSymmetry"
  ],
  keywordCells[{
    "high-symmetry antiferromagnet",
    "type-IV magnetic space group",
    "spin-space group",
    "spin-orbit coupling",
    "continuous symmetry",
    "C infinity",
    "initfromrep"
  }]
},
  TaggingRules -> <|"Paclet" -> "MagneticTB"|>,
  WindowTitle ->
    "High-symmetry tetragonal antiferromagnet: spin symmetry and SOC",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"},
    "TechNotePageStylesExt.nb",
    CharacterEncoding -> "UTF-8"
  ]
];

continuousSymmetryExpectedFragments = {
  "High-symmetry tetragonal antiferromagnet: spin symmetry and SOC",
  "High-symmetry MSG 124.360 with SOC",
  "Without SOC: add continuous C-infinity spin symmetry",
  "Band splitting along G-Z",
  "bandplot[",
  "Scope of the continuous interface",
  "typeIV[{124, 360}]",
  "showCrystalStructure[]",
  "MagneticTB/tutorial/ContinuousSymmetry"
};
