(* Generate simplified-Chinese authoring notebooks from the canonical English
   authoring sources. Code boxes, symbol names, contexts, URIs, cell styles,
   and metadata labels are intentionally preserved. *)

ClearAll["Global`*"];

(* -------------------------------------------------------------------------
   Phase 0: discover the English authoring pages and select an incremental
   subset when a documentation-only environment flag is present.
   ------------------------------------------------------------------------- *)

developerRoot = DirectoryName[ExpandFileName[$InputFileName]];
projectRoot = DirectoryName[DirectoryName[developerRoot]];
englishSourceRoot = FileNameJoin[{developerRoot, "Sources"}];
sourceRoot = FileNameJoin[{developerRoot, "SourcesChineseSimplified"}];
expectedRoot = FileNameJoin[{developerRoot, "ExpectedChineseSimplified"}];
documentationRoot = FileNameJoin[{projectRoot, "Documentation", "ChineseSimplified"}];

(* The VASP import, interactive fitting, and comparison pages read the same
   real EIGENVAL sample as their English counterparts.  Keep it beside the
   translated notebooks so the executable input stays identical. *)
compareBandEigenvalSource = FileNameJoin[{projectRoot, "Examples", "EIGENVAL"}];
compareBandEigenvalTarget = FileNameJoin[{
  documentationRoot,
  "ReferencePages",
  "Symbols",
  "EIGENVAL"
}];

ClearAll[writeUTF8TextIfChanged];
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

englishSourceFiles = Sort@FileNames[
  "*.source.wl",
  englishSourceRoot,
  Infinity
];
relativePages = StringReplace[
  StringDrop[#, StringLength[englishSourceRoot] + 1],
  ".source.wl" ~~ EndOfString -> ""
] & /@ englishSourceFiles;

(* Translate one explicitly selected page during incremental documentation
   work.  Use the same slash-separated relative path as the source tree, for
   example ReferencePages/Symbols/hoppingData. *)
requestedPage = Environment["MAGNETICTB_DOC_PAGE"];
If[StringQ[requestedPage] && StringTrim[requestedPage] =!= "",
  requestedPage = FileNameJoin@FileNameSplit[requestedPage];
  relativePages = Select[relativePages, # === requestedPage &]
];

If[MatchQ[
    Environment["MAGNETICTB_DOC_GUIDE_ONLY"],
    1 | "1" | True
  ],
  relativePages = Select[
    relativePages,
    # === FileNameJoin[{"Guides", "MagneticTB"}] &
  ]
];

If[MatchQ[
    Environment["MAGNETICTB_DOC_PLOT_PROPERTY_ONLY"],
    1 | "1" | True
  ],
  relativePages = Select[
    relativePages,
    MemberQ[
      {
        FileNameJoin[{"ReferencePages", "Symbols", "buildRealSpaceHamiltonian"}],
        FileNameJoin[{"ReferencePages", "Symbols", "plotWilsonLoop"}],
        FileNameJoin[{"ReferencePages", "Symbols", "plotSurfaceSpectrum"}],
        FileNameJoin[{"ReferencePages", "Symbols", "plotBerryCurvature2D"}],
        FileNameJoin[{"ReferencePages", "Symbols", "plotBerryCurvature3D"}],
        FileNameJoin[{"ReferencePages", "Symbols", "plotRealSpaceWavefunction"}],
        FileNameJoin[{"Tutorials", "RealSpaceAndTopology"}]
      },
      #
    ] &
  ]
];

If[MatchQ[
    Environment["MAGNETICTB_DOC_RENDER_QA_TUTORIALS_ONLY"],
    1 | "1" | True
  ],
  relativePages = Select[
    relativePages,
    MemberQ[
      {
        FileNameJoin[{"Tutorials", "GeneralExamples"}],
        FileNameJoin[{"Tutorials", "WyckoffCrystalSystems"}]
      },
      #
    ] &
  ]
];

Scan[
  Function[path,
    Scan[
      If[!DirectoryQ[#], CreateDirectory[#, CreateIntermediateDirectories -> True]] &,
      {
        DirectoryName[FileNameJoin[{sourceRoot, path}]],
        DirectoryName[FileNameJoin[{expectedRoot, path}]],
        DirectoryName[FileNameJoin[{documentationRoot, path}]]
      }
    ]
  ],
  relativePages
];

If[
  (
    AnyTrue[
      {"compareBand", "vaspEig", "bandManipulateEig", "fittingTB"},
      MemberQ[
        relativePages,
        FileNameJoin[{"ReferencePages", "Symbols", #}]
      ] &
    ] ||
    MemberQ[relativePages, FileNameJoin[{"Tutorials", "BandFitting"}]]
  ) &&
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

commonTranslations = <|
  "Examples" -> "\:8303\:4f8b",
  "Basic Examples" -> "\:57fa\:672c\:8303\:4f8b",
  "More Examples" -> "\:66f4\:591a\:8303\:4f8b",
  "Examples Initialization" -> "\:8303\:4f8b\:521d\:59cb\:5316",
  "Related Functions" -> "\:76f8\:5173\:51fd\:6570",
  "Related Guides" -> "\:76f8\:5173\:6307\:5357",
  "Related Tutorials" -> "\:76f8\:5173\:6559\:7a0b",
  "Tutorials" -> "\:6559\:7a0b",
  "History" -> "\:5386\:53f2",
  "Categorization" -> "\:5206\:7c7b",
  "Keywords" -> "\:5173\:952e\:8bcd",
  "New in: " -> "\:5f15\:5165\:7248\:672c\:ff1a",
  " | Modified in: " -> " | \:4fee\:6539\:7248\:672c\:ff1a",
  " | Obsolete in: " -> " | \:5e9f\:5f03\:7248\:672c\:ff1a",
  "Getting Started with MagneticTB" -> "MagneticTB \:5165\:95e8",
  "Using the Legacy Backend" -> "\:4f7f\:7528\:65e7\:7248\:540e\:7aef",
  "Using MagneticTB 1.0" -> "使用 MagneticTB 1.0 版",
  "Validated MagneticTB Models" -> "\:7ecf\:9a8c\:8bc1\:7684 MagneticTB \:6a21\:578b",
  "MagneticTB help home" -> "MagneticTB \:5e2e\:52a9\:4e3b\:9875",
  "Fitting MagneticTB bands to VASP" ->
    "\:5c06 MagneticTB \:80fd\:5e26\:62df\:5408\:5230 VASP",
  StringJoin[
    "fit three VASP bands, extend the hopping range, and judge the result from complete ",
    "comparison plots."
  ] ->
    "\:62df\:5408\:4e09\:6761 VASP \:80fd\:5e26\:3001\:6269\:5c55 hopping \:8303\:56f4\:ff0c\:5e76\:7528\:5b8c\:6574\:5bf9\:6bd4\:56fe\:5224\:65ad\:7ed3\:679c\:3002",
  StringJoin[
    " \[LongDash] fit three VASP bands, extend the hopping range, and judge the result from complete ",
    "comparison plots."
  ] ->
    " \[LongDash] \:62df\:5408\:4e09\:6761 VASP \:80fd\:5e26\:3001\:6269\:5c55 hopping \:8303\:56f4\:ff0c\:5e76\:7528\:5b8c\:6574\:5bf9\:6bd4\:56fe\:5224\:65ad\:7ed3\:679c\:3002"
|>;

additionalTranslations = Get[
  FileNameJoin[{developerRoot, "ChineseTranslations.wl"}]
];
If[!AssociationQ[additionalTranslations],
  Print["ChineseTranslations.wl must evaluate to an Association."];
  Exit[3]
];

translations = <|
  FileNameJoin[{"ReferencePages", "Symbols", "msgop"}] -> <|
    "returns the ordered magnetic-space-group operations for the database identifier group." -> StringJoin[
      "\:8fd4\:56de\:6570\:636e\:5e93\:6807\:8bc6 group \:5bf9\:5e94\:7684\:6709\:5e8f\:78c1\:7a7a",
      "\:95f4\:7fa4\:64cd\:4f5c\:3002"
    ],
    "Typical identifiers are gray[n], typeI[n], typeIII[n], and typeIV[n]." -> "\:5e38\:7528\:6807\:8bc6\:5305\:62ec gray[n]\:3001typeI[n]\:3001typeIII[n] \:548c typeIV[n]\:3002",
    StringJoin[
      "The returned operation order is used consistently by init when constructing the finite group, ",
      "site action, and representation matrices."
    ] -> StringJoin[
      "init \:5728\:6784\:9020\:6709\:9650\:7fa4\:3001\:683c\:70b9\:4f5c\:7528\:548c\:8868\:793a",
      "\:77e9\:9635\:65f6\:59cb\:7ec8\:6cbf\:7528\:8fd9\:91cc\:8fd4\:56de\:7684\:64cd\:4f5c\:987a",
      "\:5e8f\:3002"
    ],
    StringJoin[
      "Obtain the operations of gray magnetic space group 191. Block suppresses the informational ",
      "lattice printout so the returned group order is easy to ",
      "see:"
    ] -> StringJoin[
      "\:53d6\:5f97\:7b2c 191 \:53f7\:7070\:78c1\:7a7a\:95f4\:7fa4\:7684\:64cd\:4f5c\:3002Block ",
      "\:6291\:5236\:6676\:683c\:4fe1\:606f\:8f93\:51fa\:ff0c\:4fbf\:4e8e\:76f4\:63a5\:67e5\:770b",
      "\:8fd4\:56de\:7684\:7fa4\:9636\:ff1a"
    ],
    "Inspect the spatial part of the identity operation record:" -> "\:67e5\:770b\:5355\:4f4d\:64cd\:4f5c\:8bb0\:5f55\:7684\:7a7a\:95f4\:90e8\:5206\:ff1a",
    "magnetic space group" -> "\:78c1\:7a7a\:95f4\:7fa4",
    "symmetry operations" -> "\:5bf9\:79f0\:64cd\:4f5c",
    "gray group" -> "\:7070\:7fa4"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "init"}] -> <|
    StringJoin[
      "initializes a model and prepares its symmetry representation, site data, bond shells, and ",
      "static directed-bond orbits."
    ] -> StringJoin[
      "\:521d\:59cb\:5316\:6a21\:578b\:ff0c\:5e76\:51c6\:5907\:5176\:5bf9\:79f0\:8868\:793a\:3001",
      "\:683c\:70b9\:6570\:636e\:3001\:8fd1\:90bb hopping \:5c42\:548c\:9759\:6001\:6709\:5411\:952e\:8f68",
      "\:9053\:3002"
    ],
    StringJoin[
      "The principal options are lattice, lattpar, wyckoffposition, symminformation, basisFunctions, ",
      "InitialBondShells, GenerateSymmetryGroup, RepresentationMode, and SiteLocalData."
    ] -> StringJoin[
      "\:4e3b\:8981\:9009\:9879\:5305\:62ec lattice\:3001lattpar\:3001wyckoffposition\:3001symminform",
      "ation\:3001basisFunctions\:3001InitialBondShells\:3001GenerateSymmetryGroup\:3001Representatio",
      "nMode \:548c SiteLocalData\:3002"
    ],
    "InitialBondShells has default value 10. symham never extends this range automatically." -> StringJoin[
      "InitialBondShells \:7684\:9ed8\:8ba4\:503c\:4e3a 10\:ff1bsymham ",
      "\:7edd\:4e0d\:4f1a\:81ea\:52a8\:6269\:5c55\:8fd9\:4e2a\:8303\:56f4\:3002"
    ],
    "RepresentationMode can be \"DirectProduct\" or \"Induced\"." -> "RepresentationMode \:53ef\:4ee5\:53d6 \"DirectProduct\" \:6216 \"Induced\"\:3002",
    StringJoin[
      "Every call clears the previous CurrentModelSession before compiling the new model. A failed ",
      "initialization leaves no usable old session."
    ] -> StringJoin[
      "\:6bcf\:6b21\:8c03\:7528\:90fd\:4f1a\:5148\:6e05\:7a7a\:65e7\:7684 CurrentModelSession\:ff0c",
      "\:518d\:7f16\:8bd1\:65b0\:6a21\:578b\:ff1b\:521d\:59cb\:5316\:5931\:8d25\:540e\:4e0d\:4f1a",
      "\:7559\:4e0b\:53ef\:7ee7\:7eed\:4f7f\:7528\:7684\:65e7\:4f1a\:8bdd\:3002"
    ],
    StringJoin[
      "Initialize the built-in default hexagonal model with one s orbital. Its default lattice ",
      "parameters are a -> 1 and c -> 3:"
    ] -> StringJoin[
      "\:521d\:59cb\:5316\:7a0b\:5e8f\:5185\:7f6e\:7684\:9ed8\:8ba4\:516d\:65b9\:6a21\:578b\:ff0c",
      "\:5176\:4e2d\:53ea\:6709\:4e00\:4e2a s \:8f68\:9053\:ff0c",
      "\:9ed8\:8ba4\:6676\:683c\:53c2\:6570\:4e3a a -> 1 \:548c c -> ",
      "3\:ff1a"
    ],
    StringJoin[
      "Confirm the representation mode, Hamiltonian dimension, and number of prepared ",
      "shells:"
    ] -> StringJoin[
      "\:786e\:8ba4\:8868\:793a\:6a21\:5f0f\:3001Hamiltonian \:7ef4\:6570\:548c\:5df2\:51c6\:5907",
      "\:58f3\:5c42\:6570\:ff1a"
    ],
    StringJoin[
      "A spinful induced example uses the complete magnetic group P6'/m'mm' (BNS 191.239). ",
      "One exact in-plane spinor on the reference kagome site generates the other two local ",
      "states.  The first three neighbour-hopping layers reproduce the three-band quadratic-",
      "nodal-point model of Feng et al., Phys. Rev. B 109, 045423 (2024).  In MagneticTB's ",
      "basis and directed-bond convention, the paper parameters map to the three rules shown ",
      "in bandplot:"
    ] -> StringJoin[
      "\:4e00\:4e2a\:542b\:81ea\:65cb\:7684\:8bf1\:5bfc\:8868\:793a\:8303\:4f8b\:4f7f\:7528\:5b8c\:6574\:78c1\:7fa4 ",
      "P6'/m'mm'\:ff08BNS 191.239\:ff09\:3002\:5728 Kagome \:53c2\:8003\:683c\:70b9\:4e0a\:53ea\:8f93\:5165\:4e00\:4e2a\:7cbe\:786e\:5e73\:9762\:5185\:65cb\:91cf\:6001\:ff0c",
      "\:5bf9\:79f0\:6027\:4fbf\:4f1a\:751f\:6210\:5176\:4ed6\:4e24\:4e2a\:5c40\:57df\:6001\:3002\:524d\:4e09\:4e2a\:8fd1\:90bb hopping \:5c42\:91cd\:73b0 Feng \:7b49\:4eba",
      "\:5728 Phys. Rev. B 109, 045423 (2024) \:4e2d\:7684\:4e09\:5e26\:4e8c\:6b21\:8282\:70b9\:6a21\:578b\:3002",
      "\:6309 MagneticTB \:7684\:57fa\:5e95\:987a\:5e8f\:548c\:6709\:5411 hopping \:7ea6\:5b9a\:ff0c\:8bba\:6587\:53c2\:6570\:5bf9\:5e94 bandplot \:4e2d\:663e\:793a\:7684\:4e09\:6761\:89c4\:5219\:ff1a"
    ],
    "initialize model" -> "\:521d\:59cb\:5316\:6a21\:578b",
    "Wyckoff position" -> "Wyckoff \:4f4d\:7f6e",
    "tight binding" -> "\:7d27\:675f\:7f1a",
    "bond shells" -> "\:8fd1\:90bb hopping \:5c42"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "initfromrep"}] -> <|
    StringJoin[
      "initializes a model from an already closed finite group and program-generated exact ",
      "representation matrices."
    ] -> StringJoin[
      "\:4ece\:5df2\:7ecf\:95ed\:5408\:7684\:6709\:9650\:7fa4\:548c\:7a0b\:5e8f\:751f\:6210\:7684",
      "\:7cbe\:786e\:8868\:793a\:77e9\:9635\:521d\:59cb\:5316\:6a21\:578b\:3002"
    ],
    StringJoin[
      "Unlike init, initfromrep does not accept basisFunctions, SiteLocalData, or ",
      "GenerateSymmetryGroup."
    ] -> StringJoin[
      "\:4e0e init \:4e0d\:540c\:ff0cinitfromrep \:4e0d\:63a5\:53d7 ",
      "basisFunctions\:3001SiteLocalData \:6216 GenerateSymmetryGroup\:3002"
    ],
    StringJoin[
      "The finite symmetry list must already be complete, closed, and ordered consistently with ",
      "repinformation."
    ] -> StringJoin[
      "\:6709\:9650\:5bf9\:79f0\:64cd\:4f5c\:5217\:8868\:5fc5\:987b\:5df2\:7ecf\:5b8c\:5907\:4e14",
      "\:95ed\:5408\:ff0c\:5e76\:4e0e repinformation \:4f7f\:7528\:4e00\:81f4\:7684\:987a\:5e8f\:3002"
    ],
    StringJoin[
      "A single pure-internal continuous C_infty generator is supported through the documented ",
      "Discrete and Continuous representation-information schema."
    ] -> StringJoin[
      "\:901a\:8fc7\:6587\:6863\:89c4\:5b9a\:7684 Discrete/Continuous ",
      "\:8868\:793a\:4fe1\:606f\:7ed3\:6784\:ff0c\:53ef\:4ee5\:652f\:6301\:4e00\:4e2a\:7eaf\:5185",
      "\:90e8\:7684\:8fde\:7eed C_infty \:751f\:6210\:5143\:3002"
    ],
    StringJoin[
      "Initialize a two-dimensional explicit representation of the identity group. The outer list ",
      "enumerates Wyckoff orbits and the inner list enumerates the complete ordered finite ",
      "group:"
    ] -> StringJoin[
      "\:7528\:5355\:4f4d\:7fa4\:7684\:4e8c\:7ef4\:663e\:5f0f\:8868\:793a\:521d\:59cb\:5316\:6a21",
      "\:578b\:3002\:5916\:5c42\:5217\:8868\:679a\:4e3e Wyckoff \:8f68\:9053\:ff0c",
      "\:5185\:5c42\:5217\:8868\:679a\:4e3e\:5b8c\:6574\:4e14\:6709\:5e8f\:7684\:6709\:9650\:7fa4",
      "\:ff1a"
    ],
    StringJoin[
      "Construct the onsite Hamiltonian to confirm that the supplied two-dimensional representation ",
      "fixes the Hamiltonian dimension:"
    ] -> StringJoin[
      "\:6784\:9020 onsite Hamiltonian\:ff0c\:786e\:8ba4\:8f93\:5165\:7684\:4e8c\:7ef4\:8868\:793a",
      "\:786e\:5b9a\:4e86 Hamiltonian \:7ef4\:6570\:ff1a"
    ],
    "explicit representation" -> "\:663e\:5f0f\:8868\:793a",
    "representation matrix" -> "\:8868\:793a\:77e9\:9635",
    "continuous symmetry" -> "\:8fde\:7eed\:5bf9\:79f0\:6027"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "symham"}] -> <|
    "constructs the symmetry-constrained Hamiltonian contribution from prepared bond shell n." -> StringJoin[
      "\:7531\:5df2\:7ecf\:51c6\:5907\:597d\:7684\:7b2c n \:4e2a\:8fd1\:90bb hopping \:5c42\:6784\:9020\:6ee1",
      "\:8db3\:5bf9\:79f0\:7ea6\:675f\:7684 Hamiltonian \:8d21\:732e\:3002"
    ],
    StringJoin[
      "keeps reverse directed hoppings independent and constructs a non-Hermitian Hamiltonian ",
      "contribution."
    ] -> StringJoin[
      "\:4fdd\:6301\:53cd\:5411\:6709\:5411 hopping \:76f8\:4e92\:72ec\:7acb\:ff0c",
      "\:5e76\:6784\:9020\:975e\:5384\:7c73 Hamiltonian \:8d21\:732e\:3002"
    ],
    StringJoin[
      "symham reads static bond-orbit data prepared by init; it does not call init, repeat the bond ",
      "search, or enlarge the prepared shell range."
    ] -> StringJoin[
      "symham \:53ea\:8bfb\:53d6 init \:51c6\:5907\:7684\:9759\:6001\:952e\:8f68\:9053\:6570\:636e",
      "\:ff1b\:5b83\:4e0d\:4f1a\:8c03\:7528 init\:3001\:91cd\:590d\:641c\:7d22\:952e\:6216\:6269",
      "\:5927\:5df2\:51c6\:5907\:7684\:58f3\:5c42\:8303\:56f4\:3002"
    ],
    StringJoin[
      "The default options are \"Hermitian\" -> True, \"KernelMethod\" -> \"Iterative\", ",
      "\"ValidationLevel\" -> \"Basic\", and \"CartesianCoordinates\" -> False."
    ] -> StringJoin[
      "\:9ed8\:8ba4\:9009\:9879\:4e3a \"Hermitian\" -> True\:3001\"KernelMethod\" -> ",
      "\"Iterative\"\:3001\"ValidationLevel\" -> \"Basic\" \:548c \"CartesianCoordinates\" -> ",
      "False\:3002"
    ],
    StringJoin[
      "Requesting an unprepared shell returns $Failed and instructs the user to rerun init with a ",
      "larger InitialBondShells value."
    ] -> StringJoin[
      "\:8bf7\:6c42\:5c1a\:672a\:51c6\:5907\:7684\:58f3\:5c42\:4f1a\:8fd4\:56de $Failed\:ff0c",
      "\:5e76\:63d0\:793a\:7528\:6237\:4ee5\:66f4\:5927\:7684 InitialBondShells ",
      "\:91cd\:65b0\:8fd0\:884c init\:3002"
    ],
    StringJoin[
      "Initialize the complete model used in this example. symham never performs this step ",
      "implicitly:"
    ] -> StringJoin[
      "\:5148\:521d\:59cb\:5316\:672c\:4f8b\:4f7f\:7528\:7684\:5b8c\:6574\:6a21\:578b\:3002symham ",
      "\:7edd\:4e0d\:4f1a\:9690\:5f0f\:6267\:884c\:8fd9\:4e00\:6b65\:ff1a"
    ],
    StringJoin[
      "Construct the onsite and first-neighbour contributions, then display their dimension and ",
      "parameter names:"
    ] -> StringJoin[
      "\:6784\:9020 onsite \:9879\:548c\:7b2c\:4e00\:8fd1\:90bb\:9879\:ff0c",
      "\:7136\:540e\:663e\:793a\:7ef4\:6570\:4e0e\:53c2\:6570\:540d\:ff1a"
    ],
    "Verify exact Hermiticity with real momentum variables:" -> "\:5c06\:52a8\:91cf\:53d8\:91cf\:89c6\:4e3a\:5b9e\:6570\:ff0c\:9a8c\:8bc1\:7cbe\:786e\:5384\:7c73\:6027\:ff1a",
    StringJoin[
      "Construct the first-neighbour contribution without Hermiticity ",
      "constraints:"
    ] -> StringJoin[
      "\:5728\:4e0d\:65bd\:52a0\:5384\:7c73\:7ea6\:675f\:65f6\:6784\:9020\:7b2c\:4e00\:8fd1\:90bb",
      "\:9879\:ff1a"
    ],
    "null space" -> "\:96f6\:7a7a\:95f4",
    "hopping" -> "\:8dc3\:8fc1",
    "Hermitian" -> "\:5384\:7c73",
    "symmetry constraints" -> "\:5bf9\:79f0\:7ea6\:675f"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "bandManipulate"}] -> <|
    "creates an interactive band plot of H along path with npoint subdivisions per path segment." -> StringJoin[
      "\:6cbf path \:7ed8\:5236 H \:7684\:4ea4\:4e92\:5f0f\:80fd\:5e26\:56fe\:ff0c",
      "\:6bcf\:6bb5\:8def\:5f84\:7ec6\:5206\:4e3a npoint \:4efd\:3002"
    ],
    "Symbols in H other than kx, ky, and kz are exposed as interactive parameters." -> StringJoin[
      "H \:4e2d\:9664 kx\:3001ky \:548c kz \:5916\:7684\:7b26\:53f7\:90fd\:4f1a\:4f5c\:4e3a\:4ea4",
      "\:4e92\:53c2\:6570\:663e\:793a\:3002"
    ],
    "Parameter sliders start at zero. The function does not add a hidden energy translation." -> StringJoin[
      "\:53c2\:6570\:6ed1\:5757\:4ece\:96f6\:5f00\:59cb\:ff1b\:8be5\:51fd\:6570\:4e0d\:4f1a\:52a0",
      "\:5165\:9690\:85cf\:7684\:80fd\:91cf\:5e73\:79fb\:3002"
    ],
    "The display convention for genuinely complex non-Hermitian bands is not defined." -> StringJoin[
      "\:771f\:6b63\:590d\:6570\:7684\:975e\:5384\:7c73\:80fd\:5e26\:76ee\:524d\:5c1a\:672a\:89c4",
      "\:5b9a\:663e\:793a\:7ea6\:5b9a\:3002"
    ],
    "Initialize a model and construct the Hamiltonian before plotting:" -> "\:7ed8\:56fe\:524d\:5148\:521d\:59cb\:5316\:6a21\:578b\:5e76\:6784\:9020 Hamiltonian\:ff1a",
    StringJoin[
      "Define a standard hexagonal path and open the interactive band ",
      "panel:"
    ] -> StringJoin[
      "\:5b9a\:4e49\:4e00\:6761\:6807\:51c6\:516d\:65b9\:6676\:683c\:8def\:5f84\:ff0c",
      "\:5e76\:6253\:5f00\:4ea4\:4e92\:5f0f\:80fd\:5e26\:9762\:677f\:ff1a"
    ],
    StringJoin[
      "With every Hamiltonian parameter set to one, the single Gamma-point eigenvalue is ",
      "seven:"
    ] -> StringJoin[
      "\:5c06\:6240\:6709 Hamiltonian \:53c2\:6570\:8bbe\:4e3a 1 \:65f6\:ff0c\:5355\:4e2a Gamma ",
      "\:70b9\:672c\:5f81\:503c\:4e3a 7\:ff1a"
    ],
    "band structure" -> "\:80fd\:5e26\:7ed3\:6784",
    "interactive plot" -> "\:4ea4\:4e92\:5f0f\:7ed8\:56fe",
    "tight-binding parameters" -> "\:7d27\:675f\:7f1a\:53c2\:6570"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "CurrentModelSession"}] -> <|
    "returns the current read-only in-memory model session created by init or initfromrep." -> StringJoin[
      "\:8fd4\:56de\:7531 init \:6216 initfromrep \:521b\:5efa\:7684\:5f53\:524d\:53ea\:8bfb\:5185",
      "\:5b58\:6a21\:578b\:4f1a\:8bdd\:3002"
    ],
    StringJoin[
      "The session contains the canonical model specification, compact full representation, compiled ",
      "bond shells, and solver caches."
    ] -> StringJoin[
      "\:4f1a\:8bdd\:5305\:542b\:89c4\:8303\:6a21\:578b\:63cf\:8ff0\:3001\:7d27\:51d1\:7684\:5b8c",
      "\:6574\:8868\:793a\:3001\:5df2\:7f16\:8bd1\:7684\:8fd1\:90bb hopping \:5c42\:548c\:6c42\:89e3\:5668",
      "\:7f13\:5b58\:3002"
    ],
    "Calling CurrentModelSession before successful initialization returns $Failed." -> StringJoin[
      "\:5728\:6210\:529f\:521d\:59cb\:5316\:524d\:8c03\:7528 CurrentModelSession \:4f1a\:8fd4\:56de ",
      "$Failed\:3002"
    ],
    StringJoin[
      "Display and export functions read the session but do not implicitly initialize or solve ",
      "missing shells."
    ] -> StringJoin[
      "\:663e\:793a\:548c\:5bfc\:51fa\:51fd\:6570\:53ea\:8bfb\:53d6\:4f1a\:8bdd\:ff0c",
      "\:4e0d\:4f1a\:9690\:5f0f\:521d\:59cb\:5316\:ff0c\:4e5f\:4e0d\:4f1a\:6c42\:89e3\:7f3a\:5931",
      "\:7684\:58f3\:5c42\:3002"
    ],
    "Initialize a complete model before reading its session:" -> "\:8bfb\:53d6\:4f1a\:8bdd\:524d\:5148\:521d\:59cb\:5316\:4e00\:4e2a\:5b8c\:6574\:6a21\:578b\:ff1a",
    StringJoin[
      "Summarize the representation source, representation mode, orbital dimension, and prepared ",
      "shell count:"
    ] -> StringJoin[
      "\:6c47\:603b\:8868\:793a\:6765\:6e90\:3001\:8868\:793a\:6a21\:5f0f\:3001\:8f68\:9053\:7ef4",
      "\:6570\:548c\:5df2\:51c6\:5907\:58f3\:5c42\:6570\:ff1a"
    ],
    StringJoin[
      "List the shell indices without compiling or solving anything ",
      "new:"
    ] -> StringJoin[
      "\:5728\:4e0d\:7f16\:8bd1\:6216\:6c42\:89e3\:4efb\:4f55\:65b0\:5185\:5bb9\:7684\:60c5\:51b5",
      "\:4e0b\:5217\:51fa\:58f3\:5c42\:7f16\:53f7\:ff1a"
    ],
    "model session" -> "\:6a21\:578b\:4f1a\:8bdd",
    "cache" -> "\:7f13\:5b58",
    "compiled bond shells" -> "\:5df2\:7f16\:8bd1\:8fd1\:90bb hopping \:5c42"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "orbitalTable"}] -> <|
    StringJoin[
      "returns a Dataset describing the fixed Hamiltonian orbital order prepared by init or ",
      "initfromrep."
    ] -> StringJoin[
      "\:8fd4\:56de\:4e00\:4e2a Dataset\:ff0c\:63cf\:8ff0 init \:6216 initfromrep ",
      "\:51c6\:5907\:7684\:56fa\:5b9a Hamiltonian \:8f68\:9053\:987a\:5e8f\:3002"
    ],
    StringJoin[
      "Rows include the orbital ID, Wyckoff orbit, equivalent atom, fractional and Cartesian ",
      "position, local index, basis state, and spin metadata."
    ] -> StringJoin[
      "\:6bcf\:4e00\:884c\:5305\:542b\:8f68\:9053 ID\:3001Wyckoff ",
      "\:8f68\:9053\:3001\:7b49\:4ef7\:539f\:5b50\:3001\:5206\:6570\:4e0e\:7b1b\:5361\:5c14\:5750",
      "\:6807\:3001\:5c40\:57df\:7f16\:53f7\:3001\:57fa\:6001\:548c\:81ea\:65cb\:5143\:6570\:636e",
      "\:3002"
    ],
    "The orbital order is the common order used by Hamiltonians and full representation matrices." -> StringJoin[
      "Hamiltonian \:548c\:5b8c\:6574\:8868\:793a\:77e9\:9635\:5171\:540c\:4f7f\:7528\:8fd9\:91cc",
      "\:7ed9\:51fa\:7684\:8f68\:9053\:987a\:5e8f\:3002"
    ],
    "orbitalTable does not call symham or mutate the model session." -> StringJoin[
      "orbitalTable \:4e0d\:4f1a\:8c03\:7528 symham\:ff0c\:4e5f\:4e0d\:4f1a\:4fee\:6539\:6a21\:578b",
      "\:4f1a\:8bdd\:3002"
    ],
    "Initialize a complete model before requesting its orbital table:" -> "\:8bf7\:6c42\:8f68\:9053\:8868\:524d\:5148\:521d\:59cb\:5316\:4e00\:4e2a\:5b8c\:6574\:6a21\:578b\:ff1a",
    StringJoin[
      "Display the orbital table. The result is a Dataset suitable for interactive ",
      "inspection:"
    ] -> StringJoin[
      "\:663e\:793a\:8f68\:9053\:8868\:3002\:7ed3\:679c\:662f\:4e00\:4e2a\:9002\:5408\:4ea4\:4e92",
      "\:67e5\:770b\:7684 Dataset\:ff1a"
    ],
    StringJoin[
      "Extract a compact, reproducible view of the orbital ID, fractional position, transported ",
      "basis state, and spin state:"
    ] -> StringJoin[
      "\:63d0\:53d6\:8f68\:9053 ID\:3001\:5206\:6570\:5750\:6807\:3001\:8f93\:8fd0\:540e\:57fa\:6001",
      "\:548c\:81ea\:65cb\:72b6\:6001\:7684\:7d27\:51d1\:53ef\:590d\:73b0\:89c6\:56fe\:ff1a"
    ],
    "orbital order" -> "\:8f68\:9053\:987a\:5e8f",
    "basis" -> "\:57fa\:5e95",
    "Wyckoff orbit" -> "Wyckoff \:8f68\:9053",
    "Hamiltonian indices" -> "Hamiltonian \:7f16\:53f7"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "showbonds"}] -> <|
    "displays the cached bond length and directed bonds of prepared bond shell n." -> StringJoin[
      "\:663e\:793a\:5df2\:51c6\:5907\:7684\:7b2c n \:4e2a\:8fd1\:90bb hopping \:5c42\:4e2d\:7f13\:5b58\:7684",
      "\:952e\:957f\:548c\:6709\:5411\:952e\:3002"
    ],
    "The displayed shell must already have been prepared by init." -> "\:8981\:663e\:793a\:7684\:58f3\:5c42\:5fc5\:987b\:5df2\:7ecf\:7531 init \:51c6\:5907\:3002",
    "showbonds never extends the bond search and never evaluates symham." -> StringJoin[
      "showbonds \:7edd\:4e0d\:4f1a\:6269\:5c55\:952e\:641c\:7d22\:ff0c",
      "\:4e5f\:4e0d\:4f1a\:8ba1\:7b97 symham\:3002"
    ],
    "The output identifies source and target sites and their translated target positions." -> StringJoin[
      "\:8f93\:51fa\:4f1a\:6807\:8bc6\:8d77\:70b9\:3001\:7ec8\:70b9\:683c\:70b9\:4ee5\:53ca\:5e73",
      "\:79fb\:540e\:7684\:7ec8\:70b9\:4f4d\:7f6e\:3002"
    ],
    StringJoin[
      "Initialize a complete model before displaying a prepared ",
      "shell:"
    ] -> StringJoin[
      "\:663e\:793a\:5df2\:51c6\:5907\:58f3\:5c42\:524d\:5148\:521d\:59cb\:5316\:4e00\:4e2a\:5b8c",
      "\:6574\:6a21\:578b\:ff1a"
    ],
    StringJoin[
      "Display the first nonzero-distance shell. The returned Grid is rendered directly in a ",
      "notebook:"
    ] -> StringJoin[
      "\:663e\:793a\:7b2c\:4e00\:4e2a\:975e\:96f6\:8ddd\:79bb\:58f3\:5c42\:3002\:8fd4\:56de\:7684 ",
      "Grid \:4f1a\:76f4\:63a5\:5728 notebook \:4e2d\:6e32\:67d3\:ff1a"
    ],
    StringJoin[
      "Read the already cached directed-bond count without extending the ",
      "search:"
    ] -> StringJoin[
      "\:5728\:4e0d\:6269\:5c55\:641c\:7d22\:7684\:60c5\:51b5\:4e0b\:8bfb\:53d6\:5df2\:7f13\:5b58",
      "\:7684\:6709\:5411\:952e\:6570\:ff1a"
    ],
    "bond shell" -> "\:8fd1\:90bb hopping \:5c42",
    "directed bond" -> "\:6709\:5411\:952e",
    "neighbour" -> "\:8fd1\:90bb",
    "hopping geometry" -> "\:8dc3\:8fc1\:51e0\:4f55"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "initold"}] -> <|
    "initializes a model with the fully isolated MagneticTB legacy backend." -> StringJoin[
      "\:4f7f\:7528\:5b8c\:5168\:9694\:79bb\:7684 MagneticTB \:65e7\:7248\:540e\:7aef\:521d\:59cb",
      "\:5316\:6a21\:578b\:3002"
    ],
    StringJoin[
      "The legacy backend is retained during the MagneticTB 2.x stabilization period so established ",
      "notebooks can still be evaluated."
    ] -> StringJoin[
      "\:5728 MagneticTB 2.x \:7a33\:5b9a\:671f\:95f4\:7ee7\:7eed\:4fdd\:7559\:65e7\:7248\:540e\:7aef",
      "\:ff0c\:4f7f\:5df2\:6709 notebook \:4ecd\:53ef\:8fd0\:884c\:3002"
    ],
    StringJoin[
      "MagneticTBOld` is a separate kernel context. It does not load, call, or share implementation ",
      "functions with MagneticTB`."
    ] -> StringJoin[
      "MagneticTBOld` \:662f\:72ec\:7acb\:7684\:5185\:6838 context\:ff1b",
      "\:5b83\:4e0d\:4f1a\:52a0\:8f7d\:3001\:8c03\:7528\:6216\:5171\:4eab MagneticTB` ",
      "\:7684\:5b9e\:73b0\:51fd\:6570\:3002"
    ],
    "For unambiguous comparisons, load MagneticTB` and MagneticTBOld` in separate Wolfram kernels." -> StringJoin[
      "\:4e3a\:907f\:514d\:6bd4\:8f83\:65f6\:4ea7\:751f\:6b67\:4e49\:ff0c",
      "\:5efa\:8bae\:5728\:4e0d\:540c\:7684 Wolfram \:5185\:6838\:4e2d\:5206\:522b\:52a0\:8f7d ",
      "MagneticTB` \:548c MagneticTBOld`\:3002"
    ],
    "Initialize a one-orbital legacy model:" -> "\:521d\:59cb\:5316\:4e00\:4e2a\:5355\:8f68\:9053\:65e7\:7248\:6a21\:578b\:ff1a",
    "legacy backend" -> "\:65e7\:7248\:540e\:7aef",
    "compatibility" -> "\:517c\:5bb9\:6027",
    "initialize model" -> "\:521d\:59cb\:5316\:6a21\:578b"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "symhamold"}] -> <|
    "constructs the legacy symmetry-constrained Hamiltonian contribution for bond shell n." -> StringJoin[
      "\:4f7f\:7528\:65e7\:7248\:540e\:7aef\:6784\:9020\:7b2c n ",
      "\:4e2a\:8fd1\:90bb hopping \:5c42\:7684\:5bf9\:79f0\:7ea6\:675f Hamiltonian \:8d21\:732e\:3002"
    ],
    StringJoin[
      "symhamold is the archived pre-2.0 implementation and is intentionally independent of the new ",
      "linear-algebra backend."
    ] -> StringJoin[
      "symhamold \:662f\:5f52\:6863\:7684 2.0 \:4e4b\:524d\:5b9e\:73b0\:ff0c",
      "\:5e76\:6709\:610f\:4e0e\:65b0\:7684\:7ebf\:6027\:4ee3\:6570\:540e\:7aef\:4fdd\:6301\:72ec",
      "\:7acb\:3002"
    ],
    "It is retained for compatibility and regression comparison while the 2.x series stabilizes." -> StringJoin[
      "\:5728 2.x \:7cfb\:5217\:7a33\:5b9a\:671f\:95f4\:4fdd\:7559\:5b83\:ff0c",
      "\:7528\:4e8e\:517c\:5bb9\:548c\:56de\:5f52\:6bd4\:8f83\:3002"
    ],
    StringJoin[
      "Initialize the complete legacy model used in this example. symhamold does not perform ",
      "initialization:"
    ] -> StringJoin[
      "\:5148\:521d\:59cb\:5316\:672c\:4f8b\:4f7f\:7528\:7684\:5b8c\:6574\:65e7\:7248\:6a21\:578b",
      "\:3002symhamold \:4e0d\:4f1a\:6267\:884c\:521d\:59cb\:5316\:ff1a"
    ],
    StringJoin[
      "Construct the first two legacy shell contributions and display their dimension and parameter ",
      "names:"
    ] -> StringJoin[
      "\:6784\:9020\:524d\:4e24\:4e2a\:65e7\:7248\:58f3\:5c42\:8d21\:732e\:ff0c",
      "\:5e76\:663e\:793a\:5176\:7ef4\:6570\:548c\:53c2\:6570\:540d\:ff1a"
    ],
    "legacy backend" -> "\:65e7\:7248\:540e\:7aef",
    "compatibility" -> "\:517c\:5bb9\:6027"
  |>,
  FileNameJoin[{"ReferencePages", "Symbols", "bandManipulateold"}] -> <|
    "creates the legacy interactive band plot for H along path." -> "\:6cbf path \:4e3a H \:521b\:5efa\:65e7\:7248\:4ea4\:4e92\:5f0f\:80fd\:5e26\:56fe\:3002",
    StringJoin[
      "This function belongs to the isolated MagneticTBOld` context and uses the legacy momentum ",
      "symbols kxold, kyold, and kzold."
    ] -> StringJoin[
      "\:8be5\:51fd\:6570\:5c5e\:4e8e\:9694\:79bb\:7684 MagneticTBOld` context\:ff0c",
      "\:5e76\:4f7f\:7528\:65e7\:7248\:52a8\:91cf\:7b26\:53f7 kxold\:3001kyold \:548c kzold\:3002"
    ],
    StringJoin[
      "Initialize a legacy model and construct its Hamiltonian before ",
      "plotting:"
    ] -> StringJoin[
      "\:7ed8\:56fe\:524d\:5148\:521d\:59cb\:5316\:65e7\:7248\:6a21\:578b\:5e76\:6784\:9020\:5176 ",
      "Hamiltonian\:ff1a"
    ],
    StringJoin[
      "Open the legacy interactive band panel. The returned expression has head ",
      "Manipulate:"
    ] -> StringJoin[
      "\:6253\:5f00\:65e7\:7248\:4ea4\:4e92\:5f0f\:80fd\:5e26\:9762\:677f\:3002",
      "\:8fd4\:56de\:8868\:8fbe\:5f0f\:7684\:5934\:90e8\:4e3a ",
      "Manipulate\:ff1a"
    ],
    "legacy backend" -> "\:65e7\:7248\:540e\:7aef",
    "band structure" -> "\:80fd\:5e26\:7ed3\:6784",
    "interactive plot" -> "\:4ea4\:4e92\:5f0f\:7ed8\:56fe"
  |>,
  FileNameJoin[{"Guides", "MagneticTB"}] -> <|
    "Tutorials: start here" -> "教程：从这里开始",
    "Model setup and inspection" -> "模型输入与查看",
    "MagneticTB 1.0 functions" -> "MagneticTB 1.0 版函数",
    "Install, update, and load the package, then build a graphene Hamiltonian and plot its bands." ->
      "安装、更新并加载程序包，从石墨烯开始构造 Hamiltonian 并画出能带。",
    "Follow complete model examples from symmetry operations to Hamiltonians and band structures." ->
      "从对称操作出发，跟着完整模型例子构造 Hamiltonian、计算能带。",
    "Draw a noncollinear magnetic Kagome structure, its Brillouin zone, and its bands." ->
      "以非共线磁性 Kagome 为例，画出晶体结构、磁矩、Brillouin 区和能带。",
    "Fit VASP bands with nearest- and second-neighbor hopping, or select an energy window or k-point region." ->
      "对比最近邻和第二近邻拟合，也可只拟合指定能量范围或某个 k 点附近的能带。",
    "Compare an antiferromagnet with SOC against the continuous spin symmetry required without SOC." ->
      "对比反铁磁体有 SOC 的模型与无 SOC 时需要的连续自旋对称性。",
    "Load MSGCorep explicitly and calculate single- and double-valued band corepresentations." ->
      "手动加载 MSGCorep，分别计算单值和双值能带共表示。",
    "Load version 1.0 in a separate kernel to run existing models." ->
      "在单独的内核中加载 1.0 版，运行已有模型。",
    StringJoin[
      "MagneticTB constructs symmetry-constrained tight-binding Hamiltonians for magnetic, ",
      "nonmagnetic, and discrete spin-space-group models."
    ] -> StringJoin[
      "MagneticTB \:4e3a\:78c1\:6027\:3001\:975e\:78c1\:6027\:548c\:79bb\:6563\:81ea\:65cb\:7a7a",
      "\:95f4\:7fa4\:6a21\:578b\:6784\:9020\:6ee1\:8db3\:5bf9\:79f0\:7ea6\:675f\:7684\:7d27\:675f",
      "\:7f1a Hamiltonian\:3002"
    ],
    "Symmetry and model initialization" -> "\:5bf9\:79f0\:6027\:4e0e\:6a21\:578b\:521d\:59cb\:5316",
    "Hamiltonian construction and visualization" -> "Hamiltonian \:6784\:9020\:4e0e\:53ef\:89c6\:5316",
    "Legacy compatibility backend" -> "\:65e7\:7248\:517c\:5bb9\:540e\:7aef",
    StringJoin[
      "MagneticTBOld` is shipped as a fully isolated compatibility backend during the MagneticTB 2.x ",
      "stabilization period. Load it in a separate kernel when comparing old and new results."
    ] -> StringJoin[
      "\:5728 MagneticTB 2.x \:7a33\:5b9a\:671f\:95f4\:ff0cMagneticTBOld` ",
      "\:4f5c\:4e3a\:5b8c\:5168\:9694\:79bb\:7684\:517c\:5bb9\:540e\:7aef\:4e00\:540c\:53d1\:5e03",
      "\:3002\:6bd4\:8f83\:65b0\:65e7\:7ed3\:679c\:65f6\:ff0c\:5efa\:8bae\:5728\:4e0d\:540c\:5185",
      "\:6838\:4e2d\:52a0\:8f7d\:3002"
    ],
    " \[LongDash] obtain magnetic-space-group operations from the bundled database." -> StringJoin[
      " \[LongDash] \:4ece\:968f\:5305\:6570\:636e\:5e93\:53d6\:5f97\:78c1\:7a7a\:95f4\:7fa4\:64cd",
      "\:4f5c\:3002"
    ],
    " \[LongDash] initialize a model from basis functions." -> " \[LongDash] \:4ece\:57fa\:51fd\:6570\:521d\:59cb\:5316\:6a21\:578b\:3002",
    " \[LongDash] initialize a model from exact representation matrices." -> " \[LongDash] \:4ece\:7cbe\:786e\:8868\:793a\:77e9\:9635\:521d\:59cb\:5316\:6a21\:578b\:3002",
    " \[LongDash] inspect the current compiled model session." -> " \[LongDash] \:67e5\:770b\:5f53\:524d\:5df2\:7f16\:8bd1\:6a21\:578b\:4f1a\:8bdd\:3002",
    " \[LongDash] inspect the fixed Hamiltonian orbital order." -> " \[LongDash] \:67e5\:770b\:56fa\:5b9a\:7684 Hamiltonian \:8f68\:9053\:987a\:5e8f\:3002",
    " \[LongDash] solve a prepared bond shell and construct its Hamiltonian contribution." -> StringJoin[
      " \[LongDash] \:6c42\:89e3\:5df2\:51c6\:5907\:7684\:8fd1\:90bb hopping \:5c42\:5e76\:6784\:9020\:5176 ",
      "Hamiltonian \:8d21\:732e\:3002"
    ],
    " \[LongDash] display prepared bond geometry." -> " \[LongDash] \:663e\:793a\:5df2\:51c6\:5907\:7684\:952e\:51e0\:4f55\:3002",
    " \[LongDash] explore the band structure interactively." -> " \[LongDash] \:4ea4\:4e92\:5f0f\:67e5\:770b\:80fd\:5e26\:7ed3\:6784\:3002",
    " \[LongDash] initialize a model with the archived backend." -> " \[LongDash] \:4f7f\:7528\:5f52\:6863\:540e\:7aef\:521d\:59cb\:5316\:6a21\:578b\:3002",
    " \[LongDash] construct a Hamiltonian with the archived backend." -> " \[LongDash] \:4f7f\:7528\:5f52\:6863\:540e\:7aef\:6784\:9020 Hamiltonian\:3002",
    " \[LongDash] inspect legacy bands interactively." -> " \[LongDash] \:4ea4\:4e92\:5f0f\:67e5\:770b\:65e7\:7248\:80fd\:5e26\:3002",
    "tight binding" -> "\:7d27\:675f\:7f1a",
    "magnetic symmetry" -> "\:78c1\:5bf9\:79f0\:6027",
    "band structure" -> "\:80fd\:5e26\:7ed3\:6784"
  |>,
  FileNameJoin[{"Tutorials", "GettingStarted"}] -> <|
    "View atomic positions and magnetization directions" -> "查看原子坐标和磁化方向",
    StringJoin[
      "After init, use atompos to view the generated atomic positions and magnetization ",
      "directions. The atoms are grouped by Wyckoff orbit. Each entry is {fractional position, ",
      "magnetization direction}, with both vectors expressed in the lattice-vector basis:"
    ] -> StringJoin[
      "初始化后，输入 atompos 就可以查看生成的原子坐标和磁化方向。结果按 Wyckoff ",
      "位置分组，每个原子写成 {分数坐标, 磁化方向}，两个向量都以晶格基矢为基："
    ],
    StringJoin[
      "The two entries are the carbon atoms at {1/3,2/3,0} and {2/3,1/3,0}. Their ",
      "magnetization directions are both {0,0,0}, since this model is nonmagnetic."
    ] -> StringJoin[
      "结果中的两个原子位于 {1/3,2/3,0} 和 {2/3,1/3,0}。这里不考虑磁性，所以磁化方向都",
      "是 {0,0,0}。"
    ],
    "Use atompos again to view the two atoms and their magnetization directions:" -> "同样可以用 atompos 查看这两个原子的坐标和磁化方向：",
    StringJoin[
      "The moments at {0,0,0} and {1/2,0,0} point along +z and -z. Use showCrystalStructure to ",
      "draw the atoms and their opposite moments:"
    ] -> StringJoin[
      "位于 {0,0,0} 和 {1/2,0,0} 的两个原子分别沿 +z 和 -z 方向磁化。使用 ",
      "showCrystalStructure 就可以画出原子及方向相反的磁矩："
    ],
    "Multiple orbitals and multiple Wyckoff positions" -> "考虑多个轨道和多个 Wyckoff 点",
    StringJoin[
      "For several inequivalent Wyckoff positions, list one representative of each in ",
      "wyckoffposition and place its orbitals in the corresponding entry of basisFunctions. In ",
      "the CsCl example below, the corner site carries s, and the body-center site carries px, ",
      "py, and pz. The full Hamiltonian is shown in the order {s,px,py,pz}, including the ",
      "onsite terms and hopping between the two sites."
    ] -> StringJoin[
      "考虑多个轨道和多个 Wyckoff 点时，在 wyckoffposition ",
      "中逐个列出代表位置，再按相同顺序填写 basisFunctions。下面的 CsCl ",
      "模型在顶角位置放置 s 轨道，在体心位置放置 px、py、pz 轨道。完整 Hamiltonian ",
      "的行列顺序为 {s,px,py,pz}，包含 onsite 项和两个位置之间的跃迁。"
    ],
    StringJoin[
      "This tutorial initializes a simple model, constructs a symmetry-constrained Hamiltonian, ",
      "inspects its basis and bonds, and opens an interactive band plot."
    ] -> StringJoin[
      "\:672c\:6559\:7a0b\:521d\:59cb\:5316\:4e00\:4e2a\:7b80\:5355\:6a21\:578b\:ff0c",
      "\:6784\:9020\:6ee1\:8db3\:5bf9\:79f0\:7ea6\:675f\:7684 Hamiltonian\:ff0c",
      "\:67e5\:770b\:5176\:57fa\:5e95\:4e0e\:952e\:ff0c\:5e76\:6253\:5f00\:4ea4\:4e92\:5f0f\:80fd",
      "\:5e26\:56fe\:3002"
    ],
    "Load MagneticTB" -> "\:52a0\:8f7d MagneticTB",
    "Initialize a model" -> "\:521d\:59cb\:5316\:6a21\:578b",
    StringJoin[
      "Calling init without options uses MagneticTB's built-in default hexagonal lattice, a -> 1, c ",
      "-> 3, one s orbital, and six prepared bond shells."
    ] -> StringJoin[
      "\:4e0d\:5e26\:9009\:9879\:8c03\:7528 init \:65f6\:ff0c\:7a0b\:5e8f\:4f7f\:7528 MagneticTB ",
      "\:5185\:7f6e\:7684\:9ed8\:8ba4\:516d\:65b9\:6676\:683c\:3001a -> 1\:3001c -> ",
      "3\:3001\:4e00\:4e2a s \:8f68\:9053\:548c\:516d\:4e2a\:5df2\:51c6\:5907\:8fd1\:90bb hopping \:5c42\:3002"
    ],
    "Inspect the compiled model" -> "\:67e5\:770b\:5df2\:7f16\:8bd1\:6a21\:578b",
    "Construct the Hamiltonian" -> "\:6784\:9020 Hamiltonian",
    "The onsite term is shell 1. Shell 2 is the first nonzero-distance neighbour shell." -> StringJoin[
      "onsite \:9879\:662f\:7b2c 1 \:58f3\:5c42\:ff1b\:7b2c 2 ",
      "\:58f3\:5c42\:662f\:7b2c\:4e00\:4e2a\:975e\:96f6\:8ddd\:79bb\:8fd1\:90bb\:58f3\:5c42\:3002"
    ],
    "Explore the bands" -> "\:67e5\:770b\:80fd\:5e26",
    "Important evaluation boundaries" -> "\:91cd\:8981\:7684\:8ba1\:7b97\:8fb9\:754c",
    StringJoin[
      "symham only solves shells prepared by init. If a larger shell is required, rerun init with a ",
      "larger InitialBondShells value. Display and export functions never trigger hidden ",
      "initialization or hidden shell solving."
    ] -> StringJoin[
      "symham \:53ea\:6c42\:89e3 init \:5df2\:51c6\:5907\:7684\:58f3\:5c42\:3002",
      "\:5982\:9700\:66f4\:5927\:7684\:58f3\:5c42\:ff0c\:5e94\:589e\:5927 InitialBondShells ",
      "\:5e76\:91cd\:65b0\:8fd0\:884c init\:3002\:663e\:793a\:548c\:5bfc\:51fa\:51fd\:6570\:7edd",
      "\:4e0d\:4f1a\:89e6\:53d1\:9690\:85cf\:521d\:59cb\:5316\:6216\:9690\:85cf\:7684\:58f3\:5c42",
      "\:6c42\:89e3\:3002"
    ],
    "Related functions" -> "\:76f8\:5173\:51fd\:6570",
    "getting started" -> "\:5165\:95e8",
    "tight binding" -> "\:7d27\:675f\:7f1a",
    "band structure" -> "\:80fd\:5e26\:7ed3\:6784"
  |>,
  FileNameJoin[{"Tutorials", "LegacyBackend"}] -> <|
    StringJoin[
      "MagneticTB 2.0 ships the archived implementation as a separate compatibility backend. It is ",
      "retained during the 2.x stabilization period and can be used to evaluate established ",
      "notebooks or compare results with the new implementation."
    ] -> StringJoin[
      "MagneticTB 2.0 \:5c06\:5f52\:6863\:5b9e\:73b0\:4f5c\:4e3a\:72ec\:7acb\:7684\:517c\:5bb9\:540e",
      "\:7aef\:4e00\:540c\:53d1\:5e03\:3002\:5b83\:4f1a\:5728 2.x ",
      "\:7a33\:5b9a\:671f\:95f4\:7ee7\:7eed\:4fdd\:7559\:ff0c\:53ef\:7528\:4e8e\:8fd0\:884c\:5df2",
      "\:6709 notebook \:6216\:4e0e\:65b0\:7248\:7ed3\:679c\:6bd4\:8f83\:3002"
    ],
    "Load the legacy context" -> "\:52a0\:8f7d\:65e7\:7248 context",
    StringJoin[
      "Use a fresh Wolfram kernel for the legacy backend. This keeps the old and new global states, ",
      "symbols, and implementations completely separate."
    ] -> StringJoin[
      "\:8bf7\:5728\:5168\:65b0\:7684 Wolfram \:5185\:6838\:4e2d\:4f7f\:7528\:65e7\:7248\:540e\:7aef",
      "\:ff0c\:4ece\:800c\:4f7f\:65b0\:65e7\:5168\:5c40\:72b6\:6001\:3001\:7b26\:53f7\:548c\:5b9e",
      "\:73b0\:5b8c\:5168\:5206\:79bb\:3002"
    ],
    "Initialize a legacy model" -> "\:521d\:59cb\:5316\:65e7\:7248\:6a21\:578b",
    "Construct and display the legacy Hamiltonian" -> "\:6784\:9020\:5e76\:663e\:793a\:65e7\:7248 Hamiltonian",
    "Isolation policy" -> "\:9694\:79bb\:539f\:5219",
    StringJoin[
      "MagneticTBOld` never loads or calls MagneticTB`, and MagneticTB` never loads or calls ",
      "MagneticTBOld`. The legacy backend will not be removed merely because version 2.0 exists; ",
      "removal is a later explicit decision after several stable releases and regression comparisons."
    ] -> StringJoin[
      "MagneticTBOld` \:6c38\:8fdc\:4e0d\:4f1a\:52a0\:8f7d\:6216\:8c03\:7528 MagneticTB`\:ff0c",
      "MagneticTB` \:4e5f\:6c38\:8fdc\:4e0d\:4f1a\:52a0\:8f7d\:6216\:8c03\:7528 MagneticTBOld`\:3002",
      "\:65e7\:7248\:540e\:7aef\:4e0d\:4f1a\:4ec5\:4ec5\:56e0\:4e3a 2.0 ",
      "\:5df2\:53d1\:5e03\:5c31\:88ab\:5220\:9664\:ff1b\:5728\:7ecf\:8fc7\:82e5\:5e72\:7a33\:5b9a",
      "\:7248\:672c\:548c\:56de\:5f52\:6bd4\:8f83\:540e\:ff0c\:662f\:5426\:5220\:9664\:4ecd\:9700",
      "\:53e6\:884c\:660e\:786e\:51b3\:5b9a\:3002"
    ],
    "Related functions" -> "\:76f8\:5173\:51fd\:6570",
    "legacy backend" -> "\:65e7\:7248\:540e\:7aef",
    "compatibility" -> "\:517c\:5bb9\:6027",
    "migration" -> "\:8fc1\:79fb",
    "regression comparison" -> "\:56de\:5f52\:6bd4\:8f83"
  |>,
  FileNameJoin[{"Tutorials", "ValidatedModels"}] -> <|
    StringJoin[
      "Every model on this page is a complete, independently runnable input used by the MagneticTB ",
      "regression suite. The displayed outputs were generated in fresh kernels. Each summary is ",
      "{Hamiltonian dimension, parameter counts for shells 1 and 2}."
    ] -> StringJoin[
      "\:672c\:9875\:6bcf\:4e2a\:6a21\:578b\:90fd\:662f MagneticTB ",
      "\:56de\:5f52\:6d4b\:8bd5\:4f7f\:7528\:7684\:5b8c\:6574\:3001\:53ef\:72ec\:7acb\:8fd0\:884c",
      "\:8f93\:5165\:3002\:663e\:793a\:7684\:8f93\:51fa\:5747\:7531\:5168\:65b0\:5185\:6838\:751f",
      "\:6210\:3002\:6bcf\:4e2a\:6458\:8981\:7684\:683c\:5f0f\:4e3a {Hamiltonian \:7ef4\:6570, ",
      "\:7b2c 1 \:548c\:7b2c 2 \:58f3\:5c42\:7684\:53c2\:6570\:6570}\:3002"
    ],
    "Load MagneticTB and define a compact parameter counter" -> "\:52a0\:8f7d MagneticTB \:5e76\:5b9a\:4e49\:7d27\:51d1\:7684\:53c2\:6570\:8ba1\:6570\:51fd\:6570",
    "Graphene: one pz orbital on each of two equivalent sites" -> "\:77f3\:58a8\:70ef\:ff1a\:4e24\:4e2a\:7b49\:4ef7\:683c\:70b9\:5404\:53d6\:4e00\:4e2a pz \:8f68\:9053",
    "CsCl: rectangular s-to-p hopping between different local dimensions" -> "CsCl\:ff1a\:4e0d\:540c\:5c40\:57df\:7ef4\:6570\:4e4b\:95f4\:7684\:77e9\:5f62 s-p hopping",
    "MoS2: eleven-band d-p model" -> "MoS2\:ff1a11 \:80fd\:5e26 d-p \:6a21\:578b",
    "FeSe: ten-band, two-Wyckoff model" -> "FeSe\:ff1a10 \:80fd\:5e26\:3001\:4e24\:4e2a Wyckoff \:8f68\:9053\:6a21\:578b",
    "Discrete spin-space group: pure internal C4 spin rotation" -> "\:79bb\:6563\:81ea\:65cb\:7a7a\:95f4\:7fa4\:ff1a\:7eaf\:5185\:90e8 C4 \:81ea\:65cb\:65cb\:8f6c",
    "Discrete spin-space group: half translation times time reversal" -> "\:79bb\:6563\:81ea\:65cb\:7a7a\:95f4\:7fa4\:ff1a\:534a\:5e73\:79fb\:4e58\:4ee5\:65f6\:95f4\:53cd\:6f14",
    "Direct-product and induced representations in all seven crystal systems" -> "\:4e03\:5927\:6676\:7cfb\:4e2d\:7684\:76f4\:79ef\:8868\:793a\:4e0e\:8bf1\:5bfc\:8868\:793a",
    StringJoin[
      "The complete seven-crystal-system input is provided in the following example notebook. It ",
      "evaluates both RepresentationMode values for each model and returns identical shell parameter ",
      "counts {1,2} in every direct/induced pair."
    ] -> StringJoin[
      "\:4e0b\:9762\:7684\:4f8b\:5b50 notebook \:7ed9\:51fa\:5b8c\:6574\:7684\:4e03\:6676\:7cfb\:8f93",
      "\:5165\:3002\:5b83\:5bf9\:6bcf\:4e2a\:6a21\:578b\:5206\:522b\:8ba1\:7b97\:4e24\:79cd ",
      "RepresentationMode\:ff0c\:5e76\:4e14\:6bcf\:4e00\:5bf9\:76f4\:79ef/\:8bf1\:5bfc\:7ed3\:679c",
      "\:90fd\:8fd4\:56de\:76f8\:540c\:7684\:58f3\:5c42\:53c2\:6570\:6570 {1,2}\:3002"
    ],
    "Legacy backend" -> "\:65e7\:7248\:540e\:7aef",
    StringJoin[
      "The old implementation must be evaluated in a separate fresh kernel. See Using the Legacy ",
      "Backend for a complete initold, symhamold, and bandManipulateold example with outputs."
    ] -> StringJoin[
      "\:65e7\:7248\:5b9e\:73b0\:5fc5\:987b\:5728\:72ec\:7acb\:7684\:5168\:65b0\:5185\:6838\:4e2d",
      "\:8ba1\:7b97\:3002\:5b8c\:6574\:4e14\:5e26\:8f93\:51fa\:7684 initold\:3001symhamold \:548c ",
      "bandManipulateold \:4f8b\:5b50\:8bf7\:53c2\:9605\:300a\:4f7f\:7528\:65e7\:7248\:540e\:7aef",
      "\:300b\:3002"
    ],
    "validated models" -> "\:7ecf\:9a8c\:8bc1\:6a21\:578b",
    "spin-space group" -> "\:81ea\:65cb\:7a7a\:95f4\:7fa4",
    "induced representation" -> "\:8bf1\:5bfc\:8868\:793a"
  |>
|>;

translateNotebook[notebook_Notebook, rules_List] := Module[{result},
  result = notebook /. Cell[text_String, style_, options___] :>
    Cell[Replace[text, rules], style, options];
  result = result /. ButtonBox[label_String, options___] :>
    ButtonBox[Replace[label, rules], options];
  result /. TextData[items_List] :>
    TextData[Replace[items, text_String :> Replace[text, rules], {1}]]
];

(* -------------------------------------------------------------------------
   Phase 1: page-level Chinese summaries for testing-line documentation.

   The testing-line pages are generated from one shared English specification.
   Their executable cells and outputs remain identical, while this page-level
   pass supplies concise, natural Chinese descriptions for the complete new
   API surface.  Long UTF-8 prose is split with StringJoin for source
   readability; notebook_guard.py escapes it only in the final .nb files.
   ------------------------------------------------------------------------- *)
testingChinesePageSummaries = <|
  FileNameJoin[{"ReferencePages", "Symbols", "bandplot"}] ->
    StringJoin[
      "bandplot \:6cbf\:5012\:7a7a\:95f4\:5206\:6570\:5750\:6807\:8def\:5f84\:8ba1\:7b97\:5e76\:7ed8\:5236",
      "\:6570\:503c\:80fd\:5e26\:3002\:6bcf\:6bb5\:5199\:6210 {{k1,k2},{label1,label2}}\:ff0c",
      "\:7a0b\:5e8f\:5728\:5185\:90e8\:4f7f\:7528 2 Pi k\:3002plotRange ",
      "\:53ea\:88c1\:526a\:7eb5\:8f74\:ff1byTicks \:8bbe\:7f6e\:7eb5\:8f74\:523b\:5ea6\:ff1b",
      "FontSize\:3001FontFamily \:548c ImageSize \:63a7\:5236\:6392\:7248\:3002Automatic ",
      "\:8def\:5f84\:9700\:8981\:5148\:5b8c\:6210 init\:3002\:8def\:5f84\:3001\:91c7\:6837\:6570\:3001",
      "\:53c2\:6570\:89c4\:5219\:6216\:663e\:793a\:9009\:9879\:975e\:6cd5\:65f6\:8fd4\:56de $Failed\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "showband"}] ->
    StringJoin[
      "showband[npoint,H,rules] \:4f7f\:7528\:5f53\:524d\:6676\:683c\:7684 standardKPath ",
      "\:7ed8\:5236\:9759\:6001\:80fd\:5e26\:ff0c\:5e76\:5171\:4eab bandplot \:7684 ",
      "plotRange\:3001yTicks\:3001FontSize\:3001FontFamily \:548c ImageSize\:3002\:5f53\:524d testing ",
      "\:5b9e\:73b0\:4e2d\:ff0cUsage.wl \:6240\:5199\:7684\:56db\:53c2\:6570\:663e\:5f0f\:8def\:5f84\:5f62",
      "\:5f0f\:4f1a\:88ab\:81ea\:52a8\:8def\:5f84\:91cd\:8f7d\:8bef\:5339\:914d\:ff0c",
      "\:56e0\:6b64\:672c\:9875\:53ea\:8bb0\:5f55\:5df2\:7ecf\:9a8c\:8bc1\:7684\:4e09\:53c2\:6570\:5f62",
      "\:5f0f\:ff1b\:663e\:5f0f\:8def\:5f84\:8bf7\:8c03\:7528 bandplot\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "standardKPath"}] ->
    StringJoin[
      "standardKPath \:6839\:636e\:5f53\:524d\:6570\:503c\:6676\:683c\:7684\:5ea6\:91cf\:8bc6\:522b ",
      "Bravais/BZ \:7c7b\:578b\:ff0c\:518d\:8bfb\:53d6\:5185\:7f6e\:7684 14 \:7c7b\:4e09\:7ef4 Bravais ",
      "\:6676\:683c\:300125 \:79cd BZ \:53d8\:4f53\:8def\:5f84\:8868\:3002BravaisType ",
      "\:53ef\:663e\:5f0f\:6307\:5b9a\:4e0e\:8bc6\:522b\:7ed3\:679c\:4e00\:81f4\:7684\:7c7b\:578b\:ff0c",
      "Tolerance \:53ea\:63a7\:5236\:5ea6\:91cf\:6bd4\:8f83\:3002\:5b83\:4e0d\:662f\:7a7a\:95f4\:7fa4\:5c0f",
      "\:7fa4\:6216\:4e0d\:53ef\:7ea6\:8868\:793a\:7b97\:6cd5\:ff1b",
      "\:975e\:6807\:51c6\:57fa\:77e2\:548c\:4e2d\:5fc3\:5316\:6676\:683c\:7ea6\:5b9a\:53ef\:80fd\:9700",
      "\:8981\:5148\:6807\:51c6\:5316\:57fa\:77e2\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "showCrystalStructure"}] ->
    "绘制当前模型的原胞、原子和非零磁矩。",
  FileNameJoin[{"ReferencePages", "Symbols", "showBrillouinZone"}] ->
    "绘制第一 Brillouin 区并标出当前晶格的约定高对称路径。",
  FileNameJoin[{"ReferencePages", "Symbols", "orbitalTable"}] ->
    StringJoin[
      "orbitalTable \:663e\:793a init \:6216 initfromrep \:56fa\:5b9a\:7684 Hamiltonian ",
      "\:884c\:5217\:57fa\:987a\:5e8f\:ff0c\:5305\:62ec\:683c\:70b9\:3001Wyckoff ",
      "\:8f68\:9053\:3001\:7b49\:4ef7\:539f\:5b50\:3001\:5206\:6570/\:7b1b\:5361\:5c14\:4f4d\:7f6e\:3001",
      "\:5c40\:57df\:8f68\:9053\:3001\:5b9e\:9645\:57fa\:6001\:3001\:81ea\:65cb\:7ed3\:6784\:548c\:8fd0",
      "\:8f93\:64cd\:4f5c\:3002\:516c\:5f00\:8fd4\:56de\:503c\:662f\:5c55\:793a\:7528 Grid\:ff1b",
      "\:7a0b\:5e8f\:5185\:90e8\:4ecd\:4ee5 Association/List \:4fdd\:5b58\:6570\:636e\:ff0c",
      "\:4f46\:8fd9\:4e0d\:662f\:7b2c\:4e8c\:79cd\:516c\:5f00\:8fd4\:56de\:683c\:5f0f\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "hoppingData"}] ->
    StringJoin[
      "hoppingData \:4ece\:5df2\:7ecf\:7531 symham \:6c42\:89e3\:5e76\:7f13\:5b58\:7684\:952e\:5c42\:76f4",
      "\:63a5\:8bfb\:53d6\:6570\:503c H(R)\:ff0c\:4e0d\:4f1a\:5bf9 H(k) ",
      "\:4f5c\:53cd\:5085\:91cc\:53f6\:91cd\:5efa\:3002Hermitian\:3001KernelMethod\:3001ValidationLevel ",
      "\:5fc5\:987b\:4e0e\:6c42\:89e3\:8be5\:952e\:5c42\:65f6\:5b8c\:5168\:4e00\:81f4\:3002\:8fd4\:56de ",
      "Association \:91c7\:7528 H(R)[a,b]=<0,a|H|R,b> \:4e0e H(k)=Sum_R H(R) Exp[2 Pi I k.R] ",
      "\:7ea6\:5b9a\:ff1b\:672a\:6c42\:89e3\:3001\:8d8a\:754c\:6216\:4ecd\:542b\:7b26\:53f7\:65f6\:8fd4",
      "\:56de $Failed\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "transformHoppings"}] ->
    StringJoin[
      "transformHoppings \:7528\:975e\:5947\:5f02 3\:00d73 \:6574\:6570\:77e9\:9635 T \:628a\:5b8c\:6574 ",
      "H(R) \:7cbe\:786e\:53d8\:6362\:5230\:65b0\:539f\:80de\:ff0c\:65b0\:6676\:683c\:884c\:4e3a T ",
      "\:4e58\:65e7\:6676\:683c\:884c\:3002\:7a0b\:5e8f\:7528\:7cbe\:786e Smith ",
      "\:5206\:89e3\:679a\:4e3e\:966a\:96c6\:ff0c\:4e0d\:5f15\:5165\:7a7a\:95f4\:622a\:65ad\:3002",
      "WannierCenters \:548c Lattice \:53ef\:8986\:76d6\:5143\:6570\:636e\:ff0cHermiticityTolerance ",
      "\:68c0\:67e5 H(R)=H(-R)^dagger\:ff1b\:8f68\:9053\:7ef4\:6570\:4e58\:4ee5 Abs[Det[T]]\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "buildBlochHamiltonian"}] ->
    StringJoin[
      "buildBlochHamiltonian \:6309 H(k)=Sum_R H(R) Exp[2 Pi I k.R] ",
      "\:5728\:4e09\:7ef4\:5012\:7a7a\:95f4\:5206\:6570\:5750\:6807 k \:6784\:9020\:5468\:671f Bloch ",
      "Hamiltonian\:3002CellMatrix \:53ef\:5148\:751f\:6210\:7cbe\:786e\:6574\:6570\:8d85\:80de\:ff0c",
      "HermiticityTolerance \:68c0\:67e5\:8f93\:5165\:4e0e\:8f93\:51fa\:ff0cOutput ",
      "\:53ef\:8fd4\:56de\:77e9\:9635\:6216\:5e26\:52a8\:91cf\:3001\:539f\:80de\:4ee3\:8868\:548c\:6b8b",
      "\:5dee\:7684 Data\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "buildRealSpaceHamiltonian"}] ->
    StringJoin[
      "buildRealSpaceHamiltonian 从 H(R) 组装有限体系 SparseArray。",
      "第二个参数既可以是三个正整数尺寸，也可以是非空、无重复的整数晶胞坐标列表；",
      "显式列表保留输入顺序并只支持全 Open 边界。",
      "矩阵行列依次表示（源晶胞,轨道）与（目标晶胞,轨道）。",
      "Output -> \"Data\" 返回 Shape、Cells、CellBounds，以及与矩阵严格对齐的 BasisRecords；",
      "若 hopping 数据含 Lattice 和 WannierCenters，还会给出分数与笛卡尔位置，",
      "供 plotRealSpaceWavefunction 使用。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "plotWilsonLoop"}] ->
    StringJoin[
      "plotWilsonLoop 沿参数路径平移一条倒格矢闭合的 Wilson 回路，",
      "绘制 Wilson 相位或模一 Wannier 中心分支。",
      "start 与 end 的差必须是 2 Pi 乘整数倒格矢；",
      "带标签路径按 ParameterSubdivisions 采样，每个回路按 LoopSubdivisions 计算。",
      "HermitianTolerance、GapTolerance、CovarianceTolerance 与 OverlapTolerance ",
      "原样传给严格 Wilson-loop 检查，不提供近似回退。",
      "Output -> \"Data\" 同时返回采样路径、分支数值、边界标签和 Graphics。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "plotSurfaceSpectrum"}] ->
    StringJoin[
      "plotSurfaceSpectrum 沿二维表面动量路径与均匀能量网格绘制半无限体系的 ",
      "-Im Tr[Gsurface]/Pi。晶胞轴 1、2 张成表面，轴 3 为堆垛方向；",
      "每个网格点都调用现有主层 decimation，不收敛时明确失败。",
      "MomentumSubdivisions、EnergyPoints、CellMatrix、Broadening、Tolerance、",
      "MaxIterations、Surface 与 HermiticityTolerance 控制数值计算，",
      "其余选项控制密度图。",
      "Output -> \"Data\" 返回动量路径、能量、谱权重和 Graphics。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "plotBerryCurvature2D"}] ->
    StringJoin[
      "plotBerryCurvature2D 在有向动量切片上逐点调用 berryCurvature，",
      "绘制占据子空间的有符号 Berry 曲率密度。",
      "Directions 选择两个有序坐标并决定取向，FixedCoordinates 固定其余坐标；",
      "GridSize 决定显示网格，StepSize 决定每个点的 Wilson 小平行四边形。",
      "四个 tolerance 沿用严格 Berry 几何检查。",
      "Output -> \"Data\" 返回坐标网格、曲率数组、最大绝对值和 Graphics。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "plotBerryCurvature3D"}] ->
    StringJoin[
      "plotBerryCurvature3D 绘制三维 Berry 曲率向量场 {F23,F31,F12}。",
      "每个点分别计算有向平面 {2,3}、{3,1}、{1,2}；",
      "箭头方向表示取向，长度与颜色表示模。",
      "MagnitudeThreshold 只隐藏显示中的小向量，不删除 Data 中的完整网格。",
      "图形可调 Axes、AxesLabel、ViewPoint 等选项，",
      "同时固定 Boxed -> False 和 Orthographic 平行投影。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "plotRealSpaceWavefunction"}] ->
    StringJoin[
      "plotRealSpaceWavefunction 把有限体系态映射到 ",
      "buildRealSpaceHamiltonian[...,\"Output\"->\"Data\"] 返回的矩阵对齐 BasisRecords。",
      "球体体积表示概率；Aggregation -> \"Atom\" 合并同一笛卡尔位置的轨道，",
      "\"Orbital\" 保留每个矩阵基态。",
      "PhaseColoring 可显示每个聚合记录的代表复相位。",
      "输入必须含完整位置元数据，态向量长度必须等于 Dimension 且范数非零；",
      "图形固定隐藏坐标轴和盒框，并采用 Orthographic 平行投影。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "buildSlabHamiltonian"}] ->
    StringJoin[
      "buildSlabHamiltonian \:6784\:9020\:6df7\:5408\:7a7a\:95f4 slab\:3001ribbon \:6216 rod ",
      "Hamiltonian\:3002PeriodicDirections \:6307\:5b9a\:4e00\:4e2a\:6216\:4e24\:4e2a\:4fdd\:7559\:8fde",
      "\:7eed Bloch \:52a8\:91cf\:7684\:76f4\:683c\:77e2\:8f74\:ff0c",
      "\:8fd9\:4e9b\:8f74\:7684\:5c3a\:5bf8\:5fc5\:987b\:4e3a 1\:ff1b",
      "\:5176\:4f59\:8f74\:4fdd\:6301\:6709\:9650\:5f00\:8fb9\:754c\:3002CellMatrix ",
      "\:53ef\:91cd\:5b9a\:5411\:6574\:6570\:539f\:80de\:ff0cHermiticityTolerance ",
      "\:68c0\:67e5\:7ed3\:679c\:ff0cOutput \:53ef\:8fd4\:56de\:77e9\:9635\:6216\:5b8c\:6574 Data\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "surfaceGreenFunction"}] ->
    StringJoin[
      "surfaceGreenFunction \:7528\:4e3b\:5c42 decimation \:8ba1\:7b97\:534a\:65e0\:9650\:8868\:9762\:7684",
      "\:63a8\:8fdf Green \:51fd\:6570\:3002\:539f\:80de\:8f74 1\:30012 \:5f20\:6210\:8868\:9762\:ff0c",
      "\:8f74 3 \:4e3a\:5806\:53e0\:65b9\:5411\:ff1b\:957f\:7a0b\:8f74 3 ",
      "\:8dc3\:8fc1\:81ea\:52a8\:7ec4\:6210\:4e3b\:5c42\:3002CellMatrix\:3001Broadening\:3001Tolerance",
      "\:3001MaxIterations\:3001Surface\:3001HermiticityTolerance \:548c Output ",
      "\:5206\:522b\:63a7\:5236\:539f\:80de\:3001\:5c55\:5bbd\:3001\:6536\:655b\:3001\:8868\:9762\:65b9",
      "\:5411\:3001\:5384\:7c73\:68c0\:67e5\:4e0e\:8fd4\:56de\:683c\:5f0f\:ff1b",
      "\:4e0d\:6536\:655b\:65f6\:660e\:786e\:8fd4\:56de $Failed\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "surfaceSpectralFunction"}] ->
    StringJoin[
      "surfaceSpectralFunction \:8fd4\:56de -Im Tr[Gsurface]/Pi\:3002\:5b83\:4e0e surfaceGreenFunction ",
      "\:4f7f\:7528\:540c\:4e00\:8868\:9762\:51e0\:4f55\:3001\:4e3b\:5c42\:548c\:6536\:655b\:89c4\:5219",
      "\:ff0c\:5e76\:63d0\:4f9b CellMatrix\:3001Broadening\:3001Tolerance\:3001MaxIterations\:3001Surface",
      "\:3001HermiticityTolerance \:516d\:4e2a\:9009\:9879\:ff1b\:9700\:8981\:77e9\:9635\:5757\:548c\:6b8b",
      "\:5dee\:65f6\:5e94\:6539\:7528 surfaceGreenFunction \:7684 Data \:8f93\:51fa\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "wilsonLoop"}] ->
    StringJoin[
      "wilsonLoop \:6cbf start \:5230 end \:7684\:6b63\:5411\:95ed\:5408\:5012\:683c\:77e2\:56de\:8def",
      "\:8ba1\:7b97\:5360\:636e\:5b50\:7a7a\:95f4 Wilson \:672c\:5f81\:76f8\:4f4d\:3002end-start ",
      "\:5fc5\:987b\:662f 2 Pi \:4e58\:6574\:6570\:5012\:683c\:77e2\:ff1b\:8f68\:9053\:4e2d\:5fc3 wcc ",
      "\:7ed9\:51fa\:7aef\:70b9 sewing\:3002\:7b97\:6cd5\:4f7f\:7528 SVD ",
      "\:5e7a\:6b63\:5316\:8fde\:63a5\:ff0c\:5e76\:4e25\:683c\:68c0\:67e5\:5384\:7c73\:6027\:3001\:80fd",
      "\:9699\:3001\:7aef\:70b9\:534f\:53d8\:548c\:8fde\:63a5\:5947\:5f02\:503c\:3002PathSubdivisions ",
      "\:4e0e\:56db\:4e2a tolerance \:63a7\:5236\:79bb\:6563\:548c\:95e8\:9650\:ff0cOutput ",
      "\:53ef\:8fd4\:56de Pi \:5355\:4f4d\:76f8\:4f4d\:3001Wannier ",
      "\:4e2d\:5fc3\:3001\:672c\:5f81\:503c\:6216 Data\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "berryPhase"}] ->
    StringJoin[
      "berryPhase \:5bf9\:663e\:5f0f\:91c7\:6837\:7684\:95ed\:5408\:52a8\:91cf\:8def\:5f84\:8ba1\:7b97",
      "\:5360\:636e\:5b50\:7a7a\:95f4\:603b Berry \:76f8\:4f4d\:3002\:7aef\:70b9\:53ef\:76f8\:540c\:ff0c",
      "\:4e5f\:53ef\:76f8\:5dee 2 Pi \:6574\:6570\:5012\:683c\:77e2\:ff1b",
      "\:8f68\:9053\:4e2d\:5fc3\:8d1f\:8d23\:95ed\:5408 sewing\:3002",
      "HermitianTolerance\:3001GapTolerance\:3001CovarianceTolerance \:548c OverlapTolerance ",
      "\:63a7\:5236\:4e25\:683c\:68c0\:67e5\:ff0cOutput \:53ef\:8fd4\:56de\:5f27\:5ea6\:76f8\:4f4d\:3001Pi ",
      "\:5355\:4f4d\:76f8\:4f4d\:3001Wilson \:884c\:5217\:5f0f\:6216 Data\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "berryCurvature"}] ->
    StringJoin[
      "berryCurvature \:7531 k \:70b9\:9644\:8fd1\:6709\:5411 Wilson ",
      "\:5c0f\:5e73\:884c\:56db\:8fb9\:5f62\:7684 Berry \:901a\:91cf\:6c42\:5360\:636e\:5b50\:7a7a\:95f4",
      "\:66f2\:7387\:3002Directions \:7684\:987a\:5e8f\:51b3\:5b9a\:53d6\:5411\:ff0cStepSize ",
      "\:7ed9\:4e00\:6761\:6216\:4e24\:6761\:8fb9\:957f\:ff1b\:5176\:4f59\:56db\:4e2a tolerance ",
      "\:6cbf\:7528 Berry \:76f8\:4f4d\:68c0\:67e5\:3002Output \:53ef\:8fd4\:56de\:66f2\:7387\:3001\:901a",
      "\:91cf\:6216\:542b\:91c7\:6837\:8def\:5f84\:7684 Data\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "pointChernNumber"}] ->
    StringJoin[
      "pointChernNumber \:5728\:4e09\:7ef4\:52a8\:91cf\:70b9\:5468\:56f4\:7684\:5b9a\:5411\:7acb\:65b9",
      "\:4f53\:8868\:9762\:79ef\:5206\:5360\:636e\:5b50\:7a7a\:95f4\:901a\:91cf\:ff0c",
      "\:8fd4\:56de\:6574\:6570 Chern \:7535\:8377\:3002H \:53ef\:4e3a\:6570\:503c\:51fd\:6570\:6216\:542b ",
      "MomentumSymbols \:7684\:77e9\:9635\:3002SurfaceSubdivisions \:63a7\:5236\:7f51\:683c\:ff1b",
      "HermitianTolerance\:3001SurfaceGapTolerance\:3001CenterGapTolerance\:3001OverlapTolerance \:4e0e ",
      "IntegerTolerance \:63a7\:5236\:95e8\:9650\:ff1bRequireGaplessCenter ",
      "\:51b3\:5b9a\:662f\:5426\:8981\:6c42\:4e2d\:5fc3\:95ed\:5408\:80fd\:9699\:ff1bOutput ",
      "\:53ef\:8fd4\:56de\:6574\:6570\:6216 Data\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "findGaplessPoints"}] ->
    StringJoin[
      "findGaplessPoints \:5148\:5728\:5468\:671f BrillouinZone \:7f51\:683c\:4e0a\:8ba1\:7b97\:76f4\:63a5",
      "\:80fd\:9699\:ff0c\:518d\:4ece\:4e0d\:8d85\:8fc7 CandidateCount ",
      "\:4e2a\:79cd\:5b50\:5c40\:57df\:6781\:5c0f\:5316 gap^2\:3002",
      "GridSize\:3001GapTolerance\:3001MergeTolerance\:3001HermitianTolerance\:3001MaxIterations\:3001Refin",
      "ementMethod \:548c MomentumSymbols \:63a7\:5236\:641c\:7d22\:ff0cOutput -> \"Data\" ",
      "\:7ed9\:51fa\:8bca\:65ad\:5e76\:660e\:786e CompletenessGuaranteed -> False\:ff1b",
      "\:8fd9\:662f\:6709\:754c\:641c\:7d22\:ff0c\:4e0d\:5ba3\:79f0\:6570\:5b66\:5b8c\:5907\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "getMSGElemFromMSGCorep"}] ->
    StringJoin[
      "\:7528\:6237\:5fc5\:987b\:5148\:663e\:5f0f\:6267\:884c Needs[\"MSGCorep`\"]\:3002",
      "getMSGElemFromMSGCorep \:53ea\:628a\:5df2\:8f7d\:5165 MSGCorep \:4ea7\:751f\:7684\:6709\:5e8f\:78c1\:7a7a\:95f4\:7fa4\:64cd\:4f5c",
      "\:8f6c\:6362\:4e3a MagneticTB \:8bb0\:5f55\:ff0c\:4e0d\:4f1a\:81ea\:52a8\:52a0\:8f7d MSGCorep\:3002",
      "\:8be5\:63a5\:53e3\:53ea\:9002\:7528\:4e8e MSGCorep \:771f\:5b9e\:751f\:6210\:7684\:64cd\:4f5c\:ff1b",
      "\:4f9d\:8d56\:672a\:9884\:5148\:8f7d\:5165\:65f6\:8fd4\:56de $Failed\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "getTBBandCorep"}] ->
    StringJoin[
      "getTBBandCorep \:53ea\:9002\:7528\:4e8e getMSGElemFromMSGCorep ",
      "\:751f\:6210\:5e76\:6309\:539f\:987a\:5e8f\:4f20\:7ed9 init \:7684\:5bf9\:79f0\:64cd\:4f5c\:3002",
      "\:6807\:91cf s \:57fa\:8d70\:5355\:503c\:8def\:7ebf\:ff0c{sup,sdn} ",
      "\:65cb\:91cf\:57fa\:8d70\:53cc\:503c\:8def\:7ebf\:3002\:51fd\:6570\:5728\:6bcf\:4e2a\:4e09\:5206",
      "\:91cf\:5012\:7a7a\:95f4\:5206\:6570\:5750\:6807 k \:70b9\:6253\:5370\:540d\:79f0\:548c\:4e00\:4e2a ",
      "Grid\:ff0c\:8fd4\:56de Null\:ff1b\:4f9d\:8d56\:3001\:4f1a\:8bdd\:3001H\:3001\:53c2\:6570\:89c4\:5219",
      "\:3001kset \:6216\:5916\:90e8\:8ba1\:7b97\:5931\:8d25\:65f6\:8fd4\:56de $Failed\:ff0c",
      "\:4e0d\:4f7f\:7528\:5185\:90e8\:8fd1\:4f3c\:66ff\:4ee3\:3002\:8c03\:7528\:524d\:5fc5\:987b\:7531\:7528\:6237\:663e\:5f0f\:6267\:884c Needs[\"MSGCorep`\"]\:3002"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "vaspEig"}] ->
    StringJoin[
      "vaspEig 从 VASP EIGENVAL 文件读取指定自旋通道和闭区间能带，",
      "并从每个本征值中减去给定的费米能。返回值按 k 点排列，",
      "每条记录为 {{kx,ky,kz},{energy1,energy2,...}}。",
      "示例中的 vaspBands 只是保存该列表的普通变量，不是 MagneticTB 函数。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "bandManipulateEig"}] ->
    StringJoin[
      "bandManipulateEig 把符号 MagneticTB Hamiltonian 与 vaspEig 读取的 VASP 能带",
      "放入同一个交互图中。紧束缚能带为黑色，VASP 能带为红色；",
      "每个 Hamiltonian 参数都有滑块。函数按 2 Pi k 计算 H，",
      "并返回 Manipulate。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "fittingTB"}] ->
    StringJoin[
      "fittingTB 用有限差分 Levenberg-Marquardt 算法，把关于参数为仿射线性的 ",
      "MagneticTB Hamiltonian 拟合到数值参考能带。函数缓存常数矩阵与参数系数矩阵，",
      "并返回拟合参数、收敛信息、初始与最终残差以及能带对比图。"
    ],
  FileNameJoin[{"ReferencePages", "Symbols", "compareBand"}] ->
    StringJoin[
      "compareBand 沿给定倒空间分数坐标路径，把已经固定参数的紧束缚能带与",
      "参考能带叠加显示。紧束缚能带为紫色，参考能带为蓝色；plotRange 只控制显示范围。"
    ],
  FileNameJoin[{"Tutorials", "BandFitting"}] ->
    StringJoin[
      "本教程用随页提供的 VASP EIGENVAL 中三条能带，",
      "比较最近邻与第二近邻 MagneticTB 模型，",
      "并演示如何只拟合费米能级或指定 k 点附近的数据。",
      "判断以完整能带对比图为准，不以损失数值为准。"
    ],
  FileNameJoin[{"Tutorials", "CrystalAndKPaths"}] ->
    StringJoin[
      "本教程用每个原胞含三个磁性格点的 120° 非共线反铁磁 Kagome 模型，",
      "串联晶体结构、六角 Brillouin 区、高对称路径和三带 Hamiltonian。"
    ],
  FileNameJoin[{"Tutorials", "RealSpaceAndTopology"}] ->
    StringJoin[
      "本教程把四个由对称性生成的 Hamiltonian 各构造一次。",
      "模型构造默认折叠，正文中的每幅图只保留一条直接的 MagneticTB 调用。"
    ],
  FileNameJoin[{"Tutorials", "Corepresentations"}] ->
    StringJoin[
      "本教程以磁空间群 164.87 为例，",
      "比较旋量基与标量 s 基在 Γ、K、A 点的能带共表示，并显示实际计算表。"
    ]
|>;

(* Section headings are translated separately from page summaries so the
   same English tutorial cell can remain a stable navigation anchor. *)
testingChineseHeadingTranslations = <|
  "Crystal structure, Brillouin zone, and k paths" -> "\:6676\:4f53\:7ed3\:6784\:3001Brillouin \:533a\:4e0e k \:8def\:5f84",
  "Real-space Hamiltonians, surfaces, and topology" -> "\:5b9e\:7a7a\:95f4 Hamiltonian\:3001\:8868\:9762\:4e0e\:62d3\:6251",
  "Band corepresentations with MSGCorep" -> "\:4f7f\:7528 MSGCorep \:8ba1\:7b97\:80fd\:5e26\:5171\:8868\:793a",
  "Initialize the model" -> "\:521d\:59cb\:5316\:6a21\:578b",
  "Crystal structure" -> "\:6676\:4f53\:7ed3\:6784",
  "Symmetry-generated band" -> "\:5bf9\:79f0\:6027\:751f\:6210\:7684\:80fd\:5e26",
  StringJoin[
    "This tutorial uses one complete simple-cubic model for the orbital table, crystal ",
    "structure, conventional k path, first Brillouin zone, and band structure. The Hamiltonian ",
    "is generated by symham from gray[221]; no separate hand-written dispersion is introduced."
  ] -> StringJoin[
    "\:672c\:6559\:7a0b\:7528\:540c\:4e00\:4e2a\:5b8c\:6574\:7684\:7b80\:5355\:7acb\:65b9\:6a21\:578b\:4e32\:8054\:8f68",
    "\:9053\:8868\:3001\:6676\:4f53\:7ed3\:6784\:3001\:6807\:51c6 k \:8def\:5f84\:3001\:7b2c\:4e00 Brillouin ",
    "\:533a\:548c\:80fd\:5e26\:7ed3\:6784\:3002Hamiltonian \:7531 gray[221] \:521d\:59cb\:5316\:540e\:7684 ",
    "symham \:751f\:6210\:ff0c\:4e0d\:518d\:5f15\:5165\:72ec\:7acb\:7684\:624b\:5199\:8272\:6563\:5173\:7cfb\:3002"
  ],
  "showCrystalStructure draws the initialized simple-cubic cell with a parallel orthographic projection by default." ->
    "showCrystalStructure \:7ed8\:5236\:5df2\:521d\:59cb\:5316\:7684\:7b80\:5355\:7acb\:65b9\:6676\:80de\:ff0c\:9ed8\:8ba4\:4f7f\:7528 Orthographic \:5e73\:884c\:6295\:5f71\:3002",
  "standardKPath recognizes the lattice metric and reads a built-in Bravais-lattice path table. It does not compute space-group little-group irreducible representations." ->
    "standardKPath \:8bc6\:522b\:6676\:683c\:5ea6\:91cf\:5e76\:8bfb\:53d6\:5185\:7f6e Bravais \:6676\:683c\:8def\:5f84\:8868\:ff1b\:5b83\:4e0d\:8ba1\:7b97\:7a7a\:95f4\:7fa4\:5c0f\:7fa4\:4e0d\:53ef\:7ea6\:8868\:793a\:3002",
  "Use the second prepared bond shell to generate the nearest-neighbor simple-cubic Hamiltonian, then draw it on the same conventional path." ->
    "\:4f7f\:7528\:5df2\:51c6\:5907\:7684\:7b2c 2 \:4e2a\:8fd1\:90bb hopping \:5c42\:751f\:6210\:7b80\:5355\:7acb\:65b9\:6676\:6700\:8fd1\:90bb Hamiltonian\:ff0c\:518d\:6cbf\:540c\:4e00\:6761\:6807\:51c6\:8def\:5f84\:7ed8\:5236\:80fd\:5e26\:3002",
  "Conventional path and first Brillouin zone" -> "\:7ea6\:5b9a\:8def\:5f84\:4e0e\:7b2c\:4e00 Brillouin \:533a",
  "Static bands" -> "\:9759\:6001\:80fd\:5e26",
  "Hamiltonians used below" -> "下面使用的 Hamiltonian",
  "Chern-insulator Hamiltonian" -> "Chern 绝缘体 Hamiltonian",
  "Z2-insulator Hamiltonian" -> "Z2 绝缘体 Hamiltonian",
  "C4 higher-order topological-insulator Hamiltonian" -> "C4 高阶拓扑绝缘体 Hamiltonian",
  "Magnetic Weyl Hamiltonian" -> "磁性 Weyl Hamiltonian",
  "C4 TCI surface Green function" -> "C4 TCI 表面 Green 函数",
  "Wilson-loop plots" -> "Wilson-loop 图",
  "Surface spectra" -> "表面谱",
  "Chern insulator" -> "Chern 绝缘体",
  "Z2 insulator" -> "Z2 绝缘体",
  "C4 higher-order topological-insulator corner state" -> "C4 高阶拓扑绝缘体角态",
  "Magnetic Weyl point" -> "磁性 Weyl 点",
  "Berry-curvature density" -> "Berry 曲率密度",
  StringJoin[
    "Construct each Hamiltonian once, then pass it directly to a MagneticTB plotting function. ",
    "The model-construction groups are collapsed so that the main workflow remains one or two ",
    "lines per figure; expand a group only when you want to inspect or rerun that Hamiltonian."
  ] -> StringJoin[
    "每个 Hamiltonian 只构造一次，然后直接传给 MagneticTB 绘图函数。",
    "模型构造组默认折叠，因此每幅图的主流程只需一两行；",
    "只有在需要检查或重新生成 Hamiltonian 时才展开相应组。"
  ],
  "Evaluate a model group once before running the figures that use it." ->
    "运行相应图形前，只需计算一次该模型组。",
  "One call gives the winding Chern branch; the second gives the Z2 partner exchange." ->
    "第一条调用给出绕行的 Chern 分支，第二条调用给出 Z2 配对交换。",
  "plotSurfaceSpectrum uses the H(R) data returned by hoppingData and draws the semi-infinite surface spectrum directly." ->
    "plotSurfaceSpectrum 直接使用 hoppingData 返回的 H(R) 数据绘制半无限表面谱。",
  "After the finite Hamiltonian and corner state are prepared, one call draws the real-space probability density." ->
    "准备好有限尺寸 Hamiltonian 和角态后，一条调用即可绘制实空间概率密度。",
  "The first line evaluates the node charge; the second draws its Berry-curvature monopole field." ->
    "第一行计算节点电荷，第二行绘制其 Berry 曲率单极场。",
  "The same Chern Hamiltonian produces the signed two-dimensional curvature map in one call." ->
    "同一个 Chern Hamiltonian 只需一条调用即可生成带符号的二维曲率图。",
  StringJoin[
    "This tutorial follows five symmetry-generated models from real-space hopping data to ",
    "surface and topological plots. It uses a C4 topological crystalline insulator, Chern and ",
    "Z2 insulators, a two-dimensional C4 higher-order topological insulator, and a magnetic Weyl ",
    "model. No hand-written stand-in Hamiltonian is used."
  ] -> StringJoin[
    "本教程用五个对称性生成模型串联实空间 hopping 数据、表面计算与拓扑绘图：",
    "C4 拓扑晶体绝缘体、Chern 与 Z2 绝缘体、二维 C4 高阶拓扑绝缘体和磁性 Weyl 模型。",
    "页面不使用手写替代 Hamiltonian。"
  ],
  StringJoin[
    "The four-band tetragonal model is generated from gray[75]. hoppingData supplies the H(R) ",
    "blocks used by both the semi-infinite Green function and the plotted surface spectrum."
  ] -> "四带四方模型由 gray[75] 生成；hoppingData 给出半无限 Green 函数和表面谱图共同使用的 H(R) 块。",
  StringJoin[
    "The first plot shows one winding Chern branch. The second shows the partner exchange of the ",
    "two occupied branches in a time-reversal-invariant Z2 model."
  ] -> "第一幅图显示绕行一次的 Chern 分支；第二幅图显示时间反演不变 Z2 模型两条占据分支的配对交换。",
  "The corresponding semi-infinite calculations expose the chiral Chern branch and the counterpropagating Z2 pair." ->
    "对应的半无限计算分别显示手征 Chern 分支和反向传播的 Z2 分支对。",
  StringJoin[
    "The open two-dimensional breathing-square model is generated from gray[75]. Its four nearly ",
    "degenerate corner modes are projected onto one corner orbital before the probability density ",
    "is drawn by plotRealSpaceWavefunction."
  ] -> "开边界二维呼吸方格模型由 gray[75] 生成。先把一个角轨道投影到四个近简并角模子空间，再用 plotRealSpaceWavefunction 绘制概率密度。",
  StringJoin[
    "The BNS 143.3 model is generated by symham. pointChernNumber returns the quantized charge on ",
    "a closed cube, and plotBerryCurvature3D draws the same monopole field around the node."
  ] -> "BNS 143.3 模型由 symham 生成。pointChernNumber 返回闭合立方体上的量子化电荷，plotBerryCurvature3D 绘制节点周围同一个单极场。",
  "A symmetry-generated P4 Chern model gives a signed two-dimensional Berry-curvature density plot." ->
    "对称性生成的 P4 Chern 模型给出带符号的二维 Berry 曲率密度图。",
  "One-dimensional hopping data" -> "\:4e00\:7ef4\:8dc3\:8fc1\:6570\:636e",
  "Finite, periodic, and semi-infinite forms" -> "\:6709\:9650\:3001\:5468\:671f\:4e0e\:534a\:65e0\:9650\:5f62\:5f0f",
  "Wilson and Berry geometry" -> "Wilson \:4e0e Berry \:51e0\:4f55",
  "Point charge and gapless-point search" -> "\:70b9\:7535\:8377\:4e0e\:65e0\:9699\:70b9\:641c\:7d22",
  "Arbitrary finite shapes and wavefunctions" -> "任意有限形状与波函数",
  "Topological figures at a glance" -> "拓扑图一览",
  "Chern-insulator Wilson flow" -> "Chern 绝缘体 Wilson 谱流",
  "Z2-insulator Wilson flow" -> "Z2 绝缘体 Wilson 谱流",
  "Chern-insulator surface spectrum" -> "Chern 绝缘体表面谱",
  "Z2-insulator surface spectrum" -> "Z2 绝缘体表面谱",
  "Berry-curvature density and vector plots" -> "Berry 曲率密度图与矢量图",
  "Plot the occupied-band Berry curvature of the Qi-Wu-Zhang Chern model across the full Brillouin zone:" ->
    "绘制 Qi-Wu-Zhang Chern 模型占据带在完整 Brillouin 区内的 Berry 曲率：",
  "Plot the radial {F23,F31,F12} Berry-curvature field around a two-band Weyl point:" ->
    "绘制两带 Weyl 点周围径向的 {F23,F31,F12} Berry 曲率场：",
  "Plot the probability and complex phase of a wavefunction on a three-site L-shaped finite system:" ->
    "绘制三格点 L 形有限系统上波函数的概率与复相位：",
  "Plot a representative corner-localized state of a two-dimensional BBH higher-order topological insulator:" ->
    "绘制二维 BBH 高阶拓扑绝缘体的一个代表性角局域态：",
  "BBH quadrupole-insulator corner state" -> "BBH 四极矩绝缘体角态",
  StringJoin[
    "The second example implements the two-dimensional ",
    "Benalcazar-Bernevig-Hughes (BBH) quadrupole insulator of ",
    "arXiv:1611.07987 with |gamma/lambda| = 0.2 < 1. It projects one corner ",
    "orbital onto the four lowest-absolute-energy finite-system states to ",
    "select a reproducible representative corner mode. The plot visualizes ",
    "localization; it does not by itself calculate the quadrupole invariant."
  ] -> StringJoin[
    "第二个例子实现 arXiv:1611.07987 的二维 Benalcazar-Bernevig-Hughes（BBH）",
    "四极矩绝缘体，并取 |gamma/lambda| = 0.2 < 1。例子把一个角轨道投影到",
    "有限系统中绝对能量最低的四态子空间，从而确定一个可重现的代表性角模。",
    "该图用于显示局域性；它本身不计算四极矩拓扑不变量。"
  ],
  StringJoin[
    "For |gamma/lambda| < 1 the open BBH quadrupole model has four ",
    "mid-gap corner modes. Projecting a corner orbital onto that four-state ",
    "subspace selects one representative mode without relying on an ",
    "arbitrary eigenvector inside the nearly degenerate subspace."
  ] -> StringJoin[
    "当 |gamma/lambda| < 1 时，开边界 BBH 四极矩模型有四个能隙中角模。",
    "把一个角轨道投影到该四态子空间，可选出一个代表性角模，同时避免依赖",
    "近简并子空间中任意选择的本征向量。"
  ],
  "Use one ordered MSGCorep operation list" -> "\:4f7f\:7528\:540c\:4e00\:4efd\:6709\:5e8f MSGCorep \:64cd\:4f5c\:5217\:8868",
  "Double-valued high-symmetry calculation" -> "\:53cc\:503c\:9ad8\:5bf9\:79f0\:70b9\:8ba1\:7b97",
  "Single-valued basis" -> "\:5355\:503c\:57fa",
  "Real-space, surface, and topological calculations" -> "\:5b9e\:7a7a\:95f4\:3001\:8868\:9762\:4e0e\:62d3\:6251\:8ba1\:7b97",
  "Optional MSGCorep interface" -> "\:53ef\:9009 MSGCorep \:63a5\:53e3",
  StringJoin[
    "These four plots are generated by the current MagneticTB API, not imported images. ",
    "The first row shows Wilson spectra as unjoined sampled points. The second row shows ",
    "the matching semi-infinite boundary spectral weight."
  ] -> StringJoin[
    "下面四幅图由当前 MagneticTB 正式接口直接生成，并非导入图片。",
    "前两幅以不连线采样点显示 Wilson 谱；后两幅显示与之对应的半无限边界谱权重。"
  ],
  StringJoin[
    "The following plotSurfaceSpectrum examples use exact nearest-neighbour H(R) blocks. ",
    "Their Graphics outputs are generated by the current MagneticTB plotting API."
  ] -> StringJoin[
    "下面的 plotSurfaceSpectrum 示例使用精确的最近邻 H(R) 块；",
    "图形输出由当前 MagneticTB 绘图接口直接生成。"
  ],
  StringJoin[
    "The following plotWilsonLoop examples show the sampled Wilson spectra as unjoined points. ",
    "Each Graphics output is generated directly by the current MagneticTB plotting API."
  ] -> StringJoin[
    "下面的 plotWilsonLoop 示例把采样 Wilson 谱显示为不连线的离散点；",
    "每幅图形输出都由当前 MagneticTB 绘图接口直接生成。"
  ],
  "The occupied band of the Qi-Wu-Zhang Chern model winds once across the transverse Brillouin zone." ->
    "Qi-Wu-Zhang Chern 模型的占据带在横向 Brillouin 区内绕行一次。",
  StringJoin[
    "The two occupied branches of the lattice BHZ model exchange Kramers partners between ",
    "the two time-reversal-invariant endpoints of half the Brillouin zone."
  ] -> StringJoin[
    "晶格 BHZ 模型的两条占据分支在半个 Brillouin 区的两个时间反演不变端点之间",
    "交换 Kramers 配对。"
  ],
  "The QWZ boundary spectrum contains one chiral branch crossing the bulk gap." ->
    "QWZ 边界谱含有一条穿过体能隙的手征分支。",
  "The lattice BHZ boundary spectrum contains a time-reversed pair of helical branches." ->
    "晶格 BHZ 边界谱含有一对互为时间反演的螺旋分支。",
  StringJoin[
    "The Chern-insulator example has one occupied band whose Wilson phase winds once over a full ",
    "transverse Brillouin zone. The lattice BHZ example has two occupied bands; their Wannier-center ",
    "branches exchange partners between the two time-reversal-invariant endpoints of half the zone. ",
    "plotWilsonLoop displays this flow but does not return a Z2 integer."
  ] -> StringJoin[
    "Chern 绝缘体只有一个占据带，其 Wilson 相位在横向完整 Brillouin 区内绕行一次。",
    "晶格 BHZ 模型有两个占据带；在半个 Brillouin 区的两个时间反演不变端点之间，",
    "Wannier 中心分支发生配对交换。plotWilsonLoop 显示这种谱流，但不直接返回 Z2 整数。"
  ],
  "The occupied band of the Qi-Wu-Zhang Chern insulator has one Wilson branch that winds once across the transverse Brillouin zone:" ->
    "Qi-Wu-Zhang Chern 绝缘体的占据带只有一条 Wilson 分支，它在横向 Brillouin 区内绕行一次：",
  "For the lattice BHZ quantum-spin-Hall model, the two occupied Wannier-center branches exchange Kramers partners between 0 and Pi:" ->
    "对于晶格 BHZ 量子自旋 Hall 模型，两条占据 Wannier 中心分支在 0 与 Pi 之间交换 Kramers 配对：",
  StringJoin[
    "The topological examples place the two-dimensional model in cell axes 1 and 3, so axis 1 is the ",
    "conserved edge momentum and axis 3 is semi-infinite. The Chern insulator has one chiral branch ",
    "crossing the bulk gap. The lattice BHZ model has a time-reversed pair of counterpropagating ",
    "helical edge branches."
  ] -> StringJoin[
    "拓扑例把二维模型放在晶胞轴 1、3 平面内，因此轴 1 是守恒的边缘动量，轴 3 取半无限。",
    "Chern 绝缘体有一条穿过体能隙的手征边缘分支；晶格 BHZ 模型有一对互为时间反演、",
    "传播方向相反的螺旋边缘分支。"
  ],
  "The Qi-Wu-Zhang Chern insulator has a single chiral edge branch that crosses the bulk gap:" ->
    "Qi-Wu-Zhang Chern 绝缘体有一条穿过体能隙的手征边缘分支：",
  "The lattice BHZ quantum-spin-Hall model has a time-reversed pair of helical edge branches crossing at zero momentum:" ->
    "晶格 BHZ 量子自旋 Hall 模型有一对互为时间反演的螺旋边缘分支，并在零动量处相交：",
  StringJoin[
    "The Qi-Wu-Zhang Chern model has one winding Wilson branch and one chiral edge mode. The lattice ",
    "BHZ model has two occupied bands: its Wannier centers exchange Kramers partners over half the ",
    "Brillouin zone, and its boundary has a time-reversed pair of helical modes. plotWilsonLoop shows ",
    "the branch flow but does not return a Z2 integer. plotSurfaceSpectrum calculates the matching ",
    "semi-infinite boundary spectra from the exact nearest-neighbour H(R) blocks."
  ] -> StringJoin[
    "Qi-Wu-Zhang Chern 模型有一条绕行的 Wilson 分支和一个手征边缘模。晶格 BHZ 模型有两个",
    "占据带：Wannier 中心在半个 Brillouin 区内交换 Kramers 配对，其边界出现一对互为时间反演",
    "的螺旋模。plotWilsonLoop 显示分支谱流，但不直接返回 Z2 整数；plotSurfaceSpectrum 从精确",
    "最近邻 H(R) 块计算相应的半无限边界谱。"
  ]
|>;

(* Exact prose for pages whose examples carry physical interpretation.  The
   generic testing-page summary remains a fallback, but it must not erase a
   real translation of notes, option semantics, or example captions. *)
testingChinesePageTranslations = <|
  FileNameJoin[{"Tutorials", "Corepresentations"}] -> <|
    StringJoin[
      "To identify the symmetry of a band at a high-symmetry k point, use its small ",
      "corepresentation. For magnetic space group 164.87, the examples below compare spinor ",
      "and scalar s bases at Gamma, K, and A."
    ] -> StringJoin[
      "要判断高对称 k 点处能带的对称性，可以计算它所属的小共表示。下面以磁空间群 164.87 ",
      "为例，比较旋量基和标量 s 轨道在 Gamma、K、A 点的结果。"
    ],
    StringJoin[
      "This calculation uses the external MSGCorep package and its SpaceGroupIrep dependency. ",
      "Install both first, then load MSGCorep and MagneticTB in a fresh kernel as shown below. ",
      "MSGCorep must be loaded by the user:"
    ] -> StringJoin[
      "这项计算需要外部程序包 MSGCorep 及其依赖 SpaceGroupIrep。安装好后，在新内核中按",
      "下面的顺序加载 MSGCorep 和 MagneticTB。MSGCorep 需要由用户手动加载："
    ],
    StringJoin[
      "Use getMSGElemFromMSGCorep to obtain the symmetry operations, and pass them to init ",
      "without changing their order. The printed primitive vectors give the lattice for this ",
      "example. The corepresentation calculation requires this operation list; do not replace ",
      "it with a list from msgop."
    ] -> StringJoin[
      "先用 getMSGElemFromMSGCorep 取得对称操作，再按原顺序传给 ",
      "init。函数打印的原胞基矢就是这个例子所用的晶格。共表示计算要求使用这份操作列表，不能换成 ",
      "msgop 给出的另一份列表。"
    ],
    StringJoin[
      "Place {sup,sdn} at the origin to calculate double-valued corepresentations. Here we ",
      "keep only the onsite term, with the same energy e1 for both states, to show how the ",
      "result depends on the chosen basis."
    ] -> StringJoin[
      "在原点放置 {sup,sdn} 两个旋量基态，计算双值共表示。这里仅保留 onsite ",
      "项，两个态的能量同为 e1，用来说明基函数如何决定共表示。"
    ],
    StringJoin[
      "Specify Gamma, K, and A in reciprocal fractional coordinates and call getTBBandCorep. ",
      "The function prints the k-point names and result tables. It returns Null, so no ",
      "additional Out cell is expected:"
    ] -> StringJoin[
      "用倒空间分数坐标给出 Gamma、K、A 三个点，再调用 getTBBandCorep。函数直接打印 ",
      "k 点名称和结果表，返回值为 Null，不会再出现一个额外的 Out 单元："
    ],
    StringJoin[
      "The columns list band indices, energy, degeneracy, and the small corepresentation. With ",
      "e1 = 0, the two bands have zero energy and are doubly degenerate at all three points. ",
      "Their labels are Gamma_4(2), K_4 K_5(2), and A_4(2), respectively; each label belongs ",
      "to its own k point."
    ] -> StringJoin[
      "表格的四列依次为能带编号、能量、简并度和小共表示。取 e1 = 0 ",
      "后，两个态在三个点的能量都为零，且都是二重简并；在 Gamma、K、A 点分别属于 ",
      "Γ₄(2)、K₄K₅(2)、A₄(2)。这些标签各自对应所在的 k 点。"
    ],
    StringJoin[
      "Now replace {sup,sdn} by a scalar s orbital. Use the same ordered symmetry operations, ",
      "run init again, and recalculate the Hamiltonian before finding its corepresentations:"
    ] -> StringJoin[
      "再把 {sup,sdn} 换成一个标量 s 轨道。仍使用同一份有序对称操作，重新 init 并求 ",
      "Hamiltonian，再计算共表示："
    ],
    StringJoin[
      "There is now one nondegenerate band, with labels Gamma_1(1), K_1(1), and A_1(1). Both ",
      "examples use zero onsite energy, but the spinor and scalar bases give double-valued and ",
      "single-valued corepresentations, respectively."
    ] -> StringJoin[
      "这时只有一条非简并能带，在三个点分别属于 Γ₁(1)、K₁(1)、A₁(1)。虽然两个例子的 ",
      "onsite 能量都取零，旋量基给出双值共表示，标量 s 轨道给出单值共表示。"
    ]
  |>,
  FileNameJoin[{"Tutorials", "CrystalAndKPaths"}] -> <|
    StringJoin[
      "Use a kagome antiferromagnet to draw the crystal structure, its magnetic moments, and ",
      "the first Brillouin zone. The three in-plane moments are 120 degrees apart. After ",
      "setting up this model, we obtain a k path and use it to plot the three bands."
    ] -> StringJoin[
      "下面以 Kagome 反铁磁模型为例，介绍如何画晶体结构、磁矩和第一 Brillouin ",
      "区。模型中的三个磁矩位于平面内，两两夹角为 120°。设置好模型后，再选取高对称路径，画出三条能带。"
    ],
    StringJoin[
      "Choose magnetic space group 191.239, P6'/m'mm', with a = 1 and c = 4. One ",
      "representative atom generates {1/2,0,0}, {0,1/2,0}, and {1/2,1/2,0}, in this order. ",
      "Their moments in the lattice-vector basis are {1,0,0}, {0,1,0}, and {-1,-1,0}. In ",
      "Cartesian coordinates, these are unit vectors 120 degrees apart and their sum is zero."
    ] -> StringJoin[
      "选取磁空间群 191.239，P6'/m'mm'，晶格常数取 a = 1、c = ",
      "4。输入一个代表原子，就会依次生成 {1/2,0,0}、{0,1/2,0} 和 {1/2,1/2,0} ",
      "三个位置。以晶格基矢为基，它们的磁矩分别是 {1,0,0}、{0,1,0} 和 ",
      "{-1,-1,0}；换成笛卡尔坐标后，三个单位矢量两两夹角为 120°，合磁矩为零。"
    ],
    StringJoin[
      "Here {1,-1}/Sqrt[2] is one spinor at the representative site, not two independent basis ",
      "states. Induced generates its counterparts at the other two sites. There is one state ",
      "per site, giving a three-band Hamiltonian; keeping both spin states on every site would ",
      "instead give six basis states."
    ] -> StringJoin[
      "其中 {1,-1}/Sqrt[2] 表示代表格点上的一个旋量，不是两个独立基态。Induced ",
      "会生成另外两个格点上的对应态。每个格点保留一个态，所以 Hamiltonian ",
      "有三条带；如果每个格点都保留自旋上、下两个态，才是六维基。"
    ],
    StringJoin[
      "Draw 2 by 2 primitive cells. The picture contains twelve atoms and twelve red ",
      "magnetic-moment arrows. MomentScale changes the arrow length. The gray lines mark cell ",
      "edges, not hopping bonds. The default view uses an orthographic projection."
    ] -> StringJoin[
      "先画出 2 × 2 个原胞，共有十二个原子和十二支红色磁矩箭头。MomentScale ",
      "调整箭头长度，灰线表示晶胞边界，不是跃迁键。默认视图采用平行投影。"
    ],
    StringJoin[
      "View 3 by 3 cells from above to see the repeating 120-degree magnetic pattern. There ",
      "are twenty-seven atoms in this picture. CellRange changes how many cells are drawn, not ",
      "the Hamiltonian:"
    ] -> StringJoin[
      "再从上方看 3 × 3 个原胞，图中二十七个原子清楚地显示了重复的 120° 磁序。CellRange ",
      "只改变画出的晶胞数，不改变 Hamiltonian："
    ],
    StringJoin[
      "Use standardKPath to obtain the conventional path for this hexagonal lattice. The ",
      "returned points are in reciprocal fractional coordinates. The function chooses a ",
      "built-in path from the lattice metric; it does not calculate little-group irreducible ",
      "representations."
    ] -> StringJoin[
      "使用 standardKPath 可以得到这个六方晶格的约定高对称路径，结果采用倒空间分数坐标。函数根据",
      "晶格度量选取内置路径，不计算小群的不可约表示。"
    ],
    StringJoin[
      "For the in-plane bands, keep the first three segments, Gamma-M-K-Gamma at kz = 0. Use ",
      "the same path in the Brillouin-zone and band plots:"
    ] -> StringJoin[
      "这里看平面内的能带，只取前三段，即 kz = 0 平面上的 Gamma-M-K-Gamma。下面的 ",
      "Brillouin 区图和能带图使用同一条路径："
    ],
    StringJoin[
      "The first Brillouin zone is a hexagonal prism, and the colored path is in its central ",
      "plane. In this plot the fractional k points have been converted to Cartesian reciprocal ",
      "coordinates using reciprocal vectors that include 2 Pi."
    ] -> StringJoin[
      "图中的第一 Brillouin 区是六角柱，彩色路径位于中央平面。绘图时，分数坐标已经用包含 2 Pi ",
      "的倒格矢换成了倒空间笛卡尔坐标。"
    ],
    "Hamiltonian and bands" -> "Hamiltonian 与能带",
    StringJoin[
      "Add symham[1] and symham[2] to keep onsite and nearest-neighbor terms. ComplexExpand ",
      "and Simplify rewrite the paired exponentials as cosines for real k. The three rows and ",
      "columns follow the site order given above:"
    ] -> StringJoin[
      "把 symham[1] 与 symham[2] 相加，保留 onsite ",
      "和最近邻跃迁。ComplexExpand 和 Simplify 将实数 k ",
      "下成对的指数项写成余弦。矩阵的三行三列对应前面列出的三个格点："
    ],
    StringJoin[
      "The common onsite energy is e1, and t1 is the nearest-neighbor hopping parameter. The ",
      "relative minus sign in the (2,3) entry comes from the induced spinor basis used here. ",
      "These terms contain no interlayer hopping, so the bands are independent of kz. This is ",
      "an ideal kagome model, not a fit to a particular material."
    ] -> StringJoin[
      "其中 e1 是三个格点共有的 onsite 能量，t1 是最近邻跃迁参数。(2,3) ",
      "矩阵元中的相对负号来自这里使用的诱导旋量基。这两项不含层间跃迁，所以能带与 kz ",
      "无关。这里考虑的是理想 Kagome 模型，不是某种材料的拟合结果。"
    ],
    StringJoin[
      "Set e1 = 0 and t1 = -1 and plot the bands. The flat band has energy -2. The two ",
      "dispersive bands meet at K at energy 1, while the flat band touches the lower ",
      "dispersive band at Gamma. The energy unit is the magnitude of the hopping."
    ] -> StringJoin[
      "取 e1 = 0、t1 = -1，就得到下面的能带。平带位于能量 -2，另外两条色散带在 K 点、能量 ",
      "1 处相交；平带在 Gamma 点与下方色散带相接。这里以跃迁参数的绝对值为能量单位。"
    ]
  |>,
  FileNameJoin[{"Tutorials", "BandFitting"}] -> <|
    "Fitting MagneticTB bands to VASP" ->
      "将 MagneticTB 能带拟合到 VASP",
    StringJoin[
      "To fit a tight-binding model to VASP bands, first read EIGENVAL and construct the ",
      "Hamiltonian from symmetry. Here we fit three bands with a three-orbital model, compare ",
      "nearest-neighbor and second-neighbor results, and then fit a selected energy or k-point ",
      "range. Plot the fitted and VASP bands together to see both the agreement and any ",
      "remaining differences."
    ] -> StringJoin[
      "要用 VASP 能带拟合紧束缚模型，先读取 EIGENVAL，再根据对称性构造 ",
      "Hamiltonian。下面用三轨道模型拟合三条能带，比较最近邻和第二近邻的结果，再介绍只拟合某个能量范",
      "围或 k 点附近的用法。把拟合能带和 VASP 能带画在一起，就可以看出哪些地方符合、哪些地方还有差别。"
    ],
    "Read and inspect the three VASP bands" ->
      "读取并检查三条 VASP 能带",
    StringJoin[
      "Use vaspEig to read bands 21-23 and subtract the Fermi energy. First plot all three ",
      "bands to see the reference band structure:"
    ] -> StringJoin[
      "用 vaspEig 读取第 21 至 23 条能带，并减去费米能级。先把三条带都画出来，查看要拟合的能带",
      "："
    ],
    "Build the symmetry-generated model" ->
      "构造对称性生成的模型",
    StringJoin[
      "For this three-orbital model, take the first six operations from msgop[gray[156]], ",
      "which form its complete unitary subgroup. The onsite and nearest-neighbor terms contain ",
      "e1, e2, and t1-t9. These are the parameters to be fitted:"
    ] -> StringJoin[
      "这里的三轨道模型使用 msgop[gray[156]] ",
      "的前六个操作，它们组成完整的幺正子群。onsite 和最近邻项含有 e1、e2 以及 t1 至 ",
      "t9，这些就是需要拟合的参数："
    ],
    "Fit the onsite and nearest-neighbor model" ->
      "拟合 onsite 与最近邻模型",
    "Fit the three bands at all 90 k points. Then draw them together with the VASP bands:" -> "先拟合全部 90 个 k 点处的三条能带，再与 VASP 能带画在同一张图上：",
    "Add second-neighbor hopping" ->
      "加入第二近邻 hopping",
    StringJoin[
      "To include second-neighbor hopping, keep the fitted nearest-neighbor parameters as the ",
      "initial values and set the new parameters r1-r10 to zero. Fit again and compare the two ",
      "band plots to see where the extra hoppings improve the agreement."
    ] -> StringJoin[
      "加入第二近邻时，用刚才得到的最近邻参数作为初值，新增的 r1 至 r10 ",
      "先设为零。再拟合一次，对比两张能带图，就可以看出第二近邻跃迁改善了哪些部分。"
    ],
    "Fit the energy range of interest" ->
      "只拟合关心的能量范围",
    StringJoin[
      "If only bands near the Fermi energy are of interest, set EnergyWindow to the desired ",
      "energy interval. Only reference energies in this interval enter the fit. The plot still ",
      "shows the full bands, so changes outside the interval can also be seen."
    ] -> StringJoin[
      "如果只关心费米能级附近的能带，可以用 EnergyWindow ",
      "指定能量范围，拟合时只计入这个范围内的参考能量。图中仍画出完整能带，便于同时查看范围以外的变化。"
    ],
    "Fit near the K point" ->
      "只拟合 K 点附近",
    StringJoin[
      "To fit near K, set KPointNeighborhood to K and radius 0.15 in reciprocal fractional ",
      "coordinates. Distances include periodic equivalents of k points. ComparisonPlot ",
      "displays only this neighborhood: blue points are VASP, gray dashed curves are the ",
      "initial model, and red curves are the fitted model."
    ] -> StringJoin[
      "如果只拟合 K 点附近，可以用 KPointNeighborhood 指定 K 点和半径 ",
      "0.15。这里采用倒空间分数坐标，距离计算考虑 k 点的周期等价位置。ComparisonPlot ",
      "只显示这个范围：蓝点是 VASP 能带，灰色虚线是初始模型，红线是拟合结果。"
    ],
    "Judge the plotted bands, not only the loss" ->
      "根据能带图判断，不要只看损失",
    StringJoin[
      "fittingTB compares sorted model eigenvalues with the corresponding reference energies. ",
      "It does not impose a gap or track orbital character and band connectivity. A small ",
      "energy error alone is therefore not enough: also compare gaps and crossings in the ",
      "plot. Use physically reasonable initial parameters, and retain the fitted values when ",
      "adding more neighbor hoppings."
    ] -> StringJoin[
      "fittingTB 比较排序后的模型本征值与对应的参考能量，并不额外约束带隙，也不追踪轨道成分和能带连接",
      "关系。因此，数值误差小还不够，还要看图中的带隙和交叉是否正确。初值应有合理的物理依据；增加更远的近邻跃迁",
      "时，可以沿用已有的拟合参数。"
    ],
    "Related functions" -> "相关函数"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "hoppingData"
  }] -> <|
    "returns numerical real-space hopping blocks for cached shells 1 through n." ->
      "返回已经缓存的第 1 到第 n 个近邻项的数值实空间 hopping 矩阵。",
    "returns numerical blocks for the explicitly selected cached shells." ->
      "返回列表中明确选定且已经缓存的近邻项。",
    "Call init and symham first. hoppingData reads the solved real-space shell cache directly; it never reconstructs hopping amplitudes by Fourier-inverting H(k)." ->
      "先运行 init 和 symham；hoppingData 直接读取已经求出的实空间 H(R)，不会对 H(k) 作反 Fourier 变换。",
    "The returned Association uses H(R)[a,b]=<0,a|H|R,b> and H(k)=Sum_R H(R) Exp[2 Pi I k.R]. Lattice and WannierCenters are included when the current model supplies them." ->
      "返回的 Association 包含晶格平移、H(R) 矩阵，以及模型可提供的晶格和 Wannier 中心。",
    StringJoin[
      "hoppingData returns $Failed and issues a message when the shell selection is invalid or unprepared, ",
      "parameter values remain symbolic, the requested symham option combination was not cached, or cached ",
      "real-space data are inconsistent. It does not repair the input or switch to an approximate fallback."
    ] ->
      "所选近邻尚未求解、编号越界、参数仍含符号，或者不存在匹配的 symham 缓存时，hoppingData 返回 $Failed。",
    "Must match the Hermitian setting used when the selected shells were solved by symham." ->
      "必须与 symham 求解这些近邻时使用的 Hermitian 设置一致。",
    "Must match the exact common-kernel method used by symham." ->
      "必须与 symham 使用的精确公共零空间方法一致。",
    "Must match the validation level used by symham." ->
      "必须与 symham 使用的验证等级一致。",
    "Solve the nearest-neighbor shell of a one-orbital cubic model and read its numerical real-space blocks from the cache:" ->
      "求出单轨道简单立方模型的最近邻项，再读取相应的数值 H(R)："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "getMSGElemFromMSGCorep"
  }] -> <|
    "converts the magnetic-space-group operations of an explicitly loaded MSGCorep application to MagneticTB's ordered operation records." ->
      "把用户已经载入的 MSGCorep 所给出的磁空间群操作转换为 MagneticTB 使用的有序操作记录。",
    "This interface is defined only for operations returned by MSGCorep. It is not a general converter for arbitrary lists that merely resemble MSGCorep data." ->
      "这个接口只处理 MSGCorep 返回的操作，不能用来转换形式相似的任意列表。",
    StringJoin[
      "Each returned record contains the MSGCorep label, its rotation matrix, fractional translation, and ",
      "MagneticTB's \"F\"/\"T\" unitary flag. The function also prints the primitive lattice vectors ",
      "supplied by SpaceGroupIrep."
    ] ->
      "每条记录依次给出 MSGCorep 标签、旋转矩阵、分数平移和 MagneticTB 使用的 \"F\"/\"T\" 标志；函数同时输出 SpaceGroupIrep 给出的原胞基矢。",
    "MSGCorep and SpaceGroupIrep are optional external Mathematica applications. The user must load MSGCorep explicitly with Needs[\"MSGCorep`\"] before calling this function; MagneticTB does not load it automatically." ->
      "MSGCorep 和 SpaceGroupIrep 是可选的外部 Mathematica 应用。调用前必须由用户显式执行 Needs[\"MSGCorep`\"]；MagneticTB 不会自动加载它。",
    "getMSGElemFromMSGCorep returns $Failed and issues a message when MSGCorep or its SpaceGroupIrep dependency has not already been loaded. It does not repair the input or switch to an approximate fallback." ->
      "若 MSGCorep 或其 SpaceGroupIrep 依赖尚未载入，getMSGElemFromMSGCorep 会给出消息并返回 $Failed；函数不会修补输入，也不会改用近似算法。",
    "Load the nontrivial magnetic group {164,87}:" ->
      "先显式载入 MSGCorep，再读取非平凡磁空间群 {164,87}："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "getTBBandCorep"
  }] -> <|
    "prints the MSGCorep band-corepresentation table at every supplied reciprocal-fractional k point." ->
      "在每个给定的倒空间分数坐标 k 点输出 MSGCorep 能带共表示表。",
    "Call getMSGElemFromMSGCorep[MSG], then use exactly that ordered operation list in init. getTBBandCorep is supported only for this MSGCorep-generated route." ->
      "先调用 getMSGElemFromMSGCorep[MSG]，再把它返回的有序操作列表原样传给 init；getTBBandCorep 只支持这条由 MSGCorep 生成操作的流程。",
    "The basis selected by init fixes single- or double-valued treatment: a scalar basis such as s is single valued, while the explicit {sup,sdn} spinor basis is double valued." ->
      "init 所选的局域基决定采用单值还是双值路线：标量 s 基采用单值表示，显式的 {sup,sdn} 旋量基采用双值表示。",
    StringJoin[
      "kset is a nonempty list of three-component reciprocal-fractional points. The function returns Null ",
      "and prints one k-name line and one Grid per point; the Grid, rather than an Association, is the ",
      "public presentation result."
    ] ->
      "kset 是由三分量倒空间分数坐标组成的非空列表。函数对每个 k 点输出一行名称和一个 Grid，并返回 Null；公开展示结果是 Grid，不是 Association。",
    "The optional MSGCorep and SpaceGroupIrep applications must be installed, and the user must execute Needs[\"MSGCorep`\"] before either MagneticTB wrapper is called. MagneticTB neither loads MSGCorep automatically nor replaces its calculation with an internal approximation." ->
      "必须先安装可选的 MSGCorep 与 SpaceGroupIrep，并在调用任一 MagneticTB 包装函数前由用户执行 Needs[\"MSGCorep`\"]。MagneticTB 不会自动加载 MSGCorep，也不会用内部近似替代其计算。",
    "getTBBandCorep returns $Failed and issues a message when the dependency is unavailable, the current init session is incompatible, H/rules/kset are malformed, or MSGCorep rejects the model. It does not repair the input or switch to an approximate fallback." ->
      "依赖未载入、当前 init 模型不兼容、H、rules 或 kset 格式错误，或者 MSGCorep 拒绝该模型时，getTBBandCorep 会给出消息并返回 $Failed；函数不会修补输入，也不会改用近似算法。",
    "For MSG {164,87}, initialize a double-valued spinor model and print corepresentations at Gamma, K, and A:" ->
      "对 MSG {164,87} 初始化双值旋量模型，并输出 Gamma、K 和 A 点的共表示：",
    "The same nontrivial group with one scalar s orbital follows the single-valued MSGCorep route:" ->
      "同一个非平凡磁群改用一个标量 s 轨道时，采用 MSGCorep 的单值路线："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "vaspEig"
  }] -> <|
    "reads a selected spin channel and inclusive band interval from a complete numeric VASP EIGENVAL file." ->
      "从完整的数值 VASP EIGENVAL 文件读取指定自旋通道和闭区间能带。",
    "The result contains one {{kx,ky,kz},{energy1,energy2,...}} record per k point. The reciprocal coordinates remain fractional and every energy is shifted by subtracting efermi." ->
      "每个 k 点返回一条 {{kx,ky,kz},{energy1,energy2,...}} 记录。倒空间坐标保持为分数坐标，每个能量都减去 efermi。",
    "Band numbers are 1 based. startBand and endBand are both included. spin -> 1 selects the first energy column; spin -> 2 is available only in a two-spin EIGENVAL file." ->
      "能带编号从 1 开始，startBand 与 endBand 都包含在区间内。spin -> 1 选择第一列能量；只有双自旋 EIGENVAL 才能使用 spin -> 2。",
    "The parser accepts only the complete numeric VASP layout. It does not evaluate file contents or silently keep a partial record." ->
      "解析器只接受完整的数值 VASP 格式；它不会执行文件内容，也不会悄悄保留不完整记录。",
    "vaspEig returns $Failed and issues a message when the file is unreadable or incomplete, a numeric record is malformed, spin is unavailable, or the inclusive band interval is invalid. It does not repair the input or switch to an approximate fallback." ->
      "文件不可读或不完整、数值记录格式错误、自旋通道不存在或闭区间能带编号非法时，vaspEig 给出消息并返回 $Failed；函数不会修补输入，也不会改用近似结果。",
    "Read bands 21 through 23 from the included EIGENVAL file and show the first k-point record:" ->
      "从随页面提供的 EIGENVAL 文件读取第 21 至 23 条能带，并显示第一个 k 点记录：",
    "Plot the three imported bands before fitting:" ->
      "拟合前先绘制导入的三条能带：",
    "The included sample contains 90 k points and three selected bands:" ->
      "随页面提供的样例含 90 个 k 点和三条选定能带："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "bandManipulateEig"
  }] -> <|
    "opens a Manipulate panel for tuning the symbolic parameters of H against numeric reference bands." ->
      "打开 Manipulate 面板，交互调节 H 的符号参数，使其匹配数值参考能带。",
    "referenceBands uses the {{kx,ky,kz},{energy1,energy2,...}} record format returned by vaspEig. H is evaluated at 2 Pi times each reciprocal-fractional k point." ->
      "referenceBands 使用 vaspEig 返回的 {{kx,ky,kz},{energy1,energy2,...}} 记录格式。H 在每个倒空间分数坐标的 2 Pi 倍处计算。",
    "Every ordinary symbolic parameter in H receives a slider from -1 to 1 with initial value 0. The ExportData button prints the current parameter rules; it does not write a file." ->
      "H 中每个普通符号参数都有一个从 -1 到 1、初值为 0 的滑块。ExportData 按钮输出当前参数规则，但不写文件。",
    "Tight-binding bands are black and reference bands are red. The displayed energy range is the reference range expanded by 0.2 at each end." ->
      "紧束缚能带为黑色，参考能带为红色；显示范围在参考能量范围两端各扩展 0.2。",
    "bandManipulateEig has no public options and returns a Manipulate expression." ->
      "bandManipulateEig 没有公开选项，返回 Manipulate 表达式。",
    "A malformed Hamiltonian, malformed reference data, Subscript or Indexed parameters, or unresolved function heads return a descriptive Failure; the function does not substitute a numerical fallback." ->
      "Hamiltonian 或参考数据格式错误、参数使用 Subscript 或 Indexed、或者含有未解析函数头时，函数返回带说明的 Failure，不会改用数值 fallback。",
    "Read the included VASP bands and construct the symmetry-generated three-orbital hexagonal Hamiltonian:" ->
      "读取随页面提供的 VASP 能带，并用对称性生成六角三轨道 Hamiltonian：",
    "Move the sliders until the black tight-binding bands follow the red VASP bands:" ->
      "调节滑块，使黑色紧束缚能带尽量贴合红色 VASP 能带："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "fittingTB"
  }] -> <|
    "fits every affine-linear symbolic parameter of H to selected numeric reference bands." ->
      "把 H 中所有仿射线性符号参数拟合到选定的数值参考能带。",
    "referenceBands uses the {{kx,ky,kz},{energy1,energy2,...}} format returned by vaspEig. kRange is a nonempty list of 1-based record indices, and every selected record must contain exactly Dimensions[H][[1]] energies." ->
      "referenceBands 使用 vaspEig 返回的 {{kx,ky,kz},{energy1,energy2,...}} 格式。kRange 是从 1 开始的非空记录编号列表；每条选中记录必须恰有 Dimensions[H][[1]] 个能量。",
    "The EnergyWindow, KPointNeighborhood, BandSelection, and ResidualWeights selectors are intersected with kRange. EnergyWindow masks individual reference energies, not complete k-point records." ->
      "EnergyWindow、KPointNeighborhood、BandSelection 和 ResidualWeights 选择器与 kRange 取交集。EnergyWindow 屏蔽单个参考能量，而不是整条 k 点记录。",
    "KPointNeighborhood uses fractional reciprocal coordinates. With Periodic -> True, the displacement is reduced by the componentwise minimum-image convention before applying a spherical Radius or componentwise Range." ->
      "KPointNeighborhood 使用倒空间分数坐标。Periodic -> True 时，先按各分量最小镜像约定约化位移，再应用球形 Radius 或逐分量 Range。",
    "BandSelection uses 1-based reference positions. At every retained k point, sorted model eigenvalue j is paired with reference position j; the reference bands are not dynamically reordered." ->
      "BandSelection 使用从 1 开始的参考能带位置。在每个保留的 k 点，排序后的第 j 个模型本征值与第 j 个参考位置配对；参考能带不会动态重新排序。",
    "ResidualWeights is Automatic or a nonnegative record-by-band matrix with the same shape as referenceBands[[All,2]]. A zero weight removes that residual; positive weights minimize the weighted sum of squared band differences." ->
      "ResidualWeights 可取 Automatic，或取与 referenceBands[[All,2]] 同形状的非负记录-能带矩阵。零权重移除对应残差；正权重用于最小化能带差的加权平方和。",
    "initialRules can be a rule list or Association and must assign a real numeric starting value to every fitting parameter. H is evaluated at 2 Pi times each reciprocal-fractional k point." ->
      "initialRules 可以是规则列表或 Association，且必须给每个拟合参数指定实数初值。H 在每个倒空间分数坐标的 2 Pi 倍处计算。",
    "The Hamiltonian must be affine-linear in the fitting parameters. MagneticTB caches its constant and parameter-coefficient matrices, then uses a finite-difference Levenberg-Marquardt optimizer; it does not switch to a slower nonlinear symbolic fallback." ->
      "Hamiltonian 对拟合参数必须为仿射线性。MagneticTB 先缓存常数矩阵和参数系数矩阵，再使用有限差分 Levenberg-Marquardt 优化；不会切换到更慢的非线性符号 fallback。",
    "The objective minimizes sorted energy residuals only. It does not preserve a band gap, band character, or a chosen band ordering; always inspect ComparisonPlot and reject a lower-loss fit when its band structure is physically wrong." ->
      "目标函数只最小化排序后的能量残差，并不保持带隙、能带特征或指定的带序；必须检查 ComparisonPlot，若较低损失对应的能带结构在物理上错误，就应舍弃该拟合。",
    "The result reports FittedParams, convergence diagnostics, ComparisonPlot, CandidateKPointCount, UsedKPointCount, UsedKPointIndices, ResidualCount, the effective selectors, weighting, and band pairing." ->
      "结果给出 FittedParams、收敛诊断、ComparisonPlot、CandidateKPointCount、UsedKPointCount、UsedKPointIndices、ResidualCount，以及实际使用的选择器、权重和能带配对信息。",
    "fittingTB returns $Failed and issues a message when the Hamiltonian or reference records are invalid, kRange or a selector is invalid, the combined selection is empty or underdetermined, initial values are incomplete, the model is nonlinear in its parameters, the time limit is exceeded, or the optimizer fails. It does not repair the input or switch to an approximate fallback." ->
      "Hamiltonian 或参考记录非法、kRange 或选择器非法、组合选择为空或约束不足、初值不完整、模型对参数非线性、超过时间限制或优化失败时，fittingTB 给出消息并返回 $Failed；函数不会修补输入，也不会改用近似结果。",
    "Positive integer upper bound on Levenberg-Marquardt iterations." ->
      "Levenberg-Marquardt 迭代次数的正整数上限。",
    "Positive fitting time limit in seconds, or Infinity." ->
      "以秒计的正拟合时间上限，也可取 Infinity。",
    "All keeps every reference energy; {emin,emax} keeps energies in the inclusive interval." ->
      "All 保留全部参考能量；{emin,emax} 保留闭区间内的能量。",
    "All keeps every k point; an Association selects a periodic or nonperiodic Radius or componentwise Range around Center." ->
      "All 保留全部 k 点；Association 可选择 Center 周围周期或非周期的 Radius 邻域或逐分量 Range。",
    "All keeps every reference position; an integer, increasing integer list, or positive-step Span selects 1-based band positions." ->
      "All 保留全部参考位置；整数、递增整数列表或正步长 Span 选择从 1 开始的能带位置。",
    "Automatic assigns unit weight; a nonnegative record-by-band matrix weights individual residuals, with zero removing a residual." ->
      "Automatic 赋予单位权重；非负记录-能带矩阵可逐个加权残差，其中零会移除对应残差。",
    "Positive relative step used for central finite-difference Jacobian columns." ->
      "中心有限差分 Jacobian 各列使用的正相对步长。",
    "Positive initial Levenberg-Marquardt damping factor." ->
      "Levenberg-Marquardt 的正初始阻尼因子。",
    "Positive stopping tolerance for the gradient, step size, and objective improvement." ->
      "梯度、步长与目标函数改进量使用的正停止容差。",
    "Generate a one-band cubic Hamiltonian and exact reference bands with known parameters:" ->
      "用对称性生成单带立方 Hamiltonian，并按已知参数产生精确参考能带：",
    "Fit the two parameters from zero initial values:" ->
      "从零初值拟合两个参数：",
    "Plot the reference points, the zero-parameter bands, and the fitted bands:" ->
      "绘制参考数据点、零参数能带和拟合后的能带：",
    "Read the VASP three-band data and construct the symmetry-generated hexagonal model:" ->
      "读取 VASP 三带数据，并构造由对称性生成的六角晶格模型：",
    "Fit all three VASP bands with onsite and nearest-neighbour hopping, then draw the full-range comparison:" ->
      "用 onsite 和最近邻 hopping 拟合全部三条 VASP 能带，并绘制完整能量范围的对比图：",
    "Add second-neighbour hopping, refit the same three VASP bands, and draw the full-range comparison again:" ->
      "加入第二近邻 hopping，重新拟合同样的三条 VASP 能带，并再次绘制完整能量范围的对比图：",
    "\"EnergyWindow\" - fit only reference energies from -0.5 through 0.5:" ->
      "\"EnergyWindow\" - 只拟合 -0.5 到 0.5 之间的参考能量：",
    "\"KPointNeighborhood\" - fit the periodic radius 0.35 around Gamma in fractional reciprocal coordinates:" ->
      "\"KPointNeighborhood\" - 在倒空间分数坐标中拟合 Gamma 点周围半径 0.35 的周期邻域：",
    "\"BandSelection\" - select the first sorted reference position:" ->
      "\"BandSelection\" - 选择排序后的第一个参考能带位置：",
    "\"ResidualWeights\" - give the first path segment five times the weight of the remaining points:" ->
      "\"ResidualWeights\" - 令第一段路径的权重为其余 k 点的五倍：",
    "MaxIterations - stop this example after one optimizer iteration:" ->
      "MaxIterations - 让这个示例在一次优化迭代后停止：",
    "TimeConstraint - remove the elapsed-time limit while retaining the iteration limit:" ->
      "TimeConstraint - 保留迭代上限，但取消运行时间限制：",
    "\"FiniteDifferenceStep\" - use a larger relative step for the central-difference Jacobian:" ->
      "\"FiniteDifferenceStep\" - 为中心差分 Jacobian 使用更大的相对步长：",
    "\"InitialDamping\" - start from a stronger damping factor:" ->
      "\"InitialDamping\" - 从更强的阻尼因子开始：",
    "\"FitTolerance\" - use a looser stopping tolerance:" ->
      "\"FitTolerance\" - 使用较宽松的停止容差："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "compareBand"
  }] -> <|
    "overlays numerical tight-binding bands with supplied reference bands along a labeled path." ->
      "沿带标签的路径叠加数值紧束缚能带与输入的参考能带。",
    "path is a nonempty list of numeric three-component reciprocal-fractional segment endpoints and endpoint labels. npoint is a positive integer." ->
      "path 是非空列表，每段含数值三分量倒空间分数坐标端点及其标签；npoint 必须为正整数。",
    "referenceBands uses the vaspEig record format and must contain exactly Length[path] npoint records in the same segment order. The tight-binding curve includes both endpoints of every segment." ->
      "referenceBands 使用 vaspEig 记录格式，必须按同一线段顺序恰好包含 Length[path] npoint 条记录；紧束缚曲线包含每段的两个端点。",
    "rules can be a rule list or Association. After applying rules and 2 Pi times each path momentum, H must be a numerical square matrix." ->
      "rules 可以是规则列表或 Association。应用 rules，并把每个路径动量乘以 2 Pi 后，H 必须成为数值方阵。",
    "The result is a Labeled plot: tight-binding bands are purple and reference bands are blue. plotRange changes only the visible energy interval." ->
      "返回 Labeled 图：紧束缚能带为紫色，参考能带为蓝色；plotRange 只改变可见能量区间。",
    "compareBand returns $Failed and issues a message when H, path, npoint, rules, or referenceBands are malformed, their record counts disagree, or H remains nonnumeric. It does not repair the input or switch to an approximate fallback." ->
      "H、path、npoint、rules 或 referenceBands 格式错误、记录数不一致或 H 仍非数值时，compareBand 给出消息并返回 $Failed；函数不会修补输入，也不会改用近似结果。",
    "All displays the complete overlay; {emin,emax} restricts only the visible energy interval." ->
      "All 显示完整叠加图；{emin,emax} 只限制可见能量区间。",
    "Read the VASP bands and construct the symmetry-generated model:" ->
      "读取 VASP 能带，并用对称性生成模型：",
    "Fit the nearest-neighbour model, add second-neighbour hopping, and refit all 90 records:" ->
      "先拟合最近邻模型，再加入第二近邻 hopping，并重新拟合全部 90 条记录：",
    "Display all three selected VASP bands together with the fitted second-neighbour model from -1.5 to 3 eV:" ->
      "在 -1.5 至 3 eV 范围内同时显示选取的三条 VASP 能带和拟合后的第二近邻模型：",
    "plotRange - zoom into the interval from -0.5 to 0.5 after inspecting all three bands:" ->
      "plotRange - 查看完整三条能带后，放大到 -0.5 至 0.5 的能量区间："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "showband"
  }] -> <|
    "draws a static band plot along standardKPath[] for the current initialized lattice." ->
      "沿当前已初始化晶格的 standardKPath[] 绘制静态能带图。",
    "The verified three-argument form returns the same Graphics form and shares all display options with bandplot." ->
      "已经验证的三参数形式返回与 bandplot 相同形式的 Graphics，并共享全部显示选项。",
    "A current initialized model is required because the path is selected from the lattice metric. Use bandplot[path,npoint,H,rules] when an explicit path is required." ->
      "路径由当前晶格度量自动选择，因此必须先初始化模型；需要显式路径时请使用 bandplot[path,npoint,H,rules]。",
    StringJoin[
      "Usage.wl also advertises showband[path,npoint,H,rules], but the current testing implementation ",
      "misroutes that exact four-argument form through the automatic-path overload. This page does not ",
      "present the unresolved form as working behavior."
    ] -> StringJoin[
      "Usage.wl 还列出了 showband[path,npoint,H,rules]，但当前 testing 实现会把这个四参数形式",
      "错误分派到自动路径重载；在该问题解决前，本页不把它写成可用调用。"
    ],
    "showband returns $Failed and issues a message when automatic path recognition fails or any shared bandplot input is invalid. It does not repair the input or switch to an approximate fallback." ->
      "自动路径识别失败或任一共享 bandplot 输入非法时，showband 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Initialize the one-orbital simple-cubic model used by the automatic path:" ->
      "用 gray[221]、init 与 symham 初始化自动路径示例所用的单轨道简单立方模型：",
    "Let showband choose the conventional simple-cubic path and draw the dispersion:" ->
      "让 showband 自动选择简单立方约定路径并绘制能带：",
    "Vertical energy range passed to the band Graphics." -> "传给能带 Graphics 的纵向能量范围。",
    "Automatic, None, positions, or explicit vertical tick specifications." -> "可取 Automatic、None、位置列表或显式纵向刻度设置。",
    "Positive size for path labels, frame ticks, and frame labels." -> "路径标签、边框刻度与边框标签使用的正字号。",
    "Font family used for path labels, frame ticks, and frame labels." -> "路径标签、边框刻度与边框标签使用的字体族。",
    "Automatic or the final Graphics image size." -> "Automatic 或最终 Graphics 图像尺寸。",
    "restricts only the displayed energy interval; the sampled eigenvalues are unchanged." -> "只限制显示的能量区间，不改变已经采样的本征值。",
    "accepts explicit positions and labels for the vertical frame axis." -> "为纵向边框轴设置显式位置与标签。",
    "sets the path labels and tick-label size." -> "设置路径标签与刻度标签字号。",
    "sets one font family for path labels and frame ticks." -> "为路径标签与边框刻度设置字体族。",
    "sets the final Graphics size without changing the band data." -> "设置最终 Graphics 尺寸而不改变能带数据。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "standardKPath"
  }] -> <|
    "returns a conventional high-symmetry path for the current lattice in reciprocal fractional coordinates." ->
      "以倒空间分数坐标返回当前晶格的约定高对称路径。",
    "The algorithm recognizes the lattice metric and then reads a built-in path table covering 25 Brillouin-zone variations of the 14 three-dimensional Bravais lattice types." ->
      "算法先识别晶格度量，再读取覆盖 14 类三维 Bravais 晶格、25 种 Brillouin 区变体的内置路径表。",
    StringJoin[
      "This is not a space-group little-group or irreducible-representation calculation. Nonstandard ",
      "primitive bases and centred-lattice conventions can prevent metric recognition; use an explicitly ",
      "standardized basis or the matching BravaisType value."
    ] -> StringJoin[
      "这不是空间群小群或不可约表示计算。非标准原胞基矢与中心化晶格约定可能使度量识别失败；",
      "应使用显式标准化的基矢，或给出匹配的 BravaisType。"
    ],
    "The returned path uses Bradley-Cracknell display labels where the correspondence is unambiguous, with Setyawan-Curtarolo labels retained as a fallback." ->
      "当对应关系明确时，返回路径使用 Bradley-Cracknell 显示标签；否则保留 Setyawan-Curtarolo 标签。",
    "standardKPath returns $Failed and issues a message when no current model exists, the evaluated lattice is singular or nonnumeric, or the requested BravaisType conflicts with the recognized metric. It does not repair the input or switch to an approximate fallback." ->
      "当前模型不存在、数值晶格奇异或非数值，或者指定的 BravaisType 与识别出的度量冲突时，standardKPath 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Initialize a one-orbital simple-cubic model before requesting its conventional path:" ->
      "请求约定路径前，先用 gray[221] 初始化单轨道简单立方模型：",
    "For a simple-cubic lattice, the path contains the conventional Gamma-X-M-Gamma-R-X and M-R segments:" ->
      "简单立方晶格的路径包含约定的 Gamma-X-M-Gamma-R-X 与 M-R 线段：",
    "Automatic metric recognition, or the matching recognized family/BZ string such as \"CubicPrimitive\" or \"CUB\"." ->
      "Automatic 自动识别度量；也可给出匹配的晶格族或 Brillouin 区字符串，例如 \"CubicPrimitive\" 或 \"CUB\"。",
    "Positive numerical tolerance used only in metric classification." -> "仅在度量分类中使用的正数容差。",
    "makes the recognized convention explicit without changing the cubic path." ->
      "显式写出已经识别的晶格约定，不改变立方路径。",
    "controls only metric classification and leaves this exact cubic lattice unchanged." ->
      "只控制度量分类，不改变这个精确立方晶格。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "surfaceGreenFunction"
  }] -> <|
    "computes the retarded Green function of a semi-infinite surface by principal-layer decimation." ->
      "通过主层消元计算半无限表面的推迟 Green 函数。",
    "Cell axes 1 and 2 span the surface; axis 3 is the stacking direction. Longer-range hopping along axis 3 is grouped into an exact principal layer." ->
      "晶胞轴 1、2 张成表面，轴 3 是堆垛方向；沿轴 3 的长程 hopping 会被精确归并到主层中。",
    "The default Positive surface is terminated toward increasing layer index. Output -> \"Data\" includes the H00/H01/H10 blocks, iterations, residuals, and spectral weight." ->
      "默认 Positive 表面终止于层编号增大的方向。Output -> \"Data\" 返回 H00/H01/H10 块、迭代次数、残差与谱权重。",
    "surfaceGreenFunction returns $Failed and issues a message when data, surface momentum, real energy, cell transform, decimation controls, surface side, Hermiticity, or output form is invalid; nonconvergence is reported rather than hidden. It does not repair the input or switch to an approximate fallback." ->
      "数据、表面动量、实能量、晶胞变换、消元控制、表面方向、厄米性或输出形式非法时，surfaceGreenFunction 发出消息并返回 $Failed；不收敛会明确报告，不会被隐藏，也不会切换到近似 fallback。",
    "Build the symmetry-generated four-band C4 topological crystalline insulator (TCI):" ->
      "用 gray[75]、init、symham 与 hoppingData 构造四带 C4 拓扑晶体绝缘体（TCI）：",
    "Compute the retarded surface Green matrix at the surface Gamma point:" ->
      "计算表面 Gamma 点的推迟表面 Green 矩阵：",
    "Plot the surface spectral function of the same C4 TCI:" ->
      "绘制同一 C4 TCI 的表面谱函数：",
    "Automatic keeps the input axes; a nonsingular integer matrix defines the surface plane and stacking direction first." ->
      "Automatic 保留输入晶胞轴；非奇异整数矩阵先定义表面平面与堆垛方向。",
    "Positive imaginary retarded-energy broadening eta." -> "推迟能量的正虚部展宽 eta。",
    "Positive relative coupling tolerance for principal-layer decimation." -> "主层消元的正相对耦合容差。",
    "Positive integer decimation iteration limit." -> "正整数消元迭代次数上限。",
    "Positive or Negative termination along the stacking axis." -> "沿堆垛轴选择 Positive 或 Negative 终止面。",
    "Nonnegative residual tolerance for hopping and principal-layer blocks." -> "hopping 与主层块所用的非负残差容差。",
    "GreenFunction, SpectralWeight, or a diagnostic Data Association." -> "返回 GreenFunction、SpectralWeight 或诊断 Data Association。",
    "uses an explicit exact calculation cell before constructing the principal layers." -> "构造主层前使用显式的精确计算晶胞。",
    "changes the positive retarded imaginary part of energy." -> "改变推迟能量的正虚部。",
    "sets the decimation coupling convergence threshold." -> "设置消元耦合的收敛阈值。",
    "sets an explicit upper bound while leaving sufficient iterations for this chain." -> "显式设置迭代上限，并为该模型保留足够迭代次数。",
    "selects the negative termination of the topological model." -> "选择拓扑模型的 Negative 终止面。",
    "uses a stricter residual threshold for exact hopping blocks." -> "对精确 hopping 块采用更严格的残差阈值。",
    "returns the scalar surface spectral weight directly." -> "直接返回标量表面谱权重。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "surfaceSpectralFunction"
  }] -> <|
    "returns -Im Tr[Gsurface]/Pi for the semi-infinite surface." ->
      "返回半无限表面的 -Im Tr[Gsurface]/Pi。",
    "The geometry, principal-layer construction, and convergence rules are identical to surfaceGreenFunction." ->
      "几何约定、主层构造与收敛规则和 surfaceGreenFunction 完全相同。",
    "The return value is a real spectral weight, not the Green matrix. Use surfaceGreenFunction with Output -> \"Data\" for block and convergence diagnostics." ->
      "返回值是实数谱权重，不是 Green 矩阵；需要矩阵块与收敛诊断时，请使用 surfaceGreenFunction 的 Output -> \"Data\"。",
    "surfaceSpectralFunction returns $Failed and issues a message when any shared surface Green-function input or convergence check fails. It does not repair the input or switch to an approximate fallback." ->
      "任一共享表面 Green 函数输入或收敛检查失败时，surfaceSpectralFunction 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Build the symmetry-generated four-band C4 topological crystalline insulator (TCI):" ->
      "用 gray[75]、init、symham 与 hoppingData 构造四带 C4 拓扑晶体绝缘体（TCI）：",
    "Plot the semi-infinite surface spectral weight of the C4 TCI:" ->
      "绘制该 C4 TCI 的半无限表面谱权重：",
    "Optional exact integer cell defining surface and stacking axes." -> "定义表面轴与堆垛轴的可选精确整数晶胞。",
    "Positive retarded broadening eta." -> "正的推迟展宽 eta。",
    "Positive decimation convergence tolerance." -> "正的消元收敛容差。",
    "Positive decimation iteration limit." -> "正的消元迭代次数上限。",
    "Positive or Negative termination along cell axis 3." -> "沿晶胞轴 3 选择 Positive 或 Negative 终止面。",
    "Nonnegative residual tolerance for hopping and principal-layer blocks." -> "hopping 与主层块所用的非负残差容差。",
    "uses an explicit exact calculation cell." -> "使用显式的精确计算晶胞。",
    "broadens the surface pole by eta=0.01." -> "用 eta=0.01 展宽表面极点。",
    "sets the principal-layer coupling convergence threshold." -> "设置主层耦合的收敛阈值。",
    "sets an explicit safe iteration bound for this chain." -> "为该模型设置显式且足够的迭代上限。",
    "selects the negative termination." -> "选择 Negative 终止面。",
    "uses a stricter residual threshold for exact hopping blocks." -> "对精确 hopping 块采用更严格的残差阈值。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "symmetrizationHRInit"
  }] -> <|
    "prepares the orbital centers and symmetry matrices needed to symmetrize a Wannier90 HR model." ->
      "准备对 Wannier90 HR 模型做对称化所需的轨道中心与对称矩阵。",
    StringJoin[
      "To symmetrize a Wannier90 HR model, one needs the center of every Wannier orbital and the matrix of ",
      "every symmetry operation in the same orbital order. This function constructs those data from the ",
      "lattice, Wyckoff positions, symmetry operations, and local orbitals."
    ] -> StringJoin[
      "对 Wannier90 HR 模型做对称化，需要按同一轨道顺序给出每个 Wannier 轨道中心和每个对称操作矩阵。",
      "本函数根据晶格、Wyckoff 位置、对称操作与局域轨道构造这些数据。"
    ],
    "The orbital order must agree with the external electronic-structure calculation. At present the function supports the documented VASP ordering." ->
      "轨道顺序必须与外部电子结构计算一致；当前仅支持文档规定的 VASP 顺序。",
    "The result contains wcc for the fractional Wannier centers, DR for the full symmetry matrices, and symmetry for the corresponding ordered operations." ->
      "结果包含分数 Wannier 中心 wcc、完整对称矩阵 DR，以及对应的有序操作 symmetry。",
    "Symbolic direct-lattice vectors stored by rows." -> "按行存放的符号正格矢。",
    "Numerical replacement rules for the symbolic lattice." -> "符号晶格的数值替换规则。",
    "Wyckoff seeds and magnetic moments in the same nesting used by init." -> "与 init 使用相同嵌套形式的 Wyckoff 种子与磁矩。",
    "Complete ordered symmetry records." -> "完整且有序的对称操作记录。",
    "Local basis used to build Wannier centers and point matrices." -> "构造 Wannier 中心与点操作矩阵所用的局域基。",
    "Orbital ordering convention; VASP is the only supported value." -> "轨道排序约定；当前唯一支持的值是 VASP。",
    "Prepare a full tetragonal C4 model with the complete gray-group 75 operation list and a px/py basis, then show a compact summary instead of a large raw Association:" ->
      "使用灰群 75 的完整操作列表与 px/py 基准备四方 C4 模型，并以紧凑 Grid 显示摘要，而不输出庞大的原始 Association：",
    "\"Lattice\" -> DiagonalMatrix[{2 a,2 a,3 c}] changes the tetragonal direct-lattice rows while preserving C4 symmetry. All eight operation matrices remain 2 by 2:" ->
      "\"Lattice\" -> DiagonalMatrix[{2 a,2 a,3 c}] 在保持 C4 对称性的同时改变四方正格矢；八个操作矩阵仍均为 2×2：",
    "\"LattPar\" -> {a->2,c->3} evaluates every symbolic tetragonal lattice entry before the symmetry matrices are constructed:" ->
      "\"LattPar\" -> {a->2,c->3} 在构造对称矩阵前数值化四方晶格中的全部符号：",
    "\"WyckoffPosition\" -> {{{1/2,1/2,0},{0,0,0}}} places both C4-related local orbitals on the invariant fractional center {1/2,1/2,0}:" ->
      "\"WyckoffPosition\" -> {{{1/2,1/2,0},{0,0,0}}} 把两个 C4 相关局域轨道放在不变分数中心 {1/2,1/2,0}：",
    "\"Symminformation\" may select the complete unitary C4 subgroup from gray group 75. The returned records preserve all four ordered operations:" ->
      "\"Symminformation\" 可选取灰群 75 的完整幺正 C4 子群；返回记录保留全部四个有序操作：",
    "\"BasisFunctions\" -> {\"s\",\"pz\"} changes the VASP-ordered local basis while retaining the same complete C4 symmetry list:" ->
      "\"BasisFunctions\" -> {\"s\",\"pz\"} 改变按 VASP 排序的局域基，同时保留同一完整 C4 对称列表：",
    "\"Software\" -> \"VASP\" selects the currently supported external orbital-order convention. The result has exactly the documented three keys:" ->
      "\"Software\" -> \"VASP\" 选择当前支持的外部轨道顺序约定；结果恰好含文档列出的三个键：",
    "Only the documented VASP orbital ordering is currently available:" -> "当前仅支持文档规定的 VASP 轨道顺序："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "transformHoppings"
  }] -> <|
    "rewrites complete numerical hopping blocks in the exact integer cell whose direct-lattice rows are T times the old rows." ->
      "把完整数值 hopping 块精确改写到整数晶胞中；新正格矢行等于 T 乘旧正格矢行。",
    "T must be a nonsingular 3 by 3 integer matrix. Exact Smith decomposition enumerates all cell representatives; no real-space cutoff is introduced." ->
      "T 必须是非奇异 3×3 整数矩阵。算法用精确 Smith 分解枚举全部晶胞代表，不引入实空间截断。",
    "The returned Association has Abs[Det[T]] times the original orbital dimension and, when available, transformed lattice vectors and Wannier centers." ->
      "返回 Association 的轨道维数是原来的 Abs[Det[T]] 倍；若输入提供元数据，还返回变换后的晶格矢量与 Wannier 中心。",
    "transformHoppings returns $Failed and issues a message when the hopping Association/path is incomplete, T is not a nonsingular integer matrix, center/lattice metadata are malformed, or the Hermiticity residual exceeds the requested tolerance. It does not repair the input or switch to an approximate fallback." ->
      "hopping Association 或路径不完整、T 不是非奇异整数矩阵、中心或晶格元数据非法，或者厄米残差超过容差时，transformHoppings 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Build the tetragonal gray-group 75 hopping model used below:" ->
      "用 gray[75]、init、symham 与 hoppingData 构造以下示例共用的四方 hopping 模型：",
    "Double the chain cell along axis 3. The transformed model contains two orbitals per new cell:" ->
      "沿轴 3 把晶胞加倍；变换后的每个新晶胞含两个轨道：",
    "Automatic uses centers stored in data; None omits them; an explicit list supplies one fractional three-vector per input orbital." ->
      "Automatic 使用 data 中保存的中心；None 省略中心；显式列表为每个输入轨道给出一个三分量分数中心。",
    "Automatic uses the stored row-vector lattice; None omits it; an explicit real numerical 3 by 3 matrix replaces it." ->
      "Automatic 使用保存的行矢量晶格；None 省略晶格；显式实数 3×3 矩阵替换它。",
    "Nonnegative numerical tolerance for H(R)=H(-R)^dagger before the exact cell transform." ->
      "精确晶胞变换前检查 H(R)=H(-R)^dagger 的非负数值容差。",
    "uses explicit orbital centers instead of metadata stored in the input." -> "使用显式轨道中心，而不是输入中保存的元数据。",
    "uses an explicit old direct lattice and returns the transformed row vectors." -> "使用显式旧正格子，并返回变换后的行矢量。",
    "checks the input hopping table against the stated numerical tolerance." -> "按给定数值容差检查输入 hopping 表。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "unsymham"
  }] -> <|
    "returns the translation-periodic Hamiltonian of shell n before point-group, magnetic, or Hermitian constraints are imposed." ->
      "返回近邻 hopping 层 n 在施加点群、磁对称或厄米约束之前的平移周期 Hamiltonian。",
    "Every directed hopping matrix element is independent, although the bond geometry and Bloch phases still come from the lattice prepared by init." ->
      "每个有向 hopping 矩阵元都相互独立，但键几何与 Bloch 相位仍来自 init 准备的晶格。",
    "This is mainly useful for understanding how symmetry reduces the parameter space. Use symham to construct the physical symmetry-constrained model." ->
      "该函数主要用于理解对称性如何缩减参数空间；构造物理的对称约束模型应使用 symham。",
    "Initialize the complete cubic gray-group 221 s-orbital model, then display the independent complex amplitude carried by every directed nearest-neighbour bond:" ->
      "用完整的立方灰群 221 初始化 s 轨道模型，再显示每条有向最近邻键携带的独立复振幅：",
    "For the same prepared shell, symham reduces these independent entries to the symmetry-allowed Hermitian subspace:" ->
      "对同一已准备的近邻 hopping 层，symham 把这些独立矩阵元缩减到对称性允许的厄米子空间：",
    "The shell must already have been prepared by init; the function never extends the bond search:" ->
      "近邻 hopping 层必须已经由 init 准备；该函数不会自行扩展键搜索："
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "wilsonLoop"
  }] -> <|
    "computes Wilson-loop eigenphases for an isolated occupied subspace along a forward reciprocal-lattice loop." ->
      "沿正向倒格矢闭环计算孤立占据子空间的 Wilson-loop 本征相位。",
    StringJoin[
      "start and end have the same coordinate dimension and end-start must be 2 Pi times an integer ",
      "reciprocal vector. wcc contains one fractional orbital center per Hamiltonian basis state."
    ] -> "start 与 end 的坐标维数必须相同，且 end-start 必须是 2 Pi 乘整数倒格矢；wcc 为每个 Hamiltonian 基态给出一个分数轨道中心。",
    StringJoin[
      "Occupied-subspace overlaps are unitarized by SVD and the endpoint is closed by the orbital-center ",
      "sewing matrix. The production calculation does not use an eigenvector-gauge fallback."
    ] -> "占据子空间重叠通过 SVD 幺正化，端点由轨道中心 sewing 矩阵闭合；生产计算不使用本征矢规范 fallback。",
    "The default return is a sorted list of phases divided by Pi. Alternate outputs expose Wannier centers modulo one, Wilson eigenvalues, or a diagnostic Association." ->
      "默认返回按序排列并除以 Pi 的相位列表；其他输出可返回模 1 Wannier 中心、Wilson 本征值或诊断 Association。",
    "wilsonLoop returns $Failed and issues a message when the loop is not reciprocal-lattice closed, centers/occupation are inconsistent, H is non-Hermitian, the direct gap closes, endpoint covariance fails, an overlap is singular, or an option is invalid. It does not repair the input or switch to an approximate fallback." ->
      "闭环不满足倒格矢闭合、中心或占据数不一致、H 非厄米、直接能隙闭合、端点协变失败、重叠奇异或选项非法时，wilsonLoop 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Build a gapped Chern slice of the BNS 143.3 symmetry-generated Weyl model:" ->
      "用 BNS 143.3 对称性生成的 Weyl 模型构造有隙 Chern 切片：",
    "Compute one occupied-band Wilson loop on the Chern slice:" -> "在该 Chern 切片上计算单占据带 Wilson 环：",
    "Plot the Wilson phase while the loop is translated across the transverse Brillouin zone:" ->
      "当闭环横向平移穿过 Brillouin 区时，绘制 Wilson 相位谱流：",
    "Positive number of equal forward-path intervals." -> "正整数等分的前向路径区间数。",
    "Nonnegative Hermiticity and final unitarity tolerance." -> "厄米性与最终幺正性的非负容差。",
    "Nonnegative minimum direct-gap threshold for an occupied subspace smaller than the full Hilbert space." ->
      "占据子空间小于完整 Hilbert 空间时使用的非负最小直接能隙阈值。",
    "Nonnegative endpoint sewing-covariance tolerance." -> "端点 sewing 协变性的非负容差。",
    "Nonnegative lower bound for link singular values." -> "链接奇异值的非负下界。",
    "PhaseOverPi, WannierCenters, Eigenvalues, or Data." -> "返回 PhaseOverPi、WannierCenters、Eigenvalues 或 Data。",
    "uses twelve equal forward links on the Chern slice." -> "在 Chern 切片上使用 12 条等长前向链接。",
    "uses a stricter Hermiticity and final unitarity threshold." -> "使用更严格的厄米性与最终幺正性阈值。",
    "sets the occupied-band direct-gap threshold." -> "设置占据带直接能隙阈值。",
    "sets the endpoint sewing-covariance threshold." -> "设置端点 sewing 协变阈值。",
    "sets the minimum acceptable occupied-link singular value." -> "设置可接受占据链接奇异值的下限。",
    "returns the diagnostic Wilson matrix, residuals, gap, path, and closure vector." ->
      "返回诊断 Wilson 矩阵、残差、能隙、路径与闭合矢量。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "plotRealSpaceWavefunction"
  }] -> <|
    "data must be buildRealSpaceHamiltonian[...,\"Output\"->\"Data\"] with PositionDataAvailable -> True, complete BasisRecords, a numerical Lattice, and WannierCenters metadata." ->
      "data 必须是 buildRealSpaceHamiltonian[...,\"Output\"->\"Data\"] 的返回值，并且 PositionDataAvailable -> True，含完整 BasisRecords、数值 Lattice 与 WannierCenters 元数据。",
    "state is a finite numerical vector whose length equals data[\"Dimension\"]. Sphere volume represents probability. Atom aggregation combines all basis states at the same Cartesian position; Orbital aggregation preserves the matrix basis order." ->
      "state 是有限数值向量，长度必须等于 data[\"Dimension\"]。球体体积表示概率；Atom 聚合合并同一笛卡尔位置的全部基态，Orbital 聚合保留矩阵基的顺序。",
    "PhaseColoring uses the phase of the largest-amplitude basis state in each aggregate. Weight remains the sum of squared amplitudes in that aggregate." ->
      "PhaseColoring 使用每个聚合记录中振幅最大基态的相位；Weight 始终是该记录内各振幅模平方之和。",
    "The Graphics3D result always uses Axes -> False, Boxed -> False, and an Orthographic parallel projection. Output -> \"Data\" returns normalized amplitudes, SiteRecords, visible count, maximum weight, and Graphics." ->
      "Graphics3D 结果固定采用 Axes -> False、Boxed -> False 与 Orthographic 平行投影。Output -> \"Data\" 返回归一化振幅、SiteRecords、可见记录数、最大权重和 Graphics。",
    StringJoin[
      "The second example constructs a two-dimensional C4 breathing-square ",
      "higher-order topological insulator directly with gray group 75 and symham. ",
      "With |s1/t1| = 0.2 < 1, it projects one corner ",
      "orbital onto the four lowest-absolute-energy finite-system states to ",
      "select a reproducible representative corner mode. The plot visualizes ",
      "localization; it does not by itself calculate the quadrupole invariant."
    ] -> StringJoin[
      "第二个例子用灰群 75 和 symham 直接构造二维 C4 呼吸方格高阶拓扑绝缘体。",
      "当 |s1/t1| = 0.2 < 1 时，把一个角轨道投影到有限体系中绝对能量最低的四态子空间，",
      "从而选出可重现的代表性角模。该图显示角局域性，但本身不计算四极矩拓扑不变量。"
    ],
    "plotRealSpaceWavefunction returns $Failed and issues a message when the finite-system data lack position metadata, the state has the wrong length or zero norm, an aggregation/display option is invalid, or the color function does not return colors. It does not repair the input, reduce the grid silently, or switch to an approximate fallback." ->
      "当有限体系数据缺少位置元数据、态向量长度错误或范数为零、聚合或显示选项非法，或者颜色函数不返回颜色时，plotRealSpaceWavefunction 发出消息并返回 $Failed；它不会修补输入、静默缩小网格或切换到近似 fallback。",
    "Plot the probability and complex phase of a wavefunction on a three-site L-shaped finite system:" ->
      "绘制三格点 L 形有限体系中波函数的概率与复相位：",
    "Plot a representative corner-localized state of the symmetry-generated two-dimensional C4 higher-order topological insulator:" ->
      "绘制由对称性生成的二维 C4 高阶拓扑绝缘体的代表性角局域态：",
    "Atom combines basis states at the same Cartesian position; Orbital keeps every matrix basis state separate." ->
      "Atom 合并同一笛卡尔位置的基态；Orbital 保留每个矩阵基态。",
    "Whether the input state is normalized before probabilities are calculated." ->
      "计算概率前是否归一化输入态。",
    "Nonnegative probability threshold for visible wavefunction spheres." ->
      "显示波函数球体所用的非负概率阈值。",
    "Positive Cartesian-distance tolerance used only for Atom aggregation." ->
      "仅供 Atom 聚合使用的正笛卡尔距离容差。",
    "Whether visible spheres are colored by representative complex phase instead of probability." ->
      "是否按代表性复相位而不是概率给可见球体着色。",
    "Whether faint background spheres are drawn for every site record." ->
      "是否为每条格点记录绘制浅色背景球。",
    "Whether site indices are drawn beside every site record." ->
      "是否在每条格点记录旁绘制格点编号。",
    "Automatic uses 0.055 times the shortest lattice vector; a positive value sets background-site radius." ->
      "Automatic 取最短晶格矢量长度的 0.055 倍；正数显式设置背景格点半径。",
    "Automatic uses 0.22 times the shortest lattice vector; a positive value sets maximum wavefunction-sphere radius." ->
      "Automatic 取最短晶格矢量长度的 0.22 倍；正数显式设置波函数球体最大半径。",
    "Graphics returns the real-space plot; Data also returns amplitudes, site records, weights, and Graphics." ->
      "Graphics 返回实空间图；Data 还返回振幅、格点记录、权重与 Graphics。",
    "Automatic uses SolarColors; a named ColorData scheme or color function colors probability when PhaseColoring is False." ->
      "Automatic 使用 SolarColors；PhaseColoring 为 False 时，也可用命名 ColorData 色表或颜色函数表示概率。",
    "Positive site-label font size." -> "正的格点标签字号。",
    "Graphics3D viewpoint; projection remains orthographic and axes/box remain hidden." ->
      "Graphics3D 视点；投影保持正交，坐标轴与盒框保持隐藏。",
    "Final Graphics3D image size." -> "最终 Graphics3D 图像尺寸。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "plotSurfaceSpectrum"
  }] -> <|
    "data is the numerical hopping Association accepted by surfaceGreenFunction. Cell axes 1 and 2 span the surface, and axis 3 is the semi-infinite stacking direction." ->
      "data 是 surfaceGreenFunction 接受的数值 hopping Association。晶胞轴 1、2 张成表面，轴 3 是半无限堆垛方向。",
    "path may be explicit surface-momentum points or connected labeled two-dimensional segments. The plotted value is -Im Tr[Gsurface]/Pi on a uniform energy grid." ->
      "path 可以是显式表面动量点，也可以是相连并带标签的二维线段。图中数值是在均匀能量网格上计算的 -Im Tr[Gsurface]/Pi。",
    "Every grid point uses the same principal-layer decimation and convergence checks as surfaceSpectralFunction; nonconvergence is reported instead of hidden." ->
      "每个网格点都使用与 surfaceSpectralFunction 相同的主层消元与收敛检查；不收敛会明确报告，不会被隐藏。",
    StringJoin[
      "For the P4 Chern example, transformHoppings moves the original second in-plane direction to cell ",
      "axis 3; axis 1 is then the conserved edge momentum and axis 3 is semi-infinite. The gray-P4 ",
      "four-band example uses its stored cell with axis 3 as the stacking direction. The resulting plots ",
      "show one chiral gap-crossing branch and a counterpropagating surface-state pair, respectively."
    ] -> StringJoin[
      "在 P4 Chern 例中，transformHoppings 把原来的第二个面内方向移到晶胞轴 3；",
      "于是轴 1 是守恒边缘动量，轴 3 取半无限。灰群 P4 四带例直接使用其保存的晶胞，",
      "同样令轴 3 为堆垛方向。两幅图分别显示一条穿越能隙的手征分支和一对反向传播的表面态分支。"
    ],
    "Output -> \"Data\" returns SurfaceMomentumPath, PathCoordinate, boundary labels, Energies, SpectralWeight, Broadening, Surface, and Graphics." ->
      "Output -> \"Data\" 返回 SurfaceMomentumPath、PathCoordinate、边界标签、Energies、SpectralWeight、Broadening、Surface 与 Graphics。",
    "plotSurfaceSpectrum returns $Failed and issues a message when the hopping data, path, ordered energy range, surface options, decimation controls, or display options are invalid, or a grid-point calculation does not converge. It does not repair the input, reduce the grid silently, or switch to an approximate fallback." ->
      "当 hopping 数据、路径、有序能量范围、表面选项、消元控制或显示选项非法，或者某个网格点不收敛时，plotSurfaceSpectrum 发出消息并返回 $Failed；它不会修补输入、静默缩小网格或切换到近似 fallback。",
    "The symmetry-generated P4 Chern insulator has a single chiral edge branch that crosses the bulk gap:" ->
      "由对称性生成的 P4 Chern 绝缘体有一条穿越体能隙的手征边缘分支：",
    "The symmetry-generated gray-P4 Z2 model has a pair of helical surface branches crossing at zero momentum:" ->
      "由对称性生成的灰群 P4 Z2 模型有一对在零动量处相交的螺旋表面分支：",
    "Positive number of intervals per labeled surface-momentum segment." -> "每段带标签表面动量线段中的正整数区间数。",
    "Number of uniformly spaced energies, at least two." -> "均匀能量采样点数，至少为 2。",
    "Automatic keeps the hopping cell; a nonsingular integer 3 by 3 matrix defines the surface cell first." -> "Automatic 保留 hopping 晶胞；非奇异 3×3 整数矩阵先定义表面晶胞。",
    "Positive imaginary retarded-energy broadening." -> "推迟能量的正虚部展宽。",
    "Positive relative coupling tolerance for principal-layer decimation." -> "主层消元的正相对耦合容差。",
    "Positive decimation iteration limit." -> "正的消元迭代次数上限。",
    "Positive or Negative termination along cell axis 3." -> "晶胞轴 3 正向或负向终止面。",
    "Nonnegative tolerance for hopping and principal-layer Hermiticity." -> "hopping 与主层厄米性的非负容差。",
    "Graphics returns the density plot; Data also returns the sampled path, energies, spectral weights, and Graphics." -> "Graphics 返回密度图；Data 还返回采样路径、能量、谱权重与 Graphics。",
    "Color scheme or color function passed to ListDensityPlot." -> "传给 ListDensityPlot 的配色方案或颜色函数。",
    "Whether ListDensityPlot rescales values before applying ColorFunction." -> "ListDensityPlot 在应用 ColorFunction 前是否缩放数值。",
    "Legend specification passed to ListDensityPlot." -> "传给 ListDensityPlot 的图例设置。",
    "Displayed plot range passed to ListDensityPlot." -> "传给 ListDensityPlot 的显示范围。",
    "Zero or one for the sampled spectral grid." -> "谱采样网格使用 0 阶或 1 阶插值。",
    "Mesh specification passed to ListDensityPlot." -> "传给 ListDensityPlot 的网格线设置。",
    "Two frame labels; Automatic uses surface k path and energy." -> "两个边框标签；Automatic 使用表面 k 路径与能量。",
    "Positive size for frame ticks and labels." -> "边框刻度与标签的正字号。",
    "Font family for frame ticks and labels." -> "边框刻度与标签的字体族。",
    "Height-to-width ratio of the density plot." -> "密度图的高宽比。",
    "Final Graphics image size." -> "最终 Graphics 图像尺寸。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "plotWilsonLoop"
  }] -> <|
    StringJoin[
      "start and end have the same momentum-coordinate dimension, and end-start ",
      "must be 2 Pi times an integer reciprocal vector. centers contains one ",
      "fractional orbital center per Hamiltonian basis state."
    ] -> "start 与 end 的动量坐标维数必须相同，end-start 必须是 2 Pi 乘以整数倒格矢；centers 为每个 Hamiltonian 基态给出一个分数轨道中心。",
    StringJoin[
      "parameterPath may be an explicit ordered list of distinct momentum offsets ",
      "or connected labeled segments {{{p1,p2},{label1,label2}},...}. Labeled ",
      "segments use ParameterSubdivisions intervals per segment."
    ] -> "parameterPath 可以是互不相同的动量偏移有序列表，也可以是相连的带标签线段 {{{p1,p2},{label1,label2}},...}；每段使用 ParameterSubdivisions 个区间。",
    StringJoin[
      "Every sample calls wilsonLoop with forward occupied-subspace links, ",
      "SVD unitarization, and orbital-center endpoint sewing. The plot introduces ",
      "no gauge, gap, or numerical fallback."
    ] -> "每个采样点都调用 wilsonLoop，使用占据子空间正向链接、SVD 幺正化和轨道中心端点缝合；绘图过程不引入规范、能隙或数值 fallback。",
    StringJoin[
      "The symmetry-generated P4 Chern example has one occupied band whose Wilson phase winds once over a full ",
      "transverse Brillouin zone. The gray-P4 Z2 example has two occupied bands; their Wannier-center ",
      "branches exchange partners between the two time-reversal-invariant endpoints of half the zone. ",
      "plotWilsonLoop displays this flow but does not return a Z2 integer."
    ] -> "由对称性生成的 P4 Chern 例有一个占据带，其 Wilson 相位在完整横向 Brillouin 区绕行一次。灰群 P4 Z2 例有两个占据带，其 Wannier 中心分支在半区两个时间反演不变端点间交换伙伴。plotWilsonLoop 显示该谱流，但不返回 Z2 整数。",
    "The default return is Graphics. Output -> \"Data\" returns ParameterPath, PathCoordinate, boundary labels, PhaseConvention, Values, BranchCount, and the same Graphics." ->
      "默认返回 Graphics。Output -> \"Data\" 返回 ParameterPath、PathCoordinate、边界标签、PhaseConvention、Values、BranchCount 与同一 Graphics。",
    "plotWilsonLoop returns $Failed and issues a message when the loop or parameter path is malformed, a translated Wilson loop fails its Hermiticity/gap/covariance/overlap checks, or an option is invalid. It does not repair the input, reduce the grid silently, or switch to an approximate fallback." ->
      "当闭环或参数路径格式错误、平移后的 Wilson 环未通过厄米性、能隙、协变性或重叠检查，或者选项非法时，plotWilsonLoop 发出消息并返回 $Failed；它不会修补输入、静默缩小网格或切换到近似 fallback。",
    "The occupied band of the symmetry-generated P4 Chern insulator has one Wilson branch that winds once across the transverse Brillouin zone:" ->
      "由对称性生成的 P4 Chern 绝缘体占据带有一条 Wilson 分支，在横向 Brillouin 区绕行一次：",
    "For the symmetry-generated gray-P4 Z2 insulator, the two occupied Wannier-center branches exchange partners between 0 and Pi:" ->
      "对于由对称性生成的灰群 P4 Z2 绝缘体，两条占据 Wannier 中心分支在 0 与 Pi 之间交换伙伴：",
    "Positive number of intervals per labeled parameter-path segment." -> "每段带标签参数路径中的正整数区间数。",
    "Positive number of forward links in every translated Wilson loop." -> "每个平移 Wilson 环中的正整数正向链接数。",
    "Nonnegative Hamiltonian Hermiticity and Wilson-unitarity tolerance." -> "Hamiltonian 厄米性与 Wilson 幺正性的非负容差。",
    "Nonnegative minimum direct gap for a partially occupied subspace." -> "部分占据子空间的非负最小直接能隙。",
    "Nonnegative endpoint sewing-covariance tolerance." -> "端点缝合协变性的非负容差。",
    "Nonnegative lower bound for occupied-link singular values." -> "占据链接奇异值的非负下界。",
    "PhaseOverPi gives phases on [-1,1]; WannierCenters gives centers modulo one." -> "PhaseOverPi 给出 [-1,1] 上的相位；WannierCenters 给出模 1 的中心。",
    "Graphics returns the branch plot; Data also returns the sampled path, values, labels, and Graphics." -> "Graphics 返回分支图；Data 还返回采样路径、数值、标签与 Graphics。",
    "Automatic, None, numeric positions, or standard vertical tick specifications." -> "可取 Automatic、None、数值位置或标准纵向刻度设置。",
    "Whether sampled points in each Wilson branch are connected." -> "是否连接每条 Wilson 分支中的采样点。",
    "Marker specification passed to ListPlot." -> "传给 ListPlot 的标记设置。",
    "Style or list of styles passed to ListPlot." -> "传给 ListPlot 的样式或样式列表。",
    "Displayed x/y range; Automatic uses the full path and the selected phase convention." -> "显示的 x/y 范围；Automatic 使用完整路径和所选相位约定。",
    "Grid-line specification; Automatic marks path boundaries and reference phases." -> "网格线设置；Automatic 标出路径边界与参考相位。",
    "Two frame labels; Automatic names the parameter path and phase convention." -> "两个边框标签；Automatic 使用参数路径和相位约定名称。",
    "Positive size for ticks, path labels, and frame labels." -> "刻度、路径标签与边框标签的正字号。",
    "Font family for ticks, path labels, and frame labels." -> "刻度、路径标签与边框标签的字体族。",
    "Final Graphics image size." -> "最终 Graphics 图像尺寸。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "pointChernNumber"
  }] -> <|
    "H may be a numerical function H[k] or a matrix containing the configured momentum symbols. The occupied subspace must remain gapped on the cube surface." ->
      "H 可以是数值函数 H[k]，也可以是含所配置动量符号的矩阵；占据子空间在立方体表面必须始终保持有隙。",
    "By default the center must be gapless within CenterGapTolerance. The raw surface flux is accepted only when it lies within IntegerTolerance of an integer." ->
      "默认要求中心能隙不超过 CenterGapTolerance；只有当原始表面通量对应的 Chern 数与整数之差不超过 IntegerTolerance 时才接受结果。",
    "pointChernNumber returns $Failed and issues a message when the point/radius/occupation/Hamiltonian is invalid, the required center is gapped, the enclosing surface gap closes, an overlap is singular, or flux quantization/options fail. It does not repair the input or switch to an approximate fallback." ->
      "当点、半径、占据数或 Hamiltonian 非法，所要求的中心仍有隙，包围表面发生闭隙，重叠矩阵奇异，或者通量量子化或选项检查失败时，pointChernNumber 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Build the BNS 143.3 symmetry-generated magnetic Weyl model:" ->
      "构造由 BNS 143.3 对称性生成的磁性 Weyl 模型：",
    "Compute the charge enclosed around the symmetry-generated Weyl point:" ->
      "计算由对称性生成的 Weyl 点所包围的电荷：",
    "Plot the Berry-curvature vector field around the same Weyl point:" ->
      "绘制同一 Weyl 点周围的 Berry 曲率矢量场：",
    "Positive number of subdivisions on each cube edge." -> "立方体每条边上的正整数细分数。",
    "Nonnegative Hamiltonian Hermiticity tolerance." -> "Hamiltonian 厄米性的非负容差。",
    "Positive minimum direct gap required on the enclosing cube." -> "包围立方体表面所要求的正最小直接能隙。",
    "Nonnegative maximum center gap when RequireGaplessCenter is True." -> "RequireGaplessCenter 为 True 时中心能隙的非负上限。",
    "Nonnegative occupied-overlap singular-value threshold." -> "占据子空间重叠奇异值的非负阈值。",
    "Nonnegative allowed difference between raw flux/2 Pi and the nearest integer." -> "原始通量除以 2 Pi 后与最近整数之间允许的非负差值。",
    "Whether k0 must be gapless within CenterGapTolerance." -> "k0 是否必须在 CenterGapTolerance 内无隙。",
    "Three distinct symbols used when H is supplied as a matrix." -> "H 以矩阵给出时使用的三个互异动量符号。",
    "ChernNumber or a diagnostic Data Association." -> "返回 ChernNumber 或诊断 Data Association。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "showBrillouinZone"
  }] -> <|
    StringJoin[
      "Run init with the desired lattice before plotting. The first Brillouin zone depends on the ",
      "lattice, so no Hamiltonian or hopping parameters are needed."
    ] -> "绘图前先用 init 指定晶格。第一 Brillouin 区由晶格决定，无需先计算 Hamiltonian 或指定 hopping 参数。",
    "The lattice vectors are the rows of lattice. The plot uses Cartesian reciprocal-space coordinates, with reciprocal vectors including the factor 2 Pi." ->
      "lattice 的每一行是一条实空间基矢。图形使用倒空间笛卡尔坐标，倒格矢包含 2 Pi 因子。",
    "\"KPath\" -> Automatic uses standardKPath[]. An explicit path is a list of {{kStart,kEnd},{labelStart,labelEnd}} segments in reciprocal fractional coordinates; equivalent points are folded into the first zone." ->
      "\"KPath\" -> Automatic 使用 standardKPath[]。自定义路径由 {{起点,终点},{起点标签,终点标签}} 形式的线段组成，坐标为倒格矢的分数坐标；等价点会折回第一 Brillouin 区。",
    "The result is Graphics3D with Axes -> False, Boxed -> False, and ViewProjection -> \"Orthographic\". These are fixed current presentation properties." ->
      "结果是 Graphics3D，并固定采用 Axes -> False、Boxed -> False 与 ViewProjection -> \"Orthographic\"。",
    "showBrillouinZone returns $Failed and issues a message when the reciprocal lattice cannot define a bounded first zone, the supplied path is malformed, or TranslationRange/Tolerance is invalid. It does not repair the input or switch to an approximate fallback." ->
      "当倒格子不能定义有界第一 Brillouin 区、给定路径格式错误，或者 TranslationRange 或 Tolerance 非法时，showBrillouinZone 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Hexagonal lattice: initialize two symmetry-related pz sites using gray group 191. The in-plane primitive vectors enclose 120 degrees:" ->
      "六方晶格：用灰群 191 初始化两个对称等价的 pz 格点。两条面内原胞基矢的夹角为 120°：",
    "The first Brillouin zone is a hexagonal prism. Red lines mark the automatically selected high-symmetry path; Gamma is at the zone center:" ->
      "第一 Brillouin 区呈六棱柱形。红线标出自动选取的高对称路径，Γ 点位于区中心：",
    "\"KPath\" - draw only the G-M-K-G path in the kz = 0 plane of the hexagonal zone. Here G denotes Gamma, and coordinates are fractions of the reciprocal vectors:" ->
      "\"KPath\" - 只画六方 Brillouin 区 kz = 0 平面内的 G-M-K-G 路径。这里 G 表示 Γ 点，坐标为倒格矢的分数坐标：",
    "\"ShowKPath\" - hide the path and labels to see the hexagonal prism by itself:" ->
      "\"ShowKPath\" - 隐藏路径和标签，只显示六棱柱形的 Brillouin 区：",
    "\"TranslationRange\" - include more reciprocal-lattice translations. This hexagonal zone already has the same shape at the default value 2:" ->
      "\"TranslationRange\" - 纳入更多倒格子平移矢量。本例在默认值 2 时已得到相同的六方 Brillouin 区形状：",
    "\"Tolerance\" - change the geometric comparison tolerance. This well-conditioned hexagonal example keeps the same visible shape:" ->
      "\"Tolerance\" - 调整几何比较容差，本例的六方 Brillouin 区形状没有可见变化：",
    "FontSize - enlarge the high-symmetry point labels on the hexagonal zone:" ->
      "FontSize - 放大六方 Brillouin 区的高对称点标签：",
    "ImageSize - draw a smaller hexagonal-zone image:" ->
      "ImageSize - 缩小六方 Brillouin 区图形的显示尺寸：",
    "Primitive tetragonal lattice (gray group 123): a = b and c = 3 a/2 give a square prism. Its height in reciprocal space is 2 Pi/c:" ->
      "简单四方晶格（灰群 123）：a = b、c = 3 a/2，第一 Brillouin 区是底面为正方形的棱柱，倒空间高度为 2 Pi/c：",
    "Primitive orthorhombic lattice (gray group 47): three unequal lattice constants give a rectangular box with edge lengths 2 Pi/a, 2 Pi/b, and 2 Pi/c:" ->
      "简单正交晶格（灰群 47）：三个晶格常数互不相等，第一 Brillouin 区是边长分别为 2 Pi/a、2 Pi/b、2 Pi/c 的长方体：",
    "Monoclinic lattice (gray group 10): the first and third primitive vectors enclose 70 degrees. Hide the path to inspect the zone faces; the Wigner-Seitz zone is not simply the skew reciprocal primitive cell:" ->
      "单斜晶格（灰群 10）：第一、第三条原胞基矢夹角为 70°。隐藏路径后可直接观察各个区面；Wigner-Seitz 区并不等同于倾斜的倒格子原胞：",
    "Automatic conventional path, None, or an explicit MagneticTB band path." -> "可取 Automatic 约定路径、None 或显式 MagneticTB 能带路径。",
    "Whether the selected path and point labels are drawn." -> "是否绘制所选路径与点标签。",
    "Positive integer reciprocal-lattice search range used to construct and fold the zone." -> "构造并折叠 Brillouin 区所用的正整数倒格子搜索范围。",
    "Positive geometric tolerance for zone vertices and path folding." -> "Brillouin 区顶点与路径折叠所用的正几何容差。",
    "Font size used for high-symmetry point labels." -> "高对称点标签所用字号。",
    "Final Graphics3D image size." -> "最终 Graphics3D 图像尺寸。"
  |>,
  FileNameJoin[{
    "ReferencePages", "Symbols", "showCrystalStructure"
  }] -> <|
    "Atomic positions and moments come from the model prepared by init. A moment is transformed with the direct lattice and displayed as a red arrow." ->
      "原子位置与磁矩来自 init 准备的模型。磁矩随正格子变换，并显示为红色箭头。",
    "Symmetry expands each Wyckoff seed into its equivalent atoms. CellRange repeats these atoms in translated cells. Sphere colors distinguish Wyckoff orbits; arrows show moment directions, not their magnitudes." ->
      "对称操作会从每个 Wyckoff 代表位置生成全部等价原子，CellRange 再将它们平移到邻近原胞。球体颜色区分不同的 Wyckoff 轨道；箭头表示磁矩方向，长度不代表磁矩大小。",
    "The result is Graphics3D with Axes -> False, Boxed -> False, and ViewProjection -> \"Orthographic\". These three presentation choices are fixed current behavior rather than public options." ->
      "结果是 Graphics3D，并固定采用 Axes -> False、Boxed -> False 与 ViewProjection -> \"Orthographic\"；这三项是当前固定显示行为，不是公开选项。",
    "showCrystalStructure returns $Failed and issues a message when the current model lacks valid lattice/atom data, CellRange is malformed, a symbolic moment remains unresolved, or a positive size option is invalid. It does not repair the input or switch to an approximate fallback." ->
      "当当前模型缺少有效晶格或原子数据、CellRange 格式错误、符号磁矩未被赋值，或者正尺寸选项非法时，showCrystalStructure 发出消息并返回 $Failed；它不会修补输入或切换到近似 fallback。",
    "Ferromagnetic order: draw 2 by 2 cells in the ab plane. All eight atomic moments point along +z:" ->
      "铁磁序：画出 ab 面内 2×2 个原胞，八个原子的磁矩全部指向 +z：",
    "Antiferromagnetic order: the same eight atoms form two sublattices with equal and opposite moments. The red arrows alternate between +z and -z:" ->
      "反铁磁序：同样的八个原子分属两个子晶格，磁矩大小相等、方向相反。红色箭头交替指向 +z 和 -z：",
    "Ferrimagnetic order: each cell has two +z moments and one -z moment. The four-cell picture contains twelve atoms and has a nonzero net moment:" ->
      "亚铁磁序：每个原胞有两个 +z 磁矩和一个 -z 磁矩，四个原胞共十二个原子，净磁矩不为零：",
    "Noncollinear order: P4 generates four in-plane moments from one seed. Repeating the cell shows sixteen atoms with moments along +x, -x, +y, and -y:" ->
      "非共线磁序：P4 从一个代表原子生成四个面内磁矩。平移原胞后，图中共有十六个原子，磁矩分别沿 +x、-x、+y 和 -y：",
    "\"CellRange\" - show two neighbouring cells of the noncollinear model above. Each cell contains four atoms:" ->
      "\"CellRange\" - 显示上述非共线模型的两个相邻原胞，每个原胞有四个原子：",
    "\"MomentRules\" - P4 expands the seed into four atoms in one cell. Set m = -1 to point all four magnetic moments along -z:" ->
      "\"MomentRules\" - P4 从代表位置生成一个原胞内的四个原子。令 m = -1，使四个磁矩都指向 -z：",
    "\"MomentScale\" - use the same four-atom arrangement with +z moments and shorten the arrows. This changes only their displayed length:" ->
      "\"MomentScale\" - 保留相同的四原子排列，将磁矩设为 +z 并缩短箭头；这个选项只改变箭头的显示长度：",
    "\"AtomRadius\" - reduce the sphere radius for all four atoms:" ->
      "\"AtomRadius\" - 缩小四个原子球体的半径：",
    "\"ShowAtomLabels\" - label the four equivalent atoms 1.1 through 1.4. The first number identifies the Wyckoff orbit:" ->
      "\"ShowAtomLabels\" - 将四个等价原子标为 1.1 至 1.4，第一个数字表示 Wyckoff 轨道编号：",
    "FontSize - enlarge the four atom labels:" ->
      "FontSize - 放大四个原子的标签：",
    "ImageSize - show the four-atom cell at a smaller size:" ->
      "ImageSize - 缩小四原子原胞图的显示尺寸：",
    "Integer cell bounds, or n for translations from -n through n on every axis." -> "整数晶胞边界；也可用 n 表示每个轴从 -n 到 n 的平移。",
    "Rules that make symbolic moment components numerical." -> "把符号磁矩分量数值化的规则。",
    "Automatic or a positive numerical arrow length scale." -> "Automatic 或正数箭头长度比例。",
    "Automatic or a positive numerical sphere radius." -> "Automatic 或正数球体半径。",
    "Whether to label every atom by orbit and equivalent-atom index." -> "是否用轨道编号与等价原子编号标注每个原子。",
    "Font size used when atom labels are shown." -> "显示原子标签时使用的字号。",
    "Final Graphics3D image size." -> "最终 Graphics3D 图像尺寸。"
  |>
|>;

ClearAll[testingChineseText, translateTestingNotebook];
testingChineseText[text_String, style_String, path_String] := Module[
  {summary, option},
  If[
    StringContainsQ[text, RegularExpression["[一-龥]"]],
    Return[text]
  ];
  summary = testingChinesePageSummaries[path];
  Which[
    KeyExistsQ[testingChineseHeadingTranslations, text],
      testingChineseHeadingTranslations[text],
    style === "ExampleText",
      option = StringTrim@First@StringSplit[text, " - "];
      If[StringContainsQ[text, " - "],
        option <> " \:7684\:53ef\:8fd0\:884c\:793a\:4f8b\:ff1a\:663e\:5f0f\:8bbe\:7f6e\:8be5\:9009\:9879\:5e76\:76f4\:63a5\:663e\:793a\:7ed3\:679c\:3002",
        "\:4e0b\:9762\:7ed9\:51fa\:5b8c\:6574\:53ef\:8fd0\:884c\:7684\:8c03\:7528\:ff1b\:8f93\:51fa\:7531\:5f53\:524d testing \:6b63\:5f0f\:63a5\:53e3\:751f\:6210\:3002"
      ],
    style === "TableText",
      StringJoin[
        "\:672c\:884c\:7ed9\:51fa\:8be5\:9009\:9879\:7684\:5141\:8bb8\:5f62\:5f0f\:4e0e\:8ba1\:7b97\:4f5c",
        "\:7528\:ff1bOptions \:5c0f\:8282\:63d0\:4f9b\:663e\:5f0f\:8c03\:7528\:548c\:5b9e\:9645\:8f93\:51fa",
        "\:3002"
      ],
    style === "GuideText",
      " \:2014 \:6253\:5f00\:5bf9\:5e94\:51fd\:6570\:9875\:67e5\:770b\:7528\:9014\:3001\:8fb9\:754c\:3001\:5168\:90e8\:9009\:9879\:4e0e\:53ef\:8fd0\:884c\:793a\:4f8b\:3002",
    style === "GuideTutorial",
      Lookup[testingChineseHeadingTranslations, text, "\:76f8\:5173\:4e13\:9898\:6559\:7a0b"],
    MemberQ[{"Notes", "Text", "GuideAbstract"}, style],
      summary,
    style === "Usage",
      summary,
    True,
      Lookup[testingChineseHeadingTranslations, text, text]
  ]
];

translateTestingNotebook[notebook_Notebook, path_String] := Module[{summary},
  If[!KeyExistsQ[testingChinesePageSummaries, path], Return[notebook]];
  summary = testingChinesePageSummaries[path];
  notebook /.
    Cell[text_String, style_String, options___] /;
        MemberQ[{
          "Title", "Text", "ExampleText", "Notes", "GuideAbstract",
          "GuideText", "GuideTutorial", "Section", "Subsection",
          "GuideFunctionsSection", "TableText"
        }, style] :>
      Cell[testingChineseText[text, style, path], style, options] /.
    Cell[TextData[items_List], style : "Usage", options___] :>
      Cell[
        TextData[Replace[
          items,
          text_String /; StringLength[text] >= 20 :>
            testingChineseText[text, style, path],
          {1}
        ]],
        style,
        options
      ]
];

(* -------------------------------------------------------------------------
   Phase 2: expected-text fragments that must survive translation.

   Generic title/section/URI checks are added automatically below; this table
   contains only page-specific phrases worth guarding explicitly.
   ------------------------------------------------------------------------- *)
additionalExpectedFragments = <|
  FileNameJoin[{"Tutorials", "Corepresentations"}] -> {
    "双值高对称点计算",
    "单值基",
    "getTBBandCorep[",
    "small corep",
    "paclet:MagneticTB/guide/MagneticTB"
  },
  FileNameJoin[{"Tutorials", "CrystalAndKPaths"}] -> {
    "Kagome",
    "Hamiltonian 与能带",
    "typeIII[{191, 239}]",
    "bandplot[",
    "paclet:MagneticTB/guide/MagneticTB"
  },
  FileNameJoin[{"Guides", "MagneticTB"}] -> {
    "教程：从这里开始",
    "模型输入与查看",
    "晶体、对称性与群论",
    "Hamiltonian 与能带计算",
    "能带拟合与 VASP 对比",
    "将 MagneticTB 能带拟合到 VASP",
    "绘图与倒空间路径",
    If[
      FileExistsQ[FileNameJoin[{
        englishSourceRoot,
        "Tutorials",
        "RealSpaceAndTopology.source.wl"
      }]],
      "物性与拓扑",
      Nothing
    ]
  },
  FileNameJoin[{"Tutorials", "BandFitting"}] -> {
    "读取并检查三条 VASP 能带",
    "拟合 onsite 与最近邻模型",
    "加入第二近邻 hopping",
    "只拟合关心的能量范围",
    "只拟合 K 点附近",
    "根据能带图判断，不要只看损失"
  },
  FileNameJoin[{"Tutorials", "ContinuousSymmetry"}] -> {
    "含 SOC 的高对称 MSG 124.360",
    "无 SOC：加入连续 C-infinity 自旋对称性",
    "沿 G-Z 的能带劈裂",
    "bandplot[",
    "连续接口的适用范围",
    "typeIV[{124, 360}]",
    "showCrystalStructure[]"
  },
  FileNameJoin[{"Tutorials", "RealSpaceAndTopology"}] -> {
    "下面使用的 Hamiltonian",
    "Chern 绝缘体 Hamiltonian",
    "Z2 绝缘体 Hamiltonian",
    "Wilson-loop 图",
    "Chern 绝缘体",
    "Z2 绝缘体",
    "表面谱",
    "C4 高阶拓扑绝缘体角态",
    "磁性 Weyl 点",
    "Berry 曲率密度"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "plotWilsonLoop"}] -> {
    "由对称性生成的 P4 Chern 绝缘体",
    "由对称性生成的灰群 P4 Z2 绝缘体"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "plotSurfaceSpectrum"}] -> {
    "由对称性生成的 P4 Chern 绝缘体",
    "由对称性生成的灰群 P4 Z2 模型"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "plotRealSpaceWavefunction"}] -> {
    "二维 C4 呼吸方格高阶拓扑绝缘体",
    "四极矩拓扑不变量"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "pointChernNumber"}] -> {
    "由 BNS 143.3 对称性生成的磁性 Weyl 模型",
    "绘制同一 Weyl 点周围的 Berry 曲率矢量场"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "pointMatrix"}] -> {
    "完整幺正四方 C4 点群"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "readHR"}] -> {
    "buildBlochHamiltonian",
    "对称性生成 hopping 元"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "showBrillouinZone"}] -> {
    "msgop[gray[191]]", "msgop[gray[123]]",
    "msgop[gray[47]]", "msgop[gray[10]]",
    "六棱柱", "简单四方晶格", "简单正交晶格", "单斜晶格"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "showCrystalStructure"}] -> {
    "msgop[bnsdict[{75, 1}]]",
    "铁磁序",
    "反铁磁序",
    "亚铁磁序",
    "非共线磁序",
    "红色箭头"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "showSymmetryRepresentations"}] -> {
    "showSymmetryRepresentations[{1, 3, 25}]"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "initfromrep"}] -> {
    "\:9009\:9879",
    "InitialBondShells -> 2",
    "RepresentationMode -> \\\"Induced\\\"",
    "C4z^-1"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "symham"}] -> {
    "\:9009\:9879",
    "Hermitian -> False",
    "KernelMethod -> \\\"Stacked\\\"",
    "ValidationLevel -> \\\"Full\\\"",
    "CartesianCoordinates -> True",
    "\:628a\:76f8\:540c\:7684\:7cbe\:786e\:5bf9\:79f0\:7ea6\:675f\:5408\:5e76\:4e3a\:4e00\:4e2a\:77e9\:9635\:6c42\:89e3"
  },
  FileNameJoin[{"ReferencePages", "Symbols", "initold"}] -> {
    "\:9009\:9879",
    "debugQold -> False"
  }
|>;

(* -------------------------------------------------------------------------
   Phase 3: translate each selected page and write UTF-8 authoring sources.

   Final .nb canonicalization is intentionally a separate notebook_guard.py
   step, keeping translation and Unicode escaping as independent gates.
   ------------------------------------------------------------------------- *)
Scan[
  Function[path,
    englishPath = FileNameJoin[{englishSourceRoot, path <> ".source.wl"}];
    If[!FileExistsQ[englishPath],
      Print["Missing English source: ", englishPath];
      Exit[4]
    ];
    notebook = Get[englishPath];
    If[Head[notebook] =!= Notebook,
      Print["Invalid English notebook source: ", englishPath];
      Exit[5]
    ];
    translated = translateNotebook[
      notebook,
      Normal@Join[
        commonTranslations,
        Lookup[translations, path, <||>],
        additionalTranslations,
        Lookup[testingChinesePageTranslations, path, <||>]
      ]
    ];
    translated = translateTestingNotebook[translated, path];
    writeUTF8TextIfChanged[
      FileNameJoin[{sourceRoot, path <> ".source.wl"}],
      ToString[
        translated,
        InputForm,
        PageWidth -> If[
          MemberQ[
            FileNameJoin[{
              "ReferencePages",
              "Symbols",
              #
            }] & /@ {"bandManipulate", "bandManipulateold"},
            path
          ],
          120,
          Infinity
        ]
      ]
    ];
    writeUTF8TextIfChanged[
      FileNameJoin[{expectedRoot, path <> ".expected.txt"}],
      StringRiffle[
        DeleteMissing@DeleteDuplicates@Join[
          {
            FirstCase[
              translated,
              Cell[text_String, "ObjectName" | "Title", ___] :> text,
              Missing["Title"],
              Infinity
            ],
            Which[
              StringStartsQ[path, "ReferencePages"],
                "\:8303\:4f8b",
              StringStartsQ[path, "Tutorials"],
                FirstCase[
                  translated,
                  Cell[text_String, "Section", ___] :> text,
                  Missing["Section"],
                  Infinity
                ],
              True,
                FirstCase[
                  translated,
                  Cell[text_String, "GuideFunctionsSection", ___] :> text,
                  Missing["GuideSection"],
                  Infinity
                ]
            ],
            FirstCase[
              translated,
              Cell[text_String, "Categorization", options___] /;
                  !FreeQ[{options}, CellLabel -> "URI"] :> text,
              Missing["URI"],
              Infinity
            ]
          },
          Lookup[additionalExpectedFragments, path, {}]
        ],
        "\n"
      ] <> "\n"
    ];
    Print[FileNameJoin[{sourceRoot, path <> ".source.wl"}]]
  ],
  relativePages
];
