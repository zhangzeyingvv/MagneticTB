(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  normalizeWannierShellSelection,
  normalizeWannierParameterRules,
  integerTranslationVectorQ,
  resolveHamiltonianWannierCenters,
  extractWannierHamiltonianTerm,
  compileHamiltonianWannier90HRData,
  compileWannier90TermRecord,
  compileWannier90ShellRecords,
  compileWannier90HRData,
  formatWannier90Real,
  formatWannier90HRData,
  parseWannierInteger,
  parseWannierReal,
  parseWannier90RecordChunk,
  parseWannier90HRFile,
  filterWannier90HRData,
  resolveWannier90ExportTarget,
  emitWannier90HRData
];

hop::selection =
  "The shell selection must be a positive integer n (meaning shells 1 through n) or a nonempty list of positive shell indices; received `1`.";
hop::shellrange =
  "Shell indices `1` were not prepared by init; only `2` shells are available.";
hop::unsolved =
  "No cached real-space result exists for shell(s) `1` with Hermitian -> `2`, KernelMethod -> `3`, and ValidationLevel -> `4`. Evaluate symham for those shells with the same options first.";
hop::matrix =
  "The hand-written Hamiltonian must be a nonempty square matrix.";
hop::rules =
  "Parameter values must be supplied as a list of rules or an Association; received `1`.";
hop::numeric =
  "The real-space hopping matrices remain nonnumeric after applying the supplied parameter rules. Remaining symbols include `1`.";
hop::data =
  "The cached real-space hopping data for shell `1` is missing or inconsistent with the static bond data.";
hop::centers =
  "The wcc setting must be None, Automatic, or one real fractional three-vector per Hamiltonian row; received `1`.";
hop::phase =
  "A term in Hamiltonian element (`1`,`2`) is not a numeric coefficient times exponentials with an exactly linear phase in {kx,ky,kz}: `3`.";
hop::translation =
  "Hamiltonian element (`1`,`2`) gives the noninteger Wannier translation `3`; check the phase convention or supply the correct wcc centers.";
hop::export =
  "The Wannier90 HR file could not be written to `1`. The target directory must already exist.";
hop::option = "Invalid Wannier90 export option value(s): `1`.";

readHR::file = "The Wannier90 HR file `1` could not be read.";
readHR::format =
  "The file `1` is not a valid complete wannier90_hr.dat file: `2`.";
readHR::option = "Invalid readHR option value(s): `1`.";

normalizeWannierShellSelection[selection_, available_Integer?Positive] :=
  Module[{shells},
    shells = Which[
      IntegerQ[selection] && selection > 0,
        Range[selection],
      ListQ[selection] && selection =!= {} &&
          And @@ (IntegerQ[#] && # > 0 & /@ selection),
        DeleteDuplicates[selection],
      True,
        Return[Failure[
          "InvalidShellSelection",
          <|"Selection" -> selection|>
        ]]
    ];
    If[
      !SubsetQ[Range[available], shells],
      Failure[
        "ShellOutOfRange",
        <|"Shells" -> shells, "Available" -> available|>
      ],
      shells
    ]
  ];

normalizeWannierParameterRules[rules_] := Which[
  AssociationQ[rules],
    Normal[rules],
  ListQ[rules] && And @@ (MatchQ[#, _Rule | _RuleDelayed] & /@ rules),
    rules,
  True,
    $Failed
];

integerTranslationVectorQ[translation_] :=
  ListQ[translation] && Length[translation] == 3 &&
    And @@ (IntegerQ[Simplify[#]] & /@ translation);

resolveHamiltonianWannierCenters[setting_, numWannier_Integer?Positive] :=
  Module[{centers, model, siteOrbits, localDimensions},
    centers = Which[
      setting === None,
        ConstantArray[{0, 0, 0}, numWannier],
      setting === Automatic && AssociationQ[$CurrentModelSession],
        model = Lookup[
          $CurrentModelSession,
          "ModelSpecification",
          $Failed
        ];
        If[!AssociationQ[model], Return[$Failed, Module]];
        siteOrbits = Lookup[model, "SiteOrbits", $Failed];
        localDimensions = Lookup[model, "LocalDimensions", $Failed];
        If[
          !ListQ[siteOrbits] || !ListQ[localDimensions] ||
            Length[siteOrbits] =!= Length[localDimensions],
          Return[$Failed, Module]
        ];
        Flatten[
          MapThread[
            Function[{orbit, dimension},
              Flatten[ConstantArray[#, dimension] & /@ orbit, 1]
            ],
            {siteOrbits, localDimensions}
          ],
          1
        ],
      setting === Automatic,
        ConstantArray[{0, 0, 0}, numWannier],
      ListQ[setting],
        setting,
      True,
        $Failed
    ];
    If[
      !ListQ[centers] || Length[centers] =!= numWannier ||
        !And @@ (
          VectorQ[#, NumericQ] && Length[#] == 3 &&
            Max[Abs[Im[N[#]]]] == 0. & /@ centers
        ),
      $Failed,
      centers
    ]
  ];

extractWannierHamiltonianTerm[
    term_,
    row_Integer?Positive,
    column_Integer?Positive,
    centers_List,
    momentumSymbols_List,
    tolerance_?NumericQ
  ] := Module[
  {
    factors, phaseArguments, coefficientFactors, coefficient,
    linearPhase, phaseVector, remainder, rawTranslation,
    roundedTranslation
  },
  factors = If[Head[term] === Times, List @@ term, {term}];
  phaseArguments = Cases[
    factors,
    HoldPattern[Power[E, argument_]] :> argument
  ];
  coefficientFactors = DeleteCases[
    factors,
    HoldPattern[Power[E, _]]
  ];
  coefficient = Times @@ coefficientFactors;
  If[
    !FreeQ[coefficient, HoldPattern[Power[E, _]]] ||
      !And @@ (FreeQ[coefficient, #] & /@ momentumSymbols) ||
      !NumericQ[coefficient],
    Return[Failure[
      "InvalidHamiltonianTerm",
      <|"Row" -> row, "Column" -> column, "Term" -> term|>
    ]]
  ];
  linearPhase = Expand[Total[phaseArguments]/I];
  phaseVector = Coefficient[linearPhase, #] & /@ momentumSymbols;
  remainder = Expand[linearPhase - phaseVector.momentumSymbols];
  If[
    remainder =!= 0 || !VectorQ[phaseVector, NumericQ] ||
      Max[Abs[Im[N[phaseVector]]]] > tolerance,
    Return[Failure[
      "InvalidHamiltonianTerm",
      <|"Row" -> row, "Column" -> column, "Term" -> term|>
    ]]
  ];
  rawTranslation = N[
    Re[phaseVector] + centers[[row]] - centers[[column]]
  ];
  roundedTranslation = Round[rawTranslation];
  If[Max[Abs[rawTranslation - roundedTranslation]] > tolerance,
    Return[Failure[
      "NonintegerTranslation",
      <|
        "Row" -> row,
        "Column" -> column,
        "Translation" -> rawTranslation
      |>
    ]]
  ];
  <|
    "Row" -> row,
    "Column" -> column,
    "Translation" -> roundedTranslation,
    "Coefficient" -> Chop[N[coefficient], tolerance]
  |>
];

compileHamiltonianWannier90HRData[
    hamiltonian_?MatrixQ,
    parameterRules_List,
    centers_List,
    tolerance_?NumericQ
  ] := Module[
  {
    numWannier, numericalHamiltonian, recordsByElement, records,
    failure, nonzeroRecords, groups, translations, matrices
  },
  If[
    hamiltonian === {} || Length[Dimensions[hamiltonian]] =!= 2 ||
      First[Dimensions[hamiltonian]] =!= Last[Dimensions[hamiltonian]],
    Return[Failure["InvalidHamiltonian", <||>]]
  ];
  numWannier = Length[hamiltonian];
  numericalHamiltonian = hamiltonian /. parameterRules;
  recordsByElement = Flatten@Table[
    Map[
      extractWannierHamiltonianTerm[
        #,
        row,
        column,
        centers,
        {kx, ky, kz},
        tolerance
      ] &,
      With[
        {expanded = Expand@TrigToExp[numericalHamiltonian[[row, column]]]},
        If[Head[expanded] === Plus, List @@ expanded, {expanded}]
      ]
    ],
    {row, numWannier},
    {column, numWannier}
  ];
  failure = SelectFirst[recordsByElement, FailureQ, Missing["NoFailure"]];
  If[!MissingQ[failure], Return[failure]];
  records = recordsByElement;
  nonzeroRecords = Select[
    records,
    Abs[N[Lookup[#, "Coefficient"]]] > tolerance &
  ];
  If[nonzeroRecords === {},
    Return[<|
      "Schema" -> "Wannier90HRData",
      "SchemaVersion" -> 1,
      "Source" -> "HamiltonianCoefficientExtraction",
      "NumWannier" -> numWannier,
      "NumTranslations" -> 1,
      "Translations" -> {{0, 0, 0}},
      "Degeneracies" -> {1},
      "HoppingMatrices" -> {
        ConstantArray[0., {numWannier, numWannier}]
      }
    |>]
  ];
  groups = SortBy[
    GatherBy[nonzeroRecords, Lookup[#, "Translation"] &],
    Lookup[First[#], "Translation"] &
  ];
  translations = Lookup[First[#], "Translation"] & /@ groups;
  matrices = Map[
    Function[group,
      Table[
        Total@Lookup[
          Select[
            group,
            Lookup[#, "Row"] == row &&
              Lookup[#, "Column"] == column &
          ],
          "Coefficient",
          0
        ],
        {row, numWannier},
        {column, numWannier}
      ]
    ],
    groups
  ];
  <|
    "Schema" -> "Wannier90HRData",
    "SchemaVersion" -> 1,
    "Source" -> "HamiltonianCoefficientExtraction",
    "NumWannier" -> numWannier,
    "NumTranslations" -> Length[translations],
    "Translations" -> translations,
    "Degeneracies" -> ConstantArray[1, Length[translations]],
    "HoppingMatrices" -> matrices
  |>
];

compileWannier90TermRecord[
    shell_Integer?Positive,
    bonds_List,
    sites_List,
    blockDimensions_List,
    fullDimension_Integer?Positive,
    rawTerm_Association,
    parameterRules_List
  ] := Module[
  {
    bondIndex, rowBlock, columnBlock, hopping, remainingSymbols,
    bond, translation, rowRange, columnRange, embedded
  },
  bondIndex = Lookup[rawTerm, "BondIndex", $Failed];
  rowBlock = Lookup[rawTerm, "RowBlock", $Failed];
  columnBlock = Lookup[rawTerm, "ColumnBlock", $Failed];
  hopping = Lookup[rawTerm, "Matrix", $Failed] /. parameterRules;
  If[
    !IntegerQ[bondIndex] || !Between[bondIndex, {1, Length[bonds]}] ||
      !IntegerQ[rowBlock] || !Between[rowBlock, {1, Length[sites]}] ||
      !IntegerQ[columnBlock] ||
        !Between[columnBlock, {1, Length[sites]}] ||
      !MatrixQ[hopping] ||
      Dimensions[hopping] =!= {
        blockDimensions[[rowBlock]],
        blockDimensions[[columnBlock]]
      },
    Return[Failure["InvalidShellData", <|"Shell" -> shell|>]]
  ];
  If[!MatrixQ[hopping, NumericQ],
    remainingSymbols = DeleteDuplicates@Cases[
      hopping,
      symbol_Symbol /;
        Context[Unevaluated[symbol]] =!= "System`",
      Infinity
    ];
    Return[Failure[
      "NonnumericHopping",
      <|"Symbols" -> remainingSymbols|>
    ]]
  ];
  bond = bonds[[bondIndex]];
  If[
    Lookup[bond, "RowSite", $Failed] =!= rowBlock ||
      Lookup[bond, "ColumnSite", $Failed] =!= columnBlock,
    Return[Failure["InvalidShellData", <|"Shell" -> shell|>]]
  ];
  translation = Simplify[
    Lookup[bond, "ColumnCell", $Failed] -
      Lookup[bond, "RowCell", $Failed]
  ];
  If[!integerTranslationVectorQ[translation],
    Return[Failure["InvalidShellData", <|"Shell" -> shell|>]]
  ];
  translation = Round[translation];
  rowRange = Lookup[sites[[rowBlock]], "BlockRange", $Failed];
  columnRange = Lookup[sites[[columnBlock]], "BlockRange", $Failed];
  If[
    !ListQ[rowRange] || !ListQ[columnRange] ||
      Length[rowRange] =!= blockDimensions[[rowBlock]] ||
      Length[columnRange] =!= blockDimensions[[columnBlock]],
    Return[Failure["InvalidShellData", <|"Shell" -> shell|>]]
  ];
  embedded = ConstantArray[0, {fullDimension, fullDimension}];
  embedded[[rowRange, columnRange]] = hopping;
  <|
    "Shell" -> shell,
    "BondIndex" -> bondIndex,
    "Translation" -> translation,
    "Matrix" -> embedded
  |>
];

compileWannier90ShellRecords[
    shell_Integer?Positive,
    staticShell_Association,
    shellResult_Association,
    parameterRules_List,
    referenceBlockDimensions_List,
    fullDimension_Integer?Positive
  ] := Module[{bonds, sites, blockDimensions, records, failure},
  If[Lookup[staticShell, "Shell", Missing["Shell"]] =!= shell,
    Return[Failure["InvalidShellData", <|"Shell" -> shell|>]]
  ];
  bonds = Lookup[staticShell, "Bonds", $Failed];
  sites = Lookup[staticShell, "Sites", $Failed];
  blockDimensions = Lookup[staticShell, "BlockDimensions", $Failed];
  If[
    !ListQ[bonds] || !ListQ[sites] ||
      blockDimensions =!= referenceBlockDimensions ||
      Length[sites] =!= Length[blockDimensions] ||
      !ListQ[Lookup[shellResult, "Terms", $Failed]],
    Return[Failure["InvalidShellData", <|"Shell" -> shell|>]]
  ];
  records = compileWannier90TermRecord[
    shell,
    bonds,
    sites,
    blockDimensions,
    fullDimension,
    #,
    parameterRules
  ] & /@ shellResult["Terms"];
  failure = SelectFirst[records, FailureQ, Missing["NoFailure"]];
  If[MissingQ[failure], records, failure]
];

compileWannier90HRData[
    session_Association,
    shells_List,
    parameterRules_List,
    hermitianQ : (True | False),
    method_String,
    validationLevel_String
  ] := Module[
  {
    realSpaceCache, staticShells, cacheKeys, missingShells,
    shellResults, referenceBlockDimensions, fullDimension,
    recordsByShell, recordFailure, termRecords, groups,
    translations, matrices
  },
  realSpaceCache = Lookup[session, "RealSpaceShellCache", $Failed];
  staticShells = Lookup[session, "CompiledBondShells", $Failed];
  If[!AssociationQ[realSpaceCache] || !ListQ[staticShells],
    Return[$Failed]
  ];
  cacheKeys = realSpaceShellCacheKey[
    #,
    hermitianQ,
    method,
    validationLevel
  ] & /@ shells;
  missingShells = Pick[
    shells,
    MissingQ[Lookup[realSpaceCache, #, Missing["NotCached"]]] & /@
      cacheKeys
  ];
  If[missingShells =!= {},
    Return[Failure[
      "UnsolvedShells",
      <|"Shells" -> missingShells|>
    ]]
  ];
  shellResults = Lookup[realSpaceCache, cacheKeys];
  If[
    !And @@ (
      AssociationQ[#] &&
        Lookup[#, "Schema", None] === "MagneticTBRealSpaceShellResult" &&
        ListQ[Lookup[#, "Terms", $Failed]] & /@ shellResults
    ),
    Return[$Failed]
  ];
  referenceBlockDimensions = Lookup[
    staticShells[[First[shells]]],
    "BlockDimensions",
    $Failed
  ];
  If[
    !ListQ[referenceBlockDimensions] ||
      !And @@ (IntegerQ[#] && # > 0 & /@ referenceBlockDimensions) ||
      !And @@ (
        Lookup[staticShells[[#]], "BlockDimensions", $Failed] ===
          referenceBlockDimensions & /@ shells
      ),
    Return[$Failed]
  ];
  fullDimension = Total[referenceBlockDimensions];
  recordsByShell = MapThread[
    compileWannier90ShellRecords[
      #1,
      staticShells[[#1]],
      #2,
      parameterRules,
      referenceBlockDimensions,
      fullDimension
    ] &,
    {shells, shellResults}
  ];
  recordFailure = SelectFirst[
    recordsByShell,
    FailureQ,
    Missing["NoFailure"]
  ];
  If[!MissingQ[recordFailure], Return[recordFailure]];
  termRecords = Flatten[recordsByShell, 1];
  If[termRecords === {}, Return[$Failed]];

  groups = SortBy[
    GatherBy[termRecords, Lookup[#, "Translation"] &],
    Lookup[First[#], "Translation"] &
  ];
  translations = Lookup[First[#], "Translation"] & /@ groups;
  matrices = Total[Lookup[#, "Matrix"]] & /@ groups;
  If[!And @@ (MatrixQ[#, NumericQ] & /@ matrices), Return[$Failed]];

  <|
    "Schema" -> "Wannier90HRData",
    "SchemaVersion" -> 2,
    "Source" -> "MagneticTBRealSpaceShellCache",
    "Shells" -> shells,
    "Hermitian" -> hermitianQ,
    "NumWannier" -> fullDimension,
    "NumTranslations" -> Length[translations],
    "Translations" -> translations,
    "Degeneracies" -> ConstantArray[1, Length[translations]],
    "HoppingMatrices" -> matrices
  |>
];

formatWannier90Real[value_?NumericQ, digits_Integer?Positive] := Module[
  {numericValue},
  numericValue = If[Abs[N[value]] < 10.^(-digits), 0., N[value]];
  StringTrim@ToString[
    NumberForm[
      numericValue,
      {30, digits},
      NumberPadding -> {"", "0"},
      NumberPoint -> ".",
      ExponentFunction -> (Null &)
    ],
    OutputForm
  ]
];

formatWannier90HRData[data_Association, digits_Integer?Positive] := Module[
  {
    numWannier, translations, degeneracies, matrices,
    degeneracyLines, recordLines
  },
  numWannier = Lookup[data, "NumWannier", $Failed];
  translations = Lookup[data, "Translations", $Failed];
  degeneracies = Lookup[data, "Degeneracies", $Failed];
  matrices = Lookup[data, "HoppingMatrices", $Failed];
  If[
    !IntegerQ[numWannier] || numWannier < 1 ||
      !ListQ[translations] || translations === {} ||
      !And @@ (integerTranslationVectorQ /@ translations) ||
      !ListQ[degeneracies] ||
      Length[degeneracies] =!= Length[translations] ||
      !And @@ (IntegerQ[#] && # > 0 & /@ degeneracies) ||
      !ListQ[matrices] || Length[matrices] =!= Length[translations] ||
      !And @@ (
        MatrixQ[#, NumericQ] &&
          Dimensions[#] === {numWannier, numWannier} & /@ matrices
      ),
    Return[$Failed]
  ];
  degeneracyLines = StringRiffle[
    StringRiffle[ToString /@ #, " "] & /@
      Partition[degeneracies, UpTo[15]],
    "\n"
  ];
  recordLines = Flatten@Table[
    StringRiffle[
      Join[
        ToString /@ Join[
          translations[[translationIndex]],
          {row, column}
        ],
        {
          formatWannier90Real[
            Re[matrices[[translationIndex, row, column]]],
            digits
          ],
          formatWannier90Real[
            Im[matrices[[translationIndex, row, column]]],
            digits
          ]
        }
      ],
      " "
    ],
    {translationIndex, Length[translations]},
    {column, numWannier},
    {row, numWannier}
  ];
  StringRiffle[
    Join[
      {
        "Generated by MagneticTB",
        ToString[numWannier],
        ToString[Length[translations]],
        degeneracyLines
      },
      recordLines
    ],
    "\n"
  ] <> "\n"
];

parseWannierInteger[token_String] := Module[{trimmed, sign, digits},
  trimmed = StringTrim[token];
  If[
    !StringMatchQ[trimmed, RegularExpression["[+-]?[0-9]+"]],
    Return[$Failed]
  ];
  sign = If[StringStartsQ[trimmed, "-"], -1, 1];
  digits = StringReplace[
    trimmed,
    StartOfString ~~ ("+" | "-") -> ""
  ];
  sign FromDigits[digits]
];

parseWannierReal[token_String] := Module[
  {
    trimmed, exponentParts, mantissa, exponent, sign, unsigned,
    pointPositions, pointPosition, integerDigits, fractionalDigits,
    integerPart, fractionalPart
  },
  trimmed = StringTrim[token];
  If[
    !StringMatchQ[
      trimmed,
      RegularExpression[
        "[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?"
      ]
    ],
    Return[$Failed]
  ];
  exponentParts = StringSplit[trimmed, RegularExpression["[eEdD]"]];
  If[Length[exponentParts] > 2, Return[$Failed]];
  mantissa = First[exponentParts];
  exponent = If[
    Length[exponentParts] == 2,
    parseWannierInteger[Last[exponentParts]],
    0
  ];
  If[exponent === $Failed, Return[$Failed]];
  sign = If[StringStartsQ[mantissa, "-"], -1, 1];
  unsigned = StringReplace[
    mantissa,
    StartOfString ~~ ("+" | "-") -> ""
  ];
  pointPositions = StringPosition[unsigned, "."];
  If[pointPositions === {},
    integerDigits = unsigned;
    fractionalDigits = "",
    pointPosition = pointPositions[[1, 1]];
    integerDigits = StringTake[unsigned, pointPosition - 1];
    fractionalDigits = StringDrop[unsigned, pointPosition]
  ];
  integerPart = If[integerDigits === "", 0, FromDigits[integerDigits]];
  fractionalPart = If[
    fractionalDigits === "",
    0,
    FromDigits[fractionalDigits]/10^StringLength[fractionalDigits]
  ];
  N[sign (integerPart + fractionalPart) 10^exponent]
];

parseWannier90RecordChunk[chunk_List, numWannier_Integer?Positive] :=
  Module[{translation, expectedPairs, pairList, matrix},
    translation = chunk[[1, 1 ;; 3]];
    If[
      !And @@ (#[[1 ;; 3]] === translation & /@ chunk),
      Return[Failure["MixedTranslationBlock", <||>]]
    ];
    expectedPairs = Sort@Flatten[
      Table[{row, column}, {column, numWannier}, {row, numWannier}],
      1
    ];
    pairList = Sort[chunk[[All, {4, 5}]]];
    If[pairList =!= expectedPairs,
      Return[Failure["InvalidMatrixIndices", <||>]]
    ];
    matrix = ConstantArray[0., {numWannier, numWannier}];
    Scan[
      Function[record,
        matrix[[record[[4]], record[[5]]]] =
          record[[6]] + I record[[7]]
      ],
      chunk
    ];
    <|"Translation" -> translation, "Matrix" -> matrix|>
  ];

parseWannier90HRFile[file_String] := Module[
  {
    lines, numWannier, numTranslations, degeneracyLineCount,
    degeneracyTokens, degeneracies, recordStart, expectedRecordCount,
    recordLines, parsedRecords, values, chunks, parsedChunks,
    chunkFailure, translations, matrices
  },
  lines = Quiet@Check[Import[file, "Lines"], $Failed];
  If[lines === $Failed || !ListQ[lines],
    Return[Failure["UnreadableFile", <||>]]
  ];
  While[lines =!= {} && StringTrim[Last[lines]] === "", lines = Most[lines]];
  If[Length[lines] < 4, Return[Failure["TooFewLines", <||>]]];
  numWannier = If[
    Length[StringSplit[lines[[2]]]] == 1,
    parseWannierInteger[First[StringSplit[lines[[2]]]]],
    $Failed
  ];
  numTranslations = If[
    Length[StringSplit[lines[[3]]]] == 1,
    parseWannierInteger[First[StringSplit[lines[[3]]]]],
    $Failed
  ];
  If[
    !IntegerQ[numWannier] || numWannier < 1 ||
      !IntegerQ[numTranslations] || numTranslations < 1,
    Return[Failure["InvalidHeaderCounts", <||>]]
  ];
  degeneracyLineCount = Ceiling[numTranslations/15];
  If[Length[lines] < 3 + degeneracyLineCount,
    Return[Failure["MissingDegeneracies", <||>]]
  ];
  degeneracyTokens = Flatten[
    StringSplit /@ Take[lines, {4, 3 + degeneracyLineCount}]
  ];
  If[Length[degeneracyTokens] =!= numTranslations,
    Return[Failure["InvalidDegeneracyCount", <||>]]
  ];
  degeneracies = parseWannierInteger /@ degeneracyTokens;
  If[
    MemberQ[degeneracies, $Failed] ||
      !And @@ (IntegerQ[#] && # > 0 & /@ degeneracies),
    Return[Failure["InvalidDegeneracies", <||>]]
  ];

  recordStart = 4 + degeneracyLineCount;
  expectedRecordCount = numTranslations numWannier^2;
  recordLines = Drop[lines, recordStart - 1];
  If[Length[recordLines] =!= expectedRecordCount,
    Return[Failure[
      "InvalidRecordCount",
      <|
        "Expected" -> expectedRecordCount,
        "Actual" -> Length[recordLines]
      |>
    ]]
  ];
  parsedRecords = Map[
    Function[line,
      values = StringSplit[line];
      If[Length[values] =!= 7,
        $Failed,
        Join[
          parseWannierInteger /@ Take[values, 5],
          parseWannierReal /@ Drop[values, 5]
        ]
      ]
    ],
    recordLines
  ];
  If[MemberQ[parsedRecords, $Failed, Infinity],
    Return[Failure["InvalidRecord", <||>]]
  ];
  chunks = Partition[parsedRecords, numWannier^2];
  parsedChunks = parseWannier90RecordChunk[#, numWannier] & /@ chunks;
  chunkFailure = SelectFirst[
    parsedChunks,
    FailureQ,
    Missing["NoFailure"]
  ];
  If[!MissingQ[chunkFailure], Return[chunkFailure]];
  translations = Lookup[parsedChunks, "Translation"];
  matrices = Lookup[parsedChunks, "Matrix"];
  If[DuplicateFreeQ[translations] =!= True,
    Return[Failure["DuplicateTranslations", <||>]]
  ];
  <|
    "Schema" -> "Wannier90HRData",
    "SchemaVersion" -> 1,
    "Source" -> ExpandFileName[file],
    "Header" -> First[lines],
    "NumWannier" -> numWannier,
    "NumTranslations" -> numTranslations,
    "Translations" -> translations,
    "Degeneracies" -> degeneracies,
    "HoppingMatrices" -> matrices
  |>
];

filterWannier90HRData[
    data_Association,
    cutoff_,
    precision_?NumericQ
  ] := Module[{indices, translations, matrices},
  translations = data["Translations"];
  indices = If[
    cutoff === All,
    Range[Length[translations]],
    Select[
      Range[Length[translations]],
      Max[Abs[translations[[#]]]] <= cutoff &
    ]
  ];
  matrices = Chop[data["HoppingMatrices"][[indices]], precision];
  Join[
    data,
    <|
      "NumTranslations" -> Length[indices],
      "Translations" -> data["Translations"][[indices]],
      "Degeneracies" -> data["Degeneracies"][[indices]],
      "HoppingMatrices" -> matrices,
      "CellCutoff" -> cutoff,
      "PrecisionCutoff" -> precision
    |>
  ]
];

resolveWannier90ExportTarget[target_String] := Which[
  DirectoryQ[target],
    FileNameJoin[{target, "wannier90_hr.dat"}],
  DirectoryQ[DirectoryName[ExpandFileName[target]]],
    target,
  True,
    $Failed
];

emitWannier90HRData[
    data_Association,
    exportTarget_,
    realDigits_Integer?Positive
  ] := Module[{text, resolvedTarget, exported},
  text = formatWannier90HRData[data, realDigits];
  If[text === $Failed, Return[Failure["InvalidHRData", <||>]]];
  If[exportTarget === None,
    Print[text];
    Return[text]
  ];
  resolvedTarget = resolveWannier90ExportTarget[exportTarget];
  If[resolvedTarget === $Failed,
    Return[Failure[
      "InvalidExportTarget",
      <|"Target" -> exportTarget|>
    ]]
  ];
  exported = Quiet@Check[Export[resolvedTarget, text, "Text"], $Failed];
  If[
    exported === $Failed,
    Failure["ExportFailed", <|"Target" -> resolvedTarget|>],
    resolvedTarget
  ]
];

Options[hop] = {
  "Hermitian" -> True,
  "KernelMethod" -> "Iterative",
  "ValidationLevel" -> "Basic",
  "hrExport" -> None,
  "RealDigits" -> 12,
  "wcc" -> Automatic,
  "TranslationTolerance" -> 10^-9
};

hop[hamiltonian_?MatrixQ, rules_, OptionsPattern[]] := Module[
  {
    parameterRules, exportTarget, realDigits, centerSetting,
    tolerance, centers, data, emitted
  },
  If[
    hamiltonian === {} || Length[Dimensions[hamiltonian]] =!= 2 ||
      First[Dimensions[hamiltonian]] =!= Last[Dimensions[hamiltonian]],
    Message[hop::matrix];
    Return[$Failed]
  ];
  parameterRules = normalizeWannierParameterRules[rules];
  If[parameterRules === $Failed,
    Message[hop::rules, rules];
    Return[$Failed]
  ];
  exportTarget = OptionValue["hrExport"];
  realDigits = OptionValue["RealDigits"];
  centerSetting = OptionValue["wcc"];
  tolerance = OptionValue["TranslationTolerance"];
  If[
    !(exportTarget === None || StringQ[exportTarget]) ||
      !IntegerQ[realDigits] || realDigits < 6 ||
      !NumericQ[tolerance] || !TrueQ[tolerance > 0],
    Message[hop::option, {exportTarget, realDigits, tolerance}];
    Return[$Failed]
  ];
  centers = resolveHamiltonianWannierCenters[
    centerSetting,
    Length[hamiltonian]
  ];
  If[centers === $Failed,
    Message[hop::centers, centerSetting];
    Return[$Failed]
  ];
  data = compileHamiltonianWannier90HRData[
    hamiltonian,
    parameterRules,
    centers,
    tolerance
  ];
  Which[
    MatchQ[data, Failure["InvalidHamiltonianTerm", _Association]],
      Message[
        hop::phase,
        data[[2, "Row"]],
        data[[2, "Column"]],
        data[[2, "Term"]]
      ];
      Return[$Failed],
    MatchQ[data, Failure["NonintegerTranslation", _Association]],
      Message[
        hop::translation,
        data[[2, "Row"]],
        data[[2, "Column"]],
        data[[2, "Translation"]]
      ];
      Return[$Failed],
    FailureQ[data],
      Message[hop::matrix];
      Return[$Failed]
  ];
  emitted = emitWannier90HRData[data, exportTarget, realDigits];
  If[FailureQ[emitted],
    Message[
      hop::export,
      Lookup[emitted[[2]], "Target", exportTarget]
    ];
    Return[$Failed]
  ];
  emitted
];

hop[selection_, rules_, OptionsPattern[]] := Module[
  {
    session, availableShells, shells, parameterRules, hermitianQ,
    method, validationLevel, exportTarget, realDigits, data, emitted
  },
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  availableShells = Length[session["CompiledBondShells"]];
  shells = normalizeWannierShellSelection[selection, availableShells];
  Which[
    MatchQ[shells, Failure["InvalidShellSelection", _Association]],
      Message[hop::selection, selection];
      Return[$Failed],
    MatchQ[shells, Failure["ShellOutOfRange", _Association]],
      Message[hop::shellrange, shells[[2, "Shells"]], availableShells];
      Return[$Failed]
  ];
  parameterRules = normalizeWannierParameterRules[rules];
  If[parameterRules === $Failed,
    Message[hop::rules, rules];
    Return[$Failed]
  ];
  hermitianQ = OptionValue["Hermitian"];
  method = OptionValue["KernelMethod"];
  validationLevel = OptionValue["ValidationLevel"];
  exportTarget = OptionValue["hrExport"];
  realDigits = OptionValue["RealDigits"];
  If[
    !MemberQ[{True, False}, hermitianQ] || !StringQ[method] ||
      !StringQ[validationLevel] ||
      !(exportTarget === None || StringQ[exportTarget]) ||
      !IntegerQ[realDigits] || realDigits < 6,
    Message[hop::option, {
      hermitianQ, method, validationLevel, exportTarget, realDigits
    }];
    Return[$Failed]
  ];
  data = compileWannier90HRData[
    session,
    shells,
    parameterRules,
    hermitianQ,
    method,
    validationLevel
  ];
  Which[
    MatchQ[data, Failure["UnsolvedShells", _Association]],
      Message[
        hop::unsolved,
        data[[2, "Shells"]],
        hermitianQ,
        method,
        validationLevel
      ];
      Return[$Failed],
    MatchQ[data, Failure["NonnumericHopping", _Association]],
      Message[hop::numeric, data[[2, "Symbols"]]];
      Return[$Failed],
    MatchQ[data, Failure["InvalidShellData", _Association]],
      Message[hop::data, data[[2, "Shell"]]];
      Return[$Failed],
    data === $Failed || FailureQ[data],
      Message[hop::data, shells];
      Return[$Failed]
  ];
  emitted = emitWannier90HRData[data, exportTarget, realDigits];
  If[FailureQ[emitted],
    Message[
      hop::export,
      Lookup[emitted[[2]], "Target", exportTarget]
    ];
    Return[$Failed]
  ];
  emitted
];

Options[readHR] = {
  "prec" -> 10^-10,
  "ncell" -> All
};

readHR[file_String, OptionsPattern[]] := Module[
  {precision, cellCutoff, data},
  precision = OptionValue["prec"];
  cellCutoff = OptionValue["ncell"];
  If[
    !NumericQ[precision] || !TrueQ[precision >= 0] ||
      !(cellCutoff === All ||
        (IntegerQ[cellCutoff] && cellCutoff >= 0)),
    Message[readHR::option, {precision, cellCutoff}];
    Return[$Failed]
  ];
  data = parseWannier90HRFile[file];
  Which[
    MatchQ[data, Failure["UnreadableFile", _Association]],
      Message[readHR::file, file];
      $Failed,
    FailureQ[data],
      Message[readHR::format, file, First[data]];
      $Failed,
    True,
      filterWannier90HRData[data, cellCutoff, precision]
  ]
];

readHR[file_, OptionsPattern[]] := (
  Message[readHR::file, file];
  $Failed
);

End[]
EndPackage[]
