(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  hoppingFiniteNumberQ,
  hoppingRealVectorQ,
  hoppingRealMatrixQ,
  hoppingIntegerVectorQ,
  hoppingCellMatrixQ,
  hoppingResolveData,
  hoppingDataQ,
  hoppingEffectiveMatrices,
  hoppingTranslationKey,
  hoppingSparseRules,
  hoppingMatrixNonzeroQ,
  hoppingHermiticityResidual,
  hoppingCosetKey,
  hoppingCellRepresentatives,
  hoppingEmbedBlock,
  hoppingOptionNames,
  hoppingUnknownOptions,
  finiteRealSpaceBasisData,
  hoppingData,
  transformHoppings,
  buildBlochHamiltonian,
  buildRealSpaceHamiltonian,
  buildSlabHamiltonian,
  hoppingBlochHamiltonian
];

hoppingFiniteNumberQ[value_] :=
  NumericQ[value] && Quiet@Check[
    FreeQ[N[value], Indeterminate | ComplexInfinity | DirectedInfinity],
    False
  ];

hoppingRealVectorQ[vector_, length_Integer?Positive] :=
  ListQ[vector] && Length[vector] === length &&
    VectorQ[vector, hoppingFiniteNumberQ] &&
    TrueQ[Max[Abs[Im[N[vector]]]] == 0];

hoppingRealMatrixQ[matrix_, dimensions_List] :=
  MatrixQ[matrix, hoppingFiniteNumberQ] &&
    Dimensions[matrix] === dimensions &&
    TrueQ[Max[Abs[Im[N[matrix]]]] == 0];

hoppingIntegerVectorQ[vector_] :=
  ListQ[vector] && Length[vector] === 3 && VectorQ[vector, IntegerQ];

hoppingCellMatrixQ[matrix_] :=
  MatrixQ[matrix, IntegerQ] && Dimensions[matrix] === {3, 3} &&
    Det[matrix] =!= 0;

hoppingResolveData[input_Association] := input;
hoppingResolveData[input_String] := readHR[input];
hoppingResolveData[_] := $Failed;

hoppingDataQ[data_] := Module[
  {dimension, translations, degeneracies, matrices},
  If[!AssociationQ[data], Return[False]];
  dimension = Lookup[data, "NumWannier", $Failed];
  translations = Lookup[data, "Translations", $Failed];
  degeneracies = Lookup[data, "Degeneracies", $Failed];
  matrices = Lookup[data, "HoppingMatrices", $Failed];
  IntegerQ[dimension] && dimension > 0 &&
    ListQ[translations] && translations =!= {} &&
    And @@ (hoppingIntegerVectorQ /@ translations) &&
    DuplicateFreeQ[translations] &&
    ListQ[degeneracies] &&
    Length[degeneracies] === Length[translations] &&
    And @@ (IntegerQ[#] && # > 0 & /@ degeneracies) &&
    ListQ[matrices] && Length[matrices] === Length[translations] &&
    And @@ (
      MatrixQ[#, hoppingFiniteNumberQ] &&
        Dimensions[#] === {dimension, dimension} & /@ matrices
    )
];

hoppingEffectiveMatrices[data_Association] := MapThread[
  #1/#2 &,
  {data["HoppingMatrices"], data["Degeneracies"]}
];

hoppingTranslationKey[translation_List] :=
  ToString[translation, InputForm];

hoppingSparseRules[matrix_] := Cases[
  ArrayRules[SparseArray[matrix]],
  Rule[{row_Integer, column_Integer}, value_] /;
      !TrueQ[PossibleZeroQ[value]] :>
    Rule[{row, column}, value]
];

hoppingMatrixNonzeroQ[matrix_] := hoppingSparseRules[matrix] =!= {};

hoppingHermiticityResidual[data_Association] := Module[
  {translations, matrices, dimension, matrixByTranslation, zero, residuals},
  translations = data["Translations"];
  matrices = hoppingEffectiveMatrices[data];
  dimension = data["NumWannier"];
  zero = SparseArray[{}, {dimension, dimension}];
  matrixByTranslation = AssociationThread[
    hoppingTranslationKey /@ translations,
    matrices
  ];
  residuals = MapThread[
    Function[{translation, matrix},
      Max@Abs@Flatten@N[
        matrix - ConjugateTranspose@Lookup[
          matrixByTranslation,
          hoppingTranslationKey[-translation],
          zero
        ]
      ]
    ],
    {translations, matrices}
  ];
  Max[Prepend[residuals, 0.]]
];

(* For row-vector lattice coordinates and a new-cell matrix T, Mathematica's
   Smith decomposition satisfies left.T.right = diagonal. Therefore
   Mod[oldCell.right, diagonal] is the exact coset key modulo the rows of T. *)
hoppingCosetKey[vector_List, rightSmith_List, diagonal_List] :=
  MapThread[Mod, {vector . rightSmith, diagonal}];

hoppingCellRepresentatives[cellMatrix_List] := Module[
  {
    leftSmith, smithDiagonal, rightSmith, diagonal,
    smithCoordinates, inverseRight, representatives, records
  },
  {leftSmith, smithDiagonal, rightSmith} = SmithDecomposition[cellMatrix];
  diagonal = Abs[Diagonal[smithDiagonal]];
  If[
    Times @@ diagonal =!= Abs[Det[cellMatrix]] ||
      !And @@ (IntegerQ[#] && # > 0 & /@ diagonal),
    Return[$Failed]
  ];
  smithCoordinates = Tuples[(Range[0, # - 1] &) /@ diagonal];
  inverseRight = Inverse[rightSmith];
  representatives = Simplify[# . inverseRight] & /@ smithCoordinates;
  If[!And @@ (hoppingIntegerVectorQ /@ representatives),
    Return[$Failed]
  ];
  records = SortBy[
    MapThread[
      Function[{representative, smithCoordinate},
        <|
          "Representative" -> representative,
          "FractionalCoordinate" -> Mod[
            Simplify[representative . Inverse[cellMatrix]],
            1
          ],
          "CosetKey" -> hoppingCosetKey[
            representative,
            rightSmith,
            diagonal
          ],
          "SmithCoordinate" -> smithCoordinate
        |>
      ],
      {representatives, smithCoordinates}
    ],
    Lookup[#, "FractionalCoordinate"] &
  ];
  <|
    "Records" -> records,
    "Representatives" -> Lookup[records, "Representative"],
    "FractionalCoordinates" -> Lookup[records, "FractionalCoordinate"],
    "CosetKeys" -> Lookup[records, "CosetKey"],
    "RightSmithMatrix" -> rightSmith,
    "SmithDiagonal" -> diagonal
  |>
];

hoppingEmbedBlock[
    matrix_, rowBlock_Integer, columnBlock_Integer,
    blockDimension_Integer, blockCount_Integer
  ] := Module[{rowOffset, columnOffset, rules},
  rowOffset = (rowBlock - 1) blockDimension;
  columnOffset = (columnBlock - 1) blockDimension;
  rules = hoppingSparseRules[matrix] /.
    Rule[{row_Integer, column_Integer}, value_] :>
      Rule[{rowOffset + row, columnOffset + column}, value];
  SparseArray[
    rules,
    {blockCount blockDimension, blockCount blockDimension}
  ]
];

hoppingOptionNames[function_] := First /@ Options[function];
hoppingUnknownOptions[function_, supplied_List] := Complement[
  DeleteDuplicates[
    First /@ Cases[
      supplied,
      _Rule | _RuleDelayed,
      {1}
    ]
  ],
  hoppingOptionNames[function]
];

hoppingData::selection =
  "The shell selection must be a positive integer n (meaning shells 1 through n) or a nonempty list of positive shell indices; received `1`.";
hoppingData::shellrange =
  "Shell indices `1` were not prepared by init; only `2` shells are available.";
hoppingData::rules =
  "Parameter values must be supplied as a list of rules or an Association; received `1`.";
hoppingData::unsolved =
  "No cached real-space result exists for shell(s) `1` with Hermitian -> `2`, KernelMethod -> `3`, and ValidationLevel -> `4`. Evaluate symham for those shells with the same options first.";
hoppingData::numeric =
  "The real-space hopping matrices remain nonnumeric after applying the supplied parameter rules. Remaining symbols include `1`.";
hoppingData::data =
  "The cached real-space hopping data for shell `1` is missing or inconsistent with the static bond data.";
hoppingData::option = "Invalid hoppingData option value(s): `1`.";

Options[hoppingData] = {
  "Hermitian" -> True,
  "KernelMethod" -> "Iterative",
  "ValidationLevel" -> "Basic"
};

hoppingData[
    selection_, rules_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, session, availableShells,
    shells, parameterRules, hermitianQ, method, validationLevel,
    data, model, lattice, centers
  },
  unknown = hoppingUnknownOptions[hoppingData, supplied];
  If[unknown =!= {}, Message[hoppingData::option, unknown]; Return[$Failed]];
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  availableShells = Length[session["CompiledBondShells"]];
  shells = normalizeWannierShellSelection[selection, availableShells];
  Which[
    MatchQ[shells, Failure["InvalidShellSelection", _Association]],
      Message[hoppingData::selection, selection];
      Return[$Failed],
    MatchQ[shells, Failure["ShellOutOfRange", _Association]],
      Message[hoppingData::shellrange, shells[[2, "Shells"]], availableShells];
      Return[$Failed]
  ];
  parameterRules = normalizeWannierParameterRules[rules];
  If[parameterRules === $Failed,
    Message[hoppingData::rules, rules];
    Return[$Failed]
  ];
  hermitianQ = OptionValue["Hermitian"];
  method = OptionValue["KernelMethod"];
  validationLevel = OptionValue["ValidationLevel"];
  If[
    !MemberQ[{True, False}, hermitianQ] || !StringQ[method] ||
      !StringQ[validationLevel],
    Message[hoppingData::option, {hermitianQ, method, validationLevel}];
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
        hoppingData::unsolved,
        data[[2, "Shells"]],
        hermitianQ,
        method,
        validationLevel
      ];
      Return[$Failed],
    MatchQ[data, Failure["NonnumericHopping", _Association]],
      Message[hoppingData::numeric, data[[2, "Symbols"]]];
      Return[$Failed],
    MatchQ[data, Failure["InvalidShellData", _Association]],
      Message[hoppingData::data, data[[2, "Shell"]]];
      Return[$Failed],
    data === $Failed || FailureQ[data],
      Message[hoppingData::data, shells];
      Return[$Failed]
  ];
  model = Lookup[session, "ModelSpecification", <||>];
  lattice = Lookup[model, "Lattice", None];
  centers = If[
    ListQ[Lookup[model, "SiteOrbits", $Failed]] &&
      ListQ[Lookup[model, "LocalDimensions", $Failed]],
    compatibilityWannierCenters[
      model["SiteOrbits"],
      model["LocalDimensions"]
    ],
    None
  ];
  Join[
    data,
    <|
      "Source" -> "MagneticTBRealSpaceShellCache",
      "Convention" ->
        "H(R)[a,b]=<0,a|H|R,b>; H(k)=Sum_R H(R) Exp[2 Pi I k.R]"
    |>,
    If[hoppingRealMatrixQ[lattice, {3, 3}], <|"Lattice" -> lattice|>, <||>],
    If[
      ListQ[centers] && Length[centers] === data["NumWannier"] &&
        And @@ (hoppingRealVectorQ[#, 3] & /@ centers),
      <|"WannierCenters" -> centers|>,
      <||>
    ]
  ]
];

transformHoppings::data =
  "Expected complete numeric hopping data returned by hoppingData/readHR, or a valid wannier90_hr.dat path.";
transformHoppings::cell =
  "The cell transformation must be a nonsingular 3 by 3 integer matrix whose rows are the new direct-lattice vectors in the old fractional basis; received `1`.";
transformHoppings::centers =
  "WannierCenters must be Automatic, None, or one real fractional three-vector per input orbital.";
transformHoppings::lattice =
  "Lattice must be Automatic, None, or a real numeric 3 by 3 row-vector lattice matrix.";
transformHoppings::hermitian =
  "The hopping table is not Hermitian: H(R)=H(-R)^dagger residual `1` exceeds tolerance `2`.";
transformHoppings::internal =
  "The exact integer-cell transformation failed at old translation `1` and representative `2`.";
transformHoppings::option = "Invalid transformHoppings option value(s): `1`.";

Options[transformHoppings] = {
  "WannierCenters" -> Automatic,
  "Lattice" -> Automatic,
  "HermiticityTolerance" -> 10^-9
};

transformHoppings[
    input_, cellMatrix_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, data, tolerance,
    centerSetting, latticeSetting, inputCenters, inputLattice,
    hermiticityResidual, representativeData, representatives,
    representativeCoordinates, rightSmith, smithDiagonal,
    keyToIndex, translations, matrices, oldDimension,
    representativeCount, transformedDimension, inverseCell,
    recordsOrFailure, recordGroups, transformedTranslations,
    transformedMatrices, transformedCenters, transformedLattice,
    orbitalMap, source
  },
  unknown = hoppingUnknownOptions[transformHoppings, supplied];
  If[unknown =!= {},
    Message[transformHoppings::option, unknown];
    Return[$Failed]
  ];
  tolerance = OptionValue["HermiticityTolerance"];
  centerSetting = OptionValue["WannierCenters"];
  latticeSetting = OptionValue["Lattice"];
  If[
    !hoppingFiniteNumberQ[tolerance] ||
      !TrueQ[Im[N[tolerance]] == 0] || !TrueQ[tolerance >= 0],
    Message[transformHoppings::option, tolerance];
    Return[$Failed]
  ];
  If[!hoppingCellMatrixQ[cellMatrix],
    Message[transformHoppings::cell, cellMatrix];
    Return[$Failed]
  ];
  data = hoppingResolveData[input];
  If[data === $Failed || !hoppingDataQ[data],
    Message[transformHoppings::data];
    Return[$Failed]
  ];
  hermiticityResidual = hoppingHermiticityResidual[data];
  If[TrueQ[hermiticityResidual > tolerance],
    Message[transformHoppings::hermitian, hermiticityResidual, tolerance];
    Return[$Failed]
  ];
  oldDimension = data["NumWannier"];
  inputCenters = Which[
    centerSetting === None, None,
    centerSetting === Automatic, Lookup[data, "WannierCenters", None],
    True, centerSetting
  ];
  If[
    inputCenters =!= None &&
      (!ListQ[inputCenters] || Length[inputCenters] =!= oldDimension ||
        !And @@ (hoppingRealVectorQ[#, 3] & /@ inputCenters)),
    Message[transformHoppings::centers];
    Return[$Failed]
  ];
  inputLattice = Which[
    latticeSetting === None, None,
    latticeSetting === Automatic, Lookup[data, "Lattice", None],
    True, latticeSetting
  ];
  If[
    inputLattice =!= None && !hoppingRealMatrixQ[inputLattice, {3, 3}],
    Message[transformHoppings::lattice];
    Return[$Failed]
  ];

  representativeData = hoppingCellRepresentatives[cellMatrix];
  If[representativeData === $Failed,
    Message[transformHoppings::cell, cellMatrix];
    Return[$Failed]
  ];
  representatives = representativeData["Representatives"];
  representativeCoordinates = representativeData["FractionalCoordinates"];
  rightSmith = representativeData["RightSmithMatrix"];
  smithDiagonal = representativeData["SmithDiagonal"];
  representativeCount = Length[representatives];
  transformedDimension = representativeCount oldDimension;
  inverseCell = Inverse[cellMatrix];
  keyToIndex = AssociationThread[
    hoppingTranslationKey /@ representativeData["CosetKeys"],
    Range[representativeCount]
  ];
  translations = data["Translations"];
  matrices = hoppingEffectiveMatrices[data];

  recordsOrFailure = Catch[
    Module[{harvested},
      harvested = Reap[
        Do[
          If[hoppingMatrixNonzeroQ[matrices[[translationIndex]]],
            Do[
              Module[
                {
                  targetOldCell, cosetKey, targetRepresentative,
                  transformedTranslation, embedded
                },
                targetOldCell =
                  representatives[[representativeIndex]] +
                    translations[[translationIndex]];
                cosetKey = hoppingCosetKey[
                  targetOldCell,
                  rightSmith,
                  smithDiagonal
                ];
                targetRepresentative = Lookup[
                  keyToIndex,
                  hoppingTranslationKey[cosetKey],
                  Missing["NotFound"]
                ];
                If[MissingQ[targetRepresentative],
                  Throw[
                    Failure[
                      "MissingCoset",
                      <|
                        "Translation" -> translations[[translationIndex]],
                        "Representative" -> representativeIndex
                      |>
                    ],
                    "TransformHoppingsFailure"
                  ]
                ];
                transformedTranslation = Simplify[
                  (targetOldCell -
                    representatives[[targetRepresentative]]) . inverseCell
                ];
                If[!hoppingIntegerVectorQ[transformedTranslation],
                  Throw[
                    Failure[
                      "NonintegerTranslation",
                      <|
                        "Translation" -> translations[[translationIndex]],
                        "Representative" -> representativeIndex
                      |>
                    ],
                    "TransformHoppingsFailure"
                  ]
                ];
                embedded = hoppingEmbedBlock[
                  matrices[[translationIndex]],
                  representativeIndex,
                  targetRepresentative,
                  oldDimension,
                  representativeCount
                ];
                Sow[{transformedTranslation, embedded}]
              ],
              {representativeIndex, representativeCount}
            ]
          ],
          {translationIndex, Length[translations]}
        ]
      ][[2]];
      If[harvested === {}, {}, First[harvested]]
    ],
    "TransformHoppingsFailure"
  ];
  If[FailureQ[recordsOrFailure],
    Message[
      transformHoppings::internal,
      Lookup[recordsOrFailure[[2]], "Translation", Missing["NotAvailable"]],
      Lookup[recordsOrFailure[[2]], "Representative", Missing["NotAvailable"]]
    ];
    Return[$Failed]
  ];
  If[recordsOrFailure === {},
    transformedTranslations = {{0, 0, 0}};
    transformedMatrices = {
      SparseArray[{}, {transformedDimension, transformedDimension}]
    },
    recordGroups = SortBy[
      GatherBy[recordsOrFailure, First],
      First@First[#] &
    ];
    transformedTranslations = First@First[#] & /@ recordGroups;
    transformedMatrices = Total[Last /@ #] & /@ recordGroups
  ];
  transformedCenters = If[
    inputCenters === None,
    None,
    Flatten[
      Table[
        Mod[
          Simplify[(representative + center) . inverseCell],
          1
        ],
        {representative, representatives},
        {center, inputCenters}
      ],
      1
    ]
  ];
  transformedLattice = If[
    inputLattice === None,
    None,
    cellMatrix . inputLattice
  ];
  orbitalMap = Flatten[
    Table[
      <|
        "NewOrbitalIndex" ->
          (representativeIndex - 1) oldDimension + orbital,
        "RepresentativeIndex" -> representativeIndex,
        "OriginalOrbitalIndex" -> orbital,
        "OldCellTranslation" -> representatives[[representativeIndex]],
        "NewCellFractionalCoordinate" ->
          representativeCoordinates[[representativeIndex]]
      |>,
      {representativeIndex, representativeCount},
      {orbital, oldDimension}
    ],
    1
  ];
  source = Lookup[data, "Source", Missing["NotAvailable"]];
  Join[
    <|
      "Schema" -> "MagneticTBHoppingData",
      "SchemaVersion" -> 1,
      "Source" -> source,
      "SourceSchema" -> Lookup[data, "Schema", Missing["NotAvailable"]],
      "NumWannier" -> transformedDimension,
      "NumTranslations" -> Length[transformedTranslations],
      "Translations" -> transformedTranslations,
      "Degeneracies" -> ConstantArray[1, Length[transformedTranslations]],
      "HoppingMatrices" -> transformedMatrices,
      "CellMatrix" -> cellMatrix,
      "CellVolumeFactor" -> Abs[Det[cellMatrix]],
      "CellRepresentatives" -> representatives,
      "RepresentativeCoordinates" -> representativeCoordinates,
      "OrbitalMap" -> orbitalMap,
      "HermitianResidual" -> hoppingHermiticityResidual[
        <|
          "NumWannier" -> transformedDimension,
          "Translations" -> transformedTranslations,
          "Degeneracies" -> ConstantArray[1, Length[transformedTranslations]],
          "HoppingMatrices" -> transformedMatrices
        |>
      ],
      "Convention" ->
        "H(R)[a,b]=<0,a|H|R,b>; new lattice rows equal CellMatrix.old lattice rows"
    |>,
    If[transformedCenters === None, <||>,
      <|"WannierCenters" -> transformedCenters|>
    ],
    If[transformedLattice === None, <||>,
      <|"Lattice" -> transformedLattice|>
    ]
  ]
];

finiteRealSpaceBasisData[cells_List, data_Association] := Module[
  {
    orbitalDimension, centers, lattice, positionDataAvailable,
    records
  },
  orbitalDimension = data["NumWannier"];
  centers = Lookup[data, "WannierCenters", None];
  lattice = Lookup[data, "Lattice", None];
  positionDataAvailable =
    ListQ[centers] && Length[centers] === orbitalDimension &&
      And @@ (hoppingRealVectorQ[#, 3] & /@ centers) &&
      hoppingRealMatrixQ[lattice, {3, 3}];
  records = Flatten[
    Table[
      Module[{fractionalPosition, cartesianPosition},
        fractionalPosition = If[
          positionDataAvailable,
          N[cells[[cellIndex]] + centers[[orbitalIndex]]],
          Missing["NotAvailable"]
        ];
        cartesianPosition = If[
          positionDataAvailable,
          N[fractionalPosition . lattice],
          Missing["NotAvailable"]
        ];
        <|
          "BasisIndex" ->
            (cellIndex - 1) orbitalDimension + orbitalIndex,
          "CellIndex" -> cellIndex,
          "Cell" -> cells[[cellIndex]],
          "OrbitalIndex" -> orbitalIndex,
          "WannierCenter" -> If[
            positionDataAvailable,
            centers[[orbitalIndex]],
            Missing["NotAvailable"]
          ],
          "FractionalPosition" -> fractionalPosition,
          "CartesianPosition" -> cartesianPosition
        |>
      ],
      {cellIndex, Length[cells]},
      {orbitalIndex, orbitalDimension}
    ],
    1
  ];
  Join[
    <|
      "PositionDataAvailable" -> positionDataAvailable,
      "BasisRecords" -> records
    |>,
    If[positionDataAvailable,
      <|"Lattice" -> lattice, "WannierCenters" -> centers|>,
      <||>
    ]
  ]
];

buildRealSpaceHamiltonian::data = transformHoppings::data;
buildRealSpaceHamiltonian::size =
  "The finite geometry must be either three positive sizes or a nonempty duplicate-free list of integer cell coordinates {i,j,k}; received `1`.";
buildRealSpaceHamiltonian::boundary =
  "BoundaryConditions must contain three entries, each \"Open\" or \"Periodic\". An explicit arbitrary cell set supports only open boundaries; received `1`.";
buildRealSpaceHamiltonian::hermitian = transformHoppings::hermitian;
buildRealSpaceHamiltonian::option =
  "Invalid buildRealSpaceHamiltonian option value(s): `1`.";

Options[buildRealSpaceHamiltonian] = {
  "BoundaryConditions" -> {"Open", "Open", "Open"},
  "HermiticityTolerance" -> 10^-9,
  "Output" -> "Matrix"
};

buildRealSpaceHamiltonian[
    input_, geometry_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, data, boundary, tolerance,
    output, residual, rectangularQ, explicitCellsQ, size, shape,
    cells, cellIndex, matrices, translations,
    orbitalDimension, totalDimension, entries, groups, rules,
    hamiltonian, basisData
  },
  unknown = hoppingUnknownOptions[buildRealSpaceHamiltonian, supplied];
  If[unknown =!= {},
    Message[buildRealSpaceHamiltonian::option, unknown];
    Return[$Failed]
  ];
  boundary = OptionValue["BoundaryConditions"];
  tolerance = OptionValue["HermiticityTolerance"];
  output = OptionValue["Output"];
  rectangularQ = VectorQ[geometry, IntegerQ] &&
    Length[geometry] === 3 && Min[geometry] >= 1;
  explicitCellsQ = ListQ[geometry] && geometry =!= {} &&
    And @@ (hoppingIntegerVectorQ /@ geometry) &&
    DuplicateFreeQ[geometry];
  If[!rectangularQ && !explicitCellsQ,
    Message[buildRealSpaceHamiltonian::size, geometry];
    Return[$Failed]
  ];
  If[
    !ListQ[boundary] || Length[boundary] =!= 3 ||
      !And @@ (MemberQ[{"Open", "Periodic"}, #] & /@ boundary),
    Message[buildRealSpaceHamiltonian::boundary, boundary];
    Return[$Failed]
  ];
  If[explicitCellsQ && boundary =!= {"Open", "Open", "Open"},
    Message[buildRealSpaceHamiltonian::boundary, boundary];
    Return[$Failed]
  ];
  If[
    !hoppingFiniteNumberQ[tolerance] ||
      !TrueQ[Im[N[tolerance]] == 0] || !TrueQ[tolerance >= 0] ||
      !MemberQ[{"Matrix", "Data"}, output],
    Message[buildRealSpaceHamiltonian::option, {tolerance, output}];
    Return[$Failed]
  ];
  data = hoppingResolveData[input];
  If[data === $Failed || !hoppingDataQ[data],
    Message[buildRealSpaceHamiltonian::data];
    Return[$Failed]
  ];
  residual = hoppingHermiticityResidual[data];
  If[TrueQ[residual > tolerance],
    Message[buildRealSpaceHamiltonian::hermitian, residual, tolerance];
    Return[$Failed]
  ];
  If[rectangularQ,
    size = geometry;
    shape = "Rectangular";
    cells = Tuples[(Range[0, # - 1] &) /@ size],
    size = None;
    shape = "ExplicitCells";
    cells = geometry
  ];
  cellIndex = AssociationThread[
    hoppingTranslationKey /@ cells,
    Range[Length[cells]]
  ];
  matrices = hoppingEffectiveMatrices[data];
  translations = data["Translations"];
  orbitalDimension = data["NumWannier"];
  totalDimension = Length[cells] orbitalDimension;
  entries = Reap[
    Do[
      Module[{target, targetIndex, localRules},
        target = cells[[sourceIndex]] + translations[[translationIndex]];
        If[rectangularQ,
          Do[
            If[boundary[[axis]] === "Periodic",
              target[[axis]] = Mod[target[[axis]], size[[axis]]]
            ],
            {axis, 3}
          ]
        ];
        targetIndex = Lookup[
          cellIndex,
          hoppingTranslationKey[target],
          Missing["NotFound"]
        ];
        If[!MissingQ[targetIndex],
          localRules = hoppingSparseRules[matrices[[translationIndex]]] /.
            Rule[{row_Integer, column_Integer}, value_] :>
              {
                {
                  (sourceIndex - 1) orbitalDimension + row,
                  (targetIndex - 1) orbitalDimension + column
                },
                value
              };
          Scan[Sow, localRules]
        ]
      ],
      {sourceIndex, Length[cells]},
      {translationIndex, Length[translations]}
    ]
  ][[2]];
  entries = If[entries === {}, {}, First[entries]];
  groups = GatherBy[entries, First];
  rules = Rule[First@First[#], Total[#[[All, 2]]]] & /@ groups;
  rules = Select[rules, !TrueQ[PossibleZeroQ[Last[#]]] &];
  hamiltonian = SparseArray[rules, {totalDimension, totalDimension}];
  basisData = finiteRealSpaceBasisData[cells, data];
  If[
    output === "Matrix",
    hamiltonian,
    Join[
      <|
      "Schema" -> "MagneticTBFiniteRealSpaceHamiltonian",
      "SchemaVersion" -> 2,
      "Hamiltonian" -> hamiltonian,
      "Shape" -> shape,
      "Cells" -> cells,
      "NumCells" -> Length[cells],
      "OrbitalsPerCell" -> orbitalDimension,
      "Dimension" -> totalDimension,
      "BoundaryConditions" -> boundary,
      "CellBounds" -> (MinMax /@ Transpose[cells]),
      "HermitianResidual" -> Max@Abs@Flatten@N[
        hamiltonian - ConjugateTranspose[hamiltonian]
      ],
      "Convention" ->
        "matrix row=(source cell,orbital), column=(target cell,orbital)"
      |>,
      If[rectangularQ, <|"Size" -> size|>, <||>],
      basisData
    ]
  ]
];

buildSlabHamiltonian::data = transformHoppings::data;
buildSlabHamiltonian::size =
  "The hybrid-space size must be three positive integers, with size 1 along every periodic direction; received size `1` and periodic directions `2`.";
buildSlabHamiltonian::momentum =
  "The crystal momentum must be one finite real reciprocal-fractional coordinate per periodic direction; received `1` for directions `2`.";
buildSlabHamiltonian::directions =
  "PeriodicDirections must contain one or two distinct axis indices chosen from {1,2,3}; received `1`.";
buildSlabHamiltonian::cell =
  "CellMatrix must be Automatic or a nonsingular 3 by 3 integer matrix whose rows define the calculation axes; received `1`.";
buildSlabHamiltonian::hermitian = transformHoppings::hermitian;
buildSlabHamiltonian::option =
  "Invalid buildSlabHamiltonian option value(s): `1`.";

Options[buildSlabHamiltonian] = {
  "PeriodicDirections" -> {1, 2},
  "CellMatrix" -> Automatic,
  "HermiticityTolerance" -> 10^-9,
  "Output" -> "Matrix"
};

buildSlabHamiltonian[
    input_, size_, momentum_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, periodicDirections,
    openDirections, cellMatrix, tolerance, output, data, residual,
    coordinateRanges, cells, cellIndex, matrices, translations,
    embeddedMomentum, orbitalDimension, totalDimension, entries,
    groups, rules, hamiltonian
  },
  unknown = hoppingUnknownOptions[buildSlabHamiltonian, supplied];
  If[unknown =!= {},
    Message[buildSlabHamiltonian::option, unknown];
    Return[$Failed]
  ];
  periodicDirections = OptionValue["PeriodicDirections"];
  cellMatrix = OptionValue["CellMatrix"];
  tolerance = OptionValue["HermiticityTolerance"];
  output = OptionValue["Output"];
  If[
    !ListQ[periodicDirections] ||
      !MemberQ[{1, 2}, Length[periodicDirections]] ||
      !VectorQ[periodicDirections, IntegerQ] ||
      !DuplicateFreeQ[periodicDirections] ||
      !And @@ (MemberQ[{1, 2, 3}, #] & /@ periodicDirections),
    Message[buildSlabHamiltonian::directions, periodicDirections];
    Return[$Failed]
  ];
  If[
    !VectorQ[size, IntegerQ] || Length[size] =!= 3 || Min[size] < 1 ||
      !And @@ (# === 1 & /@ size[[periodicDirections]]),
    Message[
      buildSlabHamiltonian::size,
      size,
      periodicDirections
    ];
    Return[$Failed]
  ];
  If[!hoppingRealVectorQ[momentum, Length[periodicDirections]],
    Message[
      buildSlabHamiltonian::momentum,
      momentum,
      periodicDirections
    ];
    Return[$Failed]
  ];
  If[cellMatrix =!= Automatic && !hoppingCellMatrixQ[cellMatrix],
    Message[buildSlabHamiltonian::cell, cellMatrix];
    Return[$Failed]
  ];
  If[
    !hoppingFiniteNumberQ[tolerance] ||
      !TrueQ[Im[N[tolerance]] == 0] || !TrueQ[tolerance >= 0] ||
      !MemberQ[{"Matrix", "Data"}, output],
    Message[buildSlabHamiltonian::option, {tolerance, output}];
    Return[$Failed]
  ];
  data = hoppingResolveData[input];
  If[data === $Failed || !hoppingDataQ[data],
    Message[buildSlabHamiltonian::data];
    Return[$Failed]
  ];
  If[cellMatrix =!= Automatic,
    data = transformHoppings[
      data,
      cellMatrix,
      "HermiticityTolerance" -> tolerance
    ];
    If[data === $Failed, Return[$Failed]]
  ];
  residual = hoppingHermiticityResidual[data];
  If[TrueQ[residual > tolerance],
    Message[buildSlabHamiltonian::hermitian, residual, tolerance];
    Return[$Failed]
  ];
  openDirections = Complement[{1, 2, 3}, periodicDirections];
  coordinateRanges = Table[
    If[
      MemberQ[periodicDirections, axis],
      {0},
      Range[0, size[[axis]] - 1]
    ],
    {axis, 3}
  ];
  cells = Tuples[coordinateRanges];
  cellIndex = AssociationThread[
    hoppingTranslationKey /@ cells,
    Range[Length[cells]]
  ];
  matrices = hoppingEffectiveMatrices[data];
  translations = data["Translations"];
  embeddedMomentum = ConstantArray[0, 3];
  embeddedMomentum[[periodicDirections]] = momentum;
  orbitalDimension = data["NumWannier"];
  totalDimension = Length[cells] orbitalDimension;
  entries = Reap[
    Do[
      Module[{target, phase, targetIndex, localRules},
        target = cells[[sourceIndex]] + translations[[translationIndex]];
        phase = Exp[
          2 Pi I embeddedMomentum.translations[[translationIndex]]
        ];
        target[[periodicDirections]] =
          ConstantArray[0, Length[periodicDirections]];
        If[
          And @@ (
            Between[target[[#]], {0, size[[#]] - 1}] & /@
              openDirections
          ),
          targetIndex = Lookup[
            cellIndex,
            hoppingTranslationKey[target],
            Missing["NotFound"]
          ];
          If[!MissingQ[targetIndex],
            localRules = hoppingSparseRules[matrices[[translationIndex]]] /.
              Rule[{row_Integer, column_Integer}, value_] :>
                {
                  {
                    (sourceIndex - 1) orbitalDimension + row,
                    (targetIndex - 1) orbitalDimension + column
                  },
                  value phase
                };
            Scan[Sow, localRules]
          ]
        ]
      ],
      {sourceIndex, Length[cells]},
      {translationIndex, Length[translations]}
    ]
  ][[2]];
  entries = If[entries === {}, {}, First[entries]];
  groups = GatherBy[entries, First];
  rules = Rule[First@First[#], Total[#[[All, 2]]]] & /@ groups;
  rules = Select[rules, !TrueQ[PossibleZeroQ[Last[#]]] &];
  hamiltonian = SparseArray[rules, {totalDimension, totalDimension}];
  If[
    output === "Matrix",
    hamiltonian,
    <|
      "Schema" -> "MagneticTBHybridSpaceHamiltonian",
      "SchemaVersion" -> 1,
      "Hamiltonian" -> hamiltonian,
      "Cells" -> cells,
      "NumFiniteCells" -> Length[cells],
      "OrbitalsPerCell" -> orbitalDimension,
      "Dimension" -> totalDimension,
      "Size" -> size,
      "PeriodicDirections" -> periodicDirections,
      "OpenDirections" -> openDirections,
      "CrystalMomentum" -> momentum,
      "EmbeddedCrystalMomentum" -> embeddedMomentum,
      "CellMatrix" -> cellMatrix,
      "HermitianResidual" -> Max@Abs@Flatten@N[
        hamiltonian - ConjugateTranspose[hamiltonian]
      ],
      "Convention" ->
        "H(R)[a,b]=<0,a|H|R,b>; periodic hopping carries Exp[2 Pi I k.R], while nonperiodic directions use finite open cells"
    |>
  ]
];

buildBlochHamiltonian::data = transformHoppings::data;
buildBlochHamiltonian::momentum =
  "The crystal momentum must be a finite real three-vector in reciprocal fractional coordinates; received `1`.";
buildBlochHamiltonian::cell =
  "CellMatrix must be Automatic or a nonsingular 3 by 3 integer matrix whose rows define the supercell direct-lattice vectors; received `1`.";
buildBlochHamiltonian::hermitian = transformHoppings::hermitian;
buildBlochHamiltonian::option =
  "Invalid buildBlochHamiltonian option value(s): `1`.";

Options[buildBlochHamiltonian] = {
  "CellMatrix" -> Automatic,
  "HermiticityTolerance" -> 10^-9,
  "Output" -> "Matrix"
};

buildBlochHamiltonian[
    input_, momentum_, suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    supplied = {suppliedOptions}, unknown, cellMatrix, tolerance,
    output, data, hoppingResidual, hamiltonian, matrixResidual
  },
  unknown = hoppingUnknownOptions[buildBlochHamiltonian, supplied];
  If[unknown =!= {},
    Message[buildBlochHamiltonian::option, unknown];
    Return[$Failed]
  ];
  cellMatrix = OptionValue["CellMatrix"];
  tolerance = OptionValue["HermiticityTolerance"];
  output = OptionValue["Output"];
  If[!hoppingRealVectorQ[momentum, 3],
    Message[buildBlochHamiltonian::momentum, momentum];
    Return[$Failed]
  ];
  If[cellMatrix =!= Automatic && !hoppingCellMatrixQ[cellMatrix],
    Message[buildBlochHamiltonian::cell, cellMatrix];
    Return[$Failed]
  ];
  If[
    !hoppingFiniteNumberQ[tolerance] ||
      !TrueQ[Im[N[tolerance]] == 0] || !TrueQ[tolerance >= 0] ||
      !MemberQ[{"Matrix", "Data"}, output],
    Message[buildBlochHamiltonian::option, {tolerance, output}];
    Return[$Failed]
  ];
  data = hoppingResolveData[input];
  If[data === $Failed || !hoppingDataQ[data],
    Message[buildBlochHamiltonian::data];
    Return[$Failed]
  ];
  If[cellMatrix =!= Automatic,
    data = transformHoppings[
      data,
      cellMatrix,
      "HermiticityTolerance" -> tolerance
    ];
    If[data === $Failed, Return[$Failed]]
  ];
  hoppingResidual = hoppingHermiticityResidual[data];
  If[TrueQ[hoppingResidual > tolerance],
    Message[
      buildBlochHamiltonian::hermitian,
      hoppingResidual,
      tolerance
    ];
    Return[$Failed]
  ];
  hamiltonian = hoppingBlochHamiltonian[data, momentum];
  If[hamiltonian === $Failed,
    Message[buildBlochHamiltonian::data];
    Return[$Failed]
  ];
  matrixResidual = Max@Abs@Flatten@N[
    hamiltonian - ConjugateTranspose[hamiltonian]
  ];
  If[TrueQ[matrixResidual > tolerance],
    Message[
      buildBlochHamiltonian::hermitian,
      matrixResidual,
      tolerance
    ];
    Return[$Failed]
  ];
  If[
    output === "Matrix",
    hamiltonian,
    <|
      "Schema" -> "MagneticTBBlochHamiltonian",
      "SchemaVersion" -> 1,
      "Hamiltonian" -> hamiltonian,
      "Dimension" -> data["NumWannier"],
      "CrystalMomentum" -> momentum,
      "CellMatrix" -> cellMatrix,
      "CellVolumeFactor" -> Lookup[data, "CellVolumeFactor", 1],
      "CellRepresentatives" -> Lookup[
        data,
        "CellRepresentatives",
        {{0, 0, 0}}
      ],
      "HermitianResidual" -> matrixResidual,
      "Convention" ->
        "H(k)=Sum_R H(R) Exp[2 Pi I k.R], with k in reciprocal fractional coordinates of the active cell"
    |>
  ]
];

(* This private evaluator is shared by buildBlochHamiltonian and hopping
   transform tests. It is never the source of transform or finite/surface
   assembly. *)
hoppingBlochHamiltonian[input_, momentum_] := Module[{data, matrices},
  data = hoppingResolveData[input];
  If[
    data === $Failed || !hoppingDataQ[data] ||
      !ListQ[momentum] || Length[momentum] =!= 3,
    Return[$Failed]
  ];
  matrices = hoppingEffectiveMatrices[data];
  Total@MapThread[
    #1 Exp[2 Pi I momentum.#2] &,
    {matrices, data["Translations"]}
  ]
];

End[]

EndPackage[]
