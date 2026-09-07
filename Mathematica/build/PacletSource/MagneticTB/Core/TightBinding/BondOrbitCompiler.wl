(* ::Package:: *)

BeginPackage["MagneticTBLinearAlgebra`"]

CompileDirectedBondOrbits::usage =
  "CompileDirectedBondOrbits[prepared,n] compiles one shell; CompileDirectedBondOrbits[prepared,{n1,...}] validates the prepared model once and compiles several shells.";
Begin["`Private`"]

ClearAll[
  compilerIntegerVectorQ,
  bondIndexLookupKey,
  compileBondIndexLookup,
  canonicalSiteCoordinate,
  findBondIndex,
  transformedBondData,
  computeLocalSiteAction,
  localSiteAction,
  compiledContinuousSiteGenerators,
  compilePreparedBondContext,
  compileModelSitesFromContext,
  compileBondShellFromContext,
  compileBondAction,
  compileDirectedBondOrbitFromContext
];

compilePreparedBondContext[prepared_Association] := Module[
  {siteOrbits, localDimensions, symmetry, actionData, bondClasses,
   coordinateDimension, actionAccessor, context, sites,
   localSiteActionTable, continuousSiteGenerators},
  If[!AssociationQ[prepared] || Lookup[prepared, "Schema", None] =!= "PreparedModel" ||
      Lookup[prepared, "SchemaVersion", None] =!= 1 ||
      !TrueQ[Lookup[prepared, "Verified", False]], Return[$Failed]];
  siteOrbits = Lookup[prepared, "SiteOrbits", $Failed];
  localDimensions = Lookup[prepared, "LocalDimensions", $Failed];
  symmetry = Lookup[prepared, "Symmetry", $Failed];
  actionData = Lookup[prepared, "RepresentationActionData", $Failed];
  bondClasses = Lookup[prepared, "BondClasses", $Failed];
  actionAccessor = Quiet@
    RepresentationActionData`CompileSiteActionAccessor[actionData];
  If[!SitePermutationCompiler`ValidSiteOrbitsQ[siteOrbits] ||
      !ListQ[localDimensions] ||
      Length[localDimensions] =!= Length[siteOrbits] ||
      !And @@ (IntegerQ[#] && Positive[#] & /@ localDimensions) ||
      !AssociationQ[symmetry] ||
      actionAccessor === $Failed ||
      !ListQ[bondClasses], Return[$Failed]];
  coordinateDimension = Length[siteOrbits[[1, 1]]];
  If[
    Lookup[actionData, "SiteOrbits", $Failed] =!= siteOrbits ||
      Lookup[actionData, "LocalDimensions", $Failed] =!= localDimensions ||
      Lookup[actionData, "SpatialActions", $Failed] =!=
        Lookup[symmetry, "SpatialActions", $Failed] ||
      Lookup[actionData, "GroupAlgebra", $Failed] =!=
        Lookup[symmetry, "GroupAlgebra", $Failed] ||
      Lookup[actionData, "AntiunitaryFlags", $Failed] =!=
        Lookup[symmetry, "AntiunitaryFlags", $Failed] ||
      !And @@ Map[
        Function[shell,
          ListQ[shell] && shell =!= {} && And @@ Map[
            Function[entry,
              MatchQ[entry, {_?NumericQ, _Integer?NonNegative, _List}] &&
                entry[[2]] === Length[entry[[3]]] &&
                And @@ (MatchQ[#, {_List, _List}] &&
                    Length[#[[1]]] === coordinateDimension &&
                    Length[#[[2]]] === coordinateDimension & /@ entry[[3]])
            ], shell]
        ], bondClasses],
    Return[$Failed]
  ];
  context = <|
    "PreparedModel" -> prepared,
    "SiteOrbits" -> siteOrbits,
    "LocalDimensions" -> localDimensions,
    "BondClasses" -> bondClasses,
    "SiteActionAccessor" -> actionAccessor
  |>;
  sites = compileModelSitesFromContext[context];
  If[sites === $Failed, Return[$Failed]];
  context = Join[context, <|"Sites" -> sites|>];
  localSiteActionTable = Table[
    computeLocalSiteAction[context, siteIndex, operationIndex],
    {operationIndex, Length[symmetry["SpatialActions"]]},
    {siteIndex, Length[sites]}
  ];
  If[
    MemberQ[localSiteActionTable, $Failed, Infinity] ||
      MemberQ[localSiteActionTable, _Missing, Infinity],
    Return[$Failed]
  ];
  continuousSiteGenerators = compiledContinuousSiteGenerators[prepared];
  If[continuousSiteGenerators === $Failed, Return[$Failed]];
  Join[
    context,
    <|
      "LocalSiteActionTable" -> localSiteActionTable,
      "ContinuousSiteGenerators" -> continuousSiteGenerators
    |>
  ]
];
compilePreparedBondContext[_] := $Failed;

compileModelSitesFromContext[context_Association] := Module[
  {siteOrbits, localDimensions, rawSites, dimensions, offsets},
  siteOrbits = context["SiteOrbits"];
  localDimensions = context["LocalDimensions"];
  rawSites = Flatten[
    Table[
      {<|
        "OrbitIndex" -> orbitIndex,
        "EquivalentIndex" -> equivalentIndex,
        "Position" -> siteOrbits[[orbitIndex, equivalentIndex]],
        "Dimension" -> localDimensions[[orbitIndex]]
      |>},
      {orbitIndex, Length[siteOrbits]},
      {equivalentIndex, Length[siteOrbits[[orbitIndex]]]}
    ],
    2
  ];
  dimensions = Lookup[rawSites, "Dimension"];
  offsets = FoldList[Plus, 0, dimensions];
  MapIndexed[
    Function[{site, index},
      Join[
        site,
        <|
          "SiteIndex" -> First[index],
          "BlockRange" ->
            Range[
              offsets[[First[index]]] + 1,
              offsets[[First[index] + 1]]
            ]
        |>
      ]
    ],
    rawSites
  ]
];
compilerIntegerVectorQ[vector_List] :=
  MatrixPredicates`IntegerVectorQ[vector];

bondIndexLookupKey[
    rowSite_Integer,
    columnSite_Integer,
    displacement_List
  ] := With[
  {
    row = rowSite,
    column = columnSite,
    canonicalDisplacement = Simplify[displacement]
  },
  HoldComplete[row, column, canonicalDisplacement]
];

compileBondIndexLookup[bonds_List] := Association@Map[
  Function[bond,
    bondIndexLookupKey[
      bond["RowSite"],
      bond["ColumnSite"],
      bond["Displacement"]
    ] -> bond["BondIndex"]
  ],
  bonds
];

canonicalSiteCoordinate[sites_List, coordinate_List] := Module[{site},
  site = SelectFirst[
    sites,
    compilerIntegerVectorQ[
      Simplify[coordinate - Lookup[#, "Position"]]
    ] &,
    Missing["SiteNotFound", coordinate]
  ];
  If[
    MissingQ[site],
    site,
    Join[
      site,
      <|
        "CellTranslation" ->
          Simplify[coordinate - Lookup[site, "Position"]]
      |>
    ]
  ]
];

findBondIndex[
    bondIndexLookup_Association,
    rowSite_Integer,
    columnSite_Integer,
    displacement_List
  ] := Lookup[
  bondIndexLookup,
  bondIndexLookupKey[rowSite, columnSite, displacement],
  Missing["BondNotFound", {rowSite, columnSite, displacement}]
];

compileBondShellFromContext[
    context_Association,
    shell_Integer?Positive
  ] := Module[{bondClasses, sites, endpoints, bonds},
  bondClasses = context["BondClasses"];
  If[shell > Length[bondClasses],
    Message[CompileDirectedBondOrbits::shell, shell];
    Return[$Failed]
  ];
  sites = context["Sites"];
  endpoints = Flatten[bondClasses[[shell, All, 3]], 1];
  bonds = MapIndexed[
    Function[{pair, index},
      Module[{rowEndpoint, columnEndpoint, rowSite, columnSite},
        rowEndpoint = pair[[1]];
        columnEndpoint = pair[[2]];
        rowSite = canonicalSiteCoordinate[sites, rowEndpoint];
        columnSite = canonicalSiteCoordinate[sites, columnEndpoint];
        If[MissingQ[rowSite] || MissingQ[columnSite],
          Message[CompileDirectedBondOrbits::endpoint, pair];
          Return[$Failed, Module]
        ];
        <|
          "BondIndex" -> First[index],
          "RowSite" -> Lookup[rowSite, "SiteIndex"],
          "ColumnSite" -> Lookup[columnSite, "SiteIndex"],
          "RowCell" -> Lookup[rowSite, "CellTranslation"],
          "ColumnCell" -> Lookup[columnSite, "CellTranslation"],
          "Displacement" -> Simplify[columnEndpoint - rowEndpoint],
          "Endpoints" -> pair
        |>
      ]
    ],
    endpoints
  ];
  If[MemberQ[bonds, $Failed], Return[$Failed]];
  <|
    "PreparedModel" -> context["PreparedModel"],
    "SiteActionAccessor" -> context["SiteActionAccessor"],
    "LocalSiteActionTable" -> context["LocalSiteActionTable"],
    "Shell" -> shell,
    "Sites" -> sites,
    "BlockDimensions" -> Lookup[sites, "Dimension"],
    "Bonds" -> bonds
  |>
];
transformedBondData[
    compiled_Association,
    bond_Association,
    operationIndex_Integer
  ] := Module[
  {prepared, operation, rotation, transformedRow, transformedColumn},
  prepared = compiled["PreparedModel"];
  operation = prepared["Symmetry", "SpatialActions"][[operationIndex]];
  rotation = operation["Rotation"];
  transformedRow = localSiteAction[
    compiled,
    bond["RowSite"],
    operationIndex
  ];
  transformedColumn = localSiteAction[
    compiled,
    bond["ColumnSite"],
    operationIndex
  ];
  If[MissingQ[transformedRow] || MissingQ[transformedColumn],
    Return[Missing[
      "TransformedSiteNotFound",
      {operationIndex, bond}
    ]]
  ];
  <|
    "RowSite" -> transformedRow["ImageSite"],
    "ColumnSite" -> transformedColumn["ImageSite"],
    "Displacement" -> Simplify[rotation.Lookup[bond, "Displacement"]]
  |>
];

computeLocalSiteAction[
    compiled_Association,
    siteIndex_Integer,
    operationIndex_Integer
  ] := Module[
  {sites, actionAccessor, originalSite, orbitIndex,
   equivalentIndex, action, transformedSite},
  sites = compiled["Sites"];
  actionAccessor = compiled["SiteActionAccessor"];
  originalSite = sites[[siteIndex]];
  orbitIndex = originalSite["OrbitIndex"];
  equivalentIndex = originalSite["EquivalentIndex"];
  action = RepresentationActionData`RepresentationSiteAction[
    actionAccessor,
    orbitIndex,
    equivalentIndex,
    operationIndex
  ];
  If[action === $Failed, Return[Missing["InvalidActionData"]]];
  transformedSite = SelectFirst[
    sites,
    Lookup[#, "OrbitIndex"] == orbitIndex &&
      Lookup[#, "EquivalentIndex"] == action["ImageSiteIndex"] &,
    Missing[
      "ActionImageSiteNotFound",
      {orbitIndex, action["ImageSiteIndex"]}
    ]
  ];
  If[MissingQ[transformedSite], Return[transformedSite]];
  <|
    "ImageSite" -> transformedSite["SiteIndex"],
    "Matrix" -> action["Matrix"],
    "Antiunitary" -> action["Antiunitary"]
  |>
];

localSiteAction[
    compiled_Association,
    siteIndex_Integer,
    operationIndex_Integer
  ] := Module[{table},
  table = Lookup[compiled, "LocalSiteActionTable", Missing["NotCached"]];
  If[
    ListQ[table] &&
      Between[operationIndex, {1, Length[table]}] &&
      ListQ[table[[operationIndex]]] &&
      Between[siteIndex, {1, Length[table[[operationIndex]]]}],
    table[[operationIndex, siteIndex]],
    computeLocalSiteAction[compiled, siteIndex, operationIndex]
  ]
];

compiledContinuousSiteGenerators[prepared_Association] := Module[
  {continuousSymmetry, siteGenerators},
  continuousSymmetry = Lookup[prepared, "ContinuousSymmetry", None];
  If[continuousSymmetry === None, Return[None]];
  siteGenerators = Lookup[
    continuousSymmetry,
    "SiteGenerators",
    $Failed
  ];
  If[
    !ListQ[siteGenerators] ||
      Length[siteGenerators] =!= Length[prepared["SiteOrbits"]] ||
      !And @@ MapThread[
        Length[#1] === Length[#2] &,
        {siteGenerators, prepared["SiteOrbits"]}
      ],
    $Failed,
    Flatten[siteGenerators, 1]
  ]
];

compileBondAction[
    compiled_Association,
    bond_Association,
    operationIndex_Integer,
    target_String : "Same"
  ] := Module[{left, right, antiunitary},
  If[!MemberQ[{"Same", "Reverse"}, target] ||
      !And @@ (KeyExistsQ[bond, #] & /@ {"RowSite", "ColumnSite"}) ||
      !Between[operationIndex, {
        1, Length[compiled["PreparedModel", "Symmetry", "SpatialActions"]]
      }],
    Return[$Failed]
  ];
  left = localSiteAction[compiled, bond["RowSite"], operationIndex];
  right = localSiteAction[compiled, bond["ColumnSite"], operationIndex];
  If[MissingQ[left] || MissingQ[right], Return[$Failed]];
  If[left["Antiunitary"] =!= right["Antiunitary"], Return[$Failed]];
  antiunitary = left["Antiunitary"];
  <|
    "OperationIndex" -> operationIndex,
    "Left" -> left["Matrix"],
    "Right" -> right["Matrix"],
    "Antiunitary" -> antiunitary,
    "Target" -> target
  |>
];
compileBondAction[___] := $Failed;

CompileDirectedBondOrbits::model = "Expected a validated PreparedModel.";
CompileDirectedBondOrbits::shell =
  "Bond shell `1` is not present in the PreparedModel.";
CompileDirectedBondOrbits::endpoint =
  "Could not map bond endpoint `1` to a unit-cell site.";
CompileDirectedBondOrbits::image =
  "Operation `1` maps bond `2` outside the selected shell.";
CompileDirectedBondOrbits::generators =
  "The same-bond stabilizer generator calculation failed for operation indices `1`; bond-orbit compilation stopped.";
CompileDirectedBondOrbits::coset =
  "The reverse-bond stabilizer is not the expected right coset for representative bond `1`; bond-orbit compilation stopped.";
CompileDirectedBondOrbits::action =
  "A bond action could not be compiled for representative bond `1`; bond-orbit compilation stopped.";
compileDirectedBondOrbitFromContext[
    context_Association,
    shell_Integer?Positive
  ] := Module[
  {prepared, compiled, bonds, bondIndexLookup, sites, symmetryIndices,
   group, imageTable,
   reverseIndices,
   unseen, directedOrbits, representativeIndex, directPairs, directIndices,
   memberRecords, representativeBond, dimensions,
   sameStabilizerIndices, sameStabilizerOperations,
   stabilizerGeneratorIndices,
   sameGeneratorOperations, generatorReductionCache, generatorCacheKey,
   bondToDirectedOrbit, continuousSiteGenerators, groupAction},
  prepared = context["PreparedModel"];
  compiled = compileBondShellFromContext[context, shell];
  If[compiled === $Failed, Return[$Failed]];
  continuousSiteGenerators = context["ContinuousSiteGenerators"];
  bonds = compiled["Bonds"];
  bondIndexLookup = compileBondIndexLookup[bonds];
  If[Length[bondIndexLookup] =!= Length[bonds], Return[$Failed]];
  sites = compiled["Sites"];
  group = prepared["Symmetry", "GroupAlgebra"];
  symmetryIndices = Range[Length[prepared["Symmetry", "SpatialActions"]]];
  imageTable = Table[
    Module[{transformed, imageIndex},
      transformed = transformedBondData[compiled, bond, operationIndex];
      If[MissingQ[transformed], Return[$Failed, Module]];
      imageIndex = findBondIndex[
        bondIndexLookup,
        transformed["RowSite"],
        transformed["ColumnSite"],
        transformed["Displacement"]
      ];
      If[MissingQ[imageIndex],
        Message[
          CompileDirectedBondOrbits::image,
          operationIndex,
          bond["BondIndex"]
        ];
        $Failed,
        imageIndex
      ]
    ],
    {operationIndex, symmetryIndices},
    {bond, bonds}
  ];
  If[MemberQ[imageTable, $Failed], Return[$Failed]];
  groupAction = Quiet@GroupAlgebra`CompileGroupAction[group, imageTable];
  If[groupAction === $Failed,
    Return[$Failed]
  ];
  reverseIndices = Table[
    findBondIndex[
      bondIndexLookup,
      bond["ColumnSite"],
      bond["RowSite"],
      Simplify[-bond["Displacement"]]
    ],
    {bond, bonds}
  ];
  If[AnyTrue[reverseIndices, MissingQ], Return[$Failed]];
  unseen = Range[Length[bonds]];
  directedOrbits = {};
  generatorReductionCache = <||>;
  While[unseen =!= {},
    representativeIndex = First[unseen];
    representativeBond = bonds[[representativeIndex]];
    directIndices = GroupAlgebra`ActionOrbitIndices[
      groupAction, representativeIndex];
    If[directIndices === $Failed, Return[$Failed]];
    directPairs = Table[
      {targetIndex, First@GroupAlgebra`TransporterIndices[
        groupAction, representativeIndex, targetIndex]},
      {targetIndex, directIndices}
    ];
    memberRecords = Map[
      Function[pair,
        <|
          "BondIndex" -> pair[[1]],
          "Mode" -> "Direct",
          "Operation" -> compileBondAction[
            compiled,
            representativeBond,
            pair[[2]],
            "Same"
          ]
        |>
      ],
      directPairs
    ];
    If[MemberQ[memberRecords, $Failed, Infinity],
      Message[CompileDirectedBondOrbits::action, representativeIndex];
      Return[$Failed]
    ];
    sameStabilizerIndices = GroupAlgebra`StabilizerIndices[
      groupAction, representativeIndex];
    If[sameStabilizerIndices === $Failed, Return[$Failed]];
    generatorCacheKey = ToString[sameStabilizerIndices, InputForm];
    stabilizerGeneratorIndices = Lookup[
      generatorReductionCache,
      generatorCacheKey,
      Missing["NotCached"]
    ];
    If[MissingQ[stabilizerGeneratorIndices],
      stabilizerGeneratorIndices =
        GroupAlgebra`FindGeneratorIndices[
          group,
          sameStabilizerIndices
        ];
      If[stabilizerGeneratorIndices === $Failed,
        Message[CompileDirectedBondOrbits::generators, sameStabilizerIndices];
        Return[$Failed]
      ];
      AssociateTo[
        generatorReductionCache,
        generatorCacheKey -> stabilizerGeneratorIndices
      ]
    ];
    sameStabilizerOperations = compileBondAction[
      compiled,
      representativeBond,
      #,
      "Same"
    ] & /@ sameStabilizerIndices;
    sameGeneratorOperations = compileBondAction[
      compiled,
      representativeBond,
      #,
      "Same"
    ] & /@ stabilizerGeneratorIndices;
    If[
      MemberQ[sameStabilizerOperations, $Failed, Infinity] ||
        MemberQ[sameGeneratorOperations, $Failed, Infinity],
      Message[CompileDirectedBondOrbits::action, representativeIndex];
      Return[$Failed]
    ];
    dimensions = {
      sites[[representativeBond["RowSite"], "Dimension"]],
      sites[[representativeBond["ColumnSite"], "Dimension"]]
    };
    AppendTo[
      directedOrbits,
      <|
        "DirectedOrbitIndex" -> Length[directedOrbits] + 1,
        "RepresentativeBondIndex" -> representativeIndex,
        "RepresentativeBond" -> representativeBond,
        "Dimensions" -> dimensions,
        "DirectMemberRecords" ->
          SortBy[memberRecords, Lookup[#, "BondIndex"] &],
        "SameStabilizerOperationIndices" -> sameStabilizerIndices,
        "SameStabilizerOperations" -> sameStabilizerOperations,
        "StabilizerGeneratorIndices" -> stabilizerGeneratorIndices,
        "SameGeneratorOperations" -> sameGeneratorOperations,
        "SameGeneratorReductionVerified" -> True
      |>
    ];
    unseen = Complement[unseen, directIndices]
  ];
  bondToDirectedOrbit = ConstantArray[0, Length[bonds]];
  Do[
    Do[
      bondToDirectedOrbit[[member["BondIndex"]]] = orbitIndex,
      {member, directedOrbits[[orbitIndex, "DirectMemberRecords"]]}
    ],
    {orbitIndex, Length[directedOrbits]}
  ];
  If[MemberQ[bondToDirectedOrbit, 0], Return[$Failed]];
  directedOrbits = Map[
    Function[orbit,
      Module[
        {representativeReverse, reverseOrbitIndex,
         reverseStabilizerIndices, reverseRepresentativeOperationIndex,
         reverseCosetIndices, reverseCosetVerified,
         reverseStabilizerOperations},
        representativeIndex = orbit["RepresentativeBondIndex"];
        representativeBond = orbit["RepresentativeBond"];
        representativeReverse = reverseIndices[[representativeIndex]];
        reverseOrbitIndex = bondToDirectedOrbit[[representativeReverse]];
        reverseStabilizerIndices = If[
          reverseOrbitIndex == orbit["DirectedOrbitIndex"] &&
            representativeReverse =!= representativeIndex,
          GroupAlgebra`TransporterIndices[
            groupAction, representativeIndex, representativeReverse],
          {}
        ];
        reverseRepresentativeOperationIndex = If[
          reverseStabilizerIndices === {},
          Missing["NotApplicable"],
          First[reverseStabilizerIndices]
        ];
        reverseCosetIndices = If[
          MissingQ[reverseRepresentativeOperationIndex],
          {},
          GroupAlgebra`RightCosetIndices[
            group,
            orbit["SameStabilizerOperationIndices"],
            reverseRepresentativeOperationIndex
          ]
        ];
        reverseCosetVerified = If[
          reverseStabilizerIndices === {},
          True,
          reverseCosetIndices =!= $Failed &&
            Sort[reverseCosetIndices] === Sort[reverseStabilizerIndices]
        ];
        If[!TrueQ[reverseCosetVerified],
          Message[CompileDirectedBondOrbits::coset, representativeIndex];
          Return[$Failed, Module]
        ];
        reverseStabilizerOperations = compileBondAction[
          compiled,
          representativeBond,
          #,
          "Same"
        ] & /@ reverseStabilizerIndices;
        If[MemberQ[reverseStabilizerOperations, $Failed, Infinity],
          Message[CompileDirectedBondOrbits::action, representativeIndex];
          Return[$Failed, Module]
        ];
        Join[
          orbit,
          <|
            "ReverseBondIndex" -> representativeReverse,
            "ReverseOrbitIndex" -> reverseOrbitIndex,
            "SelfReverseOrbit" ->
              TrueQ[reverseOrbitIndex == orbit["DirectedOrbitIndex"]],
            "ReverseStabilizerOperationIndices" ->
              reverseStabilizerIndices,
            "ReverseStabilizerOperations" ->
              reverseStabilizerOperations,
            "ReverseRepresentativeOperationIndex" ->
              reverseRepresentativeOperationIndex,
            "ReverseCosetOperationIndices" -> reverseCosetIndices,
            "ReverseCosetVerified" -> reverseCosetVerified
          |>
        ]
      ]
    ],
    directedOrbits
  ];
  If[MemberQ[directedOrbits, $Failed], Return[$Failed]];
  <|
    "Schema" -> "DirectedBondOrbitData",
    "SchemaVersion" -> 1,
    "Shell" -> shell,
    "Sites" -> sites,
    "BlockDimensions" -> compiled["BlockDimensions"],
    "ContinuousSiteGenerators" -> continuousSiteGenerators,
    "Bonds" -> bonds,
    "SymmetryIndices" -> symmetryIndices,
    "ImageTable" -> imageTable,
    "ReverseIndices" -> reverseIndices,
    "BondToDirectedOrbit" -> bondToDirectedOrbit,
    "DirectedOrbits" -> directedOrbits
  |>
];
CompileDirectedBondOrbits[
    prepared_Association,
    shell_Integer?Positive
  ] := Module[{context},
  context = compilePreparedBondContext[prepared];
  If[context === $Failed,
    Message[CompileDirectedBondOrbits::model];
    Return[$Failed]
  ];
  compileDirectedBondOrbitFromContext[context, shell]
];
CompileDirectedBondOrbits[
    prepared_Association,
    shells : {__Integer?Positive}
  ] := Module[{context, results},
  context = compilePreparedBondContext[prepared];
  If[context === $Failed,
    Message[CompileDirectedBondOrbits::model];
    Return[$Failed]
  ];
  results = compileDirectedBondOrbitFromContext[context, #] & /@ shells;
  If[MemberQ[results, $Failed], $Failed, results]
];
CompileDirectedBondOrbits[___] :=
  (Message[CompileDirectedBondOrbits::model]; $Failed);

End[]

EndPackage[]
