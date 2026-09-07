(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  exactHoppingParameterCoefficient,
  makeHoppingParameterElementRecord,
  resolveHoppingParameterShellData,
  representativeHoppingsFromTerms,
  buildHoppingParameterSourceData,
  buildHoppingParameterOccurrenceData,
  renderHoppingEndpoint,
  renderHoppingParameterRecords,
  showHoppingParametersImpl
];

showHoppingParameters::shell =
  "The shell index must be a positive integer; received `1`.";
showHoppingParameters::shellrange =
  "Bond shell `1` was not prepared by init; only `2` shells are available.";
showHoppingParameters::unsolved =
  "No cached real-space symham result exists for shell `1` with Hermitian -> `2`, KernelMethod -> `3`, and ValidationLevel -> `4`. Evaluate symham for that shell with the same options first.";
showHoppingParameters::parameter =
  "Parameter `1` is not one of the cached shell-`2` parameters `3`.";
showHoppingParameters::option =
  "Invalid option value(s): `1`.";
showHoppingParameters::data =
  "The cached parameter/hopping data for shell `1` is missing or inconsistent.";

exactHoppingParameterCoefficient[expression_, parameter_] := Module[
  {coefficient},
  coefficient = Expand@Coefficient[Expand[expression], parameter];
  If[
    MatrixPredicates`ExactZeroExpressionQ[coefficient],
    Missing["ZeroCoefficient"],
    coefficient
  ]
];

makeHoppingParameterElementRecord[
    shell_Integer?Positive,
    constraintOrbitIndex_Integer?Positive,
    representativeBondIndex_Integer?Positive,
    bondIndex_Integer?Positive,
    parameter_,
    coefficient_,
    rowLocalIndex_Integer?Positive,
    columnLocalIndex_Integer?Positive,
    bond_Association,
    sites_List,
    basisOrderData_Association
  ] := Module[
  {
    rowSiteIndex, columnSiteIndex, rowSite, columnSite,
    rowHamiltonianIndex, columnHamiltonianIndex,
    rowBasisRecord, columnBasisRecord, rowCell, columnCell
  },
  rowSiteIndex = Lookup[bond, "RowSite", $Failed];
  columnSiteIndex = Lookup[bond, "ColumnSite", $Failed];
  If[
    !IntegerQ[rowSiteIndex] ||
      !Between[rowSiteIndex, {1, Length[sites]}] ||
      !IntegerQ[columnSiteIndex] ||
      !Between[columnSiteIndex, {1, Length[sites]}],
    Return[$Failed]
  ];
  rowSite = sites[[rowSiteIndex]];
  columnSite = sites[[columnSiteIndex]];
  If[
    !Between[rowLocalIndex, {1, rowSite["Dimension"]}] ||
      !Between[columnLocalIndex, {1, columnSite["Dimension"]}],
    Return[$Failed]
  ];
  rowHamiltonianIndex = rowSite["BlockRange"][[rowLocalIndex]];
  columnHamiltonianIndex =
    columnSite["BlockRange"][[columnLocalIndex]];
  rowBasisRecord =
    basisOrderData["Orbitals"][[rowHamiltonianIndex]];
  columnBasisRecord =
    basisOrderData["Orbitals"][[columnHamiltonianIndex]];
  rowCell = Lookup[bond, "RowCell", $Failed];
  columnCell = Lookup[bond, "ColumnCell", $Failed];
  If[
    !ListQ[rowCell] || !ListQ[columnCell] ||
      Length[rowCell] =!= 3 || Length[columnCell] =!= 3,
    Return[$Failed]
  ];
  <|
    "Parameter" -> parameter,
    "Coefficient" -> coefficient,
    "Shell" -> shell,
    "ConstraintOrbitIndex" -> constraintOrbitIndex,
    "BondIndex" -> bondIndex,
    "RepresentativeBondIndex" -> representativeBondIndex,
    "Representative" -> (bondIndex === representativeBondIndex),
    "ColumnSiteIndex" -> columnSiteIndex,
    "ColumnCell" -> columnCell,
    "ColumnEndpoint" -> Lookup[bond, "Endpoints", {$Failed, $Failed}][[2]],
    "ColumnLocalOrbitalIndex" -> columnLocalIndex,
    "ColumnHamiltonianIndex" -> columnHamiltonianIndex,
    "ColumnBasisRecord" -> columnBasisRecord,
    "ColumnBasisState" -> displayedHamiltonianBasisState[columnBasisRecord],
    "RowSiteIndex" -> rowSiteIndex,
    "RowCell" -> rowCell,
    "RowEndpoint" -> Lookup[bond, "Endpoints", {$Failed, $Failed}][[1]],
    "RowLocalOrbitalIndex" -> rowLocalIndex,
    "RowHamiltonianIndex" -> rowHamiltonianIndex,
    "RowBasisRecord" -> rowBasisRecord,
    "RowBasisState" -> displayedHamiltonianBasisState[rowBasisRecord],
    "CellTranslation" -> Simplify[columnCell - rowCell]
  |>
];

(* Terms is the canonical real-space cache.  The representative matrix is the
   unique term whose BondIndex equals the stored representative for the same
   constraint orbit, so caching it a second time is unnecessary. *)
representativeHoppingsFromTerms[shellResult_Association] := Module[
  {representativeBondIndices, terms, representativeTerms},
  representativeBondIndices = Lookup[
    shellResult,
    "RepresentativeBondIndices",
    $Failed
  ];
  terms = Lookup[shellResult, "Terms", $Failed];
  If[!ListQ[representativeBondIndices] || !ListQ[terms], Return[$Failed]];
  representativeTerms = MapIndexed[
    Function[{bondIndex, orbitPosition},
      SelectFirst[
        terms,
        Lookup[#, "BondIndex", $Failed] === bondIndex &&
          Lookup[#, "OrbitIndex", $Failed] === First[orbitPosition] &,
        Missing["RepresentativeTermNotFound"]
      ]
    ],
    representativeBondIndices
  ];
  If[AnyTrue[representativeTerms, MissingQ], Return[$Failed]];
  representativeTerms = Lookup[representativeTerms, "Matrix", $Failed];
  If[And @@ (MatrixQ /@ representativeTerms), representativeTerms, $Failed]
];

resolveHoppingParameterShellData[
    session_Association,
    shell_Integer?Positive,
    hermitianQ : (True | False),
    method_String,
    validationLevel_String
  ] := Module[
  {compiledShells, realSpaceCache, key, shellResult, basisOrderData},
  compiledShells = Lookup[session, "CompiledBondShells", $Failed];
  realSpaceCache = Lookup[session, "RealSpaceShellCache", $Failed];
  If[
    !ListQ[compiledShells] || shell > Length[compiledShells],
    Return[Failure[
      "ShellOutOfRange",
      <|"Available" -> If[ListQ[compiledShells], Length[compiledShells], 0]|>
    ]]
  ];
  If[!AssociationQ[realSpaceCache], Return[$Failed]];
  key = realSpaceShellCacheKey[
    shell,
    hermitianQ,
    method,
    validationLevel
  ];
  shellResult = Lookup[realSpaceCache, key, Missing["NotSolved"]];
  If[MissingQ[shellResult], Return[Failure["UnsolvedShell", <||>]]];
  basisOrderData = buildHamiltonianBasisOrderData[session];
  If[
    !AssociationQ[shellResult] ||
      Lookup[shellResult, "Schema", None] =!=
        "MagneticTBRealSpaceShellResult" ||
      !AssociationQ[basisOrderData],
    Return[$Failed]
  ];
  <|
    "StaticShell" -> compiledShells[[shell]],
    "ShellResult" -> shellResult,
    "BasisOrderData" -> basisOrderData
  |>
];

buildHoppingParameterSourceData[
    shell_Integer?Positive,
    resolvedData_Association
  ] := Module[
  {
    staticShell, shellResult, basisOrderData, bonds, sites,
    parameters, parametersByOrbit, representativeHoppings,
    representativeBondIndices, records, failure
  },
  staticShell = resolvedData["StaticShell"];
  shellResult = resolvedData["ShellResult"];
  basisOrderData = resolvedData["BasisOrderData"];
  bonds = Lookup[staticShell, "Bonds", $Failed];
  sites = Lookup[staticShell, "Sites", $Failed];
  parameters = Lookup[shellResult, "Parameters", $Failed];
  parametersByOrbit = Lookup[shellResult, "ParametersByOrbit", $Failed];
  representativeHoppings = representativeHoppingsFromTerms[shellResult];
  representativeBondIndices = Lookup[
    shellResult,
    "RepresentativeBondIndices",
    $Failed
  ];
  If[
    !ListQ[bonds] || !ListQ[sites] || !ListQ[parameters] ||
      !ListQ[parametersByOrbit] ||
      !ListQ[representativeHoppings] ||
      !ListQ[representativeBondIndices] ||
      !SameQ[
        Length[parametersByOrbit],
        Length[representativeHoppings],
        Length[representativeBondIndices]
      ],
    Return[$Failed]
  ];
  records = Flatten@Table[
    Module[{bondIndex, bond, matrix, coefficient},
      bondIndex = representativeBondIndices[[orbit]];
      If[
        !IntegerQ[bondIndex] ||
          !Between[bondIndex, {1, Length[bonds]}],
        Return[$Failed, Module]
      ];
      bond = bonds[[bondIndex]];
      matrix = representativeHoppings[[orbit]];
      If[!MatrixQ[matrix], Return[$Failed, Module]];
      Flatten@Table[
        coefficient = exactHoppingParameterCoefficient[
          matrix[[row, column]],
          parameter
        ];
        If[
          MissingQ[coefficient],
          Nothing,
          makeHoppingParameterElementRecord[
            shell,
            orbit,
            bondIndex,
            bondIndex,
            parameter,
            coefficient,
            row,
            column,
            bond,
            sites,
            basisOrderData
          ]
        ],
        {parameter, parametersByOrbit[[orbit]]},
        {row, Length[matrix]},
        {column, Length[First[matrix]]}
      ]
    ],
    {orbit, Length[parametersByOrbit]}
  ];
  If[MemberQ[records, $Failed, Infinity], Return[$Failed]];
  failure = SelectFirst[records, !AssociationQ[#] &, Missing["NotFound"]];
  If[
    !MissingQ[failure] ||
      SortBy[DeleteDuplicates[Lookup[records, "Parameter"]], SymbolName] =!=
        SortBy[parameters, SymbolName],
    Return[$Failed]
  ];
  <|
    "Schema" -> "MagneticTBHoppingParameterSourceData",
    "SchemaVersion" -> 1,
    "Shell" -> shell,
    "Parameters" -> parameters,
    "Records" -> records
  |>
];

buildHoppingParameterOccurrenceData[
    shell_Integer?Positive,
    parameter_,
    resolvedData_Association
  ] := Module[
  {
    staticShell, shellResult, basisOrderData, bonds, sites,
    parameters, representativeBondIndices, terms, records
  },
  staticShell = resolvedData["StaticShell"];
  shellResult = resolvedData["ShellResult"];
  basisOrderData = resolvedData["BasisOrderData"];
  bonds = staticShell["Bonds"];
  sites = staticShell["Sites"];
  parameters = shellResult["Parameters"];
  representativeBondIndices = shellResult["RepresentativeBondIndices"];
  terms = shellResult["Terms"];
  If[!MemberQ[parameters, parameter],
    Return[Failure[
      "UnknownParameter",
      <|"Parameter" -> parameter, "Parameters" -> parameters|>
    ]]
  ];
  records = Flatten@Map[
    Function[term,
      Module[
        {bondIndex, orbitIndex, bond, matrix, coefficient},
        bondIndex = Lookup[term, "BondIndex", $Failed];
        orbitIndex = Lookup[term, "OrbitIndex", $Failed];
        matrix = Lookup[term, "Matrix", $Failed];
        If[
          !IntegerQ[bondIndex] ||
            !Between[bondIndex, {1, Length[bonds]}] ||
            !IntegerQ[orbitIndex] ||
            !Between[
              orbitIndex,
              {1, Length[representativeBondIndices]}
            ] ||
            !MatrixQ[matrix],
          Return[{$Failed}, Module]
        ];
        bond = bonds[[bondIndex]];
        Flatten@Table[
          coefficient = exactHoppingParameterCoefficient[
            matrix[[row, column]],
            parameter
          ];
          If[
            MissingQ[coefficient],
            Nothing,
            makeHoppingParameterElementRecord[
              shell,
              orbitIndex,
              representativeBondIndices[[orbitIndex]],
              bondIndex,
              parameter,
              coefficient,
              row,
              column,
              bond,
              sites,
              basisOrderData
            ]
          ],
          {row, Length[matrix]},
          {column, Length[First[matrix]]}
        ]
      ]
    ],
    terms
  ];
  If[
    MemberQ[records, $Failed, Infinity] ||
      !And @@ (AssociationQ /@ records),
    Return[$Failed]
  ];
  <|
    "Schema" -> "MagneticTBHoppingParameterOccurrenceData",
    "SchemaVersion" -> 1,
    "Shell" -> shell,
    "Parameter" -> parameter,
    "Records" -> records
  |>
];

renderHoppingEndpoint[record_Association, side_String] := Row[{
  "site ", record[side <> "SiteIndex"],
  " / H#", record[side <> "HamiltonianIndex"],
  " / ", record[side <> "BasisState"],
  " / cell ", record[side <> "Cell"],
  " / position ", record[side <> "Endpoint"]
}];

renderHoppingParameterRecords[
    data_Association,
    title_
  ] := Grid[
  Join[
    {{
      title,
      SpanFromLeft, SpanFromLeft, SpanFromLeft,
      SpanFromLeft, SpanFromLeft, SpanFromLeft,
      SpanFromLeft
    }},
    {{
      "Parameter",
      "Constraint orbit",
      "Bond",
      "Column / ket / source",
      "Row / bra / destination",
      "Cell translation",
      "Coefficient",
      "Representative"
    }},
    Map[
      {
        #["Parameter"],
        #["ConstraintOrbitIndex"],
        #["BondIndex"],
        renderHoppingEndpoint[#, "Column"],
        renderHoppingEndpoint[#, "Row"],
        #["CellTranslation"],
        #["Coefficient"],
        #["Representative"]
      } &,
      data["Records"]
    ]
  ],
  Frame -> All,
  Alignment -> Left
];

Options[showHoppingParameters] = {
  "Hermitian" -> True,
  "KernelMethod" -> "Iterative",
  "ValidationLevel" -> "Basic"
};

showHoppingParametersImpl[
    shell_Integer?Positive,
    parameterSelector_,
    hermitianQ_,
    method_,
    validationLevel_
  ] := Module[{resolvedData, sourceData, occurrenceData},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  If[
    !MemberQ[{True, False}, hermitianQ] ||
      !StringQ[method] || !StringQ[validationLevel],
    Message[
      showHoppingParameters::option,
      {hermitianQ, method, validationLevel}
    ];
    Return[$Failed]
  ];
  resolvedData = resolveHoppingParameterShellData[
    $CurrentModelSession,
    shell,
    hermitianQ,
    method,
    validationLevel
  ];
  Which[
    MatchQ[resolvedData, Failure["ShellOutOfRange", _Association]],
      Message[
        showHoppingParameters::shellrange,
        shell,
        resolvedData[[2, "Available"]]
      ];
      Return[$Failed],
    MatchQ[resolvedData, Failure["UnsolvedShell", _Association]],
      Message[
        showHoppingParameters::unsolved,
        shell,
        hermitianQ,
        method,
        validationLevel
      ];
      Return[$Failed],
    resolvedData === $Failed || FailureQ[resolvedData],
      Message[showHoppingParameters::data, shell];
      Return[$Failed]
  ];
  If[parameterSelector === All,
    sourceData = buildHoppingParameterSourceData[shell, resolvedData];
    If[sourceData === $Failed,
      Message[showHoppingParameters::data, shell];
      Return[$Failed]
    ];
    renderHoppingParameterRecords[
      sourceData,
      Row[{
        "Independent parameter sources for shell ", shell,
        " (representative hoppings)"
      }]
    ],
    occurrenceData = buildHoppingParameterOccurrenceData[
      shell,
      parameterSelector,
      resolvedData
    ];
    If[
      MatchQ[occurrenceData, Failure["UnknownParameter", _Association]],
      Message[
        showHoppingParameters::parameter,
        parameterSelector,
        shell,
        occurrenceData[[2, "Parameters"]]
      ];
      Return[$Failed]
    ];
    If[occurrenceData === $Failed || FailureQ[occurrenceData],
      Message[showHoppingParameters::data, shell];
      Return[$Failed]
    ];
    renderHoppingParameterRecords[
      occurrenceData,
      Row[{
        "All real-space occurrences of ", parameterSelector,
        " in shell ", shell
      }]
    ]
  ]
];

showHoppingParameters[
    shell_Integer?Positive,
    OptionsPattern[]
  ] := showHoppingParametersImpl[
  shell,
  All,
  OptionValue["Hermitian"],
  OptionValue["KernelMethod"],
  OptionValue["ValidationLevel"]
];

showHoppingParameters[
    shell_Integer?Positive,
    parameter_Symbol,
    OptionsPattern[]
  ] := showHoppingParametersImpl[
  shell,
  parameter,
  OptionValue["Hermitian"],
  OptionValue["KernelMethod"],
  OptionValue["ValidationLevel"]
];

showHoppingParameters[shell_, ___] := (
  Message[showHoppingParameters::shell, shell];
  $Failed
);

End[]
EndPackage[]
