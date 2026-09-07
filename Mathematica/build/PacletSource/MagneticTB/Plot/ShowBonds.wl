(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  validShowBondSiteQ,
  validShowBondRecordQ,
  buildBondShellDisplayData,
  renderBondShellGrid
];

showbonds::shell =
  "The shell index must be a positive integer; received `1`.";
showbonds::shellrange =
  "Bond shell `1` was not prepared by init; only `2` shells are available. Rerun init with InitialBondShells -> `1`.";
showbonds::data =
  "The cached bond data for shell `1` is missing or malformed.";

validShowBondSiteQ[site_] :=
  AssociationQ[site] &&
    And @@ (KeyExistsQ[site, #] & /@ {
      "SiteIndex", "OrbitIndex", "EquivalentIndex", "Position"
    }) &&
    IntegerQ[site["SiteIndex"]] && site["SiteIndex"] > 0 &&
    ListQ[site["Position"]];

validShowBondRecordQ[bond_, siteCount_Integer?Positive] :=
  AssociationQ[bond] &&
    And @@ (KeyExistsQ[bond, #] & /@ {
      "RowSite", "ColumnSite", "Endpoints"
    }) &&
    IntegerQ[bond["RowSite"]] &&
    Between[bond["RowSite"], {1, siteCount}] &&
    IntegerQ[bond["ColumnSite"]] &&
    Between[bond["ColumnSite"], {1, siteCount}] &&
    MatchQ[bond["Endpoints"], {_List, _List}];

(* Pure data layer: consume only the session prepared by init.  Keeping this
   separate from Grid construction makes the bond report testable and keeps
   Plot free of bond-search or bond-orbit compilation work. *)
buildBondShellDisplayData[
    session_Association,
    shell_Integer?Positive
  ] := Module[
  {
    compiledShells, modelSpecification, bondClasses, staticShell, shellClass,
    sites, bonds, siteCount, bondLength, rows
  },
  compiledShells = Lookup[session, "CompiledBondShells", $Failed];
  modelSpecification = Lookup[session, "ModelSpecification", $Failed];
  bondClasses = If[
    AssociationQ[modelSpecification],
    Lookup[modelSpecification, "BondClasses", $Failed],
    $Failed
  ];
  If[
    !ListQ[compiledShells] || !ListQ[bondClasses] ||
      shell > Length[compiledShells] || shell > Length[bondClasses],
    Return[$Failed]
  ];

  staticShell = compiledShells[[shell]];
  shellClass = bondClasses[[shell]];
  If[
    !AssociationQ[staticShell] ||
      Lookup[staticShell, "Shell", Missing["Shell"]] =!= shell ||
      !ListQ[shellClass] || shellClass === {},
    Return[$Failed]
  ];

  sites = Lookup[staticShell, "Sites", $Failed];
  bonds = Lookup[staticShell, "Bonds", $Failed];
  If[
    !ListQ[sites] || sites === {} || !ListQ[bonds] ||
      !And @@ (validShowBondSiteQ /@ sites),
    Return[$Failed]
  ];
  siteCount = Length[sites];
  If[!And @@ (validShowBondRecordQ[#, siteCount] & /@ bonds),
    Return[$Failed]
  ];
  If[
    !MatchQ[shellClass[[1]], {_?NumericQ, _Integer, _List}],
    Return[$Failed]
  ];
  bondLength = shellClass[[1, 1]];

  rows = Map[
    Function[site,
      Module[{siteBonds},
        siteBonds = Select[
          bonds,
          Lookup[#, "RowSite", Missing["RowSite"]] ===
            site["SiteIndex"] &
        ];
        <|
          "SiteIndex" -> site["SiteIndex"],
          "OrbitIndex" -> site["OrbitIndex"],
          "EquivalentIndex" -> site["EquivalentIndex"],
          "AtomPosition" -> site["Position"],
          "BondCount" -> Length[siteBonds],
          "TargetPositions" -> Lookup[siteBonds, "Endpoints"][[All, 2]]
        |>
      ]
    ],
    sites
  ];

  <|
    "Schema" -> "MagneticTBBondShellDisplayData",
    "SchemaVersion" -> 1,
    "Shell" -> shell,
    "NeighbourOrder" -> shell - 1,
    "BondLength" -> bondLength,
    "SiteCount" -> siteCount,
    "DirectedBondCount" -> Length[bonds],
    "DirectedOrbitCount" -> Length[
      Lookup[staticShell, "DirectedOrbits", {}]
    ],
    "Rows" -> rows
  |>
];

renderBondShellGrid[data_Association] := Module[
  {title, body},
  title = Row[{
    "Bond shell ", data["Shell"],
    " (", data["NeighbourOrder"], "-th neighbour), ",
    "bond length = ", data["BondLength"], "; ",
    data["DirectedBondCount"], " directed bonds, ",
    data["DirectedOrbitCount"], " directed orbits"
  }];
  body = Map[
    {
      #["SiteIndex"],
      #["AtomPosition"],
      #["BondCount"],
      If[#["TargetPositions"] === {}, None, #["TargetPositions"]]
    } &,
    data["Rows"]
  ];
  Grid[
    Join[
      {{title, SpanFromLeft, SpanFromLeft, SpanFromLeft}},
      {{
        "Site index",
        "Atom position (fractional)",
        "Num. of directed bonds",
        Row[{DisplayForm@SubsuperscriptBox["p", "k", "'"],
          " (fractional)"}]
      }},
      body
    ],
    Frame -> All,
    Alignment -> Left
  ]
];

showbonds[shell_Integer?Positive] := Module[
  {session, availableShells, data},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  availableShells = Length[Lookup[session, "CompiledBondShells", {}]];
  If[shell > availableShells,
    Message[showbonds::shellrange, shell, availableShells, shell];
    Return[$Failed]
  ];
  data = buildBondShellDisplayData[session, shell];
  If[data === $Failed,
    Message[showbonds::data, shell];
    Return[$Failed]
  ];
  renderBondShellGrid[data]
];

showbonds[shell_] := (
  Message[showbonds::shell, shell];
  $Failed
);

End[]
EndPackage[]
