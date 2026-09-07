(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  resolveBrokenSymmetryGeneratorIndices,
  sameBrokenSymmetryMagneticSiteQ,
  brokenSymmetryOrbitImageIndices,
  splitBrokenSymmetryOrbit,
  compileBrokenSymmetryInitRules
];

brokenSymmetryInitRules::selector =
  "The retained symmetry generators must be operation indices or unambiguous operation labels from the current init result; received `1`. Use {} for the identity-only subgroup.";
brokenSymmetryInitRules::group =
  "The retained generators `1` could not be closed to a subgroup of the symmetry group prepared by init.";
brokenSymmetryInitRules::orbit =
  "Wyckoff orbit `1` could not be split consistently under the retained subgroup; no init rules were produced.";
brokenSymmetryInitRules::representation =
  "The current model was created by initfromrep. Subgroup representation data must be regenerated explicitly; brokenSymmetryInitRules does not convert complete representation matrices into basisFunctions input.";

resolveBrokenSymmetryGeneratorIndices[
    selector_,
    labels_List,
    operationCount_Integer?Positive
  ] := Module[{items, resolved},
  If[selector === All, Return[Range[operationCount]]];
  items = If[ListQ[selector], selector, {selector}];
  resolved = Map[
    Function[item,
      If[IntegerQ[item],
        If[1 <= item <= operationCount, item, $Failed],
        With[{matches = Flatten@Position[labels, item, {1}]},
          If[Length[matches] == 1, First[matches], $Failed]
        ]
      ]
    ],
    items
  ];
  If[MemberQ[resolved, $Failed], $Failed, DeleteDuplicates[resolved]]
];

sameBrokenSymmetryMagneticSiteQ[first_, second_] :=
  MatchQ[first, {_List, _List}] &&
    MatchQ[second, {_List, _List}] &&
    Length[first[[1]]] == 3 && Length[second[[1]]] == 3 &&
    Length[first[[2]]] == 3 && Length[second[[2]]] == 3 &&
    MatrixPredicates`IntegerVectorQ[
      FullSimplify[first[[1]] - second[[1]]]
    ] &&
    And @@ (
      MatrixPredicates`ExactZeroExpressionQ /@
        FullSimplify[first[[2]] - second[[2]]]
    );

brokenSymmetryOrbitImageIndices[
    seed_,
    fullOrbit_List,
    subgroupSymmetryInformation_
  ] := Module[{images, matches},
  images = If[
    AssociationQ[subgroupSymmetryInformation],
    nativeSpinSpaceGroupAtomPositions[
      {seed},
      Values[subgroupSymmetryInformation]
    ],
    nativeAtomPositions[{seed}, subgroupSymmetryInformation]
  ];
  If[
    images === $Failed || !MatchQ[images, {_List}] ||
      MemberQ[images, $Failed, Infinity],
    Return[$Failed]
  ];
  matches = Map[
    Function[image,
      Select[
        Range[Length[fullOrbit]],
        sameBrokenSymmetryMagneticSiteQ[fullOrbit[[#]], image] &
      ]
    ],
    First[images]
  ];
  If[
    matches === {} || !And @@ (Length[#] == 1 & /@ matches),
    $Failed,
    Sort@DeleteDuplicates[First /@ matches]
  ]
];

splitBrokenSymmetryOrbit[
    fullOrbit_List,
    subgroupSymmetryInformation_
  ] := Module[{remaining, groups, seedIndex, indices},
  If[fullOrbit === {}, Return[$Failed]];
  remaining = Range[Length[fullOrbit]];
  groups = {};
  While[remaining =!= {},
    seedIndex = First[remaining];
    indices = brokenSymmetryOrbitImageIndices[
      fullOrbit[[seedIndex]],
      fullOrbit,
      subgroupSymmetryInformation
    ];
    If[
      indices === $Failed || indices === {} ||
        !MemberQ[indices, seedIndex] ||
        !SubsetQ[remaining, indices],
      Return[$Failed]
    ];
    AppendTo[groups, indices];
    remaining = Complement[remaining, indices]
  ];
  If[
    Sort[Flatten[groups]] === Range[Length[fullOrbit]],
    groups,
    $Failed
  ]
];

compileBrokenSymmetryInitRules[
    session_Association,
    selector_
  ] := Module[
  {
    initInput, modelSpecification, metadata,
    symmetryInformation, spatialActions, labels, generatorIndices,
    groupAlgebra, subgroupIndices,
    subgroupSymmetryInformation, spinSpaceGroupInformation,
    atomPositions, basisSpecification,
    representationMode, representationSource, splitGroups,
    failedOrbit, splitRecords
  },
  initInput = Lookup[session, "InitInput", $Failed];
  modelSpecification = Lookup[session, "ModelSpecification", $Failed];
  metadata = If[
    AssociationQ[modelSpecification],
    Lookup[modelSpecification, "Metadata", $Failed],
    $Failed
  ];
  If[
    !AssociationQ[initInput] || !AssociationQ[modelSpecification] ||
      !AssociationQ[metadata],
    Return[$Failed]
  ];
  symmetryInformation = Lookup[
    metadata,
    "SymmetryInformation",
    $Failed
  ];
  spinSpaceGroupInformation = Lookup[
    metadata,
    "SpinSpaceGroupInformation",
    Missing["NotSpinSpaceGroup"]
  ];
  spatialActions = Lookup[
    Lookup[modelSpecification, "Symmetry", <||>],
    "SpatialActions",
    $Failed
  ];
  atomPositions = Lookup[metadata, "AtomPositions", $Failed];
  basisSpecification = Lookup[
    metadata,
    "BasisSpecification",
    $Failed
  ];
  representationMode = Lookup[session, "RepresentationMode", $Failed];
  representationSource = Lookup[
    session,
    "RepresentationSource",
    "BasisFunctions"
  ];
  If[
    !ListQ[symmetryInformation] || symmetryInformation === {} ||
      !ListQ[spatialActions] ||
      Length[spatialActions] =!= Length[symmetryInformation] ||
      !ListQ[atomPositions] || !ListQ[basisSpecification] ||
      Length[atomPositions] =!= Length[basisSpecification],
    Return[$Failed]
  ];
  If[representationSource === "Matrices",
    Return[Failure["ExplicitRepresentationInput", <||>]]
  ];
  labels = symmetryInformation[[All, 1]];
  generatorIndices = resolveBrokenSymmetryGeneratorIndices[
    selector,
    labels,
    Length[spatialActions]
  ];
  If[generatorIndices === $Failed,
    Return[Failure["InvalidSelector", <|"Selector" -> selector|>]]
  ];
  groupAlgebra = Lookup[
    Lookup[modelSpecification, "Symmetry", <||>],
    "GroupAlgebra",
    $Failed
  ];
  If[!GroupAlgebra`GroupAlgebraQ[groupAlgebra],
    Return[Failure["InvalidParentGroup", <||>]]
  ];
  subgroupIndices = GroupAlgebra`GeneratedSubgroupIndices[
    groupAlgebra,
    generatorIndices
  ];
  If[subgroupIndices === $Failed,
    Return[Failure[
      "InvalidSubgroup",
      <|"GeneratorIndices" -> generatorIndices|>
    ]]
  ];
  subgroupSymmetryInformation = If[
    AssociationQ[spinSpaceGroupInformation],
    AssociationThread[
      labels[[subgroupIndices]],
      Values[spinSpaceGroupInformation][[subgroupIndices]]
    ],
    symmetryInformation[[subgroupIndices]]
  ];
  splitGroups = Map[
    splitBrokenSymmetryOrbit[#, subgroupSymmetryInformation] &,
    atomPositions
  ];
  failedOrbit = FirstPosition[splitGroups, $Failed, Missing["NotFound"]];
  If[!MissingQ[failedOrbit],
    Return[Failure[
      "OrbitSplitFailed",
      <|"Orbit" -> First[failedOrbit]|>
    ]]
  ];
  splitRecords = Flatten@Table[
    Table[
      <|
        "WyckoffSeed" ->
          atomPositions[[orbit, First[splitGroups[[orbit, part]]]]],
        "BasisFunctions" -> basisSpecification[[orbit]]
      |>,
      {part, Length[splitGroups[[orbit]]]}
    ],
    {orbit, Length[atomPositions]}
  ];
  {
    MagneticTB`lattice -> Lookup[initInput, "Lattice"],
    MagneticTB`lattpar -> Lookup[initInput, "LatticeParameters"],
    MagneticTB`wyckoffposition -> Lookup[splitRecords, "WyckoffSeed"],
    MagneticTB`symminformation -> subgroupSymmetryInformation,
    MagneticTB`basisFunctions -> Lookup[splitRecords, "BasisFunctions"],
    MagneticTB`debugQ -> Lookup[initInput, "Debug", False],
    MagneticTB`InitialBondShells ->
      Lookup[initInput, "InitialBondShells", 10],
    MagneticTB`GenerateSymmetryGroup -> False,
    MagneticTB`RepresentationMode -> representationMode
  }
];

brokenSymmetryInitRules[selector_] := Module[{result},
  If[!ensureCurrentModelSession[], Return[$Failed]];
  result = compileBrokenSymmetryInitRules[$CurrentModelSession, selector];
  Which[
    MatchQ[result, Failure["InvalidSelector", _Association]],
      Message[brokenSymmetryInitRules::selector, selector];
      $Failed,
    MatchQ[result, Failure["InvalidSubgroup" | "InvalidParentGroup", _Association]],
      Message[
        brokenSymmetryInitRules::group,
        Lookup[result[[2]], "GeneratorIndices", selector]
      ];
      $Failed,
    MatchQ[result, Failure["OrbitSplitFailed", _Association]],
      Message[brokenSymmetryInitRules::orbit, result[[2, "Orbit"]]];
      $Failed,
    MatchQ[result, Failure["ExplicitRepresentationInput", _Association]],
      Message[brokenSymmetryInitRules::representation];
      $Failed,
    result === $Failed || FailureQ[result],
      Message[brokenSymmetryInitRules::group, selector];
      $Failed,
    True,
      result
  ]
];

End[]
EndPackage[]
