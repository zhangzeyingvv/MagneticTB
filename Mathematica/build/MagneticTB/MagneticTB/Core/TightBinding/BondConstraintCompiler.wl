(* ::Package:: *)

BeginPackage["MagneticTBLinearAlgebra`"]

CompileBondConstraints::usage =
  "CompileBondConstraints[data] converts precompiled directed-bond data into a Hermitian or non-Hermitian constraint problem without recomputing group actions or bond orbits.";

Begin["`Private`"]

ClearAll[
  constraintRetargetAction,
  validBondActionRecordQ,
  validDirectedBondOrbitDataQ
];

constraintRetargetAction[action_Association, target_String] := Join[
  KeyDrop[action, "Target"],
  <|"Target" -> target|>
];

validBondActionRecordQ[action_, dims : {m_Integer?Positive, n_Integer?Positive}] :=
  AssociationQ[action] &&
    MatrixPredicates`ExactSquareMatrixQ[Lookup[action, "Left", $Failed]] &&
    MatrixPredicates`ExactSquareMatrixQ[Lookup[action, "Right", $Failed]] &&
    Dimensions[action["Left"]] === {m, m} &&
    Dimensions[action["Right"]] === {n, n} &&
    BooleanQ[Lookup[action, "Antiunitary", $Failed]] &&
    MemberQ[{"Same", "Reverse"}, Lookup[action, "Target", "Same"]];
validBondActionRecordQ[_, _] := False;

validDirectedBondOrbitDataQ[data_Association] := Module[
  {sites, bonds, reverseIndices, orbits, siteCount, bondCount,
   siteDimensions, continuousGenerators, memberBondIndices},
  If[!AssociationQ[data] || Lookup[data, "Schema", None] =!= "DirectedBondOrbitData" ||
      Lookup[data, "SchemaVersion", None] =!= 1 ||
      !IntegerQ[Lookup[data, "Shell", $Failed]] || data["Shell"] < 1,
    Return[False]
  ];
  sites = Lookup[data, "Sites", $Failed];
  bonds = Lookup[data, "Bonds", $Failed];
  reverseIndices = Lookup[data, "ReverseIndices", $Failed];
  orbits = Lookup[data, "DirectedOrbits", $Failed];
  If[!ListQ[sites] || sites === {} || !ListQ[bonds] || bonds === {} ||
      !ListQ[reverseIndices] || !ListQ[orbits] || orbits === {}, Return[False]];
  siteCount = Length[sites]; bondCount = Length[bonds];
  If[!And @@ MapIndexed[
      AssociationQ[#1] && Lookup[#1, "SiteIndex", $Failed] === First[#2] &&
        IntegerQ[Lookup[#1, "Dimension", $Failed]] && #1["Dimension"] > 0 &,
      sites], Return[False]];
  siteDimensions = Lookup[sites, "Dimension"];
  continuousGenerators = Lookup[data, "ContinuousSiteGenerators", None];
  If[continuousGenerators =!= None &&
      (!ListQ[continuousGenerators] || Length[continuousGenerators] =!= siteCount ||
        !And @@ MapThread[
          MatrixPredicates`ExactHermitianMatrixQ[#1] &&
            Dimensions[#1] === {#2, #2} &,
          {continuousGenerators, siteDimensions}]), Return[False]];
  If[!And @@ MapIndexed[
      Function[{bond, index},
        AssociationQ[bond] && Lookup[bond, "BondIndex", $Failed] === First[index] &&
          IntegerQ[Lookup[bond, "RowSite", $Failed]] &&
          IntegerQ[Lookup[bond, "ColumnSite", $Failed]] &&
          Between[bond["RowSite"], {1, siteCount}] &&
          Between[bond["ColumnSite"], {1, siteCount}] &&
          VectorQ[Lookup[bond, "Displacement", $Failed]] &&
          MatchQ[Lookup[bond, "Endpoints", $Failed], {_List, _List}]
      ], bonds], Return[False]];
  If[Length[reverseIndices] =!= bondCount ||
      !And @@ (IntegerQ[#] && Between[#, {1, bondCount}] & /@ reverseIndices) ||
      !And @@ Table[reverseIndices[[reverseIndices[[bond]]]] === bond,
        {bond, bondCount}] ||
      !And @@ Table[
        bonds[[reverseIndices[[bond]], "RowSite"]] === bonds[[bond, "ColumnSite"]] &&
          bonds[[reverseIndices[[bond]], "ColumnSite"]] === bonds[[bond, "RowSite"]] &&
          And @@ (MatrixPredicates`ExactZeroExpressionQ /@
            Simplify[bonds[[reverseIndices[[bond]], "Displacement"]] +
              bonds[[bond, "Displacement"]]]),
        {bond, bondCount}], Return[False]];
  If[!And @@ MapIndexed[
      Function[{orbit, index}, Module[{dims, representativeIndex, records},
        If[!AssociationQ[orbit] ||
            Lookup[orbit, "DirectedOrbitIndex", $Failed] =!= First[index],
          Return[False, Module]];
        dims = Lookup[orbit, "Dimensions", $Failed];
        representativeIndex = Lookup[orbit, "RepresentativeBondIndex", $Failed];
        records = Lookup[orbit, "DirectMemberRecords", $Failed];
        MatchQ[dims, {_Integer?Positive, _Integer?Positive}] &&
          IntegerQ[representativeIndex] && Between[representativeIndex, {1, bondCount}] &&
          Lookup[orbit, "RepresentativeBond", $Failed] === bonds[[representativeIndex]] &&
          dims === siteDimensions[[
            {bonds[[representativeIndex, "RowSite"]],
             bonds[[representativeIndex, "ColumnSite"]]}]] &&
          ListQ[records] && records =!= {} &&
          And @@ (AssociationQ[#] &&
              IntegerQ[Lookup[#, "BondIndex", $Failed]] &&
              Between[#["BondIndex"], {1, bondCount}] &&
              Lookup[#, "Mode", $Failed] === "Direct" &&
              validBondActionRecordQ[Lookup[#, "Operation", $Failed], dims] & /@
            records) &&
          And @@ (ListQ[Lookup[orbit, #, $Failed]] & /@ {
            "SameStabilizerOperations", "SameGeneratorOperations",
            "ReverseStabilizerOperations", "SameStabilizerOperationIndices",
            "StabilizerGeneratorIndices", "ReverseStabilizerOperationIndices",
            "ReverseCosetOperationIndices"}) &&
          And @@ (validBondActionRecordQ[#, dims] & /@
            Join[orbit["SameStabilizerOperations"],
              orbit["SameGeneratorOperations"],
              orbit["ReverseStabilizerOperations"]]) &&
          IntegerQ[Lookup[orbit, "ReverseBondIndex", $Failed]] &&
          orbit["ReverseBondIndex"] === reverseIndices[[representativeIndex]] &&
          IntegerQ[Lookup[orbit, "ReverseOrbitIndex", $Failed]] &&
          Between[orbit["ReverseOrbitIndex"], {1, Length[orbits]}] &&
          BooleanQ[Lookup[orbit, "SameGeneratorReductionVerified", $Failed]] &&
          BooleanQ[Lookup[orbit, "ReverseCosetVerified", $Failed]]
      ]], orbits], Return[False]];
  If[!And @@ Table[
      orbits[[orbits[[orbit, "ReverseOrbitIndex"]], "ReverseOrbitIndex"]] === orbit,
      {orbit, Length[orbits]}], Return[False]];
  memberBondIndices = Sort@Lookup[
    Flatten[Lookup[orbits, "DirectMemberRecords"], 1],
    "BondIndex",
    $Failed
  ];
  If[memberBondIndices =!= Range[bondCount], Return[False]];
  True
];
validDirectedBondOrbitDataQ[_] := False;

Options[CompileBondConstraints] = {"Hermitian" -> True};
CompileBondConstraints::data =
  "Expected data returned by CompileDirectedBondOrbits.";
CompileBondConstraints::hermitian =
  "The Hermitian option must be True or False; received `1`.";
CompileBondConstraints::reverse =
  "The directed-orbit reverse map is inconsistent at directed orbit `1`.";
CompileBondConstraints[data_Association, OptionsPattern[]] := Module[
  {hermitianQ, directedOrbits, reverseIndices, unseen, problemOrbits,
   orbitIndex, orbit, reverseOrbitIndex,
   hermiticityOperation, reverseGeneratorOperation,
   memberRecords, stabilizerOperations, constraintGeneratorOperations,
   continuousSiteGenerators, continuousGeneratorPair,
   constraintBlocks, representativeBond},
  If[!validDirectedBondOrbitDataQ[data],
    Message[CompileBondConstraints::data];
    Return[$Failed]
  ];
  hermitianQ = OptionValue["Hermitian"];
  If[!MemberQ[{True, False}, hermitianQ],
    Message[CompileBondConstraints::hermitian, hermitianQ];
    Return[$Failed]
  ];
  directedOrbits = data["DirectedOrbits"];
  reverseIndices = data["ReverseIndices"];
  continuousSiteGenerators = Lookup[
    data,
    "ContinuousSiteGenerators",
    None
  ];
  If[
    continuousSiteGenerators =!= None &&
      (!ListQ[continuousSiteGenerators] ||
        Length[continuousSiteGenerators] =!= Length[data["Sites"]]),
    Message[CompileBondConstraints::data];
    Return[$Failed]
  ];
  unseen = Range[Length[directedOrbits]];
  problemOrbits = {};
  While[unseen =!= {},
    orbitIndex = First[unseen];
    orbit = directedOrbits[[orbitIndex]];
    representativeBond = orbit["RepresentativeBond"];
    continuousGeneratorPair = If[
      continuousSiteGenerators === None,
      None,
      continuousSiteGenerators[[
        {
          representativeBond["RowSite"],
          representativeBond["ColumnSite"]
        }
      ]]
    ];
    reverseOrbitIndex = orbit["ReverseOrbitIndex"];
    If[
      !IntegerQ[reverseOrbitIndex] ||
        !Between[reverseOrbitIndex, {1, Length[directedOrbits]}] ||
        directedOrbits[[reverseOrbitIndex, "ReverseOrbitIndex"]] =!=
          orbitIndex,
      Message[CompileBondConstraints::reverse, orbitIndex];
      Return[$Failed]
    ];
    hermiticityOperation = <|
      "OperationIndex" -> "Hermiticity",
      "Left" -> IdentityMatrix[orbit["Dimensions"][[1]]],
      "Right" -> IdentityMatrix[orbit["Dimensions"][[2]]],
      "Antiunitary" -> False,
      "Target" -> "Reverse"
    |>;
    reverseGeneratorOperation = If[
      orbit["ReverseStabilizerOperations"] === {},
      Missing["NotApplicable"],
      constraintRetargetAction[
        First[orbit["ReverseStabilizerOperations"]],
        "Reverse"
      ]
    ];
    If[!hermitianQ,
      memberRecords = orbit["DirectMemberRecords"];
      stabilizerOperations = orbit["SameStabilizerOperations"];
      constraintGeneratorOperations = orbit["SameGeneratorOperations"];
      unseen = Rest[unseen],
      If[reverseOrbitIndex == orbitIndex,
        memberRecords = orbit["DirectMemberRecords"];
        stabilizerOperations = Join[
          If[
            orbit["ReverseBondIndex"] == orbit["RepresentativeBondIndex"],
            {hermiticityOperation},
            {}
          ],
          orbit["SameStabilizerOperations"],
          constraintRetargetAction[#, "Reverse"] & /@
            orbit["ReverseStabilizerOperations"]
        ];
        constraintGeneratorOperations = Join[
          If[
            orbit["ReverseBondIndex"] == orbit["RepresentativeBondIndex"],
            {hermiticityOperation},
            {}
          ],
          orbit["SameGeneratorOperations"],
          If[MissingQ[reverseGeneratorOperation], {}, {reverseGeneratorOperation}]
        ];
        unseen = DeleteCases[unseen, orbitIndex],
        memberRecords = Join[
          orbit["DirectMemberRecords"],
          Map[
            Function[member,
              <|
                "BondIndex" -> reverseIndices[[member["BondIndex"]]],
                "Mode" -> "Dagger",
                "Operation" -> member["Operation"]
              |>
            ],
            orbit["DirectMemberRecords"]
          ]
        ];
        stabilizerOperations = orbit["SameStabilizerOperations"];
        constraintGeneratorOperations = orbit["SameGeneratorOperations"];
        unseen = DeleteCases[unseen, orbitIndex | reverseOrbitIndex]
      ]
    ];
    constraintBlocks = MagneticTBLinearAlgebra`CompileRectangularConstraintBlocks[
      orbit["Dimensions"], constraintGeneratorOperations,
      continuousGeneratorPair];
    If[constraintBlocks === $Failed,
      Message[CompileBondConstraints::data]; Return[$Failed]
    ];
    AppendTo[
      problemOrbits,
      <|
        "ConstraintOrbitIndex" -> Length[problemOrbits] + 1,
        "SourceDirectedOrbitIndices" -> If[
          hermitianQ && reverseOrbitIndex =!= orbitIndex,
          {orbitIndex, reverseOrbitIndex},
          {orbitIndex}
        ],
        "RepresentativeBondIndex" -> orbit["RepresentativeBondIndex"],
        "RepresentativeBond" -> orbit["RepresentativeBond"],
        "Dimensions" -> orbit["Dimensions"],
        "CoordinateDimension" -> 2 Times @@ orbit["Dimensions"],
        "ConstraintBlocks" -> constraintBlocks,
        "ConstraintSources" -> Map[
          KeyTake[#, {"OperationIndex", "Antiunitary", "Target"}] &,
          constraintGeneratorOperations
        ],
        "HasContinuousConstraint" -> (continuousGeneratorPair =!= None),
        "MemberRecords" ->
          SortBy[memberRecords, Lookup[#, "BondIndex"] &],
        "StabilizerOperations" -> stabilizerOperations,
        "SameStabilizerOperationIndices" ->
          orbit["SameStabilizerOperationIndices"],
        "StabilizerGeneratorIndices" ->
          orbit["StabilizerGeneratorIndices"],
        "ReverseStabilizerOperationIndices" -> If[
          hermitianQ,
          orbit["ReverseStabilizerOperationIndices"],
          {}
        ],
        "ReverseRepresentativeOperationIndex" -> If[
          hermitianQ,
          orbit["ReverseRepresentativeOperationIndex"],
          Missing["NotApplicable"]
        ],
        "ReverseCosetOperationIndices" -> If[
          hermitianQ,
          orbit["ReverseCosetOperationIndices"],
          {}
        ],
        "SameGeneratorReductionVerified" ->
          orbit["SameGeneratorReductionVerified"],
        "ReverseCosetVerified" -> If[
          hermitianQ,
          orbit["ReverseCosetVerified"],
          True
        ],
        "GeneratorReductionVerified" -> TrueQ[
          orbit["SameGeneratorReductionVerified"] &&
            (!hermitianQ || orbit["ReverseCosetVerified"])
        ]
      |>
    ]
  ];
  <|
    "Schema" -> "BondConstraintProblem",
    "SchemaVersion" -> 1,
    "StaticShellKey" -> data["Shell"],
    "Shell" -> data["Shell"],
    "Hermitian" -> hermitianQ,
    "Orbits" -> problemOrbits
  |>
];
CompileBondConstraints[_, OptionsPattern[]] :=
  (Message[CompileBondConstraints::data]; $Failed);

End[]

EndPackage[]
