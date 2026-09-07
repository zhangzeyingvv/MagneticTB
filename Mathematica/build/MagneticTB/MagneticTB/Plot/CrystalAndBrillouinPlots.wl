(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  validCrystalCellRangeQ,
  crystalCellEdges,
  buildCrystalStructureData,
  reciprocalTranslationVectors,
  firstBZFractionalRepresentatives,
  foldBandPathToFirstBZ,
  brillouinZoneVertices,
  buildBrillouinZoneData,
  renderBrillouinZone
];

showCrystalStructure::data =
  "The current model does not contain valid lattice, atomic-position, and magnetic-moment data.";
showCrystalStructure::range =
  "CellRange must be an integer n >= 0 or {{i1,i2},{j1,j2},{k1,k2}} with ordered integer bounds.";
showCrystalStructure::moment =
  "Magnetic moments remain nonnumeric after applying MomentRules. Supply numerical rules for symbolic moment components.";
showCrystalStructure::option =
  "AtomRadius and MomentScale must be Automatic or positive numeric values.";

showBrillouinZone::data =
  "The first Brillouin zone could not be constructed from the current reciprocal lattice.";
showBrillouinZone::path =
  "KPath must be Automatic, None, or a valid MagneticTB band path lying inside the first Brillouin zone.";
showBrillouinZone::option =
  "TranslationRange must be a positive integer and Tolerance must be a positive number.";

validCrystalCellRangeQ[range_] := MatchQ[
  range,
  {{_Integer, _Integer}, {_Integer, _Integer}, {_Integer, _Integer}}
] && And @@ (#[[1]] <= #[[2]] & /@ range);

crystalCellEdges[translation_List, lattice_?MatrixQ] := Module[
  {corners, pairs},
  corners = Association@Table[
    corner -> N[(corner + translation).lattice],
    {corner, Tuples[{0, 1}, 3]}
  ];
  pairs = Select[
    Subsets[Keys[corners], {2}],
    Total[Abs[#[[1]] - #[[2]]]] == 1 &
  ];
  Line[{corners[#[[1]]], corners[#[[2]]]}] & /@ pairs
];

Options[showCrystalStructure] = {
  "CellRange" -> {{0, 0}, {0, 0}, {0, 0}},
  "MomentRules" -> {},
  "MomentScale" -> Automatic,
  "AtomRadius" -> Automatic,
  "ShowAtomLabels" -> False,
  FontSize -> 14,
  ImageSize -> Large
};

buildCrystalStructureData[OptionsPattern[showCrystalStructure]] := Module[
  {
    session, model, metadata, lattice, atomPositions, cellRange,
    momentRules, momentScale, atomRadius, translations, baseRecords,
    records, invalidMoment, minimumLength
  },
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  model = Lookup[session, "ModelSpecification", $Failed];
  metadata = If[AssociationQ[model], Lookup[model, "Metadata", $Failed], $Failed];
  lattice = If[AssociationQ[metadata], Lookup[metadata, "LatticePlot", $Failed], $Failed];
  atomPositions = If[AssociationQ[metadata], Lookup[metadata, "AtomPositions", $Failed], $Failed];
  If[
    !MatrixQ[lattice, NumericQ] || Dimensions[lattice] =!= {3, 3} ||
      !ListQ[atomPositions] || atomPositions === {} ||
      !And @@ Map[
        ListQ[#] && # =!= {} &&
          And @@ (MatchQ[#, {{_, _, _}, {_, _, _}}] & /@ #) &,
        atomPositions
      ],
    Message[showCrystalStructure::data];
    Return[$Failed]
  ];
  cellRange = Replace[
    OptionValue["CellRange"],
    n_Integer?NonNegative :> {{-n, n}, {-n, n}, {-n, n}}
  ];
  If[!validCrystalCellRangeQ[cellRange],
    Message[showCrystalStructure::range];
    Return[$Failed]
  ];
  momentRules = OptionValue["MomentRules"];
  If[!ListQ[momentRules],
    Message[showCrystalStructure::moment];
    Return[$Failed]
  ];
  minimumLength = Min[Norm /@ N[lattice]];
  momentScale = Replace[
    OptionValue["MomentScale"],
    Automatic -> 0.55 minimumLength
  ];
  atomRadius = Replace[
    OptionValue["AtomRadius"],
    Automatic -> 0.08 minimumLength
  ];
  If[
    !NumericQ[momentScale] || momentScale <= 0 ||
      !NumericQ[atomRadius] || atomRadius <= 0,
    Message[showCrystalStructure::option];
    Return[$Failed]
  ];
  translations = Tuples[Range @@@ cellRange];
  baseRecords = Flatten@MapIndexed[
    Function[{orbit, orbitPosition},
      MapIndexed[
        <|
          "OrbitIndex" -> First[orbitPosition],
          "EquivalentIndex" -> First[#2],
          "FractionalPosition" -> #1[[1]],
          "FractionalMoment" -> #1[[2]]
        |> &,
        orbit
      ]
    ],
    atomPositions
  ];
  records = Flatten@Table[
    Map[
      Function[record,
        Module[{cartesianPosition, cartesianMoment, magneticQ},
          cartesianPosition = N[
            (record["FractionalPosition"] + translation).lattice
          ];
          cartesianMoment = N[
            (record["FractionalMoment"] /. momentRules).lattice
          ];
          magneticQ = VectorQ[cartesianMoment, NumericQ] &&
            Norm[cartesianMoment] > 10^-12;
          Join[
            record,
            <|
              "CellTranslation" -> translation,
              "CartesianPosition" -> cartesianPosition,
              "CartesianMoment" -> cartesianMoment,
              "Magnetic" -> magneticQ,
              "ArrowEnd" -> If[
                magneticQ,
                cartesianPosition +
                  momentScale Normalize[cartesianMoment],
                None
              ]
            |>
          ]
        ]
      ],
      baseRecords
    ],
    {translation, translations}
  ];
  invalidMoment = SelectFirst[
    records,
    !VectorQ[#1["CartesianMoment"], NumericQ] &,
    Missing["NotFound"]
  ];
  If[!MissingQ[invalidMoment],
    Message[showCrystalStructure::moment];
    Return[$Failed]
  ];
  <|
    "Lattice" -> lattice,
    "CellRange" -> cellRange,
    "Translations" -> translations,
    "AtomRadius" -> atomRadius,
    "MomentScale" -> momentScale,
    "AtomRecords" -> records,
    "AtomCount" -> Length[records],
    "MagneticAtomCount" -> Count[records, record_ /; TrueQ[record["Magnetic"]]],
    "CellEdges" -> Flatten[
      crystalCellEdges[#, lattice] & /@ translations
    ]
  |>
];

showCrystalStructure[opts : OptionsPattern[]] := Module[
  {data, records, radius, showLabels, fontSize},
  data = buildCrystalStructureData[opts];
  If[!AssociationQ[data], Return[$Failed]];
  records = data["AtomRecords"];
  radius = data["AtomRadius"];
  showLabels = TrueQ[OptionValue["ShowAtomLabels"]];
  fontSize = OptionValue[FontSize];
  Graphics3D[
    {
      Directive[GrayLevel[0.35], Thin],
      data["CellEdges"],
      Map[
        {
          ColorData[97][#1["OrbitIndex"]],
          Specularity[White, 20],
          Sphere[#1["CartesianPosition"], radius]
        } &,
        records
      ],
      Directive[Red, Thick, Arrowheads[0.035]],
      Map[
        Arrow[{#1["CartesianPosition"], #1["ArrowEnd"]}] &,
        Select[records, TrueQ[#1["Magnetic"]] &]
      ],
      If[
        showLabels,
        Map[
          Text[
            Style[
              Row[{#1["OrbitIndex"], ".", #1["EquivalentIndex"]}],
              Black,
              FontSize -> fontSize
            ],
            #1["CartesianPosition"] + {0, 0, 1.5 radius}
          ] &,
          records
        ],
        {}
      ]
    },
    Axes -> False,
    Boxed -> False,
    ViewProjection -> "Orthographic",
    Lighting -> "Neutral",
    PlotRange -> All,
    ImageSize -> OptionValue[ImageSize]
  ]
];

reciprocalTranslationVectors[reciprocal_?MatrixQ, range_Integer?Positive] :=
  N[#.reciprocal] & /@ DeleteCases[
    Tuples[Range[-range, range], 3],
    {0, 0, 0}
  ];

firstBZFractionalRepresentatives[
    point_List,
    reciprocal_?MatrixQ,
    range_Integer?Positive,
    tolerance_?NumericQ
  ] := Module[{translations, candidates, norms, minimumNorm},
  translations = Tuples[Range[-range, range], 3];
  candidates = N[point - #] & /@ translations;
  norms = Norm[#.reciprocal] & /@ candidates;
  minimumNorm = Min[norms];
  Sort@Pick[
    candidates,
    Map[# <= minimumNorm + tolerance Max[1., minimumNorm] &, norms]
  ]
];

foldBandPathToFirstBZ[
    path_List,
    reciprocal_?MatrixQ,
    range_Integer?Positive,
    tolerance_?NumericQ
  ] := Map[
  Function[segment,
    Module[{left, right, pair},
      left = firstBZFractionalRepresentatives[
        segment[[1, 1]], reciprocal, range, tolerance
      ];
      right = firstBZFractionalRepresentatives[
        segment[[1, 2]], reciprocal, range, tolerance
      ];
      If[left === {} || right === {}, Return[$Failed]];
      pair = First@MinimalBy[
        Tuples[{left, right}],
        Norm[(#[[1]] - #[[2]]).reciprocal] &
      ];
      {pair, segment[[2]]}
    ]
  ],
  path
];

brillouinZoneVertices[
    reciprocal_?MatrixQ,
    range_Integer?Positive,
    tolerance_?NumericQ
  ] := Module[
  {
    allVectors, planeVectors, allBounds, triples, rawVertices,
    point, determinant, scaleTolerance, vertices
  },
  allVectors = SortBy[
    reciprocalTranslationVectors[reciprocal, range],
    Norm
  ];
  If[Length[allVectors] < 4, Return[$Failed]];
  planeVectors = Take[allVectors, UpTo[24]];
  allBounds = (Norm[#]^2/2) & /@ allVectors;
  scaleTolerance = tolerance Max[1., Max[allBounds]];
  triples = Subsets[planeVectors, {3}];
  rawVertices = Reap[
    Do[
      determinant = Quiet@Det[triple];
      If[NumericQ[determinant] && Abs[determinant] > scaleTolerance,
        point = Quiet@Check[
          LinearSolve[triple, (Norm[#]^2/2) & /@ triple],
          $Failed
        ];
        If[
          VectorQ[point, NumericQ] &&
            Max[allVectors.point - allBounds] <= 10 scaleTolerance,
          Sow[Chop[point, scaleTolerance]]
        ]
      ],
      {triple, triples}
    ]
  ][[2]];
  If[rawVertices === {}, Return[$Failed]];
  vertices = DeleteDuplicates[
    First[rawVertices],
    Norm[#1 - #2] <= 10 scaleTolerance &
  ];
  If[Length[vertices] < 4, $Failed, vertices]
];

Options[showBrillouinZone] = {
  "KPath" -> Automatic,
  "ShowKPath" -> True,
  "TranslationRange" -> 2,
  "Tolerance" -> 10^-8,
  FontSize -> 14,
  ImageSize -> Large
};

buildBrillouinZoneData[OptionsPattern[showBrillouinZone]] := Module[
  {
    session, model, metadata, reciprocal, range, tolerance, vertices,
    mesh, suppliedPath, pathData, path, displayedPath, allVectors,
    allBounds, cartesianPathPoints, scaleTolerance
  },
  If[!ensureCurrentModelSession[], Return[$Failed]];
  session = $CurrentModelSession;
  model = Lookup[session, "ModelSpecification", $Failed];
  metadata = If[AssociationQ[model], Lookup[model, "Metadata", $Failed], $Failed];
  reciprocal = If[
    AssociationQ[metadata],
    N[
      Lookup[metadata, "ReciprocalLattice", $Failed] /.
        Lookup[metadata, "LatticeParameters", {}]
    ],
    $Failed
  ];
  range = OptionValue["TranslationRange"];
  tolerance = OptionValue["Tolerance"];
  If[
    !MatrixQ[reciprocal, NumericQ] || Dimensions[reciprocal] =!= {3, 3} ||
      !IntegerQ[range] || range < 1 ||
      !NumericQ[tolerance] || tolerance <= 0,
    Message[showBrillouinZone::option];
    Return[$Failed]
  ];
  vertices = brillouinZoneVertices[reciprocal, range, tolerance];
  If[vertices === $Failed,
    Message[showBrillouinZone::data];
    Return[$Failed]
  ];
  mesh = Quiet@Check[ConvexHullMesh[vertices], $Failed];
  If[Head[mesh] =!= BoundaryMeshRegion,
    Message[showBrillouinZone::data];
    Return[$Failed]
  ];
  suppliedPath = OptionValue["KPath"];
  pathData = Which[
    suppliedPath === None,
      <|"BravaisType" -> None, "Path" -> {}|>,
    suppliedPath === Automatic,
      standardKPathData["Tolerance" -> Max[tolerance, 10^-6]],
    validBandPathQ[suppliedPath],
      <|"BravaisType" -> "UserSupplied", "Path" -> suppliedPath|>,
    True,
      $Failed
  ];
  If[!AssociationQ[pathData],
    Message[showBrillouinZone::path];
    Return[$Failed]
  ];
  path = pathData["Path"];
  displayedPath = If[
    path === {},
    {},
    foldBandPathToFirstBZ[path, reciprocal, range, tolerance]
  ];
  If[displayedPath === $Failed,
    Message[showBrillouinZone::path];
    Return[$Failed]
  ];
  cartesianPathPoints = If[
    displayedPath === {},
    {},
    N[Flatten[displayedPath[[All, 1]], 1].reciprocal]
  ];
  allVectors = reciprocalTranslationVectors[reciprocal, range];
  allBounds = (Norm[#]^2/2) & /@ allVectors;
  scaleTolerance = tolerance Max[1., Max[allBounds]];
  If[
    cartesianPathPoints =!= {} &&
      AnyTrue[
        cartesianPathPoints,
        Max[allVectors.# - allBounds] > 20 scaleTolerance &
      ],
    Message[showBrillouinZone::path];
    Return[$Failed]
  ];
  <|
    "ReciprocalLattice" -> reciprocal,
    "Vertices" -> vertices,
    "Mesh" -> mesh,
    "FacetCount" -> Length[MeshPrimitives[mesh, 2]],
    "BravaisType" -> pathData["BravaisType"],
    "BZType" -> Lookup[pathData, "BZType", None],
    "KPath" -> path,
    "DisplayedKPath" -> displayedPath,
    "CartesianKPath" -> If[
      displayedPath === {},
      {},
      N[(#.reciprocal) & /@ displayedPath[[All, 1]]]
    ]
  |>
];

renderBrillouinZone[data_Association, showPath_, fontSize_, imageSize_] := Module[
  {mesh, path, reciprocal, labeledPoints},
  mesh = data["Mesh"];
  path = data["DisplayedKPath"];
  reciprocal = data["ReciprocalLattice"];
  labeledPoints = If[
    path === {},
    {},
    DeleteDuplicates[
      Flatten[
        Map[
          {
            {#[[2, 1]], N[#[[1, 1]].reciprocal]},
            {#[[2, 2]], N[#[[1, 2]].reciprocal]}
          } &,
          path
        ],
        1
      ]
    ]
  ];
  Graphics3D[
    {
      Directive[LightBlue, Opacity[0.18], EdgeForm[None]],
      MeshPrimitives[mesh, 2],
      Directive[GrayLevel[0.25], Thick, Opacity[1]],
      MeshPrimitives[mesh, 1],
      If[
        TrueQ[showPath] && path =!= {},
        {
          Directive[Red, Thick],
          Line /@ data["CartesianKPath"],
          Directive[Red, PointSize[0.018]],
          Point[labeledPoints[[All, 2]]],
          Map[
            Text[
              Style[#[[1]], Black, Bold, FontSize -> fontSize],
              #[[2]],
              {0, -1.2}
            ] &,
            labeledPoints
          ]
        },
        {}
      ]
    },
    Boxed -> False,
    Axes -> False,
    ViewProjection -> "Orthographic",
    PlotRange -> All,
    Lighting -> "Neutral",
    ImageSize -> imageSize
  ]
];

showBrillouinZone[opts : OptionsPattern[]] := Module[{data},
  data = buildBrillouinZoneData[opts];
  If[!AssociationQ[data], Return[$Failed]];
  renderBrillouinZone[
    data,
    OptionValue["ShowKPath"],
    OptionValue[FontSize],
    OptionValue[ImageSize]
  ]
];

End[]
EndPackage[]
