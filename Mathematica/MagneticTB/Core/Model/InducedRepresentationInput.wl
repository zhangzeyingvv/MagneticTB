(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  nativeInducedLocalData,
  nativeInducedBasisOrderingData,
  nativeInducedFunctionBasisBlocks
];

nativeInducedFunctionBasisBlocks[
    transportedBases_List,
    operationRecords_List,
    imageTable_List,
    variables_List
  ] := Module[
  {operationCount, siteCount, imageSets, blocks, positions, matrices},
  operationCount = Length[operationRecords];
  siteCount = Length[transportedBases];
  If[Dimensions[imageTable] =!= {operationCount, siteCount}, Return[$Failed]];
  imageSets = Table[
    FunctionBasisRepresentation`TransformFunctionBasisForPointOperation[
      transportedBases[[source]], operationRecords[[operation]], variables
    ],
    {operation, operationCount}, {source, siteCount}
  ];
  If[MemberQ[imageSets, $Failed, Infinity], Return[$Failed]];
  blocks = ConstantArray[$Failed, {operationCount, siteCount}];
  Do[
    positions = Position[imageTable, target, {2}];
    matrices = FunctionBasisRepresentation`BasisRepresentationMatrices[
      transportedBases[[target]], Extract[imageSets, positions], variables
    ];
    If[matrices === $Failed, Return[$Failed]];
    Do[
      blocks[[Sequence @@ positions[[index]]]] = matrices[[index]],
      {index, Length[positions]}
    ],
    {target, siteCount}
  ];
  If[MemberQ[blocks, $Failed, Infinity], $Failed, blocks]
];

nativeInducedLocalData[
    siteOrbits_List,
    localBases_List,
    pointOperationRecords_List,
    sitePermutationData_Association,
    variables_List,
    group_Association,
    antiunitaryFlags_List
  ] := Table[
  Module[
    {indices, representatives, transportedBases, blocks, matrices, imageTable},
    imageTable = sitePermutationData["ImageSiteIndices"][[orbit]];
    indices = SitePermutationCompiler`SiteSymmetryOperationIndices[
      sitePermutationData, orbit, 1
    ];
    If[indices === $Failed, Return[$Failed, Module]];
    representatives =
      PhysicalRepresentation`AutomaticCosetRepresentativeIndices[
        group, imageTable, antiunitaryFlags, 1
      ];
    If[representatives === $Failed, Return[$Failed, Module]];
    transportedBases = Map[
      FunctionBasisRepresentation`TransformFunctionBasisForPointOperation[
        localBases[[orbit]],
        pointOperationRecords[[orbit, #]],
        variables
      ] &,
      representatives
    ];
    If[MemberQ[transportedBases, $Failed, Infinity], Return[$Failed, Module]];
    blocks = nativeInducedFunctionBasisBlocks[
      transportedBases,
      pointOperationRecords[[orbit]],
      imageTable,
      variables
    ];
    If[blocks === $Failed, Return[$Failed, Module]];
    matrices = blocks[[indices, 1]];
    <|
      "ReferenceSiteIndex" -> 1,
      "SiteSymmetryOperationIndices" -> indices,
      "SiteSymmetryMatrices" -> matrices,
      "CosetRepresentativeIndices" -> representatives,
      "TransportedBases" -> transportedBases,
      "LocalBlocks" -> blocks
    |>
  ],
  {orbit, Length[siteOrbits]}
];

nativeInducedBasisOrderingData[
    actionCompilation_Association,
    localBases_List,
    pointOperationRecords_List,
    variables_List
  ] := Module[{orbitData, records, operationCount},
  orbitData = Lookup[actionCompilation, "SiteSymmetryData", $Failed];
  operationCount = Length[
    Lookup[actionCompilation, "AntiunitaryFlags", {}]
  ];
  If[
    !ListQ[orbitData] || orbitData === {} || operationCount < 1 ||
      !And @@ (AssociationQ /@ orbitData),
    Return[$Failed]
  ];
  If[
    Length[localBases] =!= Length[orbitData] ||
      Length[pointOperationRecords] =!= Length[orbitData],
    Return[$Failed]
  ];
  records = MapThread[
    Function[{data, referenceBasis, operationRecords},
      Module[{representatives, transportedBases},
        representatives = Lookup[
          data,
          "CosetRepresentativeIndices",
          $Failed
        ];
        If[
          !ListQ[representatives] ||
            !And @@ (
              IntegerQ[#] && Between[#, {1, operationCount}] & /@
                representatives
            ),
          Return[$Failed, Module]
        ];
        transportedBases = Lookup[data, "TransportedBases", $Failed];
        If[transportedBases === $Failed,
          If[
            !ListQ[referenceBasis] || referenceBasis === {} ||
              !ListQ[operationRecords] ||
              Length[operationRecords] =!= operationCount,
            Return[$Failed, Module]
          ];
          transportedBases = Map[
            Function[operationIndex,
              FunctionBasisRepresentation`TransformFunctionBasisForPointOperation[
                referenceBasis,
                operationRecords[[operationIndex]],
                variables
              ]
            ],
            representatives
          ]
        ];
        If[!ListQ[transportedBases] || transportedBases === {} ||
            Length[transportedBases] =!= Length[representatives] ||
            MemberQ[transportedBases, $Failed, Infinity],
          Return[$Failed, Module]
        ];
        Join[
          KeyTake[
            data,
            {
              "ReferenceSiteIndex",
              "CosetRepresentativeIndices"
            }
          ],
          <|"TransportedBases" -> transportedBases|>
        ]
      ]
    ],
    {orbitData, localBases, pointOperationRecords}
  ];
  If[
    MemberQ[records, $Failed, Infinity],
    $Failed,
    <|"Mode" -> "Induced", "OrbitData" -> records|>
  ]
];

End[]

EndPackage[]
