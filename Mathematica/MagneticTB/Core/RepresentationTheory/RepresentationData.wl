(* ::Package:: *)

BeginPackage["RepresentationData`"]

AssembleBlockMonomialRepresentation::usage =
  "AssembleBlockMonomialRepresentation[images,blocks,dimensions,flags] assembles orbit and full matrices from operation-by-source image tables and local blocks.";
CreateRepresentationData::usage =
  "CreateRepresentationData[spec] validates exact unitary matrices aligned with a GroupAlgebra and antiunitary flags.";

Begin["`Private`"]

assembleOrbitMatrix[images_List, blocks_List, localDimension_Integer] :=
  ArrayFlatten@Table[
    If[images[[source]] === target,
      blocks[[source]],
      ConstantArray[0, {localDimension, localDimension}]
    ],
    {target, Length[images]}, {source, Length[images]}
  ];

blockDiagonalMatrix[matrices_List] := Module[{dimensions = Length /@ matrices},
  ArrayFlatten@Table[
    If[row === column, matrices[[row]],
      ConstantArray[0, {dimensions[[row]], dimensions[[column]]}]],
    {row, Length[matrices]}, {column, Length[matrices]}
  ]
];

AssembleBlockMonomialRepresentation::data =
  "Image tables, local blocks, dimensions, and flags are incompatible.";

AssembleBlockMonomialRepresentation[
    imageSiteIndices_List,
    localBlocks_List,
    localDimensions_List,
    antiunitaryFlags_List
  ] := Module[
  {orbitCount, operationCount, orbitMatrices, fullMatrices, blockDimensions},
  orbitCount = Length[localDimensions];
  operationCount = Length[antiunitaryFlags];
  If[orbitCount < 1 || operationCount < 1 ||
      Length[imageSiteIndices] =!= orbitCount ||
      Length[localBlocks] =!= orbitCount ||
      !And @@ (IntegerQ[#] && Positive[#] & /@ localDimensions) ||
      !And @@ (BooleanQ /@ antiunitaryFlags),
    Message[AssembleBlockMonomialRepresentation::data]; Return[$Failed]
  ];
  Do[
    If[Length[imageSiteIndices[[orbit]]] =!= operationCount ||
        Length[localBlocks[[orbit]]] =!= operationCount,
      Message[AssembleBlockMonomialRepresentation::data]; Return[$Failed]
    ];
    Do[
      If[!ListQ[imageSiteIndices[[orbit, operation]]] ||
          Sort[imageSiteIndices[[orbit, operation]]] =!=
            Range[Length[imageSiteIndices[[orbit, operation]]]] ||
          Length[localBlocks[[orbit, operation]]] =!=
            Length[imageSiteIndices[[orbit, operation]]] ||
          !And @@ (
            Dimensions[#] === {localDimensions[[orbit]], localDimensions[[orbit]]} &&
              MatrixPredicates`ExactUnitaryMatrixQ[#] & /@
              localBlocks[[orbit, operation]]
          ),
        Message[AssembleBlockMonomialRepresentation::data]; Return[$Failed]
      ],
      {operation, operationCount}
    ],
    {orbit, orbitCount}
  ];
  orbitMatrices = Table[
    assembleOrbitMatrix[
      imageSiteIndices[[orbit, operation]],
      localBlocks[[orbit, operation]],
      localDimensions[[orbit]]
    ],
    {orbit, orbitCount}, {operation, operationCount}
  ];
  fullMatrices = Table[
    blockDiagonalMatrix[orbitMatrices[[All, operation]]],
    {operation, operationCount}
  ];
  If[!And @@ (MatrixPredicates`ExactUnitaryMatrixQ /@ fullMatrices),
    Message[AssembleBlockMonomialRepresentation::data]; Return[$Failed]
  ];
  blockDimensions = Flatten@Table[
    ConstantArray[localDimensions[[orbit]],
      Length[imageSiteIndices[[orbit, 1]]]],
    {orbit, orbitCount}
  ];
  <|
    "OrbitRepresentationMatrices" -> orbitMatrices,
    "RepresentationMatrices" -> fullMatrices,
    "BlockDimensions" -> blockDimensions,
    "Dimension" -> Total[blockDimensions],
    "UnitaryVerified" -> True
  |>
];
AssembleBlockMonomialRepresentation[___] :=
  (Message[AssembleBlockMonomialRepresentation::data]; $Failed);

CreateRepresentationData::data =
  "Expected same-size exact unitary matrices, matching Boolean flags, and a verified GroupAlgebra.";
CreateRepresentationData[spec_Association] := Module[
  {matrices, flags, group, dimension},
  matrices = Lookup[spec, "RepresentationMatrices", $Failed];
  flags = Lookup[spec, "AntiunitaryFlags", $Failed];
  group = Lookup[spec, "GroupAlgebra", $Failed];
  If[!ListQ[matrices] || matrices === {} ||
      !GroupAlgebra`GroupAlgebraQ[group] ||
      Length[matrices] =!= group["Order"] ||
      !ListQ[flags] || Length[flags] =!= Length[matrices] ||
      !And @@ (BooleanQ /@ flags) ||
      !And @@ (MatrixPredicates`ExactSquareMatrixQ /@ matrices) ||
      !SameQ @@ (Dimensions /@ matrices) ||
      !And @@ (MatrixPredicates`ExactUnitaryMatrixQ /@ matrices),
    Message[CreateRepresentationData::data]; Return[$Failed]
  ];
  dimension = Length[First[matrices]];
  <|
    "RepresentationMatrices" -> matrices,
    "AntiunitaryFlags" -> flags,
    "GroupAlgebra" -> group,
    "Convention" -> Lookup[spec, "Convention", "columns are source states"],
    "Dimension" -> dimension,
    "UnitaryVerified" -> True,
    "Verified" -> True
  |>
];
CreateRepresentationData[_] := (Message[CreateRepresentationData::data]; $Failed);

End[]

EndPackage[]
