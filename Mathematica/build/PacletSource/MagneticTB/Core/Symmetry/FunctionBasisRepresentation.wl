(* ::Package:: *)

BeginPackage["FunctionBasisRepresentation`"]

BasisRepresentationMatrix::usage =
  "BasisRepresentationMatrix[basis,images,variables] returns the matrix whose j-th column gives the expansion of images[[j]] in basis. Repeated scalar/vector polynomials are retained as independent ordered copies automatically.";
BasisRepresentationMatrices::usage =
  "BasisRepresentationMatrices[basis,imageSets,variables] computes several representation matrices in one coefficient-space factorization.";
PointOperationRepresentationMatrix::usage =
  "PointOperationRepresentationMatrix[basis,operation,variables] builds one point-operation representation. operation contains \"CoordinateImages\" and optional \"ValueMatrix\", \"Antiunitary\", and \"AntiunitaryMatrix\" entries.";
PointOperationRepresentationMatrices::usage =
  "PointOperationRepresentationMatrices[basis,operations,variables] builds a list of point-operation representation matrices using one shared polynomial coefficient solve and requires every generated matrix to be exactly unitary.";
TransformFunctionBasisForPointOperation::usage =
  "TransformFunctionBasisForPointOperation[basis,operation,variables] returns the actual scalar or spinor basis functions after applying the same coordinate, value-space, and antiunitary transformation convention used to build representation matrices.";

Begin["`Private`"]

ClearAll[
  normalizeFunctionBasis,
  automaticMultiplicityIndices,
  conjugatePolynomialInRealVariables,
  polynomialCoefficientData,
  coefficientColumn,
  firstNonzeroPosition,
  representationMatricesFromNormalizedBasis,
  transformBasisForOperation
];

normalizeFunctionBasis[basis_List] := Module[
  {vectorValuedQ, dimensions},
  If[basis === {}, Return[$Failed]];
  vectorValuedQ = ListQ[First[basis]];
  If[!And @@ (SameQ[ListQ[#], vectorValuedQ] & /@ basis),
    Return[$Failed]
  ];
  If[!vectorValuedQ, Return[List /@ basis]];
  dimensions = Length /@ basis;
  If[First[dimensions] < 1 || !SameQ @@ dimensions, Return[$Failed]];
  basis
];

automaticMultiplicityIndices[basisVectors_List] := Module[
  {canonicalVectors},
  canonicalVectors = Simplify[Expand[#]] & /@ basisVectors;
  Table[
    Count[
      Take[canonicalVectors, index],
      _?(SameQ[#, canonicalVectors[[index]]] &)
    ],
    {index, Length[canonicalVectors]}
  ]
];

conjugatePolynomialInRealVariables[expression_, variables_List] := Module[
  {rules},
  If[!PolynomialQ[Expand[expression], variables], Return[$Failed]];
  rules = CoefficientRules[Expand[expression], variables];
  Total@Map[
    Function[rule,
      Conjugate[Last[rule]] Times @@ MapThread[
        Power,
        {variables, First[rule]}
      ]
    ],
    rules
  ]
];

polynomialCoefficientData[functionVector_List, variables_List] := Module[
  {expanded},
  expanded = Expand /@ functionVector;
  If[!And @@ (PolynomialQ[#, variables] & /@ expanded), Return[$Failed]];
  Association[CoefficientRules[#, variables]] & /@ expanded
];

coefficientColumn[
    data_List,
    exponents_List,
    columnMultiplicity_,
    multiplicities_List
  ] := Flatten@Table[
  If[
    SameQ[multiplicity, columnMultiplicity],
    Lookup[data[[component]], Key[exponent], 0],
    0
  ],
  {multiplicity, multiplicities},
  {component, Length[data]},
  {exponent, exponents}
];

firstNonzeroPosition[row_List] := SelectFirst[
  Range[Length[row]],
  !TrueQ[PossibleZeroQ[row[[#]]]] &,
  Missing["NotFound"]
];

BasisRepresentationMatrices::basis =
  "The basis must be a nonempty list of scalar/vector polynomials.";
BasisRepresentationMatrices::images =
  "Every image set must have the same number and scalar/vector shape as the basis.";
BasisRepresentationMatrices::variables =
  "The variables must be a nonempty list of distinct symbols.";
BasisRepresentationMatrices::polynomial =
  "The basis and all images must be polynomials in the declared variables.";
BasisRepresentationMatrices::dependent =
  "The supplied function basis is linearly dependent.";
BasisRepresentationMatrices::closure =
  "At least one transformed function is outside the span of the supplied basis.";

representationMatricesFromNormalizedBasis[
    basisVectors_List,
    imageSetVectors_List,
    variables_List,
    basisMultiplicities_List
  ] := Module[
  {basisData, imageDataSets, allData, exponents, basisCoefficientMatrix,
   imageCoefficientMatrices, reducedTranspose, pivotRows, squareBasis,
   matrices},
  basisData = polynomialCoefficientData[#, variables] & /@ basisVectors;
  imageDataSets = Map[
    polynomialCoefficientData[#, variables] &,
    imageSetVectors,
    {2}
  ];
  If[
    MemberQ[basisData, $Failed] || MemberQ[imageDataSets, $Failed, Infinity],
    Message[BasisRepresentationMatrices::polynomial];
    Return[$Failed]
  ];

  allData = Join[basisData, Flatten[imageDataSets, 1]];
  exponents = DeleteDuplicates@Flatten[
    Map[Keys, allData, {2}],
    2
  ];
  exponents = Sort[exponents];

  basisCoefficientMatrix = Transpose[
    MapThread[
      coefficientColumn[
        #1,
        exponents,
        #2,
        DeleteDuplicates[basisMultiplicities]
      ] &,
      {basisData, basisMultiplicities}
    ]
  ];
  imageCoefficientMatrices = Map[
    Function[imageData,
      Transpose@MapThread[
        coefficientColumn[
          #1,
          exponents,
          #2,
          DeleteDuplicates[basisMultiplicities]
        ] &,
        {imageData, basisMultiplicities}
      ]
    ],
    imageDataSets
  ];

  reducedTranspose = RowReduce[Transpose[basisCoefficientMatrix]];
  pivotRows = DeleteMissing[firstNonzeroPosition /@ reducedTranspose];
  If[Length[pivotRows] =!= Length[basisVectors],
    Message[BasisRepresentationMatrices::dependent];
    Return[$Failed]
  ];
  squareBasis = basisCoefficientMatrix[[pivotRows]];
  matrices = Map[
    Simplify@LinearSolve[squareBasis, #[[pivotRows]]] &,
    imageCoefficientMatrices
  ];
  If[
    !And @@ MapThread[
      MatrixPredicates`ExactZeroMatrixQ[
        basisCoefficientMatrix.#1 - #2
      ] &,
      {matrices, imageCoefficientMatrices}
    ],
    Message[BasisRepresentationMatrices::closure];
    Return[$Failed]
  ];
  matrices
];

BasisRepresentationMatrices[
    basis_List,
    imageSets_List,
    variables_List
  ] := Module[
  {basisMultiplicities, basisVectors,
   imageSetVectors, basisCount, valueDimension},
  If[
    variables === {} || !DuplicateFreeQ[variables] ||
      !And @@ (MatchQ[#, _Symbol] & /@ variables),
    Message[BasisRepresentationMatrices::variables];
    Return[$Failed]
  ];
  basisVectors = normalizeFunctionBasis[basis];
  If[basisVectors === $Failed,
    Message[BasisRepresentationMatrices::basis];
    Return[$Failed]
  ];
  basisMultiplicities = automaticMultiplicityIndices[basisVectors];
  basisCount = Length[basisVectors];
  valueDimension = Length[First[basisVectors]];
  imageSetVectors = normalizeFunctionBasis /@ imageSets;
  If[
    MemberQ[imageSetVectors, $Failed] ||
      !And @@ (
        Length[#] == basisCount &&
          And @@ (Length[#] == valueDimension & /@ #) & /@
        imageSetVectors
      ),
    Message[BasisRepresentationMatrices::images];
    Return[$Failed]
  ];
  representationMatricesFromNormalizedBasis[
    basisVectors,
    imageSetVectors,
    variables,
    basisMultiplicities
  ]
];

BasisRepresentationMatrix[
    basis_List,
    images_List,
    variables_List
  ] := Module[{matrices},
  matrices = BasisRepresentationMatrices[basis, {images}, variables];
  If[matrices === $Failed, $Failed, First[matrices]]
];

PointOperationRepresentationMatrices::operation =
  "Each operation needs coordinate images matching the variable count and square value-space matrices matching the scalar/vector basis dimension.";
PointOperationRepresentationMatrices::unitary =
  "The generated point-operation representation is not unitary for operation indices `1`. Supply physically normalized basis functions; the compiler will not silently change basis.";

transformBasisForOperation[
    basisVectors_List,
    operation_Association,
    variables_List
  ] := Module[
  {valueDimension, coordinateImages, valueMatrix, antiunitary,
   antiunitaryMatrix, coordinateRules, conjugatedVector},
  valueDimension = Length[First[basisVectors]];
  coordinateImages = Lookup[operation, "CoordinateImages", $Failed];
  valueMatrix = Lookup[
    operation,
    "ValueMatrix",
    IdentityMatrix[valueDimension]
  ];
  antiunitary = Lookup[operation, "Antiunitary", False];
  antiunitaryMatrix = Lookup[
    operation,
    "AntiunitaryMatrix",
    IdentityMatrix[valueDimension]
  ];
  If[
    !BooleanQ[antiunitary] ||
      !ListQ[coordinateImages] || Length[coordinateImages] =!= Length[variables] ||
      !MatrixQ[valueMatrix] ||
      Dimensions[valueMatrix] =!= {valueDimension, valueDimension} ||
      !MatrixQ[antiunitaryMatrix] ||
      Dimensions[antiunitaryMatrix] =!= {valueDimension, valueDimension},
    Return[$Failed]
  ];
  coordinateRules = Thread[variables -> coordinateImages];
  Map[
    Function[functionVector,
      conjugatedVector = If[
        antiunitary,
        conjugatePolynomialInRealVariables[#, variables] & /@
          functionVector,
        functionVector
      ];
      If[
        MemberQ[conjugatedVector, $Failed],
        $Failed,
        Simplify[
          valueMatrix . (
            If[
              antiunitary,
              antiunitaryMatrix . conjugatedVector,
              functionVector
            ] /. coordinateRules
          )
        ]
      ]
    ],
    basisVectors
  ]
];

TransformFunctionBasisForPointOperation::basis =
  "The basis must be a nonempty list of scalar polynomials or equal-size polynomial vectors.";
TransformFunctionBasisForPointOperation::operation =
  "The point-operation record is incompatible with the supplied basis and variables.";
TransformFunctionBasisForPointOperation[
    basis_List,
    operation_Association,
    variables_List
  ] := Module[{scalarQ, basisVectors, transformed},
  scalarQ = basis =!= {} && !ListQ[First[basis]];
  basisVectors = normalizeFunctionBasis[basis];
  If[basisVectors === $Failed,
    Message[TransformFunctionBasisForPointOperation::basis];
    Return[$Failed]
  ];
  transformed = transformBasisForOperation[
    basisVectors,
    operation,
    variables
  ];
  If[transformed === $Failed || MemberQ[transformed, $Failed, Infinity],
    Message[TransformFunctionBasisForPointOperation::operation];
    Return[$Failed]
  ];
  If[scalarQ, First /@ transformed, transformed]
];
TransformFunctionBasisForPointOperation[___] := (
  Message[TransformFunctionBasisForPointOperation::basis];
  $Failed
);

PointOperationRepresentationMatrices[
    basis_List,
    operations_List,
    variables_List
  ] := Module[
  {basisMultiplicities, basisVectors,
   imageSetVectors, matrices, nonunitaryIndices},
  If[
    variables === {} || !DuplicateFreeQ[variables] ||
      !And @@ (MatchQ[#, _Symbol] & /@ variables),
    Message[BasisRepresentationMatrices::variables];
    Return[$Failed]
  ];
  basisVectors = normalizeFunctionBasis[basis];
  If[basisVectors === $Failed,
    Message[BasisRepresentationMatrices::basis];
    Return[$Failed]
  ];
  basisMultiplicities = automaticMultiplicityIndices[basisVectors];
  imageSetVectors = transformBasisForOperation[
      basisVectors,
      #,
      variables
    ] & /@ operations;
  If[MemberQ[imageSetVectors, $Failed, Infinity],
    Message[PointOperationRepresentationMatrices::operation];
    Return[$Failed]
  ];
  matrices = representationMatricesFromNormalizedBasis[
    basisVectors,
    imageSetVectors,
    variables,
    basisMultiplicities
  ];
  If[matrices === $Failed, Return[$Failed]];
  nonunitaryIndices = Pick[
    Range[Length[matrices]],
    MatrixPredicates`ExactUnitaryMatrixQ /@ matrices,
    False
  ];
  If[nonunitaryIndices =!= {},
    Message[
      PointOperationRepresentationMatrices::unitary,
      nonunitaryIndices
    ];
    Return[$Failed]
  ];
  matrices
];

PointOperationRepresentationMatrix[
    basis_List,
    operation_Association,
    variables_List
  ] := Module[{matrices},
  matrices = PointOperationRepresentationMatrices[
    basis,
    {operation},
    variables
  ];
  If[matrices === $Failed, $Failed, First[matrices]]
];

End[]

EndPackage[]
