(* ::Package:: *)

BeginPackage["OperatorSpaceBasis`"]

RectangularMatrixUnits::usage =
  "RectangularMatrixUnits[{m,n}] returns the row-major matrix-unit basis of m by n matrices.";
RectangularRealBasis::usage =
  "RectangularRealBasis[{m,n}] returns {E_ab} followed by {I E_ab}, a real basis of complex m by n matrices.";
HermitianMatrixBasis::usage =
  "HermitianMatrixBasis[n] returns an orthogonal real basis of n by n Hermitian matrices.";
MatrixToRealCoordinates::usage =
  "MatrixToRealCoordinates[T] returns row-major real parts followed by row-major imaginary parts.";
RealCoordinatesToMatrix::usage =
  "RealCoordinatesToMatrix[x,{m,n}] reconstructs an m by n complex matrix from real coordinates.";
RealHSInnerProduct::usage =
  "RealHSInnerProduct[a,b] is Re Tr[a^\[Dagger] b].";
BasisGramMatrix::usage =
  "BasisGramMatrix[basis] returns the Gram matrix for RealHSInnerProduct.";
CoordinatesInBasis::usage =
  "CoordinatesInBasis[matrix,basis] returns the real Hilbert--Schmidt coordinates of matrix in basis.";
ActionMatrixInBasis::usage =
  "ActionMatrixInBasis[basis,transform] represents a real-linear matrix transform in the supplied basis.";
DaggerCoordinateMatrix::usage =
  "DaggerCoordinateMatrix[{m,n}] maps row-major real coordinates of T to those of T^\[Dagger].";

Begin["`Private`"]

ClearAll[matrixUnit, rectangularDimensionsQ, nonemptyMatrixQ];

matrixUnit[m_Integer?Positive, n_Integer?Positive, a_Integer, b_Integer] :=
  Normal@SparseArray[{{a, b} -> 1}, {m, n}];

rectangularDimensionsQ[{m_Integer?Positive, n_Integer?Positive}] := True;
rectangularDimensionsQ[_] := False;
nonemptyMatrixQ[matrix_] := MatrixQ[matrix] &&
  Length[Dimensions[matrix]] === 2 &&
    And @@ (IntegerQ[#] && # > 0 & /@ Dimensions[matrix]);

RectangularMatrixUnits::dims =
  "Expected positive matrix dimensions {m,n}; received `1`.";
RectangularMatrixUnits[dims_?rectangularDimensionsQ] := Module[{m, n},
  {m, n} = dims;
  Flatten[
    Table[{matrixUnit[m, n, a, b]}, {a, 1, m}, {b, 1, n}],
    2
  ]
];
RectangularMatrixUnits[dims_] :=
  (Message[RectangularMatrixUnits::dims, dims]; $Failed);

RectangularRealBasis[dims_?rectangularDimensionsQ] := Module[{units},
  units = RectangularMatrixUnits[dims];
  Join[units, I units]
];
RectangularRealBasis[dims_] :=
  (Message[RectangularMatrixUnits::dims, dims]; $Failed);

HermitianMatrixBasis::dim =
  "Expected a positive matrix dimension; received `1`.";
HermitianMatrixBasis[n_Integer?Positive] := Module[
  {diagonal, symmetric, antisymmetric},
  diagonal = Table[matrixUnit[n, n, a, a], {a, 1, n}];
  symmetric = Flatten[
    Table[
      {matrixUnit[n, n, a, b] + matrixUnit[n, n, b, a]},
      {a, 1, n - 1}, {b, a + 1, n}
    ],
    2
  ];
  antisymmetric = Flatten[
    Table[
      {I (matrixUnit[n, n, a, b] - matrixUnit[n, n, b, a])},
      {a, 1, n - 1}, {b, a + 1, n}
    ],
    2
  ];
  Join[diagonal, symmetric, antisymmetric]
];
HermitianMatrixBasis[n_] :=
  (Message[HermitianMatrixBasis::dim, n]; $Failed);

MatrixToRealCoordinates[matrix_?nonemptyMatrixQ] := Join[
  Flatten[Simplify@ComplexExpand@Re[matrix]],
  Flatten[Simplify@ComplexExpand@Im[matrix]]
];
MatrixToRealCoordinates::matrix = "Expected a matrix.";
MatrixToRealCoordinates[_] :=
  (Message[MatrixToRealCoordinates::matrix]; $Failed);

RealCoordinatesToMatrix::length =
  "Expected `1` real coordinates for dimensions `2`; received `3`.";
RealCoordinatesToMatrix[
    coordinates_List,
    dims_?rectangularDimensionsQ
  ] := Module[{m, n, entryCount},
  {m, n} = dims;
  entryCount = m n;
  If[Length[coordinates] =!= 2 entryCount,
    Message[
      RealCoordinatesToMatrix::length,
      2 entryCount,
      dims,
      Length[coordinates]
    ];
    Return[$Failed]
  ];
  Partition[
    Take[coordinates, entryCount] + I Take[coordinates, -entryCount],
    n
  ]
];
RealCoordinatesToMatrix[coordinates_, dims_] :=
  (Message[RealCoordinatesToMatrix::length,
    If[rectangularDimensionsQ[dims], 2 Times @@ dims, "unknown"],
    dims, If[ListQ[coordinates], Length[coordinates], "not a list"]]; $Failed);

RealHSInnerProduct::shape =
  "Both arguments must be nonempty matrices with identical dimensions.";
RealHSInnerProduct[a_?nonemptyMatrixQ, b_?nonemptyMatrixQ] /;
    Dimensions[a] === Dimensions[b] :=
  Simplify@ComplexExpand@Re@Tr[ConjugateTranspose[a].b];
RealHSInnerProduct[_, _] :=
  (Message[RealHSInnerProduct::shape]; $Failed);

BasisGramMatrix::empty = "The operator basis must be nonempty.";
BasisGramMatrix[basis_List] /; Length[basis] > 0 := Module[{dimensions},
  If[!And @@ (nonemptyMatrixQ /@ basis),
    Message[BasisGramMatrix::empty]; Return[$Failed]
  ];
  dimensions = Dimensions /@ basis;
  If[!SameQ @@ dimensions,
    Message[RealHSInnerProduct::shape]; Return[$Failed]
  ];
  Table[RealHSInnerProduct[basis[[i]], basis[[j]]],
    {i, Length[basis]}, {j, Length[basis]}]
];
BasisGramMatrix[_] := (Message[BasisGramMatrix::empty]; $Failed);

CoordinatesInBasis::shape =
  "The matrix and all basis elements must have the same dimensions.";
CoordinatesInBasis[matrix_?nonemptyMatrixQ, basis_List] /;
    Length[basis] > 0 := Module[{dimensions, gram, rhs},
  If[!And @@ (nonemptyMatrixQ /@ basis),
    Message[CoordinatesInBasis::shape]; Return[$Failed]
  ];
  dimensions = Dimensions /@ Prepend[basis, matrix];
  If[!SameQ @@ dimensions,
    Message[CoordinatesInBasis::shape];
    Return[$Failed]
  ];
  gram = BasisGramMatrix[basis];
  If[MatrixRank[gram] =!= Length[basis],
    Message[CoordinatesInBasis::shape]; Return[$Failed]
  ];
  rhs = RealHSInnerProduct[#, matrix] & /@ basis;
  Simplify[LinearSolve[gram, rhs]]
];
CoordinatesInBasis[_, _] :=
  (Message[CoordinatesInBasis::shape]; $Failed);

ActionMatrixInBasis::image =
  "The transform did not return a matrix in the supplied operator space.";
ActionMatrixInBasis[basis_List, transform_] /;
    Length[basis] > 0 := Module[{columns},
  columns = CoordinatesInBasis[transform[#], basis] & /@ basis;
  If[MemberQ[columns, $Failed],
    Message[ActionMatrixInBasis::image];
    Return[$Failed]
  ];
  Simplify[Transpose[columns]]
];
ActionMatrixInBasis[_, _] :=
  (Message[ActionMatrixInBasis::image]; $Failed);

DaggerCoordinateMatrix::dims =
  "Expected positive matrix dimensions {m,n}; received `1`.";
DaggerCoordinateMatrix[dims_?rectangularDimensionsQ] := Module[
  {m, n, permutation},
  {m, n} = dims;
  permutation = SparseArray[
    Flatten[
      Table[
        {(b - 1) m + a, (a - 1) n + b} -> 1,
        {a, 1, m}, {b, 1, n}
      ],
      1
    ],
    {m n, m n}
  ];
  ArrayFlatten[
    {{permutation, 0 permutation}, {0 permutation, -permutation}}
  ]
];
DaggerCoordinateMatrix[dims_] :=
  (Message[DaggerCoordinateMatrix::dims, dims]; $Failed);

End[]

EndPackage[]
