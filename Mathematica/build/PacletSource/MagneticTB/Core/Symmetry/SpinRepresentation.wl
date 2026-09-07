(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  nativeCoordinateImages,
  nativePointOperationRecords,
  nativeMSGSpinActions,
  nativeSpinMatrix
];

nativeCoordinateImages[rotation_?MatrixQ, lattice_?MatrixQ, variables_List] :=
  FullSimplify[
    Inverse[rotation] .
      (variables . Inverse[lattice]) . lattice
  ];

nativeCoordinateImages[operation_List, lattice_?MatrixQ, variables_List] :=
  nativeCoordinateImages[operation[[2]], lattice, variables];

nativePointOperationRecords[
    spatialActions_List,
    spinActions_List,
    antiunitaryFlags_List,
    lattice_?MatrixQ,
    variables_List,
    spinorQ_
  ] := MapThread[
  Function[{spatialAction, spinAction, antiunitary},
    Join[
      <|
        "CoordinateImages" -> nativeCoordinateImages[
          spatialAction["Rotation"],
          lattice,
          variables
        ],
        "Antiunitary" -> antiunitary
      |>,
      If[TrueQ[spinorQ],
        <|
          "ValueMatrix" -> ExpToTrig@nativeSpinMatrix[spinAction],
          "AntiunitaryMatrix" -> I PauliMatrix[2]
        |>,
        <||>
      ]
    ]
  ],
  {spatialActions, spinActions, antiunitaryFlags}
];

nativeMSGSpinActions[
    spatialActions_List,
    evaluatedLattice_?MatrixQ
  ] := Map[
  Function[spatialAction,
    Module[{cartesianRotation},
      cartesianRotation = FullSimplify[
        Transpose[evaluatedLattice] . spatialAction["Rotation"] .
          Inverse[Transpose[evaluatedLattice]]
      ];
      FullSimplify[Det[cartesianRotation] cartesianRotation]
    ]
  ],
  spatialActions
];

nativeSpinMatrix[rotation_?MatrixQ] := Module[
  {angle, axis, norm, orthogonalVector, normalVector, matrix, sign,
   weight, mixedWeight, xx, yy, zz, rotatedVector},
  matrix = If[TrueQ[FullSimplify[Det[rotation]] == -1], -rotation, rotation];
  axis = {
    matrix[[3, 2]] - matrix[[2, 3]],
    matrix[[1, 3]] - matrix[[3, 1]],
    matrix[[2, 1]] - matrix[[1, 2]]
  };
  norm = Simplify[Norm[axis]];
  If[TrueQ[norm == 0],
    angle = Pi Boole[Total[Diagonal[matrix]] < 3];
    rotatedVector = (matrix + IdentityMatrix[3])/2;
    axis = Normalize[
      Extract[rotatedVector, Ordering[Max /@ Abs[rotatedVector], -1]]
    ],
    {xx, yy, zz} = Simplify[axis/norm];
    sign = 2 UnitStep[zz] - 1;
    weight = -1/(sign + zz);
    mixedWeight = xx yy weight;
    orthogonalVector = {
      1 + sign weight xx xx,
      sign mixedWeight,
      -sign xx
    };
    normalVector = {
      mixedWeight,
      sign + weight yy yy,
      -yy
    };
    rotatedVector = matrix.orthogonalVector;
    angle = Arg@Simplify[
      rotatedVector.orthogonalVector + I rotatedVector.normalVector
    ]
  ];
  ExpToTrig@MatrixExp[
    -I FullSimplify[angle]
      Sum[
        PauliMatrix[index] FullSimplify[Normalize[axis][[index]]],
        {index, 3}
      ]/2
  ]
];

End[]

EndPackage[]
