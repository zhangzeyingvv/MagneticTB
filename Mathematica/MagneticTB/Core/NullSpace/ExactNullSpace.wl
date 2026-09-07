(* Public single-matrix compatibility wrappers.  Common-kernel computation uses
   cqRawCommonKernel; these wrappers keep the established RREF/kernel schemas
   for fixtures and cross-language validation.  They intentionally stay in the
   full cyclotomic field so their canonical bases remain unchanged. *)

CyclotomicRREF[matrix_] := Module[{context},
  If[!cqMatrixQ[matrix],
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  context = matrix["context"];
  If[context["conductor"] === 1,
    cqRationalRREFResult[matrix],
    cqCoefficientRREFResult[matrix]
  ]
];

CyclotomicNullSpace[matrix_] := Module[{context},
  If[!cqMatrixQ[matrix],
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  context = matrix["context"];
  If[context["conductor"] === 1,
    cqRationalNullSpaceResult[matrix],
    cqCoefficientNullSpaceResult[matrix]
  ]
];
