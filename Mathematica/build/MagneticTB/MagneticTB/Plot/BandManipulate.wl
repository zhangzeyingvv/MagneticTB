(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

bandParameterSymbols[hamiltonian_] := DeleteDuplicates@Cases[
  Unevaluated[hamiltonian],
  symbol_Symbol /;
    Context[Unevaluated[symbol]] =!= "System`" &&
      Context[Unevaluated[symbol]] =!= "MagneticTB`" &&
      !MemberQ[{kx, ky, kz}, Unevaluated[symbol]],
  Infinity
];

normalizeBandHamiltonian[hamiltonian_] := FixedPoint[
  Replace[
    #,
    {
      HoldPattern[MatrixForm[inner_, ___]] :> inner,
      HoldPattern[TableForm[inner_, ___]] :> inner,
      HoldPattern[TraditionalForm[inner_]] :> inner,
      HoldPattern[StandardForm[inner_]] :> inner
    }
  ] &,
  hamiltonian
];

unsupportedBandParameterExpressions[hamiltonian_] :=
  DeleteDuplicates@Cases[
    Unevaluated[hamiltonian],
    _Subscript | _Indexed,
    Infinity
  ];

unsupportedBandFunctionHeads[hamiltonian_] := DeleteDuplicates@Cases[
  Unevaluated[hamiltonian],
  HoldPattern[head_Symbol[___]] /;
    Context[Unevaluated[head]] =!= "System`",
  Infinity
];

bandManipulate::nonnumeric =
  "The Hamiltonian remains nonnumeric after assigning all band parameters and k. Remaining expression: `1`.";
bandManipulate::function =
  "The Hamiltonian contains unsupported unresolved function heads: `1`.";
bandManipulate::matrix =
  "The Hamiltonian must be a nonempty square matrix after display wrappers are removed; received head `1` with dimensions `2`.";
bandManipulate::parameter =
  "Subscript and Indexed expressions are not supported as band parameters. Use ordinary symbols instead: `1`.";

bandManipulate[pathstr_, npoint_, h_] := Module[
  {
    params,
    controls,
    controlSpecifications,
    klist,
    hamiltonian,
    matrixDimensions,
    unsupportedFunctions,
    unsupportedParameters
  },
  hamiltonian = normalizeBandHamiltonian[h];
  matrixDimensions = Quiet@Check[Dimensions[hamiltonian], Missing[]];
  If[
    !And[
      MatrixQ[hamiltonian],
      MatchQ[matrixDimensions, {_Integer, _Integer}],
      matrixDimensions[[1]] > 0,
      Equal @@ matrixDimensions
    ],
    Message[
      bandManipulate::matrix,
      Head[hamiltonian],
      matrixDimensions
    ];
    Return[Failure[
      "InvalidHamiltonian",
      <|
        "InputHead" -> Head[h],
        "NormalizedHead" -> Head[hamiltonian],
        "Dimensions" -> matrixDimensions
      |>
    ]]
  ];
  unsupportedFunctions = unsupportedBandFunctionHeads[hamiltonian];
  If[unsupportedFunctions =!= {},
    Message[bandManipulate::function, unsupportedFunctions];
    Return[Failure[
      "NonNumericHamiltonian",
      <|"UnsupportedFunctionHeads" -> unsupportedFunctions|>
    ]]
  ];
  unsupportedParameters =
    unsupportedBandParameterExpressions[hamiltonian];
  If[unsupportedParameters =!= {},
    Message[bandManipulate::parameter, unsupportedParameters];
    Return[Failure[
      "UnsupportedBandParameter",
      <|"UnsupportedParameterExpressions" -> unsupportedParameters|>
    ]]
  ];
  params = bandParameterSymbols[TrigToExp[hamiltonian]];
  controls = Table[Unique["bandParameter$"], {Length[params]}];
  controlSpecifications = MapThread[
    {{#1, 0, #2}, -1, 1} &,
    {controls, params}
  ];
  klist = 2. Pi Flatten[
    Subdivide[#[[1]], #[[2]], npoint] & /@ (Transpose[pathstr][[1]]),
    1
  ];
  Print["Number of params:", Length[params]];
  Print["params:", params];
  With[
    {
      parameterSymbols = params,
      dynamicSymbols = controls,
      specifications = controlSpecifications,
      momenta = klist,
      hamiltonian = hamiltonian
    },
    Manipulate[
      Module[{matrices, invalid},
        matrices = Table[
          N[
            hamiltonian /.
              Thread[parameterSymbols -> dynamicSymbols] /.
              {kx -> k[[1]], ky -> k[[2]], kz -> k[[3]]}
          ],
          {k, momenta}
        ];
        invalid = SelectFirst[matrices, !MatrixQ[#, NumericQ] &, Missing[]];
        If[
          MissingQ[invalid],
          ListPlot[
            Transpose[Eigenvalues /@ matrices],
            PlotRange -> All,
            PlotStyle -> Black
          ],
          Message[bandManipulate::nonnumeric, invalid];
          Style[
            "Hamiltonian is not numeric; inspect the remaining symbols.",
            Red
          ]
        ]
      ],
      Evaluate[Sequence @@ specifications],
      Button[
        "ExportData",
        Print[Thread[parameterSymbols -> dynamicSymbols]]
      ]
    ]
  ]
];

End[]

EndPackage[]
