cqVerifyCommonKernel[constraints_List, basisMatrix_] := Module[
  {rawBasis, rawConstraint, residual, constraintIndex},
  If[!cqMatrixQ[basisMatrix], Return[False]];
  rawBasis = cqCoefficientRawFromExact[basisMatrix, False];
  If[FailureQ[rawBasis], Return[False]];
  For[constraintIndex = 1, constraintIndex <= Length[constraints], constraintIndex++,
    If[!cqMatrixQ[constraints[[constraintIndex]]] ||
       constraints[[constraintIndex, "context", "context_id"]] =!=
         basisMatrix["context"]["context_id"],
      Return[False]
    ];
    rawConstraint = cqCoefficientRawFromExact[
      constraints[[constraintIndex]], False
    ];
    If[FailureQ[rawConstraint], Return[False]];
    residual = cqCoefficientRawMultiply[rawConstraint, rawBasis];
    If[FailureQ[residual] || !cqCoefficientRawZeroQ[residual], Return[False]]
  ];
  True
];

cqValidateCompiledProblem[compiled_] := Module[
  {context, coordinateDimension, constraintCount, constraints, matrix, index},
  If[!AssociationQ[compiled] ||
     Lookup[compiled, "type", None] =!= "compiled_cyclotomic_problem",
    Return[cqFailure["MalformedSerialization", <|"input" -> HoldForm[compiled]|>]]
  ];
  context = Lookup[compiled, "context", None];
  coordinateDimension = Lookup[compiled, "coordinate_dimension", None];
  constraintCount = Lookup[compiled, "constraint_count", None];
  constraints = Lookup[compiled, "constraints", None];
  If[!cqContextQ[context] ||
     Lookup[compiled, "conductor", None] =!= context["conductor"] ||
     !IntegerQ[coordinateDimension] || coordinateDimension < 0 ||
     !IntegerQ[constraintCount] || !ListQ[constraints] ||
     constraintCount =!= Length[constraints],
    Return[cqFailure["MalformedSerialization", <|"input" -> compiled|>]]
  ];
  For[index = 1, index <= Length[constraints], index++,
    matrix = constraints[[index]];
    If[!cqMatrixQ[matrix],
      Return[cqFailure["MalformedSerialization", <|"constraint_index" -> index|>]]
    ];
    If[matrix["context"]["context_id"] =!= context["context_id"],
      Return[cqFailure["ConductorMismatch", <|"constraint_index" -> index|>]]
    ];
    If[matrix["columns"] =!= coordinateDimension,
      Return[cqFailure["DimensionMismatch", <|
        "constraint_index" -> index,
        "coordinate_dimension" -> coordinateDimension,
        "columns" -> matrix["columns"]
      |>]]
    ]
  ];
  True
];

cqWorkPlanDomain[plan_Association] := Switch[
  Lookup[plan, "backend_kind", None],
  "rational", <|"kind" -> "rational"|>,
  "coefficient", <|
    "kind" -> "coefficient",
    "context" -> Lookup[plan, "field_context", None]
  |>,
  _, cqFailure["MalformedSerialization", <|"work_plan" -> plan|>]
];

cqValidateRawDomain[
  domain_, coordinateDimension_, rawConstraints_
] := Module[{kind, fieldContext, matrix, index},
  If[!AssociationQ[domain] || !IntegerQ[coordinateDimension] ||
     coordinateDimension < 0 || !ListQ[rawConstraints],
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  kind = Lookup[domain, "kind", None];
  If[kind =!= "rational" && kind =!= "coefficient",
    Return[cqFailure["MalformedSerialization", <|"domain" -> domain|>]]
  ];
  If[kind === "coefficient",
    fieldContext = Lookup[domain, "context", None];
    If[!cqCoefficientFieldContextQ[fieldContext],
      Return[cqFailure["MalformedSerialization", <|"domain" -> domain|>]]
    ]
  ];
  For[index = 1, index <= Length[rawConstraints], index++,
    matrix = rawConstraints[[index]];
    If[
      If[kind === "rational",
        !cqRationalRawMatrixStructureQ[matrix],
        !cqCoefficientRawMatrixStructureQ[matrix]
      ],
      Return[cqFailure["MalformedSerialization", <|"constraint_index" -> index|>]]
    ];
    If[matrix["columns"] =!= coordinateDimension,
      Return[cqFailure["DimensionMismatch", <|"constraint_index" -> index|>]]
    ];
    If[kind === "coefficient" &&
       cqCoefficientFieldID[matrix["context"]] =!=
         cqCoefficientFieldID[fieldContext],
      Return[cqFailure["ConductorMismatch", <|"constraint_index" -> index|>]]
    ]
  ];
  True
];

cqCommonKernelWorkPlanQ[plan_] := Module[
  {backendKind, publicContext, publicField, fieldContext, coordinateDimension,
   constraintCount, rawConstraints, liftKind, liftMatrix, method},
  If[!AssociationQ[plan] ||
     Lookup[plan, "type", None] =!= "common_kernel_work_plan",
    Return[False]
  ];
  backendKind = Lookup[plan, "backend_kind", None];
  publicContext = Lookup[plan, "public_context", None];
  fieldContext = Lookup[plan, "field_context", None];
  coordinateDimension = Lookup[plan, "coordinate_dimension", None];
  constraintCount = Lookup[plan, "constraint_count", None];
  rawConstraints = Lookup[plan, "raw_constraints", None];
  liftKind = Lookup[plan, "lift_kind", None];
  liftMatrix = Lookup[plan, "lift_matrix", None];
  method = Lookup[plan, "method", None];
  If[!cqContextQ[publicContext] ||
     !IntegerQ[coordinateDimension] || coordinateDimension < 0 ||
     !IntegerQ[constraintCount] || !ListQ[rawConstraints] ||
     constraintCount =!= Length[rawConstraints] ||
     !MemberQ[{"identity", "linear"}, liftKind] ||
     !ListQ[liftMatrix] || !StringQ[method],
    Return[False]
  ];
  Switch[backendKind,
    "rational",
      publicContext["conductor"] === 1 && fieldContext === None &&
        liftKind === "identity" && liftMatrix === {{1}},
    "coefficient",
      publicField = cqCoefficientFieldFromCyclotomic[publicContext];
      !FailureQ[publicField] && cqCoefficientFieldContextQ[fieldContext] &&
        Length[liftMatrix] === publicContext["degree"] &&
        And @@ (ListQ[#] && Length[#] === fieldContext["degree"] & /@ liftMatrix) &&
        And @@ (cqExactRationalQ /@ Flatten[liftMatrix]) &&
        (liftKind =!= "identity" ||
          (cqCoefficientFieldID[fieldContext] ===
             cqCoefficientFieldID[publicField] &&
           liftMatrix === IdentityMatrix[publicContext["degree"]])),
    _, False
  ]
];

cqRationalWorkPlan[compiled_Association, rawConstraints_List] := <|
  "type" -> "common_kernel_work_plan",
  "backend_kind" -> "rational",
  "public_context" -> compiled["context"],
  "field_context" -> None,
  "coordinate_dimension" -> compiled["coordinate_dimension"],
  "constraint_count" -> Length[rawConstraints],
  "raw_constraints" -> rawConstraints,
  "lift_kind" -> "identity",
  "lift_matrix" -> {{1}},
  "method" -> "iterative"
|>;

cqRawDomainIdentity[domain_Association, dimension_Integer] := Switch[
  domain["kind"],
  "rational", cqRationalRawIdentity[dimension],
  "coefficient", cqCoefficientRawIdentity[domain["context"], dimension],
  _, cqFailure["MalformedSerialization", <|"domain" -> domain|>]
];

cqRawDomainMatrix[
  domain_Association, rows_Integer, columns_Integer, data_List
] := Switch[
  domain["kind"],
  "rational", cqRationalRawMatrix[rows, columns, data],
  "coefficient", cqCoefficientRawMatrix[domain["context"], rows, columns, data],
  _, cqFailure["MalformedSerialization", <|"domain" -> domain|>]
];

cqRawDomainMultiply[domain_Association, left_Association, right_Association] :=
  Switch[
    domain["kind"],
    "rational", cqRationalRawMultiply[left, right],
    "coefficient", cqCoefficientRawMultiply[left, right],
    _, cqFailure["MalformedSerialization", <|"domain" -> domain|>]
  ];

cqRawDomainZeroQ[domain_Association, matrix_Association] := Switch[
  domain["kind"],
  "rational", cqRationalRawZeroQ[matrix],
  "coefficient", cqCoefficientRawZeroQ[matrix],
  _, False
];

cqRawDomainNullSpace[domain_Association, matrix_Association] := Switch[
  domain["kind"],
  "rational", cqRationalRawNullSpace[matrix, False, True],
  "coefficient", cqCoefficientRawNullSpace[matrix, False, True],
  _, cqFailure["MalformedSerialization", <|"domain" -> domain|>]
];

cqRawDomainElementZeroQ[domain_Association] := Switch[
  domain["kind"],
  "rational", Function[value, value === 0],
  "coefficient", Function[value, cqCoefficientZeroQ[value]],
  _, Function[value, False]
];

(* This is the one shared block-level common-kernel algorithm.  Scalar and
   matrix arithmetic remain specialized by domain; only control flow is shared. *)
cqRawCommonKernel[
  domain_, coordinateDimension_, rawConstraints_
] := Module[
  {basis, iterationNullities, constraint, preparedRows, rowStart,
   targetDenseRows, blockInfo, block, restricted, kernel, updatedBasis,
   previousBasisColumns, rankGain, potentialRankGain, constraintIndex,
   residual, nullity, zeroElementQ, identityBasis = True},
  residual = cqValidateRawDomain[domain, coordinateDimension, rawConstraints];
  If[FailureQ[residual], Return[residual]];
  basis = cqRawDomainIdentity[domain, coordinateDimension];
  If[FailureQ[basis], Return[basis]];
  iterationNullities = {coordinateDimension};
  zeroElementQ = cqRawDomainElementZeroQ[domain];
  For[constraintIndex = 1,
      constraintIndex <= Length[rawConstraints],
      constraintIndex++,
    If[basis["columns"] === 0,
      AppendTo[iterationNullities, 0],
      constraint = rawConstraints[[constraintIndex]];
      preparedRows = cqAdaptivePrepareRows[constraint["data"], zeroElementQ];
      rowStart = 1;
      targetDenseRows = $cqAdaptiveInitialDenseRows;
      While[rowStart <= Length[preparedRows] && basis["columns"] > 0,
        previousBasisColumns = basis["columns"];
        blockInfo = cqAdaptiveTakeRowBlock[
          preparedRows,
          rowStart,
          previousBasisColumns,
          constraint["columns"],
          targetDenseRows,
          zeroElementQ
        ];
        block = cqRawDomainMatrix[
          domain,
          blockInfo["row_count"],
          constraint["columns"],
          blockInfo["data"]
        ];
        If[FailureQ[block], Return[block]];
        restricted = If[
          identityBasis,
          block,
          cqRawDomainMultiply[domain, block, basis]
        ];
        If[FailureQ[restricted], Return[restricted]];
        If[cqRawDomainZeroQ[domain, restricted],
          rankGain = 0,
          kernel = cqRawDomainNullSpace[domain, restricted];
          If[FailureQ[kernel], Return[kernel]];
          updatedBasis = If[
            identityBasis,
            kernel["basis_matrix"],
            cqRawDomainMultiply[domain, basis, kernel["basis_matrix"]]
          ];
          If[FailureQ[updatedBasis], Return[updatedBasis]];
          If[updatedBasis["columns"] > previousBasisColumns,
            Return[cqFailure["ResidualVerificationFailed", <||>]]
          ];
          rankGain = previousBasisColumns - updatedBasis["columns"];
          basis = updatedBasis;
          If[rankGain > 0, identityBasis = False]
        ];
        potentialRankGain = Min[blockInfo["row_count"], previousBasisColumns];
        targetDenseRows = cqAdaptiveNextDenseRows[
          targetDenseRows, rankGain, potentialRankGain
        ];
        rowStart = blockInfo["next_start"]
      ];
      AppendTo[iterationNullities, basis["columns"]]
    ]
  ];
  For[constraintIndex = 1,
      constraintIndex <= Length[rawConstraints],
      constraintIndex++,
    residual = cqRawDomainMultiply[
      domain, rawConstraints[[constraintIndex]], basis
    ];
    If[FailureQ[residual] || !cqRawDomainZeroQ[domain, residual],
      Return[cqFailure["ResidualVerificationFailed", <||>]]
    ]
  ];
  nullity = basis["columns"];
  <|
    "type" -> "raw_common_kernel_result",
    "backend_kind" -> domain["kind"],
    "field_id" -> If[
      domain["kind"] === "coefficient",
      cqCoefficientFieldID[domain["context"]],
      None
    ],
    "coordinate_dimension" -> coordinateDimension,
    "constraint_count" -> Length[rawConstraints],
    "iteration_nullities" -> iterationNullities,
    "basis_matrix" -> basis,
    "rank" -> coordinateDimension - nullity,
    "nullity" -> nullity,
    "exact_residual_verified" -> True
  |>
];

cqCommonKernelFinalize[
  plan_, rawResult_, outputBasis_
] := Module[
  {context, backendKind, expectedField, outputRows, exactBasis, exactRows,
   coordinateDimension, constraintCount, nullity, rank, iterationNullities},
  If[!cqCommonKernelWorkPlanQ[plan] || !AssociationQ[rawResult] ||
     Lookup[rawResult, "type", None] =!= "raw_common_kernel_result" ||
     Lookup[rawResult, "exact_residual_verified", False] =!= True,
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  context = plan["public_context"];
  backendKind = plan["backend_kind"];
  coordinateDimension = Lookup[rawResult, "coordinate_dimension", None];
  constraintCount = Lookup[rawResult, "constraint_count", None];
  nullity = Lookup[rawResult, "nullity", None];
  rank = Lookup[rawResult, "rank", None];
  iterationNullities = Lookup[rawResult, "iteration_nullities", None];
  If[coordinateDimension =!= plan["coordinate_dimension"] ||
     constraintCount =!= plan["constraint_count"] ||
     Lookup[rawResult, "backend_kind", None] =!= backendKind ||
     Lookup[rawResult, "field_id", Missing["field_id"]] =!= If[
       backendKind === "coefficient",
       cqCoefficientFieldID[plan["field_context"]],
       None
     ] ||
     !IntegerQ[nullity] || nullity < 0 ||
     !IntegerQ[rank] || rank < 0 || rank + nullity =!= coordinateDimension ||
     !ListQ[iterationNullities] ||
     Length[iterationNullities] =!= constraintCount + 1 ||
     First[iterationNullities] =!= coordinateDimension ||
     Last[iterationNullities] =!= nullity ||
     !And @@ MapThread[GreaterEqual, {Most[iterationNullities], Rest[iterationNullities]}],
    Return[cqFailure["MalformedSerialization", <|"raw_result" -> rawResult|>]]
  ];
  If[Lookup[outputBasis, "rows", None] =!= coordinateDimension ||
     Lookup[outputBasis, "columns", None] =!= nullity,
    Return[cqFailure["DimensionMismatch", <||>]]
  ];
  Switch[backendKind,
    "rational",
      If[!cqRationalRawMatrixQ[outputBasis],
        Return[cqFailure["MalformedSerialization", <||>]]
      ];
      outputRows = cqRationalRawTranspose[outputBasis];
      If[FailureQ[outputRows], Return[outputRows]];
      exactBasis = cqRationalRawToExact[context, outputBasis];
      exactRows = cqRationalRawToExact[context, outputRows],
    "coefficient",
      expectedField = cqCoefficientFieldFromCyclotomic[context];
      If[FailureQ[expectedField], Return[expectedField]];
      If[!cqCoefficientRawMatrixQ[outputBasis] ||
         cqCoefficientFieldID[outputBasis["context"]] =!=
           cqCoefficientFieldID[expectedField],
        Return[cqFailure["ConductorMismatch", <||>]]
      ];
      outputRows = cqCoefficientRawTranspose[outputBasis];
      If[FailureQ[outputRows], Return[outputRows]];
      exactBasis = cqCoefficientRawToExact[context, outputBasis];
      exactRows = cqCoefficientRawToExact[context, outputRows],
    _, Return[cqFailure["MalformedSerialization", <||>]]
  ];
  If[AnyTrue[{exactBasis, exactRows}, FailureQ],
    Return[FirstCase[{exactBasis, exactRows}, _Failure]]
  ];
  <|
    "type" -> "common_kernel_result",
    "context" -> context,
    "conductor" -> context["conductor"],
    "coordinate_dimension" -> rawResult["coordinate_dimension"],
    "constraint_count" -> rawResult["constraint_count"],
    "method" -> plan["method"],
    "iteration_nullities" -> rawResult["iteration_nullities"],
    "basis_matrix" -> exactBasis,
    "nullspace_rows" -> exactRows,
    "rank" -> rawResult["rank"],
    "nullity" -> rawResult["nullity"],
    "exact_residual_verified" -> rawResult["exact_residual_verified"]
  |>
];

cqRationalCommonKernel[compiled_Association] := Module[
  {rawConstraints, plan, domain, rawResult},
  rawConstraints = (cqRationalRawFromExact[#, False] &) /@
    compiled["constraints"];
  If[AnyTrue[rawConstraints, FailureQ],
    Return[FirstCase[rawConstraints, _Failure]]
  ];
  plan = cqRationalWorkPlan[compiled, rawConstraints];
  If[!cqCommonKernelWorkPlanQ[plan],
    Return[cqFailure["MalformedSerialization", <|"work_plan" -> plan|>]]
  ];
  domain = cqWorkPlanDomain[plan];
  rawResult = cqRawCommonKernel[
    domain, plan["coordinate_dimension"], plan["raw_constraints"]
  ];
  If[FailureQ[rawResult], Return[rawResult]];
  cqCommonKernelFinalize[plan, rawResult, rawResult["basis_matrix"]]
];

cqCoefficientCommonKernel[compiled_Association] := Module[
  {fullRawConstraints, plan, rawResult, outputBasis},
  (* This is the sole Exact -> full coefficient conversion.  A real-subfield
     plan projects these same vectors; it never revisits Exact Matrix objects. *)
  fullRawConstraints = (cqCoefficientRawFromExact[#, False] &) /@
    compiled["constraints"];
  If[AnyTrue[fullRawConstraints, FailureQ],
    Return[FirstCase[fullRawConstraints, _Failure]]
  ];
  plan = cqCoefficientWorkPlan[
    compiled["context"], compiled["coordinate_dimension"], fullRawConstraints
  ];
  If[FailureQ[plan], Return[plan]];
  If[!cqCommonKernelWorkPlanQ[plan],
    Return[cqFailure["MalformedSerialization", <|"work_plan" -> plan|>]]
  ];
  rawResult = cqRawCommonKernel[
    cqWorkPlanDomain[plan],
    plan["coordinate_dimension"],
    plan["raw_constraints"]
  ];
  If[FailureQ[rawResult], Return[rawResult]];
  outputBasis = cqCoefficientWorkPlanLiftBasis[
    plan, rawResult["basis_matrix"]
  ];
  If[FailureQ[outputBasis], Return[outputBasis]];
  cqCommonKernelFinalize[plan, rawResult, outputBasis]
];

CyclotomicCommonKernel[compiled_Association] := Module[{validated},
  validated = cqValidateCompiledProblem[compiled];
  If[FailureQ[validated], Return[validated]];
  If[compiled["context"]["conductor"] === 1,
    cqRationalCommonKernel[compiled],
    cqCoefficientCommonKernel[compiled]
  ]
];

CyclotomicCommonKernel[input_] :=
  cqFailure["MalformedSerialization", <|"input" -> HoldForm[input]|>];
