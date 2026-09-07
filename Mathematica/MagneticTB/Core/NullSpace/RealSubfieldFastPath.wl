(* Exact adapter for the maximal real subfield Q(zeta_N + zeta_N^-1).

   The adapter never runs a null-space algorithm.  It only constructs the
   smaller coefficient field, projects full-field constraints into it, and
   lifts a resulting basis back to the public full cyclotomic context. *)

$cqUseRealSubfield = True;
$cqRealSubfieldDescriptorCache = <||>;

cqRealSubfieldDescriptor[fullContext_] := Module[
  {key, cached, conductor, fullDegree, realDegree, fullField, zeta, zetaInverse, generator,
   powers, exponent, embedding, rowSelectionMatrix, rowSelectionRREF,
   coordinateRows, coordinateMatrix, inverseAugmented, inverseRREF,
   coordinateInverse, relationCoefficients, realModulus, realContext, descriptor},
  If[!cqContextQ[fullContext],
    Return[cqFailure["ConductorMismatch", <||>]]
  ];
  conductor = fullContext["conductor"];
  fullDegree = fullContext["degree"];
  If[conductor <= 2 || !EvenQ[fullDegree],
    Return[Missing["RealSubfieldNotApplicable"]]
  ];
  key = fullContext["context_id"];
  cached = Lookup[$cqRealSubfieldDescriptorCache, key, Missing["not_cached"]];
  If[!MissingQ[cached], Return[cached]];
  realDegree = Quotient[fullDegree, 2];
  fullField = cqCoefficientFieldFromCyclotomic[fullContext];
  If[FailureQ[fullField], Return[fullField]];

  zeta = cqRootOfUnityElement[fullContext, conductor, 1];
  If[FailureQ[zeta], Return[zeta]];
  zetaInverse = cqRootOfUnityElement[fullContext, conductor, -1];
  If[FailureQ[zetaInverse], Return[zetaInverse]];
  generator = zeta["coefficients"] + zetaInverse["coefficients"];

  powers = ConstantArray[{}, realDegree + 1];
  powers[[1]] = cqCoefficientOne[fullField];
  For[exponent = 1, exponent <= realDegree, exponent++,
    powers[[exponent + 1]] = cqCoefficientMultiply[
      fullField, powers[[exponent]], generator
    ];
    If[FailureQ[powers[[exponent + 1]]], Return[powers[[exponent + 1]]]]
  ];
  embedding = Transpose[Take[powers, realDegree]];

  (* Pivot columns of Transpose[E] select r independent coefficient rows of E.
     Inverting only that r by r coordinate minor is substantially smaller than
     reducing [E | I_d].  Full-vector equality below remains the exact test. *)
  rowSelectionMatrix = cqRationalRawMatrix[
    realDegree, fullDegree, Transpose[embedding]
  ];
  If[FailureQ[rowSelectionMatrix], Return[rowSelectionMatrix]];
  rowSelectionRREF = cqRationalRawRREF[rowSelectionMatrix];
  If[FailureQ[rowSelectionRREF], Return[rowSelectionRREF]];
  If[rowSelectionRREF["rank"] =!= realDegree,
    Return[cqFailure["RealSubfieldConstructionFailed", <|
      "conductor" -> conductor, "stage" -> "coordinate_rows"
    |>]]
  ];
  coordinateRows = 1 + rowSelectionRREF["pivot_columns"];
  coordinateMatrix = embedding[[coordinateRows]];
  inverseAugmented = cqRationalRawMatrix[
    realDegree,
    2 realDegree,
    MapThread[Join, {coordinateMatrix, IdentityMatrix[realDegree]}]
  ];
  If[FailureQ[inverseAugmented], Return[inverseAugmented]];
  inverseRREF = cqRationalRawRREF[inverseAugmented];
  If[FailureQ[inverseRREF], Return[inverseRREF]];
  If[inverseRREF["rank"] =!= realDegree ||
     Take[inverseRREF["pivot_columns"], realDegree] =!= Range[0, realDegree - 1],
    Return[cqFailure["RealSubfieldConstructionFailed", <|
      "conductor" -> conductor, "stage" -> "coordinate_inverse"
    |>]]
  ];
  coordinateInverse = inverseRREF["reduced"]["data"][[
    1 ;; realDegree, realDegree + 1 ;; 2 realDegree
  ]];
  If[coordinateInverse . coordinateMatrix =!= IdentityMatrix[realDegree],
    Return[cqFailure["RealSubfieldConstructionFailed", <|
      "conductor" -> conductor, "stage" -> "coordinate_inverse_verification"
    |>]]
  ];

  relationCoefficients = coordinateInverse . (-powers[[realDegree + 1, coordinateRows]]);
  realModulus = Join[relationCoefficients, {1}];
  If[embedding . relationCoefficients =!= -powers[[realDegree + 1]],
    Return[cqFailure["RealSubfieldConstructionFailed", <|
      "conductor" -> conductor, "stage" -> "minimal_polynomial_verification"
    |>]]
  ];

  realContext = <|
    "type" -> "coefficient_field_context",
    "field_kind" -> "maximal_real_subfield",
    "field_id" -> StringJoin["real_subfield:", ToString[conductor]],
    "source_conductor" -> conductor,
    "degree" -> realDegree,
    "defining_polynomial" -> realModulus
  |>;
  If[!cqCoefficientFieldContextQ[realContext],
    Return[cqFailure["RealSubfieldConstructionFailed", <|
      "conductor" -> conductor, "stage" -> "internal_context"
    |>]]
  ];
  descriptor = <|
    "type" -> "real_subfield_descriptor",
    "full_context" -> fullContext,
    "real_context" -> realContext,
    "generator_coefficients" -> generator,
    "embedding_matrix" -> embedding,
    "coordinate_rows" -> coordinateRows,
    "coordinate_inverse" -> coordinateInverse,
    "minimal_polynomial" -> realModulus
  |>;
  $cqRealSubfieldDescriptorCache[key] = descriptor;
  descriptor
];

cqRealSubfieldProjectVector[descriptor_Association, coefficients_List] := Module[
  {coordinates},
  coordinates = descriptor["coordinate_inverse"] .
    coefficients[[descriptor["coordinate_rows"]]];
  If[descriptor["embedding_matrix"] . coordinates === coefficients,
    coordinates,
    Missing["NotInRealSubfield"]
  ]
];

cqRealSubfieldLiftVector[descriptor_Association, coordinates_List] :=
  descriptor["embedding_matrix"] . coordinates;

cqRealSubfieldProjectRawConstraints[
  rawConstraints_List, descriptor_Association
] := Module[
  {realContext, project, converted = {}, matrix, data, projected,
   constraintIndex, notApplicableTag = Unique["notReal"]},
  realContext = descriptor["real_context"];
  project[value_List] := project[value] =
    cqRealSubfieldProjectVector[descriptor, value];
  For[constraintIndex = 1,
      constraintIndex <= Length[rawConstraints],
      constraintIndex++,
    matrix = rawConstraints[[constraintIndex]];
    data = Catch[
      Map[
        Function[value,
          projected = project[value];
          If[MissingQ[projected], Throw[projected, notApplicableTag]];
          projected
        ],
        matrix["data"],
        {2}
      ],
      notApplicableTag
    ];
    If[MissingQ[data], Return[Missing["RealSubfieldNotApplicable"]]];
    If[FailureQ[data], Return[data]];
    projected = cqCoefficientRawMatrix[
      realContext, matrix["rows"], matrix["columns"], data
    ];
    If[FailureQ[projected], Return[projected]];
    AppendTo[converted, projected]
  ];
  converted
];

cqCoefficientWorkPlan[
  fullContext_, coordinateDimension_Integer, fullRawConstraints_List
] := Module[
  {fullField, descriptor, realRawConstraints, fieldContext, rawConstraints,
   liftKind, liftMatrix, method},
  fullField = cqCoefficientFieldFromCyclotomic[fullContext];
  If[FailureQ[fullField], Return[fullField]];
  fieldContext = fullField;
  rawConstraints = fullRawConstraints;
  liftKind = "identity";
  liftMatrix = IdentityMatrix[fullContext["degree"]];
  method = "iterative";
  If[TrueQ[$cqUseRealSubfield] && fullContext["conductor"] > 2,
    descriptor = cqRealSubfieldDescriptor[fullContext];
    If[FailureQ[descriptor], Return[descriptor]];
    If[!MissingQ[descriptor],
      realRawConstraints = cqRealSubfieldProjectRawConstraints[
        fullRawConstraints, descriptor
      ];
      If[FailureQ[realRawConstraints], Return[realRawConstraints]];
      If[!MissingQ[realRawConstraints],
        fieldContext = descriptor["real_context"];
        rawConstraints = realRawConstraints;
        liftKind = "linear";
        liftMatrix = descriptor["embedding_matrix"];
        method = "iterative_real_subfield"
      ]
    ]
  ];
  <|
    "type" -> "common_kernel_work_plan",
    "backend_kind" -> "coefficient",
    "public_context" -> fullContext,
    "field_context" -> fieldContext,
    "coordinate_dimension" -> coordinateDimension,
    "constraint_count" -> Length[rawConstraints],
    "raw_constraints" -> rawConstraints,
    "lift_kind" -> liftKind,
    "lift_matrix" -> liftMatrix,
    "method" -> method
  |>
];

cqCoefficientWorkPlanLiftBasis[plan_Association, basis_Association] := Module[
  {fullField, data},
  If[!cqCommonKernelWorkPlanQ[plan] ||
     !cqCoefficientRawMatrixQ[basis] ||
     cqCoefficientFieldID[basis["context"]] =!=
       cqCoefficientFieldID[plan["field_context"]],
    Return[cqFailure["MalformedSerialization", <||>]]
  ];
  fullField = cqCoefficientFieldFromCyclotomic[plan["public_context"]];
  If[FailureQ[fullField], Return[fullField]];
  If[plan["lift_kind"] === "identity" &&
     cqCoefficientFieldID[basis["context"]] === cqCoefficientFieldID[fullField],
    Return[basis]
  ];
  data = Map[plan["lift_matrix"] . # &, basis["data"], {2}];
  cqCoefficientRawMatrix[fullField, basis["rows"], basis["columns"], data]
];
