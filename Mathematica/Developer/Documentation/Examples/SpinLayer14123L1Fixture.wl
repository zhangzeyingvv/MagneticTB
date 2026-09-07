(* Exact 14.1.2.3.L.1 data generated once from the public
   SpinLayerCorepresentations API.  Keeping this small immutable snapshot in
   the documentation makes the shipped tutorial runnable without requiring a
   second paclet. *)

<|
  "Name" -> "14.1.2.3.L.1",
  "MagneticType" -> "Collinear",
  "ContinuousSpinPart" -> <|
    "Group" -> "CInfinity",
    "ParameterizedElement" -> "c_theta",
    "ParameterDomain" -> "theta in [0,2 Pi)",
    "FullElementConvention" -> "f_i c_theta",
    "InfiniteQ" -> True
  |>,
  "FiniteRepresentatives" -> {
    <|
      "SpaceMatrix" -> IdentityMatrix[3],
      "Translation" -> {0, 0, 0},
      "SpinMatrix" -> IdentityMatrix[3],
      "AntiunitaryQ" -> False
    |>,
    <|
      "SpaceMatrix" -> IdentityMatrix[3],
      "Translation" -> {0, 0, 0},
      "SpinMatrix" -> DiagonalMatrix[{1, -1, -1}],
      "AntiunitaryQ" -> True
    |>,
    <|
      "SpaceMatrix" -> DiagonalMatrix[{1, -1, -1}],
      "Translation" -> {1/2, 1/2, 0},
      "SpinMatrix" -> DiagonalMatrix[{-1, -1, 1}],
      "AntiunitaryQ" -> True
    |>,
    <|
      "SpaceMatrix" -> DiagonalMatrix[{1, -1, -1}],
      "Translation" -> {1/2, 1/2, 0},
      "SpinMatrix" -> DiagonalMatrix[{-1, 1, -1}],
      "AntiunitaryQ" -> False
    |>,
    <|
      "SpaceMatrix" -> -IdentityMatrix[3],
      "Translation" -> {0, 0, 0},
      "SpinMatrix" -> IdentityMatrix[3],
      "AntiunitaryQ" -> False
    |>,
    <|
      "SpaceMatrix" -> -IdentityMatrix[3],
      "Translation" -> {0, 0, 0},
      "SpinMatrix" -> DiagonalMatrix[{1, -1, -1}],
      "AntiunitaryQ" -> True
    |>,
    <|
      "SpaceMatrix" -> DiagonalMatrix[{-1, 1, 1}],
      "Translation" -> {1/2, 1/2, 0},
      "SpinMatrix" -> DiagonalMatrix[{-1, -1, 1}],
      "AntiunitaryQ" -> True
    |>,
    <|
      "SpaceMatrix" -> DiagonalMatrix[{-1, 1, 1}],
      "Translation" -> {1/2, 1/2, 0},
      "SpinMatrix" -> DiagonalMatrix[{-1, 1, -1}],
      "AntiunitaryQ" -> False
    |>
  },
  "FiniteSpinorMatrices" -> {
    IdentityMatrix[2],
    IdentityMatrix[2],
    {{0, -1}, {1, 0}},
    {{0, -1}, {1, 0}},
    IdentityMatrix[2],
    IdentityMatrix[2],
    {{0, 1}, {-1, 0}},
    {{0, 1}, {-1, 0}}
  }
|>
