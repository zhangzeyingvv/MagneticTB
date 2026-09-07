(* ::Package:: *)

BeginPackage["DirectProductRepresentation`"]

CompileDirectProductRepresentation::usage =
  "CompileDirectProductRepresentation[images,localMatrices,group,flags] constructs a direct-product block-monomial representation without site coordinates.";

Begin["`Private`"]

CompileDirectProductRepresentation::data =
  "Expected one operation-by-source image table and one local matrix per operation for every orbit.";

CompileDirectProductRepresentation[
    imageSiteIndices_List,
    localRepresentations_List,
    group_Association,
    antiunitaryFlags_List
  ] := Module[{orbitCount, operationCount, localDimensions, localBlocks, assembled},
  orbitCount = Length[imageSiteIndices];
  operationCount = Length[antiunitaryFlags];
  If[orbitCount < 1 || !GroupAlgebra`GroupAlgebraQ[group] ||
      group["Order"] =!= operationCount ||
      Length[localRepresentations] =!= orbitCount ||
      !And @@ (ListQ[#] && Length[#] === operationCount & /@ localRepresentations) ||
      !And @@ (GroupAlgebra`GroupActionTableQ[group, #] & /@ imageSiteIndices) ||
      !And @@ Map[MatrixPredicates`ExactSquareMatrixQ, localRepresentations, {2}] ||
      !And @@ (SameQ @@ (Dimensions /@ #) & /@ localRepresentations) ||
      !And @@ Map[MatrixPredicates`ExactUnitaryMatrixQ, localRepresentations, {2}],
    Message[CompileDirectProductRepresentation::data]; Return[$Failed]
  ];
  localDimensions = Length[First[#]] & /@ localRepresentations;
  localBlocks = Table[
    ConstantArray[localRepresentations[[orbit, operation]],
      Length[imageSiteIndices[[orbit, operation]]]],
    {orbit, orbitCount}, {operation, operationCount}
  ];
  assembled = RepresentationData`AssembleBlockMonomialRepresentation[
    imageSiteIndices, localBlocks, localDimensions, antiunitaryFlags];
  If[assembled === $Failed, Return[$Failed]];
  Join[
    <|
      "Method" -> "DirectProduct",
      "Convention" -> "orbit-major, site-major, local-index-minor",
      "GroupAlgebra" -> group,
      "AntiunitaryFlags" -> antiunitaryFlags,
      "LocalRepresentationMatrices" -> localRepresentations,
      "LocalBlocks" -> localBlocks,
      "LocalDimensions" -> localDimensions
    |>,
    assembled
  ]
];
CompileDirectProductRepresentation[___] :=
  (Message[CompileDirectProductRepresentation::data]; $Failed);

End[]

EndPackage[]
