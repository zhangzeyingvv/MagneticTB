(* ::Package:: *)

BeginPackage["MatrixPredicates`"]

ExactSquareMatrixQ::usage = "ExactSquareMatrixQ[matrix] tests whether matrix is a nonempty square matrix.";
ExactZeroExpressionQ::usage = "ExactZeroExpressionQ[expression] tests exact symbolic zero without a numerical tolerance.";
ExactZeroMatrixQ::usage = "ExactZeroMatrixQ[matrix] tests whether every matrix entry is exactly zero.";
ExactMatrixEqualQ::usage = "ExactMatrixEqualQ[first,second] tests exact symbolic matrix equality.";
ExactUnitaryMatrixQ::usage = "ExactUnitaryMatrixQ[matrix] tests exact unitarity.";
ExactHermitianMatrixQ::usage = "ExactHermitianMatrixQ[matrix] tests exact Hermiticity.";
IntegerVectorQ::usage = "IntegerVectorQ[vector] tests whether every component is exactly integral.";

Begin["`Private`"]

ExactSquareMatrixQ[matrix_] :=
  MatrixQ[matrix] && Length[matrix] > 0 &&
    Length[Dimensions[matrix]] === 2 && SameQ @@ Dimensions[matrix];

ExactZeroExpressionQ[expression_] :=
  TrueQ[expression === 0] ||
    TrueQ[PossibleZeroQ[FullSimplify[expression]]];

ExactZeroMatrixQ[matrix_?MatrixQ] :=
  And @@ (ExactZeroExpressionQ /@ Flatten[Normal[matrix]]);
ExactZeroMatrixQ[_] := False;

ExactMatrixEqualQ[first_?MatrixQ, second_?MatrixQ] :=
  Dimensions[first] === Dimensions[second] &&
    ExactZeroMatrixQ[FullSimplify[first - second]];
ExactMatrixEqualQ[_, _] := False;

ExactUnitaryMatrixQ[matrix_] := Module[{dimension},
  If[!ExactSquareMatrixQ[matrix], Return[False]];
  dimension = Length[matrix];
  ExactZeroMatrixQ[
    FullSimplify[
      ConjugateTranspose[matrix].matrix - IdentityMatrix[dimension]
    ]
  ]
];

ExactHermitianMatrixQ[matrix_] :=
  ExactSquareMatrixQ[matrix] &&
    ExactMatrixEqualQ[matrix, ConjugateTranspose[matrix]];

IntegerVectorQ[vector_List] :=
  And @@ (
    TrueQ[Element[FullSimplify[#], Integers]] ||
      IntegerQ[FullSimplify[#]] & /@ vector
  );
IntegerVectorQ[_] := False;

End[]

EndPackage[]
