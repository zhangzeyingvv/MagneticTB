(* ::Package:: *)

BeginPackage["MagneticTB`"];

CompileCyclotomicMatrices::usage =
  "CompileCyclotomicMatrices[input] compiles structured root_sum constraint matrices into one exact cyclotomic context.";
CyclotomicRREF::usage =
  "CyclotomicRREF[matrix] computes deterministic exact reduced row echelon form without System`RowReduce.";
CyclotomicNullSpace::usage =
  "CyclotomicNullSpace[matrix] computes an exact null-space basis without System`NullSpace.";
CyclotomicCommonKernel::usage =
  "CyclotomicCommonKernel[compiled] computes the iterative common kernel of all compiled constraint blocks.";
CyclotomicCommonNullSpace::usage =
  "CyclotomicCommonNullSpace[{matrix1, matrix2, ...}] accepts ordinary exact Mathematica matrices and returns row-basis vectors for their common null space.";
RestoreCyclotomicExpression::usage =
  "RestoreCyclotomicExpression[object] restores matrices/results; RestoreCyclotomicExpression[context, element] restores one cyclotomic element.";

Begin["`Private`"];

$cyclotomicNullSpaceDirectory = DirectoryName[$InputFileName];

Scan[
  Get[FileNameJoin[{$cyclotomicNullSpaceDirectory, #}]] &,
  {
    "RationalArithmetic.wl",
    "PolynomialArithmetic.wl",
    "CyclotomicField.wl",
    "CyclotomicConversion.wl",
    "ExactMatrix.wl",
    "AdaptiveBlocking.wl",
    "RationalFastPath.wl",
    "CyclotomicLinearAlgebraFastPath.wl",
    "RealSubfieldFastPath.wl",
    "CommonKernel.wl",
    "ExactNullSpace.wl",
    "ExpressionFrontend.wl"
  }
];

End[];
EndPackage[];
