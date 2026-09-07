(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

Options[symhamII] = {"wcc" -> None};

symhamII[ham_, OptionsPattern[]] := Module[{U, hII, wc},
  If[OptionValue["wcc"] === None, wc = wcc, wc = OptionValue["wcc"]];
  U = DiagonalMatrix[
    Table[Exp[-I {kx, ky, kz} . tau], {tau, wc}]
  ];
  hII = ComplexExpand[ConjugateTranspose[U]] . ham . U;
  Expand[TrigToExp@FullSimplify@hII]
];

End[]
EndPackage[]
