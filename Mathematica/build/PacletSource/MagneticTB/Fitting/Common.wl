(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[
  normalizeFittingHamiltonian,
  validFittingHamiltonianQ,
  validEigenvalueDataQ,
  fittingParameterSymbols,
  fittingNumericEigenvalues,
  fittingRulesList
];

normalizeFittingHamiltonian[hamiltonian_] :=
  normalizeBandHamiltonian[hamiltonian];

validFittingHamiltonianQ[hamiltonian_] := Module[{dimensions},
  dimensions = Quiet@Check[Dimensions[hamiltonian], Missing[]];
  MatrixQ[hamiltonian] &&
    MatchQ[dimensions, {_Integer, _Integer}] &&
    dimensions[[1]] > 0 &&
    Equal @@ dimensions
];

validEigenvalueDataQ[data_] :=
  ListQ[data] && data =!= {} &&
    And @@ Map[
      Function[record,
        MatchQ[record, {_List, _List}] &&
          Length[record[[1]]] === 3 &&
          VectorQ[record[[1]], NumericQ] &&
          record[[2]] =!= {} &&
          VectorQ[record[[2]], NumericQ]
      ],
      data
    ] &&
    SameQ @@ (Length /@ data[[All, 2]]);

fittingParameterSymbols[hamiltonian_] :=
  bandParameterSymbols[TrigToExp[hamiltonian]];

fittingNumericEigenvalues[matrix_] :=
  Sort@Chop@Eigenvalues[N[matrix]];

fittingRulesList[rules_] := Which[
  AssociationQ[rules], Normal[rules],
  MatchQ[rules, {(_Rule | _RuleDelayed) ...}], rules,
  True, $Failed
];

End[]

EndPackage[]
