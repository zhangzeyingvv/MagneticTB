(* ::Package:: *)

BeginPackage["MagneticTB`"]


Begin["`Private`"]


texOutput[mat_]:=Module[{ptex},
ptex=ToString[TeXForm@mat];
StringReplace[ToString[ptex],{
RegularExpression["\\\\text\\{(\\D)(\\d)\\}"]->"$1"<>"_"<>"{$2}",
RegularExpression["\\\\text\\{(k)(.)\\}"]->"$1"<>"_"<>"$2",
RegularExpression["\\\\text\\{(\\D)(\\d+)(\\D)(\\d+)\\}"]->"$1"<>"_"<>"{$4}"<>"^"<>"{$2}"
}]]


braLatt=<|"CubicP" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}, 
    {{a, 0, 0}, {0, a, 0}, {0, 0, a}}}, 
  "CubicF" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}, 
    {{0, a/2, a/2}, {a/2, 0, a/2}, {a/2, a/2, 0}}}, 
  "CubicI" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}, 
    {{-a/2, a/2, a/2}, {a/2, -a/2, a/2}, {a/2, a/2, -a/2}}}, 
  "TetragonalP" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, c}}, 
    {{a, 0, 0}, {0, a, 0}, {0, 0, c}}}, "TetragonalI" -> 
   {{{a, 0, 0}, {0, a, 0}, {0, 0, c}}, {{-a/2, a/2, c/2}, {a/2, -a/2, c/2}, 
     {a/2, a/2, -c/2}}}, "OrthorhombicP" -> 
   {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, {{a, 0, 0}, {0, b, 0}, {0, 0, c}}}, 
  "OrthorhombicF" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, 
    {{0, b/2, c/2}, {a/2, 0, c/2}, {a/2, b/2, 0}}}, 
  "OrthorhombicI" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, 
    {{-a/2, b/2, c/2}, {a/2, -b/2, c/2}, {a/2, b/2, -c/2}}}, 
  "OrthorhombicC" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, 
    {{a/2, b/2, 0}, {-a/2, b/2, 0}, {0, 0, c}}}, 
  "HexagonalP" -> {{{a, 0, 0}, {-a/2, (Sqrt[3]*a)/2, 0}, {0, 0, c}}, 
    {{a, 0, 0}, {-a/2, (Sqrt[3]*a)/2, 0}, {0, 0, c}}}, 
  "TrigonalR"->{{{Sqrt[3] a,0,0},{-((Sqrt[3] a)/2),(3 a)/2,0},{0,0,3 c}},
  {{(Sqrt[3] a)/2,a/2,c},{-((Sqrt[3] a)/2),a/2,c},{0,-a,c}}},
(*  "TrigonalR" -> {{{a, 0, 0}, {a*Cos[\[Alpha]], a*Sin[\[Alpha]], 0}, 
     {a*Cos[\[Alpha]], a*(Cos[\[Alpha]] - Cos[\[Alpha]]^2)*Csc[\[Alpha]], 
      a*Sqrt[1 - 3*Cos[\[Alpha]]^2 + 2*Cos[\[Alpha]]^3]*Csc[\[Alpha]]}}, 
    {{a, 0, 0}, {a*Cos[\[Alpha]], a*Sin[\[Alpha]], 0}, 
     {a*Cos[\[Alpha]], a*(Cos[\[Alpha]] - Cos[\[Alpha]]^2)*Csc[\[Alpha]], 
      a*Sqrt[1 - 3*Cos[\[Alpha]]^2 + 2*Cos[\[Alpha]]^3]*Csc[\[Alpha]]}}}, *)
  "MonoclinicP" -> {{{a, 0, 0}, {0, b, 0}, {c*Cos[\[Beta]], 0, 
      c*Sqrt[1 - Cos[\[Beta]]^2]}}, {{a, 0, 0}, {0, b, 0}, 
     {c*Cos[\[Beta]], 0, c*Sqrt[1 - Cos[\[Beta]]^2]}}}, 
  "MonoclinicB" -> {{{a, 0, 0}, {0, b, 0}, {0, c*Cos[\[Beta]], 
      c*Sin[\[Beta]]}}, {{a/2, b/2, 0}, {-a/2, b/2, 0}, 
     {0, c*Cos[\[Beta]], c*Sin[\[Beta]]}}}, 
  "TriclinicP" -> {{{a, 0, 0}, {b*Cos[\[Gamma]], b*Sin[\[Gamma]], 0}, 
     {c*Cos[\[Beta]], c*(Cos[\[Alpha]] - Cos[\[Beta]]*Cos[\[Gamma]])*
       Csc[\[Gamma]], c*Sqrt[1 - Cos[\[Alpha]]^2 - Cos[\[Beta]]^2 + 
         2*Cos[\[Alpha]]*Cos[\[Beta]]*Cos[\[Gamma]] - Cos[\[Gamma]]^2]*
       Csc[\[Gamma]]}}, {{a, 0, 0}, {b*Cos[\[Gamma]], b*Sin[\[Gamma]], 0}, 
     {c*Cos[\[Beta]], c*(Cos[\[Alpha]] - Cos[\[Beta]]*Cos[\[Gamma]])*
       Csc[\[Gamma]], c*Sqrt[1 - Cos[\[Alpha]]^2 - Cos[\[Beta]]^2 + 
         2*Cos[\[Alpha]]*Cos[\[Beta]]*Cos[\[Gamma]] - Cos[\[Gamma]]^2]*
       Csc[\[Gamma]]}}}|>;

msgop[number_]:=Module[{data},
	data=MSGOP[number];
	Print["Magnetic space group (BNS): ",data["MSG"]];
	Print["Lattice: ",data["BRAV"]];
	Print["Primitive Lattice Vactor: ", braLatt[data["BRAV"]][[2]]];
	Print["Conventional Lattice Vactor: ", braLatt[data["BRAV"]][[1]]];
	data["SymmetryOperation"]
];



symhamII[ham_]:=Module[{U,hII},
U=DiagonalMatrix[Table[Exp[-I{kx,ky,kz} . tau],{tau,wcc}]];
hII=ComplexExpand[ConjugateTranspose[U]] . ham . U;
Expand[TrigToExp@FullSimplify@hII]
];

symmetryopsII:=Module[{tmp},

Table[
tmp=symminfo[[i,2]];
symmetryops[[i]]Table[Exp[I {kx,ky,kz} . (
Inverse[tmp] . wcc[[k]]-wcc[[l]])],{k,Length[wcc]},{l,Length[wcc]}]
,{i,Length[symminfo]}]
];

compactForm:={#[[1]],#[[2]] . {"x","y","z"},#[[3]],#[[4]],MatrixForm[#[[5]]],
MatrixForm[#[[6]]]}&/@((symminfo\[Transpose]~Join~{symmetryops}~Join~{symmetryopsII})\[Transpose])


End[]
EndPackage[]

