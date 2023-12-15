(* ::Package:: *)

BeginPackage["MagneticTB`"]





Begin["`Private`"]


directSum =(*FullSimplify@*)ArrayFlatten[{{#1, 0}, {0, #2}}] &;


spinMatrix2[m_] := Module[
  (*https://mathematica.stackexchange.com/questions/29924/axis-angle-
  from-rotation-matrix*)
  {ang, axis, ovec, nn, nvec, rm, s, w, w1, wm, xx, yy, zz, mat, 
   axisAngle},
  If[FullSimplify[Det[m]] == -1, mat = -m, mat = m];
  axis = {mat[[3, 2]] - mat[[2, 3]], mat[[1, 3]] - mat[[3, 1]], 
    mat[[2, 1]] - mat[[1, 2]]};
  nn = Simplify[Norm[axis]];
  If[nn == 0,
   ang = \[Pi] Boole[Total[Diagonal[mat]] < 3];
   rm = (mat + IdentityMatrix[3])/2;
   axis = Normalize[Extract[rm, Ordering[Max /@ Abs[rm], -1]]],
   
   {xx, yy, zz} = Simplify[axis/nn];
   s = 2 UnitStep[zz] - 1; w = -1/(s + zz); w1 = xx yy w;
   ovec = {1 + s w xx xx, s w1, -s xx};
   nvec = {w1, s + w yy yy, -yy};
   wm = mat . ovec;
   ang = Arg[Simplify[wm . ovec + I wm . nvec]]];
  (*axisAngle = {FullSimplify@ang, axis};*)
  ExpToTrig@
   MatrixExp[-I FullSimplify@ang Sum[PauliMatrix[i] FullSimplify[(Normalize@axis)][[i]], {i, 3}]/2]
  
  ];


texOutput[mat_]:=Module[{ptex},
ptex=ToString[TeXForm@mat];
StringReplace[ToString[ptex],{
RegularExpression["\\\\text\\{(\\D)(\\d+)\\}"]->"$1"<>"_"<>"{$2}",
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
  "OrthorhombicA" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, 
    {{a, 0, 0}, {0,-b/2, c/2}, {0, b/2, c/2}}},
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
  "MonoclinicP" -> {{{a,0,0},{0,0,b},{c Cos[\[Beta]],c Sin[\[Beta]],0}}, 
  {{a,0,0},{0,0,b},{c Cos[\[Beta]],c Sin[\[Beta]],0}}}, 
  "MonoclinicB" -> {{{a, 0, 0}, {0, 0, b}, {c*Cos[\[Beta]], 
      c*Sin[\[Beta]],0}}, {{a/2,  0, b/2}, {-a/2, 0, b/2}, 
     {c*Cos[\[Beta]], c*Sin[\[Beta]],0}}}, 
     
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



Options[symhamII] = {"wcc"->None};
symhamII[ham_, OptionsPattern[]]:=Module[{U,hII,wc},
If[OptionValue["wcc"]===None,wc=wcc,wc=OptionValue["wcc"]];
U=DiagonalMatrix[Table[Exp[-I{kx,ky,kz} . tau],{tau,wc}]];
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
MatrixForm[#[[6]]]}&/@((symminfo\[Transpose]~Join~{symmetryops}~Join~{symmetryopsII})\[Transpose]);

showbonds[n_]:=Module[{data,bondlength,grid,hasbond,bpos,head},
data=bondclassify[[n]];

bondlength=data[[1,1]];
head=PadRight[{ToString[n-1]<>"-th neighbour, "<>"Bond length = "<> ToString@bondlength},3,SpanFromLeft];
grid={head,{"Atom position","Num of bonds",DisplayForm@SubsuperscriptBox["p","k","'"]}};
hasbond=bondclassify[[n]][[;;,3]][[;;,1,1]];
Do[

If[MemberQ[hasbond,pos],
bpos=First@FirstPosition[hasbond,pos];
grid=Append[grid,{pos,data[[bpos,2]],data[[bpos,3]][[;;,2]]}];,
grid=Append[grid,{pos,0,None}];
];
,{pos,Flatten[atompos,1][[;;,1]]}];
Grid[grid,Frame->All]
]






GenerateGroup[gens_,identityElement_,multiply_]:=Module[{i,j,ng,MAXORDER=200,orders,subs,mlist,g1,g2,g3},ng=Length[gens];orders=subs=Range[ng]*0;
For[i=1,i<=ng,i++,mlist=FoldList[multiply,Table[gens[[i]],MAXORDER]];
orders[[i]]=FirstPosition[mlist,identityElement][[1]];
subs[[i]]=mlist[[;;orders[[i]]]];];
g1=Union@@subs;
g2=Union@@Table[multiply[g1[[i]],g1[[j]]],{i,Length[g1]},{j,Length[g1]}];
g3=Complement[g2,g1];g1=g2;
While[g3!={},g2=Union@@Table[multiply[g3[[i]],g1[[j]]],{i,Length[g3]},{j,Length[g1]}];
g3=Complement[g2,g1];g1=g2;];
g2=SortBy[g2,#!=identityElement&]];
getGenerator[groupele_,identityElement_,multiply_]:=Module[
{try,group,tmptry,tmpgroup,seteq,groupeleDisorder,gereratorlist},
(*Greedy Alg, may not give the minimal Generator set*)
seteq[s1_,s2_]:=SubsetQ[s1,s2]&&SubsetQ[s2,s1];
If[groupele=={identityElement},Return[groupele]];
gereratorlist=Table[try={};
(*groupeleDisorder=RandomSample[groupele];*)
groupeleDisorder=groupele;
Do[If[ele==identityElement,Continue[]];
group=GenerateGroup[try,identityElement,multiply];
tmptry=Append[try,ele];
(*Print[group];*)
tmpgroup=GenerateGroup[tmptry,identityElement,multiply];
If[Not@seteq[group,tmpgroup],try=tmptry];
If[seteq[groupele,group],Break[]];,{ele,groupeleDisorder}];
try,1];
First[SortBy[gereratorlist,Length[#]&]]];
findgenind[sgop_]:=Module[{opm,times,identityElement,gens},
(*opm=MapAt[Mod[#,1]&,sgop[[;;,2;;]],{;;,2}];
opm=MapAt[#/.{"T"->1,"F"->0}&,opm,{;;,3}];*)
opm=MapAt[#/.{"T"->1,"F"->0}&,sgop,{;;,4}];
opm=opm[[;;,{2,4}]];
(*times=FullSimplify@{#1[[1]].#2[[1]],Mod[#1[[1]].#2[[2]]+#1[[2]],1],Mod[(#1[[3]]+#2[[3]]),2]}&;
identityElement={IdentityMatrix[3],{0,0,0},0};*)
times=(*Full*)Simplify@{#1[[1]] . #2[[1]](*,Mod[#1[[1]].#2[[2]]+#1[[2]],1]*),Mod[(#1[[2]]+#2[[2]]),2]}&;
identityElement={IdentityMatrix[3](*,{0,0,0}*),0};
gens=getGenerator[opm,identityElement,times];
Flatten[FirstPosition[opm,#]&/@gens]
]


End[]
EndPackage[]

