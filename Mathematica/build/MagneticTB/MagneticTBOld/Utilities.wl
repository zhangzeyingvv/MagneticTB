(* ::Package:: *)

BeginPackage["MagneticTBOld`"]





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


texOutputold[mat_]:=Module[{ptex},
ptex=ToString[TeXForm@mat];
StringReplace[ToString[ptex],{
RegularExpression["\\\\text\\{(\\D)(\\d+)\\}"]->"$1"<>"_"<>"{$2}",
RegularExpression["\\\\text\\{(k)(.)\\}"]->"$1"<>"_"<>"$2",
RegularExpression["\\\\text\\{(\\D)(\\d+)(\\D)(\\d+)\\}"]->"$1"<>"_"<>"{$4}"<>"^"<>"{$2}"
}]]


braLattold=<|"CubicP" -> {{{aold, 0, 0}, {0, aold, 0}, {0, 0, aold}}, 
    {{aold, 0, 0}, {0, aold, 0}, {0, 0, aold}}}, 
  "CubicF" -> {{{aold, 0, 0}, {0, aold, 0}, {0, 0, aold}}, 
    {{0, aold/2, aold/2}, {aold/2, 0, aold/2}, {aold/2, aold/2, 0}}}, 
  "CubicI" -> {{{aold, 0, 0}, {0, aold, 0}, {0, 0, aold}}, 
    {{-aold/2, aold/2, aold/2}, {aold/2, -aold/2, aold/2}, {aold/2, aold/2, -aold/2}}}, 
  "TetragonalP" -> {{{aold, 0, 0}, {0, aold, 0}, {0, 0, cold}}, 
    {{aold, 0, 0}, {0, aold, 0}, {0, 0, cold}}}, "TetragonalI" -> 
   {{{aold, 0, 0}, {0, aold, 0}, {0, 0, cold}}, {{-aold/2, aold/2, cold/2}, {aold/2, -aold/2, cold/2}, 
     {aold/2, aold/2, -cold/2}}}, "OrthorhombicP" -> 
   {{{aold, 0, 0}, {0, bold, 0}, {0, 0, cold}}, {{aold, 0, 0}, {0, bold, 0}, {0, 0, cold}}}, 
  "OrthorhombicF" -> {{{aold, 0, 0}, {0, bold, 0}, {0, 0, cold}}, 
    {{0, bold/2, cold/2}, {aold/2, 0, cold/2}, {aold/2, bold/2, 0}}}, 
  "OrthorhombicI" -> {{{aold, 0, 0}, {0, bold, 0}, {0, 0, cold}}, 
    {{-aold/2, bold/2, cold/2}, {aold/2, -bold/2, cold/2}, {aold/2, bold/2, -cold/2}}}, 
  "OrthorhombicC" -> {{{aold, 0, 0}, {0, bold, 0}, {0, 0, cold}}, 
    {{aold/2, bold/2, 0}, {-aold/2, bold/2, 0}, {0, 0, cold}}}, 
  "OrthorhombicA" -> {{{aold, 0, 0}, {0, bold, 0}, {0, 0, cold}}, 
    {{aold, 0, 0}, {0,-bold/2, cold/2}, {0, bold/2, cold/2}}},
  "HexagonalP" -> {{{aold, 0, 0}, {-aold/2, (Sqrt[3]*aold)/2, 0}, {0, 0, cold}}, 
    {{aold, 0, 0}, {-aold/2, (Sqrt[3]*aold)/2, 0}, {0, 0, cold}}}, 
  "TrigonalR"->{{{Sqrt[3] aold,0,0},{-((Sqrt[3] aold)/2),(3 aold)/2,0},{0,0,3 cold}},
  {{(Sqrt[3] aold)/2,aold/2,cold},{-((Sqrt[3] aold)/2),aold/2,cold},{0,-aold,cold}}},
(*  "TrigonalR" -> {{{a, 0, 0}, {a*Cos[\[Alpha]], a*Sin[\[Alpha]], 0}, 
     {a*Cos[\[Alpha]], a*(Cos[\[Alpha]] - Cos[\[Alpha]]^2)*Csc[\[Alpha]], 
      a*Sqrt[1 - 3*Cos[\[Alpha]]^2 + 2*Cos[\[Alpha]]^3]*Csc[\[Alpha]]}}, 
    {{a, 0, 0}, {a*Cos[\[Alpha]], a*Sin[\[Alpha]], 0}, 
     {a*Cos[\[Alpha]], a*(Cos[\[Alpha]] - Cos[\[Alpha]]^2)*Csc[\[Alpha]], 
      a*Sqrt[1 - 3*Cos[\[Alpha]]^2 + 2*Cos[\[Alpha]]^3]*Csc[\[Alpha]]}}}, *)
  "MonoclinicP" -> {{{aold,0,0},{0,0,bold},{cold Cos[betaold],cold Sin[betaold],0}}, 
  {{aold,0,0},{0,0,bold},{cold Cos[betaold],cold Sin[betaold],0}}}, 
  "MonoclinicB" -> {{{aold, 0, 0}, {0, 0, bold}, {cold*Cos[betaold], 
      cold*Sin[betaold],0}}, {{aold/2,  0, bold/2}, {-aold/2, 0, bold/2}, 
     {cold*Cos[betaold], cold*Sin[betaold],0}}}, 
     
  "TriclinicP" -> {{{aold, 0, 0}, {bold*Cos[gammaold], bold*Sin[gammaold], 0}, 
     {cold*Cos[betaold], cold*(Cos[alphaold] - Cos[betaold]*Cos[gammaold])*
       Csc[gammaold], cold*Sqrt[1 - Cos[alphaold]^2 - Cos[betaold]^2 + 
         2*Cos[alphaold]*Cos[betaold]*Cos[gammaold] - Cos[gammaold]^2]*
       Csc[gammaold]}}, {{aold, 0, 0}, {bold*Cos[gammaold], bold*Sin[gammaold], 0}, 
     {cold*Cos[betaold], cold*(Cos[alphaold] - Cos[betaold]*Cos[gammaold])*
       Csc[gammaold], cold*Sqrt[1 - Cos[alphaold]^2 - Cos[betaold]^2 + 
         2*Cos[alphaold]*Cos[betaold]*Cos[gammaold] - Cos[gammaold]^2]*
       Csc[gammaold]}}}|>;

msgopold[number_]:=Module[{data},
	data=MSGOP[number];
	Print["Magnetic space group (BNS): ",data["MSG"]];
	Print["Lattice: ",data["BRAV"]];
	Print["Primitive Lattice Vactor: ", braLattold[data["BRAV"]][[2]]];
	Print["Conventional Lattice Vactor: ", braLattold[data["BRAV"]][[1]]];
	data["SymmetryOperation"]
];
mlgopold[mlg_] := Module[{},
  Print["Magnetic layer group:", {StringRiffle[mlg, "."], 
    lgOgToSymbol[mlg]}];
  Print["Lattice:", lgToBrav@First[mlg]];
  Print["Primitive Lattice Vactor:", 
   BasicVectorsMLG[lgToBrav@First[mlg]]];
  layerop[mlg]
  ];
mrgopold[mlg_] := Module[{},
  Print["Magnetic rod group:", {StringRiffle[mlg, "."], 
    rgOgToSymbol[mlg]}];
  Print["Lattice:", rgToBrav@First[mlg]];
  Print["Primitive Lattice Vactor:", 
   BasicVectorsMRG[rgToBrav@First[mlg]]];
  rodop[mlg]
  ];

bnsdictreverse=Association[Reverse/@Normal[bnsdictold]];
showMSGWyckoffold[msg_] := Module[{nu,nw, wyck},
If[NumberQ[msg],nu=bnsdictreverse[msg],nu=msg];
  Print["MSG:",nu];
  wyck = wyckoffmsg[nu];
  wyck = Reverse[wyck];
  nw = Length[wyck];
  wyck = 
   MapAt[If[NumberQ[#], Mod[#, 1], #] &, wyck, {;; , 4, ;; , 1, ;;}];
  wyck = 
   MapAt[Grid[#(*,Background\[Rule]{{LightYellow,LightBlue},None}*)] &,
    wyck, {;; , 4}];
  wyck = Delete[#, 2] & /@ wyck;
  PrependTo[
   wyck, {"Multiplicity", "Wyckoff Letter", 
    "Atomic Positions & Magnetization directions"}];
  Grid[wyck, Frame -> All, 
   Background -> {{LightYellow, LightBlue}, None}]
  ];





Options[symhamIIold] = {"wcc"->None};
symhamIIold[ham_, OptionsPattern[]]:=Module[{U,hII,wc},
If[OptionValue["wcc"]===None,wc=wccold,wc=OptionValue["wcc"]];
U=DiagonalMatrix[Table[Exp[-I{kxold,kyold,kzold} . tau],{tau,wc}]];
hII=ComplexExpand[ConjugateTranspose[U]] . ham . U;
Expand[TrigToExp@FullSimplify@hII]
];

symmetryopsIIold:=Module[{tmp},

Table[
tmp=symminfoold[[i,2]];
symmetryopsold[[i]]Table[Exp[I {kxold,kyold,kzold} . (
Inverse[tmp] . wccold[[k]]-wccold[[l]])],{k,Length[wccold]},{l,Length[wccold]}]
,{i,Length[symminfoold]}]
];

compactFormold:={#[[1]],#[[2]] . {"x","y","z"},#[[3]],#[[4]],MatrixForm[#[[5]]],
MatrixForm[#[[6]]]}&/@((symminfoold\[Transpose]~Join~{symmetryopsold}~Join~{symmetryopsIIold})\[Transpose]);

showbondsold[n_]:=Module[{data,bondlength,grid,hasbond,bpos,head},
data=bondclassifyold[[n]];

bondlength=data[[1,1]];
head=PadRight[{ToString[n-1]<>"-th neighbour, "<>"Bond length = "<> ToString@bondlength},3,SpanFromLeft];
grid={head,{"Atom position","Num of bonds",DisplayForm@SubsuperscriptBox["p","k","'"]}};
hasbond=bondclassifyold[[n]][[;;,3]][[;;,1,1]];
Do[

If[MemberQ[hasbond,pos],
bpos=First@FirstPosition[hasbond,pos];
grid=Append[grid,{pos,data[[bpos,2]],data[[bpos,3]][[;;,2]]}];,
grid=Append[grid,{pos,0,None}];
];
,{pos,Flatten[atomposold,1][[;;,1]]}];
Grid[grid,Frame->All]
]






GenerateGroupold[gens_,identityElement_,multiply_]:=Module[{i,j,ng,MAXORDER=200,orders,subs,mlist,g1,g2,g3},ng=Length[gens];orders=subs=Range[ng]*0;
For[i=1,i<=ng,i++,mlist=FoldList[multiply,Table[gens[[i]],MAXORDER]];
orders[[i]]=FirstPosition[mlist,identityElement][[1]];
subs[[i]]=mlist[[;;orders[[i]]]];];
g1=Union@@subs;
g2=Union@@Table[multiply[g1[[i]],g1[[j]]],{i,Length[g1]},{j,Length[g1]}];
g3=Complement[g2,g1];g1=g2;
While[g3!={},g2=Union@@Table[multiply[g3[[i]],g1[[j]]],{i,Length[g3]},{j,Length[g1]}];
g3=Complement[g2,g1];g1=g2;];
g2=SortBy[g2,#!=identityElement&]];
getGeneratorold[groupele_,identityElement_,multiply_]:=Module[
{try,group,tmptry,tmpgroup,seteq,groupeleDisorder,gereratorlist},
(*Greedy Alg, may not give the minimal Generator set*)
seteq[s1_,s2_]:=SubsetQ[s1,s2]&&SubsetQ[s2,s1];
If[groupele=={identityElement},Return[groupele]];
gereratorlist=Table[try={};
(*groupeleDisorder=RandomSample[groupele];*)
groupeleDisorder=groupele;
Do[If[ele==identityElement,Continue[]];
group=GenerateGroupold[try,identityElement,multiply];
tmptry=Append[try,ele];
(*Print[group];*)
tmpgroup=GenerateGroupold[tmptry,identityElement,multiply];
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
gens=getGeneratorold[opm,identityElement,times];
Flatten[FirstPosition[opm,#]&/@gens]
]


End[]
EndPackage[]

