(* ::Package:: *)

With[{p=DirectoryName[$InputFileName]}, If[!MemberQ[$Path,p],AppendTo[$Path, p]]];



BeginPackage["MagneticTB`"]

(*
Author: Zhang Zeying
Email: zhangzeyingvv@gmail.com
*)


Begin["`Private`"]

installdir=DirectoryName[$InputFileName];
MSGDATA=Import[installdir<>"MSGData.mx"];
MSGOP=MSGDATA["MSGOP"];
gray=MSGDATA["gray"];
typeI=MSGDATA["typeI"];
typeIII=MSGDATA["typeIII"];
typeIV=MSGDATA["typeIV"];
bnsdict=MSGDATA["bnsdict"];
ogdict=MSGDATA["ogdict"];
ognumdict=MSGDATA["ognumdict"];

(*Print[$InputFileName,vv];*)
Options[init]={
  lattice->{{a,0,0},{-(a/2),(Sqrt[3] a)/2,0},{0,0,c}},
  lattpar->{a->1,c->3},
  wyckoffposition->{{{2/3,1/3,0},{0,0,1/2}}},
  symminformation->{{"1",{{1,0,0},{0,1,0},{0,0,1}},{0,0,0},"F"}},
  origin={IdentityMatrix[3],{0,0,0}};
  basisFunctions->{"s"},
  debugQ->False
  };

init[OptionsPattern[]]:=Module[{norm,repall,symopinit},
  latt=(OptionValue[lattice]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  norm=Cross[latt[[1]],latt[[2]]] . latt[[3]];
  reclatt=2 Pi/norm {
    Cross[latt[[2]],latt[[3]]],
    Cross[latt[[3]],latt[[1]]],
    Cross[latt[[1]],latt[[2]]]
    };
  latpar=(OptionValue[lattpar]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  lattplot=latt/.latpar;
  wyckoff=(OptionValue[wyckoffposition]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  symminfo=(OptionValue[symminformation]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  ops=symminfo[[;;,1]];
  Do[matrixRep[ops[[i]]]=symminfo[[i,2]],{i,Length[ops]}];
  Do[tran[ops[[i]]]=symminfo[[i,3]],{i,Length[ops]}];
  Do[timere[ops[[i]]]=symminfo[[i,4]],{i,Length[ops]}];
  atompos=DeleteDuplicates/@Table[
    {Mod[#[[1]]+#[[2]],1],#[[3]]}&/@({#[[2]] . wyck[[1]],#[[3]],
    Which[#[[4]]=="F",Det[#[[2]]]#[[2]] . wyck[[2]],
    #[[4]]=="T",-Det[#[[2]]]#[[2]] . wyck[[2]]]}&/@symminfo),{wyck,wyckoff}];
  basis=OptionValue[basisFunctions];
(*  If[Length[Intersection[VectorAngle@@@(Subsets[latt,{2}]/.latpar),{(2 \[Pi])/3,\[Pi]/3}]]>0,
  basisdict["dx2-y2"]=(x^2 - y^2);
  basisdict["dx2-y2dn"]->Reverse@{(x^2 - y^2),0};
  basisdict["dx2-y2up"]->{(x^2 - y^2),0};
];*)
  pointops=pointMatrix[symminfo,#,lattplot]&/@basis;
(*Print[pointops];*)
  symopinit=Table[symop[i],{i,Length[wyckoff]}];
(*  Print[(symopinit)];*)
  symmetryops=Table[
    ArrayFlatten[Table[times[KroneckerDelta[i,j],symopinit[[i]][[k]]],
      {i,Length[wyckoff]},
        {j,Length[wyckoff]}]],
          {k,Length[ops]}
      ];
(*Print[symmetryops];*)
(*Print[atompos];*)
  symmcompile=Table[{i,symminfo[[i]],symmetryops[[i]],
    (Inverse[Transpose[reclatt]] . ((Transpose[latt] . symminfo[[i,2]]) . Inverse@Transpose[latt]) . Transpose[reclatt])/.latpar
      },{i,Length[symminfo]}];
  bondclassify=Split[SortBy[Flatten[findn[atompos[[;;,;;,1]],1],4],#[[1]]&],#1[[1]]==#2[[1]]&];
  bondclassify=Table[{#[[1,1]],Total@#[[;;,2]],Flatten[#[[;;,3]],1]}&/@Values[GroupBy[neigh,#[[3,1,1]]&]],{neigh,bondclassify}];
  wcc=Module[{natomperwyck,nbasesperwyck},
  natomperwyck=Length/@atompos;
  nbasesperwyck=Length/@basis;
  Flatten[Table[Table[Table[atompos[[i,j,1]],nbasesperwyck[[i]]],{j,natomperwyck[[i]]}],{i,Length[natomperwyck]}],2]
  ];
  genindex=findgenind[symminfo];
  Print["Generators:" ,symminfo[[;;,{1,4}]][[genindex]]];

];


Options[symmetrizationHRInit]={
  "Lattice"->{{a,0,0},{-(a/2),(Sqrt[3] a)/2,0},{0,0,c}},
  "LattPar"->{a->1,c->3},
  "WyckoffPosition"->{{{2/3,1/3,0},{0,0,1/2}}},
  "Symminformation"->{{"1",{{1,0,0},{0,1,0},{0,0,1}},{0,0,0},"F"}},
  "BasisFunctions"->{"s"},
  "Software"->"VASP"
  };

symmetrizationHRInit[OptionsPattern[]]:=Module[{norm,repall,symopinit},
  latt=(OptionValue["Lattice"]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  latpar=(OptionValue["LattPar"]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  lattplot=latt/.latpar;
  wyckoff=(OptionValue["WyckoffPosition"]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  symminfo=(OptionValue["Symminformation"]/.{x_?InexactNumberQ:>Rationalize@Round[x,.001]});
  software=OptionValue["Software"];
  ops=symminfo[[;;,1]];
  
  Do[matrixRep[ops[[i]]]=symminfo[[i,2]],{i,Length[ops]}];
  Do[tran[ops[[i]]]=symminfo[[i,3]],{i,Length[ops]}];
  Do[timere[ops[[i]]]=symminfo[[i,4]],{i,Length[ops]}];
  atompos=DeleteDuplicates/@Table[
    {Mod[#[[1]]+#[[2]],1],#[[3]]}&/@({#[[2]] . wyck[[1]],#[[3]],
    Which[#[[4]]=="F",Det[#[[2]]]#[[2]] . wyck[[2]],
    #[[4]]=="T",-Det[#[[2]]]#[[2]] . wyck[[2]]]}&/@symminfo),{wyck,wyckoff}];
  basis=OptionValue["BasisFunctions"];
  spinQ=Which[MissingQ[basisdict[#]],ListQ[#],
              True,ListQ[basisdict[#]]]&@basis[[1,1]];
 (* Print[basisdict[basis[[1,1]]],spinQ];*)
(*  If[Length[Intersection[VectorAngle@@@(Subsets[latt,{2}]/.latpar),{(2 \[Pi])/3,\[Pi]/3}]]>0,
  basisdict["dx2-y2"]=(x^2 - y^2);
  basisdict["dx2-y2dn"]->Reverse@{(x^2 - y^2),0};
  basisdict["dx2-y2up"]->{(x^2 - y^2),0};
];*)
  pointops=pointMatrix[symminfo,#,lattplot]&/@basis;
  (*Print[pointops];*)
  symopinit=Table[symopreal[i],{i,Length[wyckoff]}];
  (*Print[(symopinit)];*)
  symmetryops=Table[
    ArrayFlatten[Table[times[KroneckerDelta[i,j],symopinit[[i]][[k]]],
      {i,Length[wyckoff]},
        {j,Length[wyckoff]}]],
          {k,Length[ops]}
      ];
      
  wcc=Module[{natomperwyck,nbasesperwyck},
  natomperwyck=Length/@atompos;
  nbasesperwyck=Length/@basis;
  Flatten[Table[Table[Table[atompos[[i,j,1]],nbasesperwyck[[i]]],{j,natomperwyck[[i]]}],{i,Length[natomperwyck]}],2]
  ];
  Which[software=="VASP"&&spinQ,
  reordre=Join[Table[2i-1,{i,Length[wcc]/2}],Table[2i,{i,Length[wcc]/2}]];
  wcchr=wcc[[reordre]];
  symmetryopshr=Map[N@Table[#[[i,j]],{i,reordre},{j,reordre}]&,symmetryops,1],
  software=="VASP"&&Not[spinQ],
  wcchr=wcc;
  symmetryopshr=symmetryops
  ];
  (* Print[software=="VASP"&&Not[spinQ],spinQ,wcc,symmetryops];*)
  Association[{"wcc"->wcchr,"DR"->symmetryopshr,"symmetry"->symminfo}]

];


basisdict=<|
  "s"->1,
  "px"->x,
  "py"->y,
  "pz"->z,
  "px+ipy"->x+I y,
  "px-ipy"->x-I y,
  "dx2-y2"->(*Sqrt[3]*)(x^2 - y^2),
  "dz2"->2z^2-x^2-y^2,
  "dxy" -> 2 x y,
  "dyz"-> 2 y z,
  "dxz" ->2 x z,
  
  "sup"->{1,0},
  "pxup"->{x,0},
  "pyup"->{y,0},
  "pzup"->{z,0},
  "px+ipy up"->{x+I y,0},
  "px-ipy up"->{x-I y,0},
  "dx2-y2up"->{(*Sqrt[3]*)(x^2 - y^2),0},
  "dz2up"->{2z^2-x^2-y^2,0},
  "dxyup" ->{ 2 x y,0},
  "dyzup"->{2 y z,0},
  "dxzup" ->{2 x z,0},
  "sdn"->Reverse@{1,0},
  "pxdn"->Reverse@{x,0},
  "pydn"->Reverse@{y,0},
  "pzdn"->Reverse@{z,0},
  "px+ipy dn"->Reverse@{x+I y,0},
  "px-ipy dn"->Reverse@{x-I y,0},
  "dx2-y2dn"->Reverse@{(*Sqrt[3]*)(x^2 - y^2),0},
  "dz2dn"->Reverse@{2z^2-x^2-y^2,0},
  "dxydn" ->Reverse@{ 2 x y,0},
  "dyzdn"-> Reverse@{2 y z,0},
  "dxzdn" ->Reverse@{2 x z,0},
  "ptest3"->Reverse[1/Sqrt[2] {x+ I y,0}],
  "ptest4"-> Reverse[1/Sqrt[2] {0,x- I y}]

(*  "s1/2+1/2"->{1,0},
  "s1/2-1/2"->{0,-1},

  "p1/2+1/2"-> {-z,x+ I y},
  "p1/2-1/2"-> {-x+ I y,z},

  "ptest1"-> {z,0},
  "ptest2"-> {0, z},

  "ptest3"->1/Sqrt[2] {x+ I y,0},
  "ptest4"-> 1/Sqrt[2] {0,x- I y},
  "ptest5"-> -1/Sqrt[2] {0,x+ I y},
  "ptest6"->-1/Sqrt[2] {x- I y,0},

  "ptest11"-> {(I (x^2-y^2))/Sqrt[2],(x^2-y^2)/Sqrt[2]},
  "ptest12"-> {0,x},
  "ptest13"-> {y,0},
  "ptest14"->{0,y},


"p3/2+3/2"->{x+ I y,0},
"p3/2+1/2"-> {2z,x+ I y},
"p3/2-1/2"-> {x- I y,2z},
"p3/2-3/2"->{0,x- I y},

"d3/2+3/2"-> {z(x+I y),2I x y +(x^2-y^2)},
"d3/2+1/2"-> {x^2+y^2+z^2-3z^2,-3z(x+I y)},
"d3/2-1/2"-> {-3z(x-I y),-(x^2+y^2+z^2-3z^2)},
"d3/2-3/2"-> {2I x y -(x^2-y^2),z(x-I y)},

"d5/2+5/2"-> {0,0},
"d5/2+3/2"-> {0,0},
"d5/2+1/2"-> {0,0},
"d5/2-1/2"-> {0,0},
"d5/2-3/2"-> {0,0},
"d5/2-5/2"-> {0,0}*)

|>;





spinMatrix[op_] := 
    Module[{\[Alpha], \[Beta], \[Gamma], nop,m,sx,order},
    If[Det[op] == -1, nop = -op, nop = op];
      {\[Alpha], \[Beta], \[Gamma]} = -EulerAngles[nop];
      m=Transpose@{{E^(I \[Gamma]/2) Cos[\[Beta]/2] E^(I \[Alpha] /2), E^(I \[Gamma]/2) Sin[\[Beta]/2] E^(-I \[Alpha] /2)},
      {-E^(-I \[Gamma]/2) Sin[\[Beta]/2] E^(I \[Alpha] /2), E^(-I \[Gamma]/2) Cos[\[Beta]/2] E^(-I \[Alpha] /2)}};
      If[order==1,Return[IdentityMatrix[2]]];
      If[FullSimplify[MatrixPower[m,order]]==-IdentityMatrix[2],m,-m]
      ];




pointMatrix[ops_,orbs_,latt_]:=
Block[{(*x,y,z,*)xi,yi,zi,$Assumptions,dict,opsrule,
nbases,bases,pbase,pbases,coeff,opmatrix,opmatrixs,tranU,spinrule,solve
},
  $Assumptions={{x,y,z}\[Element]Reals};
  bases=basisdict[#]&/@orbs;
  bases=If[Not@MissingQ[#],#,#[[2]]]&/@bases;
  nbases=Length[orbs];
(*Print[bases];*)
  Which[
    ListQ[bases[[1]]]==False,
(*    Print["Single Value Rep."];*)
    opsrule=FullSimplify@Table[{MapThread[Rule,{{x,y,z},(Inverse[symm[[2]]] . ({x, y, z} . Inverse[latt]) . latt)}],symm[[4]]},{symm,ops}];
    (*Print[opsrule];*)
    opmatrixs=Table[
      pbases=If[rule[[2]]=="T",
      FullSimplify@Conjugate[(bases)/.First[rule]],
      (bases)/.First[rule]];
    opmatrix=Table[coeff[i,j],{i,nbases},{j,nbases}];
    solve=First@SolveAlways[pbases==opmatrix . bases,{x,y,z}];
    (*Print[solve,nbases];*)
    If[Length[solve]!=nbases nbases,Print["Error! Your basis(es) cannot be the rep of point group, I refuse to do this job."];Abort[]];
    opmatrix=opmatrix/.solve;
    opmatrix=Transpose[opmatrix]
    ,{rule,opsrule}]
  ,Length[bases[[1]]]==2,
  Print["Double Value Rep."];
  opsrule=FullSimplify@Table[{MapThread[Rule,{{x,y,z},(Inverse[symm[[2]]] . ({x, y, z} . Inverse[latt]) . latt)}],symm[[4]]},{symm,ops}];
  (*Table[Print[{Transpose[latt] . symm[[2]] . Inverse@Transpose@latt,spinMatrix2[Transpose[latt] . symm[[2]] . Inverse@Transpose@latt]}],{symm,ops}];*)
  spinrule=Table[ExpToTrig@spinMatrix2[Transpose[latt] . symm[[2]] . Inverse@Transpose@latt],{symm,ops}];
  opmatrixs=Table[
    pbases=If[opsrule[[i]][[2]]=="T",(*Print[Simplify@Conjugate[(dict[#]&/@orbs)/.First[rule]]];*)
      pbases=ComplexExpand[I PauliMatrix[2] . Conjugate[#]]&/@bases;
      pbases=(pbases)/.First[opsrule[[i]]];
      pbases=spinrule[[i]] . #&/@pbases
     ,
      pbases=(bases)/.First[opsrule[[i]]];
      pbases=spinrule[[i]] . #&/@pbases
    ];
  pbases=Simplify@pbases;
(*Print[pbases];*)

  opmatrix=Table[coeff[i,j],{i,nbases},{j,nbases}];
(*Print[MatrixForm@FullSimplify@opmatrix/.Flatten[Table[SolveAlways[(pbases[[j]])\[Equal]Total@Table[opmatrix[[j,k]]bases[[k]],{k,nbases}],{x,y,z}],{j,nbases}]]];*)
  solve=Flatten[Table[SolveAlways[(pbases[[j]])==Total@Table[opmatrix[[j,k]]bases[[k]],{k,nbases}],{x,y,z}],{j,nbases}]];
  If[Length[solve]!=nbases nbases,Print["Error! Your basis(es) cannot be the rep of point group, I refuse to do this job."];Abort[]];
  opmatrix=FullSimplify[opmatrix/.solve];
  opmatrix=Transpose[opmatrix]
,{i,Length[ops]}]
];

If[And@@(FullSimplify[ConjugateTranspose[#] . #]==IdentityMatrix[nbases]&/@opmatrixs),opmatrixs,
   tranU=DiagonalMatrix[Sqrt[#[[1]]]] . #[[2]]&[Eigensystem[Total[ConjugateTranspose[#] . #&/@opmatrixs]]];
   tranU = SortBy[tranU, First@First@Position[#, _?(# != 0 &), 1] &];
   If[DiagonalMatrixQ[tranU],
   (*Print["Note: the basis functions cannot carry a unitary representation, I will multiply the basis functions by an array:", (#/(#[[1]]))&@Diagonal[tranU]];*)
   opmatrixs=tranU . # . Inverse[tranU]&/@opmatrixs,
   Print["Warning: your basis set cannot carry a unitary representation, I will do a linear combination of your basis functions:", tranU];
   Print["However, I strongly recommend re-selecting the basis functions which can carry a unitary representation, because the order of basis functions in Hamiltonian may change"];
   opmatrixs=tranU . # . Inverse[tranU]&/@opmatrixs]]

];





(*
Matrix Rep for Space Group Symmetry Opeartors
*)
symop[n_]:=(*symop[n]=*)Module[{wyckn,wyck,aftertr,irp,pos,qop},
wyckn=Length[atompos[[n]]];
wyck=atompos[[n]];
pos=0;
Table[pos+=1;
ArrayFlatten[Table[
(*  Print[symm];*)
  qop=symm;
  qop={#[[1]],Inverse[#[[2]]],-Inverse[#[[2]]] . #[[3]],#[[4]]}&@qop;
  If[
    aftertr=({Mod[#[[1]]+#[[2]],1],#[[3]]}&@({qop[[2]] . wyck[[i,1]],qop[[3]],
      Which[
        qop[[4]]=="F",Det[qop[[2]]]qop[[2]] . wyck[[i,2]],
        qop[[4]]=="T",-Det[qop[[2]]]qop[[2]] . wyck[[i,2]]
           ]}));
(*      Print[aftertr];*)
      aftertr==wyck[[j]],
      pointops[[n]][[pos]],0],
      {i,wyckn},
        {j,wyckn}]],{symm,symminfo}]
      ];
      
      
symopreal[n_]:=(*symop[n]=*)Module[{wyckn,wyck,aftertr,irp,pos,qop},
wyckn=Length[atompos[[n]]];
wyck=atompos[[n]];
pos=0;
Table[pos+=1;
ArrayFlatten[Table[
(*  Print[symm];*)
  qop=symm;
  qop={#[[1]],Inverse[#[[2]]],-Inverse[#[[2]]] . #[[3]],#[[4]]}&@qop;
  If[
    i==j,
      pointops[[n]][[pos]],0],
      {i,wyckn},
        {j,wyckn}]],{symm,symminfo}]
      ];
times[x_,y_]:=If[x==0,0,x y];





(*Find neighberhood*)
find[bond_, bondends_] := 
  Module[{start, end, l, lbond, bondtran, findeqv, postive, negtive},
   l = Length[ops];
   lbond = Length[bondends];
   start =#[[1]]+#[[2]]&/@({#[[2]] . bond[[1]],#[[3]]}&/@symminfo);
   end = #[[1]]+#[[2]]&/@({#[[2]] . bond[[2]],#[[3]]}&/@symminfo);
   bondtran = Transpose[{start, end}];
(*Print[bondtran];*)
   (*findeqv[x_]:=Do[If[And@@IntegerQ/@(Flatten[x-bondends[[i]]]),
   debugPrint[i]],{i,lbond}];*)
   findeqv[x_] := Do[If[
      (x[[2]] - 
          x[[1]] - (bondends[[i]][[2]] - bondends[[i]][[1]]) == {0, 0,
           0}) &&
       (Mod[x[[2]] - bondends[[i]][[2]], 1] == {0, 
          0, 0}), Return[i]], {i, lbond}];
   postive = Table[findeqv[bondtran[[i]]], {i, l}] /. Null -> \[Pi] ;
   bondtran = Transpose[{end, start}];
   negtive = 
    Table[-findeqv[bondtran[[i]]], {i, l}] /. Null -> -\[Pi];
   Transpose[{postive, negtive}] /. \[Pi] -> Sequence[](*;
   Transpose[{postive}] /. \[Pi] -> Sequence[]*)
   ];
tobond=#[[2]]-#[[1]]&;

findn[c_, lbond_] := 
 Module[{atom, natom, nwyck, mksc, sizesc, sc, split},
  sizesc = {3, 3, 3};
  nwyck = Length[c];
  natom = Length /@ c;
  mksc = Flatten[
     Table[# + {i, j, k}, {i, -sizesc[[1]], 
       sizesc[[1]]}, {j, -sizesc[[2]], sizesc[[2]]}, {k, -sizesc[[3]],
        sizesc[[3]]}], 2] &;
  sc = Map[mksc, c, {2}];
  split[list_] := {N@#[[1, 1]], Length[#[[;; , 2]]], #[[;; , 2]]} & /@
     Split[list, First[#1] == First[#2] &];
  Table[split@
    Sort[Table[{EuclideanDistance[c[[wkp, wkpeqvi]] . lattplot, 
        sc[[scwkp, scwkpeqvi, i]] . lattplot], {c[[wkp, wkpeqvi]], 
        sc[[scwkp, scwkpeqvi, i]]}},
      {i, Length[sc[[scwkp, scwkpeqvi]]]}], N@#1[[1]] < N@#2[[1]] &],
   {wkp, nwyck},
   {wkpeqvi, Length[c[[wkp]]]},
   {scwkp, nwyck},
   {scwkpeqvi, Length[c[[scwkp]]]}]
  ];


(*Get Ham*)
unsymham[n_]:=Module[{bondends,bandind,hop,bonds,toband,nbh,hami},
  bondends=Flatten[bondclassify[[n]][[;;,3]],1];
  bandind=Map[FirstPosition[Map[Flatten[Table[#+{i,j,k},{i,-3,3},{j,-3,3},{k,-3,3}],2]&,atompos[[;;,;;,1]],{2}],#]&,bondends,{2}];
(*Print[bandind,Length[bandind],pointops];*)
  Do[hop[bandind[[i]]]=Array[
    ToExpression["tr"<>ToString[i]<>ToString[#1]<>ToString[#2]<>"+I ti"<>ToString[i]<>ToString[#1]<>ToString[#2]]&,
    {Length[pointops[[bandind[[i,1,1]]]][[1]]],Length[pointops[[bandind[[i,2,1]]]][[1]]]}];
(*Print[i,hop[bandind[[i]]]];*)
  ,
  {i,Length[bandind]}];

  bonds=tobond/@bondends;
(*Print[bonds];*)
  toband=Total[((Length/@atompos) (*(Length/@repMatrix[[;;,1]])*))[[;;#[[1]]-1]]]+#[[2]]&;
  nbh=(Length/@atompos);
(*Print[nbh];*)
  Do[
    h[toband[{i,j}],toband[{k,l}]]=
      Table[0,{Length[pointops[[i]][[1]]]},{Length[pointops[[k]][[1]]]}],
    {i,Length[nbh]},{j,nbh[[i]]},{k,Length[nbh]},{l,nbh[[k]]}];

  Do[
    h[toband[bandind[[i,1]]],toband[bandind[[i,2]]]]+=
     Exp[I ({kx,ky,kz}) . (bonds[[i]])] hop[bandind[[i]]];
   ,{i,Length[bondends]}];
  hami=ArrayFlatten[Table[h[i,j],{i,Total@nbh},{j,Total@nbh}]]
];


Options[symham]={symmetryset->genindex,
"CartesianCoordinates" -> False
};
(*SetOptions[symham,symmetryset\[Rule]All];*)
symham[n_,OptionsPattern[]]:=Module[{conjh,h,phpmh,recR,para,sset,opset,cart},
  sset=OptionValue["symmetryset"];
  cart=OptionValue["CartesianCoordinates"];
  If[sset===All,opset=Range[Length[symminfo]],
  opset=sset,opset=Range[Length[symminfo]]
(*Print["Use those symmetry constrains: ",TableForm[symminfo[[opset]],TableDepth\[Rule]2]]*)
];
(*Print[opset,Range[Length[symminfo]]];*)
  conjh=(ComplexExpand@ConjugateTranspose@unsymham[n])-unsymham[n];

  If[n==1,h=unsymham[n]/.ToRules@Reduce[conjh==0],

    h=unsymham[n]/.ToRules@Reduce[
      Flatten[Map[Total/@(Values[GroupBy[Cases[#,x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,TrigToExp[conjh],{2}]]==0];
      
    h=h/.ToRules@Reduce[
      Flatten[Map[Total/@(Values[GroupBy[Cases[{#},x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,TrigToExp[conjh],{2}]]==0]
  ];
(*Print[MatrixForm@ExpToTrig@h];*)

  phpmh=Table[
  (*recR= Inverse[Transpose[reclatt]] . ((Transpose[latt] . symminfo[[i,2]]) . Inverse@Transpose[latt]) . Transpose[reclatt];*)
  recR=symmcompile[[i,4]];
    recR=Inverse[recR];
  Which[

    symminfo[[i,4]]=="F",

    (Inverse[symmetryops[[i]]] . h . symmetryops[[i]])-(h/.Thread[{kx,ky,kz}->recR . {kx,ky,kz}])


    ,
    symminfo[[i,4]]=="T",

    (Expand@TrigToExp@(Inverse[symmetryops[[i]]] . (h) . symmetryops[[i]]))-Expand@TrigToExp@(ComplexExpand@Conjugate[h/.Thread[
    {kx,ky,kz}->recR . {-kx,-ky,-kz}]])
  ]
,{i,opset}];
phpmh=Expand[phpmh];
  If[n==1,
(*Print[First@Solve[phpmh\[Equal]0]];*)
    h=h/.ToRules@Reduce[phpmh==0],
    phpmh=Join[Map[Total@Cases[#,xx_ E^yy_:>xx E^Simplify@yy]&,phpmh,{3}],
    Map[Total@Cases[{#},xx_ E^yy_:>xx E^Simplify@yy]&,phpmh,{3}]
  ];
h=h/.FullSimplify@ToRules@Reduce[DeleteDuplicates@Flatten[Map[Total/@(Values[GroupBy[Cases[#,x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,phpmh,{3}]]==0];
(*Print[h];*)
h=h/.FullSimplify@ToRules@Reduce[DeleteDuplicates@Flatten[Map[Total/@(Values[GroupBy[Cases[{#},x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,phpmh,{3}]]==0];
(*Print[h];*)
];

(*Print["params:"Variables[h]];*)
para=Thread[Which[
n==1,#->Table[ToExpression["e"<>ToString[i]],{i,Length[#]}]&@Variables[h],
n==2,#->Table[ToExpression["t"<>ToString[i]],{i,Length[#]}]&@Variables[h],
n==3,#->Table[ToExpression["r"<>ToString[i]],{i,Length[#]}]&@Variables[h],
n==4,#->Table[ToExpression["s"<>ToString[i]],{i,Length[#]}]&@Variables[h],
True,#->Table[ToExpression["p"<>ToString[n]<>"n"<>ToString[i]],{i,Length[#]}]&@Variables[h]]
];

h=h/.para;
If[cart, h=h/.Thread[{kx,ky,kz}->({kx,ky,kz}(2Pi)) . Inverse@reclatt]];
(*Print[para];*)
Print["params:",Variables[h]];
h

];





End[]


EndPackage[]

