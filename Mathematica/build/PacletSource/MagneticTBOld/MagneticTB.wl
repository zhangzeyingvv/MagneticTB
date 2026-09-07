(* ::Package:: *)

BeginPackage["MagneticTBOld`"]

(*
Author: Zhang Zeying
Email: zhangzeyingvv@gmail.com
*)


Begin["`Private`"]

installdir=DirectoryName[$InputFileName];
MSGDATA=Import[FileNameJoin[{installdir,"Data","MSGData.mx"}]];
wyckoffmsg=Import[FileNameJoin[{installdir,"Data","wyckoffMSG.mx"}]];
{rgOgToSymbol,rgToBrav,BasicVectorsMRG,rodop,grayrodold}=
  Import[FileNameJoin[{installdir,"Data","rod.mx"}]];
{lgOgToSymbol,lgToBrav,BasicVectorsMLG,layerop,graylayerold}=
  Import[FileNameJoin[{installdir,"Data","layer.mx"}]];

MSGOP=MSGDATA["MSGOP"];
grayold=MSGDATA["gray"];
typeIold=MSGDATA["typeI"];
typeIIIold=MSGDATA["typeIII"];
typeIVold=MSGDATA["typeIV"];
bnsdictold=MSGDATA["bnsdict"];
ogdictold=MSGDATA["ogdict"];
ognumdictold=MSGDATA["ognumdict"];




(*Print[$InputFileName,vv];*)
Options[initold]={
  latticeold->{{aold,0,0},{-(aold/2),(Sqrt[3] aold)/2,0},{0,0,cold}},
  lattparold->{aold->1,cold->3},
  wyckoffpositionold->{{{2/3,1/3,0},{0,0,1/2}}},
  symminformationold->{{"1",{{1,0,0},{0,1,0},{0,0,1}},{0,0,0},"F"}},
  origin={IdentityMatrix[3],{0,0,0}};
  basisFunctionsold->{"s"},
  debugQold->False
  };

initold[OptionsPattern[]]:=Module[{norm,repall,symopinit},
  lattold=(OptionValue[latticeold]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  norm=Cross[lattold[[1]],lattold[[2]]] . lattold[[3]];
  reclattold=2 Pi/norm {
    Cross[lattold[[2]],lattold[[3]]],
    Cross[lattold[[3]],lattold[[1]]],
    Cross[lattold[[1]],lattold[[2]]]
    };
  latpar=(OptionValue[lattparold]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  lattplot=lattold/.latpar;
  wyckoff=(OptionValue[wyckoffpositionold]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  symminfoold=(OptionValue[symminformationold]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  opsold=symminfoold[[;;,1]];
  Do[matrixRep[opsold[[i]]]=symminfoold[[i,2]],{i,Length[opsold]}];
  Do[tran[opsold[[i]]]=symminfoold[[i,3]],{i,Length[opsold]}];
  Do[timere[opsold[[i]]]=symminfoold[[i,4]],{i,Length[opsold]}];
  atomposold=DeleteDuplicates/@Table[
    {Mod[#[[1]]+#[[2]],1],#[[3]]}&/@({#[[2]] . wyck[[1]],#[[3]],
    Which[#[[4]]=="F",Det[#[[2]]]#[[2]] . wyck[[2]],
    #[[4]]=="T",-Det[#[[2]]]#[[2]] . wyck[[2]]]}&/@symminfoold),{wyck,wyckoff}];
  basisold=OptionValue[basisFunctionsold];
(*  If[Length[Intersection[VectorAngle@@@(Subsets[latt,{2}]/.latpar),{(2 \[Pi])/3,\[Pi]/3}]]>0,
  basisdict["dx2-y2"]=(x^2 - y^2);
  basisdict["dx2-y2dn"]->Reverse@{(x^2 - y^2),0};
  basisdict["dx2-y2up"]->{(x^2 - y^2),0};
];*)
  pointops=pointMatrixold[symminfoold,#,lattplot]&/@basisold;
(*Print[pointops];*)
  symopinit=Table[symop[i],{i,Length[wyckoff]}];
(*  Print[(symopinit)];*)
  symmetryopsold=Table[
    ArrayFlatten[Table[times[KroneckerDelta[i,j],symopinit[[i]][[k]]],
      {i,Length[wyckoff]},
        {j,Length[wyckoff]}]],
          {k,Length[opsold]}
      ];
(*Print[symmetryops];*)
(*Print[atompos];*)
  symmcompileold=Table[{i,symminfoold[[i]],symmetryopsold[[i]],
    (Inverse[Transpose[reclattold]] . ((Transpose[lattold] . symminfoold[[i,2]]) . Inverse@Transpose[lattold]) . Transpose[reclattold])/.latpar
      },{i,Length[symminfoold]}];
  bondclassifyold=Split[SortBy[Flatten[findn[atomposold[[;;,;;,1]],1],4],#[[1]]&],#1[[1]]==#2[[1]]&];
  bondclassifyold=Table[{#[[1,1]],Total@#[[;;,2]],Flatten[#[[;;,3]],1]}&/@Values[GroupBy[neigh,#[[3,1,1]]&]],{neigh,bondclassifyold}];
  wccold=Module[{natomperwyck,nbasesperwyck},
  natomperwyck=Length/@atomposold;
  nbasesperwyck=Length/@basisold;
  Flatten[Table[Table[Table[atomposold[[i,j,1]],nbasesperwyck[[i]]],{j,natomperwyck[[i]]}],{i,Length[natomperwyck]}],2]
  ];
  genindex=findgenind[symminfoold];
  Print["Generators:" ,symminfoold[[;;,{1,4}]][[genindex]]];

];


Options[symmetrizationHRInitold]={
  "Lattice"->{{aold,0,0},{-(aold/2),(Sqrt[3] aold)/2,0},{0,0,cold}},
  "LattPar"->{aold->1,cold->3},
  "WyckoffPosition"->{{{2/3,1/3,0},{0,0,1/2}}},
  "Symminformation"->{{"1",{{1,0,0},{0,1,0},{0,0,1}},{0,0,0},"F"}},
  "BasisFunctions"->{"s"},
  "Software"->"VASP"
  };

symmetrizationHRInitold[OptionsPattern[]]:=Module[{norm,repall,symopinit},
  lattold=(OptionValue["Lattice"]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  latpar=(OptionValue["LattPar"]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  lattplot=lattold/.latpar;
  wyckoff=(OptionValue["WyckoffPosition"]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  symminfoold=(OptionValue["Symminformation"]/.{xold_?InexactNumberQ:>Rationalize@Round[xold,.001]});
  software=OptionValue["Software"];
  opsold=symminfoold[[;;,1]];
  
  Do[matrixRep[opsold[[i]]]=symminfoold[[i,2]],{i,Length[opsold]}];
  Do[tran[opsold[[i]]]=symminfoold[[i,3]],{i,Length[opsold]}];
  Do[timere[opsold[[i]]]=symminfoold[[i,4]],{i,Length[opsold]}];
  atomposold=DeleteDuplicates/@Table[
    {Mod[#[[1]]+#[[2]],1],#[[3]]}&/@({#[[2]] . wyck[[1]],#[[3]],
    Which[#[[4]]=="F",Det[#[[2]]]#[[2]] . wyck[[2]],
    #[[4]]=="T",-Det[#[[2]]]#[[2]] . wyck[[2]]]}&/@symminfoold),{wyck,wyckoff}];
  basisold=OptionValue["BasisFunctions"];
  spinQ=Which[MissingQ[basisdictold[#]],ListQ[#],
              True,ListQ[basisdictold[#]]]&@basisold[[1,1]];
 (* Print[basisdict[basis[[1,1]]],spinQ];*)
(*  If[Length[Intersection[VectorAngle@@@(Subsets[latt,{2}]/.latpar),{(2 \[Pi])/3,\[Pi]/3}]]>0,
  basisdict["dx2-y2"]=(x^2 - y^2);
  basisdict["dx2-y2dn"]->Reverse@{(x^2 - y^2),0};
  basisdict["dx2-y2up"]->{(x^2 - y^2),0};
];*)
  pointops=pointMatrixold[symminfoold,#,lattplot]&/@basisold;
  (*Print[pointops];*)
  symopinit=Table[symopreal[i],{i,Length[wyckoff]}];
  (*Print[(symopinit)];*)
  symmetryopsold=Table[
    ArrayFlatten[Table[times[KroneckerDelta[i,j],symopinit[[i]][[k]]],
      {i,Length[wyckoff]},
        {j,Length[wyckoff]}]],
          {k,Length[opsold]}
      ];
      
  wccold=Module[{natomperwyck,nbasesperwyck},
  natomperwyck=Length/@atomposold;
  nbasesperwyck=Length/@basisold;
  Flatten[Table[Table[Table[atomposold[[i,j,1]],nbasesperwyck[[i]]],{j,natomperwyck[[i]]}],{i,Length[natomperwyck]}],2]
  ];
  Which[software=="VASP"&&spinQ,
  reordre=Join[Table[2i-1,{i,Length[wccold]/2}],Table[2i,{i,Length[wccold]/2}]];
  wcchr=wccold[[reordre]];
  symmetryopshr=Map[N@Table[#[[i,j]],{i,reordre},{j,reordre}]&,symmetryopsold,1],
  software=="VASP"&&Not[spinQ],
  wcchr=wccold;
  symmetryopshr=symmetryopsold
  ];
  (* Print[software=="VASP"&&Not[spinQ],spinQ,wcc,symmetryops];*)
  Association[{"wcc"->wcchr,"DR"->symmetryopshr,"symmetry"->symminfoold}]

];


basisdictold=<|
  "s"->1,
  "px"->xold,
  "py"->yold,
  "pz"->zold,
  "px+ipy"->xold+I yold,
  "px-ipy"->xold-I yold,
  "dx2-y2"->(*Sqrt[3]*)(xold^2 - yold^2),
  "dz2"->2zold^2-xold^2-yold^2,
  "dxy" -> 2 xold yold,
  "dyz"-> 2 yold zold,
  "dxz" ->2 xold zold,
  
  "sup"->{1,0},
  "pxup"->{xold,0},
  "pyup"->{yold,0},
  "pzup"->{zold,0},
  "px+ipy up"->{xold+I yold,0},
  "px-ipy up"->{xold-I yold,0},
  "dx2-y2up"->{(*Sqrt[3]*)(xold^2 - yold^2),0},
  "dz2up"->{2zold^2-xold^2-yold^2,0},
  "dxyup" ->{ 2 xold yold,0},
  "dyzup"->{2 yold zold,0},
  "dxzup" ->{2 xold zold,0},
  "sdn"->Reverse@{1,0},
  "pxdn"->Reverse@{xold,0},
  "pydn"->Reverse@{yold,0},
  "pzdn"->Reverse@{zold,0},
  "px+ipy dn"->Reverse@{xold+I yold,0},
  "px-ipy dn"->Reverse@{xold-I yold,0},
  "dx2-y2dn"->Reverse@{(*Sqrt[3]*)(xold^2 - yold^2),0},
  "dz2dn"->Reverse@{2zold^2-xold^2-yold^2,0},
  "dxydn" ->Reverse@{ 2 xold yold,0},
  "dyzdn"-> Reverse@{2 yold zold,0},
  "dxzdn" ->Reverse@{2 xold zold,0},
  "ptest3"->Reverse[1/Sqrt[2] {xold+ I yold,0}],
  "ptest4"-> Reverse[1/Sqrt[2] {0,xold- I yold}]

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
    Module[{alphaold, betaold, gammaold, nop,m,sx,order},
    If[Det[op] == -1, nop = -op, nop = op];
      {alphaold, betaold, gammaold} = -EulerAngles[nop];
      m=Transpose@{{E^(I gammaold/2) Cos[betaold/2] E^(I alphaold /2), E^(I gammaold/2) Sin[betaold/2] E^(-I alphaold /2)},
      {-E^(-I gammaold/2) Sin[betaold/2] E^(I alphaold /2), E^(-I gammaold/2) Cos[betaold/2] E^(-I alphaold /2)}};
      If[order==1,Return[IdentityMatrix[2]]];
      If[FullSimplify[MatrixPower[m,order]]==-IdentityMatrix[2],m,-m]
      ];




pointMatrixold[opsold_,orbs_,lattold_]:=
Block[{(*x,y,z,*)xi,yi,zi,$Assumptions,dict,opsrule,
nbases,bases,pbase,pbases,coeff,opmatrix,opmatrixs,tranU,spinrule,solve
},
  $Assumptions={{xold,yold,zold}\[Element]Reals};
  bases=basisdictold[#]&/@orbs;
  bases=If[Not@MissingQ[#],#,#[[2]]]&/@bases;
  nbases=Length[orbs];
(*Print[bases];*)
  Which[
    ListQ[bases[[1]]]==False,
(*    Print["Single Value Rep."];*)
    opsrule=FullSimplify@Table[{MapThread[Rule,{{xold,yold,zold},(Inverse[symm[[2]]] . ({xold, yold, zold} . Inverse[lattold]) . lattold)}],symm[[4]]},{symm,opsold}];
    (*Print[opsrule];*)
    opmatrixs=Table[
      pbases=If[rule[[2]]=="T",
      FullSimplify@Conjugate[(bases)/.First[rule]],
      (bases)/.First[rule]];
    opmatrix=Table[coeff[i,j],{i,nbases},{j,nbases}];
    solve=First@SolveAlways[pbases==opmatrix . bases,{xold,yold,zold}];
    (*Print[solve,nbases];*)
    If[Length[solve]!=nbases nbases,Print["Error! Your basis(es) cannot be the rep of point group, I refuse to do this job."];Abort[]];
    opmatrix=opmatrix/.solve;
    opmatrix=Transpose[opmatrix]
    ,{rule,opsrule}]
  ,Length[bases[[1]]]==2,
  Print["Double Value Rep."];
  opsrule=FullSimplify@Table[{MapThread[Rule,{{xold,yold,zold},(Inverse[symm[[2]]] . ({xold, yold, zold} . Inverse[lattold]) . lattold)}],symm[[4]]},{symm,opsold}];
  (*Table[Print[{Transpose[latt] . symm[[2]] . Inverse@Transpose@latt,spinMatrix2[Transpose[latt] . symm[[2]] . Inverse@Transpose@latt]}],{symm,ops}];*)
  spinrule=Table[ExpToTrig@spinMatrix2[Transpose[lattold] . symm[[2]] . Inverse@Transpose@lattold],{symm,opsold}];
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
  solve=Flatten[Table[SolveAlways[(pbases[[j]])==Total@Table[opmatrix[[j,k]]bases[[k]],{k,nbases}],{xold,yold,zold}],{j,nbases}]];
  If[Length[solve]!=nbases nbases,Print["Error! Your basis(es) cannot be the rep of point group, I refuse to do this job."];Abort[]];
  opmatrix=FullSimplify[opmatrix/.solve];
  opmatrix=Transpose[opmatrix]
,{i,Length[opsold]}]
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
wyckn=Length[atomposold[[n]]];
wyck=atomposold[[n]];
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
        {j,wyckn}]],{symm,symminfoold}]
      ];
      
      
symopreal[n_]:=(*symop[n]=*)Module[{wyckn,wyck,aftertr,irp,pos,qop},
wyckn=Length[atomposold[[n]]];
wyck=atomposold[[n]];
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
        {j,wyckn}]],{symm,symminfoold}]
      ];
times[xold_,yold_]:=If[xold==0,0,xold yold];





(*Find neighberhood*)
find[bond_, bondends_] := 
  Module[{start, end, l, lbond, bondtran, findeqv, postive, negtive},
   l = Length[opsold];
   lbond = Length[bondends];
   start =#[[1]]+#[[2]]&/@({#[[2]] . bond[[1]],#[[3]]}&/@symminfoold);
   end = #[[1]]+#[[2]]&/@({#[[2]] . bond[[2]],#[[3]]}&/@symminfoold);
   bondtran = Transpose[{start, end}];
(*Print[bondtran];*)
   (*findeqv[x_]:=Do[If[And@@IntegerQ/@(Flatten[x-bondends[[i]]]),
   debugPrint[i]],{i,lbond}];*)
   findeqv[xold_] := Do[If[
      (xold[[2]] - 
          xold[[1]] - (bondends[[i]][[2]] - bondends[[i]][[1]]) == {0, 0,
           0}) &&
       (Mod[xold[[2]] - bondends[[i]][[2]], 1] == {0, 
          0, 0}), Return[i]], {i, lbond}];
   postive = Table[findeqv[bondtran[[i]]], {i, l}] /. Null -> \[Pi] ;
   bondtran = Transpose[{end, start}];
   negtive = 
    Table[-findeqv[bondtran[[i]]], {i, l}] /. Null -> -\[Pi];
   Transpose[{postive, negtive}] /. \[Pi] -> Sequence[](*;
   Transpose[{postive}] /. \[Pi] -> Sequence[]*)
   ];
tobond=#[[2]]-#[[1]]&;

findn[cold_, lbond_] := 
 Module[{atom, natom, nwyck, mksc, sizesc, sc, split},
  sizesc = {3, 3, 3};
  nwyck = Length[cold];
  natom = Length /@ cold;
  mksc = Flatten[
     Table[# + {i, j, k}, {i, -sizesc[[1]], 
       sizesc[[1]]}, {j, -sizesc[[2]], sizesc[[2]]}, {k, -sizesc[[3]],
        sizesc[[3]]}], 2] &;
  sc = Map[mksc, cold, {2}];
  split[list_] := {N@#[[1, 1]], Length[#[[;; , 2]]], #[[;; , 2]]} & /@
     Split[list, First[#1] == First[#2] &];
  Table[split@
    Sort[Table[{EuclideanDistance[cold[[wkp, wkpeqvi]] . lattplot, 
        sc[[scwkp, scwkpeqvi, i]] . lattplot], {cold[[wkp, wkpeqvi]], 
        sc[[scwkp, scwkpeqvi, i]]}},
      {i, Length[sc[[scwkp, scwkpeqvi]]]}], N@#1[[1]] < N@#2[[1]] &],
   {wkp, nwyck},
   {wkpeqvi, Length[cold[[wkp]]]},
   {scwkp, nwyck},
   {scwkpeqvi, Length[cold[[scwkp]]]}]
  ];


(*Get Ham*)
unsymhamold[n_]:=Module[{bondends,bandind,hopold,bonds,toband,nbh,hami},
  bondends=Flatten[bondclassifyold[[n]][[;;,3]],1];
  bandind=Map[FirstPosition[Map[Flatten[Table[#+{i,j,k},{i,-3,3},{j,-3,3},{k,-3,3}],2]&,atomposold[[;;,;;,1]],{2}],#]&,bondends,{2}];
(*Print[bandind,Length[bandind],pointops];*)
  Do[hopold[bandind[[i]]]=Array[
    ToExpression["tr"<>ToString[i]<>ToString[#1]<>ToString[#2]<>"old+I ti"<>ToString[i]<>ToString[#1]<>ToString[#2]<>"old"]&,
    {Length[pointops[[bandind[[i,1,1]]]][[1]]],Length[pointops[[bandind[[i,2,1]]]][[1]]]}];
(*Print[i,hop[bandind[[i]]]];*)
  ,
  {i,Length[bandind]}];

  bonds=tobond/@bondends;
(*Print[bonds];*)
  toband=Total[((Length/@atomposold) (*(Length/@repMatrix[[;;,1]])*))[[;;#[[1]]-1]]]+#[[2]]&;
  nbh=(Length/@atomposold);
(*Print[nbh];*)
  Do[
    h[toband[{i,j}],toband[{k,l}]]=
      Table[0,{Length[pointops[[i]][[1]]]},{Length[pointops[[k]][[1]]]}],
    {i,Length[nbh]},{j,nbh[[i]]},{k,Length[nbh]},{l,nbh[[k]]}];

  Do[
    h[toband[bandind[[i,1]]],toband[bandind[[i,2]]]]+=
     Exp[I ({kxold,kyold,kzold}) . (bonds[[i]])] hopold[bandind[[i]]];
   ,{i,Length[bondends]}];
  hami=ArrayFlatten[Table[h[i,j],{i,Total@nbh},{j,Total@nbh}]]
];


Options[symhamold]={symmetrysetold->genindex,
"CartesianCoordinates" -> False
};
(*SetOptions[symham,symmetryset\[Rule]All];*)
symhamold[n_,OptionsPattern[]]:=Module[{conjh,h,phpmh,recR,para,sset,opset,cart},
  sset=OptionValue[symmetrysetold];
  cart=OptionValue["CartesianCoordinates"];
  If[sset===All,opset=Range[Length[symminfoold]],
  opset=sset,opset=Range[Length[symminfoold]]
(*Print["Use those symmetry constrains: ",TableForm[symminfo[[opset]],TableDepth\[Rule]2]]*)
];
(*Print[opset,Range[Length[symminfo]]];*)
  conjh=(ComplexExpand@ConjugateTranspose@unsymhamold[n])-unsymhamold[n];

  If[n==1,h=unsymhamold[n]/.ToRules@Reduce[conjh==0],

    h=unsymhamold[n]/.ToRules@Reduce[
      Flatten[Map[Total/@(Values[GroupBy[Cases[#,xold_ E^yold_:>{xold,Expand[yold]}],Last]][[;;,;;,1]])&,TrigToExp[conjh],{2}]]==0];
      
    h=h/.ToRules@Reduce[
      Flatten[Map[Total/@(Values[GroupBy[Cases[{#},xold_ E^yold_:>{xold,Expand[yold]}],Last]][[;;,;;,1]])&,TrigToExp[conjh],{2}]]==0]
  ];
(*Print[MatrixForm@ExpToTrig@h];*)

  phpmh=Table[
  (*recR= Inverse[Transpose[reclatt]] . ((Transpose[latt] . symminfo[[i,2]]) . Inverse@Transpose[latt]) . Transpose[reclatt];*)
  recR=symmcompileold[[i,4]];
    recR=Inverse[recR];
  Which[

    symminfoold[[i,4]]=="F",

    (Inverse[symmetryopsold[[i]]] . h . symmetryopsold[[i]])-(h/.Thread[{kxold,kyold,kzold}->recR . {kxold,kyold,kzold}])


    ,
    symminfoold[[i,4]]=="T",

    (Expand@TrigToExp@(Inverse[symmetryopsold[[i]]] . (h) . symmetryopsold[[i]]))-Expand@TrigToExp@(ComplexExpand@Conjugate[h/.Thread[
    {kxold,kyold,kzold}->recR . {-kxold,-kyold,-kzold}]])
  ]
,{i,opset}];
phpmh=Expand[phpmh];
  If[n==1,
(*Print[First@Solve[phpmh\[Equal]0]];*)
    h=h/.ToRules@Reduce[phpmh==0],
    phpmh=Join[Map[Total@Cases[#,xx_ E^yy_:>xx E^Simplify@yy]&,phpmh,{3}],
    Map[Total@Cases[{#},xx_ E^yy_:>xx E^Simplify@yy]&,phpmh,{3}]
  ];
h=h/.FullSimplify@ToRules@Reduce[DeleteDuplicates@Flatten[Map[Total/@(Values[GroupBy[Cases[#,xold_ E^yold_:>{xold,Expand[yold]}],Last]][[;;,;;,1]])&,phpmh,{3}]]==0];
(*Print[h];*)
h=h/.FullSimplify@ToRules@Reduce[DeleteDuplicates@Flatten[Map[Total/@(Values[GroupBy[Cases[{#},xold_ E^yold_:>{xold,Expand[yold]}],Last]][[;;,;;,1]])&,phpmh,{3}]]==0];
(*Print[h];*)
];

(*Print["params:"Variables[h]];*)
para=Thread[Which[
n==1,#->Table[ToExpression["e"<>ToString[i]<>"old"],{i,Length[#]}]&@Variables[h],
n==2,#->Table[ToExpression["t"<>ToString[i]<>"old"],{i,Length[#]}]&@Variables[h],
n==3,#->Table[ToExpression["r"<>ToString[i]<>"old"],{i,Length[#]}]&@Variables[h],
n==4,#->Table[ToExpression["s"<>ToString[i]<>"old"],{i,Length[#]}]&@Variables[h],
True,#->Table[ToExpression["p"<>ToString[n]<>"n"<>ToString[i]<>"old"],{i,Length[#]}]&@Variables[h]]
];

h=h/.para;
If[cart, h=h/.Thread[{kxold,kyold,kzold}->({kxold,kyold,kzold}(2Pi)) . Inverse@reclattold]];
(*Print[para];*)
Print["params:",Variables[h]];
h

];





End[]


EndPackage[]
