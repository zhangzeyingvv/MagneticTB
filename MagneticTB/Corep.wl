(* ::Package:: *)

BeginPackage["MagneticTB`"]


Begin["`Private`"]


getSGLatt=SpaceGroupIrep`getSGLatt;
BasicVectors=SpaceGroupIrep`BasicVectors;
getRotMat=SpaceGroupIrep`getRotMat;
getSpinRotOp=SpaceGroupIrep`getSpinRotOp;
getRotName=SpaceGroupIrep`getRotName


MSGSymStd=MSGCorep`MSGSymStd;
getMSGElem=MSGCorep`getMSGElem;
getBandCorep=MSGCorep`getBandCorep;
showBandCorep=MSGCorep`showBandCorep;
spinMatrix[op_] := 
    N@Module[{\[Alpha], \[Beta], \[Gamma], nop,m,sx},
    sx=PauliMatrix[2];
    If[Det[op] == -1, nop = -op, nop = op];
      {\[Alpha], \[Beta], \[Gamma]} = -EulerAngles[nop];
      m={{E^(I \[Gamma]/2) Cos[\[Beta]/2] E^(I \[Alpha] /2), E^(I \[Gamma]/2) Sin[\[Beta]/2] E^(-I \[Alpha] /2)},
      {-E^(-I \[Gamma]/2) Sin[\[Beta]/2] E^(I \[Alpha] /2), E^(-I \[Gamma]/2) Cos[\[Beta]/2] E^(-I \[Alpha] /2)}};
      Transpose@m
      ];
FindOrder[mat_] := Module[{i,matn}, i = 1;
  If[mat == IdentityMatrix[3], Return[1]];
  matn=Det[mat] mat;
  (*matn=mat;*)
  While[
   MatrixPower[matn, i] != IdentityMatrix[3],
   i++
   ];
  i
  ];
transymmops[symmops_]:=Module[{order,sign},
  order=FindOrder/@(symminfo[[;;,2]]);
  sign=Table[First@First@MatrixPower[spinMatrix[Transpose[latt] . symminfo[[i,2]] . Inverse@Transpose@latt],order[[i]]],{i,Length[order]}];
Table[Which[
    order[[i]]==1,symmops[[i]],
    OddQ[order[[i]]](*True*), -sign[[i]] (*Det[symminfo[[i,2]]]*)(*(-1)^(1/order[[i]])*) symmops[[i]],
    EvenQ[order[[i]]], sign[[i]] (*Det[symminfo[[i,2]]]*) symmops[[i]]
    ]
  ,{i,Length[order]}]
]


getMSGElemFromMSGCorep[MSG_]:=Module[
{SG,brav,msgelem,latt},



If[Not@MemberQ[$Packages,"MSGCorep`"],Print["Please install and import MSGCorep package (https://github.com/goodluck1982/MSGCorep)";Abort[]];];
msgelem=MapAt[If[#==0,"F","T"]&,getMSGElem[MSG],{;;,-1}];
SG=MSG[[1]];
brav=getSGLatt[SG];
latt=BasicVectors[brav]/.{
SpaceGroupIrep`a->MagneticTB`a,
SpaceGroupIrep`b->MagneticTB`b,SpaceGroupIrep`c->MagneticTB`c,
SpaceGroupIrep`\[Alpha]->MagneticTB`\[Alpha],
SpaceGroupIrep`\[Beta]->MagneticTB`\[Beta],
SpaceGroupIrep`\[Gamma]->MagneticTB`\[Gamma]};
(*Print["Magnetic space group (BNS):",MSGSymStd[MSG]," No. ",StringRiffle[MSG,"."]];
Print["Lattice: ",brav];*)
Print["Primitive Lattice Vactor: ",latt];
msgelem=Insert[#,getRotMat[brav,#[[1]]],2]&/@msgelem(*;
SortBy[msgelem,#[[-1]]&]*)

];

findLittleGroupOfK[k_]:=Module[{rk,k0},
k0=Mod[k,1,0];

rk=Association@MapIndexed[First@#2->#1&,Insert[#,Inverse@Transpose[#[[2]]],5]&/@symminfo];
(*Print[rk];*)
rk=MapAt[Mod[# . k0,1,0]&,rk,{;;,5}];
(*Print[rk];*)
rk=Select[rk,#[[5]]==k0&&#[[4]]=="F"&];
Keys[rk]
];



getTBBandCorep[MSG_, ham_, param_, kset_] := Module[
  {rot,tmp, ops, wc, sym, tr, brav, msgele,mpgele,srot,srottb, eiv, little, trace, cr, coeff,
   U, opII, actk,
   order
   },
  tr = Association[];
  ops = N@symmetryops;
  (*ops=transymmops[ops];*)
  (*Print[ops];*)
  wc = N@wcc;
  sym = symminfo;
  tr["nelec"] = Length[ham];
  (*Print[basisdict[#]&@basis[[1,1]]];*)
  tr["soc"] =If[ListQ[basisdict[#]&@basis[[1,1]]],1,0];
  If[tr["soc"]==1,ops=transymmops[ops];];
  (*Print[tr["soc"]];*)
  (*basis*)
  tr["nsym"] = Length[sym];
  brav = getSGLatt[MSG[[1]]];
(*  msgele = getMSGElem[MSG];
  rot=getRotMat[brav, #[[1]]] & /@ msgele;*)
  rot=sym[[;;,2]];
  tr["rot"] = rot;
  tr["trans"] = N@sym[[;;,3]];
  (*tr["srot"] = N@getSpinRotOp[#[[1]]][[1]] & /@ msgele;*)
  (*Print[rot];*)
  (*latt={{0,-a,0},{(Sqrt[3] a)/2,a/2,0},{0,0,c}};*)
  srottb= N@spinMatrix[Transpose[latt] . # . Inverse@Transpose@latt] & /@ rot;
  
  (*msgele = getMSGElem[MSG];
  tr["srot"] = PauliMatrix[1].N@getSpinRotOp[#[[1]]][[1]].PauliMatrix[1] & /@ msgele;
  Print[getRotName[brav,#]&/@rot];*)
  
  mpgele=getRotName[brav,#]&/@rot;
  srot=PauliMatrix[1] . N@getSpinRotOp[#][[1]] . PauliMatrix[1] & /@ mpgele;
  tr["srot"] = srot;
  (*Print[MatrixForm/@srot];*)
  (*Print[Table[MatrixForm/@{Chop@srot[[i]],PauliMatrix[1].srottb[[i]].PauliMatrix[1]},{i,Length[rot]}]];*)
  (*Print[MatrixForm/@tr["srot"]];*)
  (*Abort[];*)
  tr["unitary"] = If[# == "F", 1, -1] & /@ sym[[;; , -1]];
  tr["nk"] = Length@kset;
  tr["kpt"] = kset;
  tr["nband"] = Length[ham];
  
  Table[tr[i] = {}, {i, {"ene", "deg", "knsym", "kisym", "trace"}}];
  Do[
   U = N@DiagonalMatrix[Table[Exp[-2 Pi I kpoint . tau], {tau, wc}]];
   (*Print[MatrixForm@U];*)
   eiv = Eigensystem[
     ConjugateTranspose[
       U] . (ham /. 
        Join[param, Thread[{kx, ky, kz} -> 2 Pi  kpoint]]) . U];
   eiv = Transpose@SortBy[Transpose[eiv], #[[1]] &];
   eiv = Transpose /@ SplitBy[Transpose[eiv], Round[#[[1]], 0.0001] &];
   (*Print[eiv];*)
   AppendTo[tr["ene"], Flatten@eiv[[;; , 1]]];
   AppendTo[tr["deg"], 
    Flatten[Table[#, {i, #}] & /@ Length /@ eiv[[;; , 1]]]];
   little = findLittleGroupOfK[kpoint];
   AppendTo[tr["knsym"], Length@little];
   AppendTo[tr["kisym"], little];
   
   trace =
    Flatten[Table[
      Table[
       Table[
        actk = Inverse[Transpose[symminfo[[i, 2]]]];
        (*opII is from GB Liu's note*)
        opII = 
         Exp[-2 Pi I symminfo[[i, 3]] . (actk . (kpoint))]  Table[
           Exp[2 Pi I actk . 
              kpoint . (wc[[m]] - (symminfo[[i, 2]] . wc[[l]]))], {m, 
            Length[wc]}, {l, Length[wc]}] ops[[i]];
        Chop@Tr[Conjugate[e[[2]]] . (opII) . Transpose[e[[2]]]]
        , {i, little}]
       , {nr, Length[e[[1]]]}]
      , {e, eiv}]];
   
   AppendTo[tr["trace"], Partition[trace, Length[little]]];
   
   , {kpoint, kset}];
  (*Print[tr];*)
  cr = getBandCorep[MSG, tr];
(*  Print["Magnetic space group (BNS): ", MSGSymStd[MSG]," No. ",StringRiffle[MSG,"."]];
  Print["Lattice: ",brav];*)
  (*Print[cr];*)
  Do[
  Print["k-name: "<>#2<>", k-point: "<>ToString[InputForm@#1]<>", little co-group: "<>ToString@#4 &@@cr["kinfo"][[k]]];
  Print@showBandCorep[cr, k]  
  
  ,{k,Length[kset]}]
  ]


End[]
EndPackage[]

