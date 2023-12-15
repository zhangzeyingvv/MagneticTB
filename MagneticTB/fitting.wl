(* ::Package:: *)

BeginPackage["MagneticTB`"]




Begin["`Private`"]



bandManipulateEig[h_, klist_, eigdata_] := Module[
     {h0, params, mparams, m, rule, mpstr, plotrange},
     plotrange = {Min[#] - .2, Max[#] + .2} &@
     Flatten[Transpose@eigdata[[;; , 2]]];
     h0 = TrigToExp[h];
     params = Variables[h0];
     mparams = ToExpression["p" <> ToString[#] & /@ params];
     Print["Number of params:", Length@params];
     Print["params:", params];
     Print["mparams:", mparams];
     rule = ToString[Thread[params -> mparams], StandardForm](*~
   Join~{kx\[Rule]k[[1]],ky\[Rule]k[[2]],kz\[Rule]k[[3]]}*);
     m = StringTake[ToString[{{#, 0}, -1, 1} & /@ mparams], {2 ;; -2}];
     m = 
    StringTake[
     ToString[{{#1, 0, #2}, -1, 1} & @@@ 
       Transpose[{mparams, params}]], {2 ;; -2}];
     mpstr = 
    "Manipulate[\[IndentingNewLine]Show[ListPlot[Transpose@Table[\
Eigenvalues[Evaluate[" <> ToString[h, StandardForm] <> "/." <> rule <>
      "~Join~{kx\[Rule]k\[LeftDoubleBracket]1\[RightDoubleBracket],ky\
\[Rule]k\[LeftDoubleBracket]2\[RightDoubleBracket],kz\[Rule]k\
\[LeftDoubleBracket]3\[RightDoubleBracket]}]],{k," <> 
     ToString[klist] <> "}],
PlotRange->All,PlotStyle->Black],ListPlot[Transpose@" <> 
     ToString[eigdata[[;; , 2]]] <> ",PlotStyle->Red],PlotRange->" <> 
     ToString[plotrange] <> "]," <> m <> "
,Button[\"ExportData\",Print[" <> rule <> "]]

]";
   
   (*Print[mpstr];*)
     ToExpression[mpstr]
   ];


tokpathvasp[pathstr_,npoint_]:=2Pi Flatten[Subdivide[#[[1]],#[[2]],npoint]&/@ToExpression@Partition[Partition[StringSplit[pathstr],5][[;;,1;;3]],2],1];
vasptoPath={Rationalize[ToExpression@Partition[Partition[StringSplit[#],5][[;;,{1,2,3}]],2]],Partition[Partition[StringSplit[#],5][[;;,-1]],2]}\[Transpose]&;


fittingTB[h_, eigdata_, krange_, initparms_] := Module[
  {params, mparams, klist, hp, modelhams, nk, objparams, res, lsq},
  params = Complement[Variables[h], {kx, ky, kz}];
  klist = eigdata[[;; , 1]];
  nk = Length[krange];
  mparams = params;
  hp = h /. Thread[params -> mparams];
  modelhams = 
   Table[hp /. {kx -> k[[1]], ky -> k[[2]], kz -> k[[3]]}, {k, 
     2 Pi klist[[krange]]}];
  res = FindMinimum[
    Total@Table[
      EuclideanDistance[Sort@Eigenvalues[modelhams[[i]]], 
       eigdata[[;; , 2]][[krange[[i]]]]], {i, nk}], 
    Transpose[{mparams, params /. initparms}]];
  objparams = res[[2]] /. Thread[mparams -> params];
  
  lsq = {Total@
      Table[EuclideanDistance[Sort@Eigenvalues[modelhams[[i]]], 
        eigdata[[;; , 2]][[krange[[i]]]]], {i, nk}] /. (initparms /. 
       Thread[params -> mparams]), res[[1]]};
  Association["FittedParams" -> objparams, 
   "LSQForInitParams" -> lsq[[1]], "LSQForFittedParams" -> lsq[[2]]]
  ]


vaspEig[filename_,efermi_,spin_,startband_,endband_]:=Module[
  {nband,nk,fermi,eigfile,todata,band},
  eigfile=StringSplit[Import[filename],"\n"];
  eigfile=StringReplace[#,{"E+":>"*^","E-":>"*^-"}]&[eigfile];
  fermi=efermi;
  {nk,nband}=ToExpression[StringSplit[eigfile[[6]]][[-2;;-1]]];
  todata=ToExpression[StringSplit/@#[[2;;]]]&;
  band=Table[{#[[1]][[1;;3]],#[[2;;]][[;;,spin+1]]-fermi}&@todata@Partition[eigfile[[7;;-1]],nband+2][[i]],{i,1,nk,1}];
{#[[1]],#[[2]][[startband;;endband]]}&/@band
];
Options[compareBand]= {plotRange->All};
compareBand[pathstr_,npoint_,ham_,rule_,vaspband_,OptionsPattern[]]:=
  Module[{tmp,kps,hkps,hkp,xticks,yticks,font,data,oriband,norband,vaspplotdata,tbband},
  oriband=Partition[vaspband[[;;,2]],npoint];
  (*Print[oriband];*)
  norband=Table[If[i==1,oriband[[i]],Join[{oriband[[i-1]][[-1]]},oriband[[i]]]],{i,Length[oriband]}];
  (*Print[norband];*)
  If[Length[norband]>1,PrependTo[norband[[1]],norband[[1,1]]]];
  vaspplotdata=Flatten[MapIndexed[{(npoint)(#2[[1]]-1)+#2[[3]]-1,#1}&,Transpose/@norband,{3}],1];
  (*Print[vaspplotdata];*)
  hkp={};
  (*hkps=Partition[Partition[StringSplit[pathstr],5][[;;,-1]],2];*)
  hkps=Transpose[pathstr][[2]];
  font="Times";
  Do[If[i==1,AppendTo[hkp,hkps[[i]][[1]]];AppendTo[hkp,hkps[[i]][[2]]],
    If[hkps[[i]][[1]]==hkps[[i-1]][[2]],AppendTo[hkp,hkps[[i]][[2]]],
      AppendTo[hkp,hkps[[i]][[2]]];hkp[[i]]=hkp[[i]]<>"|"<>hkps[[i]][[1]];
      ]
    ]
  ,{i,Length[hkps]}];
  hkp=StringReplace[#,{"\\Gamma"->"\[CapitalGamma]"}]&/@hkp;
  (*  tmp=ToExpression@Partition[Partition[StringSplit[pathstr],5][[;;,1;;3]],2];*)
  tmp=Transpose[pathstr][[1]];
  data=Table[kps=2Pi Subdivide[tmp[[i,1]],tmp[[i,2]],npoint];
  Transpose@Table[Sort@Eigenvalues[Evaluate[ham/.rule/.{kx->kps[[k,1]],ky->kps[[k,2]],kz->kps[[k,3]]}]],{k,Length[kps]}]
  ,{i,Length[tmp]}];
  data=Chop@Flatten[MapIndexed[{npoint(#2[[1]]-1)+#2[[3]]-1,#1}&,data,{3}],1];
  xticks=Transpose@{(npoint) Range[0,Length[hkp]-1],Style[#,Black,FontFamily->font,24]&/@hkp};
  {maxenergy,minenergy}={Max[Flatten[data[[;;,;;,2]]]],Min[Flatten[data[[;;,;;,2]]]]};
  (*{maxenergy,minenergy}={.6,.1};*)
  yticks={#,
    Style[Round[#,.01],Black,FontFamily->font,24]}&/@Subdivide[-maxenergy,0,2];
  yticks=Join[yticks,Drop[{#,
      Style[Round[#,.01],Black,FontFamily->font,24]}&/@Subdivide[maxenergy,0,2],-1]] ;
  tbband=ListLinePlot[data,
      PlotRange->{{0,(npoint)(Length[hkp]-1)},OptionValue[plotRange](*{-.02,.02}*)},
      PlotStyle->Purple,
      GridLines->{(npoint) Range[Length[hkp]],{0}},
      Frame->{{True,True},{True,True}},
      FrameStyle->True,
      FrameTicks->{{yticks,None},{xticks,None}},
      FrameTicksStyle->Directive[Black,24],
      GridLinesStyle->Directive[Black]
    ];
  Labeled[Show[tbband,ListLinePlot[vaspplotdata,PlotStyle-> Lighter[Blue,.5]]],{"Red: TB, Blue: VASP"},{Top}]
];


End[]
EndPackage[]
