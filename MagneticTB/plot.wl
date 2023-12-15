(* ::Package:: *)

BeginPackage["MagneticTB`"]




Begin["`Private`"]






bandManipulate[pathstr_,npoint_,h_]:=Module[
  {h0,params,mparams,m,rule,mpstr,klist},
  h0=TrigToExp[h];
  params=Variables[h0];
  mparams=ToExpression["p"<>ToString@Hash[#,"SHA"]<>ToString[#]&/@params];
  klist=2.Pi Flatten[Subdivide[#[[1]],#[[2]],npoint]&/@(Transpose[pathstr][[1]]),1];
  Print["Number of params:",Length@params];
  Print["params:",params];
(*  Print["mparams:",mparams];*)
  rule=ToString[Thread[params->mparams],StandardForm](*~Join~{kx\[Rule]k[[1]],ky\[Rule]k[[2]],kz\[Rule]k[[3]]}*);
(*  m=StringTake[ToString[{{#,0,StringTake[ToString[#],2;;]},-1,1}&/@mparams],{2;;-2}];*)
  m=StringTake[ToString[{{#1,0,#2},-1,1}&@@@Transpose[{mparams,params}]],{2;;-2}];
(*  Print[mparams,m];*)
  mpstr="Manipulate[\[IndentingNewLine]ListPlot[Transpose@Table[Eigenvalues[Evaluate["<>ToString[h,StandardForm]<>"/."<>rule<>"~Join~{kx\[Rule]k\[LeftDoubleBracket]1\[RightDoubleBracket],ky\[Rule]k\[LeftDoubleBracket]2\[RightDoubleBracket],kz\[Rule]k\[LeftDoubleBracket]3\[RightDoubleBracket]}]],{k,"<>ToString[klist]<>"}],
PlotRange->All,PlotStyle->Black],"<>m<>"
,Button[\"ExportData\",Print["<>rule<>"]]

]";
(*Print[m];
Print[mpstr];*)
  ToExpression[mpstr]
];





Options[bandplot]= {plotRange->All};
bandplot[pathstr_,npoint_,ham_,rule_,OptionsPattern[]]:=Module[{tmp,kps,hkps,hkp,xticks,yticks,data},

  hkp={};
(*  hkps=Partition[Partition[StringSplit[pathstr],5][[;;,-1]],2];*)
  hkps=Transpose[pathstr][[2]];
(*Print[hkps];*)
  font="KOZGOPRO-EXTRALIGHT";
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
  Transpose@Table[Sort@Eigenvalues[Evaluate[(N@ham)/.rule/.{kx->kps[[k,1]],ky->kps[[k,2]],kz->kps[[k,3]]}]],{k,Length[kps]}]
    ,{i,Length[tmp]}];
  data=Chop@Flatten[MapIndexed[{npoint(#2[[1]]-1)+#2[[3]]-1,#1}&,data,{3}],1];
(*Print[data];*)
(*data=Transpose@Table[Sort@Eigenvalues[Evaluate[ham/.rule/.{kx->k[[1]],ky->k[[2]],kz->k[[3]]}]],{k,kps}];*)
  xticks=Transpose@{(npoint) Range[0,Length[hkp]-1],Style[#,Black,FontFamily->font,24]&/@hkp};
  {maxenergy,minenergy}={Max[Flatten[data[[;;,;;,2]]]],Min[Flatten[data[[;;,;;,2]]]]};
  maxenergy=Max[Abs/@{maxenergy,minenergy}];
(*  {maxenergy,minenergy}={.6,-.3};*)
   (*Print[minenergy,maxenergy];*)
  yticks={#,
    Style[Round[#,.1],Black,FontFamily->font,24]}&/@Subdivide[-maxenergy,0,3];
  yticks=Join[yticks,Drop[{#,
    Style[Round[#,.1],Black,FontFamily->font,24]}&/@Subdivide[maxenergy,0,3],-1]] ;
  ListLinePlot[data,
(*PlotRange->OptionValue[PlotRange],*)
    PlotRange->{{0,(npoint)(Length[hkp]-1)},OptionValue[plotRange](*{-.02,.02}*)},
    PlotStyle->Black,
    GridLines->{(npoint) Range[Length[hkp]],{0}},
    Frame->{{True,True},{True,True}},
    FrameStyle->True,
    FrameTicks->{{yticks,None},{xticks,None}},
    FrameTicksStyle->Directive[Black,24],
    GridLinesStyle->Directive[Black]
    ]
  ];




End[]
EndPackage[]

