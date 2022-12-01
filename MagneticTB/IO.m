(* ::Package:: *)

BeginPackage["MagneticTB`"]




Begin["`Private`"]



pythtbIO2d[f_,wc_,rule_,latt_]:=Module[{wcc,
nband,toexp,onsiteham,onsitehopping,neighhopping,tmp1,tmp2,a,b,
hopping,realh,pytb,head
},
wcc=wc[[;;,;;2]];
nband=Length[f];
toexp=Expand@TrigToExp[f];
onsiteham=toexp/.E^__->0;
onsitehopping=Flatten[Table[{0,0,i,j,onsiteham[[i,j]]},{i,nband},{j,nband}],1];
(*Print[onsitehopping];*)
onsitehopping=Select[onsitehopping,Not[#[[-1]]===0]&];
(*Print[onsitehopping];*)
neighhopping=Map[Values[(GroupBy[Cases[#,x_ E^y_:>{x,Expand[y]}],Last])]&,toexp,{2}];
(*Print[neighhopping];*)
neighhopping=Map[{Total[#[[;;,1]]],#[[1,2]]}&,neighhopping,{3}];
neighhopping=Table[
tmp1=neighhopping[[i,j]];
Table[
tmp2=Values[First@SolveAlways[({a,b}-wcc[[i]]+wcc[[j]]) . {I kx,I ky}==bond[[2]],{kx,ky}]];
Flatten[{tmp2,i,j,bond[[1]]}]
,{bond,tmp1}]
,{i,nband},{j,nband}];
neighhopping=Flatten[neighhopping,2];
hopping=Join[onsitehopping,neighhopping];
head=StringRiffle[(StringReplace[ToString[N@#],"->"->"="]&/@rule),";"]<>"; I = 1.j; \n"<>
"magtb_model.set_onsite(["<>StringRiffle[Diagonal[onsiteham],","]<>"])\n";

pytb=
If[#[[3]]!=#[[4]]||#[[;;2]]!={0,0},"magtb_model.set_hop("<>ToString[N@#[[-1]],InputForm]<>","<>ToString[#[[3]]-1]<>","
<>ToString[#[[4]]-1]<>",["<>ToString[#[[1]]]<>","<>ToString[#[[2]]]<>"],allow_conjugate_pair=True)\n",""]&/@hopping;
pytb=head<>StringRiffle[pytb,""]
];

plotslab[wc_,latt_,shape_]:=Module[{wcc,lattice,atompos,pos,nfind,
shapepos,plotpos,posham,
nband,toexp,onsiteham,onsitehopping,neighhopping,tmp1,tmp2,a,b,
hopping,realh
},
nfind=30;
wcc=wc[[;;,;;2]];
lattice=latt[[;;2,;;2]];
(*atompos=Flatten[Table[{{i,j},#+{i,j}&/@wcc},{i,-nfind,nfind},{j,-nfind,nfind}],1];*)
atompos=Flatten[Table[{{i,j,#[[1]],#[[2]]}&/@({#,wcc[[#]]}&/@Range[Length[wcc]])},{i,-nfind,nfind},{j,-nfind,nfind}],3];
shapepos=Select[atompos,shape@@((#[[;;2]]+#[[4]]) . lattice)&];
plotpos=((#[[;;2]]+#[[4]]) . lattice)&/@N@shapepos

]


(*(*wannier90_hr.dat
*)
hop[f_,wc_,rule_]:=Module[{toexp,hopp,hoppm,onsite,select,hr,band,a,b,c,hopset,hopsetM,nums,sstr,result,
hrsim,rehrsim,hk,pij},
band=Length[f];(*Print[band];*)
toexp=Expand@TrigToExp[f];
(*Print[toexp[[1,2]]];*)
hopp=Flatten[Table[Cases[Collect[toexp[[i,j]],E^__],
xx_ E^yy_:>{Sequence@@First[SolveAlways[yy/I==(-wc[[i]]+wc[[j]]+{a,b,c}).{kx,ky,kz},{kx,ky,kz}]][[;;,2]],i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
(*Print[hopp];*)
hoppm=Flatten[Table[Cases[{Collect[toexp[[i,j]],E^__]},
xx_ E^yy_:>{Sequence@@First[SolveAlways[yy/I==(-wc[[i]]+wc[[j]]+{a,b,c}).{kx,ky,kz},{kx,ky,kz}]][[;;,2]],i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
(*hopp=SortBy[hopp,{#[[{1,2,3,5,4}]]&}]*)
onsite=Flatten[Table[Cases[{toexp[[i,j]]/.{E^_->0}},xx_Symbol:>{0,0,0,i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
(*Print[onsite];*)
(*Print[onsite];*)
hopp=SortBy[Join[hopp,hoppm,onsite],{#[[{1,2,3,4,5}]]&}];

hopset=Union[hopp[[;;,1;;3]]];
hopsetM=Union[hopp[[;;,1;;5]]];
Do[select[Sequence@@i,m,n]={Sequence@@i,m,n,0,0},{i,hopset},{m,band},{n,band}];
Do[select[Sequence@@i[[1;;5]]]=i/.N@rule,{i,hopp}];
result=Flatten[Table[select[Sequence@@i,m,n],{i,hopset},{n,band},{m,band}],2];
(*Print[result];*)
pij[i_,j_]:=-wc[[i]]+wc[[j]];
hrsim={Take[#,3],Take[#,{4,5}],#[[6]]+I #[[7]]}&/@result;
rehrsim={#[[1]]+pij[#[[2]][[1]],#[[2]][[2]]],#[[2]],#[[3]]}&/@hrsim;
hk=Normal[Total[SparseArray[#[[2]]->#[[3]]Exp[I #[[1]].{kx,ky,kz}],{band,band}]&/@rehrsim]];
(*If[Normalize[hk-f]\[NotEqual]0,Print[hk-f]];*)
If[Norm@Simplify[hk-(f/.rule)]!=0,Print["Error"];Abort[],Print["pass"],
Print[(*Norm@*)Norm@Chop@Simplify[hk-(f/.rule)]]];
nums=Length[hopset];
sstr=If[Mod[nums,15]==0,"\n"<>StringJoin[Table[StringJoin[Table["    1",15]]<>"\n",Quotient[nums,15]]],
"\n"<>StringJoin[Table[StringJoin[Table["    1",15]]<>"\n",Quotient[nums,15]],
StringJoin[Table["    1",Mod[nums,15]]]]<>"\n"];

"WannierHR\n"<>ToString[band]<>"\n"<>ToString[Length[hopset]]<>sstr
<>StringJoin[Table[StringRiffle[ToString[#,InputForm]&/@i,"    "]<>"\n",{i,result}]]

];*)


Options[hop] = {"wcc"->None,hrExport->None};
hop[f_,rule_, OptionsPattern[]]:=Module[{toexp,hopp,hoppm,onsite,select,hr,band,a,b,c,hopset,hopsetM,nums,sstr,result,
hrsim,rehrsim,hk,pij,
integerstr,
realstr,hrstr,
wc,export},
(*Print[OptionValue[atomPos],wcc];*)
If[OptionValue["wcc"]===None,wc=wcc,wc=OptionValue["wcc"]];
(*Print[wc];*)
(*wc=wcc;*)
band=Length[f];(*Print[band];*)
(*toexp=(*Expand@*)TrigToExp[f];*)
toexp=f;
(*Print[toexp];*)
(*Print[toexp/.{E^_->0}];*)
(*hopp=Flatten[Table[Cases[Collect[toexp[[i,j]],E^__],*)
hopp=Flatten[Table[Cases[toexp[[i,j]],
xx_ E^yy_:>{Sequence@@First[SolveAlways[yy/I==(-wc[[i]]+wc[[j]]+{a,b,c}) . {kx,ky,kz},{kx,ky,kz}]][[;;,2]],i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
(*Print[hopp];*)
hoppm=Flatten[Table[Cases[{Collect[toexp[[i,j]],E^__]},
xx_ E^yy_:>{Sequence@@First[SolveAlways[yy/I==(-wc[[i]]+wc[[j]]+{a,b,c}) . {kx,ky,kz},{kx,ky,kz}]][[;;,2]],i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
(*hopp=SortBy[hopp,{#[[{1,2,3,5,4}]]&}]*)
onsite=Flatten[Table[Cases[{toexp[[i,j]]/.{E^_->0}},xx_Symbol:>{0,0,0,i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
onsite=Flatten[Table[Cases[{toexp[[i,j]]/.{E^_->0}},xx_:>{0,0,0,i,j,ComplexExpand@Re@xx,ComplexExpand@Im@xx}],{i,band},{j,band}],2];
(*Print[onsite];*)
onsite=Select[onsite,Norm[#[[-2;;]]]=!=0&];
(*Print[onsite];*)
(*Print[onsite];*)
hopp=SortBy[Join[hopp,hoppm,onsite],{#[[{1,2,3,4,5}]]&}];

hopset=Union[hopp[[;;,1;;3]]];
hopsetM=Union[hopp[[;;,1;;5]]];
Do[select[Sequence@@i,m,n]={Sequence@@i,m,n,0,0},{i,hopset},{m,band},{n,band}];
Do[select[Sequence@@i[[1;;5]]]=i/.N@rule,{i,hopp}];
result=Flatten[Table[select[Sequence@@i,m,n],{i,hopset},{n,band},{m,band}],2];
(*Print[result];*)
(*Print[result];*)
pij[i_,j_]:=-wc[[i]]+wc[[j]];
hrsim={Take[#,3],Take[#,{4,5}],#[[6]]+I #[[7]]}&/@result;
rehrsim={#[[1]]+pij[#[[2]][[1]],#[[2]][[2]]],#[[2]],#[[3]]}&/@hrsim;
hk=Normal[Total[SparseArray[#[[2]]->#[[3]]Exp[I #[[1]] . {kx,ky,kz}],{band,band}]&/@rehrsim]];
(*If[Normalize[hk-f]\[NotEqual]0,Print[hk-f]];*)
(*If[Norm@Simplify[hk-(f/.rule)]!=0,Print["Error"];Abort[],Print["pass"],
Print[(*Norm@*)Norm@Chop@Simplify[hk-(f/.rule)]]];*)
nums=Length[hopset];
sstr=If[Mod[nums,15]==0,"\n"<>StringJoin[Table[StringJoin[Table["    1",15]]<>"\n",Quotient[nums,15]]],
"\n"<>StringJoin[Table[StringJoin[Table["    1",15]]<>"\n",Quotient[nums,15]],
StringJoin[Table["    1",Mod[nums,15]]]]<>"\n"];
integerstr[n_]:=ToString@NumberForm[n,5,NumberPadding->{" ",""}];
realstr[n_]:=ToString@NumberForm[N@n,{12,8},NumberPadding->{" ","0"},ExponentFunction->(Null&)];
res="Generated by MagneticTB\n"<>ToString[band]<>"\n"<>ToString[Length[hopset]]<>sstr
<>StringRiffle[Table[StringJoin@Flatten[{integerstr/@(#[[;;5]]),realstr/@(#[[6;;]])}&@line],
{line,result}],"\n"];
(*StringJoin[Table[StringRiffle[ToString[#,InputForm]&/@i,"    "]<>"\n",{i,result}]]*)
If[OptionValue[hrExport]=!=None,Export[FileNameJoin[{OptionValue[hrExport],"Wannier90_hr.dat"}],res,"Text"],Print[res]]
];


End[]
EndPackage[]



