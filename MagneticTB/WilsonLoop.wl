(* ::Package:: *)

(*BeginPackage["MagneticTBDevelop`"]
*)



(*Begin["`Private`"]*)



wilsonLoop[h_,wcc_,occ_,start_,end_]:=Module[
{rnk1,path,evals,evecs,eigvec,hk,hkdag,eig0,evec0,
V,P,W,Wij
},
rnk1=50;
V=DiagonalMatrix[Table[Exp[-I #1 . tau],{tau,wcc}]]&;
path=Subdivide[N@end,N@start,rnk1];
(*Print[path];*)
{eig0,evec0}=Transpose[SortBy[Transpose[Eigensystem[N@h[path[[1]]]]],#[[1]]&]];
path=Drop[path,1];
eigvec=Table[
hk=N@h[path[[ik1]]];
{evals,evecs}=Transpose[SortBy[Transpose[Eigensystem[hk]],#[[1]]&]];
evecs
,{ik1,rnk1-1}];
P[n_]:=Sum[Transpose[{eigvec[[n,i]]}] . Conjugate@{eigvec[[n,i]]},{i,occ}];
W=Chop@(V[end-start] . Dot@@Table[P[rnk1-i],{i,rnk1-1}])

;
Wij=Chop@Table[(Conjugate@{evec0[[i]]} . W . Transpose[{evec0[[j]]}])[[1,1]],
{i,occ},{j,occ}];
(*Wij=W[[;;occ,;;occ]];*)
(*Wij=W;*)
(*Mod[-Arg@Eigenvalues[Wij]/(2Pi),1]*)
Sort@Arg@Eigenvalues[Wij]/(Pi)

];
tokpathvasp[pathstr_,npoint_]:=2Pi Flatten[Subdivide[#[[1]],#[[2]],npoint]&/@ToExpression@Partition[Partition[StringSplit[pathstr],5][[;;,1;;3]],2],1];
plotWilsonLoop[h_,wcc_,occ_,start_,end_,pathstr_,npoint_]:=Module[{hkp,hkps,xticks,data},

hkp={};
hkps=Partition[Partition[StringSplit[pathstr],5][[;;,-1]],2];
(*Print[hkps];*)
Do[If[i==1,AppendTo[hkp,hkps[[i]][[1]]];AppendTo[hkp,hkps[[i]][[2]]],
If[hkps[[i]][[1]]==hkps[[i-1]][[2]],AppendTo[hkp,hkps[[i]][[2]]],
AppendTo[hkp,hkps[[i]][[2]]];hkp[[i]]=hkp[[i]]<>"|"<>hkps[[i]][[1]];
]
]
,{i,Length[hkps]}];
hkp=StringReplace[#,{"\\Gamma"->"\[CapitalGamma]"}]&/@hkp;

(*Print[(Partition[tokpathvasp[pathstr,#],#+1]&@npoint)];*)
data=Transpose@Table[wilsonLoop[h,wcc,occ,start+i,end+i],
{i,#}]&/@(Partition[tokpathvasp[pathstr,#],#+1]&@npoint);
data=Flatten[MapIndexed[{npoint(#2[[1]]-1)+#2[[3]]-1,#1}&,data,{3}],1];
xticks=Transpose@{(npoint) Range[0,Length[hkp]-1],Style[#,Black,Italic,FontFamily->"Times",24]&/@hkp};
(*yticks=Min[data] 1.05;*)
(*Print[Max@data];*)
ListPlot[data,
PlotRange->{{0,(npoint)(Length[hkp]-1)},{-1-.1,1+.1}},
PlotStyle->Black,
GridLines->{(npoint) Range[Length[hkp]],{0}},
Frame->{{True,True},{True,True}},
FrameStyle->True,
FrameTicks->{{All,None},{xticks,None}},
FrameTicksStyle->Directive[Black,24],
GridLinesStyle->Directive[Black]
]

]


z2path[h_,occ_,path1_]:=Module[
{rnk1,evals,evecs,eigvec,hk,hkdag,
lambda0,overlap,u,s,v,lambda
},
rnk1=Length[path1]-1;

eigvec=Table[
hk=N@h[path1[[ik1]]];
{evals,evecs}=Transpose[SortBy[Transpose[Eigensystem[hk]],#[[1]]&]];
evecs
,{ik1,rnk1}];

lambda0= IdentityMatrix[occ];
Do[hkdag=eigvec[[ik1]];
If[ik1==rnk1,
hk=eigvec[[1]];
,
hk=eigvec[[ik1+1]];
];
overlap=Conjugate[hkdag[[;;occ]]] . Transpose[hk[[;;occ]]];
{u,s,v}=SingularValueDecomposition[overlap];
overlap=v . ConjugateTranspose[u];
lambda=overlap . lambda0;
lambda0=lambda
,{ik1,rnk1}];
Mod[-Arg@Eigenvalues[lambda]/(2Pi),1]

];
(*
ListPlot[Table[z2path[hw,2,Subdivide[{i,-N@Pi ,0},{i,N@Pi,0},20]],{i,Subdivide[-Pi,0,200]}]\[Transpose],PlotRange\[Rule]All,PlotStyle\[Rule]Black]
*)

(*Total[SparseArray[#[[4;;5]]\[Rule](#[[6]]+I#[[7]])Exp[I #[[;;3]].{kx,ky,kz}],{4,4}] *)


(*rnk1=Length[path1];
innereig=Table[
hk=N@h[kpoint1];
{evals,evecs}=Transpose[SortBy[Transpose[Eigensystem[hk]],#[[1]]&]]

,{kpoint1,path1}];
prd=IdentityMatrix[occ];
Do[
overlap=Conjugate[innereig[[neiv,2,;;occ]]].Transpose[innereig[[neiv+1,2,;;occ]]];
{v,s,w}=SingularValueDecomposition@overlap;
(*prd=prd.(svd[[1]].ConjugateTranspose[svd[[3]]])*)
(*prd=(v.ConjugateTranspose[u].overlap).prd*)
Print[MatrixForm@(v.ConjugateTranspose[w])];
prd=prd.v.s.ConjugateTranspose[v].(v.ConjugateTranspose[w])
,{neiv,rnk1-1}];
(*Chop@Total@Im@Log@Eigenvalues@(Conjugate[innereig[[1,2,;;occ]]].Transpose[innereig[[rnk1,2,;;occ]]]) /Pi*)
-Arg@Eigenvalues@(prd)(*/(2Pi)*)
];*)








berryph[h_,occ_,path1_]:=Module[
{rnk1,evals,evecs,hk,nband,
innereig,overlap,svd,prd},
rnk1=Length[path1];
innereig=Table[
hk=N@h[kpoint1];
{evals,evecs}=Transpose[SortBy[Transpose[Eigensystem[hk]],#[[1]]&]]

,{kpoint1,path1}];
Do[
overlap=Conjugate[innereig[[neiv,2,;;occ]]] . Transpose[innereig[[neiv+1,2,;;occ]]];
svd=SingularValueDecomposition@overlap;
innereig[[neiv+1,2,;;occ]]=svd[[3]] . ConjugateTranspose[svd[[1]]] . innereig[[neiv+1,2,;;occ]];
,{neiv,rnk1-1}];

(*Chop@Total@Im@Log@Eigenvalues@(Conjugate[innereig[[1,2,;;occ]]].Transpose[innereig[[rnk1,2,;;occ]]]) /Pi*)
Chop@Im@Log@Eigenvalues@(Conjugate[innereig[[1,2,;;occ]]] . Transpose[innereig[[rnk1,2,;;occ]]]) /Pi
];

(*
dsg51E=tmp=(Subscript[c, 1]+Subscript[c, 2] Subscript[k, x]) Subscript[\[CapitalGamma], 0,0]+Subscript[c, 3] Subscript[k, y] Subscript[\[CapitalGamma], 0,1]+Subscript[k, z] Subscript[c, 1,1] Subscript[\[CapitalGamma], 1,3]+Subscript[k, z] Subscript[c, 2,1] Subscript[\[CapitalGamma], 2,3]+Subscript[k, z] Subscript[c, 3,1] Subscript[\[CapitalGamma], 3,3]/.repall;
ham[{kx_,ky_,kz_}]:=Evaluate[#/.(MapThread[Rule,{#,RandomReal[{-10,10},Length[#]]}&@Complement[Variables[#],{kx,ky,kz}]])&[tmp]];
cp={1,0,0};\[Epsilon]=1/100;
cirpara=FindInstance[{{\[Alpha]1,\[Beta]1,\[Gamma]1}.cp\[Equal]0,{\[Alpha]2,\[Beta]2,\[Gamma]2}.cp\[Equal]0,{\[Alpha]1,\[Beta]1,\[Gamma]1}.{\[Alpha]2,\[Beta]2,\[Gamma]2}\[Equal]0,Norm[{\[Alpha]1,\[Beta]1,\[Gamma]1}]\[Equal]1,Norm[{\[Alpha]2,\[Beta]2,\[Gamma]2}]\[Equal]1(*,\[Alpha]1\[Equal]1./Pi*)},{\[Alpha]1,\[Beta]1,\[Gamma]1,\[Alpha]2,\[Beta]2,\[Gamma]2},Reals][[1]];
cirtable=Table[ \[Epsilon] Cos[\[Theta]]{\[Alpha]1,\[Beta]1,\[Gamma]1}+\[Epsilon] Sin[\[Theta]]{\[Alpha]2,\[Beta]2,\[Gamma]2}/.cirpara,{\[Theta],Subdivide[0,2.Pi,300]}];
berryf[ham,2,cirtable]
*)

chernsph[h_,occ_,paths_,coors_]:=Module[
{k,hk,hks,pathhams,overlap,
npoint,svd,chern},
hks={};
npoint=Length[paths[[1]]];
Do[
k=coor;
AppendTo[hks,Transpose[SortBy[Transpose[Eigensystem[h[k]]],#[[1]]&]]]
,{coor,coors}];
pathhams=hks[[#]]&/@paths;
chern=0;
Do[
Do[
overlap=Conjugate[pathham[[neiv,2,;;occ]]] . Transpose[pathham[[neiv+1,2,;;occ]]];
svd=SingularValueDecomposition@overlap;
pathham[[neiv+1,2,;;occ]]=svd[[3]] . ConjugateTranspose[svd[[1]]] . pathham[[neiv+1,2,;;occ]];
,{neiv,npoint-1}];
chern=chern+Total@Im@Log@Eigenvalues@(Conjugate[pathham[[1,2,;;occ]]] . Transpose[pathham[[npoint,2,;;occ]]]) ;
(*Print[chern];*)
,{pathham,pathhams}];
chern/(2Pi)
];

chirality[h_,occ_,point_,r_]:=Module[
{n\[Phi],n\[Theta],k,hk,\[Theta]phase,
\[Phi]Eigsys,overlap,svd},
n\[Phi]=10;
n\[Theta]=1000;
\[Theta]phase={};
Do[
\[Phi]Eigsys={};
Do[
k=N@point+{r Sin[\[Theta]]Cos[\[Phi]],r Sin[\[Theta]]Sin[\[Phi]], r Cos[\[Theta]]};
hk=h[k];
AppendTo[\[Phi]Eigsys,Transpose[SortBy[Transpose[Eigensystem[hk]],#[[1]]&]]];
,{\[Phi],Subdivide[0,2.Pi,n\[Phi]]}];
Do[
overlap=Conjugate[\[Phi]Eigsys[[neiv,2,;;occ]]] . Transpose[\[Phi]Eigsys[[neiv+1,2,;;occ]]];
svd=SingularValueDecomposition@overlap;
\[Phi]Eigsys[[neiv+1,2,;;occ]]=svd[[3]] . ConjugateTranspose[svd[[1]]] . \[Phi]Eigsys[[neiv+1,2,;;occ]];
,{neiv,n\[Phi]}];
AppendTo[\[Theta]phase,Mod[Total@Im@Log@Eigenvalues@(Conjugate[\[Phi]Eigsys[[1,2,;;occ]]] . Transpose[\[Phi]Eigsys[[n\[Phi]+1,2,;;occ]]]) ,2.Pi]];

,{\[Theta],Subdivide[0,1.Pi,n\[Theta]]}];
\[Theta]phase
];

(*
<<NDSolve`FEM`;
mesh=ToBoundaryMesh[Sphere[{0,0,0},1],MaxCellMeasure\[Rule].2,AccuracyGoal\[Rule]2];
paths=Append[#,#[[1]]]&/@(mesh["BoundaryElements"][[1,1]]);
coors=mesh["Coordinates"];

mesh=ToBoundaryMesh[Cylinder[{{0,0,-1},{0,0,1}},1],MaxCellMeasure\[Rule].2,AccuracyGoal\[Rule]2];
paths=Append[#,#[[1]]]&/@(mesh["BoundaryElements"][[1,1]]);
coors=mesh["Coordinates"];

ham[{x_,y_,z_}]:=(25-(x-0)^2+y^2+z^2)PauliMatrix[3]+z PauliMatrix[1];
torus=ImplicitRegion[((Sqrt[(x)^2+y^2]-5.0)/1)^2+z^2\[LessEqual]1,{x,y,z}];
mesh=ToBoundaryMesh[torus,MaxCellMeasure\[Rule]2,AccuracyGoal\[Rule]2];
paths=Append[#,#[[1]]]&/@(mesh["BoundaryElements"][[1,1]]);
coors=mesh["Coordinates"];
mesh["Wireframe"]

*)


wLoop[h_,occ_,path1_,path2_]:=Module[
{nk1=20,rnk1,nk2=1020,rnk2,evals,evecs,hk,path1k,path2k,kp,
innereig,overlap,wilson},
path1k=Join@@Join[{{path1[[1,1]]}},(Drop[Subdivide[#[[1]],#[[2]],nk1],1]&/@path1)];
rnk1=Length[path1k];
path2k=Join@@Join[{{path1[[1,1]]}},(Drop[Subdivide[#[[1]],#[[2]],nk2],1]&/@path2)];
rnk2=Length[path2k];
(*Print[path2k];*)
Table[
innereig=Table[
kp=kpoint1+kpoint2;
hk=N@h[kp];
{evals,evecs}=Transpose[SortBy[Transpose[Eigensystem[hk]],#[[1]]&]],{kpoint1,path1k}];
overlap=Table[Conjugate[innereig[[neiv,2,i]]] . innereig[[neiv+1,2,j]],{neiv,rnk1-1},{i,1,occ},{j,1,occ}];
(*Print[overlap];*)
wilson=IdentityMatrix[occ];

Do[wilson=wilson . overlap[[i]],{i,rnk1-1}];
Arg@Eigenvalues[wilson]
(*Arg@Det[wilson]*)
(*Arg@Eigenvalues[#[[1]].ConjugateTranspose[#[[3]]]&@SingularValueDecomposition[wilson]]*)
,
{kpoint2,path2k}]];


wLoopmirr[h_,occ_,path1_,path2_,mirr_]:=Module[
{nk1=20,rnk1,nk2=200,rnk2,hk,path1k,path2k,kp,
innereig,overlap,wilsonup,wilsondn,
evalm,evecm,hkup,hkdn,dim,evalsup,evecsup,evalsdn,evecsdn},
{evalm,evecm}=Transpose[SortBy[Transpose[Eigensystem[mirr]],#[[1]]&]];
dim=Length[mirr]/2;
(*Print[evalm,evecm];
Print[#[[;;occ,;;occ]],#[[occ+1;;,occ+1;;]]&@FullSimplify[evecm.h[{a,b,c}].Inverse[evecm]]];*)
path1k=Join@@Join[{{path1[[1,1]]}},(Drop[Subdivide[#[[1]],#[[2]],nk1],1]&/@path1)];
rnk1=Length[path1k];
path2k=Join@@Join[{{path1[[1,1]]}},(Drop[Subdivide[#[[1]],#[[2]],nk2],1]&/@path2)];
rnk2=Length[path2k];
(*Print[path2k];*)
Table[
innereig=Table[
kp=kpoint1+kpoint2;
hk=evecm . N@h[kp] . Inverse[evecm];
hkup=hk[[;;dim,;;dim]];
hkdn=hk[[dim+1;;,dim+1;;]];
{evalsup,evecsup}=Transpose[SortBy[Transpose[Eigensystem[hkup]],#[[1]]&]];
{evalsdn,evecsdn}=Transpose[SortBy[Transpose[Eigensystem[hkdn]],#[[1]]&]];
{{evalsup,evecsup},{evalsdn,evecsdn}},{kpoint1,path1k}];
overlap=Table[Conjugate[innereig[[neiv,1,2,i]]] . innereig[[neiv+1,1,2,j]],{neiv,rnk1-1},{i,1,occ/2},{j,1,occ/2}];
(*Print[overlap];*)
wilsonup=IdentityMatrix[occ/2];

Do[wilsonup=wilsonup . overlap[[i]],{i,rnk1-1}];
wilsonup=Arg@Eigenvalues[wilsonup];

overlap=Table[Conjugate[innereig[[neiv,2,2,i]]] . innereig[[neiv+1,2,2,j]],{neiv,rnk1-1},{i,1,occ/2},{j,1,occ/2}];
(*Print[overlap];*)
wilsondn=IdentityMatrix[occ/2];
Do[wilsondn=wilsondn . overlap[[i]],{i,rnk1-1}];
wilsondn=Arg@Eigenvalues[wilsondn];
{wilsonup,wilsondn}
,
{kpoint2,path2k}]];


bandPlot[h_,path1_,nk1_]:=Module[{path1k,rnk1},
path1k=Join@@Join[{{path1[[1,1]]}},(Drop[Subdivide[#[[1]],#[[2]],nk1],1]&/@path1)];
rnk1=Length[path1k];
Table[Eigenvalues[N@h[kp]],{kp,path1k}]\[Transpose]
];


(*End[]
EndPackage[]
*)


