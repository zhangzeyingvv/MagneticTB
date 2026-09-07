(* ::Package:: *)

(*BeginPackage["MagneticTBDevelop`"]
*)



(*Begin["`Private`"]*)



Begin["MagneticTB`Private`"]

ClearAll[
  wilsonLoopFiniteNumberQ,
  wilsonLoopRealVectorQ,
  wilsonLoopNumericMatrix,
  wilsonLoopSortedFrame,
  wilsonLoopUnitaryPart,
  wilsonLoop
];

wilsonLoopFiniteNumberQ[value_] :=
  NumericQ[value] &&
    FreeQ[N[value], Indeterminate | ComplexInfinity | DirectedInfinity];

wilsonLoopRealVectorQ[vector_] :=
  VectorQ[vector, wilsonLoopFiniteNumberQ] &&
    TrueQ[Max[Abs[Im[N[vector]]]] == 0];

wilsonLoop::path =
  "The start and end points must be numeric vectors of equal positive length, and end-start must be a reciprocal-lattice vector (2 Pi times integers).";
wilsonLoop::centers =
  "Wannier centers must contain one real coordinate vector per Hamiltonian row, with the same coordinate dimension as the path.";
wilsonLoop::occ =
  "The occupied-band count must be an integer between 1 and the Hamiltonian dimension; received `1`.";
wilsonLoop::ham =
  "The Hamiltonian at path point `1` must be a finite numeric square Hermitian matrix of constant dimension.";
wilsonLoop::gap =
  "The occupied subspace is not isolated at path point `1`; the direct gap `2` does not exceed the tolerance `3`.";
wilsonLoop::covariance =
  "The endpoint Hamiltonians do not obey the tight-binding reciprocal covariance h(k+G)=V(G)^(-1).h(k).V(G); residual `1` exceeds tolerance `2`.";
wilsonLoop::overlap =
  "Occupied subspaces at adjacent/closure path positions have a singular overlap (minimum singular value `1` at link `2`); increase PathSubdivisions or check the isolated-band input.";
wilsonLoop::option = "Invalid Wilson-loop option value(s): `1`.";

Options[wilsonLoop] = {
  "PathSubdivisions" -> 50,
  "HermitianTolerance" -> 10^-10,
  "GapTolerance" -> 10^-9,
  "CovarianceTolerance" -> 10^-8,
  "OverlapTolerance" -> 10^-10,
  "Output" -> "PhaseOverPi"
};

wilsonLoopNumericMatrix[hamiltonian_, point_] := Module[{matrix},
  matrix = Quiet@Check[N[hamiltonian[point]], $Failed];
  If[
    matrix === $Failed || !MatrixQ[matrix, NumericQ] ||
      Length[Dimensions[matrix]] =!= 2 ||
      Dimensions[matrix][[1]] =!= Dimensions[matrix][[2]] ||
      Dimensions[matrix][[1]] < 1 ||
      !FreeQ[matrix, Indeterminate | ComplexInfinity | DirectedInfinity],
    $Failed,
    Normal[matrix]
  ]
];

wilsonLoopSortedFrame[matrix_, occupied_Integer] := Module[
  {values, vectors, order},
  {values, vectors} = Quiet@Check[Eigensystem[matrix], {$Failed, $Failed}];
  If[values === $Failed || !VectorQ[values, NumericQ] || !MatrixQ[vectors, NumericQ],
    Return[$Failed]
  ];
  order = Ordering[Re@values];
  {values[[order]], vectors[[order, ;;]][[;; occupied]]}
];

wilsonLoopUnitaryPart[matrix_, tolerance_] := Module[
  {left, singular, right, singularValues},
  {left, singular, right} = Quiet@Check[
    SingularValueDecomposition[matrix],
    {$Failed, $Failed, $Failed}
  ];
  If[left === $Failed, Return[$Failed]];
  singularValues = Diagonal[singular];
  If[singularValues === {} || Min[singularValues] <= tolerance,
    Return[{"Singular", If[singularValues === {}, 0, Min[singularValues]]}]
  ];
  {left . ConjugateTranspose[right], Min[singularValues]}
];

wilsonLoop[
    h_, centers_, occupied_, start_, end_,
    suppliedOptions : OptionsPattern[]
  ] := Module[
  {
    suppliedOptionList = {suppliedOptions}, allowedOptionNames,
    unknownOptionNames,
    subdivisions, hermitianTolerance, gapTolerance,
    covarianceTolerance, overlapTolerance, output,
    path, closureVector, reciprocalCoordinates, matrices, dimension,
    hermitianResiduals, endpointSewing, covarianceResidual,
    spectraAndFrames, spectra, frames, gaps, minimumGap,
    overlaps, closureOverlap, unitaryLinks, linkData,
    minimumSingularValue, wilsonMatrix, unitarityResidual,
    eigenvalues, phases, order, centersModuloOne, data
  },
  allowedOptionNames = First /@ Options[wilsonLoop];
  unknownOptionNames = Complement[
    First /@ suppliedOptionList,
    allowedOptionNames
  ];
  If[unknownOptionNames =!= {},
    Message[wilsonLoop::option, unknownOptionNames];
    Return[$Failed]
  ];
  subdivisions = OptionValue["PathSubdivisions"];
  hermitianTolerance = OptionValue["HermitianTolerance"];
  gapTolerance = OptionValue["GapTolerance"];
  covarianceTolerance = OptionValue["CovarianceTolerance"];
  overlapTolerance = OptionValue["OverlapTolerance"];
  output = OptionValue["Output"];
  If[
    !IntegerQ[subdivisions] || subdivisions < 1 ||
      !And @@ (wilsonLoopFiniteNumberQ[#] &&
        TrueQ[Im[N[#]] == 0] && TrueQ[Re[N[#]] >= 0] & /@ {
        hermitianTolerance, gapTolerance, covarianceTolerance,
        overlapTolerance
      }) ||
      !MemberQ[{"PhaseOverPi", "WannierCenters", "Eigenvalues", "Data"}, output],
    Message[wilsonLoop::option, {
      subdivisions, hermitianTolerance, gapTolerance,
      covarianceTolerance, overlapTolerance, output
    }];
    Return[$Failed]
  ];
  If[
    !wilsonLoopRealVectorQ[start] || !wilsonLoopRealVectorQ[end] ||
      Length[start] < 1 || Length[start] =!= Length[end],
    Message[wilsonLoop::path];
    Return[$Failed]
  ];
  closureVector = N[end - start];
  reciprocalCoordinates = closureVector/(2 Pi);
  If[Max[Abs[reciprocalCoordinates - Round[reciprocalCoordinates]]] > covarianceTolerance,
    Message[wilsonLoop::path];
    Return[$Failed]
  ];
  If[
    !ListQ[centers] || centers === {} ||
      !And @@ (wilsonLoopRealVectorQ[#] &&
        Length[#] === Length[start] & /@ centers),
    Message[wilsonLoop::centers];
    Return[$Failed]
  ];
  path = Subdivide[N[start], N[end], subdivisions];
  matrices = wilsonLoopNumericMatrix[h, #] & /@ path;
  If[MemberQ[matrices, $Failed],
    Message[
      wilsonLoop::ham,
      path[[First@FirstPosition[matrices, $Failed]]]
    ];
    Return[$Failed]
  ];
  dimension = Length[First[matrices]];
  If[Length[centers] =!= dimension,
    Message[wilsonLoop::centers];
    Return[$Failed]
  ];
  If[!IntegerQ[occupied] || !Between[occupied, {1, dimension}],
    Message[wilsonLoop::occ, occupied];
    Return[$Failed]
  ];
  If[!And @@ (Dimensions[#] === {dimension, dimension} & /@ matrices),
    Message[
      wilsonLoop::ham,
      path[[First@FirstPosition[
        Dimensions /@ matrices,
        dimensions_ /; dimensions =!= {dimension, dimension}
      ]]]
    ];
    Return[$Failed]
  ];
  hermitianResiduals = Max[Abs[Flatten[# - ConjugateTranspose[#]]]] & /@ matrices;
  If[Max[hermitianResiduals] > hermitianTolerance,
    Message[wilsonLoop::ham, path[[First@Ordering[hermitianResiduals, -1]]]];
    Return[$Failed]
  ];

  (* Alexandradinata-Dai-Bernevig use V(G)_aa=Exp[I G.r_a]. *)
  endpointSewing = DiagonalMatrix[Exp[I closureVector.#] & /@ N[centers]];
  covarianceResidual = Max@Abs@Flatten[
    Last[matrices] -
      ConjugateTranspose[endpointSewing] . First[matrices] . endpointSewing
  ];
  If[covarianceResidual > covarianceTolerance,
    Message[wilsonLoop::covariance, covarianceResidual, covarianceTolerance];
    Return[$Failed]
  ];

  spectraAndFrames = wilsonLoopSortedFrame[#, occupied] & /@ matrices;
  If[MemberQ[spectraAndFrames, $Failed],
    Message[
      wilsonLoop::ham,
      path[[First@FirstPosition[spectraAndFrames, $Failed]]]
    ];
    Return[$Failed]
  ];
  spectra = spectraAndFrames[[All, 1]];
  frames = spectraAndFrames[[All, 2]];
  gaps = If[
    occupied < dimension,
    Re[#[[occupied + 1]] - #[[occupied]]] & /@ spectra,
    ConstantArray[Infinity, Length[path]]
  ];
  minimumGap = Min[gaps];
  If[minimumGap <= gapTolerance,
    Message[
      wilsonLoop::gap,
      path[[First@FirstPosition[gaps, minimumGap]]],
      minimumGap,
      gapTolerance
    ];
    Return[$Failed]
  ];

  overlaps = MapThread[
    Conjugate[#2] . Transpose[#1] &,
    {Most[frames], Rest[frames]}
  ];
  closureOverlap =
    Conjugate[First[frames]] . endpointSewing . Transpose[Last[frames]];
  linkData = Map[
    wilsonLoopUnitaryPart[#, overlapTolerance] &,
    Append[overlaps, closureOverlap]
  ];
  If[MemberQ[linkData, $Failed],
    Message[wilsonLoop::ham, "SVD"];
    Return[$Failed]
  ];
  If[MemberQ[linkData, {"Singular", _}],
    With[{position = First@FirstPosition[linkData, {"Singular", _}]},
      Message[
        wilsonLoop::overlap,
        linkData[[position, 2]],
        If[position <= Length[overlaps], position, "closure"]
      ]
    ];
    Return[$Failed]
  ];
  unitaryLinks = linkData[[All, 1]];
  minimumSingularValue = Min[linkData[[All, 2]]];
  wilsonMatrix = Last[unitaryLinks] .
    (Dot @@ Reverse[Most[unitaryLinks]]);
  wilsonMatrix = Chop[wilsonMatrix, hermitianTolerance];
  unitarityResidual = Max@Abs@Flatten[
    ConjugateTranspose[wilsonMatrix] . wilsonMatrix -
      IdentityMatrix[occupied]
  ];
  eigenvalues = Eigenvalues[wilsonMatrix];
  phases = Arg[eigenvalues]/Pi;
  order = Ordering[phases];
  eigenvalues = eigenvalues[[order]];
  phases = Chop[phases[[order]], hermitianTolerance];
  centersModuloOne = Mod[Chop[phases/2, hermitianTolerance], 1];
  data = <|
    "Eigenvalues" -> eigenvalues,
    "PhasesOverPi" -> phases,
    "WannierCenters" -> centersModuloOne,
    "WilsonMatrix" -> wilsonMatrix,
    "UnitarityResidual" -> unitarityResidual,
    "MinimumLinkSingularValue" -> minimumSingularValue,
    "MinimumDirectGap" -> minimumGap,
    "Path" -> path,
    "PathSubdivisions" -> subdivisions,
    "OccupiedBands" -> occupied,
    "ClosureVector" -> closureVector,
    "ReciprocalCoordinates" -> Round[reciprocalCoordinates],
    "EndpointCovarianceResidual" -> covarianceResidual
  |>;
  Switch[output,
    "PhaseOverPi", data["PhasesOverPi"],
    "WannierCenters", data["WannierCenters"],
    "Eigenvalues", data["Eigenvalues"],
    "Data", data
  ]
];

End[]

tokpathvasp[pathstr_,npoint_]:=2Pi Flatten[Subdivide[#[[1]],#[[2]],npoint]&/@ToExpression@Partition[Partition[StringSplit[pathstr],5][[;;,1;;3]],2],1];
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


