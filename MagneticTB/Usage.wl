(* ::Package:: *)

(* Forward declarations to prevent shadowing messages *)


BeginPackage["MagneticTB`"]


(* ::Subsubsection:: *)
(*core*)


(*set input data*)
init ::usage = "Init the program"
unsymham ::usage = "Get the unsymmetried Hamiltonian"
symham ::usage = "Get the symmetried Hamiltonian in convenition I"

kx ::usage = ""
ky ::usage = ""
kz ::usage = ""
a ::usage = ""
b ::usage = ""
c ::usage = ""
\[Alpha] ::usage = ""
\[Beta] ::usage = ""
\[Gamma] ::usage = ""
x ::usage = ""
y ::usage = ""
z ::usage = ""


(*output data*)
ops  ::usage ="ops"
braLatt ::usage = ""
gray ::usage = ""
typeI ::usage = ""
typeIII ::usage = ""
typeIV ::usage = ""
bnsdict ::usage = ""
ogdict ::usage = ""
ognumdict ::usage = ""
bondclassify  ::usage = ""
symminfo ::usage = ""
symmetryops  ::usage = ""
symmetryopsII   ::usage = ""
atompos ::usage=""
reclatt ::usage=""
latt ::usage=""
wcc ::usage = "Wannier center"
symmcompile ::usage = ""
compactForm ::usage = ""
lattice::usage = "Options for init"
wyckoffposition::usage = "Options for init"
symminformation::usage = "Options for init"
basisFunctions::usage = "Options for init"
symmetryset::usage = "Options for symham"
lattpar::usage = "Options for init"


(* ::Subsubsection:: *)
(*Plot*)


(*Plot*)
bandManipulate ::usage=""
bandManipulateEig ::usage=""
bandplot ::usage=""
vaspEig ::usage=""
compareBand ::usage=""


(*Functions*)
readMsgData ::usage = "Read the magnetic space group data"
pointMatrix  ::usage = ""
tokpathvasp ::usage = ""


realham2d ::usage=""



(* ::Subsubsection:: *)
(*IO*)


(*
Interface
*)
pythtbIO2d ::usage = "init"
plotslab ::usage = ""
hop  ::usage =""


(*Wilson loop (testing)*)
wLoop ::usage = "init"
berryph ::usage = ""
z2path ::usage =""
wilsonLoop ::usage =""
plotWilsonLoop ::usage =""



msgop ::usage =""
texOutput ::usage =""





(* ::Subsubsection:: *)
(*Utilities*)


symhamII ::usage = "Get the symmetried Hamiltonian in convenition II"
showbonds ::usage = "Show the information of bond lengths and the classification of the bonds"


EndPackage[]
