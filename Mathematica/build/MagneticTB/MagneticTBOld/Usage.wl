(* ::Package:: *)

(* Forward declarations to prevent shadowing messages *)


BeginPackage["MagneticTBOld`"]


(* ::Subsubsection:: *)
(*core*)


(*set input data*)
initold ::usage = "Init the program"
unsymhamold ::usage = "Get the unsymmetried Hamiltonian"
symhamold ::usage = "Get the symmetried Hamiltonian in convenition I"
basisold::usage = "basis"
basisdictold::usage = "basis"
kxold ::usage = ""
kyold ::usage = ""
kzold ::usage = ""
aold ::usage = ""
bold ::usage = ""
cold ::usage = ""
alphaold ::usage = ""
betaold ::usage = ""
gammaold ::usage = ""
xold ::usage = ""
yold ::usage = ""
zold ::usage = ""
mxold ::usage = ""
myold ::usage = ""
mzold ::usage = ""
thetaold ::usage = ""


(*output data*)
opsold  ::usage ="ops"
braLattold ::usage = ""
grayold ::usage = ""
grayrodold::usage = ""
graylayerold::usage = ""
typeIold ::usage = ""
typeIIIold ::usage = ""
typeIVold ::usage = ""
bnsdictold ::usage = ""
ogdictold ::usage = ""
ognumdictold ::usage = ""
bondclassifyold  ::usage = ""
symminfoold ::usage = ""
symmetryopsold  ::usage = ""
symmetryopsIIold   ::usage = ""
atomposold ::usage=""
reclattold ::usage=""
lattold ::usage=""
wccold ::usage = "Wannier center"
symmcompileold ::usage = ""
compactFormold ::usage = ""
latticeold::usage = "Options for init"
wyckoffpositionold::usage = "Options for init"
symminformationold::usage = "Options for init"
basisFunctionsold::usage = "Options for init"
symmetrysetold::usage = "Options for symham"
lattparold::usage = "Options for init"
symmetrizationHRInitold::usage = "Options for init"
debugQold::usage = "Options for initold"
plotRangeold::usage = "Options for bandplotold"


(* ::Subsubsection:: *)
(*Plot*)


(*Plot*)
bandManipulateold ::usage=""

bandplotold ::usage=""
banddataold ::usage=""


(*Functions*)
readMsgDataold ::usage = "Read the magnetic space group data"
pointMatrixold  ::usage = ""
tokpathvaspold ::usage = ""


realham2dold ::usage=""



(* ::Subsubsection:: *)
(*IO*)


(*
Interface
*)
pythtbIO2dold ::usage = "init"
plotslabold ::usage = ""
hopold  ::usage =""
readHRold ::usage =""


(* ::Subsubsection:: *)
(*Wilson*)


(*Wilson loop (testing)*)
wLoopold ::usage = "init"
berryphold ::usage = ""
z2pathold ::usage =""
wilsonLoopold ::usage =""
plotWilsonLoopold ::usage =""



(* ::Subsubsection:: *)
(*Corep*)


getMSGElemFromMSGCorepold::usage = ""
getTBBandCorepold::usage = ""


(* ::Subsubsection:: *)
(*Utilities*)


symhamIIold ::usage = "Get the symmetried Hamiltonian in convenition II"
showbondsold ::usage = "Show the information of bond lengths and the classification of the bonds"
GenerateGroupold::usage =""
getGeneratorold::usage =""
showMSGWyckoffold::usage=""


msgopold ::usage =""
mlgopold::usage=""
mrgopold::usage=""
texOutputold ::usage =""


(* ::Subsubsection:: *)
(*Fitting*)


bandManipulateEigold ::usage=""
vaspEigold ::usage=""
compareBandold ::usage=""


(* ::Subsubsection:: *)
(*SSG*)


initSSGold ::usage=""
symhamSSGold ::usage=""


EndPackage[]
