(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

magneticTBPackageDirectory = DirectoryName[DirectoryName[$InputFileName]];

MSGDATA = Import[
  FileNameJoin[{magneticTBPackageDirectory, "Data", "MSGData.mx"}]
];
wyckoffmsg = Import[
  FileNameJoin[{magneticTBPackageDirectory, "Data", "wyckoffMSG.mx"}]
];
{rgOgToSymbol, rgToBrav, BasicVectorsMRG, rodop, grayrod} = Import[
  FileNameJoin[{magneticTBPackageDirectory, "Data", "rod.mx"}]
];
{lgOgToSymbol, lgToBrav, BasicVectorsMLG, layerop, graylayer} = Import[
  FileNameJoin[{magneticTBPackageDirectory, "Data", "layer.mx"}]
];

MSGOP = MSGDATA["MSGOP"];
gray = MSGDATA["gray"];
typeI = MSGDATA["typeI"];
typeIII = MSGDATA["typeIII"];
typeIV = MSGDATA["typeIV"];
bnsdict = MSGDATA["bnsdict"];
ogdict = MSGDATA["ogdict"];
ognumdict = MSGDATA["ognumdict"];

End[]

EndPackage[]
