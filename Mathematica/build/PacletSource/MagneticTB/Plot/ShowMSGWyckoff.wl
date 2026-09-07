(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

bnsdictreverse = Association[Reverse /@ Normal[bnsdict]];

showMSGWyckoff[msg_] := Module[{nu, nw, wyck},
  If[NumberQ[msg], nu = bnsdictreverse[msg], nu = msg];
  Print["MSG:", nu];
  wyck = wyckoffmsg[nu];
  wyck = Reverse[wyck];
  nw = Length[wyck];
  wyck = MapAt[If[NumberQ[#], Mod[#, 1], #] &, wyck, {;;, 4, ;;, 1, ;;}];
  wyck = MapAt[Grid[#] &, wyck, {;;, 4}];
  wyck = Delete[#, 2] & /@ wyck;
  PrependTo[
    wyck,
    {"Multiplicity", "Wyckoff Letter", "Atomic Positions & Magnetization directions"}
  ];
  Grid[wyck, Frame -> All, Background -> {{LightYellow, LightBlue}, None}]
];

End[]
EndPackage[]
