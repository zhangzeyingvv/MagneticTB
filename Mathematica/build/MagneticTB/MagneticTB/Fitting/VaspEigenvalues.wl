(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[vaspNumberStringQ, parseVaspNumber, parseVaspNumericLine];

vaspNumberStringQ[value_String] := StringMatchQ[
  value,
  RegularExpression[
    "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[Ee][+-]?[0-9]+)?$"
  ]
];

parseVaspNumber[value_String] /; vaspNumberStringQ[value] :=
  ToExpression[StringReplace[value, {"E" -> "*^", "e" -> "*^"}]];

parseVaspNumericLine[line_String] := Module[{tokens},
  tokens = StringSplit[StringTrim[line]];
  If[tokens === {} || !And @@ (vaspNumberStringQ /@ tokens),
    Return[$Failed]
  ];
  parseVaspNumber /@ tokens
];

vaspEig::file = "The EIGENVAL file could not be read: `1`.";
vaspEig::format =
  "The file does not have a complete numeric VASP EIGENVAL structure.";
vaspEig::spin =
  "spin must select an available eigenvalue column; received `1`.";
vaspEig::bands =
  "The requested 1-based band interval `1` through `2` is outside 1 through `3`.";

vaspEig[filename_, efermi_, spin_, startband_, endband_] := Module[
  {
    lines,
    spinHeader,
    nspin,
    header,
    nk,
    nband,
    dataLines,
    blocks,
    parsedBlocks,
    availableColumns
  },
  If[!StringQ[filename] || !FileExistsQ[filename],
    Message[vaspEig::file, filename];
    Return[$Failed]
  ];
  lines = Quiet@Check[Import[filename, "Lines"], $Failed];
  If[lines === $Failed || !ListQ[lines] || Length[lines] < 6,
    Message[vaspEig::file, filename];
    Return[$Failed]
  ];
  spinHeader = parseVaspNumericLine[lines[[1]]];
  If[
    spinHeader === $Failed || spinHeader === {} ||
      !IntegerQ[Last[spinHeader]] ||
      !MemberQ[{1, 2}, Last[spinHeader]],
    Message[vaspEig::format];
    Return[$Failed]
  ];
  nspin = Last[spinHeader];
  header = parseVaspNumericLine[lines[[6]]];
  If[
    header === $Failed || Length[header] < 2 ||
      !VectorQ[header[[-2 ;; -1]], IntegerQ],
    Message[vaspEig::format];
    Return[$Failed]
  ];
  {nk, nband} = header[[-2 ;; -1]];
  If[
    nk <= 0 || nband <= 0 || !NumericQ[efermi] ||
      !TrueQ[Chop[Im[N[efermi]]] === 0],
    Message[vaspEig::format];
    Return[$Failed]
  ];
  If[!IntegerQ[spin] || !Between[spin, {1, nspin}],
    Message[vaspEig::spin, spin];
    Return[$Failed]
  ];
  If[
    !IntegerQ[startband] || !IntegerQ[endband] ||
      startband < 1 || endband < startband || endband > nband,
    Message[vaspEig::bands, startband, endband, nband];
    Return[$Failed]
  ];

  dataLines = Select[
    Drop[lines, 6],
    StringTrim[#] =!= "" &
  ];
  If[Length[dataLines] =!= nk (nband + 1),
    Message[vaspEig::format];
    Return[$Failed]
  ];
  blocks = Partition[dataLines, nband + 1];
  parsedBlocks = Map[parseVaspNumericLine, blocks, {2}];
  If[!FreeQ[parsedBlocks, $Failed],
    Message[vaspEig::format];
    Return[$Failed]
  ];
  availableColumns = Min[Length /@ Flatten[parsedBlocks[[All, 2 ;;]], 1]];
  If[spin + 1 > availableColumns,
    Message[vaspEig::spin, spin];
    Return[$Failed]
  ];

  Map[
    Function[block,
      {
        block[[1, 1 ;; 3]],
        block[[2 ;;, spin + 1]][[startband ;; endband]] - efermi
      }
    ],
    parsedBlocks
  ]
];

End[]

EndPackage[]
