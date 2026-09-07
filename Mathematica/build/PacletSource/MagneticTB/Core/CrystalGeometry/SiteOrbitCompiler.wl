(* ::Package:: *)

BeginPackage["MagneticTB`"]

Begin["`Private`"]

ClearAll[nativeAtomPositions, nativeSpinSpaceGroupAtomPositions];

nativeAtomPositions[wyckoffSeeds_List, symmetryInformation_List] :=
  DeleteDuplicates /@ Map[
    Function[seed,
      Map[
        Function[operation,
          {
            Mod[operation[[2]].seed[[1]] + operation[[3]], 1],
            Which[
              operation[[4]] === "F",
                Det[operation[[2]]] operation[[2]].seed[[2]],
              operation[[4]] === "T",
                -Det[operation[[2]]] operation[[2]].seed[[2]],
              True,
                Return[$Failed, Function]
            ]
          }
        ],
        symmetryInformation
      ]
    ],
    wyckoffSeeds
  ];

nativeSpinSpaceGroupAtomPositions[
    wyckoffSeeds_List,
    elements_List
  ] := Module[{rawPositions, uniquePositions, positionGroups},
  rawPositions = Map[
    Function[seed,
      Map[
        Function[element,
          {
            Mod[element["space"][[1]].seed[[1]] +
              element["space"][[2]], 1],
            If[element["spin"][[2]] === 1, -1, 1] *
              element["spin"][[1]].seed[[2]]
          }
        ],
        elements
      ]
    ],
    wyckoffSeeds
  ];
  uniquePositions = DeleteDuplicates /@ rawPositions;
  positionGroups = GroupBy[#, First] & /@ uniquePositions;
  If[AnyTrue[positionGroups, AnyTrue[Values[#], Length[#] > 1 &] &],
    Return[$Failed]
  ];
  uniquePositions
];

End[]

EndPackage[]
