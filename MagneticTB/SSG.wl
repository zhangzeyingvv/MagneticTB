(* ::Package:: *)

BeginPackage["MagneticTB`"]
(*MagneticTB for spin space groups beta version 0.1*)


Begin["`Private`"]



findnSSG[c_, lbond_, lattplot_] := 
  Module[{atom, natom, nwyck, mksc, sizesc, sc, split}, 
   sizesc = {3, 3, 3};
   nwyck = Length[c];
   natom = Length /@ c;
   mksc = 
    Flatten[Table[# + {i, j, k}, {i, -sizesc[[1]], 
        sizesc[[1]]}, {j, -sizesc[[2]], 
        sizesc[[2]]}, {k, -sizesc[[3]], sizesc[[3]]}], 2] &;
   sc = Map[mksc, c, {2}];
   split[list_] := {N@#[[1, 1]], Length[#[[;; , 2]]], #[[;; , 2]]} & /@
      Split[list, First[#1] == First[#2] &];
   Table[
    split@Sort[
      Table[{EuclideanDistance[c[[wkp, wkpeqvi]] . lattplot, 
         sc[[scwkp, scwkpeqvi, i]] . lattplot], {c[[wkp, wkpeqvi]], 
         sc[[scwkp, scwkpeqvi, i]]}}, {i, 
        Length[sc[[scwkp, scwkpeqvi]]]}], 
      N@#1[[1]] < N@#2[[1]] &], {wkp, nwyck}, {wkpeqvi, 
     Length[c[[wkp]]]}, {scwkp, nwyck}, {scwkpeqvi, 
     Length[c[[scwkp]]]}]];


initSSG[inputDict_] := Module[
   {inputAssociation, x1, y1, z1,
    spinMatrix, spinMatrix2, directSum, timesSSG, identityElement,
    lattice, lattpar, wyckoffposition, symminformation, basisFunctions,
    lattplot, reclatt, reclattplot, volume, volumeplot, symmop, 
    Allsymmop, Gensymmop,
    atominfo, atompos, pos, magdirection, spinrot, rotorb, rotspin, 
    rule,
    ddim, repmats,
    base, pbase, nbases,
    coeff, opmatrix, sol,
    symopwyck, symmetryops,
    bondclassify, wcc
    },
   
   spinMatrix[op_] := 
    FullSimplify[
     Module[{\[Alpha], \[Beta], \[Gamma], nop}, 
      If[Det[op] == -1, nop = -op, nop = op];
      {\[Alpha], \[Beta], \[Gamma]} = -EulerAngles[nop];
      Transpose@{{E^(I \[Gamma]/2) Cos[\[Beta]/2] E^(I \[Alpha]/2), 
         E^(I \[Gamma]/2) Sin[\[Beta]/
            2] E^(-I \[Alpha]/2)}, {-E^(-I \[Gamma]/2) Sin[\[Beta]/
            2] E^(I \[Alpha]/2), 
         E^(-I \[Gamma]/2) Cos[\[Beta]/2] E^(-I \[Alpha]/2)}}]];
   
   spinMatrix2[m_] := 
    Module[{ang, axis, ovec, nn, nvec, rm, s, w, w1, wm, xx, yy, zz, 
      mat, axisAngle}, 
     If[FullSimplify[Det[m]] == -1, mat = -m, mat = m];
     axis = {mat[[3, 2]] - mat[[2, 3]], mat[[1, 3]] - mat[[3, 1]], 
       mat[[2, 1]] - mat[[1, 2]]};
     nn = Simplify[Norm[axis]];
     If[nn == 0, ang = \[Pi] Boole[Total[Diagonal[mat]] < 3];
      rm = (mat + IdentityMatrix[3])/2;
      axis = 
       Normalize[Extract[rm, Ordering[Max /@ Abs[rm], -1]]], {xx, yy, 
        zz} = Simplify[axis/nn];
      s = 2 UnitStep[zz] - 1; w = -1/(s + zz); w1 = xx yy w;
      ovec = {1 + s w xx xx, s w1, -s xx};
      nvec = {w1, s + w yy yy, -yy};
      wm = mat . ovec;
      ang = Arg[Simplify[wm . ovec + I wm . nvec]]];
     ExpToTrig@
      MatrixExp[-I FullSimplify@
         ang Sum[PauliMatrix[i] FullSimplify[(Normalize@axis)][[
            i]], {i, 3}]/2]];
   
   
   directSum = Fold[ArrayFlatten[{{#, 0}, {0, #2}}] &, #1, {##2}] &;
   timesSSG = 
    Association[
      "space" -> {First@#1["space"] . First@#2["space"], 
        Mod[First@#1["space"] . #2["space"][[2]] + #1["space"][[2]], 
         1]},
      "spin" -> {First@#1["spin"] . First@#2[["spin"]], 
        Mod[#1["spin"][[2]] + #2[["spin"]][[2]] , 2]}
      ] &;
   identityElement = <|"space" -> {{{1, 0, 0}, {0, 1, 0}, {0, 0, 
         1}}, {0, 0, 0}}, 
     "spin" -> {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 0}|>;
   
   
   inputAssociation = KeyMap[ToLowerCase, inputDict];
   lattice = 
    inputAssociation[
      "lattice"] /. {x_?InexactNumberQ :> 
       Rationalize@Round[x, .001]};
   lattpar = 
    inputAssociation[
      "lattpar"] /. {x_?InexactNumberQ :> 
       Rationalize@Round[x, .001]};
   wyckoffposition = 
    inputAssociation[
       "wyckoffposition"] /. {x_?InexactNumberQ :> 
        Rationalize@Round[x, .001]} /. {"x" -> x, "y" -> y, 
      "z" -> z};
   
   symminformation = 
    inputAssociation[
      "symminformation"] /. {x_?InexactNumberQ :> 
       Rationalize@Round[x, .001]};
   
   
   
   lattplot = lattice /. lattpar;
   volume = Cross[lattice[[1]], lattice[[2]]] . lattice[[3]];
   reclatt = 2  Pi/volume  {
      Cross[lattice[[2]], lattice[[3]]], 
      Cross[lattice[[3]], lattice[[1]]],
      Cross[lattice[[1]], lattice[[2]]]
      };
   reclattplot = reclatt /. lattpar;
   volumeplot = volume /. lattpar;
   
   symmop = Values[symminformation];
   Allsymmop = GenerateGroup[symmop, identityElement, timesSSG];
   
   Gensymmop = getGenerator[Allsymmop, identityElement, timesSSG];
   
   atompos = Table[
     Table[
      pos = 
       Mod[#1 . wyck["FractionalCoordinates"] + #2, 1] & @@ (sym[
          "space"]);
      magdirection = (If[#2 == 1, -1, 1] Det[#1] #1 . 
            wyck["MagnetizationDirection"]) & @@ (sym["spin"]);
      
      
      {pos, magdirection}
      , {sym, Allsymmop}],
     {wyck, wyckoffposition}];
   atompos = DeleteDuplicates /@ atompos;
   
   
   
   If[Max[Length /@ Values[GroupBy[Flatten[atompos, 1], #[[1]] &]]] > 
     1,
    Print[
     "Wyckoff position(s) are not correct! Because one positios have \
more than one atoms.\n", atompos]; Abort[]];
   
   repmats = Association[];
   Do[repmats[wyck["FractionalCoordinates"]] = Association[], {wyck, 
     wyckoffposition}];
   Do[
    base = wyck["BasisFunction"];
    nbases = Length[base];
    
    Do[
     
     Which[
      ListQ[base[[1]]] == False,
      
      
      rule = 
       FullSimplify@
        MapThread[
         Rule, {{x1, y1, z1}, 
          FullSimplify@(Inverse[
              sym["space"][[1]]] . ({x1, y1, z1} . 
               Inverse[lattplot]) . lattplot)}];
      
      pbase = 
       If[sym["spin"][[2]] == 1, 
        FullSimplify@
         ComplexExpand@Conjugate[(base) /. rule], (base) /. rule];
      
      opmatrix = Table[coeff[i, j], {i, nbases}, {j, nbases}];
      sol = First@SolveAlways[pbase == opmatrix . base, {x, y, z}];
      opmatrix = opmatrix /. sol;
      opmatrix = Transpose[opmatrix];
      
      repmats[wyck["FractionalCoordinates"]][sym] = opmatrix,
      
      ListQ[base[[1]]] == True,
      
      rule = 
       FullSimplify@
        MapThread[
         Rule, {{x, y, z}, 
          FullSimplify@(Inverse[
              sym["space"][[1]]] . ({x, y, z} . Inverse[lattplot]) . 
             lattplot)}];
      If[sym["spin"][[2]] == 1,
       pbase = 
        ComplexExpand[I PauliMatrix[2] . Conjugate[#]] & /@ base;
       pbase = (pbase) /. rule;
       pbase = 
        spinMatrix2[
            Transpose[lattplot] . sym["spin"][[1]] . 
             Inverse@Transpose@lattplot] . # & /@ pbase,
       pbase = (base) /. rule;
       pbase = 
        spinMatrix2[
            Transpose[lattplot] . sym["spin"][[1]] . 
             Inverse@Transpose@lattplot] . # & /@ pbase
       ];
      pbase = Simplify@pbase;
      opmatrix = Table[coeff[i, j], {i, nbases}, {j, nbases}];
      sol = 
       Flatten[Table[
         SolveAlways[(pbase[[j]]) == 
           Total@Table[opmatrix[[j, k]] base[[k]], {k, nbases}], {x, 
           y, z}], {j, nbases}]];
      opmatrix = FullSimplify[opmatrix /. sol];
      opmatrix = Transpose[opmatrix];
      
      repmats[wyck["FractionalCoordinates"]][sym] = opmatrix
      
      ]
     , {sym, Allsymmop}];
    , {wyck, wyckoffposition}];
   (*Print[repmats];*)
   symopwyck[wyck_] := 
    Module[{inputatom, pos1, wyckn, atom1, atom2, spaceop},
     pos1 = 0;
   (*  Print[wyck,repmats];*)
     wyckn = Length[wyck];
     inputatom = wyck[[1, 1]];
     Table[
      pos1 += 1;
      spaceop = sym["space"];
      spaceop = {Inverse[#[[1]]], -Inverse[#[[1]]] . #[[2]]} &@spaceop;
      ArrayFlatten[Table[
        atom1 = Mod[#1 . wyck[[i, 1]] + #2, 1] & @@ spaceop;
        atom2 = wyck[[j, 1]];
        (*Print[(*{atom1,wyck[[j, 1]]},*){wyck[[j, 1]],atom2,atom1==atom2}];*)
        If[atom1 == atom2, repmats[inputatom][sym], 0]
        , {i, wyckn}, {j, wyckn}]]
      , {sym, Allsymmop}]
     ];
   symmetryops = 
    Association[
     Thread[Allsymmop -> 
       Table[directSum @@ i, {i, Transpose[symopwyck /@ atompos]}]]];
   bondclassify = 
    Split[SortBy[
      Flatten[findnSSG[atompos[[;; , ;; , 1]], 1, lattplot], 
       4], #[[1]] &], #1[[1]] == #2[[1]] &];
   bondclassify = 
    Table[{#[[1, 1]], Total@#[[;; , 2]], Flatten[#[[;; , 3]], 1]} & /@
       Values[GroupBy[neigh, #[[3, 1, 1]] &]], {neigh, 
      bondclassify}];
   wcc = 
    Module[{natomperwyck, nbasesperwyck}, 
     natomperwyck = Length /@ atompos;
     nbasesperwyck = Length /@ wyckoffposition;
     Flatten[
      Table[Table[
        Table[atompos[[i, j, 1]], nbasesperwyck[[i]]], {j, 
         natomperwyck[[i]]}], {i, Length[natomperwyck]}], 2]];
   
   <|"atompos" -> atompos, "symmetryops" -> symmetryops, "wcc" -> wcc,
     "bondclassify" -> bondclassify, "input" -> inputAssociation, 
    "lattplot" -> lattplot, "reclatt" -> reclatt, 
    "reclattplot" -> reclattplot, "volume" -> {volume, volumeplot}, 
    "Gens" -> Gensymmop, "pointrepmats" -> repmats, 
    "pointops" -> Values /@ Values[repmats]|>
   
   ];


unsymhamSSG[n_, bondclassify_, atompos_, pointops_] := Module[
   {bondends, bandind, hop, bonds, toband, nbh, hami, tobond},
   tobond = #[[2]] - #[[1]] &;
   bondends = Flatten[bondclassify[[n]][[;; , 3]], 1]; 
   bandind = 
    Map[FirstPosition[
       Map[Flatten[
          Table[# + {i, j, k}, {i, -3, 3}, {j, -3, 3}, {k, -3, 3}], 
          2] &, atompos[[;; , ;; , 1]], {2}], #] &, bondends, {2}]; 
   Do[hop[bandind[[i]]] = 
      Array[ToExpression[
         "tr" <> ToString[i] <> ToString[#1] <> ToString[#2] <> 
          "+I ti" <> ToString[i] <> ToString[#1] <> 
          ToString[#2]] &, {Length[
         pointops[[bandind[[i, 1, 1]]]][[1]]], 
        Length[pointops[[bandind[[i, 2, 1]]]][[1]]]}];, {i, 
     Length[bandind]}]; bonds = tobond /@ bondends;
   
   toband = 
    Total[((Length /@ atompos) (*(Length/@repMatrix[[;;,
         1]])*))[[;; #[[1]] - 1]]] + #[[2]] &; 
   nbh = (Length /@ atompos);
   Do[
    h[toband[{i, j}], toband[{k, l}]] =
     Table[
      0, {Length[pointops[[i]][[1]]]}, {Length[pointops[[k]][[1]]]}]
    , {i, Length[nbh]}, {j, nbh[[i]]}, {k, Length[nbh]}, {l, 
     nbh[[k]]}]; 
   Do[h[toband[bandind[[i, 1]]], toband[bandind[[i, 2]]]] += 
      Exp[I ({kx, ky, kz}) . (bonds[[i]])] hop[bandind[[i]]];, {i, 
     Length[bondends]}];
   hami = 
    ArrayFlatten[Table[h[i, j], {i, Total@nbh}, {j, Total@nbh}]]];


Options[symhamSSG]={"symmetryset"->All};
(*SetOptions[symham,symmetryset\[Rule]All];*)
symhamSSG[n_,symmetry_,OptionsPattern[]]:=Module[
{conjh,h,phpmh,recR,para,sset,opset,symmetryinfo,unham,
ops},
ops=KeyTake[symmetry["symmetryops"],symmetry["Gens"]];
sset=OptionValue["symmetryset"];
If[sset===All,opset=Range[Length[ops]],opset=sset,opset=Range[Length[ops]]];
symmetryinfo=Normal[ops];
unham=unsymhamSSG[n,symmetry["bondclassify"],symmetry["atompos"],symmetry["pointops"]];
conjh=(ComplexExpand@ConjugateTranspose@unham)-unham;If[n==1,h=unham/. ToRules@Reduce[conjh==0],h=unham/. ToRules@Reduce[Flatten[Map[Total/@(Values[GroupBy[Cases[#,x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,TrigToExp[conjh],{2}]]==0];h=h/. ToRules@Reduce[Flatten[Map[Total/@(Values[GroupBy[Cases[{#},x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,TrigToExp[conjh],{2}]]==0]];phpmh=Table[
recR=Transpose[symmetryinfo[[i,1]]["space"][[1]]];Which[symmetryinfo[[i,1]]["spin"][[2]]==0,(Inverse[symmetryinfo[[i,2]]] . h . symmetryinfo[[i,2]])-(h/. Thread[{kx,ky,kz}->recR . {kx,ky,kz}]),symmetryinfo[[i,1]]["spin"][[2]]==1,(Expand@TrigToExp@(Inverse[symmetryinfo[[i,2]]] . (h) . symmetryinfo[[i,2]]))-Expand@TrigToExp@(ComplexExpand@Conjugate[h/. Thread[{kx,ky,kz}->recR . {-kx,-ky,-kz}]])],{i,opset}];

phpmh=Expand[phpmh];If[n==1,(*Print[First@Solve[phpmh\[Equal]0]];*)h=h/. ToRules@Reduce[phpmh==0],phpmh=Join[Map[Total@Cases[#,xx_ E^yy_:>xx E^Simplify@yy]&,phpmh,{3}],Map[Total@Cases[{#},xx_ E^yy_:>xx E^Simplify@yy]&,phpmh,{3}]];h=h/. FullSimplify@ToRules@Reduce[DeleteDuplicates@Flatten[Map[Total/@(Values[GroupBy[Cases[#,x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,phpmh,{3}]]==0];h=h/. FullSimplify@ToRules@Reduce[DeleteDuplicates@Flatten[Map[Total/@(Values[GroupBy[Cases[{#},x_ E^y_:>{x,Expand[y]}],Last]][[;;,;;,1]])&,phpmh,{3}]]==0];];para=Thread[Which[n==1,#->Table[ToExpression["e"<>ToString[i]],{i,Length[#]}]&@Variables[h],n==2,#->Table[ToExpression["t"<>ToString[i]],{i,Length[#]}]&@Variables[h],n==3,#->Table[ToExpression["r"<>ToString[i]],{i,Length[#]}]&@Variables[h],n==4,#->Table[ToExpression["s"<>ToString[i]],{i,Length[#]}]&@Variables[h],True,#->Table[ToExpression["p"<>ToString[n]<>"n"<>ToString[i]],{i,Length[#]}]&@Variables[h]]];h=h/. para;Print["params:",Variables[h]];
h];



End[]
EndPackage[]
