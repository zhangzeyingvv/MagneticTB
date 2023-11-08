

fittingTB[h_, eigdata_, krange_, initparms_] := Module[
  {params, mparams, klist, hp, modelhams, nk, objparams, res, lsq},
  params = Complement[Variables[h], {kx, ky, kz}];
  klist = eigdata[[;; , 1]];
  nk = Length[krange];
  mparams = params;
  hp = h /. Thread[params -> mparams];
  modelhams = 
   Table[hp /. {kx -> k[[1]], ky -> k[[2]], kz -> k[[3]]}, {k, 
     2 Pi klist[[krange]]}];
  res = FindMinimum[
    Total@Table[
      EuclideanDistance[Sort@Eigenvalues[modelhams[[i]]], 
       eigdata[[;; , 2]][[krange[[i]]]]], {i, nk}], 
    Transpose[{mparams, params /. initparms}]];
  objparams = res[[2]] /. Thread[mparams -> params];
  
  lsq = {Total@
      Table[EuclideanDistance[Sort@Eigenvalues[modelhams[[i]]], 
        eigdata[[;; , 2]][[krange[[i]]]]], {i, nk}] /. (initparms /. 
       Thread[params -> mparams]), res[[1]]};
  Association["FittedParams" -> objparams, 
   "LSQForInitParams" -> lsq[[1]], "LSQForFittedParams" -> lsq[[2]]]
  ]
