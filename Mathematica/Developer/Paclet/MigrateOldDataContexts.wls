(* One-time source-data migration for the archived backend.  The original MX
   files were serialized with MagneticTB` coordinate symbols; the old backend
   must own those symbols directly so importing its data is context-isolated. *)

ClearAll["Global`*"];

projectRoot = DirectoryName[DirectoryName[DirectoryName[
  ExpandFileName[$InputFileName]
]]];
dataRoot = FileNameJoin[{
  projectRoot, "MagneticTB", "Old", "MagneticTBOld", "Data"
}];

sourceNames = {
  "MagneticTB`a", "MagneticTB`b", "MagneticTB`c",
  "MagneticTB`x", "MagneticTB`y", "MagneticTB`z",
  "MagneticTB`mx", "MagneticTB`my", "MagneticTB`mz",
  "MagneticTB`\[Alpha]", "MagneticTB`\[Beta]",
  "MagneticTB`\[Gamma]", "MagneticTB`\[Theta]"
};
targetNames = {
  "MagneticTBOld`aold", "MagneticTBOld`bold", "MagneticTBOld`cold",
  "MagneticTBOld`xold", "MagneticTBOld`yold", "MagneticTBOld`zold",
  "MagneticTBOld`mxold", "MagneticTBOld`myold", "MagneticTBOld`mzold",
  "MagneticTBOld`alphaold", "MagneticTBOld`betaold",
  "MagneticTBOld`gammaold", "MagneticTBOld`thetaold"
};
sourceSymbols = Symbol /@ sourceNames;
targetSymbols = Symbol /@ targetNames;
replacementRules = Thread[sourceSymbols -> targetSymbols];

newContextSymbolQ[expression_] := !FreeQ[
  expression,
  symbol_Symbol /; Context[Unevaluated[symbol]] === "MagneticTB`",
  Infinity
];

reports = Map[
  Function[file,
    path = FileNameJoin[{dataRoot, file}];
    data = Import[path];
    migrated = data /. replacementRules;
    If[newContextSymbolQ[migrated],
      Print["Migration did not remove every MagneticTB` symbol from ", path];
      Exit[2]
    ];
    Export[path, migrated, "MX"];
    <|"File" -> file, "Bytes" -> FileByteCount[path]|>
  ],
  {"wyckoffMSG.mx", "rod.mx", "layer.mx"}
];

Print[InputForm[<|"MigratedOldData" -> reports, "Passed" -> True|>]];
Exit[0];
