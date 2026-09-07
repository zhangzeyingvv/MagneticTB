(* Build the official MagneticTB release layout without changing the source tree.

   The repository keeps the archived backend at MagneticTB/Old/MagneticTBOld,
   while a paclet requires non-overlapping Kernel extension roots. This script
   stages the new and old backends as sibling roots before calling PacletBuild.
*)

ClearAll["Global`*"];

(* The standalone macOS WolframScript launcher can expose a UTF-8 working
   directory as Latin-1 mojibake.  Do not trust either Directory[] or
   $InputFileName blindly: repair both forms, walk upwards, and accept a root
   only when the repository markers are actually present. *)
ClearAll[
  repairUTF8Mojibake,
  pathVariants,
  ancestorDirectories,
  magneticTBProjectRootQ,
  locateProjectRoot
];

repairUTF8Mojibake[path_String] := Module[{codes, repaired},
  codes = ToCharacterCode[path, "ISO8859-1"];
  (* A correctly decoded Unicode path contains characters outside Latin-1.
     In that case there is nothing to repair, and passing None entries to
     FromCharacterCode would itself emit a misleading build warning. *)
  If[!VectorQ[codes, IntegerQ], Return[path]];
  repaired = Check[FromCharacterCode[codes, "UTF8"], $Failed];
  If[StringQ[repaired], repaired, path]
];

pathVariants[path_String] := DeleteDuplicates@Select[
  DeleteCases[
    {
      path,
      repairUTF8Mojibake[path],
      Check[ExpandFileName[path], $Failed],
      Check[repairUTF8Mojibake[ExpandFileName[path]], $Failed],
      Check[ExpandFileName[repairUTF8Mojibake[path]], $Failed]
    },
    $Failed
  ],
  StringQ
];

ancestorDirectories[path_String] := DeleteDuplicates@FixedPointList[
  DirectoryName,
  path
];

magneticTBProjectRootQ[path_String] :=
  DirectoryQ[path] &&
    FileExistsQ[FileNameJoin[{path, "PacletInfo.wl"}]] &&
    DirectoryQ[FileNameJoin[{path, "MagneticTB"}]] &&
    DirectoryQ[FileNameJoin[{path, "Developer", "Documentation"}]];

locateProjectRoot[] := Module[
  {rawCandidates, variants, startingDirectories, candidates},
  rawCandidates = Select[
    Join[
      {$InputFileName, Directory[]},
      If[ListQ[$ScriptCommandLine], $ScriptCommandLine, {}]
    ],
    StringQ
  ];
  variants = DeleteDuplicates@Flatten[pathVariants /@ rawCandidates];
  startingDirectories = DeleteDuplicates@Flatten[
    {variants, DirectoryName /@ variants}
  ];
  candidates = DeleteDuplicates@Flatten[
    ancestorDirectories /@ startingDirectories
  ];
  SelectFirst[candidates, magneticTBProjectRootQ, Missing["NotFound"]]
];

projectRoot = locateProjectRoot[];
If[MissingQ[projectRoot],
  Print[
    "Unable to locate the MagneticTB project root. ",
    "Run this script from the repository or set a valid script path."
  ];
  Exit[2]
];

ClearAll[
  releaseLog,
  releaseElapsedSeconds
];

releaseLog[parts___] := Print[
  "[",
  DateString[{
    "Year", "-", "Month", "-", "Day", " ",
    "Hour", ":", "Minute", ":", "Second"
  }],
  "] ",
  parts
];

releaseElapsedSeconds[start_] :=
  Round[10 N[AbsoluteTime[] - start]]/10.;

releaseStartedAt = AbsoluteTime[];
releaseLog["MagneticTB full release build started."];

newSource = FileNameJoin[{projectRoot, "MagneticTB"}];
oldSource = FileNameJoin[{newSource, "Old", "MagneticTBOld"}];
documentationSource = FileNameJoin[{projectRoot, "Documentation"}];
chineseDocumentationSource = FileNameJoin[{
  documentationSource, "ChineseSimplified"
}];
pacletInfoSource = FileNameJoin[{projectRoot, "PacletInfo.wl"}];
licenseSource = FileNameJoin[{DirectoryName[projectRoot], "LICENSE"}];
buildRoot = FileNameJoin[{projectRoot, "build"}];
stageRoot = FileNameJoin[{buildRoot, "PacletSource"}];
chineseStageRoot = FileNameJoin[{buildRoot, "ChineseDocumentationSource"}];
chineseBuildRoot = FileNameJoin[{buildRoot, "ChineseDocumentationBuild"}];
builtPacletRoot = FileNameJoin[{buildRoot, "MagneticTB"}];

requiredPaths = {
  newSource,
  oldSource,
  documentationSource,
  chineseDocumentationSource,
  pacletInfoSource,
  licenseSource
};

If[!And @@ (FileExistsQ /@ requiredPaths),
  Print["Missing release input: ", Select[requiredPaths, Not@*FileExistsQ]];
  Exit[3]
];

pacletObject = Get[pacletInfoSource];
If[!MatchQ[pacletObject, _PacletObject],
  Print["PacletInfo.wl did not evaluate to a PacletObject."];
  Exit[11]
];
pacletVersion = pacletObject["Version"];

expectedEnglishPageCount = Length@FileNames[
  "*.nb",
  FileNameJoin[{documentationSource, "English"}],
  Infinity
];
expectedChinesePageCount = Length@FileNames[
  "*.nb",
  chineseDocumentationSource,
  Infinity
];
If[
  expectedEnglishPageCount === 0 ||
    expectedEnglishPageCount =!= expectedChinesePageCount,
  Print[
    "The bilingual authoring trees are incomplete or out of sync: ",
    {expectedEnglishPageCount, expectedChinesePageCount}
  ];
  Exit[10]
];

stageStartedAt = AbsoluteTime[];
releaseLog["[1/6] START staging source, both backends, and documentation."];
Scan[
  If[DirectoryQ[#], DeleteDirectory[#, DeleteContents -> True]] &,
  {
    stageRoot,
    builtPacletRoot,
    FileNameJoin[{buildRoot, "incremental"}]
  }
];
CreateDirectory[stageRoot, CreateIntermediateDirectories -> True];

CopyFile[pacletInfoSource, FileNameJoin[{stageRoot, "PacletInfo.wl"}]];
CopyDirectory[newSource, FileNameJoin[{stageRoot, "MagneticTB"}]];

(* The archived source remains nested in the repository, but not inside the
   new Kernel extension in the release layout. *)
DeleteDirectory[
  FileNameJoin[{stageRoot, "MagneticTB", "Old"}],
  DeleteContents -> True
];
CopyDirectory[oldSource, FileNameJoin[{stageRoot, "MagneticTBOld"}]];
CopyDirectory[documentationSource, FileNameJoin[{stageRoot, "Documentation"}]];
DeleteDirectory[
  FileNameJoin[{stageRoot, "Documentation", "ChineseSimplified"}],
  DeleteContents -> True
];

Scan[DeleteFile, FileNames[".DS_Store", stageRoot, Infinity]];
releaseLog[
  "[1/6] DONE staging source (",
  releaseElapsedSeconds[stageStartedAt],
  " s)."
];

Needs["PacletTools`"];

(* Wolfram 15 initializes the PNG/JPEG exporters by looking for the optional
   ImageMetadataTools paclet.  When that paclet is unavailable the raster
   exporters still work, but RegisterFormat emits an internal-error message
   on their first use.  Prime the two formats before the strict build-message
   gate so that missing optional metadata support cannot invalidate an
   otherwise correct documentation build. *)
Quiet[
  {
    ExportByteArray[Graphics[{Black, Point[{0, 0}]}], "PNG"],
    ExportByteArray[Graphics[{Black, Point[{0, 0}]}], "JPEG"]
  },
  {
    ImportExport`RegisterFormat::interr,
    CloudObject::cloudnf
  }
];

(* Mathematica 14.0's documentation skeletonizer tries to query
   WolframLanguageData for third-party symbol pages by replacing their names
   with the sentinel "notAnEntity".  That query emits a cascade of
   KeyTake/Lookup messages and contributes no useful ranking data to a custom
   paclet.  Disable only that system-symbol enrichment step for this build;
   notebook skeletonization and the local documentation search indexes remain
   enabled. *)
Needs["DocumentationSearch`Skeletonizer`"];
SetOptions[
  DocumentationSearch`Skeletonizer`Skeletonize,
  "UseWolframLanguageData" -> False
];
SetOptions[
  DocumentationSearch`Skeletonizer`Private`skeletonizeFiles,
  "UseWolframLanguageData" -> False
];

stageStartedAt = AbsoluteTime[];
releaseLog[
  "[2/6] START PacletBuild for the package and English documentation; ",
  "the official builder may be quiet until this stage finishes."
];
{result, buildMessages} = Block[{$MessageList = {}},
  {
    PacletTools`PacletBuild[stageRoot, buildRoot],
    $MessageList
  }
];
If[buildMessages =!= {},
  Print["Unexpected messages during the main paclet build: ", buildMessages];
  Exit[8]
];
If[!MatchQ[result, _Success],
  Print[InputForm[result]];
  Exit[1]
];
releaseLog[
  "[2/6] DONE package and English documentation (",
  releaseElapsedSeconds[stageStartedAt],
  " s)."
];

(* PacletDocumentationBuild 14.0 builds one language tree at a time and uses
   Documentation/English as its authoring input. Build the canonical Chinese
   pages through an isolated documentation-only stage, then place the built
   result under the standard ChineseSimplified language directory. *)
Scan[
  If[DirectoryQ[#], DeleteDirectory[#, DeleteContents -> True]] &,
  {chineseStageRoot, chineseBuildRoot}
];
CreateDirectory[
  FileNameJoin[{chineseStageRoot, "Documentation"}],
  CreateIntermediateDirectories -> True
];
CopyDirectory[
  chineseDocumentationSource,
  FileNameJoin[{chineseStageRoot, "Documentation", "English"}]
];
Put[
  PacletObject[<|
    "Name" -> "MagneticTB",
    "Version" -> pacletVersion,
    "WolframVersion" -> "12.1+",
    "Extensions" -> {
      {"Documentation", "LinkBase" -> "MagneticTB"}
    }
  |>],
  FileNameJoin[{chineseStageRoot, "PacletInfo.wl"}]
];

stageStartedAt = AbsoluteTime[];
releaseLog[
  "[3/6] START ChineseSimplified documentation build; ",
  "the official builder may be quiet until this stage finishes."
];
{chineseResult, chineseBuildMessages} = Block[{$MessageList = {}},
  {
    PacletTools`PacletDocumentationBuild[
      chineseStageRoot,
      chineseBuildRoot
    ],
    $MessageList
  }
];
If[chineseBuildMessages =!= {},
  Print[
    "Unexpected messages during the Chinese documentation build: ",
    chineseBuildMessages
  ];
  Exit[9]
];
If[!MatchQ[chineseResult, _Success],
  Print[InputForm[chineseResult]];
  Exit[4]
];
releaseLog[
  "[3/6] DONE ChineseSimplified documentation (",
  releaseElapsedSeconds[stageStartedAt],
  " s)."
];

stageStartedAt = AbsoluteTime[];
releaseLog["[4/6] START installing and checking the bilingual page trees."];
builtChineseEnglish = FileNameJoin[{
  chineseBuildRoot, "MagneticTB", "Documentation", "English"
}];
builtChineseTarget = FileNameJoin[{
  builtPacletRoot, "Documentation", "ChineseSimplified"
}];
If[!DirectoryQ[builtChineseEnglish],
  Print["Chinese documentation build directory was not created."];
  Exit[5]
];
If[DirectoryQ[builtChineseTarget],
  DeleteDirectory[builtChineseTarget, DeleteContents -> True]
];
CopyDirectory[builtChineseEnglish, builtChineseTarget];

englishPageCount = Length@FileNames[
  "*.nb",
  FileNameJoin[{builtPacletRoot, "Documentation", "English"}],
  Infinity
];
chinesePageCount = Length@FileNames["*.nb", builtChineseTarget, Infinity];
If[
  englishPageCount =!= expectedEnglishPageCount ||
    chinesePageCount =!= expectedChinesePageCount,
  Print["Unexpected documentation page counts: ", {
    englishPageCount, chinesePageCount
  }];
  Exit[6]
];
releaseLog[
  "[4/6] DONE bilingual page trees: English=",
  englishPageCount,
  ", ChineseSimplified=",
  chinesePageCount,
  " (",
  releaseElapsedSeconds[stageStartedAt],
  " s)."
];

(* PacletBuild writes PacletManifest.wl before the separately built Chinese
   documentation tree is installed.  Rebuild the manifest from the final
   staged paclet so every shipped file is covered.  PacletInfo.wl follows
   PacletTools' own convention and is excluded; PacletManifest.wl cannot hash
   itself. *)
stageStartedAt = AbsoluteTime[];
releaseLog["[5/6] START rebuilding and verifying the final manifest."];
CopyFile[licenseSource, FileNameJoin[{builtPacletRoot, "LICENSE"}]];
manifestExcludedFiles = {"PacletInfo.wl", "PacletManifest.wl"};
manifestRelativeFiles = Sort@Map[
  FileNameDrop[#, FileNameDepth[builtPacletRoot]] &,
  Select[
    FileNames["*", builtPacletRoot, Infinity],
    FileType[#] === File &
  ]
];
manifestRelativeFiles = Complement[
  manifestRelativeFiles,
  manifestExcludedFiles
];
completeManifest = <|
  "Files" -> Map[
    Function[path,
      <|
        "File" -> File[path],
        "Hash" -> FileHash[
          FileNameJoin[{builtPacletRoot, path}],
          "SHA256",
          "HexString"
        ]
      |>
    ],
    manifestRelativeFiles
  ]
|>;
Export[
  FileNameJoin[{builtPacletRoot, "PacletManifest.wl"}],
  CodeFormatter`CodeFormat[
    ToString[InputForm[completeManifest]],
    CodeFormatter`Airiness -> 0.95
  ],
  "String"
];

writtenManifest = Check[
  Get[FileNameJoin[{builtPacletRoot, "PacletManifest.wl"}]],
  $Failed
];
writtenManifestRecords = If[
  AssociationQ[writtenManifest],
  Lookup[writtenManifest, "Files", $Failed],
  $Failed
];
manifestVerified = ListQ[writtenManifestRecords] &&
  Map[
    Replace[Lookup[#, "File", Missing["File"]], File[path_] :> path] &,
    writtenManifestRecords
  ] === manifestRelativeFiles &&
  And @@ Map[
    Function[record,
      With[{
        path = Replace[
          Lookup[record, "File", Missing["File"]],
          File[value_] :> value
        ]
      },
        StringQ[path] &&
          Lookup[record, "Hash", Missing["Hash"]] === FileHash[
            FileNameJoin[{builtPacletRoot, path}],
            "SHA256",
            "HexString"
          ]
      ]
    ],
    writtenManifestRecords
  ];
If[!TrueQ[manifestVerified],
  Print["The final bilingual PacletManifest.wl is incomplete or stale."];
  Exit[14]
];
releaseLog[
  "[5/6] DONE manifest: ",
  Length[manifestRelativeFiles],
  " files (",
  releaseElapsedSeconds[stageStartedAt],
  " s)."
];

stageStartedAt = AbsoluteTime[];
releaseLog["[6/6] START creating the final .paclet archive."];
archivePath = result["PacletArchive"];
If[FileExistsQ[archivePath], DeleteFile[archivePath]];
archivePath = CreatePacletArchive[builtPacletRoot, buildRoot];
If[archivePath === $Failed || !FileExistsQ[archivePath],
  Print["Unable to create the bilingual paclet archive."];
  Exit[7]
];
releaseLog[
  "[6/6] DONE archive (",
  releaseElapsedSeconds[stageStartedAt],
  " s)."
];

Print[InputForm[<|
  "Status" -> "Success",
  "PacletArchive" -> archivePath,
  "EnglishDocumentationPages" -> englishPageCount,
  "ChineseSimplifiedDocumentationPages" -> chinesePageCount,
  "ManifestFiles" -> Length[manifestRelativeFiles]
|>]];
releaseLog[
  "MagneticTB full release build finished in ",
  releaseElapsedSeconds[releaseStartedAt],
  " s."
];
Exit[0];
