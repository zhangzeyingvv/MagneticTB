(* Seconds-only regression tests for the BuildRelease child-process runner.
   This file loads only the release logging and process-runner definitions; it
   never evaluates the release stages or the documentation gate. *)

ClearAll["Global`*"];
Print["START | BuildRelease process lifecycle tests"];

sourcePath = FileNameJoin[{
  DirectoryName[ExpandFileName[$InputFileName]],
  "BuildRelease.wls"
}];
heldExpressions = ReadList[sourcePath, Hold[Expression]];
expressionForms = ToString[#, InputForm] & /@ heldExpressions;
sectionStart = FirstPosition[
  expressionForms,
  form_String /; StringContainsQ[form, "ClearAll[releaseLog,"]
];
sectionEnd = FirstPosition[
  expressionForms,
  form_String /; StringContainsQ[
    form,
    "runStreamingShellScript[script_String, directory_String, OptionsPattern[]]"
  ]
];
If[
  !MatchQ[sectionStart, {_Integer}] ||
    !MatchQ[sectionEnd, {_Integer}] ||
    First[sectionStart] > First[sectionEnd],
  Print["Unable to locate the BuildRelease process-runner section."];
  Exit[2]
];
Scan[
  ReleaseHold,
  Take[heldExpressions, {First[sectionStart], First[sectionEnd]}]
];

ClearAll[recordTest, writeShellScript, readPID, ensureStopped];
passed = 0;
failed = 0;

recordTest[name_String, condition_] := If[TrueQ[condition],
  passed++;
  Print["PASS | ", name],
  failed++;
  Print["FAIL | ", name, " | actual=", InputForm[condition]]
];

writeShellScript[path_String, content_String] := Module[{stream},
  stream = OpenWrite[path];
  If[!MatchQ[stream, _OutputStream], Return[$Failed]];
  WriteString[stream, content];
  Close[stream];
  path
];

readPID[path_String] := If[FileExistsQ[path],
  With[{text = StringTrim[ReadString[path]]},
    If[StringMatchQ[text, DigitCharacter ..], FromDigits[text], $Failed]
  ],
  $Failed
];

ensureStopped[pids_List] := Scan[
  If[IntegerQ[#] && releaseProcessAliveQ[#],
    signalReleaseProcessPIDs[{#}, "KILL"]
  ] &,
  pids
];

testRoot = CreateDirectory[FileNameJoin[{
  $TemporaryDirectory,
  "MagneticTB-BuildReleaseProcessTests-" <> CreateUUID[]
}]];

naturalScript = FileNameJoin[{testRoot, "natural.sh"}];
writeShellScript[naturalScript, "printf '%s\\n' NATURAL_OUTPUT\nexit 7\n"];
naturalResult = runStreamingShellScript[
  naturalScript,
  testRoot,
  "TimeoutSeconds" -> 2.,
  "PollIntervalSeconds" -> 0.01,
  "TerminationGraceSeconds" -> 0.05
];
recordTest[
  "naturally finished child preserves exit code",
  AssociationQ[naturalResult] && naturalResult["ExitCode"] === 7
];

silentScript = FileNameJoin[{testRoot, "silent.sh"}];
writeShellScript[silentScript, "sleep 0.15\nexit 0\n"];
{silentElapsed, silentResult} = AbsoluteTiming@runStreamingShellScript[
  silentScript,
  testRoot,
  "TimeoutSeconds" -> 2.,
  "PollIntervalSeconds" -> 0.01,
  "TerminationGraceSeconds" -> 0.05
];
recordTest[
  "silent child finishes without a ReadString deadlock",
  AssociationQ[silentResult] &&
    silentResult["ExitCode"] === 0 &&
    silentElapsed < 2.
];

streamingScript = FileNameJoin[{testRoot, "streaming.sh"}];
writeShellScript[
  streamingScript,
  "printf '%s\\n' STREAM_FIRST\nsleep 0.15\nprintf '%s\\n' STREAM_LAST\nexit 0\n"
];
{streamingElapsed, streamingResult} = AbsoluteTiming@runStreamingShellScript[
  streamingScript,
  testRoot,
  "TimeoutSeconds" -> 2.,
  "PollIntervalSeconds" -> 0.01,
  "TerminationGraceSeconds" -> 0.05
];
recordTest[
  "multiple output chunks finish without a buffer spin",
  AssociationQ[streamingResult] &&
    streamingResult["ExitCode"] === 0 &&
    streamingElapsed < 2.
];

stdioScript = FileNameJoin[{testRoot, "stdio-boundary.sh"}];
writeShellScript[
  stdioScript,
  "if [ -t 1 ] || [ -t 2 ]; then exit 91; fi\nprintf '%s\\n' CHILD_STDIO_IS_PIPE\nexit 0\n"
];
stdioResult = runStreamingShellScript[
  stdioScript,
  testRoot,
  "TimeoutSeconds" -> 2.,
  "PollIntervalSeconds" -> 0.01,
  "TerminationGraceSeconds" -> 0.05
];
recordTest[
  "StartProcess child uses captured non-PTY stdout and stderr",
  AssociationQ[stdioResult] && stdioResult["ExitCode"] === 0
];

internalParentFile = FileNameJoin[{testRoot, "internal-parent.pid"}];
internalChildFile = FileNameJoin[{testRoot, "internal-child.pid"}];
internalTreeScript = FileNameJoin[{testRoot, "internal-tree.sh"}];
writeShellScript[
  internalTreeScript,
  StringJoin[
    "printf '%s' \"$$\" > \"", internalParentFile, "\"\n",
    "sleep 30 &\n",
    "child=$!\n",
    "printf '%s' \"$child\" > \"", internalChildFile, "\"\n",
    "wait \"$child\"\n"
  ]
];
internalTimeoutResult = runStreamingShellScript[
  internalTreeScript,
  testRoot,
  "TimeoutSeconds" -> 0.4,
  "PollIntervalSeconds" -> 0.01,
  "TerminationGraceSeconds" -> 0.05
];
internalParentPID = readPID[internalParentFile];
internalChildPID = readPID[internalChildFile];
Pause[0.1];
recordTest[
  "internal timeout returns an explicit Failure",
  MatchQ[internalTimeoutResult, Failure["ReleaseGateTimeout", _Association]]
];
recordTest[
  "internal timeout reaps both parent and descendant",
  IntegerQ[internalParentPID] &&
    IntegerQ[internalChildPID] &&
    !releaseProcessAliveQ[internalParentPID] &&
    !releaseProcessAliveQ[internalChildPID]
];
ensureStopped[{internalParentPID, internalChildPID}];

externalParentFile = FileNameJoin[{testRoot, "external-parent.pid"}];
externalChildFile = FileNameJoin[{testRoot, "external-child.pid"}];
externalTreeScript = FileNameJoin[{testRoot, "external-tree.sh"}];
writeShellScript[
  externalTreeScript,
  StringJoin[
    "printf '%s' \"$$\" > \"", externalParentFile, "\"\n",
    "sleep 30 &\n",
    "child=$!\n",
    "printf '%s' \"$child\" > \"", externalChildFile, "\"\n",
    "wait \"$child\"\n"
  ]
];
externalTimeoutResult = TimeConstrained[
  runStreamingShellScript[
    externalTreeScript,
    testRoot,
    "TimeoutSeconds" -> 5.,
    "PollIntervalSeconds" -> 0.01,
    "TerminationGraceSeconds" -> 0.05
  ],
  0.4,
  $OuterTimedOut
];
externalParentPID = readPID[externalParentFile];
externalChildPID = readPID[externalChildFile];
Pause[0.1];
recordTest[
  "external TimeConstrained abort returns its fallback or cleanup Failure",
  externalTimeoutResult === $OuterTimedOut ||
    MatchQ[externalTimeoutResult, Failure["ReleaseGateAborted", _Association]]
];
recordTest[
  "external timeout reaps both parent and descendant",
  IntegerQ[externalParentPID] &&
    IntegerQ[externalChildPID] &&
    !releaseProcessAliveQ[externalParentPID] &&
    !releaseProcessAliveQ[externalChildPID]
];
ensureStopped[{externalParentPID, externalChildPID}];

logResult = TimeConstrained[
  releaseLog["release log probe one"];
  releaseLog["[1/7] START complete documentation release gate."];
  True,
  1.,
  $TimedOut
];
recordTest[
  "successive releaseLog calls finish",
  TrueQ[logResult]
];

DeleteDirectory[testRoot, DeleteContents -> True];
Print["SUMMARY | passed=", passed, " | failed=", failed];
Exit[If[failed === 0, 0, 1]];
