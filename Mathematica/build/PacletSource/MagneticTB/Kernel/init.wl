(* ::Package:: *)

(* Mathematica package entry point. *)
MagneticTB`PackageLoader::corecontext =
  "MagneticTB was not loaded because conflicting core contexts are already present: `1`. Quit the kernel before loading this package; implementations are never mixed or replaced automatically.";

Module[
  {packageDirectory, coreContexts, loadedContexts},
  packageDirectory = DirectoryName[
    DirectoryName[ExpandFileName[$InputFileName]]
  ];
    coreContexts = {
      "MagneticTB`",
      "OperatorSpaceBasis`",
      "MatrixPredicates`",
      "LinearConstraintKernel`",
      "FunctionBasisRepresentation`",
      "GroupEnumeration`",
      "GroupAlgebra`",
      "RepresentationData`",
      "SymmetryAlgebra`",
      "SitePermutationCompiler`",
      "PhysicalRepresentation`",
      "RepresentationActionData`",
      "RepresentationValidation`",
      "DirectProductRepresentation`",
      "InducedRepresentation`",
      "ModelSchema`",
      "ModelPreparation`",
      "MagneticTBLinearAlgebra`"
    };
  loadedContexts = Select[coreContexts, MemberQ[$Packages, #] &];
  If[loadedContexts =!= {},
    Message[MagneticTB`PackageLoader::corecontext, loadedContexts];
    $Failed,
    CompoundExpression[
      Get[FileNameJoin[{packageDirectory, "Usage.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Data", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Core", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Interface", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Validation", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Plot", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "IO", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Utilities", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Corep", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Fitting", "init.wl"}]];
      Get[FileNameJoin[{packageDirectory, "Properties", "init.wl"}]]
    ]
  ]
]
