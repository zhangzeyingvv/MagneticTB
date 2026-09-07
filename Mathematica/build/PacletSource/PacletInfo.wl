(* BuildRelease.wls stages the repository's MagneticTB/Old/MagneticTBOld
   directory as the sibling release root MagneticTBOld. This keeps the two
   Kernel extensions non-overlapping without changing the source archive. *)
PacletObject[<|
  "Name" -> "MagneticTB",
  "Version" -> "2.0.10",
  "WolframVersion" -> "12.1+",
  "Description" ->
    "Construct symmetry-constrained tight-binding Hamiltonians for magnetic, nonmagnetic, and spin-space-group models.",
  "Creator" -> "Zeying Zhang",
  "URL" -> "https://github.com/zhangzeyingvv/MagneticTB",
  "Extensions" -> {
    {
      "Kernel",
      "Root" -> "MagneticTB",
      "Context" -> {"MagneticTB`"}
    },
    {
      "Kernel",
      "Root" -> "MagneticTBOld",
      "Context" -> {"MagneticTBOld`"}
    },
    {
      "Documentation",
      "Language" -> All,
      "MainPage" -> "Guides/MagneticTB"
    }
  }
|>]
