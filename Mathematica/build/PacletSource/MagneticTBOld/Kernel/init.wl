(* ::Package:: *)

(* Mathematica Init File *)
Scan[
  Get[FileNameJoin[{
    DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]], #
  }]] &,
  {
    "Usage.wl", "MagneticTB.wl", "Plot.wl", "IO.m", "Utilities.wl",
    "Corep.wl", "SSG.wl", "Fitting.wl", "WilsonLoop.wl"
  }
]
