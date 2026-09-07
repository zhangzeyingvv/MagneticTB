# Deferred Mathematica help pages

These pages have been removed from the active English and Simplified Chinese
help trees. Their existing authoring sources, expected text, and generated
notebooks are preserved here so that they can be rewritten and restored one
at a time.

- `buildSlabHamiltonian`
- `buildRealSpaceHamiltonian`
- `buildBlochHamiltonian`
- `plotBerryCurvature3D`
- `plotBerryCurvature2D`
- `plotRealSpaceWavefunction`
- `plotSurfaceSpectrum`
- `findGaplessPoints`
- `pointChernNumber`
- `plotWilsonLoop`
- `symmetrizationHRInit`
- `berryPhase`
- `berryCurvature`
- `surfaceGreenFunction`
- `surfaceSpectralFunction`
- `wilsonLoop`
- tutorial `ValidatedModels` (`MagneticTB 典型模型`)
- tutorial `WyckoffCrystalSystems` (`七晶系 Wyckoff 数据库模型`)
- tutorial `RealSpaceAndTopology` (`实空间 Hamiltonian、表面与拓扑`)

The `Exact mathematics and cyclotomic fields` section is also temporarily
hidden from the main guide. Its reference pages remain in the active source
tree because only the guide category was requested to be hidden.

The relative directory layout mirrors the active project layout. To restore a
page, first rewrite its English and Chinese sources here. Then move the four
source/expected files back to the matching active directories, remove its name
from `deferredDocumentationPageNames` in
`GenerateDocumentationSources.wls`, and generate the two notebooks again. The
old notebooks stored here are references, not files to publish unchanged.

This folder is not scanned by `BuildNotebooks.sh` and is not part of the
installed MagneticTB help system.
