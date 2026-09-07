#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]
// This crate intentionally mirrors Mathematica's machine-real numerical layer.
// Integer translations and bounded mesh indices are validated before these
// explicit f64/isize conversions; the conversions are not exact-core fallbacks.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::missing_panics_doc,
    clippy::too_many_lines
)]

mod berry;
mod brillouin;
mod error;
mod hopping;
mod kpath;
mod legacy;
mod spectral;
mod topology;
mod wannier_io;
mod wilson;

pub use berry::{
    BerryCurvatureData, BerryCurvatureOptions, berry_curvature, berry_curvature_from_samples,
    berry_phase, berry_plaquette_path,
};
pub use brillouin::{
    BrillouinZoneData, FoldedBandSegment, first_brillouin_zone, fold_band_path_to_first_bz,
    reciprocal_lattice,
};
pub use error::{PropertiesError, PropertiesResult};
pub use hopping::{
    BlochHamiltonianData, CellMatrix, FiniteGeometry, FiniteHamiltonianData, HoppingData,
    PrincipalLayerBlocks, SlabHamiltonianData, SurfaceGreenData, SurfaceOptions, SurfaceSide,
    Translation, build_bloch_hamiltonian, build_real_space_hamiltonian, build_slab_hamiltonian,
    principal_layer_blocks, surface_green_function, transform_hoppings,
};
pub use kpath::{BandPathSegment, NamedKPoint, StandardKPathData, standard_k_path};
pub use legacy::{legacy_berryph, legacy_wloop, legacy_z2_path};
pub use spectral::hermitian_eigenvalues;
pub use topology::{
    CubeSurfaceMesh, GaplessPointRecord, GaplessSearchData, GaplessSearchOptions, PointChernData,
    PointChernOptions, RefinementMethod, cube_surface_mesh, find_gapless_points,
    point_chern_number_from_samples,
};
pub use wannier_io::{ParsedWannier90Hr, format_wannier90_hr, parse_wannier90_hr};
pub use wilson::{BerryPhaseData, WilsonLoopData, WilsonLoopOptions, wilson_loop};

pub use nalgebra::{Complex, DMatrix};
