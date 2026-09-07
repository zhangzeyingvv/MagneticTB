// The boundary remains one Rust module so private helpers and error semantics
// are unchanged.  The included files follow the Mathematica package domains.
include!("geometry_boundary/shared.rs");
include!("geometry_boundary/tight_binding.rs");
include!("geometry_boundary/properties.rs");
include!("geometry_boundary/physical_representation_support.rs");
include!("geometry_boundary/crystal_geometry_support.rs");
include!("geometry_boundary/tight_binding_support.rs");
include!("geometry_boundary/symmetry_support.rs");
include!("geometry_boundary/crystal_geometry.rs");
include!("geometry_boundary/lattice_input.rs");
include!("geometry_boundary/model_compiler.rs");
include!("geometry_boundary/dispatch.rs");
