#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::too_many_lines)]

mod error;
mod fitting;
mod vasp;

pub use error::{FittingError, FittingResult};
pub use fitting::{
    AffineBandModel, FitOptions, FitOutput, KPointNeighborhood, SelectionSummary, fit_bands,
};
pub use vasp::{VaspEigenvalueRecord, parse_vasp_eigenval};

pub use nalgebra::{Complex, DMatrix};
