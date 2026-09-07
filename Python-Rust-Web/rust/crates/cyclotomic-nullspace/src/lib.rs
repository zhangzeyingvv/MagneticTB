#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

pub mod error;
pub mod exact;
pub mod fixture;
pub mod matrix;
mod native_backend;

pub use error::{ExactError, ExactResult};
pub use exact::{
    CyclotomicContext, Element, Polynomial, Rational, cyclotomic_polynomial, embed_element,
    euler_phi, make_context,
};
pub use matrix::{
    CommonKernelResult, CompiledProblem, ExactMatrix, KernelResult, RrefResult, common_kernel,
    multiply, null_space, rref, verify_common_kernel, vertical_stack,
};
