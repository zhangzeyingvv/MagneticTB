use std::sync::Arc;

use cyclotomic_nullspace::{Element, ExactMatrix, multiply};
use magnetictb_abstract_group::{GroupAction, GroupAlgebra};
use magnetictb_linear_algebra::{conjugate, exact_matrix_equal, exact_unitary_matrix};

use crate::{RepresentationError, RepresentationResult};

include!("representation/representation_data.rs");
include!("representation/direct_product_representation.rs");
include!("representation/induced_representation.rs");
