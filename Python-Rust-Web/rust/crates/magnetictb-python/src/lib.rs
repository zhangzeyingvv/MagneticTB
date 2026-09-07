use cyclotomic_nullspace::ExactError as CoreExactError;
use cyclotomic_nullspace::fixture::{
    common_result_to_json, context_from_json, context_to_json, element_from_json, element_to_json,
    kernel_result_to_json, matrix_from_json, matrix_to_json, object, problem_from_json, required,
};
use cyclotomic_nullspace::{common_kernel, make_context, null_space};
use magnetictb_abstract_group::{
    GroupAction as CoreGroupAction, GroupAlgebra as CoreGroupAlgebra, GroupError as CoreGroupError,
    OrderedElement as CoreOrderedElement, OrderedFiniteGroup as CoreOrderedFiniteGroup,
    find_concrete_generators_by, generate_group_by,
};
use magnetictb_crystal_geometry::compile_msg_wyckoff_sites;
use magnetictb_data::{
    DataCatalog as CoreDataCatalog, DataError as CoreDataError, MsgGroup as CoreMsgGroup,
};
use magnetictb_fitting::FittingError as CoreFittingError;
use magnetictb_properties::{
    Complex as NumericComplex, DMatrix as NumericMatrix, PropertiesError as CorePropertiesError,
    find_gapless_points as core_find_gapless_points,
};
use magnetictb_symmetry::{SymmetryError as CoreSymmetryError, compile_msg_group};
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use serde::Serialize;
use serde_json::Value;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::Arc;

mod fitting_boundary;
mod geometry_boundary;
mod linear_algebra_boundary;
mod physical_boundary;
mod properties_boundary;
mod representation_boundary;
mod tight_binding_boundary;

pyo3::create_exception!(_core, ExactError, PyException);
pyo3::create_exception!(_core, GroupError, PyException);
pyo3::create_exception!(_core, DataError, PyException);
pyo3::create_exception!(_core, SymmetryError, PyException);
pyo3::create_exception!(_core, PropertiesError, PyException);
pyo3::create_exception!(_core, FittingError, PyException);

include!("bindings/shared.rs");
include!("bindings/data.rs");
include!("bindings/abstract_group_types.rs");
include!("bindings/core_dispatch.rs");
include!("bindings/abstract_group_callbacks.rs");

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("ExactError", module.py().get_type::<ExactError>())?;
    module.add("GroupError", module.py().get_type::<GroupError>())?;
    module.add("DataError", module.py().get_type::<DataError>())?;
    module.add("SymmetryError", module.py().get_type::<SymmetryError>())?;
    module.add("PropertiesError", module.py().get_type::<PropertiesError>())?;
    module.add("FittingError", module.py().get_type::<FittingError>())?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    module.add_function(wrap_pyfunction!(null_space_json, module)?)?;
    module.add_function(wrap_pyfunction!(common_kernel_json, module)?)?;
    module.add_function(wrap_pyfunction!(cyclotomic_context_json, module)?)?;
    module.add_function(wrap_pyfunction!(linear_algebra_json, module)?)?;
    module.add_function(wrap_pyfunction!(representation_json, module)?)?;
    module.add_function(wrap_pyfunction!(tight_binding_json, module)?)?;
    module.add_function(wrap_pyfunction!(geometry_json, module)?)?;
    module.add_function(wrap_pyfunction!(properties_json, module)?)?;
    module.add_function(wrap_pyfunction!(fitting_json, module)?)?;
    module.add_function(wrap_pyfunction!(gapless_points_json, module)?)?;
    module.add_function(wrap_pyfunction!(generate_group_objects, module)?)?;
    module.add_function(wrap_pyfunction!(find_generator_objects, module)?)?;
    module.add_class::<PyGroupAlgebra>()?;
    module.add_class::<PyGroupAction>()?;
    module.add_class::<PyOrderedFiniteGroup>()?;
    module.add_class::<PyDataCatalog>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{compute_common_kernel, compute_null_space};
    use serde_json::json;

    #[test]
    fn malformed_json_fails_explicitly() {
        let error = compute_null_space("{").expect_err("malformed JSON must fail");
        assert_eq!(error.tag(), "MalformedSerialization");
    }

    #[test]
    fn empty_common_kernel_uses_rust_core() {
        let request = json!({
            "context": {
                "conductor": 1,
                "degree": 1,
                "cyclotomic_polynomial": [
                    {"numerator": "-1", "denominator": "1"},
                    {"numerator": "1", "denominator": "1"}
                ]
            },
            "coordinate_dimension": 0,
            "constraints": []
        });
        let output = compute_common_kernel(&request.to_string()).expect("exact common kernel");
        let value: serde_json::Value = serde_json::from_str(&output).expect("result JSON");
        assert_eq!(value["rank"], 0);
        assert_eq!(value["nullity"], 0);
        assert_eq!(value["iteration_nullities"], json!([0]));
        assert_eq!(value["exact_residual_verified"], true);
    }
}
