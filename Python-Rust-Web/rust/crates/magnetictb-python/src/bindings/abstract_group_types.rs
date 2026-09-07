#[pyclass(name = "GroupAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyGroupAlgebra {
    inner: CoreGroupAlgebra,
}
#[pymethods]
#[allow(clippy::needless_pass_by_value)]
impl PyGroupAlgebra {
    #[staticmethod]
    fn is_valid_table(multiplication_table: Vec<Vec<usize>>) -> bool {
        CoreGroupAlgebra::is_valid_table(&multiplication_table)
    }

    #[new]
    fn new(multiplication_table: Vec<Vec<usize>>) -> PyResult<Self> {
        CoreGroupAlgebra::new(multiplication_table)
            .map(|inner| Self { inner })
            .map_err(|error| python_group_error(&error))
    }

    #[getter]
    fn multiplication_table(&self) -> Vec<Vec<usize>> {
        self.inner.multiplication_table().to_vec()
    }

    #[getter]
    fn order(&self) -> usize {
        self.inner.order()
    }

    #[getter]
    fn identity(&self) -> usize {
        self.inner.identity()
    }

    #[getter]
    fn inverse_indices(&self) -> Vec<usize> {
        self.inner.inverse_indices().to_vec()
    }

    fn product(&self, left: usize, right: usize) -> PyResult<usize> {
        self.inner
            .product(left, right)
            .map_err(|error| python_group_error(&error))
    }

    fn generated_subgroup(&self, generators: Vec<usize>) -> PyResult<Vec<usize>> {
        self.inner
            .generated_subgroup(&generators)
            .map_err(|error| python_group_error(&error))
    }

    #[pyo3(signature = (subgroup=None))]
    fn find_generator_indices(&self, subgroup: Option<Vec<usize>>) -> PyResult<Vec<usize>> {
        self.inner
            .find_generator_indices(subgroup.as_deref())
            .map_err(|error| python_group_error(&error))
    }

    fn right_coset(&self, subgroup: Vec<usize>, representative: usize) -> PyResult<Vec<usize>> {
        self.inner
            .right_coset(&subgroup, representative)
            .map_err(|error| python_group_error(&error))
    }

    fn left_coset(&self, subgroup: Vec<usize>, representative: usize) -> PyResult<Vec<usize>> {
        self.inner
            .left_coset(&subgroup, representative)
            .map_err(|error| python_group_error(&error))
    }

    fn schreier_decomposition(
        &self,
        subgroup: Vec<usize>,
        representatives: Vec<usize>,
        operation: usize,
        source: usize,
    ) -> PyResult<(usize, usize)> {
        self.inner
            .schreier_decomposition(&subgroup, &representatives, operation, source)
            .map_err(|error| python_group_error(&error))
    }

    fn compile_action(&self, action_table: Vec<Vec<usize>>) -> PyResult<PyGroupAction> {
        CoreGroupAction::compile(&self.inner, action_table)
            .map(|inner| PyGroupAction { inner })
            .map_err(|error| python_group_error(&error))
    }

    fn action_table_is_valid(&self, action_table: Vec<Vec<usize>>) -> bool {
        CoreGroupAction::is_valid_table(&self.inner, &action_table)
    }
}

#[pyclass(name = "GroupAction", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyGroupAction {
    inner: CoreGroupAction,
}

#[pymethods]
impl PyGroupAction {
    #[getter]
    fn action_table(&self) -> Vec<Vec<usize>> {
        self.inner.action_table().to_vec()
    }

    #[getter]
    fn object_count(&self) -> usize {
        self.inner.object_count()
    }

    #[getter]
    fn faithful(&self) -> bool {
        self.inner.is_faithful()
    }

    fn image(&self, operation: usize, source: usize) -> PyResult<usize> {
        self.inner
            .image(operation, source)
            .map_err(|error| python_group_error(&error))
    }

    fn orbit(&self, source: usize) -> PyResult<Vec<usize>> {
        self.inner
            .orbit(source)
            .map_err(|error| python_group_error(&error))
    }

    fn orbits(&self) -> PyResult<Vec<Vec<usize>>> {
        self.inner
            .orbits()
            .map_err(|error| python_group_error(&error))
    }

    fn stabilizer(&self, source: usize) -> PyResult<Vec<usize>> {
        self.inner
            .stabilizer(source)
            .map_err(|error| python_group_error(&error))
    }

    fn transporters(&self, source: usize, target: usize) -> PyResult<Vec<usize>> {
        self.inner
            .transporters(source, target)
            .map_err(|error| python_group_error(&error))
    }
}

#[pyclass(name = "OrderedFiniteGroup", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyOrderedFiniteGroup {
    inner: CoreOrderedFiniteGroup,
}

#[pymethods]
#[allow(clippy::needless_pass_by_value)]
impl PyOrderedFiniteGroup {
    #[new]
    fn new(
        multiplication_table: Vec<Vec<usize>>,
        ordered_elements: Vec<(String, String, bool)>,
    ) -> PyResult<Self> {
        let algebra = CoreGroupAlgebra::new(multiplication_table)
            .map_err(|error| python_group_error(&error))?;
        let elements = ordered_elements
            .into_iter()
            .map(|(stable_id, label, antiunitary)| {
                CoreOrderedElement::new(stable_id, label, antiunitary)
            })
            .collect();
        CoreOrderedFiniteGroup::new(algebra, elements)
            .map(|inner| Self { inner })
            .map_err(|error| python_group_error(&error))
    }

    #[getter]
    fn algebra(&self) -> PyGroupAlgebra {
        PyGroupAlgebra {
            inner: self.inner.algebra().clone(),
        }
    }

    #[getter]
    fn ordered_elements(&self) -> Vec<(String, String, bool)> {
        self.inner
            .elements()
            .iter()
            .map(|element| {
                (
                    element.stable_id().to_owned(),
                    element.label().to_owned(),
                    element.antiunitary(),
                )
            })
            .collect()
    }
}
