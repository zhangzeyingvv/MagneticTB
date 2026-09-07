#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PolynomialExpression {
    Scalar(Element),
    Variable(usize),
    Add(Vec<Self>),
    Multiply(Vec<Self>),
    Power(Box<Self>, u8),
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledBasisAction {
    local_matrices: Vec<ExactMatrix>,
    spinor: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledInducedBasisAction {
    transported_basis_states: Vec<Vec<FunctionBasisState>>,
    local_blocks: Vec<Vec<ExactMatrix>>,
    spinor: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PolynomialTerm {
    exponents: [u8; 3],
    coefficient: Element,
}

impl PolynomialTerm {
    #[must_use]
    pub const fn exponents(&self) -> &[u8; 3] {
        &self.exponents
    }

    #[must_use]
    pub const fn coefficient(&self) -> &Element {
        &self.coefficient
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FunctionBasisState {
    components: Vec<Vec<PolynomialTerm>>,
}

impl FunctionBasisState {
    #[must_use]
    pub fn components(&self) -> &[Vec<PolynomialTerm>] {
        &self.components
    }
}

impl CompiledBasisAction {
    #[must_use]
    pub fn local_matrices(&self) -> &[ExactMatrix] {
        &self.local_matrices
    }

    #[must_use]
    pub const fn spinor(&self) -> bool {
        self.spinor
    }
}

impl CompiledInducedBasisAction {
    #[must_use]
    pub fn transported_basis_states(&self) -> &[Vec<FunctionBasisState>] {
        &self.transported_basis_states
    }

    #[must_use]
    pub fn local_blocks(&self) -> &[Vec<ExactMatrix>] {
        &self.local_blocks
    }

    #[must_use]
    pub const fn spinor(&self) -> bool {
        self.spinor
    }
}

fn basis_error(tag: &str, detail: impl Into<String>) -> SymmetryError {
    SymmetryError::new(tag, detail)
}

fn scalar(context: &Arc<CyclotomicContext>, value: i64) -> SymmetryResult<Element> {
    Element::from_polynomial(Arc::clone(context), &[Rational::from_i64(value)]).map_err(Into::into)
}

fn rational(
    context: &Arc<CyclotomicContext>,
    numerator: i64,
    denominator: i64,
) -> SymmetryResult<Element> {
    Element::from_polynomial(
        Arc::clone(context),
        &[Rational::parse(
            &numerator.to_string(),
            &denominator.to_string(),
        )?],
    )
    .map_err(Into::into)
}

fn monomial(coefficient: Element, exponents: Exponents) -> Polynomial {
    BTreeMap::from([(exponents, coefficient)])
}

fn constant(context: &Arc<CyclotomicContext>, value: i64) -> SymmetryResult<Polynomial> {
    Ok(monomial(scalar(context, value)?, [0, 0, 0]))
}

fn variable(context: &Arc<CyclotomicContext>, index: usize) -> SymmetryResult<Polynomial> {
    let mut exponents = [0, 0, 0];
    exponents[index] = 1;
    Ok(monomial(Element::one(context)?, exponents))
}

fn polynomial_add(left: &Polynomial, right: &Polynomial) -> SymmetryResult<Polynomial> {
    let mut result = left.clone();
    for (exponents, coefficient) in right {
        let updated = match result.get(exponents) {
            Some(current) => current.add(coefficient)?,
            None => coefficient.clone(),
        };
        if updated.is_zero() {
            result.remove(exponents);
        } else {
            result.insert(*exponents, updated);
        }
    }
    Ok(result)
}

fn polynomial_scale(polynomial: &Polynomial, coefficient: &Element) -> SymmetryResult<Polynomial> {
    let mut result = Polynomial::new();
    for (exponents, value) in polynomial {
        let product = value.multiply(coefficient)?;
        if !product.is_zero() {
            result.insert(*exponents, product);
        }
    }
    Ok(result)
}

fn polynomial_multiply(left: &Polynomial, right: &Polynomial) -> SymmetryResult<Polynomial> {
    let mut result = Polynomial::new();
    for (left_exponents, left_value) in left {
        for (right_exponents, right_value) in right {
            let exponents = [
                left_exponents[0] + right_exponents[0],
                left_exponents[1] + right_exponents[1],
                left_exponents[2] + right_exponents[2],
            ];
            let product = left_value.multiply(right_value)?;
            let updated = match result.get(&exponents) {
                Some(current) => current.add(&product)?,
                None => product,
            };
            if updated.is_zero() {
                result.remove(&exponents);
            } else {
                result.insert(exponents, updated);
            }
        }
    }
    Ok(result)
}

fn expression_polynomial(
    context: &Arc<CyclotomicContext>,
    expression: &PolynomialExpression,
) -> SymmetryResult<Polynomial> {
    match expression {
        PolynomialExpression::Scalar(value) => {
            if value.context() != context {
                return Err(basis_error(
                    "InvalidFunctionBasis",
                    "polynomial scalar uses a different exact field",
                ));
            }
            Ok(monomial(value.clone(), [0, 0, 0]))
        }
        PolynomialExpression::Variable(index) if *index < 3 => variable(context, *index),
        PolynomialExpression::Variable(_) => Err(basis_error(
            "InvalidFunctionBasis",
            "only the ordered Cartesian variables x, y, z are supported",
        )),
        PolynomialExpression::Add(terms) => {
            terms.iter().try_fold(Polynomial::new(), |total, term| {
                polynomial_add(&total, &expression_polynomial(context, term)?)
            })
        }
        PolynomialExpression::Multiply(factors) => {
            factors
                .iter()
                .try_fold(constant(context, 1)?, |product, factor| {
                    polynomial_multiply(&product, &expression_polynomial(context, factor)?)
                })
        }
        PolynomialExpression::Power(base, exponent) => {
            let base = expression_polynomial(context, base)?;
            (0..*exponent).try_fold(constant(context, 1)?, |product, _| {
                polynomial_multiply(&product, &base)
            })
        }
    }
}

fn expression_basis(
    context: &Arc<CyclotomicContext>,
    expressions: &[Vec<PolynomialExpression>],
) -> SymmetryResult<(Vec<BasisFunction>, bool)> {
    let component_count = expressions.first().map_or(0, Vec::len);
    if expressions.is_empty()
        || !matches!(component_count, 1 | 2)
        || expressions
            .iter()
            .any(|function| function.len() != component_count)
    {
        return Err(basis_error(
            "InvalidFunctionBasis",
            "an exact function basis must be a nonempty rectangular scalar or two-component list",
        ));
    }
    let basis = expressions
        .iter()
        .map(|function| {
            function
                .iter()
                .map(|expression| expression_polynomial(context, expression))
                .collect()
        })
        .collect::<SymmetryResult<Vec<_>>>()?;
    Ok((basis, component_count == 2))
}

fn polynomial_conjugate(polynomial: &Polynomial) -> SymmetryResult<Polynomial> {
    polynomial
        .iter()
        .map(|(exponents, coefficient)| Ok((*exponents, coefficient.conjugate()?)))
        .collect()
}

fn linear_form(matrix: &ExactMatrix, row: usize) -> SymmetryResult<Polynomial> {
    let mut result = Polynomial::new();
    for column in 0..3 {
        if !matrix.entry(row, column)?.is_zero() {
            let mut exponents = [0, 0, 0];
            exponents[column] = 1;
            result.insert(exponents, matrix.entry(row, column)?.clone());
        }
    }
    Ok(result)
}

fn substitute(
    polynomial: &Polynomial,
    coordinate_images: &ExactMatrix,
) -> SymmetryResult<Polynomial> {
    let context = coordinate_images.context();
    let forms = [
        linear_form(coordinate_images, 0)?,
        linear_form(coordinate_images, 1)?,
        linear_form(coordinate_images, 2)?,
    ];
    let mut result = Polynomial::new();
    for (exponents, coefficient) in polynomial {
        let mut term = constant(context, 1)?;
        for variable_index in 0..3 {
            for _ in 0..exponents[variable_index] {
                term = polynomial_multiply(&term, &forms[variable_index])?;
            }
        }
        result = polynomial_add(&result, &polynomial_scale(&term, coefficient)?)?;
    }
    Ok(result)
}

fn matrix_inverse(matrix: &ExactMatrix) -> SymmetryResult<ExactMatrix> {
    if matrix.rows() == 0 || matrix.rows() != matrix.columns() {
        return Err(basis_error(
            "BasisCoordinateTransform",
            "matrix inverse requires a nonempty square exact matrix",
        ));
    }
    let dimension = matrix.rows();
    let zero = Element::zero(matrix.context())?;
    let one = Element::one(matrix.context())?;
    let mut entries = Vec::with_capacity(dimension * dimension * 2);
    for row in 0..dimension {
        for column in 0..dimension {
            entries.push(matrix.entry(row, column)?.clone());
        }
        for column in 0..dimension {
            entries.push(if row == column {
                one.clone()
            } else {
                zero.clone()
            });
        }
    }
    let augmented = ExactMatrix::new(
        Arc::clone(matrix.context()),
        dimension,
        2 * dimension,
        entries,
    )?;
    let reduced = rref(&augmented)?;
    if reduced.rank != dimension || reduced.pivot_columns != (0..dimension).collect::<Vec<_>>() {
        return Err(basis_error(
            "BasisCoordinateTransform",
            "exact matrix is singular",
        ));
    }
    let mut inverse_entries = Vec::with_capacity(dimension * dimension);
    for row in 0..dimension {
        for column in 0..dimension {
            inverse_entries.push(
                reduced
                    .reduced_matrix
                    .entry(row, dimension + column)?
                    .clone(),
            );
        }
    }
    ExactMatrix::new(
        Arc::clone(matrix.context()),
        dimension,
        dimension,
        inverse_entries,
    )
    .map_err(Into::into)
}

fn coordinate_images(rotation: &ExactMatrix, lattice: &ExactMatrix) -> SymmetryResult<ExactMatrix> {
    if rotation.rows() != 3
        || rotation.columns() != 3
        || lattice.rows() != 3
        || lattice.columns() != 3
        || rotation.context() != lattice.context()
    {
        return Err(basis_error(
            "BasisCoordinateTransform",
            "rotation and lattice must be exact 3 by 3 matrices in one field",
        ));
    }
    let lattice_transpose = lattice.transpose()?;
    let inverse_rotation = matrix_inverse(rotation)?;
    let inverse_lattice_transpose = matrix_inverse(&lattice_transpose)?;
    Ok(multiply(
        &multiply(&lattice_transpose, &inverse_rotation)?,
        &inverse_lattice_transpose,
    )?)
}

fn apply_value_matrix(
    function: &BasisFunction,
    matrix: &ExactMatrix,
) -> SymmetryResult<BasisFunction> {
    if matrix.rows() != function.len() || matrix.columns() != function.len() {
        return Err(basis_error(
            "InvalidFunctionBasisOperation",
            "value-space matrix does not match the basis value dimension",
        ));
    }
    let mut result = Vec::with_capacity(function.len());
    for row in 0..matrix.rows() {
        let mut component = Polynomial::new();
        for (column, polynomial) in function.iter().enumerate() {
            component = polynomial_add(
                &component,
                &polynomial_scale(polynomial, matrix.entry(row, column)?)?,
            )?;
        }
        result.push(component);
    }
    Ok(result)
}

fn transform_basis(
    basis: &[BasisFunction],
    coordinate_images: &ExactMatrix,
    value_matrix: &ExactMatrix,
    antiunitary: bool,
) -> SymmetryResult<Vec<BasisFunction>> {
    let value_dimension = basis[0].len();
    let time_reversal = if antiunitary && value_dimension == 2 {
        let zero = Element::zero(value_matrix.context())?;
        let one = Element::one(value_matrix.context())?;
        let minus_one = one.negate()?;
        Some(ExactMatrix::new(
            Arc::clone(value_matrix.context()),
            2,
            2,
            vec![zero.clone(), one, minus_one, zero],
        )?)
    } else {
        None
    };
    let combined = match time_reversal {
        Some(matrix) => multiply(value_matrix, &matrix)?,
        None => value_matrix.clone(),
    };
    basis
        .iter()
        .map(|function| {
            let conjugated = if antiunitary {
                function
                    .iter()
                    .map(polynomial_conjugate)
                    .collect::<SymmetryResult<Vec<_>>>()?
            } else {
                function.clone()
            };
            let valued = apply_value_matrix(&conjugated, &combined)?;
            valued
                .iter()
                .map(|component| substitute(component, coordinate_images))
                .collect()
        })
        .collect()
}

fn snapshot_basis_function(function: &BasisFunction) -> FunctionBasisState {
    FunctionBasisState {
        components: function
            .iter()
            .map(|component| {
                component
                    .iter()
                    .map(|(exponents, coefficient)| PolynomialTerm {
                        exponents: *exponents,
                        coefficient: coefficient.clone(),
                    })
                    .collect()
            })
            .collect(),
    }
}

fn snapshot_basis(basis: &[BasisFunction]) -> Vec<FunctionBasisState> {
    basis.iter().map(snapshot_basis_function).collect()
}

fn transformed_basis_images(
    basis: &[BasisFunction],
    spinor: bool,
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> SymmetryResult<Vec<Vec<BasisFunction>>> {
    if spatial_operations.is_empty()
        || spatial_operations.len() != antiunitary_flags.len()
        || spin_rotations.is_some_and(|rotations| rotations.len() != spatial_operations.len())
    {
        return Err(basis_error(
            "InvalidFunctionBasisOperation",
            "spatial, spin, and antiunitary operation lists must align",
        ));
    }
    spatial_operations
        .iter()
        .enumerate()
        .map(|(index, spatial)| {
            transform_basis_for_point_operation(
                basis,
                spinor,
                spatial,
                spin_rotations.map(|rotations| &rotations[index]),
                antiunitary_flags[index],
                lattice,
            )
        })
        .collect()
}

fn transform_basis_for_point_operation(
    basis: &[BasisFunction],
    spinor: bool,
    spatial: &SpatialOperation,
    spin_rotation: Option<&ExactMatrix>,
    antiunitary: bool,
    lattice: &ExactMatrix,
) -> SymmetryResult<Vec<BasisFunction>> {
    let coordinate_images = coordinate_images(spatial.rotation(), lattice)?;
    let value_matrix = if spinor {
        if let Some(rotation) = spin_rotation {
            spin_matrix(rotation)?
        } else {
            compile_spatial_spinor_matrix(spatial.rotation(), lattice)?
        }
    } else {
        ExactMatrix::identity(lattice.context(), 1)?
    };
    transform_basis(basis, &coordinate_images, &value_matrix, antiunitary)
}

fn same_function(left: &BasisFunction, right: &BasisFunction) -> bool {
    left == right
}

#[allow(clippy::too_many_lines)]
fn representation_matrix(
    basis: &[BasisFunction],
    images: &[BasisFunction],
) -> SymmetryResult<ExactMatrix> {
    let context = basis[0][0]
        .values()
        .next()
        .or_else(|| {
            basis
                .iter()
                .flat_map(|function| function.iter())
                .find_map(|p| p.values().next())
        })
        .ok_or_else(|| basis_error("InvalidFunctionBasis", "basis contains only zero functions"))?
        .context()
        .clone();
    let mut multiplicities = Vec::with_capacity(basis.len());
    for index in 0..basis.len() {
        multiplicities.push(
            basis[..=index]
                .iter()
                .filter(|candidate| same_function(candidate, &basis[index]))
                .count(),
        );
    }
    let unique_multiplicities = {
        let mut values = multiplicities.clone();
        values.sort_unstable();
        values.dedup();
        values
    };
    let mut exponents = basis
        .iter()
        .chain(images)
        .flat_map(|function| function.iter())
        .flat_map(|polynomial| polynomial.keys().copied())
        .collect::<Vec<_>>();
    exponents.sort_unstable();
    exponents.dedup();
    let value_dimension = basis[0].len();
    let mut row_keys = Vec::new();
    for multiplicity in unique_multiplicities {
        for component in 0..value_dimension {
            for exponent in &exponents {
                row_keys.push((multiplicity, component, *exponent));
            }
        }
    }
    let coefficient_matrix = |functions: &[BasisFunction]| -> SymmetryResult<ExactMatrix> {
        let mut entries = Vec::with_capacity(row_keys.len() * functions.len());
        for (multiplicity, component, exponent) in &row_keys {
            for (column, function) in functions.iter().enumerate() {
                entries.push(if multiplicities[column] == *multiplicity {
                    function[*component]
                        .get(exponent)
                        .cloned()
                        .unwrap_or(Element::zero(&context)?)
                } else {
                    Element::zero(&context)?
                });
            }
        }
        ExactMatrix::new(
            Arc::clone(&context),
            row_keys.len(),
            functions.len(),
            entries,
        )
        .map_err(Into::into)
    };
    let basis_coefficients = coefficient_matrix(basis)?;
    let image_coefficients = coefficient_matrix(images)?;
    let reduced_transpose = rref(&basis_coefficients.transpose()?)?;
    let pivot_rows = reduced_transpose.pivot_columns;
    if pivot_rows.len() != basis.len() {
        return Err(basis_error(
            "DependentFunctionBasis",
            "ordered function basis is linearly dependent",
        ));
    }
    let square_rows = |matrix: &ExactMatrix| -> SymmetryResult<ExactMatrix> {
        let mut entries = Vec::with_capacity(pivot_rows.len() * matrix.columns());
        for row in &pivot_rows {
            for column in 0..matrix.columns() {
                entries.push(matrix.entry(*row, column)?.clone());
            }
        }
        ExactMatrix::new(
            Arc::clone(&context),
            pivot_rows.len(),
            matrix.columns(),
            entries,
        )
        .map_err(Into::into)
    };
    let result = multiply(
        &matrix_inverse(&square_rows(&basis_coefficients)?)?,
        &square_rows(&image_coefficients)?,
    )?;
    if !exact_matrix_equal(
        &multiply(&basis_coefficients, &result)?,
        &image_coefficients,
    ) {
        return Err(basis_error(
            "FunctionBasisNotClosed",
            "a transformed function is outside the ordered basis span",
        ));
    }
    Ok(result)
}

/// Compiles the Mathematica `BasisCatalog` and `FunctionBasisRepresentation`
/// convention into exact ordered local matrices.  `spin_rotations` is supplied
/// for spin-space groups; otherwise the axial part of the Cartesian spatial
/// operation is used for spinors.
pub fn compile_catalog_basis_action(
    labels: &[String],
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> SymmetryResult<CompiledBasisAction> {
    let (basis, spinor) = catalog_basis(lattice.context(), labels)?;
    compile_basis_action(
        &basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        lattice,
    )
}

/// Compile the stable Mathematica function-basis Induced path.  The ordered
/// transported bases are fixed by `coset_representatives`; each exact local
/// block maps the source transported basis to the target selected by the
/// operation-by-source site action.  This deliberately does not impose an
/// ordinary group representation law on double-valued spinor bases.
pub fn compile_catalog_induced_basis_action(
    labels: &[String],
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    image_site_indices: &[Vec<usize>],
    coset_representatives: &[usize],
    lattice: &ExactMatrix,
) -> SymmetryResult<CompiledInducedBasisAction> {
    let (basis, spinor) = catalog_basis(lattice.context(), labels)?;
    compile_induced_basis_action(
        &basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        image_site_indices,
        coset_representatives,
        lattice,
    )
}

/// Resolve ordered catalog labels to their exact scalar or spinor polynomial
/// states without applying a symmetry operation.
pub fn resolve_catalog_basis_states(
    labels: &[String],
    lattice: &ExactMatrix,
) -> SymmetryResult<Vec<FunctionBasisState>> {
    let (basis, _) = catalog_basis(lattice.context(), labels)?;
    Ok(snapshot_basis(&basis))
}

/// Apply the supplied ordered point operations to catalog basis functions and
/// return the exact transported states. The outer order is operation order and
/// the inner order is basis-function order.
pub fn transform_catalog_basis_states(
    labels: &[String],
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> SymmetryResult<Vec<Vec<FunctionBasisState>>> {
    let (basis, spinor) = catalog_basis(lattice.context(), labels)?;
    Ok(transformed_basis_images(
        &basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        lattice,
    )?
    .iter()
    .map(|images| snapshot_basis(images))
    .collect())
}

/// Compiles an explicitly encoded exact polynomial function basis.  The
/// outer order is the basis-function order and the inner order is the scalar
/// or two-component spinor value order used by stable Mathematica.
pub fn compile_function_basis_action(
    expressions: &[Vec<PolynomialExpression>],
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> SymmetryResult<CompiledBasisAction> {
    let (basis, spinor) = expression_basis(lattice.context(), expressions)?;
    compile_basis_action(
        &basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        lattice,
    )
}

/// Compile an explicitly encoded exact polynomial basis through the stable
/// transported-basis/local-block Induced path.
pub fn compile_function_induced_basis_action(
    expressions: &[Vec<PolynomialExpression>],
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    image_site_indices: &[Vec<usize>],
    coset_representatives: &[usize],
    lattice: &ExactMatrix,
) -> SymmetryResult<CompiledInducedBasisAction> {
    let (basis, spinor) = expression_basis(lattice.context(), expressions)?;
    compile_induced_basis_action(
        &basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        image_site_indices,
        coset_representatives,
        lattice,
    )
}

/// Resolve an explicitly encoded exact polynomial basis without applying a
/// symmetry operation.
pub fn resolve_function_basis_states(
    expressions: &[Vec<PolynomialExpression>],
    lattice: &ExactMatrix,
) -> SymmetryResult<Vec<FunctionBasisState>> {
    let (basis, _) = expression_basis(lattice.context(), expressions)?;
    Ok(snapshot_basis(&basis))
}

/// Apply the supplied ordered point operations to an explicitly encoded exact
/// polynomial basis and return exact transported states.
pub fn transform_function_basis_states(
    expressions: &[Vec<PolynomialExpression>],
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> SymmetryResult<Vec<Vec<FunctionBasisState>>> {
    let (basis, spinor) = expression_basis(lattice.context(), expressions)?;
    Ok(transformed_basis_images(
        &basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        lattice,
    )?
    .iter()
    .map(|images| snapshot_basis(images))
    .collect())
}

fn compile_basis_action(
    basis: &[BasisFunction],
    spinor: bool,
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> SymmetryResult<CompiledBasisAction> {
    let mut local_matrices = Vec::with_capacity(spatial_operations.len());
    for (index, images) in transformed_basis_images(
        basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        lattice,
    )?
    .iter()
    .enumerate()
    {
        let matrix = representation_matrix(basis, images)?;
        if !magnetictb_linear_algebra::exact_unitary_matrix(&matrix)? {
            return Err(basis_error(
                "NonunitaryFunctionBasis",
                format!("basis action {index} is not exactly unitary"),
            ));
        }
        local_matrices.push(matrix);
    }
    Ok(CompiledBasisAction {
        local_matrices,
        spinor,
    })
}

#[allow(clippy::too_many_arguments)]
fn compile_induced_basis_action(
    basis: &[BasisFunction],
    spinor: bool,
    spatial_operations: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    image_site_indices: &[Vec<usize>],
    coset_representatives: &[usize],
    lattice: &ExactMatrix,
) -> SymmetryResult<CompiledInducedBasisAction> {
    let site_count = coset_representatives.len();
    if site_count == 0
        || spatial_operations.len() != image_site_indices.len()
        || spatial_operations.len() != antiunitary_flags.len()
        || spin_rotations.is_some_and(|rotations| rotations.len() != spatial_operations.len())
        || coset_representatives
            .iter()
            .any(|&operation| operation >= spatial_operations.len())
    {
        return Err(basis_error(
            "InvalidInducedFunctionBasis",
            "operations, site action, and coset representatives must align",
        ));
    }
    for images in image_site_indices {
        let mut ordered = images.clone();
        ordered.sort_unstable();
        if ordered != (0..site_count).collect::<Vec<_>>() {
            return Err(basis_error(
                "InvalidInducedFunctionBasis",
                "every operation must permute all generated sites exactly once",
            ));
        }
    }

    let transformed_reference_bases = transformed_basis_images(
        basis,
        spinor,
        spatial_operations,
        spin_rotations,
        antiunitary_flags,
        lattice,
    )?;
    let transported_bases = coset_representatives
        .iter()
        .map(|&operation| transformed_reference_bases[operation].clone())
        .collect::<Vec<_>>();
    let mut local_blocks = Vec::with_capacity(spatial_operations.len());
    for (operation, spatial) in spatial_operations.iter().enumerate() {
        let mut operation_blocks = Vec::with_capacity(site_count);
        for (source, source_basis) in transported_bases.iter().enumerate() {
            let target = image_site_indices[operation][source];
            let image = transform_basis_for_point_operation(
                source_basis,
                spinor,
                spatial,
                spin_rotations.map(|rotations| &rotations[operation]),
                antiunitary_flags[operation],
                lattice,
            )?;
            let block = representation_matrix(&transported_bases[target], &image)?;
            if !magnetictb_linear_algebra::exact_unitary_matrix(&block)? {
                return Err(basis_error(
                    "NonunitaryFunctionBasis",
                    format!(
                        "induced basis block for operation {operation}, source {source} is not exactly unitary"
                    ),
                ));
            }
            operation_blocks.push(block);
        }
        local_blocks.push(operation_blocks);
    }
    Ok(CompiledInducedBasisAction {
        transported_basis_states: transported_bases
            .iter()
            .map(|basis| snapshot_basis(basis))
            .collect(),
        local_blocks,
        spinor,
    })
}
