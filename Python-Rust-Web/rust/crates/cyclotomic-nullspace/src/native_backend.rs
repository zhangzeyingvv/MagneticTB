//! Native raw backends for the Mathematica coefficient/rational work-plan contract.

use crate::error::{ExactResult, fail};
use crate::exact::{
    CyclotomicContext, Element, Polynomial, Rational, polynomial_extended_gcd, polynomial_pad,
};
use crate::matrix::ExactMatrix;
use num_bigint::BigInt;
use num_traits::Zero;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

const ADAPTIVE_INITIAL_DENSE_ROWS: usize = 20;
const ADAPTIVE_MINIMUM_DENSE_ROWS: usize = 8;
const ADAPTIVE_MAXIMUM_DENSE_ROWS: usize = 40;
const ADAPTIVE_MAXIMUM_BLOCK_ROWS: usize = 40;
const INVERSE_CACHE_MAXIMUM_ENTRIES: usize = 4096;

pub(crate) type Coefficient = Vec<Rational>;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum CoefficientFieldKind {
    FullCyclotomicField,
    MaximalRealSubfield,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct CoefficientFieldContext {
    pub field_kind: CoefficientFieldKind,
    pub field_id: String,
    pub degree: usize,
    pub defining_polynomial: Polynomial,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct RawMatrix<S> {
    pub rows: usize,
    pub columns: usize,
    pub data: Vec<S>,
}

impl<S> RawMatrix<S> {
    fn new(rows: usize, columns: usize, data: Vec<S>) -> ExactResult<Self> {
        let expected = rows.checked_mul(columns).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "raw matrix shape overflows")
        })?;
        if data.len() != expected {
            return fail(
                "MalformedSerialization",
                "raw matrix entry count differs from shape",
            );
        }
        Ok(Self {
            rows,
            columns,
            data,
        })
    }

    fn at(&self, row: usize, column: usize) -> &S {
        &self.data[row * self.columns + column]
    }
}

pub(crate) type RationalRawMatrix = RawMatrix<Rational>;
pub(crate) type CoefficientRawMatrix = RawMatrix<Coefficient>;

#[derive(Clone, Debug)]
pub(crate) struct RawRrefResult<S> {
    pub reduced: RawMatrix<S>,
    pub pivot_columns: Vec<usize>,
    pub free_columns: Vec<usize>,
    pub rank: usize,
    pub nullity: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct RawKernelResult<S> {
    pub rref: RawRrefResult<S>,
    pub nullspace_rows: RawMatrix<S>,
    pub basis_matrix: RawMatrix<S>,
}

#[derive(Clone, Debug)]
pub(crate) struct RawCommonKernelResult<S> {
    pub basis_matrix: RawMatrix<S>,
    pub iteration_nullities: Vec<usize>,
    pub rank: usize,
    pub nullity: usize,
    pub exact_residual_verified: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct CoefficientWorkPlan {
    pub public_context: Arc<CyclotomicContext>,
    pub field_context: Arc<CoefficientFieldContext>,
    pub coordinate_dimension: usize,
    pub raw_constraints: Vec<CoefficientRawMatrix>,
    pub lift_matrix: Vec<Rational>,
    pub lift_rows: usize,
    pub lift_columns: usize,
    pub method: String,
}

trait RawDomain {
    type Scalar: Clone + Eq;

    fn zero(&self) -> Self::Scalar;
    fn one(&self) -> Self::Scalar;
    fn is_zero(&self, value: &Self::Scalar) -> bool;
    fn add(&mut self, left: &Self::Scalar, right: &Self::Scalar) -> ExactResult<Self::Scalar>;
    fn subtract(&mut self, left: &Self::Scalar, right: &Self::Scalar) -> ExactResult<Self::Scalar>;
    fn multiply(&mut self, left: &Self::Scalar, right: &Self::Scalar) -> ExactResult<Self::Scalar>;
    fn add_product(
        &mut self,
        target: &Self::Scalar,
        left: &Self::Scalar,
        right: &Self::Scalar,
    ) -> ExactResult<Self::Scalar> {
        let product = self.multiply(left, right)?;
        self.add(target, &product)
    }
    fn subtract_product(
        &mut self,
        target: &Self::Scalar,
        left: &Self::Scalar,
        right: &Self::Scalar,
    ) -> ExactResult<Self::Scalar> {
        let product = self.multiply(left, right)?;
        self.subtract(target, &product)
    }
    fn inverse(&mut self, value: &Self::Scalar) -> ExactResult<Self::Scalar>;
    fn negate(&self, value: &Self::Scalar) -> Self::Scalar;
}

struct RationalDomain;

impl RawDomain for RationalDomain {
    type Scalar = Rational;

    fn zero(&self) -> Rational {
        Rational::zero()
    }
    fn one(&self) -> Rational {
        Rational::one()
    }
    fn is_zero(&self, value: &Rational) -> bool {
        value.is_zero()
    }
    fn add(&mut self, left: &Rational, right: &Rational) -> ExactResult<Rational> {
        Ok(left.add(right))
    }
    fn subtract(&mut self, left: &Rational, right: &Rational) -> ExactResult<Rational> {
        Ok(left.subtract(right))
    }
    fn multiply(&mut self, left: &Rational, right: &Rational) -> ExactResult<Rational> {
        Ok(left.multiply(right))
    }
    fn add_product(
        &mut self,
        target: &Rational,
        left: &Rational,
        right: &Rational,
    ) -> ExactResult<Rational> {
        let product_denominator = left.denominator() * right.denominator();
        Rational::new(
            target.numerator() * &product_denominator
                + left.numerator() * right.numerator() * target.denominator(),
            target.denominator() * product_denominator,
        )
    }
    fn subtract_product(
        &mut self,
        target: &Rational,
        left: &Rational,
        right: &Rational,
    ) -> ExactResult<Rational> {
        let product_denominator = left.denominator() * right.denominator();
        Rational::new(
            target.numerator() * &product_denominator
                - left.numerator() * right.numerator() * target.denominator(),
            target.denominator() * product_denominator,
        )
    }
    fn inverse(&mut self, value: &Rational) -> ExactResult<Rational> {
        Rational::one().divide(value)
    }
    fn negate(&self, value: &Rational) -> Rational {
        value.negate()
    }
}

type InverseCache = HashMap<String, HashMap<String, Coefficient>>;

fn inverse_cache() -> &'static Mutex<InverseCache> {
    static CACHE: OnceLock<Mutex<InverseCache>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

fn coefficient_cache_key(value: &[Rational]) -> String {
    let mut key = String::new();
    for coefficient in value {
        key.push_str(&coefficient.numerator_string());
        key.push('/');
        key.push_str(&coefficient.denominator_string());
        key.push(';');
    }
    key
}

struct CoefficientDomain {
    context: Arc<CoefficientFieldContext>,
    integer_modulus: Option<Vec<BigInt>>,
}

impl CoefficientDomain {
    fn new(context: Arc<CoefficientFieldContext>) -> Self {
        let integer_modulus = context
            .defining_polynomial
            .iter()
            .take(context.degree)
            .map(|coefficient| {
                (coefficient.denominator() == &BigInt::from(1))
                    .then(|| coefficient.numerator().clone())
            })
            .collect();
        Self {
            context,
            integer_modulus,
        }
    }
}

fn bigint_absolute(value: BigInt) -> BigInt {
    if value < BigInt::zero() {
        -value
    } else {
        value
    }
}

fn bigint_gcd(mut left: BigInt, mut right: BigInt) -> BigInt {
    left = bigint_absolute(left);
    right = bigint_absolute(right);
    while !right.is_zero() {
        let remainder = &left % &right;
        left = right;
        right = remainder;
    }
    left
}

fn coefficient_common_denominator(value: &[Rational]) -> BigInt {
    let mut denominator = BigInt::from(1);
    for coefficient in value {
        let gcd = bigint_gcd(denominator.clone(), coefficient.denominator().clone());
        denominator = (denominator / gcd) * coefficient.denominator();
    }
    denominator
}

fn coefficient_integer_numerators(value: &[Rational], denominator: &BigInt) -> Vec<BigInt> {
    value
        .iter()
        .map(|coefficient| coefficient.numerator() * (denominator / coefficient.denominator()))
        .collect()
}

struct IntegerCoefficient {
    numerators: Vec<BigInt>,
    denominator: BigInt,
}

fn integer_coefficient(value: &[Rational]) -> IntegerCoefficient {
    let denominator = coefficient_common_denominator(value);
    IntegerCoefficient {
        numerators: coefficient_integer_numerators(value, &denominator),
        denominator,
    }
}

fn integer_coefficient_product(
    left: &[Rational],
    right: &[Rational],
    integer_modulus: &[BigInt],
) -> IntegerCoefficient {
    let degree = left.len();
    let left_integer = integer_coefficient(left);
    let right_integer = integer_coefficient(right);
    let mut convolution = vec![BigInt::zero(); 2 * degree - 1];
    for (left_index, left_value) in left_integer.numerators.iter().enumerate() {
        if left_value.is_zero() {
            continue;
        }
        for (right_index, right_value) in right_integer.numerators.iter().enumerate() {
            if !right_value.is_zero() {
                convolution[left_index + right_index] += left_value * right_value;
            }
        }
    }
    for exponent in (degree..convolution.len()).rev() {
        let factor = convolution[exponent].clone();
        if factor.is_zero() {
            continue;
        }
        let shift = exponent - degree;
        for index in 0..degree {
            convolution[shift + index] -= &factor * &integer_modulus[index];
        }
    }
    convolution.truncate(degree);
    IntegerCoefficient {
        numerators: convolution,
        denominator: left_integer.denominator * right_integer.denominator,
    }
}

fn rational_coefficient(value: IntegerCoefficient) -> ExactResult<Coefficient> {
    value
        .numerators
        .into_iter()
        .map(|numerator| Rational::new(numerator, value.denominator.clone()))
        .collect()
}

fn combine_integer_product(
    target: &[Rational],
    left: &[Rational],
    right: &[Rational],
    integer_modulus: &[BigInt],
    subtract_product: bool,
) -> ExactResult<Coefficient> {
    let mut target_integer = integer_coefficient(target);
    let product = integer_coefficient_product(left, right, integer_modulus);
    let gcd = bigint_gcd(
        target_integer.denominator.clone(),
        product.denominator.clone(),
    );
    let target_scale = &product.denominator / &gcd;
    let product_scale = &target_integer.denominator / &gcd;
    let denominator = &target_integer.denominator * &target_scale;
    for (target_value, product_value) in
        target_integer.numerators.iter_mut().zip(product.numerators)
    {
        let mut product_term = product_value * &product_scale;
        if subtract_product {
            product_term = -product_term;
        }
        *target_value = &*target_value * &target_scale + product_term;
    }
    target_integer.denominator = denominator;
    rational_coefficient(target_integer)
}

impl RawDomain for CoefficientDomain {
    type Scalar = Coefficient;

    fn zero(&self) -> Coefficient {
        vec![Rational::zero(); self.context.degree]
    }

    fn one(&self) -> Coefficient {
        let mut result = self.zero();
        result[0] = Rational::one();
        result
    }

    fn is_zero(&self, value: &Coefficient) -> bool {
        value.iter().all(Rational::is_zero)
    }

    fn add(&mut self, left: &Coefficient, right: &Coefficient) -> ExactResult<Coefficient> {
        Ok(left
            .iter()
            .zip(right)
            .map(|(left_value, right_value)| left_value.add(right_value))
            .collect())
    }

    fn subtract(&mut self, left: &Coefficient, right: &Coefficient) -> ExactResult<Coefficient> {
        Ok(left
            .iter()
            .zip(right)
            .map(|(left_value, right_value)| left_value.subtract(right_value))
            .collect())
    }

    fn multiply(&mut self, left: &Coefficient, right: &Coefficient) -> ExactResult<Coefficient> {
        if self.is_zero(left) || self.is_zero(right) {
            return Ok(self.zero());
        }
        let modulus = self.integer_modulus.as_ref().ok_or_else(|| {
            crate::ExactError::new(
                "MalformedSerialization",
                "coefficient field modulus is not integral",
            )
        })?;
        rational_coefficient(integer_coefficient_product(left, right, modulus))
    }

    fn add_product(
        &mut self,
        target: &Coefficient,
        left: &Coefficient,
        right: &Coefficient,
    ) -> ExactResult<Coefficient> {
        let modulus = self.integer_modulus.as_ref().ok_or_else(|| {
            crate::ExactError::new(
                "MalformedSerialization",
                "coefficient field modulus is not integral",
            )
        })?;
        combine_integer_product(target, left, right, modulus, false)
    }

    fn subtract_product(
        &mut self,
        target: &Coefficient,
        left: &Coefficient,
        right: &Coefficient,
    ) -> ExactResult<Coefficient> {
        let modulus = self.integer_modulus.as_ref().ok_or_else(|| {
            crate::ExactError::new(
                "MalformedSerialization",
                "coefficient field modulus is not integral",
            )
        })?;
        combine_integer_product(target, left, right, modulus, true)
    }

    fn inverse(&mut self, value: &Coefficient) -> ExactResult<Coefficient> {
        if self.is_zero(value) {
            return fail("DivisionByZero", "cannot invert zero coefficient value");
        }
        if value.iter().skip(1).all(Rational::is_zero) {
            let mut result = self.zero();
            result[0] = Rational::one().divide(&value[0])?;
            return Ok(result);
        }
        let key = coefficient_cache_key(value);
        if let Some(cached) = inverse_cache()
            .lock()
            .expect("inverse cache mutex poisoned")
            .get(&self.context.field_id)
            .and_then(|field| field.get(&key))
            .cloned()
        {
            return Ok(cached);
        }
        let extended = polynomial_extended_gcd(value, &self.context.defining_polynomial)?;
        if extended.gcd != vec![Rational::one()] {
            return fail(
                "SingularPolynomialElement",
                "coefficient value is not invertible",
            );
        }
        let result = polynomial_pad(&extended.left_coefficient, self.context.degree)?;
        let mut cache = inverse_cache()
            .lock()
            .expect("inverse cache mutex poisoned");
        let field = cache.entry(self.context.field_id.clone()).or_default();
        if field.len() < INVERSE_CACHE_MAXIMUM_ENTRIES {
            field.insert(key, result.clone());
        }
        Ok(result)
    }

    fn negate(&self, value: &Coefficient) -> Coefficient {
        value.iter().map(Rational::negate).collect()
    }
}

fn raw_identity<D: RawDomain>(domain: &D, dimension: usize) -> ExactResult<RawMatrix<D::Scalar>> {
    let mut data = vec![
        domain.zero();
        dimension.checked_mul(dimension).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "identity shape overflows")
        })?
    ];
    for index in 0..dimension {
        data[index * dimension + index] = domain.one();
    }
    RawMatrix::new(dimension, dimension, data)
}

fn raw_transpose<D: RawDomain>(
    _domain: &D,
    matrix: &RawMatrix<D::Scalar>,
) -> ExactResult<RawMatrix<D::Scalar>> {
    let mut data = Vec::with_capacity(matrix.data.len());
    for column in 0..matrix.columns {
        for row in 0..matrix.rows {
            data.push(matrix.at(row, column).clone());
        }
    }
    RawMatrix::new(matrix.columns, matrix.rows, data)
}

fn raw_multiply<D: RawDomain>(
    domain: &mut D,
    left: &RawMatrix<D::Scalar>,
    right: &RawMatrix<D::Scalar>,
) -> ExactResult<RawMatrix<D::Scalar>> {
    if left.columns != right.rows {
        return fail(
            "DimensionMismatch",
            "raw matrix multiplication dimensions differ",
        );
    }
    let count = left.rows.checked_mul(right.columns).ok_or_else(|| {
        crate::ExactError::new("MalformedSerialization", "raw product shape overflows")
    })?;
    let mut data = vec![domain.zero(); count];
    for row in 0..left.rows {
        for inner in 0..left.columns {
            let left_value = left.at(row, inner);
            if domain.is_zero(left_value) {
                continue;
            }
            for column in 0..right.columns {
                let right_value = right.at(inner, column);
                if domain.is_zero(right_value) {
                    continue;
                }
                let offset = row * right.columns + column;
                data[offset] = domain.add_product(&data[offset], left_value, right_value)?;
            }
        }
    }
    RawMatrix::new(left.rows, right.columns, data)
}

fn raw_zero_q<D: RawDomain>(domain: &D, matrix: &RawMatrix<D::Scalar>) -> bool {
    matrix.data.iter().all(|value| domain.is_zero(value))
}

fn prepared_unique_rows<D: RawDomain>(
    domain: &D,
    matrix: &RawMatrix<D::Scalar>,
    sort_by_sparsity: bool,
) -> Vec<Vec<D::Scalar>> {
    let mut rows = Vec::new();
    for row in 0..matrix.rows {
        let current: Vec<D::Scalar> = (0..matrix.columns)
            .map(|column| matrix.at(row, column).clone())
            .collect();
        if current.iter().all(|value| domain.is_zero(value)) || rows.contains(&current) {
            continue;
        }
        rows.push(current);
    }
    if sort_by_sparsity {
        rows.sort_by_key(|row| row.iter().filter(|value| !domain.is_zero(value)).count());
    }
    rows
}

fn raw_rref<D: RawDomain>(
    domain: &mut D,
    matrix: &RawMatrix<D::Scalar>,
    prepared_rows: bool,
) -> ExactResult<RawRrefResult<D::Scalar>> {
    let mut rows: Vec<Vec<D::Scalar>> = if prepared_rows {
        (0..matrix.rows)
            .map(|row| {
                (0..matrix.columns)
                    .map(|column| matrix.at(row, column).clone())
                    .collect()
            })
            .collect()
    } else {
        prepared_unique_rows(domain, matrix, false)
    };
    let mut pivot_row = 0usize;
    let mut pivot_columns = Vec::new();
    for column in 0..matrix.columns {
        if pivot_row >= rows.len() {
            break;
        }
        let Some(found_row) =
            (pivot_row..rows.len()).find(|&candidate| !domain.is_zero(&rows[candidate][column]))
        else {
            continue;
        };
        rows.swap(pivot_row, found_row);
        let inverse = domain.inverse(&rows[pivot_row][column].clone())?;
        for value in rows[pivot_row].iter_mut().skip(column) {
            *value = domain.multiply(&value.clone(), &inverse)?;
        }
        rows[pivot_row][column] = domain.one();
        let pivot_values = rows[pivot_row].clone();
        for (candidate, candidate_values) in rows.iter_mut().enumerate() {
            if candidate == pivot_row {
                continue;
            }
            let factor = candidate_values[column].clone();
            if domain.is_zero(&factor) {
                continue;
            }
            for (current, value) in candidate_values.iter_mut().enumerate().skip(column) {
                *value =
                    domain.subtract_product(&value.clone(), &factor, &pivot_values[current])?;
            }
            candidate_values[column] = domain.zero();
        }
        pivot_columns.push(column);
        pivot_row += 1;
    }
    let free_columns: Vec<usize> = (0..matrix.columns)
        .filter(|column| !pivot_columns.contains(column))
        .collect();
    let mut data: Vec<D::Scalar> = rows.into_iter().flatten().collect();
    data.resize(matrix.rows * matrix.columns, domain.zero());
    let reduced = RawMatrix::new(matrix.rows, matrix.columns, data)?;
    let rank = pivot_columns.len();
    let nullity = free_columns.len();
    Ok(RawRrefResult {
        reduced,
        pivot_columns,
        free_columns,
        rank,
        nullity,
    })
}

fn raw_null_space<D: RawDomain>(
    domain: &mut D,
    matrix: &RawMatrix<D::Scalar>,
    verify_residual: bool,
    prepared_rows: bool,
) -> ExactResult<RawKernelResult<D::Scalar>> {
    let rref = raw_rref(domain, matrix, prepared_rows)?;
    let mut row_data = Vec::with_capacity(rref.free_columns.len() * matrix.columns);
    for &free_column in &rref.free_columns {
        let mut vector = vec![domain.zero(); matrix.columns];
        vector[free_column] = domain.one();
        for (pivot_index, &pivot_column) in rref.pivot_columns.iter().enumerate() {
            vector[pivot_column] = domain.negate(rref.reduced.at(pivot_index, free_column));
        }
        row_data.extend(vector);
    }
    let nullspace_rows = RawMatrix::new(rref.free_columns.len(), matrix.columns, row_data)?;
    let basis_matrix = raw_transpose(domain, &nullspace_rows)?;
    if verify_residual {
        let residual = raw_multiply(domain, matrix, &basis_matrix)?;
        if !raw_zero_q(domain, &residual) {
            return fail(
                "ResidualVerificationFailed",
                "raw null-space residual is nonzero",
            );
        }
    }
    Ok(RawKernelResult {
        rref,
        nullspace_rows,
        basis_matrix,
    })
}

fn adaptive_next_dense_rows(current: usize, rank_gain: usize, potential_rank_gain: usize) -> usize {
    if potential_rank_gain == 0 {
        current
    } else if rank_gain == 0 {
        ADAPTIVE_MAXIMUM_DENSE_ROWS.min(2 * current)
    } else if 3 * rank_gain >= 2 * potential_rank_gain {
        ADAPTIVE_MINIMUM_DENSE_ROWS.max(3 * current / 4)
    } else if 4 * rank_gain <= potential_rank_gain {
        ADAPTIVE_MAXIMUM_DENSE_ROWS.min((3 * current).div_ceil(2))
    } else {
        current
    }
}

fn raw_common_kernel<D: RawDomain>(
    domain: &mut D,
    coordinate_dimension: usize,
    constraints: &[RawMatrix<D::Scalar>],
) -> ExactResult<RawCommonKernelResult<D::Scalar>> {
    let mut basis = raw_identity(domain, coordinate_dimension)?;
    let mut iteration_nullities = vec![coordinate_dimension];
    let mut identity_basis = true;
    for constraint in constraints {
        if constraint.columns != coordinate_dimension {
            return fail("DimensionMismatch", "raw constraint column count differs");
        }
        if basis.columns == 0 {
            iteration_nullities.push(0);
            continue;
        }
        let rows = prepared_unique_rows(domain, constraint, true);
        let mut row_start = 0usize;
        let mut target_dense_rows = ADAPTIVE_INITIAL_DENSE_ROWS;
        while row_start < rows.len() && basis.columns > 0 {
            let previous_columns = basis.columns;
            let remaining = rows.len() - row_start;
            let dimension_limit = ADAPTIVE_MINIMUM_DENSE_ROWS.max(2 * previous_columns);
            let max_rows = ADAPTIVE_MAXIMUM_BLOCK_ROWS
                .min(dimension_limit)
                .min(remaining);
            let effective_dense_rows = target_dense_rows.min(dimension_limit);
            let target_weight = (effective_dense_rows * constraint.columns.max(1)).max(1);
            let mut row_count = 0usize;
            let mut nonzero_weight = 0usize;
            let mut block_data = Vec::new();
            while row_count < max_rows && (row_count == 0 || nonzero_weight < target_weight) {
                let row = &rows[row_start + row_count];
                nonzero_weight += row.iter().filter(|value| !domain.is_zero(value)).count();
                block_data.extend(row.iter().cloned());
                row_count += 1;
            }
            let block = RawMatrix::new(row_count, constraint.columns, block_data)?;
            let restricted = if identity_basis {
                block
            } else {
                raw_multiply(domain, &block, &basis)?
            };
            let mut rank_gain = 0usize;
            if !raw_zero_q(domain, &restricted) {
                let kernel = raw_null_space(domain, &restricted, false, true)?;
                let updated = if identity_basis {
                    kernel.basis_matrix
                } else {
                    raw_multiply(domain, &basis, &kernel.basis_matrix)?
                };
                if updated.columns > previous_columns {
                    return fail(
                        "ResidualVerificationFailed",
                        "common-kernel nullity increased",
                    );
                }
                rank_gain = previous_columns - updated.columns;
                basis = updated;
                if rank_gain > 0 {
                    identity_basis = false;
                }
            }
            target_dense_rows = adaptive_next_dense_rows(
                target_dense_rows,
                rank_gain,
                row_count.min(previous_columns),
            );
            row_start += row_count;
        }
        iteration_nullities.push(basis.columns);
    }
    for constraint in constraints {
        let residual = raw_multiply(domain, constraint, &basis)?;
        if !raw_zero_q(domain, &residual) {
            return fail(
                "ResidualVerificationFailed",
                "raw common-kernel residual is nonzero",
            );
        }
    }
    let nullity = basis.columns;
    Ok(RawCommonKernelResult {
        basis_matrix: basis,
        iteration_nullities,
        rank: coordinate_dimension - nullity,
        nullity,
        exact_residual_verified: true,
    })
}

pub(crate) fn coefficient_field_from_cyclotomic(
    context: &Arc<CyclotomicContext>,
) -> Arc<CoefficientFieldContext> {
    Arc::new(CoefficientFieldContext {
        field_kind: CoefficientFieldKind::FullCyclotomicField,
        field_id: format!("coefficient:cyclotomic:{}", context.conductor),
        degree: context.degree,
        defining_polynomial: context.cyclotomic_polynomial.clone(),
    })
}

pub(crate) fn rational_raw_from_exact(matrix: &ExactMatrix) -> ExactResult<RationalRawMatrix> {
    if matrix.context().conductor != 1 {
        return fail("ConductorMismatch", "rational backend requires conductor 1");
    }
    RawMatrix::new(
        matrix.rows(),
        matrix.columns(),
        matrix
            .entries()
            .iter()
            .map(|entry| entry.coefficients()[0].clone())
            .collect(),
    )
}

pub(crate) fn rational_raw_to_exact(
    context: &Arc<CyclotomicContext>,
    matrix: &RationalRawMatrix,
) -> ExactResult<ExactMatrix> {
    let entries = matrix
        .data
        .iter()
        .map(|value| Element::from_polynomial(Arc::clone(context), std::slice::from_ref(value)))
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(Arc::clone(context), matrix.rows, matrix.columns, entries)
}

pub(crate) fn coefficient_raw_from_exact(
    matrix: &ExactMatrix,
) -> ExactResult<CoefficientRawMatrix> {
    RawMatrix::new(
        matrix.rows(),
        matrix.columns(),
        matrix
            .entries()
            .iter()
            .map(|entry| entry.coefficients().to_vec())
            .collect(),
    )
}

pub(crate) fn coefficient_raw_to_exact(
    context: &Arc<CyclotomicContext>,
    matrix: &CoefficientRawMatrix,
) -> ExactResult<ExactMatrix> {
    let entries = matrix
        .data
        .iter()
        .map(|value| Element::from_polynomial(Arc::clone(context), value))
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(Arc::clone(context), matrix.rows, matrix.columns, entries)
}

pub(crate) fn rational_raw_multiply(
    left: &RationalRawMatrix,
    right: &RationalRawMatrix,
) -> ExactResult<RationalRawMatrix> {
    raw_multiply(&mut RationalDomain, left, right)
}

pub(crate) fn coefficient_raw_multiply(
    context: Arc<CoefficientFieldContext>,
    left: &CoefficientRawMatrix,
    right: &CoefficientRawMatrix,
) -> ExactResult<CoefficientRawMatrix> {
    raw_multiply(&mut CoefficientDomain::new(context), left, right)
}

pub(crate) fn rational_raw_rref(
    matrix: &RationalRawMatrix,
) -> ExactResult<RawRrefResult<Rational>> {
    raw_rref(&mut RationalDomain, matrix, false)
}

pub(crate) fn coefficient_raw_rref(
    context: Arc<CoefficientFieldContext>,
    matrix: &CoefficientRawMatrix,
) -> ExactResult<RawRrefResult<Coefficient>> {
    raw_rref(&mut CoefficientDomain::new(context), matrix, false)
}

pub(crate) fn rational_raw_null_space(
    matrix: &RationalRawMatrix,
) -> ExactResult<RawKernelResult<Rational>> {
    raw_null_space(&mut RationalDomain, matrix, true, false)
}

pub(crate) fn coefficient_raw_null_space(
    context: Arc<CoefficientFieldContext>,
    matrix: &CoefficientRawMatrix,
) -> ExactResult<RawKernelResult<Coefficient>> {
    raw_null_space(&mut CoefficientDomain::new(context), matrix, true, false)
}

pub(crate) fn rational_raw_common_kernel(
    coordinate_dimension: usize,
    constraints: &[RationalRawMatrix],
) -> ExactResult<RawCommonKernelResult<Rational>> {
    raw_common_kernel(&mut RationalDomain, coordinate_dimension, constraints)
}

pub(crate) fn coefficient_raw_common_kernel(
    coordinate_dimension: usize,
    constraints: &[CoefficientRawMatrix],
    context: Arc<CoefficientFieldContext>,
) -> ExactResult<RawCommonKernelResult<Coefficient>> {
    raw_common_kernel(
        &mut CoefficientDomain::new(context),
        coordinate_dimension,
        constraints,
    )
}

fn rational_matrix_multiply(
    left: &[Rational],
    left_rows: usize,
    inner: usize,
    right: &[Rational],
    right_columns: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::zero(); left_rows * right_columns];
    for row in 0..left_rows {
        for index in 0..inner {
            if left[row * inner + index].is_zero() {
                continue;
            }
            for column in 0..right_columns {
                let product =
                    left[row * inner + index].multiply(&right[index * right_columns + column]);
                result[row * right_columns + column] =
                    result[row * right_columns + column].add(&product);
            }
        }
    }
    result
}

#[derive(Clone, Debug)]
struct RealSubfieldDescriptor {
    real_context: Arc<CoefficientFieldContext>,
    embedding: Vec<Rational>,
    coordinate_rows: Vec<usize>,
    coordinate_inverse: Vec<Rational>,
    full_degree: usize,
    real_degree: usize,
}

fn real_coordinate_system(
    embedding: &[Rational],
    full_degree: usize,
    real_degree: usize,
) -> ExactResult<(Vec<usize>, Vec<Rational>)> {
    let mut transposed = vec![Rational::zero(); real_degree * full_degree];
    for row in 0..real_degree {
        for column in 0..full_degree {
            transposed[row * full_degree + column] = embedding[column * real_degree + row].clone();
        }
    }
    let selection = raw_rref(
        &mut RationalDomain,
        &RawMatrix::new(real_degree, full_degree, transposed)?,
        false,
    )?;
    if selection.rank != real_degree {
        return fail(
            "RealSubfieldConstructionFailed",
            "cannot select coordinate rows",
        );
    }
    let coordinate_rows = selection.pivot_columns;
    let mut coordinate_matrix = vec![Rational::zero(); real_degree * real_degree];
    let mut augmented = vec![Rational::zero(); real_degree * 2 * real_degree];
    for row in 0..real_degree {
        for column in 0..real_degree {
            let value = embedding[coordinate_rows[row] * real_degree + column].clone();
            coordinate_matrix[row * real_degree + column] = value.clone();
            augmented[row * 2 * real_degree + column] = value;
        }
        augmented[row * 2 * real_degree + real_degree + row] = Rational::one();
    }
    let inverse_rref = raw_rref(
        &mut RationalDomain,
        &RawMatrix::new(real_degree, 2 * real_degree, augmented)?,
        false,
    )?;
    if inverse_rref.rank != real_degree {
        return fail(
            "RealSubfieldConstructionFailed",
            "coordinate minor is singular",
        );
    }
    let mut coordinate_inverse = vec![Rational::zero(); real_degree * real_degree];
    for row in 0..real_degree {
        for column in 0..real_degree {
            coordinate_inverse[row * real_degree + column] =
                inverse_rref.reduced.at(row, real_degree + column).clone();
        }
    }
    let inverse_check = rational_matrix_multiply(
        &coordinate_inverse,
        real_degree,
        real_degree,
        &coordinate_matrix,
        real_degree,
    );
    for row in 0..real_degree {
        for column in 0..real_degree {
            let expected = if row == column {
                Rational::one()
            } else {
                Rational::zero()
            };
            if inverse_check[row * real_degree + column] != expected {
                return fail(
                    "RealSubfieldConstructionFailed",
                    "coordinate inverse verification failed",
                );
            }
        }
    }
    Ok((coordinate_rows, coordinate_inverse))
}

fn build_real_subfield_descriptor(
    full_context: &Arc<CyclotomicContext>,
) -> ExactResult<Option<RealSubfieldDescriptor>> {
    let full_degree = full_context.degree;
    if full_context.conductor <= 2 || full_degree % 2 != 0 {
        return Ok(None);
    }
    let real_degree = full_degree / 2;
    let full_field = coefficient_field_from_cyclotomic(full_context);
    let mut full_domain = CoefficientDomain::new(full_field);
    let zeta = Element::root_of_unity(full_context, full_context.conductor, 1)?;
    let zeta_inverse = Element::root_of_unity(full_context, full_context.conductor, -1)?;
    let generator = full_domain.add(
        &zeta.coefficients().to_vec(),
        &zeta_inverse.coefficients().to_vec(),
    )?;
    let mut powers = vec![full_domain.zero(); real_degree + 1];
    powers[0] = full_domain.one();
    for exponent in 1..=real_degree {
        powers[exponent] = full_domain.multiply(&powers[exponent - 1].clone(), &generator)?;
    }
    let mut embedding = vec![Rational::zero(); full_degree * real_degree];
    for column in 0..real_degree {
        for row in 0..full_degree {
            embedding[row * real_degree + column] = powers[column][row].clone();
        }
    }
    let (coordinate_rows, coordinate_inverse) =
        real_coordinate_system(&embedding, full_degree, real_degree)?;
    let selected_power: Vec<Rational> = coordinate_rows
        .iter()
        .map(|&row| powers[real_degree][row].negate())
        .collect();
    let relation = rational_matrix_multiply(
        &coordinate_inverse,
        real_degree,
        real_degree,
        &selected_power,
        1,
    );
    let relation_check =
        rational_matrix_multiply(&embedding, full_degree, real_degree, &relation, 1);
    for row in 0..full_degree {
        if relation_check[row] != powers[real_degree][row].negate() {
            return fail(
                "RealSubfieldConstructionFailed",
                "minimal polynomial verification failed",
            );
        }
    }
    let mut real_modulus = relation;
    real_modulus.push(Rational::one());
    let real_context = Arc::new(CoefficientFieldContext {
        field_kind: CoefficientFieldKind::MaximalRealSubfield,
        field_id: format!("real_subfield:{}", full_context.conductor),
        degree: real_degree,
        defining_polynomial: real_modulus,
    });
    Ok(Some(RealSubfieldDescriptor {
        real_context,
        embedding,
        coordinate_rows,
        coordinate_inverse,
        full_degree,
        real_degree,
    }))
}

fn real_descriptor_cache() -> &'static Mutex<HashMap<usize, Option<RealSubfieldDescriptor>>> {
    static CACHE: OnceLock<Mutex<HashMap<usize, Option<RealSubfieldDescriptor>>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

fn real_subfield_descriptor(
    context: &Arc<CyclotomicContext>,
) -> ExactResult<Option<RealSubfieldDescriptor>> {
    if let Some(cached) = real_descriptor_cache()
        .lock()
        .expect("real descriptor cache mutex poisoned")
        .get(&context.conductor)
        .cloned()
    {
        return Ok(cached);
    }
    let descriptor = build_real_subfield_descriptor(context)?;
    real_descriptor_cache()
        .lock()
        .expect("real descriptor cache mutex poisoned")
        .insert(context.conductor, descriptor.clone());
    Ok(descriptor)
}

fn project_real_value(
    descriptor: &RealSubfieldDescriptor,
    value: &Coefficient,
) -> Option<Coefficient> {
    let selected: Vec<Rational> = descriptor
        .coordinate_rows
        .iter()
        .map(|&row| value[row].clone())
        .collect();
    let coordinates = rational_matrix_multiply(
        &descriptor.coordinate_inverse,
        descriptor.real_degree,
        descriptor.real_degree,
        &selected,
        1,
    );
    let lifted = rational_matrix_multiply(
        &descriptor.embedding,
        descriptor.full_degree,
        descriptor.real_degree,
        &coordinates,
        1,
    );
    (lifted == *value).then_some(coordinates)
}

pub(crate) fn coefficient_work_plan(
    full_context: &Arc<CyclotomicContext>,
    coordinate_dimension: usize,
    full_constraints: &[CoefficientRawMatrix],
) -> ExactResult<CoefficientWorkPlan> {
    let full_field = coefficient_field_from_cyclotomic(full_context);
    let mut lift_matrix = vec![Rational::zero(); full_context.degree * full_context.degree];
    for index in 0..full_context.degree {
        lift_matrix[index * full_context.degree + index] = Rational::one();
    }
    let mut plan = CoefficientWorkPlan {
        public_context: Arc::clone(full_context),
        field_context: full_field,
        coordinate_dimension,
        raw_constraints: full_constraints.to_vec(),
        lift_matrix,
        lift_rows: full_context.degree,
        lift_columns: full_context.degree,
        method: "iterative".to_owned(),
    };
    let Some(descriptor) = real_subfield_descriptor(full_context)? else {
        return Ok(plan);
    };
    let mut projected_constraints = Vec::with_capacity(full_constraints.len());
    for matrix in full_constraints {
        let mut data = Vec::with_capacity(matrix.data.len());
        for value in &matrix.data {
            let Some(projected) = project_real_value(&descriptor, value) else {
                return Ok(plan);
            };
            data.push(projected);
        }
        projected_constraints.push(RawMatrix::new(matrix.rows, matrix.columns, data)?);
    }
    plan.field_context = descriptor.real_context;
    plan.raw_constraints = projected_constraints;
    plan.lift_matrix = descriptor.embedding;
    plan.lift_rows = descriptor.full_degree;
    plan.lift_columns = descriptor.real_degree;
    "iterative_real_subfield".clone_into(&mut plan.method);
    Ok(plan)
}

pub(crate) fn coefficient_work_plan_lift_basis(
    plan: &CoefficientWorkPlan,
    basis: &CoefficientRawMatrix,
) -> ExactResult<CoefficientRawMatrix> {
    if plan.lift_rows != plan.public_context.degree {
        return fail(
            "MalformedSerialization",
            "work-plan lift rows differ from public field degree",
        );
    }
    if plan.lift_rows == plan.lift_columns {
        return Ok(basis.clone());
    }
    let data = basis
        .data
        .iter()
        .map(|value| {
            rational_matrix_multiply(
                &plan.lift_matrix,
                plan.lift_rows,
                plan.lift_columns,
                value,
                1,
            )
        })
        .collect();
    RawMatrix::new(basis.rows, basis.columns, data)
}
