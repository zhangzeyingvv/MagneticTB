fn root_sum(context: &Arc<CyclotomicContext>, order: usize) -> SymmetryResult<Element> {
    let root = Element::root_of_unity(context, order, 1).map_err(|_| {
        basis_error(
            "FieldRepresentationInsufficient",
            format!("cyclotomic field cannot represent an order-{order} root"),
        )
    })?;
    Ok(root.add(&root.conjugate()?)?)
}
fn exact_matrix_power(matrix: &ExactMatrix, exponent: usize) -> SymmetryResult<ExactMatrix> {
    let mut result = ExactMatrix::identity(matrix.context(), matrix.rows())?;
    for _ in 0..exponent {
        result = multiply(&result, matrix)?;
    }
    Ok(result)
}

fn rotation_order(rotation: &ExactMatrix) -> SymmetryResult<usize> {
    let identity = ExactMatrix::identity(rotation.context(), 3)?;
    for order in [1, 2, 3, 4, 6, 8, 12] {
        if exact_matrix_equal(&exact_matrix_power(rotation, order)?, &identity) {
            return Ok(order);
        }
    }
    Err(basis_error(
        "UnsupportedSpinRotation",
        "proper spin rotation order is not crystallographic",
    ))
}

fn norm_square(vector: &[Element]) -> SymmetryResult<Element> {
    let mut result = Element::zero(vector[0].context())?;
    for value in vector {
        result = result.add(&value.multiply(value)?)?;
    }
    Ok(result)
}

fn known_square_root(value: &Element) -> SymmetryResult<Element> {
    let context = value.context();
    let one = Element::one(context)?;
    let two = scalar(context, 2)?;
    let three = scalar(context, 3)?;
    let sqrt_two = root_sum(context, 8).ok();
    let sqrt_three = root_sum(context, 12).ok();
    let mut candidates = vec![one.clone(), two.clone(), three.clone()];
    if let Some(root) = &sqrt_two {
        candidates.extend([root.clone(), root.multiply(&two)?]);
    }
    if let Some(root) = &sqrt_three {
        candidates.extend([root.clone(), root.multiply(&two)?]);
    }
    if let (Some(left), Some(right)) = (&sqrt_two, &sqrt_three) {
        candidates.push(left.multiply(right)?);
    }
    let fractions = [rational(context, 1, 2)?, rational(context, 1, 3)?];
    let original = candidates.clone();
    for candidate in original {
        for fraction in &fractions {
            candidates.push(candidate.multiply(fraction)?);
        }
    }
    candidates
        .into_iter()
        .find(|candidate| {
            candidate
                .multiply(candidate)
                .is_ok_and(|square| &square == value)
        })
        .ok_or_else(|| {
            basis_error(
                "FieldRepresentationInsufficient",
                "exact square root required by a two-fold spin rotation is unavailable",
            )
        })
}

fn proper_rotation(spin_rotation: &ExactMatrix) -> SymmetryResult<ExactMatrix> {
    if spin_rotation.rows() != 3 || spin_rotation.columns() != 3 {
        return Err(basis_error(
            "InvalidSpinRotation",
            "spin action must be an exact 3 by 3 matrix",
        ));
    }
    let inverse = matrix_inverse(spin_rotation)?;
    if !exact_matrix_equal(&inverse, &spin_rotation.transpose()?) {
        return Err(basis_error(
            "InvalidSpinRotation",
            "spin action must be exact orthogonal",
        ));
    }
    let determinant = determinant_three(spin_rotation)?;
    if determinant == Element::one(spin_rotation.context())? {
        Ok(spin_rotation.clone())
    } else if determinant == scalar(spin_rotation.context(), -1)? {
        let minus_one = scalar(spin_rotation.context(), -1)?;
        let entries = spin_rotation
            .entries()
            .iter()
            .map(|entry| entry.multiply(&minus_one))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(ExactMatrix::new(
            Arc::clone(spin_rotation.context()),
            3,
            3,
            entries,
        )?)
    } else {
        Err(basis_error(
            "InvalidSpinRotation",
            "spin action determinant must be exactly plus or minus one",
        ))
    }
}

fn determinant_three(matrix: &ExactMatrix) -> SymmetryResult<Element> {
    let positive = matrix
        .entry(0, 0)?
        .multiply(matrix.entry(1, 1)?)?
        .multiply(matrix.entry(2, 2)?)?
        .add(
            &matrix
                .entry(0, 1)?
                .multiply(matrix.entry(1, 2)?)?
                .multiply(matrix.entry(2, 0)?)?,
        )?
        .add(
            &matrix
                .entry(0, 2)?
                .multiply(matrix.entry(1, 0)?)?
                .multiply(matrix.entry(2, 1)?)?,
        )?;
    let negative = matrix
        .entry(0, 2)?
        .multiply(matrix.entry(1, 1)?)?
        .multiply(matrix.entry(2, 0)?)?
        .add(
            &matrix
                .entry(0, 1)?
                .multiply(matrix.entry(1, 0)?)?
                .multiply(matrix.entry(2, 2)?)?,
        )?
        .add(
            &matrix
                .entry(0, 0)?
                .multiply(matrix.entry(1, 2)?)?
                .multiply(matrix.entry(2, 1)?)?,
        )?;
    Ok(positive.subtract(&negative)?)
}

#[allow(clippy::too_many_lines)]
fn spin_matrix(spin_rotation: &ExactMatrix) -> SymmetryResult<ExactMatrix> {
    let rotation = proper_rotation(spin_rotation)?;
    let context = rotation.context();
    let order = rotation_order(&rotation)?;
    let (q0, qx, qy, qz) = if order == 1 {
        let zero = Element::zero(context)?;
        (Element::one(context)?, zero.clone(), zero.clone(), zero)
    } else if order == 2 {
        let mut axis_row = None;
        let mut largest_row_key: Option<(bool, Rational)> = None;
        let identity = ExactMatrix::identity(context, 3)?;
        let sum = magnetictb_linear_algebra::add(&rotation, &identity)?;
        // Mathematica's nativeSpinMatrix first takes the numerical Max of
        // Abs within each projector row, then applies Ordering[..., -1] to
        // those exact expressions.  Ordering is canonical-expression order:
        // a radical sorts after a rational even when its numerical magnitude
        // is smaller.  Stable crystallographic rows contain rational values
        // or one canonical radical whose square is rational.  The pair below
        // reproduces that order exactly: radical before rational, then exact
        // squared magnitude, with equal keys resolved to the later row.
        for row in 0..3 {
            let vector = (0..3)
                .map(|column| sum.entry(row, column).cloned())
                .collect::<Result<Vec<_>, _>>()?;
            let mut row_score: Option<(bool, Rational)> = None;
            for entry in &vector {
                let square = entry.multiply(entry)?;
                if square
                    .coefficients()
                    .iter()
                    .skip(1)
                    .any(|coefficient| !coefficient.is_zero())
                {
                    return Err(basis_error(
                        "UnsupportedSpinRotationGauge",
                        "a two-fold spin lift requires projector entries with rational squares",
                    ));
                }
                let square = square.coefficients()[0].clone();
                let radical = entry
                    .coefficients()
                    .iter()
                    .skip(1)
                    .any(|coefficient| !coefficient.is_zero());
                if row_score
                    .as_ref()
                    .is_none_or(|(_, largest_square)| square > *largest_square)
                {
                    row_score = Some((radical, square));
                }
            }
            let Some(key) = row_score.filter(|(_, square)| !square.is_zero()) else {
                continue;
            };
            if largest_row_key
                .as_ref()
                .is_none_or(|largest| key >= *largest)
            {
                axis_row = Some(vector);
                largest_row_key = Some(key);
            }
        }
        let axis_row = axis_row.ok_or_else(|| {
            basis_error(
                "InvalidSpinRotation",
                "two-fold rotation has no invariant axis",
            )
        })?;
        let norm = known_square_root(&norm_square(&axis_row)?)?;
        (
            Element::zero(context)?,
            axis_row[0].divide(&norm)?,
            axis_row[1].divide(&norm)?,
            axis_row[2].divide(&norm)?,
        )
    } else {
        let root = Element::root_of_unity(context, 2 * order, 1).map_err(|_| {
            basis_error(
                "FieldRepresentationInsufficient",
                format!("spin lift needs an order-{} root of unity", 2 * order),
            )
        })?;
        let q0 = root
            .add(&root.conjugate()?)?
            .scale(&Rational::parse("1", "2")?)?;
        let four_q0 = q0.scale(&Rational::from_i64(4))?;
        (
            q0,
            rotation
                .entry(2, 1)?
                .subtract(rotation.entry(1, 2)?)?
                .divide(&four_q0)?,
            rotation
                .entry(0, 2)?
                .subtract(rotation.entry(2, 0)?)?
                .divide(&four_q0)?,
            rotation
                .entry(1, 0)?
                .subtract(rotation.entry(0, 1)?)?
                .divide(&four_q0)?,
        )
    };
    let imaginary = Element::root_of_unity(context, 4, 1)?;
    let iqx = imaginary.multiply(&qx)?;
    let iqz = imaginary.multiply(&qz)?;
    ExactMatrix::new(
        Arc::clone(context),
        2,
        2,
        vec![
            q0.subtract(&iqz)?,
            qy.negate()?.subtract(&iqx)?,
            qy.subtract(&iqx)?,
            q0.add(&iqz)?,
        ],
    )
    .map_err(Into::into)
}

/// Lift one fractional-coordinate spatial operation to the exact two-component
/// spin matrix used by the stable `MagneticTB` basis compiler.
pub fn compile_spatial_spin_action(
    spatial_rotation: &ExactMatrix,
    lattice: &ExactMatrix,
) -> SymmetryResult<ExactMatrix> {
    if spatial_rotation.rows() != 3
        || spatial_rotation.columns() != 3
        || lattice.rows() != 3
        || lattice.columns() != 3
        || spatial_rotation.context() != lattice.context()
    {
        return Err(basis_error(
            "InvalidSpinRotation",
            "spatial rotation and lattice must be aligned exact 3 by 3 matrices",
        ));
    }
    let lattice_transpose = lattice.transpose()?;
    let inverse_lattice_transpose = matrix_inverse(&lattice_transpose)?;
    let cartesian = multiply(
        &multiply(&lattice_transpose, spatial_rotation)?,
        &inverse_lattice_transpose,
    )?;
    let determinant = determinant_three(&cartesian)?;
    let entries = cartesian
        .entries()
        .iter()
        .map(|entry| entry.multiply(&determinant))
        .collect::<Result<Vec<_>, _>>()?;
    ExactMatrix::new(Arc::clone(lattice.context()), 3, 3, entries).map_err(Into::into)
}

/// Lift one fractional-coordinate spatial operation to the exact two-component
/// spin matrix used by the stable `MagneticTB` basis compiler.
pub fn compile_spatial_spinor_matrix(
    spatial_rotation: &ExactMatrix,
    lattice: &ExactMatrix,
) -> SymmetryResult<ExactMatrix> {
    spin_matrix(&compile_spatial_spin_action(spatial_rotation, lattice)?)
}
