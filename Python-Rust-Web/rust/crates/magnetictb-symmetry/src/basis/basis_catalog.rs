fn orbital_polynomial(
    context: &Arc<CyclotomicContext>,
    label: &str,
) -> SymmetryResult<(Polynomial, Option<usize>)> {
    let minus_one = scalar(context, -1)?;
    let two = scalar(context, 2)?;
    let imaginary = Element::root_of_unity(context, 4, 1)?;
    let minus_imaginary = imaginary.negate()?;
    let x = variable(context, 0)?;
    let y = variable(context, 1)?;
    let z = variable(context, 2)?;
    let x2 = polynomial_multiply(&x, &x)?;
    let y2 = polynomial_multiply(&y, &y)?;
    let z2 = polynomial_multiply(&z, &z)?;
    let polynomial = match label {
        "s" => constant(context, 1)?,
        "px" => x,
        "py" => y,
        "pz" => z,
        "px+ipy" => polynomial_add(&x, &polynomial_scale(&y, &imaginary)?)?,
        "px-ipy" => polynomial_add(&x, &polynomial_scale(&y, &minus_imaginary)?)?,
        "dx2-y2" => polynomial_add(&x2, &polynomial_scale(&y2, &minus_one)?)?,
        "dz2" => polynomial_add(
            &polynomial_scale(&z2, &two)?,
            &polynomial_scale(&polynomial_add(&x2, &y2)?, &minus_one)?,
        )?,
        "dxy" => polynomial_scale(&polynomial_multiply(&x, &y)?, &two)?,
        "dyz" => polynomial_scale(&polynomial_multiply(&y, &z)?, &two)?,
        "dxz" => polynomial_scale(&polynomial_multiply(&x, &z)?, &two)?,
        _ => {
            let (orbital, spin) = split_spin_label(label).ok_or_else(|| {
                basis_error(
                    "UnknownBasisFunction",
                    format!("basis label {label:?} is not in the MagneticTB catalog"),
                )
            })?;
            let (mut polynomial, _) = orbital_polynomial(context, orbital)?;
            if matches!(label, "ptest3" | "ptest4") {
                let sqrt_two = root_sum(context, 8)?;
                polynomial = polynomial_scale(&polynomial, &sqrt_two.inverse()?)?;
            }
            return Ok((polynomial, Some(spin)));
        }
    };
    Ok((polynomial, None))
}
fn split_spin_label(label: &str) -> Option<(&str, usize)> {
    match label {
        "sup" => Some(("s", 0)),
        "sdn" => Some(("s", 1)),
        "pxup" => Some(("px", 0)),
        "pxdn" => Some(("px", 1)),
        "pyup" => Some(("py", 0)),
        "pydn" => Some(("py", 1)),
        "pzup" => Some(("pz", 0)),
        "pzdn" => Some(("pz", 1)),
        "px+ipy up" => Some(("px+ipy", 0)),
        "px+ipy dn" | "ptest3" => Some(("px+ipy", 1)),
        "px-ipy up" | "ptest4" => Some(("px-ipy", 0)),
        "px-ipy dn" => Some(("px-ipy", 1)),
        "dx2-y2up" => Some(("dx2-y2", 0)),
        "dx2-y2dn" => Some(("dx2-y2", 1)),
        "dz2up" => Some(("dz2", 0)),
        "dz2dn" => Some(("dz2", 1)),
        "dxyup" => Some(("dxy", 0)),
        "dxydn" => Some(("dxy", 1)),
        "dyzup" => Some(("dyz", 0)),
        "dyzdn" => Some(("dyz", 1)),
        "dxzup" => Some(("dxz", 0)),
        "dxzdn" => Some(("dxz", 1)),
        _ => None,
    }
}

fn catalog_basis(
    context: &Arc<CyclotomicContext>,
    labels: &[String],
) -> SymmetryResult<(Vec<BasisFunction>, bool)> {
    if labels.is_empty() {
        return Err(basis_error(
            "InvalidFunctionBasis",
            "basis must contain at least one ordered function",
        ));
    }
    let resolved = labels
        .iter()
        .map(|label| orbital_polynomial(context, label))
        .collect::<SymmetryResult<Vec<_>>>()?;
    let spinor = resolved[0].1.is_some();
    if resolved.iter().any(|entry| entry.1.is_some() != spinor) {
        return Err(basis_error(
            "InvalidFunctionBasis",
            "scalar and spinor basis functions cannot be mixed",
        ));
    }
    let basis = resolved
        .into_iter()
        .map(|(polynomial, spin)| {
            if let Some(component) = spin {
                let mut value = vec![Polynomial::new(), Polynomial::new()];
                value[component] = polynomial;
                value
            } else {
                vec![polynomial]
            }
        })
        .collect();
    Ok((basis, spinor))
}
