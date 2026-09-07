use std::collections::BTreeMap;

use nalgebra::Matrix3;
use serde::Serialize;

use crate::{PropertiesError, PropertiesResult};

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct NamedKPoint {
    pub label: String,
    pub display_label: String,
    pub coordinates: [f64; 3],
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct BandPathSegment {
    pub start: [f64; 3],
    pub end: [f64; 3],
    pub start_label: String,
    pub end_label: String,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct StandardKPathData {
    pub bravais_type: String,
    pub bz_type: String,
    pub parameters: BTreeMap<String, f64>,
    pub gram_matrix: [[f64; 3]; 3],
    pub reciprocal_cosines: [f64; 3],
    pub points: Vec<NamedKPoint>,
    pub sequences: Vec<Vec<String>>,
    pub path: Vec<BandPathSegment>,
}

#[derive(Clone, Debug)]
struct Classification {
    family: &'static str,
    bz_type: &'static str,
    parameters: BTreeMap<String, f64>,
    gram: Matrix3<f64>,
    reciprocal_cosines: [f64; 3],
}

type Point = (&'static str, [f64; 3]);

fn invalid_lattice(detail: impl Into<String>) -> PropertiesError {
    PropertiesError::new("InvalidStandardKPathLattice", detail)
}

fn parameter_map(values: &[(&str, f64)]) -> BTreeMap<String, f64> {
    values
        .iter()
        .map(|(name, value)| ((*name).to_owned(), *value))
        .collect()
}

fn family(bz_type: &str) -> &'static str {
    match bz_type {
        "CUB" => "CubicPrimitive",
        "FCC" => "CubicFaceCentered",
        "BCC" => "CubicBodyCentered",
        "TET" => "TetragonalPrimitive",
        "BCT1" | "BCT2" => "TetragonalBodyCentered",
        "ORC" => "OrthorhombicPrimitive",
        "ORCF1" | "ORCF2" | "ORCF3" => "OrthorhombicFaceCentered",
        "ORCI" => "OrthorhombicBodyCentered",
        "ORCC" => "OrthorhombicBaseCentered",
        "HEX" => "HexagonalPrimitive",
        "RHL1" | "RHL2" => "RhombohedralPrimitive",
        "MCL" => "MonoclinicPrimitive",
        "MCLC1" | "MCLC2" | "MCLC3" | "MCLC4" | "MCLC5" => "MonoclinicBaseCentered",
        _ => "TriclinicPrimitive",
    }
}

fn classify(lattice: Matrix3<f64>, tolerance: f64) -> PropertiesResult<Classification> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(PropertiesError::new(
            "InvalidPropertiesOption",
            "Tolerance must be a finite positive real number",
        ));
    }
    if lattice.iter().any(|value| !value.is_finite()) || lattice.determinant().abs() <= tolerance {
        return Err(invalid_lattice(
            "lattice must be a finite nonsingular 3 by 3 numerical matrix",
        ));
    }
    let gram = lattice * lattice.transpose();
    let diagonal = [gram[(0, 0)], gram[(1, 1)], gram[(2, 2)]];
    let off = [gram[(0, 1)], gram[(0, 2)], gram[(1, 2)]];
    let scale = gram.iter().map(|value| value.abs()).fold(1.0, f64::max);
    let close = |left: f64, right: f64| (left - right).abs() <= tolerance * scale;
    let zero = |value: f64| close(value, 0.0);
    let equal_diagonal = close(diagonal[0], diagonal[1]) && close(diagonal[1], diagonal[2]);
    let reciprocal = lattice
        .try_inverse()
        .ok_or_else(|| invalid_lattice("lattice inverse does not exist"))?
        .transpose()
        * std::f64::consts::TAU;
    let reciprocal_gram = reciprocal * reciprocal.transpose();
    let reciprocal_cosines = [
        reciprocal_gram[(1, 2)] / (reciprocal_gram[(1, 1)] * reciprocal_gram[(2, 2)]).sqrt(),
        reciprocal_gram[(0, 2)] / (reciprocal_gram[(0, 0)] * reciprocal_gram[(2, 2)]).sqrt(),
        reciprocal_gram[(0, 1)] / (reciprocal_gram[(0, 0)] * reciprocal_gram[(1, 1)]).sqrt(),
    ];

    let (bz_type, parameters) = if off.iter().copied().all(zero) {
        let kind = if equal_diagonal {
            "CUB"
        } else if close(diagonal[0], diagonal[1]) {
            "TET"
        } else {
            "ORC"
        };
        (
            kind,
            parameter_map(&[
                ("A", diagonal[0].sqrt()),
                ("B", diagonal[1].sqrt()),
                ("C", diagonal[2].sqrt()),
            ]),
        )
    } else if equal_diagonal && off.iter().all(|value| close(*value, diagonal[0] / 2.0)) {
        ("FCC", parameter_map(&[("A", (2.0 * diagonal[0]).sqrt())]))
    } else if equal_diagonal && off.iter().all(|value| close(*value, -diagonal[0] / 3.0)) {
        (
            "BCC",
            parameter_map(&[("A", (4.0 * diagonal[0] / 3.0).sqrt())]),
        )
    } else if close(diagonal[0], diagonal[1])
        && close(off[0], -diagonal[0] / 2.0)
        && zero(off[1])
        && zero(off[2])
    {
        (
            "HEX",
            parameter_map(&[("A", diagonal[0].sqrt()), ("C", diagonal[2].sqrt())]),
        )
    } else if equal_diagonal && close(off[0], off[1]) && close(off[1], off[2]) {
        let cosine = off[0] / diagonal[0];
        (
            if cosine > 0.0 { "RHL1" } else { "RHL2" },
            parameter_map(&[
                ("A", diagonal[0].sqrt()),
                ("CosAlpha", cosine),
                ("Alpha", cosine.acos()),
            ]),
        )
    } else if equal_diagonal
        && close(off[1], off[2])
        && !close(off[0], off[1])
        && off[1] < -tolerance * scale
    {
        let squared_c = -4.0 * off[1];
        let squared_a = squared_c / 2.0 - 2.0 * off[0];
        if squared_a <= 0.0 || squared_c <= 0.0 {
            return Err(invalid_lattice("invalid body-centered tetragonal metric"));
        }
        let a = squared_a.sqrt();
        let c = squared_c.sqrt();
        (
            if c < a { "BCT1" } else { "BCT2" },
            parameter_map(&[("A", a), ("C", c)]),
        )
    } else if equal_diagonal {
        let squared_a = -2.0 * (off[0] + off[1]);
        let squared_b = -2.0 * (off[0] + off[2]);
        let squared_c = -2.0 * (off[1] + off[2]);
        if squared_a <= 0.0 || squared_b <= 0.0 || squared_c <= 0.0 {
            return Err(invalid_lattice("invalid body-centered orthorhombic metric"));
        }
        (
            "ORCI",
            parameter_map(&[
                ("A", squared_a.sqrt()),
                ("B", squared_b.sqrt()),
                ("C", squared_c.sqrt()),
            ]),
        )
    } else if off.iter().all(|value| *value > tolerance * scale)
        && close(diagonal[0], off[0] + off[1])
        && close(diagonal[1], off[0] + off[2])
        && close(diagonal[2], off[1] + off[2])
    {
        let squared_a = 4.0 * off[2];
        let squared_b = 4.0 * off[1];
        let squared_c = 4.0 * off[0];
        let criterion = 1.0 - squared_a / squared_b - squared_a / squared_c;
        let kind = if criterion.abs() <= tolerance {
            "ORCF3"
        } else if criterion > 0.0 {
            "ORCF1"
        } else {
            "ORCF2"
        };
        (
            kind,
            parameter_map(&[
                ("A", squared_a.sqrt()),
                ("B", squared_b.sqrt()),
                ("C", squared_c.sqrt()),
            ]),
        )
    } else if close(diagonal[0], diagonal[1]) && close(off[1], off[2]) {
        if zero(off[1]) {
            let squared_a = 2.0 * (diagonal[0] + off[0]);
            let squared_b = 2.0 * (diagonal[0] - off[0]);
            let squared_c = diagonal[2];
            if squared_a <= 0.0 || squared_b <= 0.0 || squared_c <= 0.0 {
                return Err(invalid_lattice("invalid base-centered orthorhombic metric"));
            }
            (
                "ORCC",
                parameter_map(&[
                    ("A", squared_a.sqrt()),
                    ("B", squared_b.sqrt()),
                    ("C", squared_c.sqrt()),
                ]),
            )
        } else {
            let squared_a = 2.0 * (diagonal[0] - off[0]);
            let squared_b = 2.0 * (diagonal[0] + off[0]);
            let squared_c = diagonal[2];
            if squared_a <= 0.0 || squared_b <= 0.0 || squared_c <= 0.0 {
                return Err(invalid_lattice("invalid base-centered monoclinic metric"));
            }
            let a = squared_a.sqrt();
            let b = squared_b.sqrt();
            let c = squared_c.sqrt();
            let cosine = 2.0 * off[1] / (b * c);
            let sine = (1.0 - cosine * cosine).max(0.0).sqrt();
            let cosine_gamma = reciprocal_cosines[2];
            let criterion = b * cosine / c + b * b * sine * sine / (a * a);
            let kind = if cosine_gamma < -tolerance {
                "MCLC1"
            } else if cosine_gamma.abs() <= tolerance {
                "MCLC2"
            } else if (criterion - 1.0).abs() <= tolerance {
                "MCLC4"
            } else if criterion < 1.0 {
                "MCLC3"
            } else {
                "MCLC5"
            };
            (
                kind,
                parameter_map(&[
                    ("A", a),
                    ("B", b),
                    ("C", c),
                    ("CosAlpha", cosine),
                    ("SinAlpha", sine),
                    ("Alpha", cosine.acos()),
                ]),
            )
        }
    } else if zero(off[0]) && zero(off[1]) && !zero(off[2]) {
        let a = diagonal[0].sqrt();
        let b = diagonal[1].sqrt();
        let c = diagonal[2].sqrt();
        let cosine = off[2] / (b * c);
        (
            "MCL",
            parameter_map(&[
                ("A", a),
                ("B", b),
                ("C", c),
                ("CosAlpha", cosine),
                ("SinAlpha", (1.0 - cosine * cosine).max(0.0).sqrt()),
                ("Alpha", cosine.acos()),
            ]),
        )
    } else {
        let kind = if reciprocal_cosines.iter().all(|value| *value < -tolerance)
            && reciprocal_cosines[2] >= reciprocal_cosines[0].max(reciprocal_cosines[1]) - tolerance
        {
            "TRI1a"
        } else if reciprocal_cosines.iter().all(|value| *value > tolerance)
            && reciprocal_cosines[2] <= reciprocal_cosines[0].min(reciprocal_cosines[1]) + tolerance
        {
            "TRI1b"
        } else if reciprocal_cosines[0] < -tolerance
            && reciprocal_cosines[1] < -tolerance
            && reciprocal_cosines[2].abs() <= tolerance
        {
            "TRI2a"
        } else if reciprocal_cosines[0] > tolerance
            && reciprocal_cosines[1] > tolerance
            && reciprocal_cosines[2].abs() <= tolerance
        {
            "TRI2b"
        } else {
            "TRI"
        };
        (
            kind,
            parameter_map(&[
                ("ReciprocalCosineAlpha", reciprocal_cosines[0]),
                ("ReciprocalCosineBeta", reciprocal_cosines[1]),
                ("ReciprocalCosineGamma", reciprocal_cosines[2]),
            ]),
        )
    };

    Ok(Classification {
        family: family(bz_type),
        bz_type,
        parameters,
        gram,
        reciprocal_cosines,
    })
}

fn get(parameters: &BTreeMap<String, f64>, name: &str) -> f64 {
    parameters[name]
}

fn sequences(values: &[&[&str]]) -> Vec<Vec<String>> {
    values
        .iter()
        .map(|sequence| sequence.iter().map(|label| (*label).to_owned()).collect())
        .collect()
}

#[allow(clippy::many_single_char_names)]
fn definition(kind: &str, p: &BTreeMap<String, f64>) -> Option<(Vec<Point>, Vec<Vec<String>>)> {
    let g = "Γ";
    let result = match kind {
        "CUB" => (
            vec![
                (g, [0., 0., 0.]),
                ("X", [0., 0.5, 0.]),
                ("M", [0.5, 0.5, 0.]),
                ("R", [0.5, 0.5, 0.5]),
            ],
            sequences(&[&[g, "X", "M", g, "R", "X"], &["M", "R"]]),
        ),
        "FCC" => (
            vec![
                (g, [0., 0., 0.]),
                ("K", [0.375, 0.375, 0.75]),
                ("L", [0.5, 0.5, 0.5]),
                ("U", [0.625, 0.25, 0.625]),
                ("W", [0.5, 0.25, 0.75]),
                ("X", [0.5, 0., 0.5]),
            ],
            sequences(&[&[g, "X", "W", "K", g, "L", "U", "W", "L", "K"], &["U", "X"]]),
        ),
        "BCC" => (
            vec![
                (g, [0., 0., 0.]),
                ("H", [0.5, -0.5, 0.5]),
                ("P", [0.25, 0.25, 0.25]),
                ("N", [0., 0., 0.5]),
            ],
            sequences(&[&[g, "H", "N", g, "P", "H"], &["P", "N"]]),
        ),
        "TET" => (
            vec![
                (g, [0., 0., 0.]),
                ("A", [0.5, 0.5, 0.5]),
                ("M", [0.5, 0.5, 0.]),
                ("R", [0., 0.5, 0.5]),
                ("X", [0., 0.5, 0.]),
                ("Z", [0., 0., 0.5]),
            ],
            sequences(&[
                &[g, "X", "M", g, "Z", "R", "A", "Z"],
                &["X", "R"],
                &["M", "A"],
            ]),
        ),
        "BCT1" => {
            let a = get(p, "A");
            let c = get(p, "C");
            let eta = (1. + c * c / (a * a)) / 4.;
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("M", [-0.5, 0.5, 0.5]),
                    ("N", [0., 0.5, 0.]),
                    ("P", [0.25, 0.25, 0.25]),
                    ("X", [0., 0., 0.5]),
                    ("Z", [eta, eta, -eta]),
                    ("Z1", [-eta, 1. - eta, eta]),
                ],
                sequences(&[&[g, "X", "M", g, "Z", "P", "N", "Z1", "M"], &["X", "P"]]),
            )
        }
        "BCT2" => {
            let a = get(p, "A");
            let c = get(p, "C");
            let eta = (1. + a * a / (c * c)) / 4.;
            let zeta = a * a / (2. * c * c);
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("N", [0., 0.5, 0.]),
                    ("P", [0.25, 0.25, 0.25]),
                    ("Σ", [-eta, eta, eta]),
                    ("Σ1", [eta, 1. - eta, -eta]),
                    ("X", [0., 0., 0.5]),
                    ("Y", [-zeta, zeta, 0.5]),
                    ("Y1", [0.5, 0.5, -zeta]),
                    ("Z", [0.5, 0.5, -0.5]),
                ],
                sequences(&[
                    &[g, "X", "Y", "Σ", g, "Z", "Σ1", "N", "P", "Y1", "Z"],
                    &["X", "P"],
                ]),
            )
        }
        "ORC" => (
            vec![
                (g, [0., 0., 0.]),
                ("R", [0.5, 0.5, 0.5]),
                ("S", [0.5, 0.5, 0.]),
                ("T", [0., 0.5, 0.5]),
                ("U", [0.5, 0., 0.5]),
                ("X", [0.5, 0., 0.]),
                ("Y", [0., 0.5, 0.]),
                ("Z", [0., 0., 0.5]),
            ],
            sequences(&[
                &[g, "X", "S", "Y", g, "Z", "U", "R", "T", "Z"],
                &["Y", "T"],
                &["U", "X"],
                &["S", "R"],
            ]),
        ),
        "ORCF1" | "ORCF3" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let c = get(p, "C");
            let zeta = (1. + a * a / (b * b) - a * a / (c * c)) / 4.;
            let eta = (1. + a * a / (b * b) + a * a / (c * c)) / 4.;
            let seqs = if kind == "ORCF1" {
                sequences(&[
                    &[g, "Y", "T", "Z", g, "X", "A1", "Y"],
                    &["T", "X1"],
                    &["X", "A", "Z"],
                    &["L", g],
                ])
            } else {
                sequences(&[
                    &[g, "Y", "T", "Z", g, "X", "A1", "Y"],
                    &["X", "A", "Z"],
                    &["L", g],
                ])
            };
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("A", [0.5, 0.5 + zeta, zeta]),
                    ("A1", [0.5, 0.5 - zeta, 1. - zeta]),
                    ("L", [0.5, 0.5, 0.5]),
                    ("T", [1., 0.5, 0.5]),
                    ("X", [0., eta, eta]),
                    ("X1", [1., 1. - eta, 1. - eta]),
                    ("Y", [0.5, 0., 0.5]),
                    ("Z", [0.5, 0.5, 0.]),
                ],
                seqs,
            )
        }
        "ORCF2" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let c = get(p, "C");
            let eta = (1. + a * a / (b * b) - a * a / (c * c)) / 4.;
            let delta = (1. + b * b / (a * a) - b * b / (c * c)) / 4.;
            let phi = (1. + c * c / (b * b) - c * c / (a * a)) / 4.;
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("C", [0.5, 0.5 - eta, 1. - eta]),
                    ("C1", [0.5, 0.5 + eta, eta]),
                    ("D", [0.5 - delta, 0.5, 1. - delta]),
                    ("D1", [0.5 + delta, 0.5, delta]),
                    ("L", [0.5, 0.5, 0.5]),
                    ("H", [1. - phi, 0.5 - phi, 0.5]),
                    ("H1", [phi, 0.5 + phi, 0.5]),
                    ("X", [0., 0.5, 0.5]),
                    ("Y", [0.5, 0., 0.5]),
                    ("Z", [0.5, 0.5, 0.]),
                ],
                sequences(&[
                    &[g, "Y", "C", "D", "X", g, "Z", "D1", "H", "C"],
                    &["C1", "Z"],
                    &["X", "H1"],
                    &["H", "Y"],
                    &["L", g],
                ]),
            )
        }
        "ORCI" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let c = get(p, "C");
            let zeta = (1. + a * a / (c * c)) / 4.;
            let delta = (b * b - a * a) / (4. * c * c);
            let eta = (1. + b * b / (c * c)) / 4.;
            let mu = (a * a + b * b) / (4. * c * c);
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("L", [-mu, mu, 0.5 - delta]),
                    ("L1", [mu, -mu, 0.5 + delta]),
                    ("L2", [0.5 - delta, 0.5 + delta, -mu]),
                    ("R", [0., 0.5, 0.]),
                    ("S", [0.5, 0., 0.]),
                    ("T", [0., 0., 0.5]),
                    ("W", [0.25, 0.25, 0.25]),
                    ("X", [-zeta, zeta, zeta]),
                    ("X1", [zeta, 1. - zeta, -zeta]),
                    ("Y", [eta, -eta, eta]),
                    ("Y1", [1. - eta, eta, -eta]),
                    ("Z", [0.5, 0.5, -0.5]),
                ],
                sequences(&[
                    &[g, "X", "L", "T", "W", "R", "X1", "Z", g, "Y", "S", "W"],
                    &["L1", "Y"],
                    &["Y1", "Z"],
                ]),
            )
        }
        "ORCC" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let zeta = (1. + a * a / (b * b)) / 4.;
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("A", [zeta, zeta, 0.5]),
                    ("A1", [-zeta, 1. - zeta, 0.5]),
                    ("R", [0., 0.5, 0.5]),
                    ("S", [0., 0.5, 0.]),
                    ("T", [-0.5, 0.5, 0.5]),
                    ("X", [zeta, zeta, 0.]),
                    ("X1", [-zeta, 1. - zeta, 0.]),
                    ("Y", [-0.5, 0.5, 0.]),
                    ("Z", [0., 0., 0.5]),
                ],
                sequences(&[
                    &[g, "X", "S", "R", "A", "Z", g, "Y", "X1", "A1", "T", "Y"],
                    &["Z", "T"],
                ]),
            )
        }
        "HEX" => (
            vec![
                (g, [0., 0., 0.]),
                ("A", [0., 0., 0.5]),
                ("H", [1. / 3., 1. / 3., 0.5]),
                ("K", [1. / 3., 1. / 3., 0.]),
                ("L", [0.5, 0., 0.5]),
                ("M", [0.5, 0., 0.]),
            ],
            sequences(&[
                &[g, "M", "K", g, "A", "L", "H", "A"],
                &["L", "M"],
                &["K", "H"],
            ]),
        ),
        "RHL1" => {
            let cosine = get(p, "CosAlpha");
            let eta = (1. + 4. * cosine) / (2. + 4. * cosine);
            let nu = 0.75 - eta / 2.;
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("B", [eta, 0.5, 1. - eta]),
                    ("B1", [0.5, 1. - eta, eta - 1.]),
                    ("F", [0.5, 0.5, 0.]),
                    ("L", [0.5, 0., 0.]),
                    ("L1", [0., 0., -0.5]),
                    ("P", [eta, nu, nu]),
                    ("P1", [1. - nu, 1. - nu, 1. - eta]),
                    ("P2", [nu, nu, eta - 1.]),
                    ("Q", [1. - nu, nu, 0.]),
                    ("X", [nu, 0., -nu]),
                    ("Z", [0.5, 0.5, 0.5]),
                ],
                sequences(&[
                    &[g, "L", "B1"],
                    &["B", "Z", g, "X"],
                    &["Q", "F", "P1", "Z"],
                    &["L", "P"],
                ]),
            )
        }
        "RHL2" => {
            let alpha = get(p, "Alpha");
            let eta = 1. / (2. * (alpha / 2.).tan().powi(2));
            let nu = 0.75 - eta / 2.;
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("F", [0.5, -0.5, 0.]),
                    ("L", [0.5, 0., 0.]),
                    ("P", [1. - nu, -nu, 1. - nu]),
                    ("P1", [nu, nu - 1., nu - 1.]),
                    ("Q", [eta, eta, eta]),
                    ("Q1", [1. - eta, -eta, -eta]),
                    ("Z", [0.5, -0.5, 0.5]),
                ],
                sequences(&[&[g, "P", "Z", "Q", g, "F", "P1", "Q1", "L", "Z"]]),
            )
        }
        "MCL" => {
            let b = get(p, "B");
            let c = get(p, "C");
            let cosine = get(p, "CosAlpha");
            let sine = get(p, "SinAlpha");
            let eta = (1. - b * cosine / c) / (2. * sine * sine);
            let nu = 0.5 - eta * c * cosine / b;
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("A", [0.5, 0.5, 0.]),
                    ("C", [0., 0.5, 0.5]),
                    ("D", [0.5, 0., 0.5]),
                    ("D1", [0.5, 0., -0.5]),
                    ("E", [0.5, 0.5, 0.5]),
                    ("H", [0., eta, 1. - nu]),
                    ("H1", [0., 1. - eta, nu]),
                    ("H2", [0., eta, -nu]),
                    ("M", [0.5, eta, 1. - nu]),
                    ("M1", [0.5, 1. - eta, nu]),
                    ("M2", [0.5, eta, -nu]),
                    ("X", [0., 0.5, 0.]),
                    ("Y", [0., 0., 0.5]),
                    ("Y1", [0., 0., -0.5]),
                    ("Z", [0.5, 0., 0.]),
                ],
                sequences(&[
                    &[g, "Y", "H", "C", "E", "M1", "A", "X", "H1"],
                    &["M", "D", "Z"],
                    &["Y", "D"],
                ]),
            )
        }
        "MCLC1" | "MCLC2" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let c = get(p, "C");
            let cosine = get(p, "CosAlpha");
            let sine = get(p, "SinAlpha");
            let zeta = (2. - b * cosine / c) / (4. * sine * sine);
            let eta = 0.5 + 2. * zeta * c * cosine / b;
            let psi = 0.75 - a * a / (4. * b * b * sine * sine);
            let phi = psi + (0.75 - psi) * b * cosine / c;
            let seqs = if kind == "MCLC1" {
                sequences(&[
                    &[g, "Y", "F", "L", "I"],
                    &["I1", "Z", "F1"],
                    &["Y", "X1"],
                    &["X", g, "N"],
                    &["M", g],
                ])
            } else {
                sequences(&[&[g, "Y", "F", "L", "I"], &["I1", "Z", "F1"], &["N", g, "M"]])
            };
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("N", [0.5, 0., 0.]),
                    ("N1", [0., -0.5, 0.]),
                    ("F", [1. - zeta, 1. - zeta, 1. - eta]),
                    ("F1", [zeta, zeta, eta]),
                    ("F2", [-zeta, -zeta, 1. - eta]),
                    ("F3", [1. - zeta, -zeta, 1. - eta]),
                    ("I", [phi, 1. - phi, 0.5]),
                    ("I1", [1. - phi, phi - 1., 0.5]),
                    ("L", [0.5, 0.5, 0.5]),
                    ("M", [0.5, 0., 0.5]),
                    ("X", [1. - psi, psi - 1., 0.]),
                    ("X1", [psi, 1. - psi, 0.]),
                    ("X2", [psi - 1., -psi, 0.]),
                    ("Y", [0.5, 0.5, 0.]),
                    ("Y1", [-0.5, -0.5, 0.]),
                    ("Z", [0., 0., 0.5]),
                ],
                seqs,
            )
        }
        "MCLC3" | "MCLC4" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let c = get(p, "C");
            let cosine = get(p, "CosAlpha");
            let sine = get(p, "SinAlpha");
            let mu = (1. + b * b / (a * a)) / 4.;
            let delta = b * c * cosine / (2. * a * a);
            let zeta = mu - 0.25 + (1. - b * cosine / c) / (4. * sine * sine);
            let eta = 0.5 + 2. * zeta * c * cosine / b;
            let phi = 1. + zeta - 2. * mu;
            let psi = eta - 2. * delta;
            let seqs = if kind == "MCLC3" {
                sequences(&[
                    &[g, "Y", "F", "H", "Z", "I", "F1"],
                    &["H1", "Y1", "X", g, "N"],
                    &["M", g],
                ])
            } else {
                sequences(&[
                    &[g, "Y", "F", "H", "Z", "I"],
                    &["H1", "Y1", "X", g, "N"],
                    &["M", g],
                ])
            };
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("F", [1. - phi, 1. - phi, 1. - psi]),
                    ("F1", [phi, phi, psi]),
                    ("F2", [1. - phi, -phi, 1. - psi]),
                    ("H", [zeta, zeta, eta]),
                    ("H1", [1. - zeta, -zeta, 1. - eta]),
                    ("H2", [-zeta, -zeta, 1. - eta]),
                    ("I", [0.5, -0.5, 0.5]),
                    ("M", [0.5, 0., 0.5]),
                    ("N", [0.5, 0., 0.]),
                    ("N1", [0., -0.5, 0.]),
                    ("X", [0.5, -0.5, 0.]),
                    ("Y", [mu, mu, delta]),
                    ("Y1", [1. - mu, -mu, -delta]),
                    ("Y2", [-mu, -mu, -delta]),
                    ("Y3", [mu, mu - 1., delta]),
                    ("Z", [0., 0., 0.5]),
                ],
                seqs,
            )
        }
        "MCLC5" => {
            let a = get(p, "A");
            let b = get(p, "B");
            let c = get(p, "C");
            let cosine = get(p, "CosAlpha");
            let sine = get(p, "SinAlpha");
            let zeta = (b * b / (a * a) + (1. - b * cosine / c) / (sine * sine)) / 4.;
            let eta = 0.5 + 2. * zeta * c * cosine / b;
            let mu = eta / 2. + b * b / (4. * a * a) - b * c * cosine / (2. * a * a);
            let nu = 2. * mu - zeta;
            let omega = (4. * nu - 1. - b * b * sine * sine / (a * a)) * c / (2. * b * cosine);
            let delta = zeta * c * cosine / b + omega / 2. - 0.25;
            let rho = 1. - zeta * a * a / (b * b);
            (
                vec![
                    (g, [0., 0., 0.]),
                    ("F", [nu, nu, omega]),
                    ("F1", [1. - nu, 1. - nu, 1. - omega]),
                    ("F2", [nu, nu - 1., omega]),
                    ("H", [zeta, zeta, eta]),
                    ("H1", [1. - zeta, -zeta, 1. - eta]),
                    ("H2", [-zeta, -zeta, 1. - eta]),
                    ("I", [rho, 1. - rho, 0.5]),
                    ("I1", [1. - rho, rho - 1., 0.5]),
                    ("L", [0.5, 0.5, 0.5]),
                    ("M", [0.5, 0., 0.5]),
                    ("N", [0.5, 0., 0.]),
                    ("N1", [0., -0.5, 0.]),
                    ("X", [0.5, -0.5, 0.]),
                    ("Y", [mu, mu, delta]),
                    ("Y1", [1. - mu, -mu, -delta]),
                    ("Y2", [-mu, -mu, -delta]),
                    ("Y3", [mu, mu - 1., delta]),
                    ("Z", [0., 0., 0.5]),
                ],
                sequences(&[
                    &[g, "Y", "F", "L", "I"],
                    &["I1", "Z", "H", "F1"],
                    &["H1", "Y1", "X", g, "N"],
                    &["M", g],
                ]),
            )
        }
        "TRI1a" | "TRI2a" => (
            vec![
                (g, [0., 0., 0.]),
                ("L", [0.5, 0.5, 0.]),
                ("M", [0., 0.5, 0.5]),
                ("N", [0.5, 0., 0.5]),
                ("R", [0.5, 0.5, 0.5]),
                ("X", [0.5, 0., 0.]),
                ("Y", [0., 0.5, 0.]),
                ("Z", [0., 0., 0.5]),
            ],
            sequences(&[&["X", g, "Y"], &["L", g, "Z"], &["N", g, "M"], &["R", g]]),
        ),
        "TRI1b" | "TRI2b" => (
            vec![
                (g, [0., 0., 0.]),
                ("L", [0.5, -0.5, 0.]),
                ("M", [0., 0., 0.5]),
                ("N", [-0.5, -0.5, 0.5]),
                ("R", [0., -0.5, 0.5]),
                ("X", [0., -0.5, 0.]),
                ("Y", [0.5, 0., 0.]),
                ("Z", [-0.5, 0., 0.5]),
            ],
            sequences(&[&["X", g, "Y"], &["L", g, "Z"], &["N", g, "M"], &["R", g]]),
        ),
        _ => return None,
    };
    Some(result)
}

#[allow(clippy::match_same_arms)]
fn display_label(kind: &str, label: &str) -> String {
    let replacement = match (kind, label) {
        ("BCT1", "M") => "Z",
        ("BCT1", "Z") => "Λ",
        ("BCT1", "Z1") => "V",
        ("BCT2", "Y1") => "U",
        ("ORC", "T") => "U",
        ("ORC", "U") => "T",
        ("ORC", "X") => "Y",
        ("ORC", "Y") => "X",
        ("ORCF1", "A" | "A1") => "B",
        ("ORCF1", "T") => "Y",
        ("ORCF1", "X" | "X1") => "Δ",
        ("ORCF1", "Y") => "X",
        ("ORCF2", "C") => "B",
        ("ORCF2", "C1") => "D",
        ("ORCF2", "D") => "C",
        ("ORCF2", "D1") => "A",
        ("ORCF2", "H") => "G",
        ("ORCF2", "H1") => "H",
        ("ORCF2", "X") => "Y",
        ("ORCF2", "Y") => "X",
        ("ORCI", "R") => "S",
        ("ORCI", "S") => "R",
        ("ORCI", "X") => "U",
        ("ORCI", "X1") => "Δ",
        ("ORCI", "Y") => "Σ",
        ("ORCI", "Y1") => "F",
        ("ORCI", "Z") => "X",
        ("ORCC", "A1") => "E",
        ("ORCC", "X") => "Σ",
        ("ORCC", "X1") => "C",
        ("MCL", "X") => "Y",
        ("MCL", "Y" | "Y1") => "Z",
        ("MCL", "Z") => "B",
        ("MCLC1" | "MCLC2", "N") => "A",
        ("MCLC1" | "MCLC2", "N1") => "V",
        ("MCLC1" | "MCLC2", "L") => "M",
        ("MCLC1" | "MCLC2", "M") => "L",
        ("MCLC1" | "MCLC2", "Y" | "Y1") => "L",
        ("MCLC1" | "MCLC2", "Z") => "V",
        ("MCLC3" | "MCLC4", "I") => "M",
        ("MCLC3" | "MCLC4", "M") => "L",
        ("MCLC3" | "MCLC4", "N") => "A",
        ("MCLC3" | "MCLC4", "N1") => "V",
        ("MCLC3" | "MCLC4", "X") => "L",
        ("MCLC3" | "MCLC4", "Z") => "V",
        ("MCLC5", "L") => "M",
        ("MCLC5", "M") => "L",
        ("MCLC5", "N") => "A",
        ("MCLC5", "N1") => "V",
        ("MCLC5", "X") => "L",
        ("MCLC5", "Z") => "V",
        ("TRI1a" | "TRI2a", "X") => "B",
        ("TRI1a" | "TRI2a", "Y") => "F",
        ("TRI1a" | "TRI2a", "Z") => "G",
        ("TRI1b" | "TRI2b", "M") => "G",
        ("TRI1b" | "TRI2b", "X") => "F",
        ("TRI1b" | "TRI2b", "Y") => "B",
        _ => label,
    };
    replacement.to_owned()
}

#[allow(clippy::too_many_lines)]
pub fn standard_k_path(
    lattice: [[f64; 3]; 3],
    requested: Option<&str>,
    tolerance: f64,
) -> PropertiesResult<StandardKPathData> {
    let matrix = Matrix3::from_row_slice(&lattice.concat());
    let classification = classify(matrix, tolerance)?;
    if let Some(requested) = requested
        && requested != classification.bz_type
        && requested != classification.family
    {
        return Err(PropertiesError::new(
            "UnsupportedStandardBravaisType",
            format!(
                "requested {requested}, but the lattice classifies as {} ({})",
                classification.bz_type, classification.family
            ),
        ));
    }
    let (raw_points, sequences) = definition(classification.bz_type, &classification.parameters)
        .ok_or_else(|| {
            PropertiesError::new(
                "UnsupportedStandardBravaisType",
                format!(
                    "Bravais/BZ type {} is not supported",
                    classification.bz_type
                ),
            )
        })?;
    let points = raw_points
        .iter()
        .map(|(label, coordinates)| NamedKPoint {
            label: (*label).to_owned(),
            display_label: display_label(classification.bz_type, label),
            coordinates: *coordinates,
        })
        .collect::<Vec<_>>();
    let mut path = Vec::new();
    for pair in sequences.iter().flat_map(|sequence| sequence.windows(2)) {
        let start = points
            .iter()
            .find(|point| point.label == pair[0])
            .ok_or_else(|| {
                PropertiesError::new(
                    "InvalidStandardKPathDefinition",
                    format!("path point {} is absent from the catalogue", pair[0]),
                )
            })?;
        let end = points
            .iter()
            .find(|point| point.label == pair[1])
            .ok_or_else(|| {
                PropertiesError::new(
                    "InvalidStandardKPathDefinition",
                    format!("path point {} is absent from the catalogue", pair[1]),
                )
            })?;
        path.push(BandPathSegment {
            start: start.coordinates,
            end: end.coordinates,
            start_label: start.display_label.clone(),
            end_label: end.display_label.clone(),
        });
    }
    let gram_matrix =
        std::array::from_fn(|row| std::array::from_fn(|column| classification.gram[(row, column)]));
    Ok(StandardKPathData {
        bravais_type: classification.family.to_owned(),
        bz_type: classification.bz_type.to_owned(),
        parameters: classification.parameters,
        gram_matrix,
        reciprocal_cosines: classification.reciprocal_cosines,
        points,
        sequences,
        path,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_and_hexagonal_paths_match_stable_segment_counts() {
        let cubic = standard_k_path([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], None, 1e-6)
            .expect("cubic path");
        assert_eq!(cubic.bz_type, "CUB");
        assert_eq!(cubic.path.len(), 6);
        assert_eq!(cubic.path[0].start_label, "Γ");
        assert_eq!(cubic.path[0].end_label, "X");

        let hexagonal = standard_k_path(
            [
                [0.5, -3.0_f64.sqrt() / 2.0, 0.0],
                [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
                [0.0, 0.0, 2.0],
            ],
            None,
            1e-6,
        )
        .expect("hexagonal path");
        assert_eq!(hexagonal.bz_type, "HEX");
        assert_eq!(hexagonal.path.len(), 9);
    }
}
