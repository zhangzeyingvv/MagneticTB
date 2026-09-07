use nalgebra::{Complex, DMatrix};
use serde::Serialize;

use crate::wilson::berry_phase_from_samples;
use crate::{BerryPhaseData, PropertiesError, PropertiesResult, WilsonLoopOptions};

#[derive(Clone, Copy, Debug)]
pub struct BerryCurvatureOptions {
    pub directions: [usize; 2],
    pub step_size: [f64; 2],
    pub wilson: WilsonLoopOptions,
}

impl Default for BerryCurvatureOptions {
    fn default() -> Self {
        Self {
            directions: [0, 1],
            step_size: [1.0e-3, 1.0e-3],
            wilson: WilsonLoopOptions::default(),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct BerryCurvatureData {
    pub curvature: f64,
    pub flux: f64,
    pub area: f64,
    pub point: Vec<f64>,
    pub directions: [usize; 2],
    pub step_size: [f64; 2],
    pub plaquette_path: Vec<Vec<f64>>,
    pub phase_data: BerryPhaseData,
}

pub fn berry_phase(
    matrices: &[DMatrix<Complex<f64>>],
    path: &[Vec<f64>],
    centers: &[Vec<f64>],
    occupied: usize,
    options: WilsonLoopOptions,
) -> PropertiesResult<BerryPhaseData> {
    berry_phase_from_samples(matrices, path, centers, occupied, options)
}

pub fn berry_curvature<F>(
    mut hamiltonian: F,
    centers: &[Vec<f64>],
    occupied: usize,
    point: &[f64],
    options: BerryCurvatureOptions,
) -> PropertiesResult<BerryCurvatureData>
where
    F: FnMut(&[f64]) -> PropertiesResult<DMatrix<Complex<f64>>>,
{
    let path = berry_plaquette_path(point, options)?;
    let matrices = path
        .iter()
        .map(|coordinate| hamiltonian(coordinate))
        .collect::<PropertiesResult<Vec<_>>>()?;
    berry_curvature_from_samples(&matrices, centers, occupied, point, options)
}

pub fn berry_plaquette_path(
    point: &[f64],
    options: BerryCurvatureOptions,
) -> PropertiesResult<Vec<Vec<f64>>> {
    if point.len() < 2 || point.iter().any(|value| !value.is_finite()) {
        return Err(PropertiesError::new(
            "InvalidBerryCurvaturePoint",
            "point must be a finite real vector of dimension at least two",
        ));
    }
    if options.directions[0] == options.directions[1]
        || options.directions.iter().any(|&index| index >= point.len())
    {
        return Err(PropertiesError::new(
            "InvalidBerryCurvatureDirections",
            "directions must be distinct zero-based indices inside the point dimension",
        ));
    }
    if options
        .step_size
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(PropertiesError::new(
            "InvalidBerryCurvatureStep",
            "step sizes must be finite positive real numbers",
        ));
    }
    let mut first = vec![0.0; point.len()];
    let mut second = vec![0.0; point.len()];
    first[options.directions[0]] = options.step_size[0] / 2.0;
    second[options.directions[1]] = options.step_size[1] / 2.0;
    let signed = |first_sign: f64, second_sign: f64| {
        point
            .iter()
            .zip(&first)
            .zip(&second)
            .map(|((&value, &left), &right)| value + first_sign * left + second_sign * right)
            .collect::<Vec<_>>()
    };
    Ok(vec![
        signed(-1.0, -1.0),
        signed(1.0, -1.0),
        signed(1.0, 1.0),
        signed(-1.0, 1.0),
        signed(-1.0, -1.0),
    ])
}

pub fn berry_curvature_from_samples(
    matrices: &[DMatrix<Complex<f64>>],
    centers: &[Vec<f64>],
    occupied: usize,
    point: &[f64],
    options: BerryCurvatureOptions,
) -> PropertiesResult<BerryCurvatureData> {
    let path = berry_plaquette_path(point, options)?;
    if matrices.len() != path.len() {
        return Err(PropertiesError::new(
            "InvalidHamiltonianSamples",
            "berryCurvature requires one Hamiltonian matrix at each of five plaquette points",
        ));
    }
    let phase = berry_phase(matrices, &path, centers, occupied, options.wilson)?;
    let area = options.step_size[0] * options.step_size[1];
    Ok(BerryCurvatureData {
        curvature: phase.phase / area,
        flux: phase.phase,
        area,
        point: point.to_vec(),
        directions: options.directions,
        step_size: options.step_size,
        plaquette_path: path,
        phase_data: phase,
    })
}

#[cfg(test)]
mod tests {
    use nalgebra::{Complex, DMatrix};

    use super::{BerryCurvatureOptions, berry_curvature, berry_phase};
    use crate::WilsonLoopOptions;

    fn chern_hamiltonian(point: &[f64]) -> DMatrix<Complex<f64>> {
        let (x, y) = (point[0], point[1]);
        let dz = -1.0 + x.cos() + y.cos();
        DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex::new(dz, 0.0),
                Complex::new(x.sin(), -y.sin()),
                Complex::new(x.sin(), y.sin()),
                Complex::new(-dz, 0.0),
            ],
        )
    }

    #[test]
    fn atomic_berry_phase_matches_stable_orientation() {
        let path = (0..=40)
            .map(|index| vec![std::f64::consts::TAU * f64::from(index) / 40.0])
            .collect::<Vec<_>>();
        let matrices = path
            .iter()
            .map(|_| DMatrix::from_element(1, 1, Complex::new(0.0, 0.0)))
            .collect::<Vec<_>>();
        let data = berry_phase(
            &matrices,
            &path,
            &[vec![0.25]],
            1,
            WilsonLoopOptions::default(),
        )
        .expect("atomic Berry phase");
        assert!((data.phase - std::f64::consts::PI / 2.0).abs() <= 1.0e-10);
    }

    #[test]
    fn chern_model_curvature_at_gamma_is_one_half() {
        let data = berry_curvature(
            |point| Ok(chern_hamiltonian(point)),
            &[vec![0.0, 0.0], vec![0.0, 0.0]],
            1,
            &[0.0, 0.0],
            BerryCurvatureOptions::default(),
        )
        .expect("Berry curvature");
        assert!((data.curvature - 0.5).abs() <= 1.0e-5);
    }
}
