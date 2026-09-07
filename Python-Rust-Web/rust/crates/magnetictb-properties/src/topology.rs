use std::collections::HashMap;
use std::f64::consts::TAU;

use nalgebra::linalg::SVD;
use nalgebra::{Complex, DMatrix};

use crate::spectral::occupied_frame;
use crate::{PropertiesError, PropertiesResult};

#[derive(Clone, Copy, Debug)]
pub struct PointChernOptions {
    pub subdivisions: usize,
    pub hermitian_tolerance: f64,
    pub surface_gap_tolerance: f64,
    pub center_gap_tolerance: f64,
    pub overlap_tolerance: f64,
    pub integer_tolerance: f64,
    pub require_gapless_center: bool,
}

impl Default for PointChernOptions {
    fn default() -> Self {
        Self {
            subdivisions: 8,
            hermitian_tolerance: 1.0e-10,
            surface_gap_tolerance: 1.0e-8,
            center_gap_tolerance: 1.0e-6,
            overlap_tolerance: 1.0e-10,
            integer_tolerance: 5.0e-3,
            require_gapless_center: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CubeSurfaceMesh {
    pub points: Vec<[f64; 3]>,
    pub triangles: Vec<[usize; 3]>,
}

#[derive(Clone, Debug)]
pub struct PointChernData {
    pub chern_number: i64,
    pub raw_chern_number: f64,
    pub quantization_error: f64,
    pub point: [f64; 3],
    pub radius: f64,
    pub occupied_bands: usize,
    pub hamiltonian_dimension: usize,
    pub center_gap: f64,
    pub require_gapless_center: bool,
    pub minimum_surface_gap: f64,
    pub minimum_overlap_singular_value: f64,
    pub surface_subdivisions: usize,
    pub triangle_count: usize,
    pub unique_vertex_count: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RefinementMethod {
    PrincipalAxis,
    QuasiNewton,
}

#[derive(Clone, Debug)]
pub struct GaplessSearchOptions {
    pub brillouin_zone: Vec<[f64; 2]>,
    pub grid_size: Vec<usize>,
    pub candidate_count: usize,
    pub gap_tolerance: f64,
    pub merge_tolerance: f64,
    pub hermitian_tolerance: f64,
    pub max_iterations: usize,
    pub refinement_method: RefinementMethod,
}

impl Default for GaplessSearchOptions {
    fn default() -> Self {
        Self {
            brillouin_zone: vec![[-std::f64::consts::PI, std::f64::consts::PI]; 3],
            grid_size: vec![15; 3],
            candidate_count: 32,
            gap_tolerance: 1.0e-7,
            merge_tolerance: 1.0e-4,
            hermitian_tolerance: 1.0e-10,
            max_iterations: 500,
            refinement_method: RefinementMethod::PrincipalAxis,
        }
    }
}

#[derive(Clone, Debug)]
pub struct GaplessPointRecord {
    pub point: Vec<f64>,
    pub gap: f64,
    pub source: &'static str,
    pub grid_seed: Vec<f64>,
    pub grid_gap: f64,
    pub refined: bool,
    pub refinement_status: &'static str,
    pub objective_minimum: Option<f64>,
    pub eigenvalues: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct GaplessSearchData {
    pub points: Vec<Vec<f64>>,
    pub point_records: Vec<GaplessPointRecord>,
    pub occupied_bands: usize,
    pub hamiltonian_dimension: usize,
    pub brillouin_zone: Vec<[f64; 2]>,
    pub grid_size: Vec<usize>,
    pub grid_point_count: usize,
    pub local_minimum_count: usize,
    pub seed_count: usize,
    pub failed_refinements: usize,
    pub gap_tolerance: f64,
    pub merge_tolerance: f64,
    pub refinement_method: RefinementMethod,
    pub minimum_gap: f64,
}

fn invalid(tag: &str, detail: impl Into<String>) -> PropertiesError {
    PropertiesError::new(tag, detail)
}

fn point_key(point: [f64; 3]) -> [u64; 3] {
    point.map(f64::to_bits)
}

pub fn cube_surface_mesh(
    center: [f64; 3],
    radius: f64,
    subdivisions: usize,
) -> PropertiesResult<CubeSurfaceMesh> {
    if center.iter().any(|value| !value.is_finite())
        || !radius.is_finite()
        || radius <= 0.0
        || subdivisions < 2
    {
        return Err(invalid(
            "InvalidPointChernMesh",
            "center, positive radius, and at least two subdivisions are required",
        ));
    }
    let coordinates = (0..=subdivisions)
        .map(|index| -radius + 2.0 * radius * index as f64 / subdivisions as f64)
        .collect::<Vec<_>>();
    let faces = [
        (0usize, 1.0, 1i8),
        (0usize, -1.0, -1i8),
        (1usize, 1.0, -1i8),
        (1usize, -1.0, 1i8),
        (2usize, 1.0, 1i8),
        (2usize, -1.0, -1i8),
    ];
    let mut points = Vec::new();
    let mut point_indices = HashMap::new();
    let mut triangles = Vec::with_capacity(12 * subdivisions * subdivisions);
    let mut insert = |point: [f64; 3]| {
        let key = point_key(point);
        *point_indices.entry(key).or_insert_with(|| {
            let index = points.len();
            points.push(point);
            index
        })
    };
    for (axis, side, orientation) in faces {
        let other = (0..3).filter(|&index| index != axis).collect::<Vec<_>>();
        for first in 0..subdivisions {
            for second in 0..subdivisions {
                let make = |u: f64, v: f64| {
                    let mut point = center;
                    point[axis] += side * radius;
                    point[other[0]] += u;
                    point[other[1]] += v;
                    point
                };
                let vertices = [
                    insert(make(coordinates[first], coordinates[second])),
                    insert(make(coordinates[first + 1], coordinates[second])),
                    insert(make(coordinates[first + 1], coordinates[second + 1])),
                    insert(make(coordinates[first], coordinates[second + 1])),
                ];
                let mut first_triangle = [vertices[0], vertices[1], vertices[2]];
                let mut second_triangle = [vertices[0], vertices[2], vertices[3]];
                if orientation < 0 {
                    first_triangle.reverse();
                    second_triangle.reverse();
                }
                triangles.push(first_triangle);
                triangles.push(second_triangle);
            }
        }
    }
    Ok(CubeSurfaceMesh { points, triangles })
}

fn occupied_link(
    first: &DMatrix<Complex<f64>>,
    second: &DMatrix<Complex<f64>>,
    tolerance: f64,
) -> PropertiesResult<(Complex<f64>, f64)> {
    let overlap = first.adjoint() * second;
    let decomposition = SVD::new(overlap.clone(), false, false);
    let minimum = decomposition
        .singular_values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let determinant = overlap.determinant();
    if !minimum.is_finite() || minimum <= tolerance || determinant.norm() <= tolerance {
        return Err(invalid(
            "SingularOccupiedOverlap",
            format!("minimum occupied overlap singular value is {minimum}"),
        ));
    }
    Ok((determinant / determinant.norm(), minimum))
}

pub fn point_chern_number_from_samples(
    center_matrix: &DMatrix<Complex<f64>>,
    surface_matrices: &[DMatrix<Complex<f64>>],
    mesh: &CubeSurfaceMesh,
    center: [f64; 3],
    radius: f64,
    occupied: usize,
    options: PointChernOptions,
) -> PropertiesResult<PointChernData> {
    if surface_matrices.len() != mesh.points.len() {
        return Err(invalid(
            "InvalidHamiltonianSamples",
            "one Hamiltonian matrix is required per unique cube vertex",
        ));
    }
    for (name, value) in [
        ("HermitianTolerance", options.hermitian_tolerance),
        ("SurfaceGapTolerance", options.surface_gap_tolerance),
        ("CenterGapTolerance", options.center_gap_tolerance),
        ("OverlapTolerance", options.overlap_tolerance),
        ("IntegerTolerance", options.integer_tolerance),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(invalid(
                "InvalidPropertiesOption",
                format!("{name} must be finite and nonnegative"),
            ));
        }
    }
    if center_matrix.nrows() < 2 || occupied == 0 || occupied >= center_matrix.nrows() {
        return Err(invalid(
            "InvalidOccupiedBands",
            "pointChernNumber requires 1 through dimension-1 occupied bands",
        ));
    }
    let center_frame = occupied_frame(center_matrix, occupied, options.hermitian_tolerance)?;
    let center_gap = center_frame.eigenvalues[occupied] - center_frame.eigenvalues[occupied - 1];
    if options.require_gapless_center && center_gap > options.center_gap_tolerance {
        return Err(invalid(
            "PointChernCenterGapped",
            format!(
                "center direct gap {center_gap} exceeds tolerance {}",
                options.center_gap_tolerance
            ),
        ));
    }
    let dimension = center_matrix.nrows();
    let mut frames = Vec::with_capacity(surface_matrices.len());
    let mut minimum_gap = f64::INFINITY;
    for (index, matrix) in surface_matrices.iter().enumerate() {
        if matrix.nrows() != dimension || matrix.ncols() != dimension {
            return Err(invalid(
                "InvalidHamiltonian",
                format!("Hamiltonian dimension changed at surface vertex {index}"),
            ));
        }
        let frame = occupied_frame(matrix, occupied, options.hermitian_tolerance)?;
        let gap = frame.eigenvalues[occupied] - frame.eigenvalues[occupied - 1];
        if !gap.is_finite() || gap <= options.surface_gap_tolerance {
            return Err(invalid(
                "PointChernSurfaceGapClosed",
                format!(
                    "surface gap {gap} does not exceed tolerance {} at vertex {index}",
                    options.surface_gap_tolerance
                ),
            ));
        }
        minimum_gap = minimum_gap.min(gap);
        frames.push(frame.vectors);
    }
    let mut total_flux = 0.0;
    let mut minimum_singular = f64::INFINITY;
    for (triangle_index, triangle) in mesh.triangles.iter().enumerate() {
        let (first, first_singular) = occupied_link(
            &frames[triangle[0]],
            &frames[triangle[1]],
            options.overlap_tolerance,
        )
        .map_err(|error| {
            invalid(
                error.tag(),
                format!("triangle {triangle_index}: {}", error.detail()),
            )
        })?;
        let (second, second_singular) = occupied_link(
            &frames[triangle[1]],
            &frames[triangle[2]],
            options.overlap_tolerance,
        )?;
        let (third, third_singular) = occupied_link(
            &frames[triangle[2]],
            &frames[triangle[0]],
            options.overlap_tolerance,
        )?;
        total_flux += (first * second * third).arg();
        minimum_singular = minimum_singular
            .min(first_singular)
            .min(second_singular)
            .min(third_singular);
    }
    let raw = total_flux / TAU;
    if !raw.is_finite() || raw < i64::MIN as f64 || raw > i64::MAX as f64 {
        return Err(invalid(
            "PointChernQuantizationFailure",
            "surface flux is not finite",
        ));
    }
    let integer = raw.round() as i64;
    let error = (raw - integer as f64).abs();
    if error > options.integer_tolerance {
        return Err(invalid(
            "PointChernQuantizationFailure",
            format!(
                "raw Chern number {raw} is not within {} of an integer",
                options.integer_tolerance
            ),
        ));
    }
    Ok(PointChernData {
        chern_number: integer,
        raw_chern_number: raw,
        quantization_error: error,
        point: center,
        radius,
        occupied_bands: occupied,
        hamiltonian_dimension: dimension,
        center_gap,
        require_gapless_center: options.require_gapless_center,
        minimum_surface_gap: minimum_gap,
        minimum_overlap_singular_value: minimum_singular,
        surface_subdivisions: options.subdivisions,
        triangle_count: mesh.triangles.len(),
        unique_vertex_count: mesh.points.len(),
    })
}

fn canonical_point(point: &[f64], zone: &[[f64; 2]]) -> Vec<f64> {
    point
        .iter()
        .zip(zone)
        .map(|(&value, interval)| {
            interval[0] + (value - interval[0]).rem_euclid(interval[1] - interval[0])
        })
        .collect()
}

fn periodic_distance(first: &[f64], second: &[f64], zone: &[[f64; 2]]) -> f64 {
    first
        .iter()
        .zip(second)
        .zip(zone)
        .map(|((&left, &right), interval)| {
            let width = interval[1] - interval[0];
            let difference = (left - right).abs();
            difference.min(width - difference).powi(2)
        })
        .sum::<f64>()
        .sqrt()
}

fn cartesian_product<T: Clone>(axes: &[Vec<T>]) -> Vec<Vec<T>> {
    let mut result = vec![Vec::new()];
    for axis in axes {
        let mut next = Vec::with_capacity(result.len() * axis.len());
        for prefix in &result {
            for value in axis {
                let mut point = prefix.clone();
                point.push(value.clone());
                next.push(point);
            }
        }
        result = next;
    }
    result
}

fn flat_index(indices: &[usize], shape: &[usize]) -> usize {
    indices
        .iter()
        .zip(shape)
        .fold(0, |offset, (&index, &length)| offset * length + index)
}

fn local_minimum_indices(gaps: &[f64], shape: &[usize]) -> Vec<Vec<usize>> {
    let axes = shape
        .iter()
        .map(|&length| (0..length).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let offsets = cartesian_product(
        &(0..shape.len())
            .map(|_| vec![-1i8, 0, 1])
            .collect::<Vec<_>>(),
    );
    cartesian_product(&axes)
        .into_iter()
        .filter(|index| {
            let value = gaps[flat_index(index, shape)];
            offsets.iter().all(|offset| {
                if offset.iter().all(|&value| value == 0) {
                    return true;
                }
                let neighbor = index
                    .iter()
                    .zip(offset)
                    .zip(shape)
                    .map(|((&value, &delta), &length)| {
                        (value as isize + delta as isize).rem_euclid(length as isize) as usize
                    })
                    .collect::<Vec<_>>();
                value <= gaps[flat_index(&neighbor, shape)]
            })
        })
        .collect()
}

fn principal_axis_refinement<F>(
    start: &[f64],
    zone: &[[f64; 2]],
    grid_size: &[usize],
    max_iterations: usize,
    gap_tolerance: f64,
    objective: &mut F,
) -> PropertiesResult<(Vec<f64>, f64)>
where
    F: FnMut(&[f64]) -> PropertiesResult<f64>,
{
    let mut point = canonical_point(start, zone);
    let mut value = objective(&point)?;
    let mut steps = zone
        .iter()
        .zip(grid_size)
        .map(|(interval, &size)| (interval[1] - interval[0]) / size as f64)
        .collect::<Vec<_>>();
    for _ in 0..max_iterations {
        if value.sqrt() <= gap_tolerance {
            break;
        }
        let mut improved = false;
        for axis in 0..point.len() {
            let mut minus = point.clone();
            minus[axis] -= steps[axis];
            minus = canonical_point(&minus, zone);
            let minus_value = objective(&minus)?;
            let mut plus = point.clone();
            plus[axis] += steps[axis];
            plus = canonical_point(&plus, zone);
            let plus_value = objective(&plus)?;
            let denominator = minus_value - 2.0 * value + plus_value;
            let mut candidates = vec![(minus, minus_value), (plus, plus_value)];
            if denominator.is_finite() && denominator.abs() > f64::EPSILON {
                let offset =
                    (0.5 * (minus_value - plus_value) / denominator).clamp(-1.0, 1.0) * steps[axis];
                let mut parabolic = point.clone();
                parabolic[axis] += offset;
                parabolic = canonical_point(&parabolic, zone);
                let parabolic_value = objective(&parabolic)?;
                candidates.push((parabolic, parabolic_value));
            }
            candidates.sort_by(|left, right| left.1.total_cmp(&right.1));
            if candidates[0].1 < value {
                point = candidates.remove(0).0;
                value = objective(&point)?;
                improved = true;
            }
        }
        if !improved {
            for step in &mut steps {
                *step *= 0.5;
            }
        }
        if steps.iter().copied().fold(0.0, f64::max) <= 1.0e-12 {
            break;
        }
    }
    Ok((point, value))
}

fn finite_difference_gradient<F>(
    point: &[f64],
    zone: &[[f64; 2]],
    objective: &mut F,
) -> PropertiesResult<Vec<f64>>
where
    F: FnMut(&[f64]) -> PropertiesResult<f64>,
{
    (0..point.len())
        .map(|axis| {
            let step = 1.0e-6 * (zone[axis][1] - zone[axis][0]).max(1.0);
            let mut plus = point.to_vec();
            plus[axis] += step;
            plus = canonical_point(&plus, zone);
            let mut minus = point.to_vec();
            minus[axis] -= step;
            minus = canonical_point(&minus, zone);
            Ok((objective(&plus)? - objective(&minus)?) / (2.0 * step))
        })
        .collect()
}

fn quasi_newton_refinement<F>(
    start: &[f64],
    zone: &[[f64; 2]],
    max_iterations: usize,
    gap_tolerance: f64,
    objective: &mut F,
) -> PropertiesResult<(Vec<f64>, f64)>
where
    F: FnMut(&[f64]) -> PropertiesResult<f64>,
{
    let dimension = start.len();
    let mut point = canonical_point(start, zone);
    let mut value = objective(&point)?;
    let mut inverse_hessian = DMatrix::<f64>::identity(dimension, dimension);
    let mut gradient = finite_difference_gradient(&point, zone, objective)?;
    for _ in 0..max_iterations {
        if value.sqrt() <= gap_tolerance {
            break;
        }
        let gradient_vector = nalgebra::DVector::from_vec(gradient.clone());
        let direction = -&inverse_hessian * gradient_vector;
        let mut scale = 1.0;
        let mut candidate = point.clone();
        let mut candidate_value = value;
        for _ in 0..30 {
            candidate = canonical_point(
                &(0..dimension)
                    .map(|axis| point[axis] + scale * direction[axis])
                    .collect::<Vec<_>>(),
                zone,
            );
            candidate_value = objective(&candidate)?;
            if candidate_value < value {
                break;
            }
            scale *= 0.5;
        }
        if candidate_value >= value {
            break;
        }
        let next_gradient = finite_difference_gradient(&candidate, zone, objective)?;
        let s = nalgebra::DVector::from_iterator(
            dimension,
            candidate
                .iter()
                .zip(&point)
                .map(|(&next, &current)| next - current),
        );
        let y = nalgebra::DVector::from_iterator(
            dimension,
            next_gradient
                .iter()
                .zip(&gradient)
                .map(|(&next, &current)| next - current),
        );
        let ys = y.dot(&s);
        if ys > 1.0e-14 {
            let rho = 1.0 / ys;
            let identity = DMatrix::<f64>::identity(dimension, dimension);
            inverse_hessian = (&identity - rho * &s * y.transpose())
                * inverse_hessian
                * (&identity - rho * &y * s.transpose())
                + rho * &s * s.transpose();
        }
        point = candidate;
        value = candidate_value;
        gradient = next_gradient;
    }
    Ok((point, value))
}

pub fn find_gapless_points<F>(
    mut hamiltonian: F,
    occupied: usize,
    options: GaplessSearchOptions,
) -> PropertiesResult<GaplessSearchData>
where
    F: FnMut(&[f64]) -> PropertiesResult<DMatrix<Complex<f64>>>,
{
    let dimension = options.brillouin_zone.len();
    if !(1..=3).contains(&dimension)
        || options.grid_size.len() != dimension
        || options.brillouin_zone.iter().any(|interval| {
            !interval[0].is_finite() || !interval[1].is_finite() || interval[0] >= interval[1]
        })
        || options.grid_size.iter().any(|&size| size < 2)
        || options.candidate_count == 0
        || options.max_iterations == 0
        || [
            options.gap_tolerance,
            options.merge_tolerance,
            options.hermitian_tolerance,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(invalid(
            "InvalidGaplessSearchOption",
            "gapless-search options are invalid",
        ));
    }
    let axes = options
        .brillouin_zone
        .iter()
        .zip(&options.grid_size)
        .map(|(interval, &size)| {
            (0..size)
                .map(|index| interval[0] + (interval[1] - interval[0]) * index as f64 / size as f64)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let grid_points = cartesian_product(&axes);
    let mut cache = HashMap::<Vec<u64>, (Vec<f64>, f64)>::new();
    let mut spectrum = |point: &[f64]| -> PropertiesResult<(Vec<f64>, f64)> {
        let canonical = canonical_point(point, &options.brillouin_zone);
        let key = canonical
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>();
        if let Some(record) = cache.get(&key) {
            return Ok(record.clone());
        }
        let matrix = hamiltonian(&canonical)?;
        if matrix.nrows() < 2 || occupied == 0 || occupied >= matrix.nrows() {
            return Err(invalid(
                "InvalidOccupiedBands",
                "findGaplessPoints requires 1 through dimension-1 occupied bands",
            ));
        }
        let frame = occupied_frame(&matrix, occupied, options.hermitian_tolerance)?;
        let gap = frame.eigenvalues[occupied] - frame.eigenvalues[occupied - 1];
        let record = (frame.eigenvalues, gap);
        cache.insert(key, record.clone());
        Ok(record)
    };
    let first = spectrum(&grid_points[0])?;
    let hamiltonian_dimension = first.0.len();
    let grid_gaps = grid_points
        .iter()
        .map(|point| spectrum(point).map(|record| record.1))
        .collect::<PropertiesResult<Vec<_>>>()?;
    let local_indices = local_minimum_indices(&grid_gaps, &options.grid_size);
    let mut local_records = local_indices
        .iter()
        .map(|index| {
            let point = index
                .iter()
                .zip(&axes)
                .map(|(&coordinate, axis)| axis[coordinate])
                .collect::<Vec<_>>();
            (
                point,
                grid_gaps[flat_index(index, &options.grid_size)],
                "LocalGridMinimum",
            )
        })
        .collect::<Vec<_>>();
    local_records.sort_by(|left, right| left.1.total_cmp(&right.1));
    let mut all_grid = grid_points
        .iter()
        .cloned()
        .zip(grid_gaps.iter().copied())
        .map(|(point, gap)| (point, gap, "SmallGridGap"))
        .collect::<Vec<_>>();
    all_grid.sort_by(|left, right| left.1.total_cmp(&right.1));
    let mut seeds = local_records
        .into_iter()
        .take(options.candidate_count)
        .collect::<Vec<_>>();
    for record in all_grid {
        if seeds.len() == options.candidate_count {
            break;
        }
        if !seeds.iter().any(|existing| existing.0 == record.0) {
            seeds.push(record);
        }
    }
    let mut refinement_records = Vec::with_capacity(seeds.len());
    for (seed, grid_gap, source) in &seeds {
        if *grid_gap <= options.gap_tolerance {
            let eigenvalues = spectrum(seed)?.0;
            refinement_records.push(GaplessPointRecord {
                point: seed.clone(),
                gap: *grid_gap,
                source,
                grid_seed: seed.clone(),
                grid_gap: *grid_gap,
                refined: false,
                refinement_status: "NotNeeded",
                objective_minimum: None,
                eigenvalues,
            });
            continue;
        }
        let mut objective = |point: &[f64]| spectrum(point).map(|record| record.1.powi(2));
        let (point, minimum) = match options.refinement_method {
            RefinementMethod::PrincipalAxis => principal_axis_refinement(
                seed,
                &options.brillouin_zone,
                &options.grid_size,
                options.max_iterations,
                options.gap_tolerance,
                &mut objective,
            )?,
            RefinementMethod::QuasiNewton => quasi_newton_refinement(
                seed,
                &options.brillouin_zone,
                options.max_iterations,
                options.gap_tolerance,
                &mut objective,
            )?,
        };
        let (eigenvalues, gap) = spectrum(&point)?;
        refinement_records.push(GaplessPointRecord {
            point,
            gap,
            source,
            grid_seed: seed.clone(),
            grid_gap: *grid_gap,
            refined: true,
            refinement_status: "Succeeded",
            objective_minimum: Some(minimum),
            eigenvalues,
        });
    }
    refinement_records.sort_by(|left, right| left.gap.total_cmp(&right.gap));
    let mut final_records = Vec::<GaplessPointRecord>::new();
    for record in refinement_records
        .iter()
        .filter(|record| record.gap <= options.gap_tolerance)
    {
        if final_records.iter().all(|accepted| {
            periodic_distance(&accepted.point, &record.point, &options.brillouin_zone)
                > options.merge_tolerance
        }) {
            final_records.push(record.clone());
        }
    }
    let minimum_gap = grid_gaps
        .iter()
        .copied()
        .chain(refinement_records.iter().map(|record| record.gap))
        .fold(f64::INFINITY, f64::min);
    Ok(GaplessSearchData {
        points: final_records
            .iter()
            .map(|record| record.point.clone())
            .collect(),
        point_records: final_records,
        occupied_bands: occupied,
        hamiltonian_dimension,
        brillouin_zone: options.brillouin_zone,
        grid_size: options.grid_size,
        grid_point_count: grid_points.len(),
        local_minimum_count: local_indices.len(),
        seed_count: seeds.len(),
        failed_refinements: 0,
        gap_tolerance: options.gap_tolerance,
        merge_tolerance: options.merge_tolerance,
        refinement_method: options.refinement_method,
        minimum_gap,
    })
}

#[cfg(test)]
mod tests {
    use nalgebra::{Complex, DMatrix};

    use super::{
        GaplessSearchOptions, PointChernOptions, cube_surface_mesh, find_gapless_points,
        point_chern_number_from_samples,
    };

    fn weyl(point: [f64; 3], reversed: bool) -> DMatrix<Complex<f64>> {
        let z = if reversed { -point[2] } else { point[2] };
        DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex::new(z, 0.0),
                Complex::new(point[0], -point[1]),
                Complex::new(point[0], point[1]),
                Complex::new(-z, 0.0),
            ],
        )
    }

    #[test]
    fn single_weyl_charge_matches_stable_orientation() {
        let options = PointChernOptions {
            subdivisions: 6,
            ..PointChernOptions::default()
        };
        let mesh = cube_surface_mesh([0.0; 3], 0.2, options.subdivisions).expect("mesh");
        let matrices = mesh
            .points
            .iter()
            .map(|&point| weyl(point, false))
            .collect::<Vec<_>>();
        let data = point_chern_number_from_samples(
            &weyl([0.0; 3], false),
            &matrices,
            &mesh,
            [0.0; 3],
            0.2,
            1,
            options,
        )
        .expect("Chern number");
        assert_eq!(data.chern_number, -1);
        assert_eq!(data.triangle_count, 432);
    }

    #[test]
    fn periodic_grid_finds_stable_lattice_weyl_pair() {
        let node = std::f64::consts::PI / 3.0;
        let options = GaplessSearchOptions {
            grid_size: vec![8, 8, 12],
            candidate_count: 12,
            ..GaplessSearchOptions::default()
        };
        let result = find_gapless_points(
            |point| {
                let dz = point[2].cos() - node.cos() + 2.0 - point[0].cos() - point[1].cos();
                Ok(DMatrix::from_row_slice(
                    2,
                    2,
                    &[
                        Complex::new(dz, 0.0),
                        Complex::new(point[0].sin(), -point[1].sin()),
                        Complex::new(point[0].sin(), point[1].sin()),
                        Complex::new(-dz, 0.0),
                    ],
                ))
            },
            1,
            options,
        )
        .expect("gapless search");
        assert_eq!(result.points.len(), 2);
        let mut z = result
            .points
            .iter()
            .map(|point| point[2])
            .collect::<Vec<_>>();
        z.sort_by(f64::total_cmp);
        assert!((z[0] + node).abs() <= 1.0e-6);
        assert!((z[1] - node).abs() <= 1.0e-6);
    }
}
