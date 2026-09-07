use nalgebra::{Matrix3, Vector3};
use serde::Serialize;

use crate::{PropertiesError, PropertiesResult};

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct BrillouinZoneData {
    pub reciprocal_lattice: [[f64; 3]; 3],
    pub vertices: Vec<[f64; 3]>,
    pub facets: Vec<Vec<usize>>,
    pub translation_range: usize,
    pub tolerance: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct FoldedBandSegment {
    pub start: [f64; 3],
    pub end: [f64; 3],
    pub start_label: String,
    pub end_label: String,
    pub cartesian_start: [f64; 3],
    pub cartesian_end: [f64; 3],
}

fn finite_matrix(values: [[f64; 3]; 3]) -> bool {
    values.iter().flatten().all(|value| value.is_finite())
}

pub fn reciprocal_lattice(lattice: [[f64; 3]; 3]) -> PropertiesResult<[[f64; 3]; 3]> {
    if !finite_matrix(lattice) {
        return Err(PropertiesError::new(
            "InvalidBrillouinZoneInput",
            "the direct lattice must be a finite 3 by 3 matrix",
        ));
    }
    let direct = Matrix3::from_row_slice(&lattice.concat());
    let inverse = direct.try_inverse().ok_or_else(|| {
        PropertiesError::new(
            "InvalidBrillouinZoneInput",
            "the direct lattice must be nonsingular",
        )
    })?;
    let reciprocal = std::f64::consts::TAU * inverse.transpose();
    Ok(std::array::from_fn(|row| {
        std::array::from_fn(|column| reciprocal[(row, column)])
    }))
}

fn reciprocal_vectors(reciprocal: &Matrix3<f64>, range: i32) -> Vec<Vector3<f64>> {
    let mut values = Vec::new();
    for first in -range..=range {
        for second in -range..=range {
            for third in -range..=range {
                if first == 0 && second == 0 && third == 0 {
                    continue;
                }
                let fractional =
                    Vector3::new(f64::from(first), f64::from(second), f64::from(third));
                values.push(reciprocal.transpose() * fractional);
            }
        }
    }
    values.sort_by(|left, right| left.norm_squared().total_cmp(&right.norm_squared()));
    values
}

fn fractional_representatives(
    point: [f64; 3],
    reciprocal: &Matrix3<f64>,
    range: i32,
    tolerance: f64,
) -> Vec<Vector3<f64>> {
    let point = Vector3::from_row_slice(&point);
    let mut candidates = Vec::new();
    for first in -range..=range {
        for second in -range..=range {
            for third in -range..=range {
                let translation =
                    Vector3::new(f64::from(first), f64::from(second), f64::from(third));
                let candidate = point - translation;
                let norm = (reciprocal.transpose() * candidate).norm();
                candidates.push((candidate, norm));
            }
        }
    }
    let minimum = candidates
        .iter()
        .map(|(_, norm)| *norm)
        .fold(f64::INFINITY, f64::min);
    let mut representatives = candidates
        .into_iter()
        .filter_map(|(candidate, norm)| {
            (norm <= minimum + tolerance * minimum.max(1.0)).then_some(candidate)
        })
        .collect::<Vec<_>>();
    representatives.sort_by(|left, right| {
        left[0]
            .total_cmp(&right[0])
            .then_with(|| left[1].total_cmp(&right[1]))
            .then_with(|| left[2].total_cmp(&right[2]))
    });
    representatives
}

pub fn fold_band_path_to_first_bz(
    reciprocal_lattice: [[f64; 3]; 3],
    path: &[([f64; 3], [f64; 3], String, String)],
    translation_range: usize,
    tolerance: f64,
) -> PropertiesResult<Vec<FoldedBandSegment>> {
    if !finite_matrix(reciprocal_lattice)
        || translation_range == 0
        || !tolerance.is_finite()
        || tolerance <= 0.0
    {
        return Err(PropertiesError::new(
            "InvalidBrillouinZoneOption",
            "reciprocal lattice, TranslationRange, and Tolerance must be finite and valid",
        ));
    }
    let range = i32::try_from(translation_range).map_err(|_| {
        PropertiesError::new(
            "InvalidBrillouinZoneOption",
            "TranslationRange exceeds the supported integer range",
        )
    })?;
    let reciprocal = Matrix3::from_row_slice(&reciprocal_lattice.concat());
    let all_vectors = reciprocal_vectors(&reciprocal, range);
    let all_bounds = all_vectors
        .iter()
        .map(|vector| vector.norm_squared() / 2.0)
        .collect::<Vec<_>>();
    let scale_tolerance = tolerance * all_bounds.iter().copied().fold(1.0_f64, f64::max);
    path.iter()
        .map(|(start, end, start_label, end_label)| {
            let starts = fractional_representatives(*start, &reciprocal, range, tolerance);
            let ends = fractional_representatives(*end, &reciprocal, range, tolerance);
            let mut best: Option<(Vector3<f64>, Vector3<f64>, f64)> = None;
            for left in &starts {
                for right in &ends {
                    let distance = (reciprocal.transpose() * (left - right)).norm();
                    if best
                        .as_ref()
                        .is_none_or(|(_, _, current)| distance < *current)
                    {
                        best = Some((*left, *right, distance));
                    }
                }
            }
            let (left, right, _) = best.ok_or_else(|| {
                PropertiesError::new(
                    "InvalidBrillouinZonePath",
                    "a path endpoint has no representative in the reciprocal search range",
                )
            })?;
            let cartesian_left = reciprocal.transpose() * left;
            let cartesian_right = reciprocal.transpose() * right;
            for point in [&cartesian_left, &cartesian_right] {
                if all_vectors
                    .iter()
                    .zip(&all_bounds)
                    .any(|(vector, bound)| vector.dot(point) - bound > 20.0 * scale_tolerance)
                {
                    return Err(PropertiesError::new(
                        "InvalidBrillouinZonePath",
                        "a folded path endpoint lies outside the first Brillouin zone",
                    ));
                }
            }
            Ok(FoldedBandSegment {
                start: [left[0], left[1], left[2]],
                end: [right[0], right[1], right[2]],
                start_label: start_label.clone(),
                end_label: end_label.clone(),
                cartesian_start: [cartesian_left[0], cartesian_left[1], cartesian_left[2]],
                cartesian_end: [cartesian_right[0], cartesian_right[1], cartesian_right[2]],
            })
        })
        .collect()
}

fn deduplicate_vertices(values: Vec<Vector3<f64>>, tolerance: f64) -> Vec<Vector3<f64>> {
    values.into_iter().fold(Vec::new(), |mut result, value| {
        if result
            .iter()
            .all(|candidate: &Vector3<f64>| (candidate - value).norm() > tolerance)
        {
            result.push(value);
        }
        result
    })
}

fn ordered_facet(
    indices: &[usize],
    vertices: &[Vector3<f64>],
    normal: &Vector3<f64>,
) -> Vec<usize> {
    let center = indices
        .iter()
        .fold(Vector3::zeros(), |sum, &index| sum + vertices[index])
        / indices.len() as f64;
    let normal = normal.normalize();
    let reference = if normal[0].abs() <= normal[1].abs() && normal[0].abs() <= normal[2].abs() {
        Vector3::x()
    } else if normal[1].abs() <= normal[2].abs() {
        Vector3::y()
    } else {
        Vector3::z()
    };
    let first_axis = normal.cross(&reference).normalize();
    let second_axis = normal.cross(&first_axis);
    let mut ordered = indices
        .iter()
        .map(|&index| {
            let offset = vertices[index] - center;
            (
                index,
                offset.dot(&second_axis).atan2(offset.dot(&first_axis)),
            )
        })
        .collect::<Vec<_>>();
    ordered.sort_by(|left, right| left.1.total_cmp(&right.1));
    ordered.into_iter().map(|(index, _)| index).collect()
}

pub fn first_brillouin_zone(
    lattice: [[f64; 3]; 3],
    translation_range: usize,
    tolerance: f64,
) -> PropertiesResult<BrillouinZoneData> {
    if translation_range == 0 || !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(PropertiesError::new(
            "InvalidBrillouinZoneOption",
            "TranslationRange and Tolerance must be positive",
        ));
    }
    let signed_range = i32::try_from(translation_range).map_err(|_| {
        PropertiesError::new(
            "InvalidBrillouinZoneOption",
            "TranslationRange exceeds the supported integer range",
        )
    })?;
    let reciprocal_values = reciprocal_lattice(lattice)?;
    let reciprocal = Matrix3::from_row_slice(&reciprocal_values.concat());
    let all_vectors = reciprocal_vectors(&reciprocal, signed_range);
    if all_vectors.len() < 4 {
        return Err(PropertiesError::new(
            "BrillouinZoneConstructionFailed",
            "fewer than four reciprocal translations were generated",
        ));
    }
    let bounds = all_vectors
        .iter()
        .map(|vector| vector.norm_squared() / 2.0)
        .collect::<Vec<_>>();
    let scale_tolerance = tolerance * bounds.iter().copied().fold(1.0_f64, f64::max);
    let plane_count = all_vectors.len().min(24);
    let mut raw_vertices = Vec::new();
    for first in 0..plane_count.saturating_sub(2) {
        for second in first + 1..plane_count.saturating_sub(1) {
            for third in second + 1..plane_count {
                let rows = [
                    &all_vectors[first],
                    &all_vectors[second],
                    &all_vectors[third],
                ];
                let system = Matrix3::from_row_slice(&[
                    rows[0][0], rows[0][1], rows[0][2], rows[1][0], rows[1][1], rows[1][2],
                    rows[2][0], rows[2][1], rows[2][2],
                ]);
                if system.determinant().abs() <= scale_tolerance {
                    continue;
                }
                let right = Vector3::new(bounds[first], bounds[second], bounds[third]);
                let Some(point) = system.lu().solve(&right) else {
                    continue;
                };
                if point.iter().all(|value| value.is_finite())
                    && all_vectors
                        .iter()
                        .zip(&bounds)
                        .all(|(vector, bound)| vector.dot(&point) - bound <= 10.0 * scale_tolerance)
                {
                    raw_vertices.push(point);
                }
            }
        }
    }
    let vertices = deduplicate_vertices(raw_vertices, 10.0 * scale_tolerance);
    if vertices.len() < 4 {
        return Err(PropertiesError::new(
            "BrillouinZoneConstructionFailed",
            "the reciprocal half spaces produced fewer than four vertices",
        ));
    }
    let mut facet_keys = Vec::<Vec<usize>>::new();
    let mut facets = Vec::new();
    for (vector, bound) in all_vectors.iter().zip(&bounds) {
        let mut indices = vertices
            .iter()
            .enumerate()
            .filter_map(|(index, point)| {
                ((vector.dot(point) - bound).abs() <= 10.0 * scale_tolerance).then_some(index)
            })
            .collect::<Vec<_>>();
        if indices.len() < 3 {
            continue;
        }
        indices.sort_unstable();
        if facet_keys.contains(&indices) {
            continue;
        }
        facet_keys.push(indices.clone());
        facets.push(ordered_facet(&indices, &vertices, vector));
    }
    if facets.len() < 4 {
        return Err(PropertiesError::new(
            "BrillouinZoneConstructionFailed",
            "the reciprocal half spaces produced too few facets",
        ));
    }
    Ok(BrillouinZoneData {
        reciprocal_lattice: reciprocal_values,
        vertices: vertices
            .iter()
            .map(|point| [point[0], point[1], point[2]])
            .collect(),
        facets,
        translation_range,
        tolerance,
    })
}

#[cfg(test)]
mod tests {
    use super::{first_brillouin_zone, fold_band_path_to_first_bz};

    #[test]
    fn cubic_first_bz_has_stable_vertex_and_facet_counts() {
        let data =
            first_brillouin_zone([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 2, 1e-8)
                .expect("cubic BZ");
        assert_eq!(data.vertices.len(), 8);
        assert_eq!(data.facets.len(), 6);
    }

    #[test]
    fn cubic_path_is_folded_to_nearest_fractional_representatives() {
        let reciprocal =
            first_brillouin_zone([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 2, 1e-8)
                .expect("cubic BZ")
                .reciprocal_lattice;
        let path = vec![([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "G".into(), "G2".into())];
        let folded =
            fold_band_path_to_first_bz(reciprocal, &path, 2, 1e-8).expect("folded cubic path");
        assert_eq!(folded.len(), 1);
        assert!(folded[0].start.iter().all(|value| value.abs() <= 1e-12));
        assert!(folded[0].end.iter().all(|value| value.abs() <= 1e-12));
    }
}
