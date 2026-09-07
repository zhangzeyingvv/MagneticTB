use std::collections::{HashMap, VecDeque};
use std::f64::consts::{PI, TAU};

use nalgebra::{Complex, DMatrix};

use crate::{PropertiesError, PropertiesResult};

pub type Translation = [i64; 3];
pub type CellMatrix = [[i64; 3]; 3];

#[derive(Clone, Debug)]
pub struct HoppingData {
    pub num_wannier: usize,
    pub translations: Vec<Translation>,
    pub degeneracies: Vec<u64>,
    pub matrices: Vec<DMatrix<Complex<f64>>>,
    pub lattice: Option<[[f64; 3]; 3]>,
    pub wannier_centers: Option<Vec<[f64; 3]>>,
    pub source: Option<String>,
    pub cell_matrix: Option<CellMatrix>,
    pub cell_representatives: Vec<Translation>,
}

#[derive(Clone, Debug)]
pub struct BlochHamiltonianData {
    pub hamiltonian: DMatrix<Complex<f64>>,
    pub momentum: [f64; 3],
    pub hermitian_residual: f64,
    pub cell_volume_factor: u64,
    pub cell_representatives: Vec<Translation>,
}

#[derive(Clone, Debug)]
pub enum FiniteGeometry {
    Rectangular([usize; 3]),
    ExplicitCells(Vec<Translation>),
}

#[derive(Clone, Debug)]
pub struct FiniteBasisRecord {
    pub basis_index: usize,
    pub cell_index: usize,
    pub cell: Translation,
    pub orbital_index: usize,
    pub wannier_center: Option<[f64; 3]>,
    pub fractional_position: Option<[f64; 3]>,
    pub cartesian_position: Option<[f64; 3]>,
}

#[derive(Clone, Debug)]
pub struct FiniteHamiltonianData {
    pub hamiltonian: DMatrix<Complex<f64>>,
    pub shape: &'static str,
    pub cells: Vec<Translation>,
    pub size: Option<[usize; 3]>,
    pub boundary_conditions: [bool; 3],
    pub orbitals_per_cell: usize,
    pub hermitian_residual: f64,
    pub cell_bounds: [[i64; 2]; 3],
    pub position_data_available: bool,
    pub basis_records: Vec<FiniteBasisRecord>,
    pub lattice: Option<[[f64; 3]; 3]>,
    pub wannier_centers: Option<Vec<[f64; 3]>>,
}

#[derive(Clone, Debug)]
pub struct SlabHamiltonianData {
    pub hamiltonian: DMatrix<Complex<f64>>,
    pub cells: Vec<Translation>,
    pub size: [usize; 3],
    pub periodic_directions: Vec<usize>,
    pub open_directions: Vec<usize>,
    pub momentum: Vec<f64>,
    pub embedded_momentum: [f64; 3],
    pub orbitals_per_cell: usize,
    pub hermitian_residual: f64,
}

#[derive(Clone, Debug)]
pub struct PrincipalLayerBlocks {
    pub h00: DMatrix<Complex<f64>>,
    pub h01: DMatrix<Complex<f64>>,
    pub h10: DMatrix<Complex<f64>>,
    pub thickness: usize,
    pub orbitals_per_cell: usize,
}

#[derive(Clone, Copy, Debug)]
pub enum SurfaceSide {
    Positive,
    Negative,
}

#[derive(Clone, Copy, Debug)]
pub struct SurfaceOptions {
    pub broadening: f64,
    pub tolerance: f64,
    pub max_iterations: usize,
    pub side: SurfaceSide,
    pub hermiticity_tolerance: f64,
}

impl Default for SurfaceOptions {
    fn default() -> Self {
        Self {
            broadening: 1.0e-3,
            tolerance: 1.0e-10,
            max_iterations: 200,
            side: SurfaceSide::Positive,
            hermiticity_tolerance: 1.0e-9,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SurfaceGreenData {
    pub green_function: DMatrix<Complex<f64>>,
    pub spectral_weight: f64,
    pub energy: f64,
    pub momentum: [f64; 2],
    pub blocks: PrincipalLayerBlocks,
    pub iterations: usize,
    pub coupling_residual: f64,
    pub dyson_residual: f64,
    pub hopping_hermitian_residual: f64,
    pub block_hermitian_residual: f64,
}

fn invalid(tag: &str, detail: impl Into<String>) -> PropertiesError {
    PropertiesError::new(tag, detail)
}

fn finite_complex(value: Complex<f64>) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn matrix_max_norm(matrix: &DMatrix<Complex<f64>>) -> f64 {
    matrix.iter().map(|value| value.norm()).fold(0.0, f64::max)
}

fn matrix_frobenius_norm(matrix: &DMatrix<Complex<f64>>) -> f64 {
    matrix.iter().map(Complex::norm_sqr).sum::<f64>().sqrt()
}

fn checked_translation_add(left: Translation, right: Translation) -> PropertiesResult<Translation> {
    let mut result = [0; 3];
    for index in 0..3 {
        result[index] = left[index].checked_add(right[index]).ok_or_else(|| {
            invalid(
                "IntegerOverflow",
                "hopping translation addition overflows i64",
            )
        })?;
    }
    Ok(result)
}

impl HoppingData {
    pub fn validate(&self) -> PropertiesResult<()> {
        if self.num_wannier == 0
            || self.translations.is_empty()
            || self.translations.len() != self.degeneracies.len()
            || self.translations.len() != self.matrices.len()
            || self.degeneracies.contains(&0)
        {
            return Err(invalid(
                "InvalidHoppingData",
                "NumWannier, translations, degeneracies, and hopping matrices are inconsistent",
            ));
        }
        let mut unique = self.translations.clone();
        unique.sort_unstable();
        unique.dedup();
        if unique.len() != self.translations.len() {
            return Err(invalid(
                "InvalidHoppingData",
                "hopping translations must be duplicate-free",
            ));
        }
        if self.matrices.iter().any(|matrix| {
            matrix.nrows() != self.num_wannier
                || matrix.ncols() != self.num_wannier
                || matrix.iter().any(|&value| !finite_complex(value))
        }) {
            return Err(invalid(
                "InvalidHoppingData",
                "every hopping matrix must be finite and NumWannier by NumWannier",
            ));
        }
        if let Some(lattice) = self.lattice
            && lattice.iter().flatten().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "InvalidHoppingData",
                "lattice entries must be finite",
            ));
        }
        if let Some(centers) = &self.wannier_centers
            && (centers.len() != self.num_wannier
                || centers.iter().flatten().any(|value| !value.is_finite()))
        {
            return Err(invalid(
                "InvalidHoppingData",
                "one finite Wannier center is required per orbital",
            ));
        }
        Ok(())
    }

    pub fn effective_matrices(&self) -> PropertiesResult<Vec<DMatrix<Complex<f64>>>> {
        self.validate()?;
        Ok(self
            .matrices
            .iter()
            .zip(&self.degeneracies)
            .map(|(matrix, &degeneracy)| matrix / Complex::new(degeneracy as f64, 0.0))
            .collect())
    }

    pub fn hermiticity_residual(&self) -> PropertiesResult<f64> {
        let matrices = self.effective_matrices()?;
        let by_translation = self
            .translations
            .iter()
            .copied()
            .zip(matrices.iter())
            .collect::<HashMap<_, _>>();
        let zero = DMatrix::zeros(self.num_wannier, self.num_wannier);
        Ok(self
            .translations
            .iter()
            .zip(&matrices)
            .map(|(translation, matrix)| {
                let reverse = [-translation[0], -translation[1], -translation[2]];
                let opposite = by_translation.get(&reverse).copied().unwrap_or(&zero);
                matrix_max_norm(&(matrix - opposite.adjoint()))
            })
            .fold(0.0, f64::max))
    }
}

fn require_hermitian(data: &HoppingData, tolerance: f64) -> PropertiesResult<f64> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(invalid(
            "InvalidPropertiesOption",
            "HermiticityTolerance must be finite and nonnegative",
        ));
    }
    let residual = data.hermiticity_residual()?;
    if residual > tolerance {
        return Err(invalid(
            "NonHermitianHoppingData",
            format!("hopping Hermiticity residual {residual} exceeds tolerance {tolerance}"),
        ));
    }
    Ok(residual)
}

pub fn build_bloch_hamiltonian(
    data: &HoppingData,
    momentum: [f64; 3],
    hermiticity_tolerance: f64,
) -> PropertiesResult<BlochHamiltonianData> {
    if momentum.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "InvalidBlochMomentum",
            "crystal momentum must be a finite reciprocal-fractional three-vector",
        ));
    }
    require_hermitian(data, hermiticity_tolerance)?;
    let matrices = data.effective_matrices()?;
    let mut hamiltonian = DMatrix::zeros(data.num_wannier, data.num_wannier);
    for (matrix, translation) in matrices.iter().zip(&data.translations) {
        let angle = TAU
            * momentum
                .iter()
                .zip(translation)
                .map(|(left, &right)| left * right as f64)
                .sum::<f64>();
        hamiltonian += matrix * Complex::from_polar(1.0, angle);
    }
    let residual = matrix_max_norm(&(&hamiltonian - hamiltonian.adjoint()));
    if residual > hermiticity_tolerance {
        return Err(invalid(
            "NonHermitianBlochHamiltonian",
            format!(
                "Bloch Hermiticity residual {residual} exceeds tolerance {hermiticity_tolerance}"
            ),
        ));
    }
    let volume = data.cell_matrix.map_or(Ok(1), |matrix| {
        u64::try_from(determinant(matrix).unsigned_abs()).map_err(|_| {
            invalid(
                "IntegerOverflow",
                "cell volume factor is outside the supported u64 range",
            )
        })
    })?;
    Ok(BlochHamiltonianData {
        hamiltonian,
        momentum,
        hermitian_residual: residual,
        cell_volume_factor: volume,
        cell_representatives: if data.cell_representatives.is_empty() {
            vec![[0, 0, 0]]
        } else {
            data.cell_representatives.clone()
        },
    })
}

fn rectangular_cells(size: [usize; 3]) -> PropertiesResult<Vec<Translation>> {
    if size.contains(&0) {
        return Err(invalid(
            "InvalidFiniteGeometry",
            "rectangular sizes must be positive",
        ));
    }
    let mut cells = Vec::with_capacity(size.iter().product());
    for x in 0..size[0] {
        for y in 0..size[1] {
            for z in 0..size[2] {
                cells.push([
                    i64::try_from(x)
                        .map_err(|_| invalid("IntegerOverflow", "cell index overflows i64"))?,
                    i64::try_from(y)
                        .map_err(|_| invalid("IntegerOverflow", "cell index overflows i64"))?,
                    i64::try_from(z)
                        .map_err(|_| invalid("IntegerOverflow", "cell index overflows i64"))?,
                ]);
            }
        }
    }
    Ok(cells)
}

pub fn build_real_space_hamiltonian(
    data: &HoppingData,
    geometry: FiniteGeometry,
    boundary_conditions: [bool; 3],
    hermiticity_tolerance: f64,
) -> PropertiesResult<FiniteHamiltonianData> {
    require_hermitian(data, hermiticity_tolerance)?;
    let (shape, cells, size) = match geometry {
        FiniteGeometry::Rectangular(size) => ("Rectangular", rectangular_cells(size)?, Some(size)),
        FiniteGeometry::ExplicitCells(cells) => {
            if cells.is_empty() {
                return Err(invalid(
                    "InvalidFiniteGeometry",
                    "explicit cell set must be nonempty",
                ));
            }
            let mut unique = cells.clone();
            unique.sort_unstable();
            unique.dedup();
            if unique.len() != cells.len() {
                return Err(invalid(
                    "InvalidFiniteGeometry",
                    "explicit cells must be duplicate-free",
                ));
            }
            if boundary_conditions.iter().any(|&periodic| periodic) {
                return Err(invalid(
                    "InvalidBoundaryConditions",
                    "explicit cell sets support open boundaries only",
                ));
            }
            ("ExplicitCells", cells, None)
        }
    };
    let cell_indices = cells
        .iter()
        .copied()
        .enumerate()
        .map(|(index, cell)| (cell, index))
        .collect::<HashMap<_, _>>();
    let matrices = data.effective_matrices()?;
    let dimension = cells.len().checked_mul(data.num_wannier).ok_or_else(|| {
        invalid(
            "DimensionOverflow",
            "finite Hamiltonian dimension overflows usize",
        )
    })?;
    let mut hamiltonian = DMatrix::zeros(dimension, dimension);
    for (source_index, &source) in cells.iter().enumerate() {
        for (translation, matrix) in data.translations.iter().zip(&matrices) {
            let mut target = checked_translation_add(source, *translation)?;
            if let Some(size) = size {
                for axis in 0..3 {
                    if boundary_conditions[axis] {
                        let modulus = i64::try_from(size[axis]).map_err(|_| {
                            invalid("IntegerOverflow", "periodic size overflows i64")
                        })?;
                        target[axis] = target[axis].rem_euclid(modulus);
                    }
                }
            }
            if let Some(&target_index) = cell_indices.get(&target) {
                for row in 0..data.num_wannier {
                    for column in 0..data.num_wannier {
                        hamiltonian[(
                            source_index * data.num_wannier + row,
                            target_index * data.num_wannier + column,
                        )] += matrix[(row, column)];
                    }
                }
            }
        }
    }
    let residual = matrix_max_norm(&(&hamiltonian - hamiltonian.adjoint()));
    let cell_bounds = std::array::from_fn(|axis| {
        let mut values = cells.iter().map(|cell| cell[axis]);
        let first = values.next().unwrap_or(0);
        values.fold([first, first], |[minimum, maximum], value| {
            [minimum.min(value), maximum.max(value)]
        })
    });
    let position_data_available = data.lattice.is_some() && data.wannier_centers.is_some();
    let basis_records = cells
        .iter()
        .enumerate()
        .flat_map(|(cell_index, &cell)| {
            (0..data.num_wannier).map(move |orbital_index| {
                let wannier_center = data
                    .wannier_centers
                    .as_ref()
                    .map(|centers| centers[orbital_index]);
                let fractional_position = wannier_center
                    .map(|center| std::array::from_fn(|axis| cell[axis] as f64 + center[axis]));
                let cartesian_position =
                    fractional_position
                        .zip(data.lattice)
                        .map(|(position, lattice)| {
                            std::array::from_fn(|column| {
                                (0..3).map(|row| position[row] * lattice[row][column]).sum()
                            })
                        });
                FiniteBasisRecord {
                    basis_index: cell_index * data.num_wannier + orbital_index + 1,
                    cell_index: cell_index + 1,
                    cell,
                    orbital_index: orbital_index + 1,
                    wannier_center,
                    fractional_position,
                    cartesian_position,
                }
            })
        })
        .collect();
    Ok(FiniteHamiltonianData {
        hamiltonian,
        shape,
        cells,
        size,
        boundary_conditions,
        orbitals_per_cell: data.num_wannier,
        hermitian_residual: residual,
        cell_bounds,
        position_data_available,
        basis_records,
        lattice: data.lattice,
        wannier_centers: data.wannier_centers.clone(),
    })
}

pub fn build_slab_hamiltonian(
    data: &HoppingData,
    size: [usize; 3],
    periodic_directions: &[usize],
    momentum: &[f64],
    hermiticity_tolerance: f64,
) -> PropertiesResult<SlabHamiltonianData> {
    require_hermitian(data, hermiticity_tolerance)?;
    if !matches!(periodic_directions.len(), 1 | 2)
        || periodic_directions.iter().any(|&axis| axis >= 3)
        || {
            let mut directions = periodic_directions.to_vec();
            directions.sort_unstable();
            directions.dedup();
            directions.len() != periodic_directions.len()
        }
    {
        return Err(invalid(
            "InvalidPeriodicDirections",
            "PeriodicDirections must contain one or two distinct axes",
        ));
    }
    if size.contains(&0) || periodic_directions.iter().any(|&axis| size[axis] != 1) {
        return Err(invalid(
            "InvalidSlabSize",
            "periodic directions must have size one and all sizes must be positive",
        ));
    }
    if momentum.len() != periodic_directions.len()
        || momentum.iter().any(|value| !value.is_finite())
    {
        return Err(invalid(
            "InvalidSurfaceMomentum",
            "one finite reciprocal-fractional momentum is required per periodic direction",
        ));
    }
    let open_directions = (0..3)
        .filter(|axis| !periodic_directions.contains(axis))
        .collect::<Vec<_>>();
    let cells = rectangular_cells(size)?;
    let indices = cells
        .iter()
        .copied()
        .enumerate()
        .map(|(index, cell)| (cell, index))
        .collect::<HashMap<_, _>>();
    let mut embedded_momentum = [0.0; 3];
    for (&axis, &value) in periodic_directions.iter().zip(momentum) {
        embedded_momentum[axis] = value;
    }
    let matrices = data.effective_matrices()?;
    let dimension = cells.len() * data.num_wannier;
    let mut hamiltonian = DMatrix::zeros(dimension, dimension);
    for (source_index, &source) in cells.iter().enumerate() {
        for (translation, matrix) in data.translations.iter().zip(&matrices) {
            let mut target = checked_translation_add(source, *translation)?;
            if open_directions.iter().any(|&axis| {
                target[axis] < 0
                    || usize::try_from(target[axis]).map_or(true, |value| value >= size[axis])
            }) {
                continue;
            }
            let angle = TAU
                * embedded_momentum
                    .iter()
                    .zip(translation)
                    .map(|(left, &right)| left * right as f64)
                    .sum::<f64>();
            let phase = Complex::from_polar(1.0, angle);
            for &axis in periodic_directions {
                target[axis] = 0;
            }
            if let Some(&target_index) = indices.get(&target) {
                for row in 0..data.num_wannier {
                    for column in 0..data.num_wannier {
                        hamiltonian[(
                            source_index * data.num_wannier + row,
                            target_index * data.num_wannier + column,
                        )] += phase * matrix[(row, column)];
                    }
                }
            }
        }
    }
    let residual = matrix_max_norm(&(&hamiltonian - hamiltonian.adjoint()));
    Ok(SlabHamiltonianData {
        hamiltonian,
        cells,
        size,
        periodic_directions: periodic_directions.to_vec(),
        open_directions,
        momentum: momentum.to_vec(),
        embedded_momentum,
        orbitals_per_cell: data.num_wannier,
        hermitian_residual: residual,
    })
}

fn determinant(matrix: CellMatrix) -> i128 {
    let a = matrix.map(|row| row.map(i128::from));
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

fn adjugate(matrix: CellMatrix) -> [[i128; 3]; 3] {
    let a = matrix.map(|row| row.map(i128::from));
    [
        [
            a[1][1] * a[2][2] - a[1][2] * a[2][1],
            a[0][2] * a[2][1] - a[0][1] * a[2][2],
            a[0][1] * a[1][2] - a[0][2] * a[1][1],
        ],
        [
            a[1][2] * a[2][0] - a[1][0] * a[2][2],
            a[0][0] * a[2][2] - a[0][2] * a[2][0],
            a[0][2] * a[1][0] - a[0][0] * a[1][2],
        ],
        [
            a[1][0] * a[2][1] - a[1][1] * a[2][0],
            a[0][1] * a[2][0] - a[0][0] * a[2][1],
            a[0][0] * a[1][1] - a[0][1] * a[1][0],
        ],
    ]
}

fn row_times_i128(vector: Translation, matrix: [[i128; 3]; 3]) -> [i128; 3] {
    let vector = vector.map(i128::from);
    std::array::from_fn(|column| (0..3).map(|row| vector[row] * matrix[row][column]).sum())
}

fn quotient_key(vector: Translation, adjugate: [[i128; 3]; 3], modulus: i128) -> [i128; 3] {
    row_times_i128(vector, adjugate).map(|value| value.rem_euclid(modulus))
}

fn cell_representatives(matrix: CellMatrix) -> PropertiesResult<Vec<Translation>> {
    let determinant = determinant(matrix);
    if determinant == 0 {
        return Err(invalid(
            "InvalidCellMatrix",
            "cell matrix must be nonsingular",
        ));
    }
    let volume = determinant.unsigned_abs();
    let count = usize::try_from(volume)
        .map_err(|_| invalid("CellVolumeTooLarge", "cell volume does not fit in memory"))?;
    let adj = adjugate(matrix);
    let modulus = i128::try_from(volume).expect("u128 from i128 magnitude fits i128");
    let origin = [0, 0, 0];
    let mut records = HashMap::from([(quotient_key(origin, adj, modulus), origin)]);
    let mut queue = VecDeque::from([origin]);
    let steps = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    while records.len() < count {
        let Some(current) = queue.pop_front() else {
            return Err(invalid(
                "CellTransformFailure",
                "could not enumerate cell cosets",
            ));
        };
        for step in steps {
            let candidate = checked_translation_add(current, step)?;
            let key = quotient_key(candidate, adj, modulus);
            if let std::collections::hash_map::Entry::Vacant(entry) = records.entry(key) {
                entry.insert(candidate);
                queue.push_back(candidate);
                if records.len() == count {
                    break;
                }
            }
        }
    }
    let sign = determinant.signum();
    let mut representatives = records.into_values().collect::<Vec<_>>();
    representatives.sort_by_key(|representative| {
        row_times_i128(*representative, adj).map(|value| (value * sign).rem_euclid(modulus))
    });
    Ok(representatives)
}

fn transformed_translation(
    difference: Translation,
    adjugate: [[i128; 3]; 3],
    determinant: i128,
) -> PropertiesResult<Translation> {
    let numerator = row_times_i128(difference, adjugate);
    if numerator.iter().any(|value| value % determinant != 0) {
        return Err(invalid(
            "CellTransformFailure",
            "transformed hopping translation is not integral",
        ));
    }
    let mut result = [0; 3];
    for index in 0..3 {
        result[index] = i64::try_from(numerator[index] / determinant).map_err(|_| {
            invalid(
                "IntegerOverflow",
                "transformed hopping translation overflows i64",
            )
        })?;
    }
    Ok(result)
}

pub fn transform_hoppings(
    data: &HoppingData,
    cell: CellMatrix,
    tolerance: f64,
) -> PropertiesResult<HoppingData> {
    require_hermitian(data, tolerance)?;
    let det = determinant(cell);
    if det == 0 {
        return Err(invalid(
            "InvalidCellMatrix",
            "cell matrix must be nonsingular",
        ));
    }
    let adj = adjugate(cell);
    let modulus = det.unsigned_abs() as i128;
    let representatives = cell_representatives(cell)?;
    let key_to_index = representatives
        .iter()
        .copied()
        .enumerate()
        .map(|(index, representative)| (quotient_key(representative, adj, modulus), index))
        .collect::<HashMap<_, _>>();
    let old_dimension = data.num_wannier;
    let new_dimension = old_dimension
        .checked_mul(representatives.len())
        .ok_or_else(|| {
            invalid(
                "DimensionOverflow",
                "transformed orbital dimension overflows usize",
            )
        })?;
    let effective = data.effective_matrices()?;
    let mut output = HashMap::<Translation, DMatrix<Complex<f64>>>::new();
    for (&translation, matrix) in data.translations.iter().zip(&effective) {
        if matrix.iter().all(|value| value.norm() == 0.0) {
            continue;
        }
        for (source_index, &representative) in representatives.iter().enumerate() {
            let target_old = checked_translation_add(representative, translation)?;
            let target_index = *key_to_index
                .get(&quotient_key(target_old, adj, modulus))
                .ok_or_else(|| invalid("CellTransformFailure", "target coset is missing"))?;
            let difference = [
                target_old[0] - representatives[target_index][0],
                target_old[1] - representatives[target_index][1],
                target_old[2] - representatives[target_index][2],
            ];
            let transformed = transformed_translation(difference, adj, det)?;
            let combined = output
                .entry(transformed)
                .or_insert_with(|| DMatrix::zeros(new_dimension, new_dimension));
            for row in 0..old_dimension {
                for column in 0..old_dimension {
                    combined[(
                        source_index * old_dimension + row,
                        target_index * old_dimension + column,
                    )] += matrix[(row, column)];
                }
            }
        }
    }
    if output.is_empty() {
        output.insert([0, 0, 0], DMatrix::zeros(new_dimension, new_dimension));
    }
    let mut translations = output.keys().copied().collect::<Vec<_>>();
    translations.sort_unstable();
    let matrices = translations
        .iter()
        .map(|translation| {
            output.remove(translation).ok_or_else(|| {
                invalid(
                    "CellTransformFailure",
                    "a transformed translation disappeared during stable ordering",
                )
            })
        })
        .collect::<PropertiesResult<Vec<_>>>()?;
    let wannier_centers = data.wannier_centers.as_ref().map(|centers| {
        representatives
            .iter()
            .flat_map(|&representative| {
                centers.iter().map(move |center| {
                    let old = [
                        representative[0] as f64 + center[0],
                        representative[1] as f64 + center[1],
                        representative[2] as f64 + center[2],
                    ];
                    let inverse = adj.map(|row| row.map(|value| value as f64 / det as f64));
                    std::array::from_fn(|column| {
                        (0..3)
                            .map(|row| old[row] * inverse[row][column])
                            .sum::<f64>()
                            .rem_euclid(1.0)
                    })
                })
            })
            .collect()
    });
    let lattice = data.lattice.map(|lattice| {
        std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                (0..3)
                    .map(|inner| cell[row][inner] as f64 * lattice[inner][column])
                    .sum()
            })
        })
    });
    let result = HoppingData {
        num_wannier: new_dimension,
        translations,
        degeneracies: vec![1; matrices.len()],
        matrices,
        lattice,
        wannier_centers,
        source: data.source.clone(),
        cell_matrix: Some(cell),
        cell_representatives: representatives,
    };
    require_hermitian(&result, tolerance)?;
    Ok(result)
}

pub fn principal_layer_blocks(
    data: &HoppingData,
    momentum: [f64; 2],
) -> PropertiesResult<PrincipalLayerBlocks> {
    if momentum.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "InvalidSurfaceMomentum",
            "surface momentum must be finite",
        ));
    }
    let effective = data.effective_matrices()?;
    let normal_range = data
        .translations
        .iter()
        .map(|translation| translation[2].unsigned_abs() as usize)
        .max()
        .unwrap_or(0);
    let thickness = normal_range.max(1);
    let mut slices = HashMap::<i64, DMatrix<Complex<f64>>>::new();
    for (&translation, matrix) in data.translations.iter().zip(&effective) {
        let phase = Complex::from_polar(
            1.0,
            TAU * (momentum[0] * translation[0] as f64 + momentum[1] * translation[1] as f64),
        );
        *slices
            .entry(translation[2])
            .or_insert_with(|| DMatrix::zeros(data.num_wannier, data.num_wannier)) +=
            matrix * phase;
    }
    let dimension = thickness * data.num_wannier;
    let mut h00 = DMatrix::zeros(dimension, dimension);
    let mut h01 = DMatrix::zeros(dimension, dimension);
    let mut h10 = DMatrix::zeros(dimension, dimension);
    for source in 0..thickness {
        for target in 0..thickness {
            for (output, normal) in [
                (&mut h00, target as i64 - source as i64),
                (&mut h01, thickness as i64 + target as i64 - source as i64),
                (
                    &mut h10,
                    -(thickness as i64) + target as i64 - source as i64,
                ),
            ] {
                if let Some(slice) = slices.get(&normal) {
                    for row in 0..data.num_wannier {
                        for column in 0..data.num_wannier {
                            output[(
                                source * data.num_wannier + row,
                                target * data.num_wannier + column,
                            )] += slice[(row, column)];
                        }
                    }
                }
            }
        }
    }
    Ok(PrincipalLayerBlocks {
        h00,
        h01,
        h10,
        thickness,
        orbitals_per_cell: data.num_wannier,
    })
}

pub fn surface_green_function(
    data: &HoppingData,
    momentum: [f64; 2],
    energy: f64,
    options: SurfaceOptions,
) -> PropertiesResult<SurfaceGreenData> {
    if !energy.is_finite()
        || !options.broadening.is_finite()
        || options.broadening <= 0.0
        || !options.tolerance.is_finite()
        || options.tolerance <= 0.0
        || options.max_iterations == 0
    {
        return Err(invalid(
            "InvalidSurfaceOption",
            "surface options and energy are invalid",
        ));
    }
    let hopping_residual = require_hermitian(data, options.hermiticity_tolerance)?;
    let blocks = principal_layer_blocks(data, momentum)?;
    let h00_residual = matrix_frobenius_norm(&(&blocks.h00 - blocks.h00.adjoint()));
    let coupling_residual = matrix_frobenius_norm(&(&blocks.h10 - blocks.h01.adjoint()));
    let block_residual = h00_residual.max(coupling_residual);
    if block_residual > options.hermiticity_tolerance {
        return Err(invalid(
            "NonHermitianPrincipalLayer",
            format!(
                "principal-layer residual {block_residual} exceeds tolerance {}",
                options.hermiticity_tolerance
            ),
        ));
    }
    let dimension = blocks.h00.nrows();
    let identity = DMatrix::<Complex<f64>>::identity(dimension, dimension);
    let z = &identity * Complex::new(energy, options.broadening);
    let mut epsilon = blocks.h00.clone();
    let mut epsilon_surface = blocks.h00.clone();
    let mut alpha = blocks.h01.clone();
    let mut beta = blocks.h10.clone();
    let mut residual = f64::INFINITY;
    let mut iterations = 0;
    let mut converged = false;
    for current in 1..=options.max_iterations {
        let green = (&z - &epsilon).try_inverse().ok_or_else(|| {
            invalid(
                "SurfaceLinearSolveFailure",
                format!("decimation solve failed at iteration {current}"),
            )
        })?;
        let forward = &alpha * &green * &beta;
        let backward = &beta * &green * &alpha;
        let alpha_new = &alpha * &green * &alpha;
        let beta_new = &beta * &green * &beta;
        epsilon_surface += match options.side {
            SurfaceSide::Positive => &forward,
            SurfaceSide::Negative => &backward,
        };
        epsilon += forward + backward;
        residual = matrix_frobenius_norm(&alpha_new).max(matrix_frobenius_norm(&beta_new))
            / matrix_frobenius_norm(&epsilon).max(1.0);
        alpha = alpha_new;
        beta = beta_new;
        iterations = current;
        if residual <= options.tolerance {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(invalid(
            "SurfaceGreenFunctionNotConverged",
            format!(
                "decimation did not converge after {} iterations; residual {residual}",
                options.max_iterations
            ),
        ));
    }
    let inverse = &z - &epsilon_surface;
    let green = inverse
        .clone()
        .try_inverse()
        .ok_or_else(|| invalid("SurfaceLinearSolveFailure", "final surface solve failed"))?;
    let dyson = matrix_frobenius_norm(&(inverse * &green - &identity))
        / matrix_frobenius_norm(&identity).max(1.0);
    let spectral_weight = -(0..dimension)
        .map(|index| green[(index, index)].im)
        .sum::<f64>()
        / PI;
    Ok(SurfaceGreenData {
        green_function: green,
        spectral_weight,
        energy,
        momentum,
        blocks,
        iterations,
        coupling_residual: residual,
        dyson_residual: dyson,
        hopping_hermitian_residual: hopping_residual,
        block_hermitian_residual: block_residual,
    })
}

#[cfg(test)]
mod tests {
    use nalgebra::{Complex, DMatrix};

    use super::{
        FiniteGeometry, HoppingData, SurfaceOptions, build_bloch_hamiltonian,
        build_real_space_hamiltonian, surface_green_function, transform_hoppings,
    };

    fn scalar(value: f64) -> DMatrix<Complex<f64>> {
        DMatrix::from_element(1, 1, Complex::new(value, 0.0))
    }

    fn chain() -> HoppingData {
        HoppingData {
            num_wannier: 1,
            translations: vec![[0, 0, 0], [0, 0, 1], [0, 0, -1]],
            degeneracies: vec![1, 1, 1],
            matrices: vec![scalar(0.0), scalar(1.0), scalar(1.0)],
            lattice: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            wannier_centers: Some(vec![[0.0, 0.0, 0.0]]),
            source: Some("stable-test".to_owned()),
            cell_matrix: None,
            cell_representatives: Vec::new(),
        }
    }

    #[test]
    fn finite_chain_and_bloch_match_stable() {
        let data = chain();
        let finite = build_real_space_hamiltonian(
            &data,
            FiniteGeometry::Rectangular([1, 1, 5]),
            [false; 3],
            1.0e-9,
        )
        .expect("finite chain");
        assert_eq!(finite.hamiltonian.nrows(), 5);
        assert_eq!(finite.hamiltonian[(0, 1)], Complex::new(1.0, 0.0));
        assert_eq!(finite.cell_bounds, [[0, 0], [0, 0], [0, 4]]);
        assert!(finite.position_data_available);
        assert_eq!(finite.basis_records.len(), 5);
        assert_eq!(
            finite.basis_records[1].fractional_position,
            Some([0.0, 0.0, 1.0])
        );
        assert_eq!(
            finite.basis_records[1].cartesian_position,
            Some([0.0, 0.0, 1.0])
        );
        let bloch =
            build_bloch_hamiltonian(&data, [0.0, 0.0, 0.25], 1.0e-9).expect("Bloch Hamiltonian");
        assert!(bloch.hamiltonian[(0, 0)].norm() <= 1.0e-12);
    }

    #[test]
    fn doubled_cell_and_surface_match_stable() {
        let data = chain();
        let doubled = transform_hoppings(&data, [[1, 0, 0], [0, 1, 0], [0, 0, 2]], 1.0e-9)
            .expect("doubled chain");
        assert_eq!(doubled.num_wannier, 2);
        assert_eq!(doubled.cell_representatives.len(), 2);
        let surface = surface_green_function(&data, [0.0, 0.0], 0.0, SurfaceOptions::default())
            .expect("surface Green function");
        let eta = 1.0e-3;
        let analytic = (Complex::new(0.0, eta)
            - (Complex::new(0.0, eta).powi(2) - Complex::new(4.0, 0.0)).sqrt())
            / 2.0;
        assert!((surface.green_function[(0, 0)] - analytic).norm() <= 1.0e-10);
    }
}
