use crate::{FittingError, FittingResult};

#[derive(Clone, Debug, PartialEq)]
pub struct VaspEigenvalueRecord {
    pub k_point: [f64; 3],
    pub energies: Vec<f64>,
}

fn invalid(tag: &str, detail: impl Into<String>) -> FittingError {
    FittingError::new(tag, detail)
}

fn valid_number_token(token: &str) -> bool {
    let bytes = token.as_bytes();
    if bytes.is_empty() {
        return false;
    }
    let mut index = usize::from(matches!(bytes[0], b'+' | b'-'));
    let integer_start = index;
    while index < bytes.len() && bytes[index].is_ascii_digit() {
        index += 1;
    }
    let integer_digits = index - integer_start;
    let mut fractional_digits = 0;
    if index < bytes.len() && bytes[index] == b'.' {
        index += 1;
        let fractional_start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        fractional_digits = index - fractional_start;
    }
    if integer_digits + fractional_digits == 0 {
        return false;
    }
    if index < bytes.len() && matches!(bytes[index], b'E' | b'e') {
        index += 1;
        if index < bytes.len() && matches!(bytes[index], b'+' | b'-') {
            index += 1;
        }
        let exponent_start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == exponent_start {
            return false;
        }
    }
    index == bytes.len()
}

fn numeric_line(line: &str) -> Option<Vec<f64>> {
    let tokens = line.split_whitespace().collect::<Vec<_>>();
    if tokens.is_empty() || tokens.iter().any(|token| !valid_number_token(token)) {
        return None;
    }
    tokens
        .iter()
        .map(|token| token.parse::<f64>().ok().filter(|value| value.is_finite()))
        .collect()
}

fn integer_token(token: &str) -> Option<usize> {
    if token.is_empty()
        || token
            .bytes()
            .enumerate()
            .any(|(index, byte)| !(byte.is_ascii_digit() || index == 0 && byte == b'+'))
    {
        return None;
    }
    token.parse::<usize>().ok()
}

pub fn parse_vasp_eigenval(
    text: &str,
    fermi_energy: f64,
    spin: usize,
    start_band: usize,
    end_band: usize,
) -> FittingResult<Vec<VaspEigenvalueRecord>> {
    if !fermi_energy.is_finite() {
        return Err(invalid(
            "InvalidVaspEigenvalFormat",
            "Fermi energy must be finite",
        ));
    }
    let lines = text.lines().collect::<Vec<_>>();
    if lines.len() < 6 {
        return Err(invalid(
            "InvalidVaspEigenvalFormat",
            "EIGENVAL must contain at least six header lines",
        ));
    }
    let first_tokens = lines[0].split_whitespace().collect::<Vec<_>>();
    if numeric_line(lines[0]).is_none() {
        return Err(invalid(
            "InvalidVaspEigenvalFormat",
            "first header line must be fully numeric",
        ));
    }
    let spin_count = first_tokens
        .last()
        .and_then(|token| integer_token(token))
        .filter(|value| matches!(value, 1 | 2))
        .ok_or_else(|| {
            invalid(
                "InvalidVaspEigenvalFormat",
                "first header line must end in spin count 1 or 2",
            )
        })?;
    if spin == 0 || spin > spin_count {
        return Err(invalid(
            "InvalidVaspSpin",
            format!("spin must be between 1 and {spin_count}; received {spin}"),
        ));
    }
    let header_tokens = lines[5].split_whitespace().collect::<Vec<_>>();
    if header_tokens.len() < 2 || numeric_line(lines[5]).is_none() {
        return Err(invalid(
            "InvalidVaspEigenvalFormat",
            "sixth header line must end in k-point and band counts",
        ));
    }
    let k_point_count = integer_token(header_tokens[header_tokens.len() - 2])
        .filter(|value| *value > 0)
        .ok_or_else(|| invalid("InvalidVaspEigenvalFormat", "invalid k-point count"))?;
    let band_count = integer_token(header_tokens[header_tokens.len() - 1])
        .filter(|value| *value > 0)
        .ok_or_else(|| invalid("InvalidVaspEigenvalFormat", "invalid band count"))?;
    if start_band == 0 || end_band < start_band || end_band > band_count {
        return Err(invalid(
            "InvalidVaspBandRange",
            format!(
                "requested bands {start_band} through {end_band} are outside 1 through {band_count}"
            ),
        ));
    }
    let data_lines = lines[6..]
        .iter()
        .filter(|line| !line.trim().is_empty())
        .copied()
        .collect::<Vec<_>>();
    if data_lines.len() != k_point_count * (band_count + 1) {
        return Err(invalid(
            "InvalidVaspEigenvalFormat",
            "EIGENVAL data does not contain complete k-point blocks",
        ));
    }
    let mut result = Vec::with_capacity(k_point_count);
    for block in data_lines.chunks_exact(band_count + 1) {
        let point = numeric_line(block[0]).ok_or_else(|| {
            invalid(
                "InvalidVaspEigenvalFormat",
                "k-point record contains a nonnumeric token",
            )
        })?;
        if point.len() < 3 {
            return Err(invalid(
                "InvalidVaspEigenvalFormat",
                "k-point record must contain three coordinates",
            ));
        }
        let band_rows = block[1..]
            .iter()
            .map(|line| {
                numeric_line(line).ok_or_else(|| {
                    invalid(
                        "InvalidVaspEigenvalFormat",
                        "band record contains a nonnumeric token",
                    )
                })
            })
            .collect::<FittingResult<Vec<_>>>()?;
        if band_rows.iter().any(|row| row.len() <= spin) {
            return Err(invalid(
                "InvalidVaspSpin",
                format!("spin column {spin} is not available"),
            ));
        }
        let energies = band_rows[start_band - 1..end_band]
            .iter()
            .map(|row| row[spin] - fermi_energy)
            .collect();
        result.push(VaspEigenvalueRecord {
            k_point: [point[0], point[1], point[2]],
            energies,
        });
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "1 1 1 1\nheader\nheader\nheader\nheader\n1 2 2\n\n0 0 0 1\n1 -1 0\n2 2 0\n\n0.5 0 0 1\n1 -0.5 0\n2 3 0\n";

    #[test]
    fn parses_selected_bands_and_subtracts_fermi_energy() {
        let records = parse_vasp_eigenval(SAMPLE, 0.5, 1, 1, 2).expect("valid EIGENVAL");
        assert_eq!(records.len(), 2);
        assert!(records[0].k_point.iter().all(|value| value.abs() < 1.0e-15));
        assert!((records[0].energies[0] + 1.5).abs() < 1.0e-15);
        assert!((records[0].energies[1] - 1.5).abs() < 1.0e-15);
    }

    #[test]
    fn rejects_executable_or_nonfinite_tokens() {
        let malicious = SAMPLE.replace("1 -1 0", "1 0;System`Run[\"bad\"] 0");
        let error =
            parse_vasp_eigenval(&malicious, 0.0, 1, 1, 1).expect_err("executable token must fail");
        assert_eq!(error.tag(), "InvalidVaspEigenvalFormat");
        let nonfinite = SAMPLE.replace("1 -1 0", "1 NaN 0");
        assert!(parse_vasp_eigenval(&nonfinite, 0.0, 1, 1, 1).is_err());
    }
}
