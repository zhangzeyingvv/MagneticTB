use nalgebra::{Complex, DMatrix};

use crate::{HoppingData, PropertiesError, PropertiesResult, Translation};

#[derive(Clone, Debug)]
pub struct ParsedWannier90Hr {
    pub data: HoppingData,
    pub header: String,
}

fn invalid(tag: &str, detail: impl Into<String>) -> PropertiesError {
    PropertiesError::new(tag, detail)
}

fn integer(token: &str, field: &str) -> PropertiesResult<i64> {
    let trimmed = token.trim();
    let digits = trimmed
        .strip_prefix('+')
        .or_else(|| trimmed.strip_prefix('-'))
        .unwrap_or(trimmed);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(invalid(
            "InvalidWannier90HR",
            format!("{field} is not an integer"),
        ));
    }
    trimmed.parse::<i64>().map_err(|_| {
        invalid(
            "Wannier90IntegerOverflow",
            format!("{field} is outside the signed 64-bit range"),
        )
    })
}

fn real(token: &str, field: &str) -> PropertiesResult<f64> {
    let value = token
        .replace(['d', 'D'], "e")
        .parse::<f64>()
        .map_err(|_| invalid("InvalidWannier90HR", format!("{field} is not real")))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "InvalidWannier90HR",
            format!("{field} must be finite"),
        ))
    }
}

fn positive_count(line: &str, field: &str) -> PropertiesResult<usize> {
    let tokens = line.split_whitespace().collect::<Vec<_>>();
    if tokens.len() != 1 {
        return Err(invalid(
            "InvalidWannier90HR",
            format!("{field} must occupy one header field"),
        ));
    }
    usize::try_from(integer(tokens[0], field)?)
        .ok()
        .filter(|value| *value > 0)
        .ok_or_else(|| invalid("InvalidWannier90HR", format!("{field} must be positive")))
}

pub fn parse_wannier90_hr(
    text: &str,
    source: Option<String>,
    precision_cutoff: f64,
    cell_cutoff: Option<i64>,
) -> PropertiesResult<ParsedWannier90Hr> {
    if !precision_cutoff.is_finite() || precision_cutoff < 0.0 {
        return Err(invalid(
            "InvalidWannier90Option",
            "prec must be finite and nonnegative",
        ));
    }
    if cell_cutoff.is_some_and(|value| value < 0) {
        return Err(invalid(
            "InvalidWannier90Option",
            "ncell must be All or a nonnegative integer",
        ));
    }
    let mut lines = text.lines().collect::<Vec<_>>();
    while lines.last().is_some_and(|line| line.trim().is_empty()) {
        lines.pop();
    }
    if lines.len() < 4 {
        return Err(invalid(
            "InvalidWannier90HR",
            "the file contains too few lines",
        ));
    }
    let num_wannier = positive_count(lines[1], "NumWannier")?;
    let num_translations = positive_count(lines[2], "NumTranslations")?;
    let degeneracy_line_count = num_translations.div_ceil(15);
    if lines.len() < 3 + degeneracy_line_count {
        return Err(invalid(
            "InvalidWannier90HR",
            "the degeneracy table is incomplete",
        ));
    }
    let degeneracy_tokens = lines[3..3 + degeneracy_line_count]
        .iter()
        .flat_map(|line| line.split_whitespace())
        .collect::<Vec<_>>();
    if degeneracy_tokens.len() != num_translations {
        return Err(invalid(
            "InvalidWannier90HR",
            "the degeneracy count does not match NumTranslations",
        ));
    }
    let input_degeneracies = degeneracy_tokens
        .iter()
        .map(|token| {
            u64::try_from(integer(token, "degeneracy")?)
                .ok()
                .filter(|value| *value > 0)
                .ok_or_else(|| invalid("InvalidWannier90HR", "all degeneracies must be positive"))
        })
        .collect::<PropertiesResult<Vec<_>>>()?;
    let block_size = num_wannier.checked_mul(num_wannier).ok_or_else(|| {
        invalid(
            "DimensionOverflow",
            "NumWannier squared overflows the platform size",
        )
    })?;
    let expected_records = num_translations.checked_mul(block_size).ok_or_else(|| {
        invalid(
            "DimensionOverflow",
            "wannier90 record count overflows the platform size",
        )
    })?;
    let records = &lines[3 + degeneracy_line_count..];
    if records.len() != expected_records {
        return Err(invalid(
            "InvalidWannier90HR",
            format!(
                "expected {expected_records} hopping records but found {}",
                records.len()
            ),
        ));
    }
    let mut translations = Vec::with_capacity(num_translations);
    let mut degeneracies = Vec::with_capacity(num_translations);
    let mut matrices = Vec::with_capacity(num_translations);
    for (block_index, block) in records.chunks_exact(block_size).enumerate() {
        let mut translation: Option<Translation> = None;
        let mut seen = vec![false; block_size];
        let mut matrix = DMatrix::zeros(num_wannier, num_wannier);
        for line in block {
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            if tokens.len() != 7 {
                return Err(invalid(
                    "InvalidWannier90HR",
                    "each hopping record must contain seven fields",
                ));
            }
            let current = [
                integer(tokens[0], "R1")?,
                integer(tokens[1], "R2")?,
                integer(tokens[2], "R3")?,
            ];
            if translation.is_some_and(|value| value != current) {
                return Err(invalid(
                    "InvalidWannier90HR",
                    "one translation block mixes distinct translations",
                ));
            }
            translation = Some(current);
            let row = integer(tokens[3], "row")?
                .checked_sub(1)
                .and_then(|value| usize::try_from(value).ok())
                .filter(|value| *value < num_wannier)
                .ok_or_else(|| invalid("InvalidWannier90HR", "row index is out of range"))?;
            let column = integer(tokens[4], "column")?
                .checked_sub(1)
                .and_then(|value| usize::try_from(value).ok())
                .filter(|value| *value < num_wannier)
                .ok_or_else(|| invalid("InvalidWannier90HR", "column index is out of range"))?;
            let pair = row * num_wannier + column;
            if std::mem::replace(&mut seen[pair], true) {
                return Err(invalid(
                    "InvalidWannier90HR",
                    "a translation block contains duplicate matrix indices",
                ));
            }
            let mut re = real(tokens[5], "real hopping")?;
            let mut im = real(tokens[6], "imaginary hopping")?;
            if re.abs() < precision_cutoff {
                re = 0.0;
            }
            if im.abs() < precision_cutoff {
                im = 0.0;
            }
            matrix[(row, column)] = Complex::new(re, im);
        }
        let translation = translation.ok_or_else(|| {
            invalid(
                "InvalidWannier90HR",
                "a nonempty hopping block did not contain a translation",
            )
        })?;
        let keep = cell_cutoff.map(i64::unsigned_abs).is_none_or(|cutoff| {
            translation
                .iter()
                .all(|coordinate| coordinate.unsigned_abs() <= cutoff)
        });
        if keep {
            translations.push(translation);
            degeneracies.push(input_degeneracies[block_index]);
            matrices.push(matrix);
        }
    }
    if translations.is_empty() {
        return Err(invalid(
            "InvalidWannier90HR",
            "ncell removed every translation block",
        ));
    }
    let mut unique = translations.clone();
    unique.sort_unstable();
    unique.dedup();
    if unique.len() != translations.len() {
        return Err(invalid(
            "InvalidWannier90HR",
            "translation blocks must be duplicate-free",
        ));
    }
    let data = HoppingData {
        num_wannier,
        translations,
        degeneracies,
        matrices,
        lattice: None,
        wannier_centers: None,
        source,
        cell_matrix: None,
        cell_representatives: Vec::new(),
    };
    data.validate()?;
    Ok(ParsedWannier90Hr {
        data,
        header: lines[0].to_owned(),
    })
}

fn formatted_real(value: f64, digits: usize) -> PropertiesResult<String> {
    if !value.is_finite() {
        return Err(invalid(
            "InvalidHoppingData",
            "Wannier90 export requires finite hopping values",
        ));
    }
    let exponent = -i32::try_from(digits).unwrap_or(i32::MAX);
    let value = if value.abs() < 10.0_f64.powi(exponent) {
        0.0
    } else {
        value
    };
    Ok(format!("{value:.digits$}"))
}

pub fn format_wannier90_hr(data: &HoppingData, digits: usize) -> PropertiesResult<String> {
    if digits < 6 {
        return Err(invalid(
            "InvalidWannier90Option",
            "RealDigits must be at least six",
        ));
    }
    data.validate()?;
    let mut lines = vec![
        "Generated by MagneticTB".to_owned(),
        data.num_wannier.to_string(),
        data.translations.len().to_string(),
    ];
    for chunk in data.degeneracies.chunks(15) {
        lines.push(
            chunk
                .iter()
                .map(u64::to_string)
                .collect::<Vec<_>>()
                .join(" "),
        );
    }
    for (translation, matrix) in data.translations.iter().zip(&data.matrices) {
        for column in 0..data.num_wannier {
            for row in 0..data.num_wannier {
                let value = matrix[(row, column)];
                lines.push(format!(
                    "{} {} {} {} {} {} {}",
                    translation[0],
                    translation[1],
                    translation[2],
                    row + 1,
                    column + 1,
                    formatted_real(value.re, digits)?,
                    formatted_real(value.im, digits)?,
                ));
            }
        }
    }
    Ok(lines.join("\n") + "\n")
}

#[cfg(test)]
mod tests {
    use super::{format_wannier90_hr, parse_wannier90_hr};
    use crate::HoppingData;
    use nalgebra::{Complex, DMatrix};

    #[test]
    fn stable_hr_round_trip_and_filter() {
        let data = HoppingData {
            num_wannier: 1,
            translations: vec![[-1, 0, 0], [0, 0, 0], [1, 0, 0]],
            degeneracies: vec![1, 2, 1],
            matrices: vec![
                DMatrix::from_element(1, 1, Complex::new(0.25, 0.0)),
                DMatrix::from_element(1, 1, Complex::new(1.0, 0.0)),
                DMatrix::from_element(1, 1, Complex::new(0.25, 0.0)),
            ],
            lattice: None,
            wannier_centers: None,
            source: None,
            cell_matrix: None,
            cell_representatives: Vec::new(),
        };
        let text = format_wannier90_hr(&data, 12).expect("format");
        let parsed =
            parse_wannier90_hr(&text, Some("memory".to_owned()), 1.0e-10, None).expect("parse");
        assert_eq!(parsed.header, "Generated by MagneticTB");
        assert_eq!(parsed.data.translations, data.translations);
        assert_eq!(parsed.data.degeneracies, data.degeneracies);
        let filtered = parse_wannier90_hr(&text, None, 1.0e-10, Some(0)).expect("filter");
        assert_eq!(filtered.data.translations, vec![[0, 0, 0]]);
        assert_eq!(filtered.data.degeneracies, vec![2]);
    }
}
