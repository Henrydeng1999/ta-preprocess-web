use crate::TaData;

pub fn parse_csv(text: &str) -> TaData {
    let mut valid_lines: Vec<Vec<f64>> = Vec::new();

    for line in text.lines() {
        let stripped = line.trim();
        if stripped.is_empty() {
            continue;
        }
        let first_val = stripped.split(',').next().unwrap_or("").trim();
        if first_val.parse::<f64>().is_err() {
            break;
        }
        let row: Vec<f64> = stripped
            .split(',')
            .map(|v| {
                let v = v.trim();
                if v.is_empty() {
                    f64::NAN
                } else {
                    v.parse().unwrap_or(f64::NAN)
                }
            })
            .collect();
        valid_lines.push(row);
    }

    if valid_lines.is_empty() || valid_lines[0].len() < 2 {
        return TaData {
            time_array: vec![],
            wavelength_array: vec![],
            ta_2d: vec![],
        };
    }

    let n_cols = valid_lines[0].len();
    let time_array: Vec<f64> = valid_lines[0][1..].to_vec();
    let mut wavelength_array = Vec::new();
    let mut ta_2d = Vec::new();

    for i in 1..valid_lines.len() {
        wavelength_array.push(valid_lines[i][0]);
        let mut row = valid_lines[i][1..].to_vec();
        if row.len() < n_cols - 1 {
            row.resize(n_cols - 1, f64::NAN);
        } else if row.len() > n_cols - 1 {
            row.truncate(n_cols - 1);
        }
        ta_2d.push(row);
    }

    TaData {
        time_array,
        wavelength_array,
        ta_2d,
    }
}
