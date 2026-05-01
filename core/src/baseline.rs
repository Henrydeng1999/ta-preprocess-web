pub struct CropResult {
    pub wavelength_array: Vec<f64>,
    pub ta_2d: Vec<Vec<f64>>,
}

pub fn crop_wavelength(wl: &[f64], ta: &[Vec<f64>], wl_min: f64, wl_max: f64) -> CropResult {
    let mut new_wl = Vec::new();
    let mut new_ta = Vec::new();

    for (i, &w) in wl.iter().enumerate() {
        if w >= wl_min && w <= wl_max {
            let row_all_nan = ta[i].iter().all(|v| v.is_nan());
            if !row_all_nan {
                new_wl.push(w);
                new_ta.push(ta[i].clone());
            }
        }
    }

    CropResult {
        wavelength_array: new_wl,
        ta_2d: new_ta,
    }
}

pub fn baseline_subtraction(time: &[f64], ta: &[Vec<f64>], n_baseline: usize) -> Vec<Vec<f64>> {
    let mut t0_idx = 0;
    let mut min_abs = f64::INFINITY;
    for (i, &t) in time.iter().enumerate() {
        if t.abs() < min_abs {
            min_abs = t.abs();
            t0_idx = i;
        }
    }

    let start_idx = if t0_idx >= n_baseline {
        t0_idx - n_baseline
    } else if t0_idx > 0 {
        0
    } else {
        1
    };

    let end_idx = if t0_idx > 0 { t0_idx } else { 1 };

    ta.iter()
        .map(|row| {
            let mut sum = 0.0;
            let mut count = 0;
            for j in start_idx..end_idx {
                if j < row.len() && !row[j].is_nan() {
                    sum += row[j];
                    count += 1;
                }
            }
            let baseline = if count > 0 { sum / count as f64 } else { 0.0 };
            row.iter().map(|v| if v.is_nan() { f64::NAN } else { v - baseline }).collect()
        })
        .collect()
}

pub fn polyfit(x: &[f64], y: &[f64], order: usize) -> Vec<f64> {
    let n = x.len();
    let m = order + 1;

    if n < m {
        return vec![0.0; m];
    }

    let mut a = vec![vec![0.0; m]; m];
    let mut b = vec![0.0; m];

    for i in 0..m {
        for j in 0..m {
            let mut sum = 0.0;
            for k in 0..n {
                sum += x[k].powi((i + j) as i32);
            }
            a[i][j] = sum;
        }
        let mut sum = 0.0;
        for k in 0..n {
            sum += y[k] * x[k].powi(i as i32);
        }
        b[i] = sum;
    }

    for col in 0..m {
        let mut max_row = col;
        for row in (col + 1)..m {
            if a[row][col].abs() > a[max_row][col].abs() {
                max_row = row;
            }
        }
        a.swap(col, max_row);
        b.swap(col, max_row);

        let pivot = a[col][col];
        if pivot.abs() < 1e-14 {
            continue;
        }

        for j in col..m {
            a[col][j] /= pivot;
        }
        b[col] /= pivot;

        for row in 0..m {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            for j in col..m {
                a[row][j] -= factor * a[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    for i in 0..m {
        if a[i][i].abs() < 1e-14 {
            b[i] = 0.0;
        }
    }

    b
}

pub fn polyval(coeffs: &[f64], x: &[f64]) -> Vec<f64> {
    if coeffs.is_empty() {
        return vec![0.0; x.len()];
    }
    x.iter()
        .map(|&xi| {
            let mut val = 0.0;
            for k in (0..coeffs.len()).rev() {
                val = val * xi + coeffs[k];
            }
            val
        })
        .collect()
}

pub fn linspace(min: f64, max: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![min];
    }
    let step = (max - min) / (n - 1) as f64;
    (0..n).map(|i| min + i as f64 * step).collect()
}

pub fn interp1d(x_data: &[f64], y_data: &[f64], x_new: &[f64]) -> Vec<f64> {
    if x_data.len() < 2 {
        return vec![f64::NAN; x_new.len()];
    }
    x_new
        .iter()
        .map(|&xn| {
            if xn.is_nan() || xn < x_data[0] || xn > *x_data.last().unwrap() {
                return f64::NAN;
            }
            let mut lo = 0usize;
            let mut hi = x_data.len() - 1;
            while hi - lo > 1 {
                let mid = lo + (hi - lo) / 2;
                if x_data[mid] <= xn {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let denom = x_data[hi] - x_data[lo];
            if denom.abs() < 1e-15 {
                y_data[lo]
            } else {
                let t = (xn - x_data[lo]) / denom;
                y_data[lo] + t * (y_data[hi] - y_data[lo])
            }
        })
        .collect()
}
