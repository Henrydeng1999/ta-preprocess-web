#[allow(unused_imports)]
use crate::baseline::interp1d;
use crate::fitting::nelder_mead;

/// Fitted CPM parameters for a single wavelength.
#[derive(Debug, Clone)]
pub struct CpmParams {
    pub amplitude: f64,       // Gaussian amplitude (A)
    pub sigma: f64,           // Gaussian width in ps (σ)
    pub deriv_amplitude: f64, // Derivative-Gaussian amplitude (B)
    pub baseline: f64,        // Local baseline offset (C)
}

/// Result of CPM removal across all wavelengths.
#[derive(Debug, Clone)]
pub struct CpmResult {
    pub ta_2d: Vec<Vec<f64>>,      // TA data after CPM removal
    pub cpm_params: Vec<CpmParams>, // fitted CPM parameters per wavelength
}

/// Evaluate the CPM model: y(t) = A * exp(-t²/(2σ²)) + B * (-t/σ²) * exp(-t²/(2σ²)) + C
fn cpm_model(params: &[f64], t: f64) -> f64 {
    let a = params[0];
    let sigma = params[1];
    let b = params[2];
    let c = params[3];

    if sigma <= 0.0 {
        return c;
    }

    let t2 = t * t;
    let s2 = sigma * sigma;
    let gauss = (-t2 / (2.0 * s2)).exp();
    let deriv_gauss = (-t / s2) * gauss;

    a * gauss + b * deriv_gauss + c
}

/// Remove CPM coherent artifact from TA data.
///
/// For each wavelength, fits the signal near t=0 with:
///   y(t) = A·exp(-t²/(2σ²)) + B·(-t/σ²)·exp(-t²/(2σ²)) + C
///
/// Then subtracts only the CPM component (Gaussian + derivative Gaussian),
/// leaving the baseline C intact.
///
/// # Arguments
/// * `time` - time delay array (ps), assumed sorted
/// * `wl` - probe wavelength array
/// * `ta` - 2D TA data, ta[i][j] = ΔA at wavelength wl[i], time time[j]
/// * `fit_window` - half-width of fitting range around t=0 (ps), e.g. 0.5
/// * `n_baseline` - number of pre-t=0 points used for baseline estimate
pub fn remove_cpm(
    time: &[f64],
    _wl: &[f64],
    ta: &[Vec<f64>],
    fit_window: f64,
    n_baseline: usize,
) -> CpmResult {
    let n_time = time.len();
    let n_wl = ta.len();

    // Find the index closest to t=0
    let mut t0_idx = 0;
    let mut min_abs_t = f64::INFINITY;
    for (i, &t) in time.iter().enumerate() {
        if t.abs() < min_abs_t {
            min_abs_t = t;
            t0_idx = i;
        }
    }

    // Collect indices within [-fit_window, fit_window]
    let fit_indices: Vec<usize> = (0..n_time)
        .filter(|&j| time[j].abs() <= fit_window)
        .collect();

    // Collect indices for baseline estimation (t < 0 region)
    let baseline_indices: Vec<usize> = (0..n_time)
        .filter(|&j| time[j] < 0.0)
        .collect();

    let mut result_ta = Vec::with_capacity(n_wl);
    let mut result_params = Vec::with_capacity(n_wl);

    for i in 0..n_wl {
        let row = &ta[i];

        // Skip if row is all NaN
        let all_nan = row.iter().all(|v| v.is_nan());
        if all_nan {
            result_ta.push(row.clone());
            result_params.push(CpmParams {
                amplitude: 0.0,
                sigma: 0.0,
                deriv_amplitude: 0.0,
                baseline: 0.0,
            });
            continue;
        }

        // Estimate baseline from t < 0 region
        let baseline_est: f64 = {
            let mut sum = 0.0;
            let mut count = 0;
            // Use up to n_baseline points from the t<0 region
            for &j in &baseline_indices {
                if count >= n_baseline {
                    break;
                }
                if !row[j].is_nan() {
                    sum += row[j];
                    count += 1;
                }
            }
            if count > 0 { sum / count as f64 } else { 0.0 }
        };

        // Estimate peak amplitude at t=0
        let peak_est = if !row[t0_idx].is_nan() {
            row[t0_idx] - baseline_est
        } else {
            // Try nearby points
            let mut found = 0.0;
            for delta in 1..5 {
                if t0_idx + delta < n_time && !row[t0_idx + delta].is_nan() {
                    found = row[t0_idx + delta] - baseline_est;
                    break;
                }
                if t0_idx >= delta && !row[t0_idx - delta].is_nan() {
                    found = row[t0_idx - delta] - baseline_est;
                    break;
                }
            }
            found
        };

        // Extract fitting data within the window (skip NaN values)
        let fit_t_clean: Vec<f64> = fit_indices.iter()
            .zip(fit_indices.iter().map(|&j| row[j]))
            .filter(|(_, v)| !v.is_nan())
            .map(|(&j, _)| time[j])
            .collect();
        let fit_y_clean: Vec<f64> = fit_indices.iter()
            .map(|&j| row[j])
            .filter(|v| !v.is_nan())
            .collect();

        if fit_t_clean.len() < 4 {
            // Not enough points to fit 4 parameters, leave unchanged
            result_ta.push(row.clone());
            result_params.push(CpmParams {
                amplitude: 0.0,
                sigma: 0.0,
                deriv_amplitude: 0.0,
                baseline: 0.0,
            });
            continue;
        }

        // Initial guess: [A, σ, B, C]
        let x0 = vec![peak_est, 0.05, 0.0, baseline_est];

        // Cost function: sum of squared residuals
        let cost = |params: &[f64]| -> f64 {
            let mut ssr = 0.0;
            for k in 0..fit_t_clean.len() {
                let y_pred = cpm_model(params, fit_t_clean[k]);
                let residual = y_pred - fit_y_clean[k];
                ssr += residual * residual;
            }
            // Penalty for unphysical sigma
            if params[1] <= 0.0 {
                ssr += 1e10;
            }
            ssr
        };

        let nm_result = nelder_mead(&cost, &x0, 2000);

        let fitted = nm_result.x;

        // Validate fitted sigma
        if fitted[1] <= 0.0 || fitted[1].is_nan() || fitted[1].is_infinite() {
            // Fit failed, leave unchanged
            result_ta.push(row.clone());
            result_params.push(CpmParams {
                amplitude: 0.0,
                sigma: 0.0,
                deriv_amplitude: 0.0,
                baseline: 0.0,
            });
            continue;
        }

        let a_fit = fitted[0];
        let sigma_fit = fitted[1];
        let b_fit = fitted[2];
        let c_fit = fitted[3];

        // Subtract CPM component from the full signal
        let new_row: Vec<f64> = row.iter().enumerate().map(|(j, &v)| {
            if v.is_nan() {
                f64::NAN
            } else {
                let cpm_component = a_fit * (-(time[j] * time[j]) / (2.0 * sigma_fit * sigma_fit)).exp()
                    + b_fit * (-time[j] / (sigma_fit * sigma_fit))
                        * (-(time[j] * time[j]) / (2.0 * sigma_fit * sigma_fit)).exp();
                v - cpm_component
            }
        }).collect();

        result_ta.push(new_row);
        result_params.push(CpmParams {
            amplitude: a_fit,
            sigma: sigma_fit,
            deriv_amplitude: b_fit,
            baseline: c_fit,
        });
    }

    CpmResult {
        ta_2d: result_ta,
        cpm_params: result_params,
    }
}
