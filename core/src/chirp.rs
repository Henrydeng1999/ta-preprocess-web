use crate::baseline::{interp1d, polyfit, polyval};

pub struct ChirpResult {
    pub ta_2d: Vec<Vec<f64>>,
    pub coeffs: Option<Vec<f64>>,
    pub t0_per_wl: Vec<f64>,
}

pub struct ChirpShiftResult {
    pub ta_2d: Vec<Vec<f64>>,
    pub t0_fitted: Vec<f64>,
    pub ref_t0: f64,
}

pub struct ChirpOpts {
    pub search_range: (f64, f64),
    pub poly_order: usize,
    pub snr_threshold: f64,
    pub n_iter: usize,
    pub n_sigma: f64,
    pub n_baseline: usize,
}

impl Default for ChirpOpts {
    fn default() -> Self {
        ChirpOpts {
            search_range: (-2.0, 2.0),
            poly_order: 2,
            snr_threshold: 3.0,
            n_iter: 3,
            n_sigma: 2.5,
            n_baseline: 10,
        }
    }
}

pub fn find_t0_half_height(time: &[f64], ta: &[Vec<f64>], search_range: (f64, f64)) -> Vec<f64> {
    let mut t_search_mask = Vec::new();
    let mut t_search = Vec::new();
    for (j, &t) in time.iter().enumerate() {
        if t >= search_range.0 && t <= search_range.1 {
            t_search_mask.push(j);
            t_search.push(t);
        }
    }

    let mut t0_list = Vec::with_capacity(ta.len());

    for row in ta {
        let sig: Vec<f64> = t_search_mask.iter().map(|&j| row[j]).collect();
        let valid_idx: Vec<usize> = sig
            .iter()
            .enumerate()
            .filter(|(_, &v)| !v.is_nan())
            .map(|(i, _)| i)
            .collect();

        if valid_idx.len() < 5 {
            t0_list.push(f64::NAN);
            continue;
        }

        let t_valid: Vec<f64> = valid_idx.iter().map(|&i| t_search[i]).collect();
        let s_valid: Vec<f64> = valid_idx.iter().map(|&i| sig[i]).collect();
        let abs_sig: Vec<f64> = s_valid.iter().map(|v| v.abs()).collect();

        let mut peak_idx = 0;
        let mut peak_val = 0.0;
        for (k, &v) in abs_sig.iter().enumerate() {
            if v > peak_val {
                peak_val = v;
                peak_idx = k;
            }
        }

        if peak_val < 1e-6 {
            t0_list.push(f64::NAN);
            continue;
        }

        let half_val = peak_val * 0.5;
        let search_start = peak_idx.saturating_sub(20);
        let mut cross_idx: Option<usize> = None;

        for k in search_start..peak_idx {
            if abs_sig[k] >= half_val {
                cross_idx = Some(k);
                break;
            }
        }

        let t0 = match cross_idx {
            None => t_valid[peak_idx],
            Some(ci) => {
                if ci > 0 && ci < t_valid.len() {
                    let t1 = t_valid[ci - 1];
                    let t2 = t_valid[ci];
                    let v1 = abs_sig[ci - 1];
                    let v2 = abs_sig[ci];
                    if (v2 - v1).abs() > 1e-10 {
                        t1 + (half_val - v1) / (v2 - v1) * (t2 - t1)
                    } else {
                        t_valid[ci]
                    }
                } else {
                    t_valid[ci]
                }
            }
        };
        t0_list.push(t0);
    }

    t0_list
}

/// Compute SNR per wavelength (same logic as JS fallback)
pub fn compute_snr_per_wl(time: &[f64], ta: &[Vec<f64>], n_baseline: usize) -> Vec<f64> {
    let mut t0_idx = 0;
    let mut min_abs = f64::INFINITY;
    for (i, &t) in time.iter().enumerate() {
        if t.abs() < min_abs {
            min_abs = t.abs();
            t0_idx = i;
        }
    }
    let baseline_start = if t0_idx >= n_baseline { t0_idx - n_baseline } else { 0 };

    let mut snr_list = Vec::with_capacity(ta.len());
    for row in ta {
        let mut noise_sum = 0.0;
        let mut noise_count = 0;
        for j in baseline_start..t0_idx {
            if j < row.len() && !row[j].is_nan() {
                noise_sum += row[j];
                noise_count += 1;
            }
        }
        let noise_mean = if noise_count > 0 { noise_sum / noise_count as f64 } else { 0.0 };
        let mut noise_var_sum = 0.0;
        for j in baseline_start..t0_idx {
            if j < row.len() && !row[j].is_nan() {
                noise_var_sum += (row[j] - noise_mean).powi(2);
            }
        }
        let noise_std = if noise_count > 1 { (noise_var_sum / (noise_count - 1) as f64).sqrt() } else { 1e-10 };
        let mut peak_val = 0.0;
        for v in row {
            if !v.is_nan() && v.abs() > peak_val {
                peak_val = v.abs();
            }
        }
        snr_list.push(if noise_std > 1e-15 { peak_val / noise_std } else { 0.0 });
    }
    snr_list
}

/// Sigma-clipped polyfit: iteratively remove outliers beyond n_sigma
pub fn sigma_clip_polyfit(
    x: &[f64],
    y: &[f64],
    order: usize,
    n_iter: usize,
    n_sigma: f64,
) -> (Option<Vec<f64>>, Vec<usize>, Vec<usize>) {
    let mut current_idx: Vec<usize> = (0..x.len())
        .filter(|&i| !x[i].is_nan() && !y[i].is_nan())
        .collect();

    if current_idx.len() < order + 1 {
        return (None, current_idx, vec![]);
    }

    let mut all_rejected: Vec<usize> = Vec::new();

    for _ in 0..n_iter {
        let cx: Vec<f64> = current_idx.iter().map(|&i| x[i]).collect();
        let cy: Vec<f64> = current_idx.iter().map(|&i| y[i]).collect();
        let iter_coeffs = polyfit(&cx, &cy, order);
        let y_fit = polyval(&iter_coeffs, &cx);
        let residuals: Vec<f64> = cy.iter().zip(y_fit.iter()).map(|(a, b)| a - b).collect();
        let res_mean = residuals.iter().sum::<f64>() / residuals.len() as f64;
        let res_std = {
            let var: f64 = residuals.iter().map(|r| (r - res_mean).powi(2)).sum();
            (var / (residuals.len().max(1) - 1).max(1) as f64).sqrt()
        };

        let mut new_idx = Vec::new();
        let mut iter_rejected = Vec::new();
        for (k, &idx) in current_idx.iter().enumerate() {
            if (residuals[k] - res_mean).abs() > n_sigma * res_std {
                iter_rejected.push(idx);
            } else {
                new_idx.push(idx);
            }
        }
        all_rejected.extend(iter_rejected);
        current_idx = new_idx;
        if current_idx.len() < order + 1 {
            break;
        }
    }

    let coeffs = if current_idx.len() >= order + 1 {
        let cx: Vec<f64> = current_idx.iter().map(|&i| x[i]).collect();
        let cy: Vec<f64> = current_idx.iter().map(|&i| y[i]).collect();
        Some(polyfit(&cx, &cy, order))
    } else {
        None
    };

    (coeffs, current_idx, all_rejected)
}

pub fn apply_chirp_shift(time: &[f64], wl: &[f64], ta: &[Vec<f64>], coeffs: &[f64]) -> ChirpShiftResult {
    // Use λ⁻² as x-axis for chirp (physical dispersion: t₀ ∝ λ⁻²)
    let inv_l2: Vec<f64> = wl.iter().map(|&w| 1.0 / (w * w)).collect();
    let t0_fitted = polyval(coeffs, &inv_l2);

    // Shift each wavelength so its t0 aligns to t=0
    let corrected: Vec<Vec<f64>> = ta
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let dt = t0_fitted[i];
            let time_shifted: Vec<f64> = time.iter().map(|&t| t - dt).collect();

            let mut valid_x = Vec::new();
            let mut valid_y = Vec::new();
            for (j, &v) in row.iter().enumerate() {
                if !v.is_nan() {
                    valid_x.push(time_shifted[j]);
                    valid_y.push(v);
                }
            }

            if valid_x.len() < 3 {
                return vec![f64::NAN; time.len()];
            }

            interp1d(&valid_x, &valid_y, time)
        })
        .collect();

    ChirpShiftResult {
        ta_2d: corrected,
        t0_fitted,
        ref_t0: 0.0,
    }
}

pub fn chirp_correction_half_height(time: &[f64], wl: &[f64], ta: &[Vec<f64>], opts: &ChirpOpts) -> ChirpResult {
    let t0_per_wl = find_t0_half_height(time, ta, opts.search_range);
    let snr_per_wl = compute_snr_per_wl(time, ta, opts.n_baseline);

    // Filter by SNR: use λ⁻² for fitting
    let inv_l2: Vec<f64> = wl.iter().map(|&w| 1.0 / (w * w)).collect();
    let mut valid_x = Vec::new();
    let mut valid_y = Vec::new();
    for (i, _) in wl.iter().enumerate() {
        if !t0_per_wl[i].is_nan() && snr_per_wl[i] >= opts.snr_threshold {
            valid_x.push(inv_l2[i]);
            valid_y.push(t0_per_wl[i]);
        }
    }

    if valid_x.len() < opts.poly_order + 1 {
        return ChirpResult {
            ta_2d: ta.to_vec(),
            coeffs: None,
            t0_per_wl,
        };
    }

    let (coeffs_opt, _kept, _rejected) = sigma_clip_polyfit(
        &valid_x, &valid_y, opts.poly_order, opts.n_iter, opts.n_sigma
    );

    match coeffs_opt {
        Some(coeffs) => {
            let result = apply_chirp_shift(time, wl, ta, &coeffs);
            ChirpResult {
                ta_2d: result.ta_2d,
                coeffs: Some(coeffs),
                t0_per_wl,
            }
        }
        None => ChirpResult {
            ta_2d: ta.to_vec(),
            coeffs: None,
            t0_per_wl,
        },
    }
}

pub fn chirp_correction_global(time: &[f64], wl: &[f64], ta: &[Vec<f64>], opts: &ChirpOpts) -> ChirpResult {
    let t0_per_wl = find_t0_half_height(time, ta, opts.search_range);
    let snr_per_wl = compute_snr_per_wl(time, ta, opts.n_baseline);

    // Filter by SNR: use λ⁻² for fitting
    let inv_l2: Vec<f64> = wl.iter().map(|&w| 1.0 / (w * w)).collect();
    let mut valid_x = Vec::new();
    let mut valid_y = Vec::new();
    for (i, _) in wl.iter().enumerate() {
        if !t0_per_wl[i].is_nan() && snr_per_wl[i] >= opts.snr_threshold {
            valid_x.push(inv_l2[i]);
            valid_y.push(t0_per_wl[i]);
        }
    }

    if valid_x.len() < opts.poly_order + 1 {
        return ChirpResult {
            ta_2d: ta.to_vec(),
            coeffs: None,
            t0_per_wl,
        };
    }

    let (initial_coeffs_opt, _kept, _rejected) = sigma_clip_polyfit(
        &valid_x, &valid_y, opts.poly_order, opts.n_iter, opts.n_sigma
    );

    let initial_coeffs = match initial_coeffs_opt {
        Some(c) => c,
        None => {
            return ChirpResult {
                ta_2d: ta.to_vec(),
                coeffs: None,
                t0_per_wl,
            };
        }
    };

    let n_wl = wl.len();
    let n_coeffs = initial_coeffs.len();

    // Precompute valid data per wavelength (skip SNR-filtered)
    let precomp_valid: Vec<Option<PrecompData>> = (0..n_wl)
        .map(|i| {
            if snr_per_wl[i] < opts.snr_threshold {
                return None;
            }
            let row = &ta[i];
            let mut vx = Vec::new();
            let mut vy = Vec::new();
            for (j, &v) in row.iter().enumerate() {
                if !v.is_nan() {
                    vx.push(time[j]);
                    vy.push(v);
                }
            }
            if vx.len() >= 5 {
                Some(PrecompData { vx, vy })
            } else {
                None
            }
        })
        .collect();

    let t_eval: Vec<f64> = time.iter().filter(|&&t| t >= -0.5 && t <= 0.5).copied().collect();
    let t_eval_len = t_eval.len();
    let dt_step = if t_eval_len >= 2 {
        (t_eval[t_eval_len - 1] - t_eval[0]) / (t_eval_len - 1) as f64
    } else {
        0.01
    };
    let inv_dt_step = 1.0 / dt_step;

    let cost_func = |coeffs: &[f64]| -> f64 {
        let dt_arr: Vec<f64> = wl
            .iter()
            .map(|&xi| {
                let inv_l2 = 1.0 / (xi * xi);
                let mut t0 = 0.0;
                let mut xi_pow = 1.0;
                for k in 0..n_coeffs {
                    t0 += coeffs[k] * xi_pow;
                    xi_pow *= inv_l2;
                }
                t0
            })
            .collect();

        let max_shift = dt_arr
            .iter()
            .map(|d| d.abs())
            .fold(0.0f64, f64::max);
        if max_shift > 5.0 {
            return 1e10 * (1.0 + max_shift);
        }

        let mut total_sharpness = 0.0;
        let mut n_valid_wl = 0;

        for (i, pv) in precomp_valid.iter().enumerate() {
            let pv = match pv {
                Some(p) => p,
                None => continue,
            };

            let dt = dt_arr[i];
            let vx = &pv.vx;
            let vy = &pv.vy;
            let vx_len = vx.len();
            let vx0 = vx[0];
            let vx_last = vx[vx_len - 1];

            let mut max_grad = 0.0;
            let mut prev_sig = f64::NAN;
            let mut valid_count = 0;

            for &te in &t_eval {
                let xn = te + dt;
                if xn < vx0 || xn > vx_last {
                    continue;
                }

                let mut lo = 0usize;
                let mut hi = vx_len - 1;
                while hi - lo > 1 {
                    let mid = lo + (hi - lo) / 2;
                    if vx[mid] <= xn {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }

                let denom = vx[hi] - vx[lo];
                let sig = if denom.abs() > 1e-15 {
                    vy[lo] + (xn - vx[lo]) / denom * (vy[hi] - vy[lo])
                } else {
                    vy[lo]
                };

                if valid_count > 0 {
                    let grad = (sig - prev_sig).abs() * inv_dt_step;
                    if grad > max_grad {
                        max_grad = grad;
                    }
                }
                prev_sig = sig;
                valid_count += 1;
            }

            if valid_count >= 2 {
                total_sharpness += max_grad;
                n_valid_wl += 1;
            }
        }

        if n_valid_wl == 0 {
            return 1e10;
        }

        let mut reg = 0.0;
        for k in 0..n_coeffs {
            reg += (coeffs[k] - initial_coeffs[k]).powi(2);
        }
        reg = 0.01 * reg / n_coeffs as f64;

        -total_sharpness / n_valid_wl as f64 + reg
    };

    let result = crate::fitting::nelder_mead(&cost_func, &initial_coeffs, 3000);
    let optimal_coeffs = result.x;

    let shift_result = apply_chirp_shift(time, wl, ta, &optimal_coeffs);

    ChirpResult {
        ta_2d: shift_result.ta_2d,
        coeffs: Some(optimal_coeffs),
        t0_per_wl,
    }
}

struct PrecompData {
    vx: Vec<f64>,
    vy: Vec<f64>,
}
