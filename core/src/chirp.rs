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

pub fn apply_chirp_shift(time: &[f64], wl: &[f64], ta: &[Vec<f64>], coeffs: &[f64]) -> ChirpShiftResult {
    let t0_fitted = polyval(coeffs, wl);
    let ref_t0 = t0_fitted.iter().sum::<f64>() / t0_fitted.len() as f64;

    let corrected: Vec<Vec<f64>> = ta
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let dt = t0_fitted[i] - ref_t0;
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
        ref_t0,
    }
}

pub fn chirp_correction_half_height(time: &[f64], wl: &[f64], ta: &[Vec<f64>]) -> ChirpResult {
    let t0_per_wl = find_t0_half_height(time, ta, (-2.0, 2.0));

    let mut valid_x = Vec::new();
    let mut valid_y = Vec::new();
    for (i, &w) in wl.iter().enumerate() {
        if !t0_per_wl[i].is_nan() {
            valid_x.push(w);
            valid_y.push(t0_per_wl[i]);
        }
    }

    if valid_x.len() < 3 {
        return ChirpResult {
            ta_2d: ta.to_vec(),
            coeffs: None,
            t0_per_wl,
        };
    }

    let coeffs = polyfit(&valid_x, &valid_y, 2);
    let result = apply_chirp_shift(time, wl, ta, &coeffs);

    ChirpResult {
        ta_2d: result.ta_2d,
        coeffs: Some(coeffs),
        t0_per_wl,
    }
}

pub fn chirp_correction_global(time: &[f64], wl: &[f64], ta: &[Vec<f64>]) -> ChirpResult {
    let t0_per_wl = find_t0_half_height(time, ta, (-2.0, 2.0));

    let mut valid_x = Vec::new();
    let mut valid_y = Vec::new();
    for (i, &w) in wl.iter().enumerate() {
        if !t0_per_wl[i].is_nan() {
            valid_x.push(w);
            valid_y.push(t0_per_wl[i]);
        }
    }

    if valid_x.len() < 3 {
        return ChirpResult {
            ta_2d: ta.to_vec(),
            coeffs: None,
            t0_per_wl,
        };
    }

    let initial_coeffs = polyfit(&valid_x, &valid_y, 2);

    let n_wl = wl.len();

    let precomp_valid: Vec<Option<PrecompData>> = (0..n_wl)
        .map(|i| {
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
        let c0 = coeffs[0];
        let c1 = coeffs[1];
        let c2 = coeffs[2];

        let mut ref_t0 = 0.0;
        let dt_arr: Vec<f64> = wl
            .iter()
            .map(|&xi| {
                let t0 = c0 + c1 * xi + c2 * xi * xi;
                ref_t0 += t0;
                t0
            })
            .collect();
        ref_t0 /= n_wl as f64;

        let max_shift = dt_arr
            .iter()
            .map(|d| (d - ref_t0).abs())
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

            let dt = dt_arr[i] - ref_t0;
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

        let d0 = coeffs[0] - initial_coeffs[0];
        let d1 = coeffs[1] - initial_coeffs[1];
        let d2 = coeffs[2] - initial_coeffs[2];
        let reg = 0.01 * (d0 * d0 + d1 * d1 + d2 * d2) / 3.0;

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
