use crate::baseline::linspace;
use crate::FitResult;

// ==================== Nelder-Mead (for chirp correction) ====================

pub struct NelderMeadResult {
    pub x: Vec<f64>,
    pub f_val: f64,
    pub iterations: usize,
}

pub fn nelder_mead(cost_func: &dyn Fn(&[f64]) -> f64, x0: &[f64], max_iter: usize) -> NelderMeadResult {
    let n = x0.len();
    let alpha = 1.0;
    let gamma = 2.0;
    let rho = 0.5;
    let sigma = 0.5;

    let mut simplex: Vec<Vec<f64>> = vec![x0.to_vec()];
    let mut f_vals = vec![cost_func(x0)];

    for i in 0..n {
        let mut xi = x0.to_vec();
        let step = if xi[i].abs() > 1e-4 {
            xi[i].abs() * 0.05
        } else {
            0.1
        };
        xi[i] += step;
        simplex.push(xi.clone());
        f_vals.push(cost_func(&xi));
    }

    let mut iter = 0;
    while iter < max_iter {
        let mut order: Vec<usize> = (0..=n).collect();
        order.sort_by(|&a, &b| {
            f_vals[a].partial_cmp(&f_vals[b]).unwrap_or(std::cmp::Ordering::Equal)
        });

        let best = order[0];
        let worst = order[n];
        let second_worst = order[n - 1];

        let x_best = simplex[best].clone();
        let f_best = f_vals[best];
        let f_second = f_vals[second_worst];
        let f_worst = f_vals[worst];

        let centroid: Vec<f64> = (0..n)
            .map(|j| {
                let sum: f64 = simplex
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != worst)
                    .map(|(_, s)| s[j])
                    .sum();
                sum / n as f64
            })
            .collect();

        let x_ref: Vec<f64> = centroid
            .iter()
            .zip(simplex[worst].iter())
            .map(|(&c, &w)| c + alpha * (c - w))
            .collect();
        let f_ref = cost_func(&x_ref);

        let mut shrink = false;

        if f_ref < f_second && f_ref >= f_best {
            simplex[worst] = x_ref;
            f_vals[worst] = f_ref;
        } else if f_ref < f_best {
            let x_exp: Vec<f64> = centroid
                .iter()
                .zip(x_ref.iter())
                .map(|(&c, &r)| c + gamma * (r - c))
                .collect();
            let f_exp = cost_func(&x_exp);
            if f_exp < f_ref {
                simplex[worst] = x_exp;
                f_vals[worst] = f_exp;
            } else {
                simplex[worst] = x_ref;
                f_vals[worst] = f_ref;
            }
        } else {
            let x_con: Vec<f64> = centroid
                .iter()
                .zip(simplex[worst].iter())
                .map(|(&c, &w)| c + rho * (w - c))
                .collect();
            let f_con = cost_func(&x_con);
            if f_con < f_worst {
                simplex[worst] = x_con;
                f_vals[worst] = f_con;
            } else {
                shrink = true;
            }
        }

        if shrink {
            for i in 1..=n {
                simplex[i] = x_best
                    .iter()
                    .zip(simplex[i].iter())
                    .map(|(&b, &s)| b + sigma * (s - b))
                    .collect();
                f_vals[i] = cost_func(&simplex[i]);
            }
        }

        let f_max = f_vals.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let f_min = f_vals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let f_range = (f_max - f_min).abs();
        let f_scale = f_max.abs().max(f_min.abs()).max(1.0);
        if f_range / f_scale < 1e-10 {
            break;
        }

        iter += 1;
    }

    let mut best_idx = 0;
    for i in 1..=n {
        if f_vals[i] < f_vals[best_idx] {
            best_idx = i;
        }
    }

    NelderMeadResult {
        x: simplex[best_idx].clone(),
        f_val: f_vals[best_idx],
        iterations: iter,
    }
}

// ==================== Model Functions ====================

/// Evaluate model with actual tau values: y(t) = y0 + sum(Ai * exp(-t/taui))
fn eval_model(params: &[f64], t: f64, n_exp: usize) -> f64 {
    let mut y = params[params.len() - 1];
    for k in 0..n_exp {
        let amp = params[k * 2];
        let tau = params[k * 2 + 1];
        if t >= 0.0 && tau > 1e-12 {
            y += amp * (-t / tau).exp();
        } else if t >= 0.0 {
            y += amp;
        }
    }
    y
}

/// Evaluate model with log(tau) parameterization for optimization
/// params layout: [A1, log_tau1, A2, log_tau2, ..., y0]
fn eval_model_log(params: &[f64], t: f64, n_exp: usize) -> f64 {
    let mut y = params[params.len() - 1];
    for k in 0..n_exp {
        let amp = params[k * 2];
        let tau = params[k * 2 + 1].exp();
        if t >= 0.0 {
            y += amp * (-t / tau).exp();
        }
    }
    y
}

fn compute_residuals_log(params: &[f64], t: &[f64], y: &[f64], n_exp: usize) -> Vec<f64> {
    t.iter()
        .zip(y.iter())
        .map(|(&ti, &yi)| eval_model_log(params, ti, n_exp) - yi)
        .collect()
}

fn compute_jacobian_log(params: &[f64], t: &[f64], n_exp: usize, n_p: usize) -> Vec<Vec<f64>> {
    let n = t.len();
    let mut jac = vec![vec![0.0; n_p]; n];
    let h = 1e-7;

    for j in 0..n_p {
        let mut p_plus = params.to_vec();
        p_plus[j] += h;

        for i in 0..n {
            let f_plus = eval_model_log(&p_plus, t[i], n_exp);
            let f_base = eval_model_log(params, t[i], n_exp);
            jac[i][j] = (f_plus - f_base) / h;
        }
    }

    jac
}

// ==================== Matrix Operations ====================

fn transpose(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if m.is_empty() {
        return Vec::new();
    }
    let rows = m.len();
    let cols = m[0].len();
    (0..cols)
        .map(|j| (0..rows).map(|i| m[i][j]).collect())
        .collect()
}

fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = b[0].len();
    let p = b.len();
    let mut result = vec![vec![0.0; m]; n];
    for i in 0..n {
        for j in 0..m {
            let mut sum = 0.0;
            for k in 0..p {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }
    result
}

fn mat_vec_mul(a: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
        .collect()
}

fn add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn scale(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn invert_matrix(m: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = m.len();
    if n == 0 {
        return None;
    }
    let mut a: Vec<Vec<f64>> = m
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .enumerate()
                .map(|(j, &v)| v + if i == j { 1e-10 } else { 0.0 })
                .collect()
        })
        .collect();

    let mut inv: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0; n];
            row[i] = 1.0;
            row
        })
        .collect();

    for col in 0..n {
        let mut max_row = col;
        for row in (col + 1)..n {
            if a[row][col].abs() > a[max_row][col].abs() {
                max_row = row;
            }
        }
        a.swap(col, max_row);
        inv.swap(col, max_row);

        if a[col][col].abs() < 1e-14 {
            return None;
        }

        let pivot = a[col][col];
        for j in 0..n {
            a[col][j] /= pivot;
            inv[col][j] /= pivot;
        }

        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            for j in 0..n {
                a[row][j] -= factor * a[col][j];
                inv[row][j] -= factor * inv[col][j];
            }
        }
    }

    Some(inv)
}

// ==================== Levenberg-Marquardt ====================

/// LM optimizer using log(tau) parameterization
/// Returns (params_in_log_space, chi2, std_errs_in_log_space)
fn levenberg_marquardt_log(
    t: &[f64],
    y: &[f64],
    n_exp: usize,
    initial_params: &[f64],
    max_iter: usize,
) -> Option<(Vec<f64>, f64, Vec<f64>)> {
    let n_p = n_exp * 2 + 1;
    let n_d = t.len();

    if n_d < n_p {
        return None;
    }

    let mut params = initial_params.to_vec();
    let mut lambda = 0.001;
    let lambda_up = 10.0;
    let lambda_down = 0.1;
    let max_lambda = 1e12;

    let mut residuals = compute_residuals_log(&params, t, y, n_exp);
    let mut chi2 = residuals.iter().map(|r| r * r).sum::<f64>();

    for _iter in 0..max_iter {
        let jac = compute_jacobian_log(&params, t, n_exp, n_p);
        let jac_t = transpose(&jac);
        let jtj = mat_mul(&jac_t, &jac);
        let jtr = mat_vec_mul(&jac_t, &residuals);

        let mut jtj_damped = jtj.clone();
        for i in 0..n_p {
            jtj_damped[i][i] += lambda;
        }

        let inv = match invert_matrix(&jtj_damped) {
            Some(m) => m,
            None => {
                lambda *= lambda_up;
                continue;
            }
        };
        let delta = mat_vec_mul(&inv, &scale(&jtr, -1.0));

        let new_params = add(&params, &delta);

        let new_residuals = compute_residuals_log(&new_params, t, y, n_exp);
        let new_chi2 = new_residuals.iter().map(|r| r * r).sum::<f64>();

        if new_chi2 < chi2 {
            let improvement = (chi2 - new_chi2) / (chi2 + 1e-30);
            params = new_params;
            residuals = new_residuals;
            chi2 = new_chi2;

            if improvement > 0.01 {
                lambda *= lambda_down;
            }
        } else {
            lambda *= lambda_up;
        }

        let delta_norm = norm(&delta);
        if delta_norm < 1e-12 || lambda > max_lambda {
            break;
        }
    }

    let final_jac = compute_jacobian_log(&params, t, n_exp, n_p);
    let jac_t = transpose(&final_jac);
    let jtj = mat_mul(&jac_t, &final_jac);
    let cov = invert_matrix(&jtj);

    let dof = (n_d as f64) - (n_p as f64);
    let sigma2 = if dof > 0.0 { chi2 / dof } else { 0.0 };

    let std_errs = if let Some(cov) = cov {
        (0..n_p)
            .map(|i| (sigma2 * cov[i][i]).max(0.0).sqrt())
            .collect()
    } else {
        vec![0.0; n_p]
    };

    Some((params, chi2, std_errs))
}

// ==================== Initial Guess Strategy ====================

/// Fit single exponential and return (A, log_tau, y0)
fn fit_single_exp_initial(t: &[f64], y: &[f64]) -> Vec<f64> {
    let t_max = t[t.len() - 1];
    let t_min = t[0].max(0.01);
    let log_t_min = t_min.ln();
    let log_t_max = t_max.ln();
    let log_range = log_t_max - log_t_min;

    let y_first = y[0];
    let y_last = y[y.len() - 1];
    let y_range = y.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
        - y.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let amp_est = y_first - y_last;

    // Strategy A: First/last points
    let guess_a = vec![amp_est, log_t_min + 0.5 * log_range, y_last * 0.5];

    // Strategy B: Area-based tau estimate
    let mut area_sum = 0.0;
    for i in 1..t.len() {
        area_sum += (y[i] - y_last) * (t[i] - t[i - 1]);
    }
    let tau_from_area = if amp_est.abs() > 1e-15 {
        (area_sum / amp_est).abs().clamp(t_min, t_max)
    } else {
        (t_max - t_min) * 0.3
    };
    let guess_b = vec![amp_est, tau_from_area.ln(), y_last * 0.5];

    // Strategy C: Half-decay point
    let half_val = (y_first + y_last) / 2.0;
    let mut t_half = t_min + (t_max - t_min) * 0.3;
    for i in 1..t.len() {
        if (y[i - 1] - half_val) * (y[i] - half_val) <= 0.0 {
            let denom = (y[i - 1] - half_val).abs() + (y[i] - half_val).abs() + 1e-30;
            let frac = (y[i - 1] - half_val).abs() / denom;
            t_half = t[i - 1] + frac * (t[i] - t[i - 1]);
            break;
        }
    }
    let tau_from_half = (t_half / std::f64::consts::LN_2).max(t_min);
    let guess_c = vec![amp_est, tau_from_half.ln(), y_last * 0.5];

    // Strategy D: Short tau
    let guess_d = vec![y_range * 0.8, log_t_min + 0.3 * log_range, y_last * 0.5];

    // Try all guesses with quick LM (NM too slow just for initial guess)
    let guesses = vec![guess_a, guess_b, guess_c, guess_d];
    let mut best = guesses[0].clone();
    let mut best_chi2 = f64::INFINITY;

    for x0 in &guesses {
        if let Some((params, chi2, _)) = levenberg_marquardt_log(t, y, 1, x0, 200) {
            if chi2 < best_chi2 {
                best_chi2 = chi2;
                best = params;
            }
        }
    }

    best
}

/// Returns (progressive_guesses, grid_guesses).
/// Progressive guesses are derived from sequential residual decomposition.
/// Grid guesses are uniform tau-spacing backups, including negated-amplitude variants.
fn generate_initial_guesses(
    t: &[f64],
    y: &[f64],
    n_exp: usize,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let t_max = t[t.len() - 1];
    let t_min = t[0].max(0.01);
    let log_t_min = t_min.ln();
    let log_t_max = t_max.ln();
    let log_range = log_t_max - log_t_min;

    let y_first = y[0];
    let y_last = y[y.len() - 1];
    let y_range = y.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
        - y.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let amp_est = y_first - y_last;

    let mut progressive = Vec::new();

    // ===== Progressive (Origin-style) guesses =====
    let single_params = fit_single_exp_initial(t, y);

    if n_exp == 1 {
        progressive.push(single_params.clone());
    } else {
        let residuals: Vec<f64> = t.iter()
            .enumerate()
            .map(|(i, &ti)| y[i] - eval_model_log(&single_params, ti, 1))
            .collect();

        let res_params = fit_single_exp_initial(t, &residuals);

        if n_exp >= 2 {
            let a1 = single_params[0];
            let lt1 = single_params[1];
            let a2 = res_params[0];
            let lt2 = res_params[1];
            let y0 = single_params[2] + res_params[2];

            if lt1 <= lt2 {
                progressive.push(vec![a1, lt1, a2, lt2, y0]);
            } else {
                progressive.push(vec![a2, lt2, a1, lt1, y0]);
            }
            progressive.push(vec![a2, lt2, a1, lt1, y0]);
        }

        if n_exp >= 3 {
            let residuals2: Vec<f64> = t.iter()
                .enumerate()
                .map(|(i, &ti)| residuals[i] - eval_model_log(&res_params, ti, 1))
                .collect();
            let res2_range = residuals2.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
                - residuals2.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let res2_amp = residuals2.first().copied().unwrap_or(0.0)
                - residuals2.last().copied().unwrap_or(0.0);

            let third_log_taus = [
                log_t_min + 0.1 * log_range,
                log_t_min + 0.3 * log_range,
                log_t_min + 0.7 * log_range,
            ];

            for log_tau3 in &third_log_taus {
                let a1 = single_params[0];
                let lt1 = single_params[1];
                let a2 = res_params[0];
                let lt2 = res_params[1];
                // Accumulate y0 from all residual stages
                let y0 = single_params[2] + res_params[2];
                let a3 = if res2_amp.abs() > 1e-15 { res2_amp } else { res2_range * 0.3 };

                let mut components = vec![(a1, lt1), (a2, lt2), (a3, *log_tau3)];
                components.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

                let mut x0 = Vec::new();
                for (a, lt) in &components {
                    x0.push(*a);
                    x0.push(*lt);
                }
                x0.push(y0);
                progressive.push(x0);
            }
        }
    }

    // ===== Grid-based guesses as backup =====
    let tau_strategies: Vec<Vec<f64>> = match n_exp {
        1 => vec![vec![0.5]],
        2 => vec![
            vec![0.2, 0.8],
            vec![0.15, 0.65],
            vec![0.3, 0.85],
            vec![0.1, 0.5],
        ],
        3 => vec![
            vec![0.1, 0.4, 0.8],
            vec![0.15, 0.45, 0.75],
            vec![0.05, 0.3, 0.7],
        ],
        4 => vec![
            vec![0.08, 0.25, 0.55, 0.85],
            vec![0.1, 0.3, 0.6, 0.9],
        ],
        5 => vec![
            vec![0.05, 0.2, 0.4, 0.65, 0.9],
            vec![0.08, 0.25, 0.5, 0.75, 0.95],
        ],
        _ => {
            let fracs = (0..n_exp).map(|k| (k as f64 + 0.5) / n_exp as f64).collect();
            vec![fracs]
        }
    };

    let amp_strategies: Vec<Box<dyn Fn(usize) -> f64>> = vec![
        Box::new(move |k: usize| -> f64 {
            if k == 0 { amp_est } else { y_range / n_exp as f64 }
        }),
        Box::new(move |_: usize| -> f64 { y_range / n_exp as f64 }),
        Box::new(move |k: usize| -> f64 { y_range * (0.5_f64.powi(k as i32 + 1)) }),
    ];

    let mut grid = Vec::new();
    for tau_fracs in &tau_strategies {
        for amp_fn in &amp_strategies {
            let mut x0 = Vec::new();
            for k in 0..n_exp {
                x0.push(amp_fn(k));
                x0.push(log_t_min + tau_fracs[k] * log_range);
            }
            x0.push(y_last * 0.5);
            // Also push negated-amplitude variant for negative-dominant signals (GSB)
            let mut x0_neg = x0.clone();
            for k in 0..n_exp {
                x0_neg[k * 2] = -x0_neg[k * 2];
            }
            grid.push(x0);
            grid.push(x0_neg);
        }
    }

    (progressive, grid)
}

// ==================== Multi-Exponential Fitting ====================

/// Log-spaced time grid from t_min to t_max.
/// t_fit_min is always the signal peak (positive), so log spacing naturally
/// concentrates points near the start where dynamics are fastest.
fn make_fit_tgrid(t_min: f64, t_max: f64, n: usize) -> Vec<f64> {
    if n == 0 || t_max <= t_min || t_min <= 0.0 {
        return linspace(t_min, t_max, n);
    }
    let (log_lo, log_hi) = (t_min.ln(), t_max.ln());
    if log_hi <= log_lo {
        return linspace(t_min, t_max, n);
    }
    (0..n).map(|i| {
        (log_lo + (log_hi - log_lo) * i as f64 / (n - 1).max(1) as f64).exp()
    }).collect()
}

pub fn fit_multi_exp(
    time: &[f64],
    signal: &[f64],
    n_exp: usize,
    t_fit_min: f64,
    t_fit_max: f64,
) -> Option<FitResult> {
    let fit_idx: Vec<usize> = time
        .iter()
        .enumerate()
        .filter(|(j, &t)| t >= t_fit_min && t <= t_fit_max && !signal[*j].is_nan())
        .map(|(j, _)| j)
        .collect();

    let n_p = n_exp * 2 + 1;
    if fit_idx.len() <= n_p {
        return None;
    }

    let t_fit: Vec<f64> = fit_idx.iter().map(|&j| time[j]).collect();
    let s_fit: Vec<f64> = fit_idx.iter().map(|&j| signal[j]).collect();

    let s_max = s_fit.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    if s_max < 1e-15 {
        return None;
    }

    let (progressive, grid) = generate_initial_guesses(&t_fit, &s_fit, n_exp);

    let cost_fn = |p: &[f64]| -> f64 {
        let mut ss = 0.0;
        for i in 0..t_fit.len() {
            let r = eval_model_log(p, t_fit[i], n_exp) - s_fit[i];
            ss += r * r;
        }
        ss
    };

    let mut best_coarse: Option<Vec<f64>> = None;
    let mut best_chi2 = f64::INFINITY;

    // Progressive guesses: NM coarse search → LM refinement
    for x0 in &progressive {
        if x0.len() != n_p { continue; }
        let nm = nelder_mead(&cost_fn, x0, 1000);
        if let Some((params, chi2, _)) = levenberg_marquardt_log(&t_fit, &s_fit, n_exp, &nm.x, 300) {
            if chi2 < best_chi2 {
                best_chi2 = chi2;
                best_coarse = Some(params);
            }
        }
    }

    // Grid guesses: lighter NM → LM (fallback)
    for x0 in &grid {
        if x0.len() != n_p { continue; }
        let nm = nelder_mead(&cost_fn, x0, 400);
        if let Some((params, chi2, _)) = levenberg_marquardt_log(&t_fit, &s_fit, n_exp, &nm.x, 200) {
            if chi2 < best_chi2 {
                best_chi2 = chi2;
                best_coarse = Some(params);
            }
        }
    }

    // Final polish with full LM to get accurate std_errs
    let coarse = best_coarse?;
    let (log_params, final_chi2, std_errs_log) =
        levenberg_marquardt_log(&t_fit, &s_fit, n_exp, &coarse, 1000)?;

    // Convert from log(tau) to actual tau for output
    let mut sorted: Vec<(f64, f64, f64, f64)> = (0..n_exp)
        .map(|k| {
            let amp = log_params[k * 2];
            let tau = log_params[k * 2 + 1].exp();
            let amp_err = std_errs_log[k * 2];
            let tau_err = tau * std_errs_log[k * 2 + 1]; // d(tau)/d(log_tau) = tau
            (amp, tau, amp_err, tau_err)
        })
        .collect();
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut final_params = Vec::new();
    let mut final_errs = Vec::new();
    for k in 0..n_exp {
        final_params.push(sorted[k].0);
        final_params.push(sorted[k].1);
        final_errs.push(sorted[k].2);
        final_errs.push(sorted[k].3);
    }
    final_params.push(log_params[log_params.len() - 1]); // y0
    final_errs.push(std_errs_log[std_errs_log.len() - 1]);

    // Compute R² using actual tau model
    let s_mean: f64 = s_fit.iter().sum::<f64>() / s_fit.len() as f64;
    let ss_tot: f64 = s_fit.iter().map(|v| (v - s_mean).powi(2)).sum();
    let r2 = if ss_tot > 1e-15 {
        1.0 - final_chi2 / ss_tot
    } else {
        0.0
    };

    // Compute fit curve using actual tau model (NOT log model!)
    let t_fine = make_fit_tgrid(t_fit_min, t_fit_max, 500);
    let y_fine: Vec<f64> = t_fine
        .iter()
        .map(|&t| eval_model(&final_params, t, n_exp))
        .collect();

    Some(FitResult {
        wavelength: 0.0,
        n_exp,
        params: final_params,
        std_errs: final_errs,
        r2,
        t_data: t_fit,
        y_data: s_fit,
        t_fit: t_fine,
        y_fit: y_fine,
    })
}
