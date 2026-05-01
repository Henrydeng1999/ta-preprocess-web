use crate::baseline::linspace;
use crate::FitResult;

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

fn multi_exp_func(params: &[f64], t: f64, n_exp: usize) -> f64 {
    let mut y = params[params.len() - 1];
    for k in 0..n_exp {
        let tau = params[k * 2 + 1];
        if tau.abs() < 1e-12 {
            y += params[k * 2];
        } else {
            y += params[k * 2] * (-t / tau).exp();
        }
    }
    y
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

    let tau_guesses = [1.0, 0.3, 5.0];

    let mut x0 = Vec::new();
    for k in 0..n_exp {
        x0.push(s_max / n_exp as f64 * if k == 0 { 1.0 } else { 0.5 });
        x0.push(tau_guesses[k]);
    }
    x0.push(0.0);

    let cost_func = |params: &[f64]| -> f64 {
        let mut ss = 0.0;
        for i in 0..t_fit.len() {
            let y_pred = multi_exp_func(params, t_fit[i], n_exp);
            ss += (y_pred - s_fit[i]).powi(2);
        }
        for k in 0..n_exp {
            if params[k * 2 + 1] < 0.001 {
                ss += 1e6;
            }
        }
        ss
    };

    let result = nelder_mead(&cost_func, &x0, 5000);
    let mut best_params = result.x;

    for k in 0..n_exp {
        if best_params[k * 2 + 1] < 0.0 {
            best_params[k * 2 + 1] = best_params[k * 2 + 1].abs();
            best_params[k * 2] = -best_params[k * 2];
        }
    }

    let mut sorted: Vec<(f64, f64)> = (0..n_exp)
        .map(|k| (best_params[k * 2], best_params[k * 2 + 1]))
        .collect();
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut final_params = Vec::new();
    for (a, tau) in &sorted {
        final_params.push(*a);
        final_params.push(*tau);
    }
    final_params.push(best_params[best_params.len() - 1]);

    let s_mean: f64 = s_fit.iter().sum::<f64>() / s_fit.len() as f64;
    let ss_tot: f64 = s_fit.iter().map(|v| (v - s_mean).powi(2)).sum();
    let r2 = if ss_tot > 1e-15 {
        1.0 - result.f_val / ss_tot
    } else {
        0.0
    };

    let n_d = t_fit.len();
    let dof = (n_d as f64) - (n_p as f64);
    let sigma2 = if dof > 0.0 {
        result.f_val / dof
    } else {
        0.0
    };

    let eps = 1e-8;
    let jacobian: Vec<Vec<f64>> = (0..n_d)
        .map(|i| {
            (0..n_p)
                .map(|j| {
                    let mut p_plus = final_params.clone();
                    p_plus[j] += eps;
                    let f_plus = multi_exp_func(&p_plus, t_fit[i], n_exp);
                    let f_orig = multi_exp_func(&final_params, t_fit[i], n_exp);
                    (f_plus - f_orig) / eps
                })
                .collect()
        })
        .collect();

    let jtj: Vec<Vec<f64>> = (0..n_p)
        .map(|i| {
            (0..n_p)
                .map(|j| {
                    let s: f64 = (0..n_d).map(|k| jacobian[k][i] * jacobian[k][j]).sum();
                    s
                })
                .collect()
        })
        .collect();

    let cov = invert_matrix(&jtj);
    let std_errs = if let Some(cov) = cov {
        (0..n_p)
            .map(|i| (sigma2 * cov[i][i]).max(0.0).sqrt())
            .collect()
    } else {
        vec![0.0; n_p]
    };

    let t_fine = linspace(t_fit_min, t_fit_max, 500);
    let y_fine: Vec<f64> = t_fine.iter().map(|&t| multi_exp_func(&final_params, t, n_exp)).collect();

    Some(FitResult {
        wavelength: 0.0,
        n_exp,
        params: final_params,
        std_errs,
        r2,
        t_data: t_fit,
        y_data: s_fit,
        t_fit: t_fine,
        y_fit: y_fine,
    })
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
