pub struct IrfResult {
    pub ta_2d: Vec<Vec<f64>>,
    pub irf_fwhm: f64,
}

pub struct IrfEstimate {
    pub fwhm: f64,
    pub sigma: f64,
}

/// Estimate IRF width from the coherent artifact near t=0.
/// Averages signal across central wavelengths (strongest SNR),
/// then fits the derivative of the rising edge to a Gaussian.
pub fn estimate_irf(
    time: &[f64],
    ta: &[Vec<f64>],
    n_wl_avg: usize,
) -> Option<IrfEstimate> {
    if time.len() < 5 || ta.is_empty() {
        return None;
    }

    let n_wl = ta.len();
    let n_time = time.len();

    // Pick central wavelengths (strongest SNR region)
    let half = n_wl_avg / 2;
    let wl_start = if n_wl / 2 > half { n_wl / 2 - half } else { 0 };
    let wl_end = (n_wl / 2 + half).min(n_wl);

    // Average signal across selected wavelengths
    let mut avg_signal = vec![0.0; n_time];
    let mut counts = vec![0usize; n_time];
    for i in wl_start..wl_end {
        for j in 0..n_time {
            if !ta[i][j].is_nan() {
                avg_signal[j] += ta[i][j];
                counts[j] += 1;
            }
        }
    }
    for j in 0..n_time {
        if counts[j] > 0 {
            avg_signal[j] /= counts[j] as f64;
        } else {
            avg_signal[j] = f64::NAN;
        }
    }

    // Find t=0 index
    let mut t0_idx = 0;
    let mut min_abs = f64::INFINITY;
    for (i, &t) in time.iter().enumerate() {
        if t.abs() < min_abs {
            min_abs = t.abs();
            t0_idx = i;
        }
    }

    // Compute numerical derivative of the averaged signal near t=0
    // Use a window around t=0 to capture the rising edge
    let window = if n_time > 40 { 20 } else { n_time / 2 };
    let start = if t0_idx > window { t0_idx - window } else { 0 };
    let end = (t0_idx + window).min(n_time);

    if end - start < 5 {
        return None;
    }

    let mut deriv_t = Vec::new();
    let mut deriv_y = Vec::new();
    for j in start..end - 1 {
        let dt = time[j + 1] - time[j];
        if dt.abs() < 1e-15 {
            continue;
        }
        let dy = avg_signal[j + 1] - avg_signal[j];
        if dy.is_nan() {
            continue;
        }
        deriv_t.push((time[j] + time[j + 1]) / 2.0);
        deriv_y.push(dy / dt);
    }

    if deriv_t.len() < 5 {
        return None;
    }

    // The derivative of the rising edge should look like a Gaussian.
    // Find the peak of the absolute derivative
    let mut peak_idx = 0;
    let mut peak_val = 0.0;
    for (k, &v) in deriv_y.iter().enumerate() {
        if v.abs() > peak_val {
            peak_val = v.abs();
            peak_idx = k;
        }
    }

    if peak_val < 1e-10 {
        return None;
    }

    // Fit the derivative to a Gaussian: g(t) = A * exp(-(t - mu)^2 / (2*sigma^2))
    // Use log-space linear fit on the half-max region
    let mu = deriv_t[peak_idx];
    let half_max = peak_val * peak_val / 2.0;

    // Collect points above half-max for sigma estimation
    let mut sum_inv_sq = 0.0;
    let mut sum_ln = 0.0;
    let mut n_pts = 0;

    for k in 0..deriv_y.len() {
        let val_sq = deriv_y[k] * deriv_y[k];
        if val_sq > half_max && val_sq > 0.0 {
            let dt = deriv_t[k] - mu;
            sum_inv_sq += dt * dt;
            sum_ln += (-val_sq / (peak_val * peak_val)).ln();
            n_pts += 1;
        }
    }

    if n_pts < 2 || sum_inv_sq < 1e-30 {
        // Fallback: estimate sigma from the width at half-max of the derivative
        let mut left_idx = peak_idx;
        while left_idx > 0 {
            if deriv_y[left_idx].abs() < peak_val * 0.5 {
                break;
            }
            left_idx -= 1;
        }
        let mut right_idx = peak_idx;
        while right_idx < deriv_y.len() - 1 {
            if deriv_y[right_idx].abs() < peak_val * 0.5 {
                break;
            }
            right_idx += 1;
        }
        let fwhm_deriv = (deriv_t[right_idx] - deriv_t[left_idx]).abs();
        if fwhm_deriv < 1e-15 {
            return None;
        }
        // Derivative of Gaussian has FWHM = sigma * 2.355 / sqrt(2)
        // because derivative is sigma * gaussian, so FWHM_deriv = FWHM_gauss
        // Actually: derivative of Gaussian is proportional to Gaussian with same sigma
        // So FWHM of derivative = 2.355 * sigma
        let sigma = fwhm_deriv / 2.355;
        let fwhm = 2.355 * sigma;
        return Some(IrfEstimate { fwhm, sigma });
    }

    // sigma from linear fit: ln(val^2/A^2) = -dt^2 / sigma^2
    // so sigma^2 = -sum(dt^2) / sum(ln(val^2/A^2))
    let sigma_sq = -sum_inv_sq / sum_ln;
    if sigma_sq <= 0.0 {
        return None;
    }
    let sigma = sigma_sq.sqrt();
    let fwhm = 2.355 * sigma;

    Some(IrfEstimate { fwhm, sigma })
}

/// 1D discrete convolution with zero-padding.
/// Kernel is normalized so it sums to 1 before convolution.
fn convolve_1d(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let k = kernel.len();
    if n == 0 || k == 0 {
        return vec![0.0; n];
    }

    // Normalize kernel
    let k_sum: f64 = kernel.iter().sum();
    let norm_kernel: Vec<f64> = if k_sum.abs() > 1e-15 {
        kernel.iter().map(|&v| v / k_sum).collect()
    } else {
        kernel.to_vec()
    };

    let half_k = k / 2;
    let mut result = vec![0.0; n];

    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..k {
            let si = i as isize + j as isize - half_k as isize;
            if si >= 0 && si < n as isize {
                sum += signal[si as usize] * norm_kernel[k - 1 - j];
            }
        }
        result[i] = sum;
    }

    result
}

/// Generate a Gaussian IRF kernel on the given time grid.
/// sigma = fwhm / 2.355
fn make_gaussian_kernel(time: &[f64], sigma: f64) -> Vec<f64> {
    let dt = if time.len() >= 2 {
        (time[time.len() - 1] - time[0]) / (time.len() - 1) as f64
    } else {
        1.0
    };

    // Kernel spans ±4*sigma
    let n_kernel = (4.0 * sigma / dt).ceil() as usize * 2 + 1;
    let n_kernel = n_kernel.max(3);

    let mut kernel = Vec::with_capacity(n_kernel);
    let center = (n_kernel as f64 - 1.0) / 2.0;
    for i in 0..n_kernel {
        let x = (i as f64 - center) * dt;
        kernel.push((-x * x / (2.0 * sigma * sigma)).exp())
    }

    kernel
}

/// Deconvolve IRF using Lucy-Richardson iterative method.
/// For each wavelength, iterates:
///   f_{n+1} = f_n * [h(-t) ⊗ (measured / (h ⊗ f_n))]
/// where h is the Gaussian IRF kernel.
pub fn deconvolve_irf(
    time: &[f64],
    _wl: &[f64],
    ta: &[Vec<f64>],
    irf_fwhm: f64,
    n_iter: usize,
) -> IrfResult {
    let n_time = time.len();
    if n_time == 0 || ta.is_empty() {
        return IrfResult {
            ta_2d: ta.to_vec(),
            irf_fwhm,
        };
    }

    let dt = if n_time >= 2 {
        (time[n_time - 1] - time[0]) / (n_time - 1) as f64
    } else {
        1.0
    };

    // If IRF FWHM is smaller than time step, deconvolution is meaningless
    if irf_fwhm < dt * 0.5 {
        return IrfResult {
            ta_2d: ta.to_vec(),
            irf_fwhm,
        };
    }

    let sigma = irf_fwhm / 2.355;
    let kernel = make_gaussian_kernel(time, sigma);

    // Flipped kernel for correlation (h(-t))
    let kernel_flipped: Vec<f64> = kernel.iter().rev().copied().collect();

    let n_iter = if n_iter == 0 { 15 } else { n_iter };

    let mut result = Vec::with_capacity(ta.len());

    for row in ta {
        // Check if row has enough valid data
        let valid_count = row.iter().filter(|v| !v.is_nan()).count();
        if valid_count < 5 {
            result.push(row.clone());
            continue;
        }

        // Replace NaN with 0 for convolution, track NaN mask
        let mut signal = vec![0.0; n_time];
        let mut nan_mask = vec![false; n_time];
        for j in 0..n_time {
            if row[j].is_nan() {
                nan_mask[j] = true;
            } else {
                signal[j] = row[j];
            }
        }

        // Initial estimate: the measured signal itself
        let mut f = signal.clone();

        // Lucy-Richardson iteration
        let mut diverged = false;
        for _ in 0..n_iter {
            // h ⊗ f_n
            let hf = convolve_1d(&f, &kernel);

            // measured / (h ⊗ f_n), with protection against division by zero
            let mut ratio = vec![0.0; n_time];
            for j in 0..n_time {
                if nan_mask[j] {
                    ratio[j] = 1.0; // neutral element for multiplication
                } else if hf[j].abs() > 1e-15 {
                    ratio[j] = signal[j] / hf[j];
                } else {
                    ratio[j] = 1.0;
                }
            }

            // h(-t) ⊗ ratio
            let corr = convolve_1d(&ratio, &kernel_flipped);

            // f_{n+1} = f_n * corr
            let mut f_new = vec![0.0; n_time];
            for j in 0..n_time {
                f_new[j] = f[j] * corr[j];
            }

            // Check for divergence (values blowing up)
            let max_val = f_new.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
            if max_val > 1e10 || max_val.is_nan() || max_val.is_infinite() {
                diverged = true;
                break;
            }

            f = f_new;
        }

        if diverged {
            result.push(row.clone());
        } else {
            // Restore NaN positions
            let mut deconvolved = f;
            for j in 0..n_time {
                if nan_mask[j] {
                    deconvolved[j] = f64::NAN;
                }
            }
            result.push(deconvolved);
        }
    }

    IrfResult {
        ta_2d: result,
        irf_fwhm,
    }
}
