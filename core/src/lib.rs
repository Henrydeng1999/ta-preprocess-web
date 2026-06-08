pub mod ufs;
pub mod csv_parser;
pub mod baseline;
pub mod chirp;
pub mod fitting;
pub mod irf;
pub mod cpm;

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaData {
    pub time_array: Vec<f64>,
    pub wavelength_array: Vec<f64>,
    #[serde(rename = "ta_2d")]
    pub ta_2d: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessResult {
    pub wavelength_array: Vec<f64>,
    #[serde(rename = "ta_after")]
    pub ta_after: Vec<Vec<f64>>,
    #[serde(rename = "ta_before_chirp")]
    pub ta_before_chirp: Vec<Vec<f64>>,
    pub chirp_coeffs: Option<Vec<f64>>,
    #[serde(rename = "t0_per_wl")]
    pub t0_per_wl: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FitResult {
    pub wavelength: f64,
    pub n_exp: usize,
    pub params: Vec<f64>,
    pub std_errs: Vec<f64>,
    pub r2: f64,
    #[serde(rename = "t_data")]
    pub t_data: Vec<f64>,
    #[serde(rename = "y_data")]
    pub y_data: Vec<f64>,
    #[serde(rename = "t_fit")]
    pub t_fit: Vec<f64>,
    #[serde(rename = "y_fit")]
    pub y_fit: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ChirpOutput {
    #[serde(rename = "ta_2d")]
    ta_2d: Vec<f64>,
    coeffs: Option<Vec<f64>>,
    #[serde(rename = "t0_per_wl")]
    t0_per_wl: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CropOutput {
    wavelength_array: Vec<f64>,
    #[serde(rename = "ta_2d")]
    ta_2d: Vec<f64>,
    n_wl: usize,
    n_time: usize,
}

#[wasm_bindgen]
pub fn greet() -> String {
    "TA-WASM core loaded!".to_string()
}

#[wasm_bindgen]
pub fn parse_csv_wasm(text: &str) -> JsValue {
    let data = csv_parser::parse_csv(text);
    serde_wasm_bindgen::to_value(&data).unwrap()
}

#[wasm_bindgen]
pub fn parse_ufs_wasm(buffer: &[u8]) -> JsValue {
    match ufs::parse_ufs_file(buffer) {
        Ok(data) => serde_wasm_bindgen::to_value(&data).unwrap(),
        Err(e) => {
            let err = serde_wasm_bindgen::to_value(&format!("UFS parse error: {}", e)).unwrap();
            err
        }
    }
}

#[wasm_bindgen]
pub fn crop_wavelength(wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize, wl_min: f64, wl_max: f64) -> JsValue {
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let result = baseline::crop_wavelength(wl, &ta_2d, wl_min, wl_max);
    let n_wl_out = result.wavelength_array.len();
    let output = CropOutput {
        wavelength_array: result.wavelength_array,
        ta_2d: flatten_2d(&result.ta_2d),
        n_wl: n_wl_out,
        n_time,
    };
    serde_wasm_bindgen::to_value(&output).unwrap()
}

#[wasm_bindgen]
pub fn baseline_subtraction(time: &[f64], ta: &[f64], n_wl: usize, n_time: usize, n_baseline: usize) -> Vec<f64> {
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let result = baseline::baseline_subtraction(time, &ta_2d, n_baseline);
    flatten_2d(&result)
}

#[wasm_bindgen]
pub fn chirp_correction_half_height(
    time: &[f64], wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize,
    search_min: f64, search_max: f64, poly_order: usize,
    snr_threshold: f64, n_iter: usize, n_sigma: f64, n_baseline: usize,
) -> JsValue {
    if wl.len() != n_wl || ta.len() != n_wl * n_time {
        let err = serde_json::json!({"error": format!("Dimension mismatch: wl={}, n_wl={}, ta={}, expected={}", wl.len(), n_wl, ta.len(), n_wl * n_time)});
        return serde_wasm_bindgen::to_value(&err).unwrap();
    }
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let opts = chirp::ChirpOpts {
        search_range: (search_min, search_max),
        poly_order,
        snr_threshold,
        n_iter,
        n_sigma,
        n_baseline,
    };
    let result = chirp::chirp_correction_half_height(time, wl, &ta_2d, &opts);
    let output = ChirpOutput {
        ta_2d: flatten_2d(&result.ta_2d),
        coeffs: result.coeffs,
        t0_per_wl: result.t0_per_wl,
    };
    serde_wasm_bindgen::to_value(&output).unwrap()
}

#[wasm_bindgen]
pub fn chirp_correction_global(
    time: &[f64], wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize,
    search_min: f64, search_max: f64, poly_order: usize,
    snr_threshold: f64, n_iter: usize, n_sigma: f64, n_baseline: usize,
) -> JsValue {
    if wl.len() != n_wl || ta.len() != n_wl * n_time {
        let err = serde_json::json!({"error": format!("Dimension mismatch: wl={}, n_wl={}, ta={}, expected={}", wl.len(), n_wl, ta.len(), n_wl * n_time)});
        return serde_wasm_bindgen::to_value(&err).unwrap();
    }
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let opts = chirp::ChirpOpts {
        search_range: (search_min, search_max),
        poly_order,
        snr_threshold,
        n_iter,
        n_sigma,
        n_baseline,
    };
    let result = chirp::chirp_correction_global(time, wl, &ta_2d, &opts);
    let output = ChirpOutput {
        ta_2d: flatten_2d(&result.ta_2d),
        coeffs: result.coeffs,
        t0_per_wl: result.t0_per_wl,
    };
    serde_wasm_bindgen::to_value(&output).unwrap()
}

#[wasm_bindgen]
pub fn apply_chirp_with_coeffs(time: &[f64], wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize, coeffs: &[f64]) -> Vec<f64> {
    if wl.len() != n_wl || ta.len() != n_wl * n_time {
        return vec![];
    }
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let result = chirp::apply_chirp_shift(time, wl, &ta_2d, coeffs);
    flatten_2d(&result.ta_2d)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CpmOutput {
    #[serde(rename = "ta_2d")]
    ta_2d: Vec<f64>,
    amplitudes: Vec<f64>,
    sigmas: Vec<f64>,
    deriv_amplitudes: Vec<f64>,
    baselines: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct IrfOutput {
    #[serde(rename = "ta_2d")]
    ta_2d: Vec<f64>,
    irf_fwhm: f64,
}

#[wasm_bindgen]
pub fn remove_cpm_wasm(
    time: &[f64], wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize,
    fit_window: f64, n_baseline: usize,
) -> JsValue {
    if wl.len() != n_wl || ta.len() != n_wl * n_time {
        let err = serde_json::json!({"error": "Dimension mismatch"});
        return serde_wasm_bindgen::to_value(&err).unwrap();
    }
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let result = cpm::remove_cpm(time, wl, &ta_2d, fit_window, n_baseline);
    let output = CpmOutput {
        ta_2d: flatten_2d(&result.ta_2d),
        amplitudes: result.cpm_params.iter().map(|p| p.amplitude).collect(),
        sigmas: result.cpm_params.iter().map(|p| p.sigma).collect(),
        deriv_amplitudes: result.cpm_params.iter().map(|p| p.deriv_amplitude).collect(),
        baselines: result.cpm_params.iter().map(|p| p.baseline).collect(),
    };
    serde_wasm_bindgen::to_value(&output).unwrap()
}

#[wasm_bindgen]
pub fn estimate_irf_wasm(
    time: &[f64], wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize,
    n_wl_avg: usize,
) -> JsValue {
    if wl.len() != n_wl || ta.len() != n_wl * n_time {
        let err = serde_json::json!({"error": "Dimension mismatch"});
        return serde_wasm_bindgen::to_value(&err).unwrap();
    }
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    match irf::estimate_irf(time, &ta_2d, n_wl_avg) {
        Some(est) => {
            let output = serde_json::json!({"fwhm": est.fwhm, "sigma": est.sigma});
            serde_wasm_bindgen::to_value(&output).unwrap()
        }
        None => JsValue::NULL,
    }
}

#[wasm_bindgen]
pub fn deconvolve_irf_wasm(
    time: &[f64], wl: &[f64], ta: &[f64], n_wl: usize, n_time: usize,
    irf_fwhm: f64, n_iter: usize,
) -> JsValue {
    if wl.len() != n_wl || ta.len() != n_wl * n_time {
        let err = serde_json::json!({"error": "Dimension mismatch"});
        return serde_wasm_bindgen::to_value(&err).unwrap();
    }
    let ta_2d = flatten_to_2d(ta, n_wl, n_time);
    let result = irf::deconvolve_irf(time, wl, &ta_2d, irf_fwhm, n_iter);
    let output = IrfOutput {
        ta_2d: flatten_2d(&result.ta_2d),
        irf_fwhm: result.irf_fwhm,
    };
    serde_wasm_bindgen::to_value(&output).unwrap()
}

#[wasm_bindgen]
pub fn fit_multi_exp(time: &[f64], signal: &[f64], n_exp: usize, t_fit_min: f64, t_fit_max: f64) -> JsValue {
    let result = fitting::fit_multi_exp(time, signal, n_exp, t_fit_min, t_fit_max);
    match result {
        Some(r) => serde_wasm_bindgen::to_value(&r).unwrap(),
        None => JsValue::NULL,
    }
}

#[wasm_bindgen]
pub fn polyfit_wasm(x: &[f64], y: &[f64], order: usize) -> Vec<f64> {
    baseline::polyfit(x, y, order)
}

#[wasm_bindgen]
pub fn polyval_wasm(coeffs: &[f64], x: &[f64]) -> Vec<f64> {
    baseline::polyval(coeffs, x)
}

fn flatten_to_2d(flat: &[f64], n_wl: usize, n_time: usize) -> Vec<Vec<f64>> {
    let mut result = Vec::with_capacity(n_wl);
    for i in 0..n_wl {
        let row = flat[i * n_time..(i + 1) * n_time].to_vec();
        result.push(row);
    }
    result
}

fn flatten_2d(data: &[Vec<f64>]) -> Vec<f64> {
    let mut flat = Vec::with_capacity(data.len() * data.first().map_or(0, |r| r.len()));
    for row in data {
        flat.extend_from_slice(row);
    }
    flat
}
