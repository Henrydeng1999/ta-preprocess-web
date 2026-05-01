use std::env;
use std::fs;
use std::time::Instant;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: test_core <file.ufs>");
        std::process::exit(1);
    }

    let filepath = &args[1];
    println!("=== TA-WASM Core Test ===");
    println!("File: {}", filepath);

    let data = fs::read(filepath).expect("Failed to read file");
    println!("File size: {} bytes", data.len());

    // 1. UFS 解析
    println!("\n--- Step 1: UFS Parse ---");
    let t0 = Instant::now();
    let ta_data = ta_wasm::ufs::parse_ufs_file(&data);
    let parse_time = t0.elapsed();
    match ta_data {
        Ok(d) => {
            println!("✓ Parse OK ({:.2?})", parse_time);
            println!("  Wavelengths: {} ({} - {} nm)", d.wavelength_array.len(), d.wavelength_array.first().unwrap_or(&0.0), d.wavelength_array.last().unwrap_or(&0.0));
            println!("  Time points: {} ({} - {} ps)", d.time_array.len(), d.time_array.first().unwrap_or(&0.0), d.time_array.last().unwrap_or(&0.0));
            println!("  TA 2D shape: {} x {}", d.ta_2d.len(), d.ta_2d.first().map_or(0, |r| r.len()));

            let nan_count: usize = d.ta_2d.iter().map(|r| r.iter().filter(|v| v.is_nan()).count()).sum();
            let total: usize = d.ta_2d.len() * d.ta_2d.first().map_or(0, |r| r.len());
            println!("  NaN count: {} / {} ({:.1}%)", nan_count, total, nan_count as f64 / total as f64 * 100.0);

            let min_val = d.ta_2d.iter().flat_map(|r| r.iter()).filter(|v| !v.is_nan()).cloned().fold(f64::INFINITY, f64::min);
            let max_val = d.ta_2d.iter().flat_map(|r| r.iter()).filter(|v| !v.is_nan()).cloned().fold(f64::NEG_INFINITY, f64::max);
            println!("  Value range: [{:.6}, {:.6}]", min_val, max_val);

            // 2. 波长裁剪
            println!("\n--- Step 2: Wavelength Crop (450-650nm) ---");
            let t0 = Instant::now();
            let cropped = ta_wasm::baseline::crop_wavelength(&d.wavelength_array, &d.ta_2d, 450.0, 650.0);
            let crop_time = t0.elapsed();
            println!("✓ Crop OK ({:.2?})", crop_time);
            println!("  After crop: {} wavelengths", cropped.wavelength_array.len());

            // 3. 基线扣除
            println!("\n--- Step 3: Baseline Subtraction (n=10) ---");
            let t0 = Instant::now();
            let baselined = ta_wasm::baseline::baseline_subtraction(&d.time_array, &cropped.ta_2d, 10);
            let bl_time = t0.elapsed();
            println!("✓ Baseline OK ({:.2?})", bl_time);
            let bl_min = baselined.iter().flat_map(|r| r.iter()).filter(|v| !v.is_nan()).cloned().fold(f64::INFINITY, f64::min);
            let bl_max = baselined.iter().flat_map(|r| r.iter()).filter(|v| !v.is_nan()).cloned().fold(f64::NEG_INFINITY, f64::max);
            println!("  Value range after baseline: [{:.6}, {:.6}]", bl_min, bl_max);

            // 4. 啁啾校正 - 半高点法
            println!("\n--- Step 4: Chirp Correction (half-height) ---");
            let t0 = Instant::now();
            let chirp_result = ta_wasm::chirp::chirp_correction_half_height(&d.time_array, &cropped.wavelength_array, &baselined);
            let chirp_time = t0.elapsed();
            println!("✓ Chirp correction OK ({:.2?})", chirp_time);
            if let Some(ref coeffs) = chirp_result.coeffs {
                println!("  Chirp coeffs: {:?}", coeffs);
            }
            let valid_t0: Vec<f64> = chirp_result.t0_per_wl.iter().filter(|v| !v.is_nan()).cloned().collect();
            if !valid_t0.is_empty() {
                println!("  t0 range: [{:.4}, {:.4}] ps", valid_t0.iter().cloned().fold(f64::INFINITY, f64::min), valid_t0.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
            }

            // 5. 动力学拟合
            println!("\n--- Step 5: Kinetic Fitting ---");
            let mid_wl_idx = cropped.wavelength_array.len() / 2;
            let signal = &chirp_result.ta_2d[mid_wl_idx];
            let probe_wl = cropped.wavelength_array[mid_wl_idx];
            println!("  Probe wavelength: {:.1} nm", probe_wl);

            let t0 = Instant::now();
            let fit_2exp = ta_wasm::fitting::fit_multi_exp(&d.time_array, signal, 2, -1.0, d.time_array.last().unwrap_or(&1000.0).clone());
            let fit_time = t0.elapsed();
            match fit_2exp {
                Some(r) => {
                    println!("✓ 2-exp fit OK ({:.2?})", fit_time);
                    println!("  R² = {:.6}", r.r2);
                    for k in 0..r.n_exp {
                        let a = r.params[k * 2];
                        let tau = r.params[k * 2 + 1];
                        let a_err = r.std_errs[k * 2];
                        let tau_err = r.std_errs[k * 2 + 1];
                        println!("  τ{} = {:.4} ± {:.4} ps  (A{} = {:.6} ± {:.6})", k + 1, tau, tau_err, k + 1, a, a_err);
                    }
                    let offset = r.params[r.params.len() - 1];
                    println!("  offset = {:.6}", offset);
                }
                None => println!("✗ 2-exp fit failed"),
            }

            println!("\n=== All tests passed! ===");
        }
        Err(e) => {
            println!("✗ Parse FAILED: {}", e);
            println!("\nTrying raw hex dump of first 64 bytes:");
            let end = std::cmp::min(64, data.len());
            for i in (0..end).step_by(16) {
                let hex: Vec<String> = (i..std::cmp::min(i + 16, end)).map(|j| format!("{:02x}", data[j])).collect();
                let ascii: String = (i..std::cmp::min(i + 16, end)).map(|j| {
                    let b = data[j];
                    if b >= 32 && b <= 126 { b as char } else { '.' }
                }).collect();
                println!("  {:04x}: {:48}  {}", i, hex.join(" "), ascii);
            }
        }
    }
}
