/* tslint:disable */
/* eslint-disable */

export function apply_chirp_with_coeffs(time: Float64Array, wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, coeffs: Float64Array): Float64Array;

export function baseline_subtraction(time: Float64Array, ta: Float64Array, n_wl: number, n_time: number, n_baseline: number): Float64Array;

export function chirp_correction_global(time: Float64Array, wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, search_min: number, search_max: number, poly_order: number, snr_threshold: number, n_iter: number, n_sigma: number, n_baseline: number): any;

export function chirp_correction_half_height(time: Float64Array, wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, search_min: number, search_max: number, poly_order: number, snr_threshold: number, n_iter: number, n_sigma: number, n_baseline: number): any;

export function crop_wavelength(wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, wl_min: number, wl_max: number): any;

export function deconvolve_irf_wasm(time: Float64Array, wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, irf_fwhm: number, n_iter: number): any;

export function estimate_irf_wasm(time: Float64Array, wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, n_wl_avg: number): any;

export function fit_multi_exp(time: Float64Array, signal: Float64Array, n_exp: number, t_fit_min: number, t_fit_max: number): any;

export function greet(): string;

export function parse_csv_wasm(text: string): any;

export function parse_ufs_wasm(buffer: Uint8Array): any;

export function polyfit_wasm(x: Float64Array, y: Float64Array, order: number): Float64Array;

export function polyval_wasm(coeffs: Float64Array, x: Float64Array): Float64Array;

export function remove_cpm_wasm(time: Float64Array, wl: Float64Array, ta: Float64Array, n_wl: number, n_time: number, fit_window: number, n_baseline: number): any;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly apply_chirp_with_coeffs: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => [number, number];
    readonly baseline_subtraction: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly chirp_correction_global: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number, n: number, o: number) => any;
    readonly chirp_correction_half_height: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number, n: number, o: number) => any;
    readonly crop_wavelength: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => any;
    readonly deconvolve_irf_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => any;
    readonly estimate_irf_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => any;
    readonly fit_multi_exp: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => any;
    readonly greet: () => [number, number];
    readonly parse_csv_wasm: (a: number, b: number) => any;
    readonly parse_ufs_wasm: (a: number, b: number) => any;
    readonly polyfit_wasm: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly polyval_wasm: (a: number, b: number, c: number, d: number) => [number, number];
    readonly remove_cpm_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => any;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
