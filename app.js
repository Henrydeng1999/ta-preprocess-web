var uploadedFiles = [];
var ufsFiles = [];

var _wasmReady = false;

function escapeHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}

function _flattenTA(ta) {
  var flat = [];
  for (var i = 0; i < ta.length; i++) flat.push.apply(flat, ta[i]);
  return flat;
}

function _unflattenTA(flat, nWl, nTime) {
  if (!flat) return [];
  var arr = Array.isArray(flat) ? flat : Array.from(flat);
  var result = [];
  for (var i = 0; i < nWl; i++) {
    result.push(arr.slice(i * nTime, (i + 1) * nTime));
  }
  return result;
}

function _wasmChirpResult(ret, nWl, nTime) {
  if (!ret || ret.error) return null;
  var ta2d = ret.ta_2d;
  if (ta2d === undefined || ta2d === null) return null;
  return { TA2D: _unflattenTA(ta2d, nWl, nTime), coeffs: ret.coeffs, t0PerWl: ret.t0_per_wl };
}

(async function initWasm() {
  try {
    const mod = await import('./pkg/ta_wasm.js');
    await mod.default();
    window.taWasm = {
      parse_csv_wasm: mod.parse_csv_wasm,
      parse_ufs_wasm: mod.parse_ufs_wasm,
      crop_wavelength: mod.crop_wavelength,
      baseline_subtraction: mod.baseline_subtraction,
      chirp_correction_half_height: mod.chirp_correction_half_height,
      chirp_correction_global: mod.chirp_correction_global,
      apply_chirp_with_coeffs: mod.apply_chirp_with_coeffs,
      fit_multi_exp: mod.fit_multi_exp,
      polyfit_wasm: mod.polyfit_wasm,
      polyval_wasm: mod.polyval_wasm,
      remove_cpm_wasm: mod.remove_cpm_wasm,
      deconvolve_irf_wasm: mod.deconvolve_irf_wasm,
      greet: mod.greet
    };
    console.log('[WASM] TA core loaded:', mod.greet());
    window.dispatchEvent(new CustomEvent('wasm-ready'));
  } catch(e) {
    console.error('[WASM] Load failed:', e);
    const errDiv = document.createElement('div');
    errDiv.style.cssText = 'position:fixed;top:0;left:0;right:0;padding:12px 20px;background:#dc3545;color:#fff;z-index:9999;text-align:center;font-family:monospace;font-size:13px;';
    errDiv.textContent = '⚠️ WASM 核心加载失败，数据处理不可用。请检查浏览器是否支持 WebAssembly 并刷新页面。';
    document.body.prepend(errDiv);
  }
})();

window.addEventListener('wasm-ready', function() {
  _wasmReady = true;
});

function $(id) { return document.getElementById(id); }

// ============================================================
// UFS (Version2) Binary Parser
// ============================================================

function readUfsString(buffer, offset, length) {
  const arr = new Uint8Array(buffer, offset, length);
  let end = 0;
  while (end < length && arr[end] >= 32 && arr[end] <= 126) end++;
  let result = '';
  for (let k = 0; k < end; k++) result += String.fromCharCode(arr[k]);
  return { value: result, consumed: length };
}

function readUfsBeUint32(buffer, offset) {
  const view = new DataView(buffer);
  return view.getUint32(offset, false);
}

function readUfsBeUint16(buffer, offset) {
  const view = new DataView(buffer);
  return view.getUint16(offset, false);
}

function readUfsBeFloat64(buffer, offset) {
  const view = new DataView(buffer);
  return view.getFloat64(offset, false);
}

function parseUfsSectionHeader(buffer, offset) {
  // [4B name_len][name][2B pad][2B dtype][2B unit][4B count]
  const nameLen = readUfsBeUint32(buffer, offset);
  offset += 4;
  const nameResult = readUfsString(buffer, offset, nameLen);
  const name = nameResult.value;
  offset += nameLen;

  const pad1 = readUfsBeUint16(buffer, offset); // skip
  const dtype = readUfsBeUint16(buffer, offset + 2);
  offset += 4;

  const unit = String.fromCharCode(
    new Uint8Array(buffer, offset, 1)[0],
    new Uint8Array(buffer, offset + 1, 1)[0]
  ).replace(/\x00/g, '');
  offset += 2;

  const count = readUfsBeUint32(buffer, offset);
  offset += 4;

  return { name, dtype, unit, count, dataOffset: offset, headerStart: offset - nameLen - 12 };
}

function parseUfsIntensityHeader(buffer, offset) {
  // [4B name_len][name][2B pad][2B dtype][2B unit][2B wl_count][2B pad][2B t_count]
  const nameLen = readUfsBeUint32(buffer, offset);
  offset += 4;
  const nameResult = readUfsString(buffer, offset, nameLen);
  const name = nameResult.value;
  offset += nameLen;

  const pad1 = readUfsBeUint16(buffer, offset);
  const dtype = readUfsBeUint16(buffer, offset + 2);
  const unit = String.fromCharCode(
    new Uint8Array(buffer, offset + 4, 1)[0],
    new Uint8Array(buffer, offset + 5, 1)[0]
  ).replace(/\x00/g, '');
  offset += 6;

  const wlCount = readUfsBeUint16(buffer, offset);
  offset += 2;
  const pad2 = readUfsBeUint16(buffer, offset);
  offset += 2;
  const tCount = readUfsBeUint16(buffer, offset);
  offset += 2;

  return { name, dtype, unit, wlCount, tCount, dataOffset: offset };
}

function parseUfsFile(arrayBuffer) {
  const data = arrayBuffer;
  let offset = 0;

  // Global header: [4B len]["Version2"]
  const verLen = readUfsBeUint32(data, offset);
  offset += 4;
  const verResult = readUfsString(data, offset, verLen);
  const version = verResult.value;
  offset += verLen;

  // Wavelength section
  const wlSection = parseUfsSectionHeader(data, offset);
  const wavelengths = [];
  for (let i = 0; i < wlSection.count; i++) {
    wavelengths.push(readUfsBeFloat64(data, wlSection.dataOffset + i * 8));
  }
  offset = wlSection.dataOffset + wlSection.count * 8;

  // Time section
  const tSection = parseUfsSectionHeader(data, offset);
  const times = [];
  for (let i = 0; i < tSection.count; i++) {
    times.push(readUfsBeFloat64(data, tSection.dataOffset + i * 8));
  }
  offset = tSection.dataOffset + tSection.count * 8;

  // Intensity section
  const intSection = parseUfsIntensityHeader(data, offset);
  const wlCount = intSection.wlCount;
  const tCount = intSection.tCount;
  const intStart = intSection.dataOffset;

  const intensity = [];
  for (let wi = 0; wi < wlCount; wi++) {
    const row = [];
    for (let ti = 0; ti < tCount; ti++) {
      const flat = wi * tCount + ti;
      const off = intStart + flat * 8;
      let val = readUfsBeFloat64(data, off);
      if (isNaN(val) || !isFinite(val) || Math.abs(val) < 1e-300) {
        val = NaN;
      }
      row.push(val);
    }
    intensity.push(row);
  }

  // Extract text metadata from file tail
  const fileLen = data.byteLength;
  let textStart = fileLen;
  const uint8 = new Uint8Array(data);
  for (let i = fileLen - 1; i >= Math.max(0, fileLen - 10000); i--) {
    if (uint8[i] >= 32 && uint8[i] <= 126) {
      let j = i;
      while (j > 0 && ((uint8[j-1] >= 32 && uint8[j-1] <= 126) || [9,10,13].includes(uint8[j-1]))) j--;
      textStart = j;
      break;
    }
  }

  const metadata = {};
  if (textStart < fileLen) {
    const textBytes = uint8.slice(textStart);
    let text = '';
    for (let i = 0; i < textBytes.length; i++) {
      if (textBytes[i] === 0) break;
      if (textBytes[i] >= 32 && textBytes[i] <= 126) text += String.fromCharCode(textBytes[i]);
      else if ([9,10,13].includes(textBytes[i])) text += ' ';
    }
    text = text.replace(/\r\n/g, '\n').replace(/\r/g, '\n').trim();
    for (const line of text.split('\n')) {
      const colonIdx = line.indexOf(':');
      if (colonIdx > 0) {
        const key = line.slice(0, colonIdx).trim();
        const value = line.slice(colonIdx + 1).trim();
        if (key) metadata[key] = value;
      }
    }
  }

  // Add units from section headers
  if (tSection.unit) metadata['Time units'] = tSection.unit;
  if (intSection.name && intSection.name !== 'DA') {
    metadata['Z axis title'] = intSection.name;
  } else if (intSection.name === 'DA') {
    metadata['Z axis title'] = 'dA';
  }

  return { wavelengths, times, intensity, metadata, version };
}



function enterApp() {
  $('coverPage').classList.add('hidden');
}

document.addEventListener('DOMContentLoaded', function() {
  if (location.protocol === 'file:') {
    alert('⚠️ 请通过 HTTP 服务器打开本页面（file:// 协议下文件读取受限）。\n\n推荐：\n  npx serve .\n  或 python -m http.server 8080\n  或 VS Code Live Server 插件');
  }
  const enterBtn = document.getElementById('enterBtn');
  if (enterBtn) enterBtn.addEventListener('click', enterApp);

  // Show/hide manual fit order selector based on chirp method
  const chirpSel = $('chirpMethod');
  const orderGroup = $('manualFitOrderGroup');
  if (chirpSel && orderGroup) {
    chirpSel.addEventListener('change', () => {
      orderGroup.style.display = chirpSel.value === 'manual' ? '' : 'none';
    });
  }
});

$('uploadZone').addEventListener('click', () => $('fileInput').click());
$('uploadZone').addEventListener('dragover', e => { e.preventDefault(); $('uploadZone').classList.add('dragover'); });
$('uploadZone').addEventListener('dragleave', () => $('uploadZone').classList.remove('dragover'));
$('uploadZone').addEventListener('drop', e => {
  e.preventDefault();
  $('uploadZone').classList.remove('dragover');
  handleFiles(e.dataTransfer.files);
});
$('fileInput').addEventListener('change', e => handleFiles(e.target.files));

async function handleFiles(files) {
  for (const f of files) {
    if (uploadedFiles.find(u => u.name === f.name)) continue;

    // Primary: file extension check; Fallback: magic bytes
    const extMatch = f.name.toLowerCase().endsWith('.ufs');
    const isUfs = extMatch || await new Promise(resolve => {
      const reader = new FileReader();
      reader.onload = e => {
        try {
          const bytes = new Uint8Array(e.target.result);
          const len = (bytes[0] << 24 | bytes[1] << 16 | bytes[2] << 8 | bytes[3]) >>> 0;
          if (len < 1 || len > 100 || 4 + len > bytes.length) { resolve(false); return; }
          let str = '';
          for (let i = 4; i < 4 + len; i++) {
            if (bytes[i] < 32 || bytes[i] > 126) { resolve(false); return; }
            str += String.fromCharCode(bytes[i]);
          }
          resolve(str === 'Version2');
        } catch { resolve(false); }
      };
      reader.onerror = () => resolve(false);
      reader.readAsArrayBuffer(f.slice(0, 128));
    });

    uploadedFiles.push(f);
    if (isUfs) ufsFiles.push(f.name);
  }
  renderFileList();
  $('processBtn').disabled = uploadedFiles.length === 0;
}

function renderFileList() {
  let html = '';
  uploadedFiles.forEach((f, i) => {
    const sizeMB = (f.size / 1024 / 1024).toFixed(2);
    html += `<div class="file-item"><span class="name">📄 ${escapeHtml(f.name)}</span><span class="size">${sizeMB} MB</span><span><button class="btn" style="padding:4px 10px;font-size:12px;background:#dc3545;color:white;" onclick="removeFile(${i})">✕</button></span></div>`;
  });
  $('fileList').innerHTML = html;
}

function removeFile(idx) {
  const name = uploadedFiles[idx].name;
  uploadedFiles.splice(idx, 1);
  ufsFiles = ufsFiles.filter(n => n !== name);
  renderFileList();
  $('processBtn').disabled = uploadedFiles.length === 0;
}

function setStatus(msg, type, progress) {
  const bar = progress !== undefined
    ? `<div class="progress-bar"><div class="fill" style="width:${Math.min(100, Math.max(0, progress)).toFixed(1)}%"></div></div>`
    : '';
  $('statusArea').innerHTML = `<div class="status status-${type}">${msg}${bar}</div>`;
}


function cropWavelength(wl, ta, wlMin, wlMax) {
  const idx = [];
  const newWl = [];
  for (let i = 0; i < wl.length; i++) {
    if (wl[i] >= wlMin && wl[i] <= wlMax) {
      const rowAllNaN = ta[i].every(v => isNaN(v));
      if (!rowAllNaN) { idx.push(i); newWl.push(wl[i]); }
    }
  }
  const newTA = idx.map(i => [...ta[i]]);
  return { wavelengthArray: newWl, TA2D: newTA };
}

function baselineSubtraction(time, ta, nBaseline) {
  let t0Idx = 0;
  let minAbs = Infinity;
  for (let i = 0; i < time.length; i++) {
    if (Math.abs(time[i]) < minAbs) { minAbs = Math.abs(time[i]); t0Idx = i; }
  }
  const startIdx = Math.max(0, t0Idx - nBaseline);
  const nTimes = time.length;
  const result = ta.map(row => {
    let sum = 0, count = 0;
    for (let j = startIdx; j < t0Idx; j++) {
      if (!isNaN(row[j])) { sum += row[j]; count++; }
    }
    const baseline = count > 0 ? sum / count : 0;
    return row.map(v => isNaN(v) ? NaN : v - baseline);
  });
  return result;
}

function computeSnrPerWl(time, ta, nBaseline) {
  let t0Idx = 0;
  let minAbs = Infinity;
  for (let i = 0; i < time.length; i++) {
    if (Math.abs(time[i]) < minAbs) { minAbs = Math.abs(time[i]); t0Idx = i; }
  }
  const baselineStart = Math.max(0, t0Idx - nBaseline);
  const snrPerWl = [];
  for (let i = 0; i < ta.length; i++) {
    const row = ta[i];
    let noiseSum = 0, noiseCount = 0;
    for (let j = baselineStart; j < t0Idx; j++) {
      if (!isNaN(row[j])) { noiseSum += row[j]; noiseCount++; }
    }
    const noiseMean = noiseCount > 0 ? noiseSum / noiseCount : 0;
    let noiseVarSum = 0;
    for (let j = baselineStart; j < t0Idx; j++) {
      if (!isNaN(row[j])) noiseVarSum += (row[j] - noiseMean) ** 2;
    }
    const noiseStd = noiseCount > 1 ? Math.sqrt(noiseVarSum / (noiseCount - 1)) : 1e-10;
    let peakVal = 0;
    for (let j = 0; j < row.length; j++) {
      if (!isNaN(row[j]) && Math.abs(row[j]) > peakVal) peakVal = Math.abs(row[j]);
    }
    snrPerWl.push(noiseStd > 1e-15 ? peakVal / noiseStd : 0);
  }
  return snrPerWl;
}

function polyfit(x, y, order) {
  return Array.from(window.taWasm.polyfit_wasm(new Float64Array(x), new Float64Array(y), order));
}

function polyval(coeffs, x) {
  return Array.from(window.taWasm.polyval_wasm(new Float64Array(coeffs), new Float64Array(x)));
}

function linspace(min, max, n) {
  const step = (max - min) / (n - 1);
  return Array.from({ length: n }, (_, i) => min + i * step);
}


function applyChirpShift(time, wl, ta, coeffs) {
  var flat = _flattenTA(ta);
  var ret = window.taWasm.apply_chirp_with_coeffs(
    new Float64Array(time), new Float64Array(wl), new Float64Array(flat), ta.length, time.length, new Float64Array(coeffs)
  );
  if (ret && ret.length > 0) {
    var result = _unflattenTA(ret, ta.length, time.length);
    if (result && result.length > 0 && result[0].length > 0) return { TA2D: result };
  }
  return { TA2D: ta };
}

function chirpCorrectionHalfHeight(time, wl, ta, opts) {
  var flat = _flattenTA(ta);
  var ret = window.taWasm.chirp_correction_half_height(
    new Float64Array(time), new Float64Array(wl), new Float64Array(flat), ta.length, time.length,
    opts.searchRange[0], opts.searchRange[1], opts.polyOrder,
    opts.snrThreshold, opts.nIter, opts.nSigma, opts.nBaseline
  );
  return _wasmChirpResult(ret, ta.length, time.length) ||
    { TA2D: ta, coeffs: null, t0PerWl: new Array(ta.length).fill(NaN), snrPerWl: new Array(ta.length).fill(0), snrFilteredOut: new Array(ta.length).fill(true), sigmaClippedOut: new Array(ta.length).fill(false) };
}

async function chirpCorrectionGlobal(time, wl, ta, opts) {
  var flat = _flattenTA(ta);
  var ret = await new Promise((resolve, reject) => setTimeout(() => {
    try {
      resolve(window.taWasm.chirp_correction_global(
        new Float64Array(time), new Float64Array(wl), new Float64Array(flat), ta.length, time.length,
        opts.searchRange[0], opts.searchRange[1], opts.polyOrder,
        opts.snrThreshold, opts.nIter, opts.nSigma, opts.nBaseline
      ));
    } catch(e) { reject(e); }
  }, 0));
  return _wasmChirpResult(ret, ta.length, time.length) ||
    { TA2D: ta, coeffs: null, t0PerWl: new Array(ta.length).fill(NaN), snrPerWl: new Array(ta.length).fill(0), snrFilteredOut: new Array(ta.length).fill(true), sigmaClippedOut: new Array(ta.length).fill(false) };
}

// ==================== CPM Removal ====================
async function removeCpm(time, wl, ta, fitWindow, nBaseline) {
  var flat = _flattenTA(ta);
  var ret = window.taWasm.remove_cpm_wasm(
    new Float64Array(time), new Float64Array(wl), new Float64Array(flat), ta.length, time.length,
    fitWindow, nBaseline
  );
  if (ret && ret.ta_2d) return { TA2D: _unflattenTA(ret.ta_2d, ta.length, time.length) };
  return { TA2D: ta };
}

// ==================== IRF Deconvolution ====================
function estimateIrf(time, ta, nWlAvg) {
  const nWl = ta.length, nTime = time.length;
  if (nTime < 5 || nWl === 0) return null;

  // Average signal across central wavelengths
  const half = Math.floor(nWlAvg / 2);
  const wlStart = Math.max(0, Math.floor(nWl / 2) - half);
  const wlEnd = Math.min(nWl, Math.floor(nWl / 2) + half);
  const avgSignal = Array(nTime).fill(0);
  const counts = Array(nTime).fill(0);
  for (let i = wlStart; i < wlEnd; i++) {
    for (let j = 0; j < nTime; j++) {
      if (!isNaN(ta[i][j])) { avgSignal[j] += ta[i][j]; counts[j]++; }
    }
  }
  for (let j = 0; j < nTime; j++) avgSignal[j] = counts[j] > 0 ? avgSignal[j] / counts[j] : NaN;

  // Find t=0 index
  let t0Idx = 0, minAbs = Infinity;
  for (let i = 0; i < nTime; i++) { if (Math.abs(time[i]) < minAbs) { minAbs = Math.abs(time[i]); t0Idx = i; } }

  // Compute derivative near t=0
  const window = Math.min(nTime > 40 ? 20 : Math.floor(nTime / 2), t0Idx, nTime - t0Idx - 1);
  const start = Math.max(0, t0Idx - window), end = Math.min(nTime, t0Idx + window);
  if (end - start < 5) return null;

  const derivT = [], derivY = [];
  for (let j = start; j < end - 1; j++) {
    const dt = time[j + 1] - time[j];
    if (Math.abs(dt) < 1e-15 || isNaN(avgSignal[j + 1]) || isNaN(avgSignal[j])) continue;
    derivT.push((time[j] + time[j + 1]) / 2);
    derivY.push((avgSignal[j + 1] - avgSignal[j]) / dt);
  }
  if (derivT.length < 5) return null;

  // Find peak of absolute derivative
  let peakIdx = 0, peakVal = 0;
  for (let k = 0; k < derivY.length; k++) { if (Math.abs(derivY[k]) > peakVal) { peakVal = Math.abs(derivY[k]); peakIdx = k; } }
  if (peakVal < 1e-10) return null;

  // Estimate FWHM from half-max width of derivative
  let leftIdx = peakIdx, rightIdx = peakIdx;
  while (leftIdx > 0 && Math.abs(derivY[leftIdx]) > peakVal * 0.5) leftIdx--;
  while (rightIdx < derivY.length - 1 && Math.abs(derivY[rightIdx]) > peakVal * 0.5) rightIdx++;
  const fwhmDeriv = Math.abs(derivT[rightIdx] - derivT[leftIdx]);
  if (fwhmDeriv < 1e-15) return null;

  const sigma = fwhmDeriv / 2.355;
  const fwhm = 2.355 * sigma;
  return { fwhm, sigma };
}

async function deconvolveIrf(time, wl, ta, irfFwhm, nIter) {
  var flat = _flattenTA(ta);
  var ret = window.taWasm.deconvolve_irf_wasm(
    new Float64Array(time), new Float64Array(wl), new Float64Array(flat), ta.length, time.length,
    irfFwhm, nIter
  );
  if (ret && ret.ta_2d) return { TA2D: _unflattenTA(ret.ta_2d, ta.length, time.length), irfFwhm: ret.irf_fwhm };
  return { TA2D: ta, irfFwhm };
}

function makeHeatmapData(time, wl, ta, tRange) {
  const tMaskIdx = [];
  const tPlot = [];
  for (let j = 0; j < time.length; j++) {
    if (time[j] >= tRange[0] && time[j] <= tRange[1]) { tMaskIdx.push(j); tPlot.push(time[j]); }
  }
  // Build z matrix: z[row][col] = z[y][x]
  // x-axis = wavelength, y-axis = time
  // So z has time rows and wavelength columns
  const z = [];
  for (let j = 0; j < tMaskIdx.length; j++) {
    const row = [];
    for (let i = 0; i < wl.length; i++) {
      row.push(isNaN(ta[i][tMaskIdx[j]]) ? null : ta[i][tMaskIdx[j]] * 1000);
    }
    z.push(row);
  }
  return { tPlot, z };
}

let cancelProcessing = false;
let cancelGlobalKineticFit = false;

function ensureGlobalFitControls() {
  if ($('globalFitControls')) return;
  const html = `
    <div class="card" id="globalFitControls" style="display:none;">
      <h2>🔬 全部动力学拟合</h2>
      <p style="font-size:13px;color:#888;margin-bottom:12px;">可先设置统一拟合选项并应用到全部文件，再执行批量拟合或导出汇总结果。</p>
      <div class="fit-controls" style="margin-bottom:12px;">
        <div class="params" style="margin-bottom:12px;">
          <div class="param-group">
            <label>探测波长 (逗号分隔, nm)</label>
            <input type="text" id="globalFitWlInput" value="">
          </div>
          <div class="param-group">
            <label>指数个数</label>
            <select id="globalFitNExp">
              <option value="1">1指数</option>
              <option value="2" selected>2指数</option>
              <option value="3">3指数</option>
              <option value="4">4指数</option>
              <option value="5">5指数</option>
            </select>
          </div>
          <div class="param-group">
            <label>时间轴标度</label>
            <select id="globalFitTimeScale">
              <option value="linear">线性</option>
              <option value="log" selected>对数</option>
            </select>
          </div>
          <div class="param-group">
            <label>拟合速度/精度</label>
            <select id="globalFitQuality">
              <option value="0" selected>快速</option>
              <option value="1">标准</option>
              <option value="2">高精度</option>
            </select>
          </div>
        </div>
        <div style="display:flex;gap:12px;align-items:center;flex-wrap:wrap;">
          <button class="btn" id="globalFitApplyBtn" onclick="applyGlobalFitOptionsToAll()">应用到全部文件</button>
          <button class="btn btn-primary" id="globalFitBtn" onclick="doAllKineticFits()">全部开始拟合</button>
          <button class="btn" id="globalFitCancelBtn" style="display:none;background:#dc3545;color:white;" onclick="cancelGlobalKineticFit=true;">取消批量拟合</button>
          <button class="btn btn-download" id="globalFitSummaryBtn" onclick="downloadAllFitSummaryCSV()">导出全部拟合汇总 CSV</button>
          <button class="btn btn-download" id="globalZipProcessedBtn" onclick="downloadAllByTypeZip('processed_csv')">打包全部处理后 CSV</button>
          <button class="btn btn-download" id="globalZipFitBtn" onclick="downloadAllByTypeZip('fit_csv')">打包全部拟合参数 CSV</button>
          <span id="globalFitStatus" style="font-size:13px;color:#999;"></span>
        </div>
      </div>
    </div>`;
  $('resultsArea').insertAdjacentHTML('afterbegin', html);
}

function updateGlobalFitControls() {
  const panel = $('globalFitControls');
  if (!panel) return;
  const count = document.querySelectorAll('#resultsArea .card[data-fit-card="1"]').length;
  panel.style.display = count > 0 ? '' : 'none';
  const status = $('globalFitStatus');
  if (status && count > 0) status.textContent = `当前可批量拟合 ${count} 个文件`;
}

async function processAll() {
  if (!_wasmReady) { setStatus('❌ WASM 核心尚未加载，请刷新页面重试', 'error'); return; }
  if (uploadedFiles.length === 0) return;
  const wlMin = parseFloat($('wlMin').value);
  const wlMax = parseFloat($('wlMax').value);
  const nBaseline = parseInt($('nBaseline').value);
  const chirpMethod = $('chirpMethod').value;
  const polyOrder = parseInt($('polyOrder').value);
  const snrThreshold = parseFloat($('snrThreshold').value);
  const t0SearchMin = parseFloat($('t0SearchMin').value);
  const t0SearchMax = parseFloat($('t0SearchMax').value);
  const tViewMin = parseFloat($('tViewMin').value);
  const tViewMax = parseFloat($('tViewMax').value);
  const probeWlStr = $('probeWl').value;
  const probeWavelengths = probeWlStr.split(',').map(s => parseFloat(s.trim())).filter(v => !isNaN(v));

  const chirpOpts = {
    snrThreshold,
    polyOrder,
    searchRange: [t0SearchMin, t0SearchMax],
    nIter: 3,
    nSigma: 2.5,
    nBaseline
  };

  $('resultsArea').innerHTML = '';
  ensureGlobalFitControls();
  if ($('globalFitWlInput')) $('globalFitWlInput').value = probeWlStr;
  updateGlobalFitControls();
  cancelProcessing = false;
  $('processBtn').disabled = true;
  $('cancelBtn').style.display = 'inline-block';

  const yieldThread = () => new Promise(r => setTimeout(r, 0));
  const totalFiles = uploadedFiles.length;
  const stepsPerFile = 5;

  try {
  for (let fi = 0; fi < totalFiles; fi++) {
    if (cancelProcessing) break;
    const file = uploadedFiles[fi];
    const base = (fi / totalFiles) * 100;
    const step = (1 / totalFiles / stepsPerFile) * 100;

    try {
    setStatus(`⏳ [${fi + 1}/${totalFiles}] 读取 ${escapeHtml(file.name)}...`, 'info', base + step * 0);
    await yieldThread();
    if (cancelProcessing) break;

    let timeArray, wavelengthArray, TA2D;
    const isUfs = ufsFiles.includes(file.name);

    if (isUfs) {
      const arrayBuffer = await new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = e => resolve(e.target.result);
        reader.onerror = () => reject(new Error('文件读取失败，请确认以 HTTP 方式打开页面（非 file://）'));
        reader.readAsArrayBuffer(file);
      });

      if (!arrayBuffer || arrayBuffer.byteLength === 0) {
        throw new Error('文件内容为空或读取失败');
      }

      setStatus(`⏳ [${fi + 1}/${totalFiles}] 解析 UFS ${escapeHtml(file.name)}...`, 'info', base + step * 1);
      await yieldThread();
      if (cancelProcessing) break;

      const ufsData = parseUfsFile(arrayBuffer);
      timeArray = ufsData.times;
      wavelengthArray = ufsData.wavelengths;
      TA2D = ufsData.intensity;
    } else {
      const text = await new Promise(resolve => {
        const reader = new FileReader();
        reader.onload = e => resolve(e.target.result);
        reader.readAsText(file, 'utf-8');
      });

      setStatus(`⏳ [${fi + 1}/${totalFiles}] 解析 ${escapeHtml(file.name)}...`, 'info', base + step * 1);
      await yieldThread();
      if (cancelProcessing) break;

      const parsed = window.taWasm.parse_csv_wasm(text);
      if (!parsed || !parsed.time_array || parsed.time_array.length === 0) {
        throw new Error('CSV 解析失败或文件为空');
      }
      timeArray = Array.from(parsed.time_array);
      wavelengthArray = Array.from(parsed.wavelength_array);
      TA2D = Array.from(parsed.ta_2d).map(row => Array.from(row));
    }

    const { wavelengthArray: wl, TA2D: taCropped } = cropWavelength(wavelengthArray, TA2D, wlMin, wlMax);
    const taBeforeChirp = baselineSubtraction(timeArray, taCropped, nBaseline);

    setStatus(`⏳ [${fi + 1}/${totalFiles}] 啁啾校正 ${escapeHtml(file.name)}（可能需要数秒）...`, 'info', base + step * 2);
    await yieldThread();
    if (cancelProcessing) break;

    let chirpResult;
    if (chirpMethod === 'global') {
      chirpResult = await chirpCorrectionGlobal(timeArray, wl, taCropped, chirpOpts);
    } else if (chirpMethod === 'manual') {
      chirpResult = chirpCorrectionHalfHeight(timeArray, wl, taCropped, chirpOpts);
      chirpResult = { TA2D: taCropped, coeffs: null, t0PerWl: chirpResult.t0PerWl, snrPerWl: chirpResult.snrPerWl, snrFilteredOut: chirpResult.snrFilteredOut, sigmaClippedOut: chirpResult.sigmaClippedOut };
    } else {
      chirpResult = chirpCorrectionHalfHeight(timeArray, wl, taCropped, chirpOpts);
    }
    if (chirpResult && typeof chirpResult.then === 'function') chirpResult = await chirpResult;

    const { TA2D: taChirped, coeffs, t0PerWl, snrPerWl, snrFilteredOut, sigmaClippedOut, initialCoeffs } = chirpResult;
    const taAfter = baselineSubtraction(timeArray, taChirped, nBaseline);

    let taFinal = taAfter;
    const doCpm = $('doCpm').checked;
    const cpmFitWindow = parseFloat($('cpmFitWindow').value);
    if (doCpm) {
      setStatus(`⏳ [${fi + 1}/${totalFiles}] CPM 去除 ${escapeHtml(file.name)}...`, 'info', base + step * 2.5);
      await yieldThread();
      const cpmResult = await removeCpm(timeArray, wl, taAfter, cpmFitWindow, nBaseline);
      taFinal = cpmResult.TA2D;
    }

    const doIrf = $('doIrf').checked;
    let irfFwhmVal = parseFloat($('irfFwhm').value);
    const irfNIter = parseInt($('irfNIter').value);
    if (doIrf) {
      if (irfFwhmVal <= 0 || isNaN(irfFwhmVal)) {
        const est = estimateIrf(timeArray, taFinal, 10);
        irfFwhmVal = est ? est.fwhm : 0;
        if (!est) console.warn('[IRF] 自动估算失败，跳过解卷积');
      }
      if (irfFwhmVal > 0) {
        setStatus(`⏳ [${fi + 1}/${totalFiles}] IRF 解卷积 ${escapeHtml(file.name)} (FWHM=${irfFwhmVal.toFixed(3)} ps)...`, 'info', base + step * 2.7);
        await yieldThread();
        const irfResult = await deconvolveIrf(timeArray, wl, taFinal, irfFwhmVal, irfNIter);
        taFinal = irfResult.TA2D;
      }
    }

    if (cancelProcessing) break;
    setStatus(`⏳ [${fi + 1}/${totalFiles}] 渲染结果 ${escapeHtml(file.name)}...`, 'info', base + step * 3);
    await yieldThread();

    renderResults(file.name, timeArray, wl, taBeforeChirp, taFinal, coeffs, t0PerWl, chirpMethod, tViewMin, tViewMax, probeWavelengths, snrPerWl, snrFilteredOut, sigmaClippedOut, initialCoeffs);

    } catch(e) {
      if (e instanceof WebAssembly.RuntimeError) {
        _wasmReady = false;
        setStatus(`❌ WASM 运行时崩溃，请刷新页面重试（${e.message}）`, 'error');
        break;
      }
      setStatus(`❌ [${fi + 1}/${totalFiles}] ${escapeHtml(file.name)} 处理失败：${escapeHtml(e.message)}`, 'error');
    }
  }
  } finally {
    $('processBtn').disabled = false;
    $('cancelBtn').style.display = 'none';
    if (cancelProcessing) {
      setStatus('⏹ 处理已取消', 'error');
      cancelProcessing = false;
    } else if (_wasmReady) {
      setStatus('✅ 全部处理完成!', 'success', 100);
    }
  }
}

function renderResults(fileName, time, wl, taBefore, taAfter, coeffs, t0PerWl, chirpMethod, tViewMin, tViewMax, probeWavelengths, snrPerWl, snrFilteredOut, sigmaClippedOut, initialCoeffs) {
  const baseName = fileName.replace(/[^a-zA-Z0-9]/g, '_');
  const divId = `result_${baseName}`;
  const methodName = chirpMethod === 'global' ? '全局优化法' : chirpMethod === 'manual' ? '手动选点法' : '半高点法';

  let maxSigTime = 0, maxSigVal = 0;
  for (let j = 0; j < time.length; j++) {
    if (time[j] < 0) continue;
    let colSum = 0;
    for (let i = 0; i < taAfter.length; i++) {
      if (!isNaN(taAfter[i][j])) colSum += Math.abs(taAfter[i][j]);
    }
    if (colSum > maxSigVal) { maxSigVal = colSum; maxSigTime = time[j]; }
  }
  const maxTime = time[time.length - 1];

  const timeDiffs = [];
  for (let j = 1; j < time.length; j++) {
    const dt = time[j] - time[j - 1];
    if (dt > 0) timeDiffs.push(dt);
  }
  timeDiffs.sort((a, b) => a - b);
  const medianDt = timeDiffs.length > 0 ? timeDiffs[Math.floor(timeDiffs.length / 2)] : 0.5;

  let html = `<div class="card" id="${divId}" data-fit-card="1" data-base-name="${baseName}" data-div-id="${divId}">
    <h2>📄 ${escapeHtml(fileName)}</h2>
    <p style="font-size:13px;color:#888;">波长: ${wl[0].toFixed(1)} ~ ${wl[wl.length-1].toFixed(1)} nm | 时间: ${time[0].toFixed(3)} ~ ${time[time.length-1].toFixed(3)} ps | 数据: ${wl.length}×${time.length}</p>
    <div class="tabs">
      <div class="tab active" onclick="switchTab(this, '${divId}', 'tabViz')">数据可视化</div>
      <div class="tab" onclick="switchTab(this, '${divId}', 'tabSlice')">光谱切片</div>
      <div class="tab" onclick="switchTab(this, '${divId}', 'tabFit')">动力学拟合</div>
      <div class="tab" onclick="switchTab(this, '${divId}', 'tabDownload')">下载数据</div>
    </div>
    <div class="tab-content active" id="${divId}_tabViz">
      <h3 style="font-size:14px;color:#ff6666;margin-bottom:8px;">啁啾校正</h3>
      ${chirpMethod === 'manual' ? `
      <div id="${divId}_manualPanel" style="background:#1a2a1a;border:1px solid #3a5a3a;border-radius:8px;padding:12px 16px;margin-bottom:12px;">
        <p style="font-size:13px;color:#8f8;margin-bottom:8px;">🖱️ 点击啁啾曲线图添加控制点（至少选 ${parseInt($('manualFitOrder')?.value || 3) + 1} 个），然后点击「应用拟合」</p>
        <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;">
          <button class="btn btn-primary" id="${divId}_applyManual" onclick="applyManualChirp('${baseName}','${divId}')" disabled>✅ 应用拟合</button>
          <button class="btn" style="background:#555;color:#ddd;" onclick="undoLastManualPoint('${divId}')">↩ 撤销上一点</button>
          <button class="btn" style="background:#555;color:#ddd;" onclick="clearManualPoints('${divId}')">🗑 清除所有点</button>
          <span id="${divId}_pointCount" style="font-size:12px;color:#aaa;margin-left:8px;">已选 0 个点</span>
        </div>
      </div>` : ''}
      <div class="plot-grid">
        <div class="plot-box"><h3>啁啾曲线${chirpMethod === 'manual' ? '（点击添加控制点）' : ''}</h3><div id="${divId}_chirpCurve" style="width:100%;height:350px;"></div></div>
        <div class="plot-box"><h3>校正前后对比</h3><div id="${divId}_chirpCompare" style="width:100%;height:350px;"></div></div>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">动力学曲线</h3>
      <div class="plot-grid">
        <div class="plot-box"><h3>校正前</h3><div id="${divId}_beforeKin" style="width:100%;height:350px;"></div></div>
        <div class="plot-box"><h3>校正后</h3><div id="${divId}_afterKin" style="width:100%;height:350px;"></div></div>
      </div>
    </div>
    <div class="tab-content" id="${divId}_tabSlice">
      <div class="fit-controls">
        <div class="params" style="margin-bottom:16px;">
          <div class="param-group" style="grid-column:1/-1;">
            <label>时间延迟 (逗号分隔, ps)</label>
            <input type="text" id="${divId}_sliceTimes" value="0.1, 1, 10, 100, 1000, ${maxTime.toFixed(1)}">
          </div>
        </div>
        <div style="text-align:center;margin-bottom:16px;">
          <button class="btn btn-primary" onclick="doSpectralSlice('${baseName}', '${divId}')">📊 绘制光谱切片</button>
        </div>
      </div>
      <div class="plot-box" style="min-height:400px;">
        <div id="${divId}_slicePlot" style="width:100%;height:400px;"></div>
      </div>
    </div>
    <div class="tab-content" id="${divId}_tabFit">
      <div class="fit-controls">
        <div class="params" style="margin-bottom:16px;">
          <div class="param-group">
            <label>探测波长 (逗号分隔, nm)</label>
            <input type="text" id="${divId}_fitWlInput" value="${probeWavelengths.join(', ')}">
          </div>
          <div class="param-group">
            <label>指数个数</label>
            <select id="${divId}_fitNExp">
              <option value="1">1指数</option>
              <option value="2" selected>2指数</option>
              <option value="3">3指数</option>
              <option value="4">4指数</option>
              <option value="5">5指数</option>
            </select>
          </div>
          <div class="param-group">
            <label>拟合时间下限 (ps)</label>
            <input type="number" id="${divId}_fitTMin" value="${maxSigTime.toFixed(2)}" step="${medianDt.toFixed(6)}" min="${time[0].toFixed(6)}" max="${time[time.length-1].toFixed(6)}">
          </div>
          <div class="param-group">
            <label>拟合时间上限 (ps)</label>
            <input type="number" id="${divId}_fitTMax" value="${maxTime.toFixed(2)}" step="${medianDt.toFixed(6)}" min="${time[0].toFixed(6)}" max="${time[time.length-1].toFixed(6)}">
          </div>
          <div class="param-group">
            <label>时间轴标度</label>
            <select id="${divId}_fitTimeScale" onchange="updateFitTimeScale('${divId}')">
              <option value="linear">线性</option>
              <option value="log" selected>对数</option>
            </select>
          </div>
          <div class="param-group">
            <label>拟合速度/精度</label>
            <select id="${divId}_fitQuality">
              <option value="0" selected>快速</option>
              <option value="1">标准</option>
              <option value="2">高精度</option>
            </select>
          </div>
        </div>
        <div style="text-align:center;margin-bottom:16px;">
          <button class="btn btn-primary" onclick="doKineticFit('${baseName}', '${divId}')">🔬 开始拟合</button>
        </div>
      </div>
      <div class="plot-box" style="min-height:400px;">
        <div id="${divId}_fitPlot" style="width:100%;height:400px;"></div>
      </div>
      <div id="${divId}_fitResult"></div>
    </div>
    <div class="tab-content" id="${divId}_tabDownload">
      <h3 style="font-size:14px;color:#ff6666;margin-bottom:8px;">预处理数据下载</h3>
      <div class="download-section">
        <button class="btn btn-download" onclick="downloadJSON('${baseName}')">⬇️ JSON</button>
        <button class="btn btn-download" onclick="downloadCSV('${baseName}')">⬇️ CSV</button>
        <button class="btn btn-download" onclick="downloadASCII('${baseName}')">⬇️ ASCII (空格分隔)</button>
        <button class="btn btn-download" onclick="downloadTSV('${baseName}')">⬇️ TSV (Tab分隔)</button>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">拟合结果下载</h3>
      <div class="download-section">
        <button class="btn btn-download" onclick="downloadFitCSV('${baseName}')">⬇️ 拟合参数 (Excel CSV)</button>
        <button class="btn btn-download" onclick="downloadFitRawCSV('${baseName}')">⬇️ 拟合原始+拟合数据 (Excel CSV)</button>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">光谱切片下载</h3>
      <div class="download-section">
        <button class="btn btn-download" onclick="downloadSliceCSV('${baseName}')">⬇️ 光谱切片 (Excel CSV)</button>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">图片导出 (PNG)</h3>
      <div class="download-section">
        <button class="btn btn-download" onclick="downloadAllPlotsPNG('${baseName}','${divId}')">⬇️ 全部图片</button>
        <button class="btn btn-download" onclick="downloadPlotPNG('${divId}','chirpCurve','${baseName}')">🖼️ 啁啾曲线</button>
        <button class="btn btn-download" onclick="downloadPlotPNG('${divId}','beforeKin','${baseName}')">🖼️ 校正前动力学</button>
        <button class="btn btn-download" onclick="downloadPlotPNG('${divId}','afterKin','${baseName}')">🖼️ 校正后动力学</button>
        <button class="btn btn-download" onclick="downloadPlotPNG('${divId}','fitPlot','${baseName}')">🖼️ 动力学拟合</button>
        <button class="btn btn-download" onclick="downloadPlotPNG('${divId}','slicePlot','${baseName}')">🖼️ 光谱切片</button>
      </div>
      <div style="margin-top:12px;font-size:13px;color:#888;" id="${divId}_info"></div>
    </div>
  </div>`;

  $('resultsArea').insertAdjacentHTML('beforeend', html);
  updateGlobalFitControls();

  const tRange = [tViewMin, tViewMax];

  const beforeData = makeHeatmapData(time, wl, taBefore, tRange);
  const afterData = makeHeatmapData(time, wl, taAfter, tRange);

  // Compute shared z range so both heatmaps use the same color scale
  let zMin = Infinity, zMax = -Infinity;
  for (const row of beforeData.z) {
    for (const v of row) {
      if (v !== null && v < zMin) zMin = v;
      if (v !== null && v > zMax) zMax = v;
    }
  }
  for (const row of afterData.z) {
    for (const v of row) {
      if (v !== null && v < zMin) zMin = v;
      if (v !== null && v > zMax) zMax = v;
    }
  }
  const zRange = Math.max(Math.abs(zMin), Math.abs(zMax));

  const tMaskIdx = [];
  const tPlot = [];
  for (let j = 0; j < time.length; j++) {
    if (time[j] >= tViewMin && time[j] <= tViewMax) { tMaskIdx.push(j); tPlot.push(time[j]); }
  }

  const colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'];
  const kinTracesBefore = [], kinTracesAfter = [];
  probeWavelengths.forEach((pw, pi) => {
    let idxWl = 0;
    let minDiff = Infinity;
    for (let i = 0; i < wl.length; i++) {
      if (Math.abs(wl[i] - pw) < minDiff) { minDiff = Math.abs(wl[i] - pw); idxWl = i; }
    }
    const color = colors[pi % colors.length];
    kinTracesBefore.push({
      x: tPlot, y: tMaskIdx.map(j => taBefore[idxWl][j] * 1000),
      name: `${pw}nm`, line: { color, dash: 'dash' }, opacity: 0.5
    });
    kinTracesAfter.push({
      x: tPlot, y: tMaskIdx.map(j => taAfter[idxWl][j] * 1000),
      name: `${pw}nm`, line: { color }
    });
  });

  const kinLayout = {
    xaxis: { title: '时间 (ps)' }, yaxis: { title: 'ΔA (mOD)' },
    margin: { l: 60, r: 20, t: 30, b: 50 },
    shapes: [{ type: 'line', x0: 0, x1: 0, y0: 0, y1: 1, yref: 'paper', line: { color: 'gray', dash: 'dot' } }],
    legend: { font: { size: 10 } }
  };

  Plotly.newPlot($(`${divId}_beforeKin`), kinTracesBefore, { ...kinLayout, title: '校正前' }, { responsive: true });
  Plotly.newPlot($(`${divId}_afterKin`), kinTracesAfter, { ...kinLayout, title: '校正后' }, { responsive: true });

  const chirpTraces = [];

  if (snrFilteredOut) {
    const snrFiltWl = [], snrFiltT0 = [];
    for (let i = 0; i < wl.length; i++) {
      if (snrFilteredOut[i] && !isNaN(t0PerWl[i])) {
        snrFiltWl.push(wl[i]);
        snrFiltT0.push(t0PerWl[i] * 1000);
      }
    }
    if (snrFiltWl.length > 0) {
      chirpTraces.push({
        x: snrFiltWl, y: snrFiltT0,
        mode: 'markers', name: 'SNR过滤',
        marker: { size: 4, color: 'gray', opacity: 0.4, symbol: 'circle' }
      });
    }
  }

  if (sigmaClippedOut) {
    const sigmaWl = [], sigmaT0 = [];
    for (let i = 0; i < wl.length; i++) {
      if (sigmaClippedOut[i] && !isNaN(t0PerWl[i])) {
        sigmaWl.push(wl[i]);
        sigmaT0.push(t0PerWl[i] * 1000);
      }
    }
    if (sigmaWl.length > 0) {
      chirpTraces.push({
        x: sigmaWl, y: sigmaT0,
        mode: 'markers', name: 'Sigma剔除',
        marker: { size: 6, color: 'orange', opacity: 0.7, symbol: 'x' }
      });
    }
  }

  const keptWl = [], keptT0 = [];
  for (let i = 0; i < wl.length; i++) {
    const isFiltered = snrFilteredOut && snrFilteredOut[i];
    const isClipped = sigmaClippedOut && sigmaClippedOut[i];
    if (!isFiltered && !isClipped && !isNaN(t0PerWl[i])) {
      keptWl.push(wl[i]);
      keptT0.push(t0PerWl[i] * 1000);
    }
  }
  chirpTraces.push({
    x: keptWl, y: keptT0,
    mode: 'markers', name: '有效t0',
    marker: { size: 4, opacity: 0.5, color: 'blue' }
  });

  if (chirpMethod === 'manual') {
    chirpTraces.push({ x: [], y: [], mode: 'markers+text', name: '手动选点', marker: { size: 10, color: '#00cc00', symbol: 'circle' }, text: [], textposition: 'top center', textfont: { size: 9, color: '#00cc00' } });
  } else if (coeffs) {
    const wlFine = linspace(wl[0], wl[wl.length - 1], 200);
    const invL2Fine = wlFine.map(w => 1.0 / (w * w));
    const t0Fit = polyval(coeffs, invL2Fine).map(v => v * 1000);
    chirpTraces.push({ x: wlFine, y: t0Fit, mode: 'lines', name: chirpMethod === 'global' ? '全局优化拟合' : '多项式拟合', line: { color: 'red', width: 2 } });
    if (chirpMethod === 'global' && initialCoeffs) {
      const t0Init = polyval(initialCoeffs, invL2Fine).map(v => v * 1000);
      chirpTraces.push({ x: wlFine, y: t0Init, mode: 'lines', name: '初始估计拟合', line: { color: 'green', dash: 'dash', width: 1.5 } });
    }
  }

  // Store initial axis ranges for manual mode to prevent auto-rescaling
  if (chirpMethod === 'manual') {
    const validT0Fs = t0PerWl.filter(v => !isNaN(v)).map(v => v * 1000);
    const yMin = validT0Fs.length > 0 ? Math.min(...validT0Fs) * 0.8 : -1000;
    const yMax = validT0Fs.length > 0 ? Math.max(...validT0Fs) * 1.2 : 1000;
    window._manualChirpAxisRange = {
      x: [wl[0], wl[wl.length - 1]],
      y: [yMin, yMax]
    };
  }

  Plotly.newPlot($(`${divId}_chirpCurve`), chirpTraces, {
    xaxis: { title: '波长 (nm)', range: chirpMethod === 'manual' ? window._manualChirpAxisRange?.x : undefined },
    yaxis: { title: 't0 (fs)', range: chirpMethod === 'manual' ? window._manualChirpAxisRange?.y : undefined },
    title: '啁啾曲线', margin: { l: 60, r: 20, t: 40, b: 50 },
    clickmode: 'event+select'
  }, { responsive: true });

  // Manual chirp: set up click-to-add-point interaction
  if (chirpMethod === 'manual') {
    if (!window._manualChirpData) window._manualChirpData = {};
    window._manualChirpData[divId] = { points: [], baseName, wl, time, taBefore };

    const chirpEl = $(`${divId}_chirpCurve`);
    // Add a transparent scatter trace covering the whole plot area to capture clicks anywhere
    Plotly.addTraces(chirpEl, {
      x: [wl[0], wl[wl.length - 1]],
      y: keptT0.length > 0 ? [Math.min(...keptT0) * 0.5, Math.max(...keptT0) * 1.5] : [-1000, 1000],
      mode: 'markers',
      marker: { size: 1, opacity: 0 },
      name: '_click_capture_',
      hoverinfo: 'skip',
      showlegend: false
    });

    chirpEl.on('plotly_click', function(eventData) {
      var xVal, yVal;
      // If clicked on an existing trace point, use that
      if (eventData.points && eventData.points.length > 0) {
        const pt = eventData.points[0];
        xVal = pt.x;
        yVal = pt.y;
      }
      if (xVal === undefined || yVal === undefined) return;
      const md = window._manualChirpData[divId];
      md.points.push({ wl: xVal, t0: yVal });
      updateManualChirpPlot(divId);
    });
  }

  Plotly.newPlot($(`${divId}_chirpCompare`), [
    { z: beforeData.z, x: wl, y: beforeData.tPlot, type: 'heatmap', colorscale: 'RdBu', reversescale: true, zmid: 0, zmin: -zRange, zmax: zRange, name: '校正前', showscale: false },
    { z: afterData.z, x: wl, y: afterData.tPlot, type: 'heatmap', colorscale: 'RdBu', reversescale: true, zmid: 0, zmin: -zRange, zmax: zRange, name: '校正后', xaxis: 'x2', yaxis: 'y2', showscale: false }
  ], {
    title: '校正前后对比',
    grid: { columns: 2, rows: 1, pattern: 'independent', xgap: 0.08 },
    margin: { l: 60, r: 20, t: 40, b: 50 },
    xaxis: { title: '波长 (nm)', domain: [0, 0.46] },
    yaxis: { title: '时间 (ps)' },
    xaxis2: { title: '波长 (nm)', anchor: 'y2' },
    yaxis2: { title: '时间 (ps)', anchor: 'x2' },
    annotations: [
      { text: '校正前', x: 0.22, y: 1.05, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 12 } },
      { text: '校正后', x: 0.78, y: 1.05, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 12 } }
    ]
  }, { responsive: true });

  if (coeffs) {
    const validT0Fs = t0PerWl.filter(v => !isNaN(v)).map(v => v * 1000);
    const refT0 = validT0Fs.reduce((a, b) => a + b, 0) / validT0Fs.length;
    const nSnrFilt = snrFilteredOut ? snrFilteredOut.filter(v => v).length : 0;
    const nSigmaClip = sigmaClippedOut ? sigmaClippedOut.filter(v => v).length : 0;
    const nKept = keptWl.length;
    $(`${divId}_info`).innerHTML = `参考t0 = ${refT0.toFixed(1)} fs | 啁啾范围: ${Math.min(...validT0Fs).toFixed(1)} ~ ${Math.max(...validT0Fs).toFixed(1)} fs | 方法: ${methodName} | SNR过滤: ${nSnrFilt}个 | Sigma剔除: ${nSigmaClip}个 | 拟合使用: ${nKept}个`;
  }

  window[`data_${baseName}`] = { timeArray: time, wavelengthArray: wl, TA2D: taAfter, taBefore: taBefore, coeffs, t0PerWl, snrPerWl, snrFilteredOut, sigmaClippedOut };
}

function doSpectralSlice(baseName, divId) {
  const data = window[`data_${baseName}`];
  if (!data) return;

  const timesStr = $(`${divId}_sliceTimes`).value;
  const delayTimes = timesStr.split(',').map(s => parseFloat(s.trim())).filter(v => !isNaN(v));

  if (delayTimes.length === 0) {
    $(`${divId}_slicePlot`).innerHTML = '<div class="status status-error">请输入有效的时间延迟</div>';
    return;
  }

  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;

  const colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'];
  const traces = [];
  const sliceData = [];

  delayTimes.forEach((dt, pi) => {
    let idxT = 0;
    let minDiff = Infinity;
    for (let j = 0; j < time.length; j++) {
      if (Math.abs(time[j] - dt) < minDiff) { minDiff = Math.abs(time[j] - dt); idxT = j; }
    }
    const actualTime = time[idxT];
    const spectrum = ta.map(row => row[idxT] * 1000);
    const color = colors[pi % colors.length];

    traces.push({
      x: wl, y: spectrum,
      mode: 'lines', name: `${actualTime.toFixed(2)} ps`,
      line: { color, width: 2 }
    });

    sliceData.push({ time: actualTime, spectrum: spectrum.map(v => isNaN(v) ? '' : v) });
  });

  window[`sliceData_${baseName}`] = { wl, slices: sliceData };

  Plotly.newPlot($(`${divId}_slicePlot`), traces, {
    xaxis: { title: '波长 (nm)' },
    yaxis: { title: 'ΔA (mOD)' },
    title: '光谱切片',
    margin: { l: 60, r: 20, t: 40, b: 50 },
    legend: { font: { size: 10 } }
  }, { responsive: true });
}

function downloadSliceCSV(baseName) {
  const data = window[`sliceData_${baseName}`];
  if (!data) return;

  const headers = ['波长(nm)'];
  for (const s of data.slices) headers.push(`${s.time.toFixed(2)}ps(mOD)`);

  let csv = headers.join(',') + '\n';
  for (let i = 0; i < data.wl.length; i++) {
    const row = [data.wl[i].toFixed(2)];
    for (const s of data.slices) row.push(s.spectrum[i]);
    csv += row.join(',') + '\n';
  }

  const bom = '\uFEFF';
  triggerDownload(bom + csv, `${baseName}_spectral_slices.csv`, 'text/csv;charset=utf-8');
}

function switchTab(el, divId, tabName) {
  const parent = el.parentElement;
  parent.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
  el.classList.add('active');
  const container = $(divId);
  container.querySelectorAll('.tab-content').forEach(tc => tc.classList.remove('active'));
  $(`${divId}_${tabName}`).classList.add('active');
  setTimeout(() => window.dispatchEvent(new Event('resize')), 100);
}

function triggerBlobDownload(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function triggerDownload(content, filename, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  triggerBlobDownload(blob, filename);
}

function buildProcessedJSON(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return null;
  return {
    filename: `${baseName}_processed.json`,
    mimeType: 'application/json',
    content: JSON.stringify({
      time_array: data.timeArray,
      wavelength_array: data.wavelengthArray,
      TA_2D_data: data.TA2D,
      chirp_coeffs: data.coeffs
    }, null, 2)
  };
}

function buildProcessedCSV(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return null;
  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;
  let csv = ',' + time.map(t => t.toFixed(6)).join(',') + '\n';
  for (let i = 0; i < wl.length; i++) {
    csv += wl[i].toFixed(2) + ',' + ta[i].map(v => isNaN(v) ? '' : v.toExponential(8)).join(',') + '\n';
  }
  return { filename: `${baseName}_processed.csv`, mimeType: 'text/csv', content: csv };
}

function buildProcessedASCII(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return null;
  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;
  let txt = 'wavelength(nm)';
  for (let j = 0; j < time.length; j++) txt += `  ${time[j].toFixed(6)}`;
  txt += '\n';
  for (let i = 0; i < wl.length; i++) {
    txt += wl[i].toFixed(2);
    for (let j = 0; j < time.length; j++) {
      txt += `  ${isNaN(ta[i][j]) ? 'NaN' : ta[i][j].toExponential(8)}`;
    }
    txt += '\n';
  }
  return { filename: `${baseName}_processed.dat`, mimeType: 'text/plain', content: txt };
}

function buildProcessedTSV(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return null;
  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;
  let tsv = '\t' + time.map(t => t.toFixed(6)).join('\t') + '\n';
  for (let i = 0; i < wl.length; i++) {
    tsv += wl[i].toFixed(2) + '\t' + ta[i].map(v => isNaN(v) ? '' : v.toExponential(8)).join('\t') + '\n';
  }
  return { filename: `${baseName}_processed.tsv`, mimeType: 'text/tab-separated-values', content: tsv };
}

function buildFitCSV(baseName) {
  const results = window[`fitResults_${baseName}`];
  if (!results || results.length === 0) return null;

  const maxNExp = Math.max(...results.map(r => r.nExp));
  const headers = ['波长(nm)', '指数个数', 'R²'];
  for (let k = 0; k < maxNExp; k++) {
    headers.push(`τ${k + 1}(ps)`, `τ${k + 1}标准差`, `A${k + 1}(mOD)`, `A${k + 1}标准差(mOD)`);
  }
  headers.push('y0(mOD)', 'y0标准差(mOD)', '异常标记');

  let csv = headers.join(',') + '\n';
  for (const r of results) {
    const row = [r.wavelength.toFixed(2), r.nExp, r.r2.toFixed(6)];
    for (let k = 0; k < maxNExp; k++) {
      if (k < r.nExp) {
        const tau = r.params[k * 2 + 1];
        const A = r.params[k * 2] * 1000;
        const hasStd = r.stdErrs && r.stdErrs.length > 0;
        const tauStd = hasStd ? r.stdErrs[k * 2 + 1] : '';
        const aStd = hasStd ? r.stdErrs[k * 2] * 1000 : '';
        row.push(tau.toFixed(6), hasStd ? tauStd.toFixed(6) : '', A.toFixed(4), hasStd ? aStd.toFixed(4) : '');
      } else {
        row.push('', '', '', '');
      }
    }
    const offset = r.params[r.params.length - 1] * 1000;
    const hasStd = r.stdErrs && r.stdErrs.length > 0;
    const offStd = hasStd ? r.stdErrs[r.stdErrs.length - 1] * 1000 : '';
    row.push(offset.toFixed(4), hasStd ? offStd.toFixed(4) : '', (r.flags || []).join('; '));
    csv += row.join(',') + '\n';
  }

  return { filename: `${baseName}_fit_results.csv`, mimeType: 'text/csv;charset=utf-8', content: '\uFEFF' + csv };
}

function buildFitRawCSV(baseName) {
  const results = window[`fitResults_${baseName}`];
  if (!results || results.length === 0) return null;

  let csv = '\uFEFF';
  for (const r of results) {
    if (!r.tData || r.tData.length === 0) continue;
    const wl = r.wavelength.toFixed(2);
    const nData = r.tData.length;
    const nFit = r.tFit ? r.tFit.length : 0;
    const maxLen = Math.max(nData, nFit);

    csv += `# ${wl} nm  R²=${r.r2.toFixed(6)}  nExp=${r.nExp}\n`;
    csv += `时间(ps),原始ΔA(mOD),拟合ΔA(mOD)\n`;
    for (let k = 0; k < maxLen; k++) {
      const t = k < nData ? r.tData[k] : '';
      const yd = k < nData ? (isNaN(r.yData[k]) ? '' : r.yData[k].toExponential(6)) : '';
      const yf = k < nFit ? (isNaN(r.yFit[k]) ? '' : r.yFit[k].toExponential(6)) : '';
      csv += `${t},${yd},${yf}\n`;
    }
    csv += '\n';
  }

  return { filename: `${baseName}_fit_raw_data.csv`, mimeType: 'text/csv;charset=utf-8', content: csv };
}

function buildSliceCSV(baseName) {
  const data = window[`sliceData_${baseName}`];
  if (!data) return null;

  const headers = ['波长(nm)'];
  for (const s of data.slices) headers.push(`${s.time.toFixed(2)}ps(mOD)`);

  let csv = headers.join(',') + '\n';
  for (let i = 0; i < data.wl.length; i++) {
    const row = [data.wl[i].toFixed(2)];
    for (const s of data.slices) row.push(s.spectrum[i]);
    csv += row.join(',') + '\n';
  }

  return { filename: `${baseName}_spectral_slices.csv`, mimeType: 'text/csv;charset=utf-8', content: '\uFEFF' + csv };
}

function downloadJSON(baseName) {
  const file = buildProcessedJSON(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
}

function downloadCSV(baseName) {
  const file = buildProcessedCSV(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
}

function downloadASCII(baseName) {
  const file = buildProcessedASCII(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
}

function downloadTSV(baseName) {
  const file = buildProcessedTSV(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
}

function downloadFitCSV(baseName) {
  const results = window[`fitResults_${baseName}`];
  if (!results || results.length === 0) return;

  const maxNExp = Math.max(...results.map(r => r.nExp));
  const headers = ['波长(nm)', '指数个数', 'R²'];
  for (let k = 0; k < maxNExp; k++) {
    headers.push(`τ${k + 1}(ps)`, `τ${k + 1}标准差`, `A${k + 1}(mOD)`, `A${k + 1}标准差(mOD)`);
  }
  headers.push('y0(mOD)', 'y0标准差(mOD)');

  let csv = headers.join(',') + '\n';
  for (const r of results) {
    const row = [r.wavelength.toFixed(2), r.nExp, r.r2.toFixed(6)];
    for (let k = 0; k < maxNExp; k++) {
      if (k < r.nExp) {
        const tau = r.params[k * 2 + 1];
        const A = r.params[k * 2] * 1000;
        const hasStd = r.stdErrs && r.stdErrs.length > 0;
        const tauStd = hasStd ? r.stdErrs[k * 2 + 1] : '';
        const aStd = hasStd ? r.stdErrs[k * 2] * 1000 : '';
        row.push(tau.toFixed(6), hasStd ? tauStd.toFixed(6) : '', A.toFixed(4), hasStd ? aStd.toFixed(4) : '');
      } else {
        row.push('', '', '', '');
      }
    }
    const offset = r.params[r.params.length - 1] * 1000;
    const hasStd = r.stdErrs && r.stdErrs.length > 0;
    const offStd = hasStd ? r.stdErrs[r.stdErrs.length - 1] * 1000 : '';
    row.push(offset.toFixed(4), hasStd ? offStd.toFixed(4) : '');
    csv += row.join(',') + '\n';
  }

  const bom = '\uFEFF';
  triggerDownload(bom + csv, `${baseName}_fit_results.csv`, 'text/csv;charset=utf-8');
}

function downloadFitRawCSV(baseName) {
  const results = window[`fitResults_${baseName}`];
  if (!results || results.length === 0) return;

  const bom = '\uFEFF';
  let csv = bom;

  for (const r of results) {
    if (!r.tData || r.tData.length === 0) continue;
    const wl = r.wavelength.toFixed(2);
    const nData = r.tData.length;
    const nFit = r.tFit ? r.tFit.length : 0;
    const maxLen = Math.max(nData, nFit);

    csv += `# ${wl} nm  R²=${r.r2.toFixed(6)}  nExp=${r.nExp}\n`;
    csv += `时间(ps),原始ΔA(mOD),拟合ΔA(mOD)\n`;
    for (let k = 0; k < maxLen; k++) {
      const t = k < nData ? r.tData[k] : '';
      const yd = k < nData ? (isNaN(r.yData[k]) ? '' : r.yData[k].toExponential(6)) : '';
      const yf = k < nFit ? (isNaN(r.yFit[k]) ? '' : r.yFit[k].toExponential(6)) : '';
      csv += `${t},${yd},${yf}\n`;
    }
    csv += '\n';
  }

  triggerDownload(csv, `${baseName}_fit_raw_data.csv`, 'text/csv;charset=utf-8');
}

function downloadPlotPNG(divId, plotSuffix, baseName) {
  const el = document.getElementById(`${divId}_${plotSuffix}`);
  if (!el) return;
  Plotly.downloadImage(el, {
    format: 'png',
    width: 1200,
    height: 600,
    filename: `${baseName}_${plotSuffix}`,
    scale: 2
  });
}

function downloadAllPlotsPNG(baseName, divId) {
  const plots = ['chirpCurve', 'beforeKin', 'afterKin', 'fitPlot', 'slicePlot'];
  plots.forEach((p, i) => {
    setTimeout(() => downloadPlotPNG(divId, p, baseName), i * 300);
  });
}

// ============================================================
// Manual Chirp Correction - Interactive Point Selection
// ============================================================

function updateManualChirpPlot(divId) {
  const md = window._manualChirpData[divId];
  if (!md) return;
  const pts = md.points;
  const fitOrder = parseInt($('manualFitOrder')?.value || 3);
  const minPts = fitOrder + 1;

  // Update point count display
  const countEl = $(`${divId}_pointCount`);
  if (countEl) countEl.textContent = `已选 ${pts.length} 个点（需 ≥ ${minPts}）`;

  // Enable/disable apply button
  const applyBtn = $(`${divId}_applyManual`);
  if (applyBtn) applyBtn.disabled = pts.length < minPts;

  // Update the manual points trace on Plotly
  const chirpEl = $(`${divId}_chirpCurve`);
  if (!chirpEl) return;

  const xPts = pts.map(p => p.wl);
  const yPts = pts.map(p => p.t0);
  const labels = pts.map((p, i) => `#${i+1}: ${p.wl.toFixed(1)}nm, ${p.t0.toFixed(0)}fs`);

  // Update or add manual points trace (keep axis range fixed)
  const existingTraces = chirpEl.data;
  const manualTraceIdx = existingTraces.findIndex(t => t.name === '手动选点');
  const axisRange = window._manualChirpAxisRange;
  const fixedLayout = axisRange ? {
    'xaxis.range': axisRange.x,
    'yaxis.range': axisRange.y
  } : {};
  if (manualTraceIdx >= 0) {
    Plotly.update(chirpEl, {
      x: [xPts],
      y: [yPts],
      text: [labels]
    }, fixedLayout, [manualTraceIdx]);
  } else {
    Plotly.addTraces(chirpEl, {
      x: xPts, y: yPts, mode: 'markers+text',
      name: '手动选点',
      marker: { size: 10, color: '#ff6600', symbol: 'diamond' },
      text: labels,
      textposition: 'top center'
    });
    if (axisRange) {
      Plotly.relayout(chirpEl, fixedLayout);
    }
  }

  // If enough points, also show the fit curve as trace 3
  // We add trace 3 dynamically or update it
  if (pts.length >= minPts) {
    const xArr = pts.map(p => p.wl);
    const yArr = pts.map(p => p.t0 / 1000); // convert back to ps for polyfit
    const actualOrder = Math.min(fitOrder, pts.length - 1);
    const coeffs = polyfit(xArr, yArr, actualOrder);
    const wlFine = linspace(md.wl[0], md.wl[md.wl.length - 1], 200);
    const t0Fit = polyval(coeffs, wlFine).map(v => v * 1000);

    // Check if fit curve trace exists
    const existingTraces = chirpEl.data;
    const fitTraceIdx = existingTraces.findIndex(t => t.name && t.name.startsWith('手动拟合'));
    if (fitTraceIdx >= 0) {
      Plotly.update(chirpEl, {
        x: [wlFine],
        y: [t0Fit]
      }, fixedLayout, [fitTraceIdx]);
    } else {
      Plotly.addTraces(chirpEl, {
        x: wlFine, y: t0Fit, mode: 'lines',
        name: `手动拟合 (${actualOrder}阶)`,
        line: { color: '#ff0000', width: 2 }
      });
      if (axisRange) {
        Plotly.relayout(chirpEl, fixedLayout);
      }
    }
  } else {
    // Remove fit curve if not enough points
    const existingTraces = chirpEl.data;
    const fitTraceIdx = existingTraces.findIndex(t => t.name && t.name.startsWith('手动拟合'));
    if (fitTraceIdx >= 0) {
      Plotly.deleteTraces(chirpEl, [fitTraceIdx]);
    }
  }
}

function undoLastManualPoint(divId) {
  const md = window._manualChirpData?.[divId];
  if (!md || md.points.length === 0) return;
  md.points.pop();

  // Also remove fit curve trace if it exists
  const chirpEl = $(`${divId}_chirpCurve`);
  if (chirpEl) {
    const fitTraceIdx = chirpEl.data.findIndex(t => t.name && t.name.startsWith('手动拟合'));
    if (fitTraceIdx >= 0) {
      Plotly.deleteTraces(chirpEl, [fitTraceIdx]);
    }
  }

  updateManualChirpPlot(divId);
}

function clearManualPoints(divId) {
  const md = window._manualChirpData?.[divId];
  if (!md) return;
  md.points = [];

  // Remove fit curve trace if it exists
  const chirpEl = $(`${divId}_chirpCurve`);
  if (chirpEl && chirpEl.data.length > 3) {
    Plotly.deleteTraces(chirpEl, [3]);
  }

  updateManualChirpPlot(divId);
}

function applyManualChirp(baseName, divId) {
  const md = window._manualChirpData?.[divId];
  if (!md || md.points.length < 3) return;
  try {

  const pts = md.points;
  const fitOrder = parseInt($('manualFitOrder')?.value || 3);
  const actualOrder = Math.min(fitOrder, pts.length - 1);

  // Fit polynomial to user-selected points
  const xArr = pts.map(p => p.wl);
  const yArr = pts.map(p => p.t0 / 1000); // fs -> ps
  const coeffs = polyfit(xArr, yArr, actualOrder);

  // Apply chirp correction using the fitted coefficients
  const data = window[`data_${baseName}`];
  if (!data) return;

  var taAfter = applyChirpShift(data.timeArray, data.wavelengthArray, data.taBefore, coeffs);
  if (!taAfter || taAfter.length === 0) {
    console.warn('[applyManualChirp] applyChirpShift returned empty, using taBefore');
    taAfter = data.taBefore;
  }

  // Update stored data
  data.TA2D = taAfter;
  data.coeffs = coeffs;

  const tViewMin = parseFloat($('tViewMin').value);
  const tViewMax = parseFloat($('tViewMax').value);
  const tRange = [tViewMin, tViewMax];

  // Update chirp compare
  const beforeHm = makeHeatmapData(data.timeArray, data.wavelengthArray, data.taBefore, tRange);
  const afterHm = makeHeatmapData(data.timeArray, data.wavelengthArray, taAfter, tRange);
  let hmZMin = Infinity, hmZMax = -Infinity;
  for (const row of beforeHm.z) { for (const v of row) { if (v !== null && v < hmZMin) hmZMin = v; if (v !== null && v > hmZMax) hmZMax = v; } }
  for (const row of afterHm.z) { for (const v of row) { if (v !== null && v < hmZMin) hmZMin = v; if (v !== null && v > hmZMax) hmZMax = v; } }
  const hmZRange = Math.max(Math.abs(hmZMin), Math.abs(hmZMax));
  Plotly.react($(`${divId}_chirpCompare`), [
    { z: beforeHm.z, x: data.wavelengthArray, y: beforeHm.tPlot, type: 'heatmap', colorscale: 'RdBu', reversescale: true, zmid: 0, zmin: -hmZRange, zmax: hmZRange, name: '校正前', showscale: false },
    { z: afterHm.z, x: data.wavelengthArray, y: afterHm.tPlot, type: 'heatmap', colorscale: 'RdBu', reversescale: true, zmid: 0, zmin: -hmZRange, zmax: hmZRange, name: '校正后', xaxis: 'x2', yaxis: 'y2', showscale: false }
  ], {
    title: '校正前后对比',
    grid: { columns: 2, rows: 1, pattern: 'independent', xgap: 0.08 },
    margin: { l: 60, r: 20, t: 40, b: 50 },
    xaxis: { title: '波长 (nm)', domain: [0, 0.46] }, yaxis: { title: '时间 (ps)' },
    xaxis2: { title: '波长 (nm)', anchor: 'y2' }, yaxis2: { title: '时间 (ps)', anchor: 'x2' },
    annotations: [
      { text: '校正前', x: 0.22, y: 1.05, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 12 } },
      { text: '校正后', x: 0.78, y: 1.05, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 12 } }
    ]
  });

  // Update kinetic traces (after)
  const probeWlStr = $('probeWl').value;
  const probeWavelengths = probeWlStr.split(',').map(s => parseFloat(s.trim())).filter(v => !isNaN(v));
  const kinColors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'];
  const kinTracesAfter = [];
  for (let pi = 0; pi < probeWavelengths.length; pi++) {
    const pw = probeWavelengths[pi];
    let idxWl = 0, minDiff = Infinity;
    for (let i = 0; i < data.wavelengthArray.length; i++) {
      if (Math.abs(data.wavelengthArray[i] - pw) < minDiff) { minDiff = Math.abs(data.wavelengthArray[i] - pw); idxWl = i; }
    }
    const sig = taAfter[idxWl];
    const sigPlot = [], tPlot = [];
    for (let j = 0; j < data.timeArray.length; j++) {
      if (data.timeArray[j] >= tViewMin && data.timeArray[j] <= tViewMax && !isNaN(sig[j])) { sigPlot.push(sig[j]); tPlot.push(data.timeArray[j]); }
    }
    kinTracesAfter.push({ x: tPlot, y: sigPlot, mode: 'lines', name: `${data.wavelengthArray[idxWl].toFixed(1)} nm`, line: { color: kinColors[pi % kinColors.length] } });
  }
  Plotly.react($(`${divId}_afterKin`), kinTracesAfter, {
    xaxis: { title: '时间 (ps)' }, yaxis: { title: 'ΔA (mOD)' },
    legend: { font: { size: 10 } }, margin: { l: 60, r: 20, t: 40, b: 50 }, title: '校正后'
  }, { responsive: true });

  // Update chirp curve fit style
  const chirpEl = $(`${divId}_chirpCurve`);
  if (chirpEl && chirpEl.data.length > 3) {
    Plotly.restyle(chirpEl, { line: [{ color: '#00cc00', width: 3 }], name: [`手动拟合 (${actualOrder}阶) ✓`] }, [3]);
  }

  // Update info
  const invL2display = data.wavelengthArray.map(w => 1.0 / (w * w));
  const t0Fitted = polyval(coeffs, invL2display);
  const t0Fs = t0Fitted.map(v => v * 1000);
  const ref = t0Fs.reduce((a, b) => a + b, 0) / t0Fs.length;
  $(`${divId}_info`).innerHTML = `参考t0 = ${ref.toFixed(1)} fs | 啁啾范围: ${Math.min(...t0Fs).toFixed(1)} ~ ${Math.max(...t0Fs).toFixed(1)} fs | 方法: 手动选点法 (${actualOrder}阶) | ✅ 已应用`;

  // Visual feedback
  const panel = $(`${divId}_manualPanel`);
  if (panel) panel.style.borderColor = '#0f0';
  } catch(e) {
    const infoEl = $(`${divId}_info`);
    if (infoEl) infoEl.innerHTML = `<span style="color:#dc3545">❌ 手动校正失败：${escapeHtml(e.message)}</span>`;
    console.error('[applyManualChirp]', e);
  }
}

async function fitMultiExp(time, signal, nExp, tFitMin, tFitMax, quality) {
  try {
    var r = window.taWasm.fit_multi_exp(
      new Float64Array(time), new Float64Array(signal), nExp, tFitMin, tFitMax, quality
    );
    if (!r) return null;
    return { params: Array.from(r.params), stdErrs: Array.from(r.std_errs), r2: r.r2, t_fit: Array.from(r.t_fit), y_fit: Array.from(r.y_fit) };
  } catch(e) {
    if (e instanceof WebAssembly.RuntimeError) {
      _wasmReady = false;
      throw new Error('WASM 运行时崩溃，请刷新页面重试（' + e.message + '）');
    }
    throw e;
  }
}

function updateFitTimeScale(divId) {
  const plotEl = $(`${divId}_fitPlot`);
  if (!plotEl) return;
  const scale = $(`${divId}_fitTimeScale`).value;
  const shapes = scale === 'linear'
    ? [{ type: 'line', x0: 0, x1: 0, y0: 0, y1: 1, yref: 'paper', line: { color: 'gray', dash: 'dot' } }]
    : [];
  try { Plotly.relayout(plotEl, { 'xaxis.type': scale, shapes }); } catch(e) {}
}

async function doKineticFit(baseName, divId) {
  if (!_wasmReady) { setStatus('❌ WASM 核心尚未加载，请刷新页面重试', 'error'); return; }
  const data = window[`data_${baseName}`];
  if (!data) return;

  const wlInput = $(`${divId}_fitWlInput`).value;
  const wavelengths = wlInput.split(',').map(s => parseFloat(s.trim())).filter(v => !isNaN(v));
  const nExp = parseInt($(`${divId}_fitNExp`).value);
  const tFitMin = Math.max(parseFloat($(`${divId}_fitTMin`).value), data.timeArray[0]);
  const tFitMax = Math.min(parseFloat($(`${divId}_fitTMax`).value), data.timeArray[data.timeArray.length - 1]);
  const timeScale = $(`${divId}_fitTimeScale`).value;
  const fitQuality = parseInt($(`${divId}_fitQuality`).value);

  if (wavelengths.length === 0) {
    $(`${divId}_fitResult`).innerHTML = '<div class="status status-error">请输入有效的探测波长</div>';
    return;
  }

  const fitResultDiv = $(`${divId}_fitResult`);
  const fitBtn = fitResultDiv.closest('.tab-content')?.querySelector('.btn-primary');
  if (fitBtn) fitBtn.disabled = true;

  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;

  const colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'];
  const traces = [];
  let resultHtml = '';
  const fitResultsStore = [];

  try {
  for (let pi = 0; pi < wavelengths.length; pi++) {
    const pw = wavelengths[pi];
    let idxWl = 0;
    let minDiff = Infinity;
    for (let i = 0; i < wl.length; i++) {
      if (Math.abs(wl[i] - pw) < minDiff) { minDiff = Math.abs(wl[i] - pw); idxWl = i; }
    }
    const actualWl = wl[idxWl];
    const signal = ta[idxWl];

    fitResultDiv.innerHTML = `<div class="status status-info">⏳ 拟合中 (${pi + 1}/${wavelengths.length})：${actualWl.toFixed(1)} nm，${nExp} 指数…</div>`;
    await new Promise(r => setTimeout(r, 0));

    const tMaskIdx = [];
    const tPlot = [];
    for (let j = 0; j < time.length; j++) {
      if (time[j] >= tFitMin && time[j] <= tFitMax) { tMaskIdx.push(j); tPlot.push(time[j]); }
    }
    const sigPlot = tMaskIdx.map(j => signal[j] * 1000);

    const color = colors[pi % colors.length];
    traces.push({
      x: tPlot, y: sigPlot,
      mode: 'markers', name: `${actualWl}nm 数据`,
      marker: { color, size: 5 }
    });

    let fitResult;
    try {
      fitResult = await fitMultiExp(time, signal, nExp, tFitMin, tFitMax, fitQuality);
    } catch(e) {
      resultHtml += `<div style="color:#dc3545;margin-bottom:8px;">${actualWl.toFixed(1)}nm: ❌ ${escapeHtml(e.message)}</div>`;
      if (!_wasmReady) break;
      continue;
    }
    if (!fitResult) {
      resultHtml += `<div style="color:#dc3545;margin-bottom:8px;">${actualWl.toFixed(1)}nm: 数据点不足，无法拟合</div>`;
      continue;
    }

    const tFine = fitResult.t_fit;
    const yFine = fitResult.y_fit.map(v => v * 1000);
    const fitFlags = getFitFlags({ ...fitResult, nExp });
    traces.push({
      x: tFine, y: yFine,
      mode: 'lines', name: `${actualWl}nm 拟合`,
      line: { color, width: 2 }
    });

    let formula = `ΔA(t) = `;
    const terms = [];
    for (let k = 0; k < nExp; k++) {
      const A = fitResult.params[k * 2];
      const tau = fitResult.params[k * 2 + 1];
      const sign = A >= 0 ? (k === 0 ? '' : ' + ') : (k === 0 ? '-' : ' - ');
      terms.push(`${sign}${Math.abs(A * 1000).toFixed(2)}·exp(-t/${tau.toFixed(3)})`);
    }
    const offset = fitResult.params[fitResult.params.length - 1];
    const offSign = offset >= 0 ? ' + ' : ' - ';
    formula += terms.join('');
    formula += `${offSign}${Math.abs(offset * 1000).toFixed(2)}`;

    const hasStd = fitResult.stdErrs && fitResult.stdErrs.length > 0;
    let paramRows = '';
    for (let k = 0; k < nExp; k++) {
      const A = fitResult.params[k * 2];
      const tau = fitResult.params[k * 2 + 1];
      let tauStr = `${tau.toFixed(4)} ps`;
      let aStr = `${(A * 1000).toFixed(3)} mOD`;
      if (hasStd) {
        const tauStd = fitResult.stdErrs[k * 2 + 1];
        const aStd = fitResult.stdErrs[k * 2];
        tauStr += ` ± ${tauStd.toFixed(4)}`;
        aStr += ` ± ${(aStd * 1000).toFixed(3)}`;
      }
      paramRows += `<tr><td>τ${k + 1}</td><td>${tauStr}</td><td>A${k + 1}</td><td>${aStr}</td></tr>`;
    }
    let offStr = `${(offset * 1000).toFixed(3)} mOD`;
    if (hasStd) {
      const offStd = fitResult.stdErrs[fitResult.stdErrs.length - 1];
      offStr += ` ± ${(offStd * 1000).toFixed(3)}`;
    }
    paramRows += `<tr><td>y0</td><td colspan="3">${offStr}</td></tr>`;

    resultHtml += `
      <div style="margin-bottom:12px;padding:12px;background:#ffffff;border-radius:6px;border-left:4px solid ${color};color:#222;">
        <strong style="color:${color};">${actualWl} nm</strong>
        <table style="width:100%;font-size:12px;margin-top:6px;border-collapse:collapse;color:#222;">
          <tr style="font-weight:600;color:#555;"><td>参数</td><td>值</td><td>参数</td><td>值</td></tr>
          ${paramRows}
        </table>
        <div style="margin-top:4px;font-size:12px;color:#333;">R² = ${fitResult.r2.toFixed(6)}</div>
        ${fitFlags.length > 0 ? `<div style="margin-top:4px;font-size:12px;color:#b26a00;">异常标记: ${fitFlags.join(', ')}</div>` : ''}
        <div style="margin-top:2px;font-size:11px;color:#555;font-family:monospace;">${formula}</div>
      </div>`;

    fitResultsStore.push({ wavelength: actualWl, nExp, params: fitResult.params, stdErrs: fitResult.stdErrs, r2: fitResult.r2, tData: tPlot, yData: sigPlot, tFit: tFine, yFit: yFine, flags: fitFlags });
  }
  } catch(e) {
    resultHtml += `<div class="status status-error">❌ 拟合中断：${escapeHtml(e.message)}</div>`;
  } finally {
    if (fitBtn) fitBtn.disabled = false;
  }

  if (traces.length > 0) {
    Plotly.newPlot($(`${divId}_fitPlot`), traces, {
      xaxis: { title: '时间 (ps)', type: timeScale },
      yaxis: { title: 'ΔA (mOD)' },
      title: `${nExp}指数拟合`,
      margin: { l: 60, r: 20, t: 40, b: 50 },
      shapes: timeScale === 'linear' ? [{ type: 'line', x0: 0, x1: 0, y0: 0, y1: 1, yref: 'paper', line: { color: 'gray', dash: 'dot' } }] : [],
      legend: { font: { size: 10 } }
    }, { responsive: true });
  }

  window[`fitResults_${baseName}`] = fitResultsStore;
  $(`${divId}_fitResult`).innerHTML = resultHtml || '<div class="status status-error">所有波长均无法拟合</div>';
}

function applyGlobalFitOptionsToAll() {
  const cards = Array.from(document.querySelectorAll('#resultsArea .card[data-fit-card="1"]'));
  const wlInput = $('globalFitWlInput')?.value ?? '';
  const nExp = $('globalFitNExp')?.value ?? '2';
  const timeScale = $('globalFitTimeScale')?.value ?? 'log';
  const fitQuality = $('globalFitQuality')?.value ?? '0';

  for (const card of cards) {
    const divId = card.getAttribute('data-div-id');
    const wlEl = $(`${divId}_fitWlInput`);
    const nExpEl = $(`${divId}_fitNExp`);
    const scaleEl = $(`${divId}_fitTimeScale`);
    const qualityEl = $(`${divId}_fitQuality`);
    if (wlEl) wlEl.value = wlInput;
    if (nExpEl) nExpEl.value = nExp;
    if (scaleEl) scaleEl.value = timeScale;
    if (qualityEl) qualityEl.value = fitQuality;
    updateFitTimeScale(divId);
  }

  const status = $('globalFitStatus');
  if (status) status.textContent = `已将全局拟合选项应用到 ${cards.length} 个文件`;
}

function getFitFlags(result) {
  const flags = [];
  if (!result) return flags;
  if (result.r2 < 0.98) flags.push('R2_LOW');
  if (result.r2 < 0.9) flags.push('R2_VERY_LOW');
  if (result.stdErrs && result.stdErrs.length === result.params.length) {
    for (let k = 0; k < result.nExp; k++) {
      const amp = Math.abs(result.params[k * 2]);
      const tau = Math.abs(result.params[k * 2 + 1]);
      const ampErr = Math.abs(result.stdErrs[k * 2] || 0);
      const tauErr = Math.abs(result.stdErrs[k * 2 + 1] || 0);
      if (amp > 1e-12 && ampErr > amp) flags.push(`A${k + 1}_UNSTABLE`);
      if (tau > 1e-12 && tauErr > tau) flags.push(`T${k + 1}_UNSTABLE`);
    }
    for (let k = 1; k < result.nExp; k++) {
      const prevTau = Math.abs(result.params[(k - 1) * 2 + 1]);
      const tau = Math.abs(result.params[k * 2 + 1]);
      if (prevTau > 1e-12 && tau / prevTau < 1.5) flags.push(`TAU_CLOSE_${k}_${k + 1}`);
    }
  }
  return Array.from(new Set(flags));
}

async function doAllKineticFits() {
  if (!_wasmReady) {
    setStatus('❌ WASM 核心尚未加载，请刷新页面重试', 'error');
    return;
  }

  const cards = Array.from(document.querySelectorAll('#resultsArea .card[data-fit-card="1"]'));
  if (cards.length === 0) return;

  const globalBtn = $('globalFitBtn');
  const globalCancelBtn = $('globalFitCancelBtn');
  const globalStatus = $('globalFitStatus');
  cancelGlobalKineticFit = false;
  if (globalBtn) globalBtn.disabled = true;
  if (globalCancelBtn) globalCancelBtn.style.display = 'inline-block';

  try {
    for (let i = 0; i < cards.length; i++) {
      if (cancelGlobalKineticFit) break;
      const card = cards[i];
      const baseName = card.getAttribute('data-base-name');
      const divId = card.getAttribute('data-div-id');
      const fileTitle = card.querySelector('h2')?.textContent?.trim() || baseName;
      if (globalStatus) globalStatus.textContent = `正在拟合 ${i + 1}/${cards.length}: ${fileTitle}`;
      await doKineticFit(baseName, divId);
    }
    if (globalStatus) {
      globalStatus.textContent = cancelGlobalKineticFit
        ? '批量拟合已取消'
        : `全部拟合完成，共 ${cards.length} 个文件`;
    }
  } catch (e) {
    if (globalStatus) globalStatus.textContent = `批量拟合中断: ${e.message}`;
  } finally {
    if (globalBtn) globalBtn.disabled = false;
    if (globalCancelBtn) globalCancelBtn.style.display = 'none';
    cancelGlobalKineticFit = false;
  }
}

function downloadAllFitSummaryCSV() {
  const cards = Array.from(document.querySelectorAll('#resultsArea .card[data-fit-card="1"]'));
  const rows = [];
  let maxNExp = 0;

  for (const card of cards) {
    const baseName = card.getAttribute('data-base-name');
    const results = window[`fitResults_${baseName}`];
    if (!results || results.length === 0) continue;
    for (const r of results) {
      maxNExp = Math.max(maxNExp, r.nExp);
      rows.push({ baseName, result: r });
    }
  }
  if (rows.length === 0) return;

  const headers = ['文件', '波长(nm)', '指数个数', 'R²'];
  for (let k = 0; k < maxNExp; k++) {
    headers.push(`τ${k + 1}(ps)`, `τ${k + 1}标准差`, `A${k + 1}(mOD)`, `A${k + 1}标准差(mOD)`);
  }
  headers.push('y0(mOD)', 'y0标准差(mOD)', '异常标记');

  let csv = headers.join(',') + '\n';
  for (const item of rows) {
    const r = item.result;
    const row = [item.baseName, r.wavelength.toFixed(2), r.nExp, r.r2.toFixed(6)];
    for (let k = 0; k < maxNExp; k++) {
      if (k < r.nExp) {
        const tau = r.params[k * 2 + 1];
        const A = r.params[k * 2] * 1000;
        const tauStd = r.stdErrs?.[k * 2 + 1];
        const aStd = r.stdErrs?.[k * 2];
        row.push(
          tau.toFixed(6),
          tauStd != null ? tauStd.toFixed(6) : '',
          A.toFixed(4),
          aStd != null ? (aStd * 1000).toFixed(4) : ''
        );
      } else {
        row.push('', '', '', '');
      }
    }
    const off = r.params[r.params.length - 1] * 1000;
    const offStd = r.stdErrs?.[r.stdErrs.length - 1];
    row.push(off.toFixed(4), offStd != null ? (offStd * 1000).toFixed(4) : '', (r.flags || []).join('; '));
    csv += row.join(',') + '\n';
  }

  triggerDownload('\uFEFF' + csv, 'all_fit_summary.csv', 'text/csv;charset=utf-8');
}

async function downloadAllByTypeZip(type) {
  if (typeof JSZip === 'undefined') {
    setStatus('❌ JSZip 未加载，无法打包 ZIP', 'error');
    return;
  }

  const cards = Array.from(document.querySelectorAll('#resultsArea .card[data-fit-card="1"]'));
  const zip = new JSZip();
  let count = 0;

  for (const card of cards) {
    const baseName = card.getAttribute('data-base-name');
    let file = null;
    if (type === 'processed_csv') file = buildProcessedCSV(baseName);
    else if (type === 'fit_csv') file = buildFitCSV(baseName);
    if (!file) continue;
    zip.file(file.filename, file.content);
    count++;
  }

  if (count === 0) {
    setStatus('❌ 当前没有可打包的数据', 'error');
    return;
  }

  const blob = await zip.generateAsync({ type: 'blob' });
  const filename = type === 'processed_csv' ? 'all_processed_csv.zip' : 'all_fit_csv.zip';
  triggerBlobDownload(blob, filename);
}

const _downloadFitCSV_original = downloadFitCSV;
downloadFitCSV = function(baseName) {
  const file = buildFitCSV(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
};

const _downloadFitRawCSV_original = downloadFitRawCSV;
downloadFitRawCSV = function(baseName) {
  const file = buildFitRawCSV(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
};

const _downloadSliceCSV_original = downloadSliceCSV;
downloadSliceCSV = function(baseName) {
  const file = buildSliceCSV(baseName);
  if (!file) return;
  triggerDownload(file.content, file.filename, file.mimeType);
};
