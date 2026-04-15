let uploadedFiles = [];

function $(id) { return document.getElementById(id); }

function enterApp() {
  $('coverPage').classList.add('hidden');
}

document.addEventListener('DOMContentLoaded', function() {
  const enterBtn = document.getElementById('enterBtn');
  if (enterBtn) enterBtn.addEventListener('click', enterApp);
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

function handleFiles(files) {
  for (let f of files) {
    if (!uploadedFiles.find(u => u.name === f.name)) {
      uploadedFiles.push(f);
    }
  }
  renderFileList();
  $('processBtn').disabled = uploadedFiles.length === 0;
}

function renderFileList() {
  let html = '';
  uploadedFiles.forEach((f, i) => {
    const sizeMB = (f.size / 1024 / 1024).toFixed(2);
    html += `<div class="file-item"><span class="name">📄 ${f.name}</span><span class="size">${sizeMB} MB</span><span><button class="btn" style="padding:4px 10px;font-size:12px;background:#dc3545;color:white;" onclick="removeFile(${i})">✕</button></span></div>`;
  });
  $('fileList').innerHTML = html;
}

function removeFile(idx) {
  uploadedFiles.splice(idx, 1);
  renderFileList();
  $('processBtn').disabled = uploadedFiles.length === 0;
}

function setStatus(msg, type, progress) {
  const bar = progress !== undefined
    ? `<div class="progress-bar"><div class="fill" style="width:${Math.min(100, Math.max(0, progress)).toFixed(1)}%"></div></div>`
    : '';
  $('statusArea').innerHTML = `<div class="status status-${type}">${msg}${bar}</div>`;
}

function parseCSV(text) {
  const lines = text.split(/\r?\n/);
  const validLines = [];
  for (let line of lines) {
    const stripped = line.trim();
    if (!stripped) continue;
    const firstVal = stripped.split(',')[0].trim();
    if (isNaN(parseFloat(firstVal))) break;
    validLines.push(stripped.split(',').map(v => v.trim() === '' ? NaN : parseFloat(v)));
  }
  const nRows = validLines.length;
  const nCols = validLines[0].length;
  const timeArray = validLines[0].slice(1);
  const wavelengthArray = [];
  const TA2D = [];
  for (let i = 1; i < nRows; i++) {
    wavelengthArray.push(validLines[i][0]);
    TA2D.push(validLines[i].slice(1));
  }
  return { timeArray, wavelengthArray, TA2D };
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

function polyfit(x, y, order) {
  const n = x.length;
  const m = order + 1;
  const A = [];
  const b = [];
  for (let i = 0; i < m; i++) {
    A[i] = [];
    for (let j = 0; j < m; j++) {
      let sum = 0;
      for (let k = 0; k < n; k++) sum += Math.pow(x[k], i + j);
      A[i][j] = sum;
    }
    let sum = 0;
    for (let k = 0; k < n; k++) sum += y[k] * Math.pow(x[k], i);
    b[i] = sum;
  }
  for (let col = 0; col < m; col++) {
    let maxRow = col;
    for (let row = col + 1; row < m; row++) {
      if (Math.abs(A[row][col]) > Math.abs(A[maxRow][col])) maxRow = row;
    }
    [A[col], A[maxRow]] = [A[maxRow], A[col]];
    [b[col], b[maxRow]] = [b[maxRow], b[col]];
    for (let row = col + 1; row < m; row++) {
      const factor = A[row][col] / A[col][col];
      for (let j = col; j < m; j++) A[row][j] -= factor * A[col][j];
      b[row] -= factor * b[col];
    }
  }
  const coeffs = new Array(m);
  for (let i = m - 1; i >= 0; i--) {
    coeffs[i] = b[i];
    for (let j = i + 1; j < m; j++) coeffs[i] -= A[i][j] * coeffs[j];
    coeffs[i] /= A[i][i];
  }
  return coeffs;
}

function polyval(coeffs, x) {
  return x.map(xi => {
    let val = 0;
    for (let i = 0; i < coeffs.length; i++) val += coeffs[i] * Math.pow(xi, i);
    return val;
  });
}

function linspace(min, max, n) {
  const step = (max - min) / (n - 1);
  return Array.from({ length: n }, (_, i) => min + i * step);
}

function interp1d(xData, yData, xNew) {
  return xNew.map(xn => {
    if (isNaN(xn)) return NaN;
    if (xn < xData[0] || xn > xData[xData.length - 1]) return NaN;
    let lo = 0, hi = xData.length - 1;
    while (hi - lo > 1) {
      const mid = Math.floor((lo + hi) / 2);
      if (xData[mid] <= xn) lo = mid; else hi = mid;
    }
    const t = (xn - xData[lo]) / (xData[hi] - xData[lo]);
    return yData[lo] + t * (yData[hi] - yData[lo]);
  });
}

function findT0HalfHeight(time, ta, searchRange) {
  const t0List = [];
  const tSearchMask = [];
  const tSearch = [];
  for (let j = 0; j < time.length; j++) {
    if (time[j] >= searchRange[0] && time[j] <= searchRange[1]) {
      tSearchMask.push(j);
      tSearch.push(time[j]);
    }
  }
  for (let i = 0; i < ta.length; i++) {
    const sig = tSearchMask.map(j => ta[i][j]);
    const validIdx = sig.map((v, idx) => !isNaN(v) ? idx : -1).filter(idx => idx >= 0);
    if (validIdx.length < 5) { t0List.push(NaN); continue; }
    const tValid = validIdx.map(idx => tSearch[idx]);
    const sValid = validIdx.map(idx => sig[idx]);
    const absSig = sValid.map(v => Math.abs(v));
    let peakIdx = 0, peakVal = 0;
    for (let k = 0; k < absSig.length; k++) {
      if (absSig[k] > peakVal) { peakVal = absSig[k]; peakIdx = k; }
    }
    if (peakVal < 1e-6) { t0List.push(NaN); continue; }
    const halfVal = peakVal * 0.5;
    let crossIdx = -1;
    const searchStart = Math.max(0, peakIdx - 20);
    for (let k = searchStart; k < peakIdx; k++) {
      if (absSig[k] >= halfVal) { crossIdx = k; break; }
    }
    if (crossIdx === -1) { t0List.push(tValid[peakIdx]); continue; }
    if (crossIdx > 0 && crossIdx < tValid.length) {
      const t1 = tValid[crossIdx - 1], t2 = tValid[crossIdx];
      const v1 = absSig[crossIdx - 1], v2 = absSig[crossIdx];
      if (Math.abs(v2 - v1) > 1e-10) {
        t0List.push(t1 + (halfVal - v1) / (v2 - v1) * (t2 - t1));
      } else {
        t0List.push(tValid[crossIdx]);
      }
    } else {
      t0List.push(tValid[crossIdx]);
    }
  }
  return t0List;
}

function applyChirpShift(time, wl, ta, coeffs) {
  const t0Fitted = polyval(coeffs, wl);
  const refT0 = t0Fitted.reduce((a, b) => a + b, 0) / t0Fitted.length;
  const corrected = ta.map((row, i) => {
    const dt = t0Fitted[i] - refT0;
    const timeShifted = time.map(t => t - dt);
    const validIdx = [];
    const validX = [], validY = [];
    for (let j = 0; j < row.length; j++) {
      if (!isNaN(row[j])) { validX.push(timeShifted[j]); validY.push(row[j]); validIdx.push(j); }
    }
    if (validX.length < 3) return row.map(v => NaN);
    return interp1d(validX, validY, time);
  });
  return { TA2D: corrected, t0Fitted, refT0 };
}

async function nelderMead(costFunc, x0, maxIter = 3000) {
  const n = x0.length;
  const alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;
  let simplex = [x0.slice()];
  const f0 = costFunc(x0);
  const fVals = [f0];
  for (let i = 0; i < n; i++) {
    const xi = x0.slice();
    xi[i] += Math.abs(xi[i]) * 0.05 + 1e-6;
    simplex.push(xi);
    fVals.push(costFunc(xi));
  }
  let iter = 0;
  while (iter < maxIter) {
    if (cancelProcessing) break;
    if (iter % 50 === 0) await new Promise(r => setTimeout(r, 0));
    let order = [...Array(n + 1).keys()].sort((a, b) => fVals[a] - fVals[b]);
    const best = order[0], worst = order[n], secondWorst = order[n - 1];
    const xBest = simplex[best], fBest = fVals[best];
    const centroid = new Array(n).fill(0);
    for (let i = 0; i < n + 1; i++) {
      if (i === worst) continue;
      for (let j = 0; j < n; j++) centroid[j] += simplex[i][j] / n;
    }
    const xRef = centroid.map((c, j) => c + alpha * (c - simplex[worst][j]));
    const fRef = costFunc(xRef);
    let shrink = false;
    if (fRef < fVals[secondWorst] && fRef >= fBest) {
      simplex[worst] = xRef; fVals[worst] = fRef;
    } else if (fRef < fBest) {
      const xExp = centroid.map((c, j) => c + gamma * (xRef[j] - c));
      const fExp = costFunc(xExp);
      if (fExp < fRef) { simplex[worst] = xExp; fVals[worst] = fExp; }
      else { simplex[worst] = xRef; fVals[worst] = fRef; }
    } else {
      const xCon = centroid.map((c, j) => c + rho * (simplex[worst][j] - c));
      const fCon = costFunc(xCon);
      if (fCon < fVals[worst]) { simplex[worst] = xCon; fVals[worst] = fCon; }
      else { shrink = true; }
    }
    if (shrink) {
      for (let i = 1; i < n + 1; i++) {
        simplex[i] = xBest.map((b, j) => b + sigma * (simplex[i][j] - b));
        fVals[i] = costFunc(simplex[i]);
      }
    }
    const fMax = Math.max(...fVals), fMin = Math.min(...fVals);
    if (Math.abs(fMax - fMin) < 1e-12) break;
    iter++;
  }
  let bestIdx = 0;
  for (let i = 1; i <= n; i++) { if (fVals[i] < fVals[bestIdx]) bestIdx = i; }
  return { x: simplex[bestIdx], fVal: fVals[bestIdx], iterations: iter };
}

function chirpCorrectionHalfHeight(time, wl, ta) {
  const t0PerWl = findT0HalfHeight(time, ta, [-2.0, 2.0]);
  const validX = [], validY = [];
  for (let i = 0; i < wl.length; i++) {
    if (!isNaN(t0PerWl[i])) { validX.push(wl[i]); validY.push(t0PerWl[i]); }
  }
  if (validX.length < 3) return { TA2D: ta, coeffs: null, t0PerWl };
  const coeffs = polyfit(validX, validY, 2);
  const { TA2D } = applyChirpShift(time, wl, ta, coeffs);
  return { TA2D, coeffs, t0PerWl };
}

async function chirpCorrectionGlobal(time, wl, ta) {
  const t0PerWl = findT0HalfHeight(time, ta, [-2.0, 2.0]);
  const validX = [], validY = [];
  for (let i = 0; i < wl.length; i++) {
    if (!isNaN(t0PerWl[i])) { validX.push(wl[i]); validY.push(t0PerWl[i]); }
  }
  if (validX.length < 3) return { TA2D: ta, coeffs: null, t0PerWl };
  const initialCoeffs = polyfit(validX, validY, 2);

  const nWl = wl.length;
  const validRows = ta.map(row => row.map(v => !isNaN(v)));
  const tEvalIdx = [];
  const tEval = [];
  for (let j = 0; j < time.length; j++) {
    if (time[j] >= -0.5 && time[j] <= 0.5) { tEvalIdx.push(j); tEval.push(time[j]); }
  }
  const dtStep = tEval.length >= 2 ? (tEval[tEval.length - 1] - tEval[0]) / (tEval.length - 1) : 0.01;

  function costFunc(coeffs) {
    const t0Fit = polyval(coeffs, wl);
    const refT0 = t0Fit.reduce((a, b) => a + b, 0) / t0Fit.length;
    const dtShifts = t0Fit.map(t => t - refT0);
    const maxShift = Math.max(...dtShifts.map(Math.abs));
    if (maxShift > 5.0) return 1e10 * (1 + maxShift);
    let totalSharpness = 0, nValidWl = 0;
    for (let i = 0; i < nWl; i++) {
      const dt = dtShifts[i];
      const timeShifted = time.map(t => t - dt);
      const vx = [], vy = [];
      for (let j = 0; j < ta[i].length; j++) {
        if (!isNaN(ta[i][j])) { vx.push(timeShifted[j]); vy.push(ta[i][j]); }
      }
      if (vx.length < 5) continue;
      const sigAtEval = interp1d(vx, vy, tEval);
      const validSig = sigAtEval.filter(v => !isNaN(v));
      if (validSig.length < 2) continue;
      let maxGrad = 0;
      for (let k = 1; k < validSig.length; k++) {
        const grad = Math.abs((validSig[k] - validSig[k - 1]) / dtStep);
        if (grad > maxGrad) maxGrad = grad;
      }
      totalSharpness += maxGrad;
      nValidWl++;
    }
    if (nValidWl === 0) return 1e10;
    let reg = 0;
    for (let k = 0; k < coeffs.length; k++) reg += (coeffs[k] - initialCoeffs[k]) ** 2;
    reg = 0.01 * reg / coeffs.length;
    return -totalSharpness / nValidWl + reg;
  }

  const result = await nelderMead(costFunc, initialCoeffs, 3000);
  const optimalCoeffs = result.x;
  const { TA2D } = applyChirpShift(time, wl, ta, optimalCoeffs);
  return { TA2D, coeffs: optimalCoeffs, t0PerWl };
}

function makeHeatmapData(time, wl, ta, tRange) {
  const tMaskIdx = [];
  const tPlot = [];
  for (let j = 0; j < time.length; j++) {
    if (time[j] >= tRange[0] && time[j] <= tRange[1]) { tMaskIdx.push(j); tPlot.push(time[j]); }
  }
  const z = ta.map(row => tMaskIdx.map(j => isNaN(row[j]) ? null : row[j] * 1000));
  return { tPlot, z };
}

let cancelProcessing = false;

async function processAll() {
  if (uploadedFiles.length === 0) return;
  const wlMin = parseFloat($('wlMin').value);
  const wlMax = parseFloat($('wlMax').value);
  const nBaseline = parseInt($('nBaseline').value);
  const chirpMethod = $('chirpMethod').value;
  const tViewMin = parseFloat($('tViewMin').value);
  const tViewMax = parseFloat($('tViewMax').value);
  const probeWlStr = $('probeWl').value;
  const probeWavelengths = probeWlStr.split(',').map(s => parseFloat(s.trim())).filter(v => !isNaN(v));

  $('resultsArea').innerHTML = '';
  cancelProcessing = false;
  $('processBtn').disabled = true;
  $('processBtn').textContent = '⏹ 取消处理';
  $('processBtn').onclick = () => { cancelProcessing = true; };

  const yieldThread = () => new Promise(r => setTimeout(r, 0));
  const totalFiles = uploadedFiles.length;
  const stepsPerFile = 4;

  for (let fi = 0; fi < totalFiles; fi++) {
    if (cancelProcessing) break;
    const file = uploadedFiles[fi];
    const base = (fi / totalFiles) * 100;
    const step = (1 / totalFiles / stepsPerFile) * 100;

    setStatus(`⏳ [${fi + 1}/${totalFiles}] 读取 ${file.name}...`, 'info', base + step * 0);
    await yieldThread();
    if (cancelProcessing) break;

    const text = await new Promise(resolve => {
      const reader = new FileReader();
      reader.onload = e => resolve(e.target.result);
      reader.readAsText(file, 'latin-1');
    });

    setStatus(`⏳ [${fi + 1}/${totalFiles}] 解析 ${file.name}...`, 'info', base + step * 1);
    await yieldThread();
    if (cancelProcessing) break;

    const { timeArray, wavelengthArray, TA2D } = parseCSV(text);
    const { wavelengthArray: wl, TA2D: taCropped } = cropWavelength(wavelengthArray, TA2D, wlMin, wlMax);
    const taBaseline = baselineSubtraction(timeArray, taCropped, nBaseline);
    const taBeforeChirp = taBaseline.map(r => [...r]);

    setStatus(`⏳ [${fi + 1}/${totalFiles}] 啁啾校正 ${file.name}（可能需要数秒）...`, 'info', base + step * 2);
    await yieldThread();
    if (cancelProcessing) break;

    let chirpResult;
    if (chirpMethod === 'global') {
      chirpResult = await chirpCorrectionGlobal(timeArray, wl, taBaseline);
    } else {
      chirpResult = chirpCorrectionHalfHeight(timeArray, wl, taBaseline);
    }
    const { TA2D: taAfter, coeffs, t0PerWl } = chirpResult;

    if (cancelProcessing) break;
    setStatus(`⏳ [${fi + 1}/${totalFiles}] 渲染结果 ${file.name}...`, 'info', base + step * 3);
    await yieldThread();

    renderResults(file.name, timeArray, wl, taBeforeChirp, taAfter, coeffs, t0PerWl, chirpMethod, tViewMin, tViewMax, probeWavelengths);
  }
  $('processBtn').disabled = false;
  $('processBtn').textContent = '🚀 开始处理';
  $('processBtn').onclick = processAll;
  if (cancelProcessing) {
    setStatus('⏹ 处理已取消', 'error');
    cancelProcessing = false;
  } else {
    setStatus('✅ 全部处理完成!', 'success', 100);
  }
}

function renderResults(fileName, time, wl, taBefore, taAfter, coeffs, t0PerWl, chirpMethod, tViewMin, tViewMax, probeWavelengths) {
  const baseName = fileName.replace('.csv', '').replace(/ /g, '_');
  const divId = `result_${baseName}`;
  const methodName = chirpMethod === 'global' ? '全局优化法' : '半高点法';

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

  let html = `<div class="card" id="${divId}">
    <h2>📄 ${fileName}</h2>
    <p style="font-size:13px;color:#888;">波长: ${wl[0].toFixed(1)} ~ ${wl[wl.length-1].toFixed(1)} nm | 时间: ${time[0].toFixed(3)} ~ ${time[time.length-1].toFixed(3)} ps | 数据: ${wl.length}×${time.length}</p>
    <div class="tabs">
      <div class="tab active" onclick="switchTab(this, '${divId}', 'tabViz')">数据可视化</div>
      <div class="tab" onclick="switchTab(this, '${divId}', 'tabSlice')">光谱切片</div>
      <div class="tab" onclick="switchTab(this, '${divId}', 'tabFit')">动力学拟合</div>
      <div class="tab" onclick="switchTab(this, '${divId}', 'tabDownload')">下载数据</div>
    </div>
    <div class="tab-content active" id="${divId}_tabViz">
      <h3 style="font-size:14px;color:#ff6666;margin-bottom:8px;">2D 伪彩图</h3>
      <div class="plot-grid">
        <div class="plot-box"><h3>校正前</h3><div id="${divId}_before2d" style="width:100%;height:350px;"></div></div>
        <div class="plot-box"><h3>校正后</h3><div id="${divId}_after2d" style="width:100%;height:350px;"></div></div>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">动力学曲线</h3>
      <div class="plot-grid">
        <div class="plot-box"><h3>校正前</h3><div id="${divId}_beforeKin" style="width:100%;height:350px;"></div></div>
        <div class="plot-box"><h3>校正后</h3><div id="${divId}_afterKin" style="width:100%;height:350px;"></div></div>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">啁啾校正</h3>
      <div class="plot-grid">
        <div class="plot-box"><h3>啁啾曲线</h3><div id="${divId}_chirpCurve" style="width:100%;height:350px;"></div></div>
        <div class="plot-box"><h3>校正前后对比</h3><div id="${divId}_chirpCompare" style="width:100%;height:350px;"></div></div>
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
            <select id="${divId}_fitTimeScale">
              <option value="linear">线性</option>
              <option value="log">对数</option>
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
        <button class="btn btn-download" onclick="downloadFitCSV('${baseName}')">⬇️ 拟合结果 (Excel CSV)</button>
      </div>
      <h3 style="font-size:14px;color:#ff6666;margin:16px 0 8px;">光谱切片下载</h3>
      <div class="download-section">
        <button class="btn btn-download" onclick="downloadSliceCSV('${baseName}')">⬇️ 光谱切片 (Excel CSV)</button>
      </div>
      <div style="margin-top:12px;font-size:13px;color:#888;" id="${divId}_info"></div>
    </div>
  </div>`;

  $('resultsArea').insertAdjacentHTML('beforeend', html);

  const tRange = [tViewMin, tViewMax];

  const beforeData = makeHeatmapData(time, wl, taBefore, tRange);
  const afterData = makeHeatmapData(time, wl, taAfter, tRange);

  const heatmapLayout = {
    xaxis: { title: '时间 (ps)' },
    yaxis: { title: '波长 (nm)' },
    margin: { l: 60, r: 20, t: 30, b: 50 },
    colorscale: 'RdBu',
  };

  Plotly.newPlot($(`${divId}_before2d`), [{
    z: beforeData.z, x: beforeData.tPlot, y: wl, type: 'heatmap',
    colorscale: 'RdBu', reversescale: true, zmid: 0,
    colorbar: { title: 'mOD' }
  }], { ...heatmapLayout, title: '校正前' }, { responsive: true });

  Plotly.newPlot($(`${divId}_after2d`), [{
    z: afterData.z, x: afterData.tPlot, y: wl, type: 'heatmap',
    colorscale: 'RdBu', reversescale: true, zmid: 0,
    colorbar: { title: 'mOD' }
  }], { ...heatmapLayout, title: '校正后' }, { responsive: true });

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
  const validWl = [], validT0 = [];
  for (let i = 0; i < wl.length; i++) {
    if (!isNaN(t0PerWl[i])) { validWl.push(wl[i]); validT0.push(t0PerWl[i] * 1000); }
  }
  chirpTraces.push({ x: validWl, y: validT0, mode: 'markers', name: '半高点法t0', marker: { size: 4, opacity: 0.5, color: 'blue' } });

  if (coeffs) {
    const wlFine = linspace(wl[0], wl[wl.length - 1], 200);
    const t0Fit = polyval(coeffs, wlFine).map(v => v * 1000);
    chirpTraces.push({ x: wlFine, y: t0Fit, mode: 'lines', name: chirpMethod === 'global' ? '全局优化拟合' : '多项式拟合', line: { color: 'red', width: 2 } });
    if (chirpMethod === 'global') {
      const initCoeffs = polyfit(validWl, validT0.map(v => v / 1000), 2);
      const t0Init = polyval(initCoeffs, wlFine).map(v => v * 1000);
      chirpTraces.push({ x: wlFine, y: t0Init, mode: 'lines', name: '初始估计拟合', line: { color: 'green', dash: 'dash', width: 1.5 } });
    }
  }

  Plotly.newPlot($(`${divId}_chirpCurve`), chirpTraces, {
    xaxis: { title: '波长 (nm)' }, yaxis: { title: 't0 (fs)' },
    title: '啁啾曲线', margin: { l: 60, r: 20, t: 40, b: 50 }
  }, { responsive: true });

  Plotly.newPlot($(`${divId}_chirpCompare`), [
    { z: beforeData.z, x: beforeData.tPlot, y: wl, type: 'heatmap', colorscale: 'RdBu', reversescale: true, zmid: 0, name: '校正前', showscale: false },
    { z: afterData.z, x: afterData.tPlot, y: wl, type: 'heatmap', colorscale: 'RdBu', reversescale: true, zmid: 0, name: '校正后', xaxis: 'x2', yaxis: 'y2', showscale: false }
  ], {
    title: '校正前后对比',
    grid: { columns: 2, rows: 1, pattern: 'independent', xgap: 0.08 },
    margin: { l: 60, r: 20, t: 40, b: 50 },
    xaxis: { title: '时间 (ps)', domain: [0, 0.46] },
    yaxis: { title: '波长 (nm)' },
    xaxis2: { title: '时间 (ps)', anchor: 'y2' },
    yaxis2: { title: '波长 (nm)', anchor: 'x2' },
    annotations: [
      { text: '校正前', x: 0.22, y: 1.05, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 12 } },
      { text: '校正后', x: 0.78, y: 1.05, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 12 } }
    ]
  }, { responsive: true });

  if (coeffs) {
    const validT0Fs = t0PerWl.filter(v => !isNaN(v)).map(v => v * 1000);
    const refT0 = validT0Fs.reduce((a, b) => a + b, 0) / validT0Fs.length;
    $(`${divId}_info`).innerHTML = `参考t0 = ${refT0.toFixed(1)} fs | 啁啾范围: ${Math.min(...validT0Fs).toFixed(1)} ~ ${Math.max(...validT0Fs).toFixed(1)} fs | 方法: ${methodName}`;
  }

  window[`data_${baseName}`] = { timeArray: time, wavelengthArray: wl, TA2D: taAfter, coeffs, t0PerWl };
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

function triggerDownload(content, filename, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  a.click(); URL.revokeObjectURL(url);
}

function downloadJSON(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return;
  const json = JSON.stringify({
    time_array: data.timeArray,
    wavelength_array: data.wavelengthArray,
    TA_2D_data: data.TA2D,
    chirp_coeffs: data.coeffs
  }, null, 2);
  triggerDownload(json, `${baseName}_processed.json`, 'application/json');
}

function downloadCSV(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return;
  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;
  let csv = ',' + time.map(t => t.toFixed(6)).join(',') + '\n';
  for (let i = 0; i < wl.length; i++) {
    csv += wl[i].toFixed(2) + ',' + ta[i].map(v => isNaN(v) ? '' : v.toExponential(8)).join(',') + '\n';
  }
  triggerDownload(csv, `${baseName}_processed.csv`, 'text/csv');
}

function downloadASCII(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return;
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
  triggerDownload(txt, `${baseName}_processed.dat`, 'text/plain');
}

function downloadTSV(baseName) {
  const data = window[`data_${baseName}`];
  if (!data) return;
  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;
  let tsv = '\t' + time.map(t => t.toFixed(6)).join('\t') + '\n';
  for (let i = 0; i < wl.length; i++) {
    tsv += wl[i].toFixed(2) + '\t' + ta[i].map(v => isNaN(v) ? '' : v.toExponential(8)).join('\t') + '\n';
  }
  triggerDownload(tsv, `${baseName}_processed.tsv`, 'text/tab-separated-values');
}

function downloadFitCSV(baseName) {
  const results = window[`fitResults_${baseName}`];
  if (!results || results.length === 0) return;

  const maxNExp = Math.max(...results.map(r => r.nExp));
  const headers = ['波长(nm)', '指数个数', 'R²'];
  for (let k = 0; k < maxNExp; k++) {
    headers.push(`τ${k + 1}(ps)`, `τ${k + 1}标准差`, `A${k + 1}(mOD)`, `A${k + 1}标准差(mOD)`);
  }
  headers.push('偏移(mOD)', '偏移标准差(mOD)');

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

function multiExpFunc(params, t, nExp) {
  let y = params[params.length - 1];
  for (let k = 0; k < nExp; k++) {
    y += params[k * 2] * Math.exp(-t / params[k * 2 + 1]);
  }
  return y;
}

async function fitMultiExp(time, signal, nExp, tFitMin, tFitMax) {
  const fitIdx = [];
  for (let j = 0; j < time.length; j++) {
    if (time[j] >= tFitMin && time[j] <= tFitMax && !isNaN(signal[j])) {
      fitIdx.push(j);
    }
  }
  if (fitIdx.length < nExp * 2 + 1) return null;

  const tFit = fitIdx.map(j => time[j]);
  const sFit = fitIdx.map(j => signal[j]);

  const sMax = Math.max(...sFit.map(Math.abs));
  const sMean = sFit.reduce((a, b) => a + b, 0) / sFit.length;

  let x0 = [];
  const tauGuesses = [1.0, 0.3, 5.0];
  for (let k = 0; k < nExp; k++) {
    x0.push(sMax / nExp * (k === 0 ? 1 : 0.5));
    x0.push(tauGuesses[k]);
  }
  x0.push(0);

  function costFunc(params) {
    let ss = 0;
    for (let i = 0; i < tFit.length; i++) {
      const yPred = multiExpFunc(params, tFit[i], nExp);
      ss += (yPred - sFit[i]) ** 2;
    }
    for (let k = 0; k < nExp; k++) {
      if (params[k * 2 + 1] < 0.001) ss += 1e6;
    }
    return ss;
  }

  const result = await nelderMead(costFunc, x0, 5000);
  const bestParams = result.x;

  for (let k = 0; k < nExp; k++) {
    if (bestParams[k * 2 + 1] < 0) {
      bestParams[k * 2 + 1] = Math.abs(bestParams[k * 2 + 1]);
      bestParams[k * 2] = -bestParams[k * 2];
    }
  }

  const sorted = [];
  for (let k = 0; k < nExp; k++) {
    sorted.push({ A: bestParams[k * 2], tau: bestParams[k * 2 + 1] });
  }
  sorted.sort((a, b) => a.tau - b.tau);

  const finalParams = [];
  for (const item of sorted) {
    finalParams.push(item.A, item.tau);
  }
  finalParams.push(bestParams[bestParams.length - 1]);

  let ssTot = 0;
  const sMeanVal = sFit.reduce((a, b) => a + b, 0) / sFit.length;
  for (let i = 0; i < sFit.length; i++) ssTot += (sFit[i] - sMeanVal) ** 2;
  const r2 = 1 - result.fVal / ssTot;

  const nP = finalParams.length;
  const nD = tFit.length;
  const sigma2 = result.fVal / (nD - nP);

  const eps = 1e-8;
  const J = [];
  for (let i = 0; i < nD; i++) {
    const row = [];
    for (let j = 0; j < nP; j++) {
      const pPlus = finalParams.slice();
      pPlus[j] += eps;
      const fPlus = multiExpFunc(pPlus, tFit[i], nExp);
      const fOrig = multiExpFunc(finalParams, tFit[i], nExp);
      row.push((fPlus - fOrig) / eps);
    }
    J.push(row);
  }

  const JtJ = [];
  for (let i = 0; i < nP; i++) {
    JtJ[i] = [];
    for (let j = 0; j < nP; j++) {
      let s = 0;
      for (let k = 0; k < nD; k++) s += J[k][i] * J[k][j];
      JtJ[i][j] = s;
    }
  }

  const covInv = JtJ;
  const cov = invertMatrix(covInv);
  const stdErrs = [];
  if (cov) {
    for (let i = 0; i < nP; i++) {
      stdErrs.push(Math.sqrt(Math.max(0, sigma2 * cov[i][i])));
    }
  }

  return { params: finalParams, r2, tFit, sFit, nExp, stdErrs };
}

function invertMatrix(M) {
  const n = M.length;
  const A = M.map((row, i) => row.map((v, j) => v + (i === j ? 1e-10 : 0)));
  const I = [];
  for (let i = 0; i < n; i++) {
    I[i] = new Array(n).fill(0);
    I[i][i] = 1;
  }
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(A[row][col]) > Math.abs(A[maxRow][col])) maxRow = row;
    }
    [A[col], A[maxRow]] = [A[maxRow], A[col]];
    [I[col], I[maxRow]] = [I[maxRow], I[col]];
    if (Math.abs(A[col][col]) < 1e-14) return null;
    const pivot = A[col][col];
    for (let j = 0; j < n; j++) { A[col][j] /= pivot; I[col][j] /= pivot; }
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = A[row][col];
      for (let j = 0; j < n; j++) { A[row][j] -= factor * A[col][j]; I[row][j] -= factor * I[col][j]; }
    }
  }
  return I;
}

async function doKineticFit(baseName, divId) {
  const data = window[`data_${baseName}`];
  if (!data) return;

  const wlInput = $(`${divId}_fitWlInput`).value;
  const wavelengths = wlInput.split(',').map(s => parseFloat(s.trim())).filter(v => !isNaN(v));
  const nExp = parseInt($(`${divId}_fitNExp`).value);
  const tFitMin = Math.max(parseFloat($(`${divId}_fitTMin`).value), data.timeArray[0]);
  const tFitMax = Math.min(parseFloat($(`${divId}_fitTMax`).value), data.timeArray[data.timeArray.length - 1]);
  const timeScale = $(`${divId}_fitTimeScale`).value;

  if (wavelengths.length === 0) {
    $(`${divId}_fitResult`).innerHTML = '<div class="status status-error">请输入有效的探测波长</div>';
    return;
  }

  const time = data.timeArray;
  const wl = data.wavelengthArray;
  const ta = data.TA2D;

  const colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'];
  const traces = [];
  let resultHtml = '';
  const fitResultsStore = [];

  for (let pi = 0; pi < wavelengths.length; pi++) {
    const pw = wavelengths[pi];
    let idxWl = 0;
    let minDiff = Infinity;
    for (let i = 0; i < wl.length; i++) {
      if (Math.abs(wl[i] - pw) < minDiff) { minDiff = Math.abs(wl[i] - pw); idxWl = i; }
    }
    const actualWl = wl[idxWl];
    const signal = ta[idxWl];

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

    const fitResult = await fitMultiExp(time, signal, nExp, tFitMin, tFitMax);
    if (!fitResult) {
      resultHtml += `<div style="color:#dc3545;margin-bottom:8px;">${actualWl}nm: 数据点不足，无法拟合</div>`;
      continue;
    }

    const tFine = linspace(tFitMin, tFitMax, 500);
    const yFine = tFine.map(t => multiExpFunc(fitResult.params, t, nExp) * 1000);
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
    paramRows += `<tr><td>偏移</td><td colspan="3">${offStr}</td></tr>`;

    resultHtml += `
      <div style="margin-bottom:12px;padding:12px;background:#ffffff;border-radius:6px;border-left:4px solid ${color};color:#222;">
        <strong style="color:${color};">${actualWl} nm</strong>
        <table style="width:100%;font-size:12px;margin-top:6px;border-collapse:collapse;color:#222;">
          <tr style="font-weight:600;color:#555;"><td>参数</td><td>值</td><td>参数</td><td>值</td></tr>
          ${paramRows}
        </table>
        <div style="margin-top:4px;font-size:12px;color:#333;">R² = ${fitResult.r2.toFixed(6)}</div>
        <div style="margin-top:2px;font-size:11px;color:#555;font-family:monospace;">${formula}</div>
      </div>`;

    fitResultsStore.push({ wavelength: actualWl, nExp, params: fitResult.params, stdErrs: fitResult.stdErrs, r2: fitResult.r2 });
  }

  Plotly.newPlot($(`${divId}_fitPlot`), traces, {
    xaxis: { title: '时间 (ps)', type: timeScale },
    yaxis: { title: 'ΔA (mOD)' },
    title: `${nExp}指数拟合`,
    margin: { l: 60, r: 20, t: 40, b: 50 },
    shapes: timeScale === 'linear' ? [{ type: 'line', x0: 0, x1: 0, y0: 0, y1: 1, yref: 'paper', line: { color: 'gray', dash: 'dot' } }] : [],
    legend: { font: { size: 10 } }
  }, { responsive: true });

  window[`fitResults_${baseName}`] = fitResultsStore;

  $(`${divId}_fitResult`).innerHTML = resultHtml;
}
