# ta-preprocess-web

Web-based 瞬态吸收 (TA) 数据预处理工具。支持 UFS/CSV 格式导入、啁啾校正、基线扣除、SNR 过滤等。

## 技术栈

- **前端**：纯 HTML/JS（无框架），Plotly.js 绘图
- **核心算法**：Rust → WASM（`core/` 目录），计算密集型运算走 WASM，启动失败则自动 fallback 到 JS 实现
- **部署**：GitHub Pages（静态站点，`index.html` + `app.js`）

## 项目结构

```
ta-preprocess-web/
├── index.html              # 主页面
├── app.js                  # 全部前端逻辑 + JS fallback 算法
├── core/                   # Rust WASM 核心
│   ├── src/
│   │   ├── lib.rs          # WASM 接口导出
│   │   ├── chirp.rs        # 啁啾校正（半高点法 + 全局优化法）
│   │   ├── cpm.rs          # CPM 相干伪影去除
│   │   ├── irf.rs          # IRF 解卷积（Lucy-Richardson）
│   │   ├── baseline.rs     # 基线扣除、polyfit/polyval 工具函数
│   │   ├── fitting.rs      # Nelder-Mead 优化、多指数拟合
│   │   ├── ufs.rs          # UFS 文件解析
│   │   └── csv_parser.rs   # CSV 解析
│   └── Cargo.toml
├── build.ps1 / build.sh    # WASM 编译脚本
├── .claude/                # 权限配置
└── 朱海明..._cleaned.md     # TA 光谱综述（参考用）
```

## WASM 编译与部署

```bash
# 一键编译（输出直接到根目录 pkg/）
.\build.ps1              # PowerShell
./build.sh               # Bash

# 编译并推送
.\build.ps1 -Push -Message "update WASM"
```

当前 GitHub Pages 上无 WASM 文件，走 JS fallback 路径。JS 与 WASM 算法逻辑保持一致。

## 啁啾校正算法

拟合公式采用物理色散模型：

```
t₀(λ) = c₀ + c₁·(1/λ²) + c₂·(1/λ⁴)
```

拟合在 λ⁻² 空间进行，符合材料群延迟色散（GDD）的物理本质。三种方法：

| 方法 | 说明 | 实现位置（Rust） | 实现位置（JS fallback） |
|------|------|-----------------|----------------------|
| 半高点法 | 基于半高宽检测 t₀，二阶多项式拟合 | `chirp.rs::chirp_correction_half_height` | `app.js::chirpCorrectionHalfHeight` |
| 全局优化法 | Nelder-Mead 最大化信号锐度 | `chirp.rs::chirp_correction_global` | `app.js::chirpCorrectionGlobal` |
| 手动选点 | 用户交互点击标记 t₀ 点 | — | `app.js` sigmaClip 互动 |

## 预处理流程

```
UFS/CSV 导入 → 探针波长裁剪 → 基线扣除 → 啁啾校正 → [CPM去除] → [IRF解卷积] → SNR 过滤 → 可选手动选点 → 结果可视化 + 下载
```

CPM 去除和 IRF 解卷积为可选步骤，默认关闭，勾选后启用。

## CPM 相干伪影去除

CPM (Cross-Phase Modulation) 信号在 t₀ 附近表现为高斯+导数高斯形状的相干尖峰，混在真实 ΔA 信号里。

- 模型：`y(t) = A·exp(-t²/(2σ²)) + B·(-t/σ²)·exp(-t²/(2σ²)) + C`
- 使用 Nelder-Mead 拟合，然后减去 CPM 分量（高斯+导数高斯），保留基线 C
- 默认拟合窗口 ±0.5 ps
- 实现位置：`cpm.rs::remove_cpm` / `app.js::removeCpm`

## IRF 解卷积

测量信号 = 真实动力学 ⊗ IRF（高斯卷积）。解卷积恢复真实动力学。

- IRF 宽度可手动输入或自动估算（从 t₀ 附近上升沿的导数拟合）
- 使用 Lucy-Richardson 迭代法（默认 15 次），比傅里叶除法更稳定
- 实现位置：`irf.rs::deconvolve_irf` / `app.js::deconvolveIrf`

## 支持的文件格式

- **UFS**：飞秒 TA 实验专有格式，多文件载入（参考、信号、延迟时间自动匹配）
- **CSV**：通用格式，包含 time/wavelength/ΔA 三列数据
- 标准数据文件：`paper/bibliography/pdf/` 下的 PDF 文献通过 `batch_convert.bat` 批量转为清洗后 Markdown

## 2026-06-09 变更

- **啁啾拟合公式修复**：从 `t₀ = c₀ + c₁·λ + c₂·λ²`（无物理意义）改为 `t₀ = c₀ + c₁·(1/λ²) + c₂·(1/λ⁴)`（材料色散模型）
- **零点对齐修复**：校正后 t₀ 对齐到 t=0（之前是对齐到平均 t₀）
- **WASM 啁啾校正参数化**：之前硬编码，现在接受用户参数（搜索范围、SNR 阈值等）
- **默认多项式阶数改为 2 阶**：λ⁻² 空间下 2 阶足够
- **新增 CPM 相干伪影去除**：可选步骤，高斯+导数高斯模型拟合扣除
- **新增 IRF 解卷积**：可选步骤，Lucy-Richardson 迭代法
- **热力图共享色标**：校正前后热力图使用同一颜色范围
- **构建脚本**：`build.ps1` / `build.sh` 一键编译+部署
