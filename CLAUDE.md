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
│   │   ├── baseline.rs     # 基线扣除、polyfit/polyval 工具函数
│   │   ├── fitting.rs      # Nelder-Mead 优化、多指数拟合
│   │   ├── ufs.rs          # UFS 文件解析
│   │   └── csv_parser.rs   # CSV 解析
│   └── Cargo.toml
├── .claude/                # 权限配置
└── 朱海明..._cleaned.md     # TA 光谱综述（参考用）
```

## WASM 编译与部署

```bash
# 编译 WASM（需要 wasm32 target）
cd core
rustup target add wasm32-unknown-unknown
cargo build --release --target wasm32-unknown-unknown

# 生成胶水文件（如需部署到 GitHub Pages 启用 WASM 加速）
wasm-pack build --target web --release
# 将 pkg/ 目录部署到 GitHub Pages 根目录
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
UFS/CSV 导入 → 探针波长裁剪 → 基线扣除 → 啁啾校正 → SNR 过滤 → 可选手动选点 → 结果可视化 + 下载
```

## 支持的文件格式

- **UFS**：飞秒 TA 实验专有格式，多文件载入（参考、信号、延迟时间自动匹配）
- **CSV**：通用格式，包含 time/wavelength/ΔA 三列数据
- 标准数据文件：`paper/bibliography/pdf/` 下的 PDF 文献通过 `batch_convert.bat` 批量转为清洗后 Markdown

## 2026-06-09 变更

- **啁啾拟合公式修复**：从 `t₀ = c₀ + c₁·λ + c₂·λ²`（无物理意义）改为 `t₀ = c₀ + c₁·(1/λ²) + c₂·(1/λ⁴)`（材料色散模型）
- Rust（`chirp.rs`）和 JS fallback（`app.js`）同步更新
- WASM 已重新编译（待部署 `pkg/` 后生效）
