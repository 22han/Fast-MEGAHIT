# 绘图脚本使用说明

本项目包含三个专门的绘图脚本，用于生成论文中的图表：

## 目录

1. [图3.2：双峰拟合效果对比](#图32双峰拟合效果对比)
2. [图3.3：硬截断 vs PESC软计数对比](#图33硬截断-vs-pesc软计数对比)
3. [图3.4：边界熵热力图](#图34边界熵热力图)

---

## 图3.2：双峰拟合效果对比

### 文件
`plot_bimodal_comparison.py`

### 功能
展示不同 k 值下频谱的单峰 vs 双峰拟合效果对比

### 使用方法

```bash
# 基本用法
python3 plot_bimodal_comparison.py --spectrum-dir <频谱目录>

# 指定 k 值
python3 plot_bimodal_comparison.py --spectrum-dir ./data --k-values "21,31,41,51"

# 同时生成详细图
python3 plot_bimodal_comparison.py --spectrum-dir ./data --individual
```

### 参数说明
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--spectrum-dir` | 频谱文件所在目录（必填） | - |
| `--k-values` | 要分析的 k 值，逗号分隔 | 21,31,41,51 |
| `--output-dir` | 输出目录 | plots |
| `--individual` | 同时绘制每个 k 的详细图 | 不启用 |

### 输出文件
- `图3.2_双峰拟合效果对比.png` - 主对比图
- `k{K}_双峰拟合详细图.png` - 各 k 值的详细图（可选）

---

## 图3.3：硬截断 vs PESC软计数对比

### 文件
`plot_hard_vs_soft.py`

### 功能
展示有效 k-mer 数随 k 值变化的曲线，对比硬截断法和 PESC 软计数法的差异

### 使用方法

```bash
# 基本用法
python3 plot_hard_vs_soft.py --spectrum-dir <频谱目录>

# 指定 k 值范围
python3 plot_hard_vs_soft.py --spectrum-dir ./data --k-start 21 --k-end 121 --k-step 10
```

### 参数说明
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--spectrum-dir` | 频谱文件所在目录（必填） | - |
| `--k-start` | 起始 k 值 | 21 |
| `--k-end` | 结束 k 值 | 121 |
| `--k-step` | k 值步长 | 10 |
| `--output-dir` | 输出目录 | plots |

### 输出文件
- `图3.3_硬截断vs软计数对比.png` - 主对比图
- `硬截断vs软计数联合对比.png` - 联合对比图
- `差分对比图.png` - 差分变化图

---

## 图3.4：边界熵热力图

### 文件
`plot_boundary_entropy.py`

### 功能
展示不同 k 值下边界附近频率点的后验概率分布，颜色代表熵值

### 使用方法

```bash
# 基本用法
python3 plot_boundary_entropy.py --spectrum-dir <频谱目录>

# 指定参数
python3 plot_boundary_entropy.py --spectrum-dir ./data --k-start 21 --k-end 81 --window 5 --individual
```

### 参数说明
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--spectrum-dir` | 频谱文件所在目录（必填） | - |
| `--k-start` | 起始 k 值 | 21 |
| `--k-end` | 结束 k 值 | 81 |
| `--k-step` | k 值步长 | 10 |
| `--window` | 边界窗口大小 | 5 |
| `--output-dir` | 输出目录 | plots |
| `--individual` | 同时绘制每个 k 的详细图 | 不启用 |

### 输出文件
- `图3.4_边界熵热力图.png` - 主热力图
- `k{K}_边界熵详细图.png` - 各 k 值的详细图（可选）
- `边界熵总结图.png` - 熵统计总结图

---

## 快速开始指南

### 1. 准备频谱文件
确保您的频谱文件位于指定目录，文件名格式支持：
- `spectrum_k{K}.txt`
- `k{K}_spectrum.txt`
- `ntcard_k{K}.txt`
- `k{K}.txt`

### 2. 运行所有绘图脚本

```bash
# 创建输出目录
mkdir -p plots

# 图3.2
python3 plot_bimodal_comparison.py --spectrum-dir ./ERR13033322_kmer --k-values "21,31,41,51" --individual

# 图3.3
python3 plot_hard_vs_soft.py --spectrum-dir ./ERR13033322_kmer --k-start 21 --k-end 121 --k-step 10

# 图3.4
python3 plot_boundary_entropy.py --spectrum-dir ./ERR13033322_kmer --k-start 21 --k-end 81 --k-step 10 --individual
```

### 3. 查看结果
所有图表将保存在 `plots` 目录中。

---

## 依赖要求

所有脚本都需要以下依赖：
- `numpy`
- `matplotlib`
- `scipy`
- `analyze_single_k.py`（或 `kmer_analyzer.py`）

这些依赖已包含在项目环境中。

---

## 注意事项

1. 确保频谱文件目录包含对应的 k 值文件
2. 如果某个 k 值的文件缺失，该 k 值会被自动跳过
3. 建议使用足够范围的 k 值以获得更好的对比效果
4. 生成的图片使用 300 DPI 高分辨率，适合论文发表

---

## 示例输出

```
============================================================
图3.2：不同 k 值下频谱的单峰 vs 双峰拟合效果对比
============================================================
频谱目录: ./ERR13033322_kmer
分析 k 值: [21, 31, 41, 51]
输出目录: plots

正在处理 k = 21...
正在处理 k = 31...
正在处理 k = 41...
正在处理 k = 51...

图片已保存到: plots/图3.2_双峰拟合效果对比.png
详细图已保存到: plots/k21_双峰拟合详细图.png
详细图已保存到: plots/k31_双峰拟合详细图.png
详细图已保存到: plots/k41_双峰拟合详细图.png
详细图已保存到: plots/k51_双峰拟合详细图.png
```

---

## 故障排除

### 问题：找不到频谱文件
**解决**：检查 `--spectrum-dir` 参数是否正确，确认文件名格式是否支持

### 问题：拟合失败
**解决**：某些 k 值可能拟合失败，脚本会自动跳过这些 k 值并继续

### 问题：内存不足
**解决**：减小分析的 k 值范围，或使用更大的步长

---

## 技术细节

### 图表配色方案

| 颜色代码 | 用途 |
|----------|------|
| `#3498db` | 真实频谱数据 |
| `#e74c3c` | 错误成分 / 硬截断法 |
| `#27ae60` | 优势种信号 / PESC软计数 |
| `#f39c12` | 稀有种信号 |
| `#8e44ad` | 总拟合曲线 |
| `#9b59b6` | 熵相关可视化 |

### 熵值颜色映射
- **深色**：低熵，边界清晰
- **浅色**：高熵，边界模糊
