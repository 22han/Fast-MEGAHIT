import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# 1. 生成理想化双峰数据 (带错误尾)
x = np.arange(1, 550, 1)
# 错误尾：指数衰减
error_tail = 2e7 * np.exp(-1.5 * (x - 1)) 
error_tail[x < 1] = 0
# 稀有种峰：高斯分布，均值 25
rare_peak = 6e6 * np.exp(-0.5 * ((x - 25) / 6)**2)
# 优势种峰：高斯分布，均值 320
dominant_peak = 8e7 * np.exp(-0.5 * ((x - 320) / 60)**2)
# 叠加 + 少量噪声
y = error_tail + rare_peak + dominant_peak + np.random.normal(0, 1e4, len(x))
y[y < 0] = 0

# 2. 绘图设置
plt.style.use('default') # 使用干净的基础风格
fig, ax = plt.subplots(figsize=(9, 5))

# 3. 分区填充 (学术半透明色)
# 错误区 (灰色)
ax.fill_between(x, y, where=(x <= 8), color='#95A5A6', alpha=0.5, label='Sequencing Errors (f≤8)')
# 稀有种区 (橙色)
ax.fill_between(x, y, where=(x > 8) & (x <= 100), color='#E67E22', alpha=0.4, label='Rare Strains')
# 优势种区 (蓝色)
ax.fill_between(x, y, where=(x > 100), color='#3498DB', alpha=0.3, label='Dominant Strains')

# 绘制主轮廓线
ax.plot(x, y, color='#2C3E50', linewidth=2, zorder=10)

# 4. 核心标注 (学术级)
# 标注 1：双峰异质性 (双向箭头)
# 找到稀有种峰顶和优势种峰顶的 x 坐标
rare_x = 25
dom_x = 320
rare_y_val = rare_peak[24]
dom_y_val = dominant_peak[319]

# 在峰顶画箭头
ax.annotate(r'Rare Peak ($\mu \approx 25$)', xy=(rare_x, rare_y_val), xytext=(rare_x - 30, rare_y_val * 1.2),
            arrowprops=dict(arrowstyle='->', color='#E67E22', lw=1.5),
            ha='center', fontsize=11, fontweight='bold', color='#D35400')

ax.annotate(r'Dominant Peak ($\mu \approx 320$)', xy=(dom_x, dom_y_val), xytext=(dom_x + 40, dom_y_val * 1.1),
            arrowprops=dict(arrowstyle='->', color='#2980B9', lw=1.5),
            ha='center', fontsize=11, fontweight='bold', color='#2980B9')

# 在底部画双向箭头表示覆盖度差异
ax.annotate('', xy=(rare_x, rare_y_val * 0.3), xytext=(dom_x, rare_y_val * 0.3),
            arrowprops=dict(arrowstyle='<->', color='#8E44AD', lw=2.0))
ax.text((rare_x + dom_x)/2, rare_y_val * 0.35, 'Coverage Heterogeneity: $10^3-10^4$ fold', 
        ha='center', fontsize=11, color='#8E44AD', fontweight='bold', 
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8))

# 5. 坐标轴美化 (去边框)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('k-mer Frequency', fontsize=12, fontweight='bold')
ax.set_ylabel('k-mer Count', fontsize=12, fontweight='bold')
ax.set_xlim(0, 500)
ax.set_ylim(bottom=1000)
ax.set_yscale('log') # 对数坐标更能看清差异
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('Academic_Bimodal_Spectrum.png', dpi=300, bbox_inches='tight')
print("✅ 已生成学术级示意图: Academic_Bimodal_Spectrum.png")