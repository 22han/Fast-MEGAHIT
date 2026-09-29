import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# 1. 配色方案
C_STAGE1 = '#2E86AB'  # 深蓝
C_STAGE2 = '#A23B72'  # 紫红
C_STAGE3 = '#F18F01'  # 橙色
C_STAGE4 = '#C73E1D'  # 砖红
C_BG = '#F5F5F5'
C_TEXT = '#333333'

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10})

fig = plt.figure(figsize=(12, 10), facecolor=C_BG)

# ==========================================
# 顶部横幅 (占全宽，顶部 10%)
# ==========================================
ax_top = fig.add_axes([0.05, 0.90, 0.90, 0.08])
ax_top.axis('off')
ax_top.text(0.5, 0.7, 'AdaptiveKmer: A Hierarchical Optimization Framework', 
            ha='center', va='center', fontsize=16, weight='bold', color=C_TEXT)
ax_top.text(0.5, 0.3, r'$\mathcal{K}^* = \mathcal{G}(\mathcal{R}(\mathcal{P}(\mathcal{E}(\mathcal{K}))))$  '
            r'(E=PESC, P=Partition, R=Ranking, G=Generation)', 
            ha='center', va='center', fontsize=12, color=C_TEXT)

# ==========================================
# 左侧主区域：严格占 2/3 宽度 (0.05 到 0.68)
# ==========================================
ax_left = fig.add_axes([0.05, 0.05, 0.63, 0.82])
ax_left.set_xlim(0, 1)
ax_left.set_ylim(0, 1)
ax_left.axis('off')
ax_left.set_title('Methodological Framework (Stages I-IV)', fontsize=12, pad=10)

# 定义四个 Stage 的垂直位置和高度
stage_y = [0.75, 0.52, 0.29, 0.06]
stage_h = 0.18
stage_colors = [C_STAGE1, C_STAGE2, C_STAGE3, C_STAGE4]
stage_titles = ['Stage I: PESC', 'Stage II: Partition', 'Stage III: CRITIC', 'Stage IV: Schedule']
stage_formulas = [
    r'$P(f|k) = w_{err}G + w_{sig}[(1-w_{rare})N_1 + w_{rare}N_2]$' + '\n' + r'$\rightarrow H_{soft}(k)$',
    r'$dH/dk \rightarrow$ Top-down segmentation' + '\n' + r'$P_1, P_2, P_3$ (MSE minimization)',
    r'$O_{ret}, O_{err}, O_{res}, O_{sharp}$' + '\n' + r'$\downarrow$ CRITIC weights $\rightarrow k^*$',
    r'$|k_{i+1} - k_i| \le 28$' + '\n' + 'Gap repair + Insertion'
]
stage_challenges = [
    'Challenge 1: Stabilize local evaluation',
    'Challenge 2: Decompose heterogeneity',
    'Objective weighting without manual tuning',
    'Assembler-compatible schedule'
]

# 绘制四个 Stage 的框和内容
for i in range(4):
    y = stage_y[i]
    color = stage_colors[i]
    
    # 绘制圆角矩形框
    box = FancyBboxPatch((0.05, y), 0.9, stage_h, 
                         boxstyle="round,pad=0.02", 
                         linewidth=2, edgecolor=color, facecolor='white', alpha=0.9)
    ax_left.add_patch(box)
    
    # 内部色块 (左侧细条)
    inner_box = FancyBboxPatch((0.05, y), 0.05, stage_h, 
                               boxstyle="round,pad=0.02", 
                               linewidth=0, facecolor=color, alpha=0.8)
    ax_left.add_patch(inner_box)
    
    # 标题
    ax_left.text(0.12, y + stage_h - 0.03, stage_titles[i], fontsize=11, 
                 weight='bold', color=color, va='top')
    
    # 公式/内容
    ax_left.text(0.5, y + stage_h/2 - 0.02, stage_formulas[i], fontsize=10, 
                 ha='center', va='center', color=C_TEXT)
    
    # 底部挑战标注
    ax_left.text(0.5, y - 0.02, stage_challenges[i], fontsize=8, 
                 ha='center', va='top', color='gray', style='italic')

# 绘制连接箭头
for i in range(3):
    arrow = FancyArrowPatch((0.5, stage_y[i]), (0.5, stage_y[i] + stage_h - 0.005),
                            arrowstyle='->', mutation_scale=20, color=C_TEXT, linewidth=1.5)
    ax_left.add_patch(arrow)

# ==========================================
# 右侧区域：严格占 1/3 宽度 (0.70 到 0.95)
# ==========================================
ax_right = fig.add_axes([0.70, 0.05, 0.25, 0.82])
ax_right.set_xlim(0, 1)
ax_right.set_ylim(0, 1)
ax_right.axis('off')
ax_right.set_title('Final Output', fontsize=12, pad=10)

# 绘制水平坐标轴
ax_right.axhline(0.6, color=C_TEXT, linewidth=2, xmin=0.1, xmax=0.9)
ax_right.text(0.9, 0.6, 'k', fontsize=12, ha='left', va='center')

# 非均匀分布的圆点
k_points = [0.15, 0.25, 0.45, 0.6, 0.85]
ax_right.scatter(k_points, [0.6]*len(k_points), color=C_STAGE3, s=80, zorder=5, edgecolors='white')

# 绘制虚线弧
for i in range(len(k_points)-1):
    start = k_points[i]
    end = k_points[i+1]
    center = (start + end) / 2
    radius = (end - start) / 2
    arc = patches.Arc((center, 0.6), width=radius*2, height=radius*2, angle=0, 
                      theta1=0, theta2=180, color=C_STAGE4, linestyle='--', linewidth=1.5)
    ax_right.add_patch(arc)
    ax_right.text(center, 0.6 + radius + 0.02, 'Gap ≤ 28', fontsize=7, ha='center', color=C_STAGE4)

# 底部图标与说明
ax_right.text(0.5, 0.4, '⚙ MEGAHIT\n(Assembler)', fontsize=10, ha='center', va='center', color=C_TEXT)
ax_right.text(0.5, 0.25, 'Compact, Non-uniform,\nDataset-specific', fontsize=9, ha='center', va='center', 
              color=C_TEXT, style='italic')

plt.savefig('AdaptiveKmer_Architecture.png', dpi=300, bbox_inches='tight', facecolor=C_BG)
plt.savefig('AdaptiveKmer_Architecture.pdf', format='pdf', bbox_inches='tight', facecolor=C_BG)
print("生成成功：AdaptiveKmer_Architecture.png / .pdf")
