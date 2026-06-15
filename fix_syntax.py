import numpy as np
import pandas as pd

df = pd.read_csv('/data2/zhangzihan/sample_0_reads/sample_k/kmer_results_metrics.csv')
df = df.sort_values('k')
k_vals = df['k'].values
soft_eff = df['soft_effective_diff'].values

# 计算融合信号（与之前相同）
abs_diff = np.abs(np.diff(soft_eff))
rel_diff = abs_diff / np.maximum(soft_eff[:-1], 1e-10)
def curvature(y):
    x = np.arange(len(y))
    dy = np.gradient(y, x)
    ddy = np.gradient(dy, x)
    return np.abs(ddy) / (1 + dy**2)**1.5
curv = curvature(soft_eff)[:-1]
def normalize(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x) + 1e-10)
fused_raw = 0.5*normalize(abs_diff) + 0.4*normalize(curv) + 0.1*normalize(rel_diff)

# 计算一阶差分（即 fused_raw 的差异）
diff_fused = np.diff(fused_raw)
max_diff_idx = np.argmax(diff_fused)  # 差分最大的位置
print(f"原始融合信号差分最大的索引: {max_diff_idx}, 对应k: {k_vals[max_diff_idx]} ~ {k_vals[max_diff_idx+1]}")