#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
直接测试 min_size=4 的 k 值选择结果
修改 kmer_batch.py 中的默认 min_size 参数
"""

import os
import sys

# 修改 kmer_batch.py 中的默认 min_size
print("修改 kmer_batch.py 中的默认 min_size 参数...")

kmer_batch_file = '/home/zhangzihan/Fast_MEGAHIT/kmer_batch.py'
with open(kmer_batch_file, 'r') as f:
    content = f.read()

# 查找并替换 fit_two_stage_em 中的 min_size=5 为 min_size=4
old_line = "self.construct_adaptive_partitions(min_size=5)"
new_line = "self.construct_adaptive_partitions(min_size=4)  # [TEST] changed from 5 to 4"

if old_line in content:
    content = content.replace(old_line, new_line)
    with open(kmer_batch_file, 'w') as f:
        f.write(content)
    print(f"✓ 已将 min_size 从 5 改为 4")
else:
    print(f"✗ 未找到目标行：{old_line}")
    sys.exit(1)

print("\n现在可以运行 run_pipeline.py 测试 min_size=4 的效果")
print("命令示例:")
print("  python run_pipeline.py -i /data2/zhangzihan/sample_16_reads/sample_kmer -o /data2/zhangzihan/sample_16_reads/sample_k_min4")
print("\n或者运行 experiment4_critic_comparison.py:")
print("  python experiment4_critic_comparison.py --spectrum-dir /data2/zhangzihan/sample_16_reads/sample_kmer --output-dir experiment4_output_min4")
