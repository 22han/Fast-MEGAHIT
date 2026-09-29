#!/usr/bin/env python3
"""
GitHub 友好版：生成极小但特征完美的双端 FASTQ 测试数据。
数据量控制在几十 KB，方便上传到 GitHub，同时保留泊松波动以确保 KMC 频谱有双峰。
"""
import random
import numpy as np

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(comp.get(base, 'N') for base in reversed(seq))

def generate_tiny_fastq(r1_path, r2_rc_path, seed=42):
    random.seed(seed)
    np.random.seed(seed)
    bases = 'ACGT'
    read_len = 150
    
    reads_r1 = []
    reads_r2_rc = []
    read_id_counter = 0
    
    # 【核心修改】：大幅缩减数据量
    # 主物种：1条 1000bp 序列，50x 覆盖度 -> 约 300 条 reads
    # 稀有种：1条 1000bp 序列，15x 覆盖度 -> 约 90 条 reads
    # 总数据量约 400 条 reads，文件大小 < 200 KB
    species_config = [
        (1, 1000, 50, "Dominant"),
        (1, 1000, 15, "Rare")
    ]
    
    for n_contigs, contig_len, target_cov, label in species_config:
        for c_idx in range(n_contigs):
            # 1. 生成随机 Contig
            contig = ''.join(random.choices(bases, k=contig_len))
            
            # 2. 模拟真实测序：将 Contig 分块，每块的覆盖度服从泊松分布
            num_fragments = contig_len // read_len
            for i in range(num_fragments):
                start = i * read_len
                # 使用泊松分布产生覆盖度波动
                n_reads_for_this_block = np.random.poisson(target_cov)
                
                for _ in range(n_reads_for_this_block):
                    # 在 block 内随机起始
                    local_start = random.randint(0, max(0, read_len - 20))
                    actual_start = start + local_start
                    if actual_start + read_len > contig_len:
                        continue
                        
                    r1_seq = contig[actual_start : actual_start + read_len]
                    r2_raw_seq = contig[actual_start : actual_start + read_len]
                    r2_rc_seq = reverse_complement(r2_raw_seq)
                    
                    # 3. 模拟 1% 的测序错误（产生低频错误峰）
                    def add_errors(seq):
                        seq_list = list(seq)
                        for j in range(len(seq_list)):
                            if random.random() < 0.01:
                                seq_list[j] = random.choice(bases)
                        return ''.join(seq_list)
                        
                    r1_seq = add_errors(r1_seq)
                    r2_rc_seq = add_errors(r2_rc_seq)
                    
                    qual = 'I' * read_len
                    rid = f"read_{read_id_counter}_{label}_{c_idx}"
                    read_id_counter += 1
                    
                    reads_r1.append(f"@{rid}/1\n{r1_seq}\n+\n{qual}\n")
                    reads_r2_rc.append(f"@{rid}/2\n{r2_rc_seq}\n+\n{qual}\n")
                    
    # 写入文件
    with open(r1_path, 'w') as f1, open(r2_rc_path, 'w') as f2:
        f1.writelines(reads_r1)
        f2.writelines(reads_r2_rc)
        
    print(f"✅ 极小 FASTQ 生成完毕 (GitHub 友好版):")
    print(f"   - {r1_path} ({len(reads_r1)} reads)")
    print(f"   - {r2_rc_path} ({len(reads_r2_rc)} reads)")
    print(f"   - 预期 KMC 频谱特征：")
    print(f"     * 1x~10x: 错误峰")
    print(f"     * ~15x:   稀有种峰")
    print(f"     * ~50x:   主物种峰")

if __name__ == '__main__':
    # 生成文件
    generate_tiny_fastq("mock_1.fastq", "mock_2.fastq")