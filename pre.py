#!/usr/bin/env python3
import sys
from pathlib import Path

def process_fastq_pair(r1_path, r2_path):
    """处理FASTQ对：过滤含N的reads，并生成3个输出文件"""
    r1_path, r2_path = Path(r1_path), Path(r2_path)
    
    out1_path = r1_path.with_name('clean_1.fastq')
    out2_path = r2_path.with_name('clean_2.fastq')
    out2_rc_path = r2_path.with_name('clean_2_rc.fastq')

    # 创建转换表
    trans_table = str.maketrans('ACGTacgt', 'TGCAtgca')
    
    record_count = 0
    with open(r1_path, 'r') as f1, \
         open(r2_path, 'r') as f2, \
         open(out1_path, 'w') as o1, \
         open(out2_path, 'w') as o2, \
         open(out2_rc_path, 'w') as o2_rc:

        while True:
            # 读取R1
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline().rstrip()
            p1 = f1.readline()
            q1 = f1.readline().rstrip()
            
            # 读取R2
            h2 = f2.readline()
            s2 = f2.readline().rstrip()
            p2 = f2.readline()
            q2 = f2.readline().rstrip()
            
            # 检查是否包含N
            if 'N' in s1 or 'n' in s1 or 'N' in s2 or 'n' in s2:
                continue
            
            # 对R2进行反向互补（用于第三个文件）
            rev_comp_s2 = s2.translate(trans_table)[::-1]
            rev_q2 = q2[::-1]
            
            # 写入结果
            # 文件1: clean_1.fastq (原始R1)
            o1.write(f"{h1}{s1}\n{p1}{q1}\n")
            # 文件2: clean_2.fastq (原始R2)
            o2.write(f"{h2}{s2}\n{p2}{q2}\n")
            # 文件3: clean_2_rc.fastq (反向互补的R2)
            o2_rc.write(f"{h2}{rev_comp_s2}\n{p2}{rev_q2}\n")
            
            record_count += 1

    print(f"处理完成！共保留 {record_count} 条reads")
    print(f"输出文件：")
    print(f"  1. {out1_path}      (clean R1)")
    print(f"  2. {out2_path}      (clean R2)")
    print(f"  3. {out2_rc_path}   (clean R2反向互补)")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python3 pre.py SFR_1.fastq SFR_2.fastq')
        sys.exit(1)
    process_fastq_pair(sys.argv[1], sys.argv[2])