# analyze_single_k.py
import argparse
import math
import os

def factorial(n):
    """计算阶乘"""
    if n == 0:
        return 1
    return math.factorial(n)

def find_error_cutoff(F):
    """泊松模型计算错误k-mer阈值"""
    lambda_val = 1.0
    max_check = 10
    min_cons = 2
    
    # 计算F[1]到F[max_check]的总和
    total = sum(F[1:max_check+1]) if len(F) > max_check else sum(F[1:])
    
    if total == 0:
        return 2
    
    # 计算泊松预测值
    pred = []
    for f in range(1, max_check + 1):
        # Poisson(f) = total * exp(-lambda) * lambda^f / f!
        poisson_val = total * math.exp(-lambda_val) * (lambda_val ** f) / factorial(f)
        pred.append(poisson_val)
    
    # 寻找第一个连续min_cons个频率都大于预测值的点
    for f in range(1, max_check - min_cons + 2):
        valid = True
        for j in range(min_cons):
            freq_index = f + j
            if freq_index >= len(F) or F[freq_index] <= pred[freq_index - 1]:
                valid = False
                break
        if valid:
            return f
    
    return 2

def calculate_effective_diff(F):
    """计算有效差异值 F[0] - sum(F[1...error_cutoff])"""
    cutoff = find_error_cutoff(F)
    
    # F[0] 是 distinct kmers
    distinct_kmers = F[0]
    
    # 计算错误k-mer的总数 (F[1]到F[cutoff-1])
    error_sum = sum(F[1:cutoff])
    
    effective_diff = distinct_kmers - error_sum
    if effective_diff < 0:
        effective_diff = 0
    
    return effective_diff, cutoff, error_sum

def parse_spectrum_file(file_path):
    """解析KmerLight输出文件 - 根据实际格式调整"""
    F = [0] * 11  # F[0]到F[10]，因为nfreq=10
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # 解析 F0
                if line.startswith('#F0 = '):
                    F[0] = int(line.split('= ')[1])
                
                # 解析 f1 到 f10
                elif line.startswith('f') and not line.startswith('f\t'):
                    parts = line.split()
                    if len(parts) >= 2:
                        # 提取频率编号，如 'f1' -> 1
                        freq_num = parts[0][1:]  # 去掉'f'
                        if freq_num.isdigit():
                            freq_idx = int(freq_num)
                            if 1 <= freq_idx <= 10:  # 根据nfreq调整
                                F[freq_idx] = int(parts[1])
        
        print(f"Parsed F array: {F}")
        
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
    
    return F

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spectrum-file', required=True)
    parser.add_argument('--k-value', type=int, required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    print(f"Analyzing k={args.k_value} from {args.spectrum_file}")
    
    # 解析频谱文件
    F = parse_spectrum_file(args.spectrum_file)
    
    # 计算有效差异值
    effective_diff, cutoff, error_sum = calculate_effective_diff(F)
    
    # 输出结果
    with open(args.output, 'w') as f:
        f.write(f"k_value\t{args.k_value}\n")
        f.write(f"effective_diff\t{effective_diff}\n")
        f.write(f"error_cutoff\t{cutoff}\n")
        f.write(f"error_sum\t{error_sum}\n")
        f.write(f"distinct_kmers\t{F[0]}\n")
    
    print(f"k={args.k_value}: F[0]={F[0]}, error_cutoff={cutoff}, error_sum={error_sum}, effective_diff={effective_diff}")

if __name__ == "__main__":
    main()