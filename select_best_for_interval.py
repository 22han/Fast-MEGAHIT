# select_best_for_interval.py
import glob
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--interval', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    print(f"Selecting best k for interval {args.interval}")
    
    interval_results = {}
    
    # 读取该区间内所有k值的分析结果
    result_files = glob.glob(f"{args.results_dir}/k*_result.txt")
    
    for result_file in result_files:
        try:
            # 从文件名提取k值
            filename = os.path.basename(result_file)
            k_value = int(filename.split('_')[0][1:])  # 从"k21_result.txt"提取21
            
            # 读取分析结果
            with open(result_file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('effective_diff'):
                        effective_diff = float(line.split('\t')[1])
                        interval_results[k_value] = effective_diff
                        break
                        
            print(f"  k={k_value}: effective_diff={effective_diff}")
            
        except Exception as e:
            print(f"Error processing {result_file}: {e}")
    
    # 选择区间内最佳k值
    if interval_results:
        best_k = max(interval_results.keys(), key=lambda k: interval_results[k])
        best_diff = interval_results[best_k]
        print(f"  Best k={best_k}, effective_diff={best_diff}")
    else:
        # 如果没有结果，使用区间起始值
        start = int(args.interval.split('-')[0])
        best_k = start
        print(f"  No results, using default k={best_k}")
    
    # 输出该区间的最佳k值
    with open(args.output, 'w') as f:
        f.write(f"{best_k}\n")
    
    print(f"Interval {args.interval} best k saved to {args.output}")

if __name__ == "__main__":
    main()