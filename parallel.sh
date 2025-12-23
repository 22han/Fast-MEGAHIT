#!/bin/bash
# run_multi_k_parallel.sh

# 记录开始时间
START_TIME=$(date +%s)
echo "Pipeline started at: $(date)"

# 使用方法: ./run_multi_k_parallel.sh <kmerlight_dir> <scripts_dir> <r1_file> <r2_rc_file>
if [ $# -ne 4 ]; then
    echo "Usage: $0 <kmerlight_dir> <scripts_dir> <r1_file> <r2_rc_file>"
    echo "Example: $0 /home/user/kmerlight-master /home/user/scripts /data/SRR33132322_1.fastq /data/SRR33132322_R2_rc.fastq"
    exit 1
fi

KMERLIGHT_DIR="$1"
SCRIPTS_DIR="$2"
R1_FILE="$3"
R2_RC_FILE="$4"    # 用于kmerlight的反向互补文件

# 检查文件是否存在
if [ ! -f "$KMERLIGHT_DIR/kmerlight" ]; then
    echo "Error: kmerlight not found at $KMERLIGHT_DIR/kmerlight"
    exit 1
fi

if [ ! -f "$SCRIPTS_DIR/analyze_single_k.py" ]; then
    echo "Error: analyze_single_k.py not found at $SCRIPTS_DIR/analyze_single_k.py"
    exit 1
fi

if [ ! -f "$SCRIPTS_DIR/select_best_for_interval.py" ]; then
    echo "Error: select_best_for_interval.py not found at $SCRIPTS_DIR/select_best_for_interval.py"
    exit 1
fi

if [ ! -f "$R1_FILE" ]; then
    echo "Error: R1 file not found: $R1_FILE"
    exit 1
fi

if [ ! -f "$R2_RC_FILE" ]; then
    echo "Error: R2 RC file not found: $R2_RC_FILE"
    exit 1
fi

echo "KmerLight directory: $KMERLIGHT_DIR"
echo "Scripts directory: $SCRIPTS_DIR"
echo "Input files: $R1_FILE, $R2_RC_FILE"

# 创建输出目录
mkdir -p spectra analysis_results interval_results

# 定义8个区间
intervals=("21-36" "37-52" "53-68" "69-84" "85-100" "101-116" "117-132" "133-148")

echo "Starting parallel k-mer spectrum computation and analysis by intervals..."

# 为每个区间创建并行任务
for interval in "${intervals[@]}"; do
    (
        start=${interval%-*}
        end=${interval#*-}
        
        echo "Processing interval [$start, $end] in parallel..."
        
        mkdir -p analysis_results/interval_${start}_${end}
        
        # 使用函数方式避免特殊字符问题
        process_k() {
            k=$1
            echo "Processing k=$k in interval [$start, $end]..."
            
            # 使用绝对路径调用kmerlight
            "$KMERLIGHT_DIR/kmerlight" -k "$k" -f 10 -i "$R1_FILE" "$R2_RC_FILE" -o "spectra/k${k}.txt"
            
            /home/zhangzihan/miniconda3/bin/python3 "$SCRIPTS_DIR/analyze_single_k.py" --spectrum-file "spectra/k${k}.txt" --k-value "$k" --output "analysis_results/interval_${start}_${end}/k${k}_result.txt"
        }
        
        # 导出函数以便parallel使用
        export -f process_k
        export KMERLIGHT_DIR R1_FILE R2_RC_FILE SCRIPTS_DIR start end
        
        # 使用函数调用parallel
        parallel -j 2 process_k ::: $(seq $start 2 $end)
        
        echo "Selecting best k for interval [$start, $end]..."
        
        /home/zhangzihan/miniconda3/bin/python3 "$SCRIPTS_DIR/select_best_for_interval.py" \
            --results-dir "analysis_results/interval_${start}_${end}" \
            --interval "$start-$end" \
            --output "interval_results/interval_${start}_${end}_best.txt"
        
        if [ -f "interval_results/interval_${start}_${end}_best.txt" ]; then
            best_k=$(cat "interval_results/interval_${start}_${end}_best.txt")
            echo "Interval [$start, $end] completed. Best k = $best_k"
        else
            echo "Interval [$start, $end] completed. No best k found."
        fi
        
    ) &
done

wait

echo "All intervals completed. Combining results..."

# 检查是否有结果文件
if ls interval_results/interval_*_best.txt >/dev/null 2>&1; then
    cat interval_results/interval_*_best.txt > selected_k.txt
    k_list=$(cat selected_k.txt | tr '\n' ',')
    echo "================================================"
    echo "OPTIMAL K-VALUES SELECTED: ${k_list%,}"
    echo "Results saved to: selected_k.txt"
    echo "================================================"
    
    # 显示选择的k值详情
    echo "Detailed results per interval:"
    for interval in "${intervals[@]}"; do
        start=${interval%-*}
        end=${interval#*-}
        if [ -f "interval_results/interval_${start}_${end}_best.txt" ]; then
            best_k=$(cat "interval_results/interval_${start}_${end}_best.txt")
            echo "  [$start-$end]: k = $best_k"
        fi
    done
else
    echo "Error: No interval results found!"
    exit 1
fi

# 计算总执行时间
END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))

# 格式化输出总时间
echo "================================================"
echo "Pipeline completed at: $(date)"
echo "Total k-value selection time: $(($TOTAL_DURATION / 3600)) hours, $((($TOTAL_DURATION % 3600) / 60)) minutes and $(($TOTAL_DURATION % 60)) seconds"
echo "================================================"
echo "Next step: Use the selected k-values with MEGAHIT:"
echo "  megahit -1 $R1_FILE -2 <your_original_R2_file> --k-list \"${k_list%,}\" -o megahit_output"
echo "================================================" 