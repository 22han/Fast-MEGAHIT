#!/bin/bash
# run_multi_k_parallel.sh (已修复输出目录参数化问题)

# 获取当前脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KMERLIGHT_DIR="${SCRIPT_DIR}/kmerlight-master"

# --- 修改点 1: 更新使用方法说明和参数检查 ---
# 现在接受 2 个或 3 个参数
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <raw_r1> <raw_r2> [output_dir]"
    echo "Example: $0 /data/R1.fastq /data/R2.fastq my_assembly_result"
    exit 1
fi

RAW_R1="$1"
RAW_R2="$2"
# --- 修改点 2: 灵活获取输出目录 ---
# 如果提供了第3个参数则使用它，否则默认为 megahit_out
OUTPUT_DIR="${3:-megahit_out}"

# 记录开始时间
START_TIME=$(date +%s)
echo "Pipeline started at: $(date)"

# 检查必需文件
if [ ! -f "$KMERLIGHT_DIR/kmerlight" ]; then
    echo "Error: kmerlight not found at $KMERLIGHT_DIR/kmerlight"
    exit 1
fi

for script in analyze_single_k.py select_best_for_interval.py pre.py; do
    if [ ! -f "$SCRIPT_DIR/$script" ]; then
        echo "Error: $script not found in $SCRIPT_DIR"
        exit 1
    fi
done

# -------------------- 预处理步骤 --------------------
echo "Step 1: Preprocessing..."
python3 "$SCRIPT_DIR/pre.py" "$RAW_R1" "$RAW_R2"
if [ $? -ne 0 ]; then echo "Preprocessing failed"; exit 1; fi

base=$(basename "$RAW_R1" | sed 's/_1\.fastq$//')
DATA_DIR=$(dirname "$RAW_R1")
R1_FILE="${DATA_DIR}/${base}clean_1.fastq"
R2_CLEAN="${DATA_DIR}/${base}clean_2.fastq"
R2_RC_FILE="${DATA_DIR}/${base}clean_2_rc.fastq"

# -------------------- k-mer 分析部分 --------------------
# -------------------- k-mer 分析部分 --------------------
# -------------------- k-mer 分析部分 --------------------
mkdir -p spectra analysis_results interval_results

# 获取当前脚本的绝对路径（锁死 py 文件位置）
ABS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
K_EXE="${ABS_DIR}/kmerlight-master/kmerlight"

intervals=("21-36" "37-52" "53-68" "69-84" "85-100" "101-116" "117-132" "133-148")

# 定义一个彻底独立的函数
do_single_k_analysis() {
    local k=$1
    local s=$2
    local e=$3
    local script_path=$4
    local kmerlight_exe=$5
    local r1=$6
    local r2rc=$7

    local out_dir="analysis_results/interval_${s}_${e}"
    
    # 1. 运行 kmerlight
    "$kmerlight_exe" -k "$k" -f 10 -i "$r1" "$r2rc" -o "spectra/k${k}.txt"
    
    # 2. 运行你的 analyze_single_k.py (使用绝对路径调用)
    python3 "${script_path}/analyze_single_k.py" \
        --spectrum-file "spectra/k${k}.txt" \
        --k-value "$k" \
        --output "${out_dir}/k${k}_result.txt"
}

export -f do_single_k_analysis

for interval in "${intervals[@]}"; do
    start=${interval%-*}
    end=${interval#*-}
    
    # 关键：在启动后台进程前，先在主进程创建好文件夹，防止 Python 写入失败
    mkdir -p "analysis_results/interval_${start}_${end}"

    (
        echo "Processing interval [$start, $end]..."
        
        # 将路径作为参数传给 parallel，这样子进程绝对不会弄丢路径
        parallel -j 2 do_single_k_analysis {} "$start" "$end" \
            "$ABS_DIR" "$K_EXE" "$R1_FILE" "$R2_RC_FILE" \
            ::: $(seq "$start" 2 "$end")
        
        # 运行选择最佳 k 的 py 文件
        python3 "${ABS_DIR}/select_best_for_interval.py" \
            --results-dir "analysis_results/interval_${start}_${end}" \
            --interval "$interval" \
            --output "interval_results/interval_${start}_${end}_best.txt"
    ) &
done

wait

# 合并 k-list
k_list=$(cat interval_results/interval_*_best.txt 2>/dev/null | tr '\n' ',' | sed 's/,$//')

# -------------------- MEGAHIT 组装步骤 --------------------
echo "Step 2: Running MEGAHIT assembly..."

if [ -z "$k_list" ]; then
    echo "Error: No k-values selected"; exit 1
fi

BUILD_DIR="./build"
mkdir -p "$BUILD_DIR"


ABS_OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")

cd "$BUILD_DIR" || exit 1

if [ ! -x "./megahit" ]; then
    echo "Error: ./megahit not found in $(pwd)"; exit 1
fi

# 检查目标目录是否存在（MEGAHIT不允许输出目录已存在）
if [ -d "$ABS_OUTPUT_DIR" ]; then
    echo "Warning: $ABS_OUTPUT_DIR already exists. Deleting it for MEGAHIT..."
    rm -rf "$ABS_OUTPUT_DIR"
fi

echo "Running: ./megahit -1 $R1_FILE -2 $R2_CLEAN --k-list $k_list -o $ABS_OUTPUT_DIR" 
./megahit -1 "$R1_FILE" -2 "$R2_CLEAN" --k-list "$k_list" -o "$ABS_OUTPUT_DIR" 

cd - > /dev/null

# 计算时间并结束
END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
echo "Pipeline completed. Results in: $OUTPUT_DIR"