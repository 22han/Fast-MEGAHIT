#!/bin/bash
# ============================================================
# 文件名: full_pipeline.sh
# 功能: 自动化处理双端测序数据（预处理 → KMC → pipeline_custom → MEGAHIT）
# 用法: bash full_pipeline.sh [工作目录]
# ============================================================

set -euo pipefail

# -------------------- 用户配置 --------------------
THREADS=16
MEMORY=32
FREQ=1200
START_K=21
END_K=141
STEP=2

# ---- pipeline 自定义参数（你可以在这里改）----
MIN_PARTITIONS=6
MAX_PARTITIONS=12
MIN_SIZE=4
SENSITIVITY=0.5

MEGAHIT="/home/zhangzihan/Fast_MEGAHIT/build/megahit"

ENV_KMC="kmc_env"
ENV_VS2="vs2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -------------------- Conda 初始化 --------------------
if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
else
    echo "错误：未找到 conda"
    exit 1
fi

# -------------------- 参数 --------------------
WORK_DIR="${1:-$PWD}"
cd "$WORK_DIR" || { echo "错误：无法进入 $WORK_DIR"; exit 1; }
echo "工作目录：$WORK_DIR"

# -------------------- 步骤1：预处理（不需要 conda 环境）--------------------
step1_preprocess() {
    echo "========== 步骤1：预处理 =========="
    for r1_raw in *_1.fastq; do
        [[ -f "$r1_raw" ]] || continue
        if [[ "$r1_raw" == *clean* ]]; then
            continue
        fi
        base="${r1_raw%_1.fastq}"
        r2_raw="${base}_2.fastq"
        if [[ ! -f "$r2_raw" ]]; then
            echo "警告：跳过 $base，缺少 $r2_raw"
            continue
        fi
        clean1="${base}clean_1.fastq"
        clean2_rc="${base}clean_2_rc.fastq"
        if [[ -f "$clean1" && -f "$clean2_rc" ]]; then
            echo "跳过 $base（clean 文件已存在）"
            continue
        fi
        echo "预处理 $base ..."
        python3 "${SCRIPT_DIR}/pre.py" "$r1_raw" "$r2_raw"
        if [[ $? -ne 0 ]]; then
            echo "错误：pre.py 失败，样本 $base"
            exit 1
        fi
    done
    echo "预处理完成。"
}

# -------------------- 步骤2：KMC（需要 kmc_env）--------------------
step2_kmc() {
    echo "========== 步骤2：KMC 多 k-mer 计数 =========="
    conda activate "$ENV_KMC"
    command -v kmc >/dev/null || { echo "错误：环境 $ENV_KMC 中缺少 kmc"; exit 1; }
    command -v kmc_tools >/dev/null || { echo "错误：环境 $ENV_KMC 中缺少 kmc_tools"; exit 1; }

    for r1_clean in *clean_1.fastq; do
        [[ -f "$r1_clean" ]] || continue
        base="${r1_clean%clean_1.fastq}"
        r2_clean="${base}clean_2_rc.fastq"
        if [[ ! -f "$r2_clean" ]]; then
            echo "警告：跳过 $base，缺少 $r2_clean"
            continue
        fi
        out_dir="${base}_kmer"
        final_check="${out_dir}/k${END_K}_spectrum.txt"
        if [[ -f "$final_check" ]]; then
            echo "跳过 $base（$final_check 已存在）"
            continue
        fi
        echo "处理样本 $base ..."
        mkdir -p "$out_dir"
        list_file="${out_dir}/fastq_files.lst"
        > "$list_file"
        echo "$r1_clean" >> "$list_file"
        echo "$r2_clean" >> "$list_file"
        tmp_base="${out_dir}/tmp"
        mkdir -p "$tmp_base"
        for (( k=START_K; k<=END_K; k+=STEP )); do
            echo "  -> k=$k"
            tmp_dir="${tmp_base}/k${k}"
            mkdir -p "$tmp_dir"
            db_prefix="${out_dir}/kmc_db_k${k}"
            kmc -k"$k" -ci1 -cs"$FREQ" -t"$THREADS" -m"$MEMORY" \
                "@${list_file}" "${db_prefix}" "$tmp_dir" || exit 1
            kmc_tools transform "${db_prefix}" histogram "${out_dir}/k${k}_spectrum.txt" || exit 1
            rm -f "${db_prefix}.kmc_pre" "${db_prefix}.kmc_suf"
        done
        rm -rf "$tmp_base" "$list_file"
        echo "样本 $base 的 KMC 部分完成"
    done
    conda deactivate
    echo "KMC 处理完成。"
}

# -------------------- 步骤3：run_pipeline_custom.py（需要 vs2）--------------------
step3_pipeline() {
    echo "========== 步骤3：选取最佳 k 值组合（自定义 pipeline）=========="
    conda activate "$ENV_VS2"

    if [[ ! -f "${SCRIPT_DIR}/run_pipeline_custom.py" ]]; then
        echo "错误：run_pipeline_custom.py 未找到，请放在 $SCRIPT_DIR"
        exit 1
    fi

    # 遍历所有 _kmer 目录
    for kmer_dir in *_kmer; do
        [[ -d "$kmer_dir" ]] || continue
        sample="${kmer_dir%_kmer}"
        out_k_dir="${sample}_k_custom"

        if [[ -d "$out_k_dir" ]]; then
            echo "跳过 $sample（$out_k_dir 已存在）"
            continue
        fi

        echo "运行 run_pipeline_custom.py 于 $sample ..."
        echo "  参数: --min_partitions $MIN_PARTITIONS --max_partitions $MAX_PARTITIONS --min_size $MIN_SIZE --sensitivity $SENSITIVITY"

        python3 "${SCRIPT_DIR}/run_pipeline_custom.py" \
            -i "$kmer_dir" \
            -o "$out_k_dir" \
            --min_partitions "$MIN_PARTITIONS" \
            --max_partitions "$MAX_PARTITIONS" \
            --min_size "$MIN_SIZE" \
            --sensitivity "$SENSITIVITY" || exit 1

        echo "$sample pipeline 完成 -> $out_k_dir"
    done

    conda deactivate
    echo "自定义 pipeline 完成。"
}

# -------------------- 步骤4：MEGAHIT 组装 --------------------
step4_megahit() {
    echo "========== 步骤4：MEGAHIT 双模式组装 =========="
    if [[ ! -x "$MEGAHIT" ]]; then
        echo "错误：MEGAHIT 不可执行：$MEGAHIT"
        exit 1
    fi
    for r1_raw in *_1.fastq; do
        [[ -f "$r1_raw" ]] || continue
        if [[ "$r1_raw" == *clean* ]]; then
            continue
        fi
        sample="${r1_raw%_1.fastq}"
        r2_raw="${sample}_2.fastq"
        if [[ ! -f "$r2_raw" ]]; then
            echo "跳过 $sample，缺少 $r2_raw"
            continue
        fi

        # 普通模式（原始数据） -> _ou1
        out_ou="${sample}_ou1"
        if [[ ! -d "$out_ou" ]]; then
            echo "普通 megahit: $sample"
            "$MEGAHIT" -1 "$r1_raw" -2 "$r2_raw" -o "$out_ou" -t "$THREADS" \
                --k-min "$START_K" --k-max "$END_K" --k-step "$STEP" || \
                echo "警告：普通 megahit 失败"
        else
            echo "跳过普通组装 $out_ou（已存在）"
        fi

        # 自定义 k‑list 模式（使用 _k_custom 的结果） -> _oli_custom
        out_o="${sample}_oli_custom"
        clean1="${sample}clean_1.fastq"
        clean2="${sample}clean_2.fastq"
        klist_file="${sample}_k_custom/kmer_results_custom_selected_ks.txt"
        if [[ -d "$out_o" ]]; then
            echo "跳过自定义组装 $out_o（已存在）"
        elif [[ ! -f "$clean1" || ! -f "$clean2" ]]; then
            echo "跳过自定义组装 $sample：缺少 clean 文件"
        elif [[ ! -f "$klist_file" ]]; then
            echo "跳过自定义组装 $sample：缺少 $klist_file"
        else
            klist=$(grep '^K_LIST=\[' "$klist_file" | sed 's/K_LIST=\[\(.*\)\]/\1/' | tr -d ' ')
            if [[ -z "$klist" ]]; then
                echo "警告：$klist_file 中未找到 K_LIST，跳过自定义组装"
            else
                echo "自定义 k‑list megahit: $sample (klist = $klist)"
                "$MEGAHIT" -1 "$clean1" -2 "$clean2" -o "$out_o" -t "$THREADS" \
                    --k-list "$klist" || echo "警告：自定义 megahit 失败"
            fi
        fi
    done
    echo "组装完成。"
}

# -------------------- 主流程 --------------------
START_TIME=$(date +%s)
step1_preprocess
step2_kmc
step3_pipeline
step4_megahit
END_TIME=$(date +%s)
echo "================================================"
echo "全流程运行完毕，总耗时：$((END_TIME - START_TIME)) 秒"
echo "================================================"