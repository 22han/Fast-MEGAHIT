#!/bin/bash

BASE_DIR="/data2/zhangzihan/sample_19_reads"

# 遍历所有以 _kmer 结尾的目录
for dir in "${BASE_DIR}"/*_kmer; do
    # 确保是目录
    if [ -d "$dir" ]; then
        # 提取样本名（去掉末尾的 _kmer）
        SAMPLE=$(basename "$dir" | sed 's/_kmer$//')
        
        # 🆕 1. 修改输出目录名，加上 _ablation 后缀
        # 防止覆盖你之前跑的正常模式结果！
        OUTPUT_DIR="${BASE_DIR}/${SAMPLE}_k_ablation"

        # 如果输出目录不存在，则运行 pipeline
        if [ ! -d "$OUTPUT_DIR" ]; then
            echo "🧪 Running ABLATION pipeline for ${SAMPLE}..."
            
            # 🆕 2. 调用 run_pipeline_nopartition.py，并传入消融参数
            python run_pipeline_nopartition.py \
                -i "$dir" \
                -o "$OUTPUT_DIR" \
                --ablation \
                --n_top_k 12 \
                --min_gap 0 \
                --result_prefix "ablation_results"
                
            echo "✅ ${SAMPLE} ablation done."
        else
            echo "⏩ Output ${OUTPUT_DIR} already exists, skipping ${SAMPLE}."
        fi
    fi
done

echo "🎉 All ablation pipeline jobs finished."