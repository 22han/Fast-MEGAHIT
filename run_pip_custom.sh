#!/bin/bash

BASE_DIR="/data2/zhangzihan/sample_19_reads"

# 遍历所有以 _kmer 结尾的目录
for dir in "${BASE_DIR}"/*_kmer; do
    # 确保是目录
    if [ -d "$dir" ]; then
        # 提取样本名（去掉末尾的 _kmer）
        SAMPLE=$(basename "$dir" | sed 's/_kmer$//')
        OUTPUT_DIR="${BASE_DIR}/${SAMPLE}_k_nopesc"

        # 如果输出目录不存在，则运行 pipeline
        if [ ! -d "$OUTPUT_DIR" ]; then
            echo "Running custom pipeline for ${SAMPLE}..."
            python run_pipeline_no_pesc.py -i "$dir" -o "$OUTPUT_DIR" \
                --min_partitions 6 \
                --max_partitions 12 \
                --min_size 4 \
                --sensitivity 0.5
            echo "${SAMPLE} done."
        else
            echo "Output ${OUTPUT_DIR} already exists, skipping ${SAMPLE}."
        fi
    fi
done

echo "All custom pipeline jobs finished."
