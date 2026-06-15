#!/bin/bash

# 样本名前缀（XXXX）
SAMPLES=("ERR16144122" "SRR36991061" "SRR37941741")

# 输入输出路径
BASE_DIR="/data2/zhangzihan1"

for SAMPLE in "${SAMPLES[@]}"
do
    echo "Processing ${SAMPLE}..."

    ./kmer.sh \
        ${BASE_DIR}/${SAMPLE}clean_1.fastq \
        ${BASE_DIR}/${SAMPLE}clean_2_rc.fastq \
        ${BASE_DIR}/${SAMPLE}_kmer

    echo "${SAMPLE} done."
done

echo "All samples finished."