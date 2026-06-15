#!/bin/bash
# batch_kmc.sh - 批量运行KMC于 /data2/zhangzihan1 下的所有双端样本
# 若检测到最终结果 k141_spectrum.txt 已存在，则跳过该样本

# ===============================
# 固定参数（可按需修改）
# ===============================
THREADS=16
MEMORY=32
FREQ=1200
START_K=21
END_K=141
STEP=2

# 输入目录
DATA_DIR="/data2/zhangzihan"

# ===============================
# 查找所有 clean_1.fastq 文件
# ===============================
for R1 in "${DATA_DIR}"/*clean_1.fastq; do

    # 若没有匹配文件，避免把通配符当字符串执行
    [ -e "$R1" ] || continue

    # 提取样本名前缀
    # 例如：
    # ERR14788902clean_1.fastq -> ERR14788902
    base=$(basename "$R1" "clean_1.fastq")

    # 对应R2文件
    R2="${DATA_DIR}/${base}clean_2_rc.fastq"

    # 检查R2是否存在
    if [ ! -f "$R2" ]; then
        echo "警告: 未找到对应R2文件，跳过样本 $base"
        echo "缺失文件: $R2"
        continue
    fi

    # 输出目录
    OUTDIR="${DATA_DIR}/${base}_kmer"
    mkdir -p "$OUTDIR"

    # ==========================================
    # 如果最终结果已存在，则整个样本跳过
    # ==========================================
    if [ -f "${OUTDIR}/k141_spectrum.txt" ]; then
        echo "========================================="
        echo "样本 $base 已完成，跳过"
        echo "检测到文件: ${OUTDIR}/k141_spectrum.txt"
        echo "========================================="
        continue
    fi

    echo "========================================="
    echo "开始处理样本: $base"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "输出目录: $OUTDIR"
    echo "========================================="

    # 输入列表文件（给 KMC 使用）
    LIST_FILE="${OUTDIR}/fastq_files.lst"
    > "$LIST_FILE"
    echo "$R1" >> "$LIST_FILE"
    echo "$R2" >> "$LIST_FILE"

    # 临时目录
    TMP_BASE="${OUTDIR}/tmp"
    mkdir -p "$TMP_BASE"

    # ==========================================
    # 循环运行不同 k 值
    # ==========================================
    for (( k=$START_K; k<=$END_K; k+=$STEP )); do

        echo "当前 k=$k"

        TMP_DIR="${TMP_BASE}/k${k}"
        mkdir -p "$TMP_DIR"

        DB_PREFIX="${OUTDIR}/kmc_db_k${k}"

        # 运行 KMC
        kmc -k${k} \
            -ci1 \
            -cs${FREQ} \
            -t${THREADS} \
            -m${MEMORY} \
            "@${LIST_FILE}" \
            "${DB_PREFIX}" \
            "${TMP_DIR}"

        if [ $? -ne 0 ]; then
            echo "错误: kmc 在 k=$k 失败，终止"
            exit 1
        fi

        # 输出 histogram
        kmc_tools transform \
            "${DB_PREFIX}" \
            histogram \
            "${OUTDIR}/k${k}_spectrum.txt"

        if [ $? -ne 0 ]; then
            echo "错误: kmc_tools 在 k=$k 失败，终止"
            exit 1
        fi

        # 删除数据库文件
        rm -f "${DB_PREFIX}.kmc_pre"
        rm -f "${DB_PREFIX}.kmc_suf"

    done

    # 删除临时文件
    rm -rf "$TMP_BASE"
    rm -f "$LIST_FILE"

    echo "样本 $base 处理完成"
    echo ""

done

echo "所有样本处理完毕！"