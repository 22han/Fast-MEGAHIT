#!/bin/bash
set -euo pipefail

# =========================
# 配置区
# =========================

THREADS=16

# 直接指向包含多个 .fna 文件的文件夹路径（这是完全正确的）
REF_DIR="/data2/zhangzihan/source_genomes"

INPUT_OLI="/data2/zhangzihan/sample_1_reads/sample_oli/final.contigs.fa"
INPUT_OU="/data2/zhangzihan/sample_1_reads/sample_ou/final.contigs.fa"
INPUT_OU1="/data2/zhangzihan/sample_1_reads/sample_ou1/final.contigs.fa"

OUT_OLI="/data2/zhangzihan/sample_1_reads/quast_result_oli"
OUT_OU="/data2/zhangzihan/sample_1_reads/quast_result_ou"
OUT_OU1="/data2/zhangzihan/sample_1_reads/quast_result_ou1"

FINAL_TABLE="/data2/zhangzihan/sample_1_reads/comparison_metrics.tsv"

# =========================
# 运行 MetaQUAST
# =========================

run_metaquast () {
    local INPUT=$1
    local OUTPUT=$2
    local LABEL=$3

    echo ">>> Running $LABEL"

    # 清理可能存在的上次中断留下的半成品目录，防止报错或读取旧数据
    if [ -d "$OUTPUT" ]; then
        echo "    Cleaning old output directory: $OUTPUT"
        rm -rf "$OUTPUT"
    fi

    # 核心改动：加入了 --fast 参数！
    # 它会跳过对巨大参考基因组的GC计算等耗时操作，极大提升速度
    metaquast.py \
        -r "$REF_DIR" \
        -o "$OUTPUT" \
        -t $THREADS \
        --min-contig 500 \
        --fast \
        "$INPUT"
}

# =========================
# 提取指标
# =========================

extract_metrics () {
    local REPORT=$1

    if [ ! -f "$REPORT" ]; then
        echo -e "NA\tNA\tNA\tNA"
        return
    fi

    # 使用 -m1 防止匹配到多行
    local GF=$(grep -m1 'Genome fraction (%)' "$REPORT" | cut -f2)
    local MIS=$(grep -m1 '# misassemblies' "$REPORT" | cut -f2)
    local N50=$(grep -m1 'N50' "$REPORT" | cut -f2)
    local CONTIGS=$(grep -m1 '# contigs' "$REPORT" | cut -f2)

    echo -e "${GF}\t${MIS}\t${N50}\t${CONTIGS}"
}

# =========================
# 统计 >90% genome recovery
# =========================

count_high_recovery () {
    local OUTDIR=$1
    local RUNS_PER_REF="$OUTDIR/runs_per_reference"

    if [ ! -d "$RUNS_PER_REF" ]; then
        echo "0"
        return
    fi

    local count=0

    for ref_dir in "$RUNS_PER_REF"/*/; do
        ref_dir=${ref_dir%/}
        local tsv="$ref_dir/report.tsv"
        
        if [ ! -f "$tsv" ]; then
            continue
        fi

        local gf=$(grep -m1 'Genome fraction (%)' "$tsv" | cut -f2)
        
        if [ -n "$gf" ]; then
            local ok=$(awk "BEGIN{print ($gf >= 90) ? 1 : 0}")
            count=$((count + ok))
        fi
    done

    echo "$count"
}

# =========================
# 主流程
# =========================

run_metaquast "$INPUT_OLI" "$OUT_OLI" "Your Method"
run_metaquast "$INPUT_OU" "$OUT_OU" "Default"
run_metaquast "$INPUT_OU1" "$OUT_OU1" "Comparison"

echo ">>> Extracting metrics..."

echo -e "Method\tGenomeFraction\tMisassemblies\tN50\tContigs\tGenome>90%" > "$FINAL_TABLE"

METRICS=$(extract_metrics "$OUT_OLI/combined_reference/report.tsv")
RECOV=$(count_high_recovery "$OUT_OLI")
echo -e "OLI\t$METRICS\t$RECOV" >> "$FINAL_TABLE"

METRICS=$(extract_metrics "$OUT_OU/combined_reference/report.tsv")
RECOV=$(count_high_recovery "$OUT_OU")
echo -e "Default\t$METRICS\t$RECOV" >> "$FINAL_TABLE"

METRICS=$(extract_metrics "$OUT_OU1/combined_reference/report.tsv")
RECOV=$(count_high_recovery "$OUT_OU1")
echo -e "Comparison\t$METRICS\t$RECOV" >> "$FINAL_TABLE"

echo "=============================="
echo "✅ Final comparison table:"
cat "$FINAL_TABLE"
echo "=============================="
