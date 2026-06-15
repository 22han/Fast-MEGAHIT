#!/bin/bash
#===============================================================================
# 真实数据集质量评估脚本 (统一 Reads 版本)
# 功能：CheckV 完整性评估 + 统一使用 Original Reads 的 Read Mapping 回贴率统计
# 特点：所有方法均映射到同一套 Original Reads，确保可比性；自动清理中间文件，详细日志记录
#===============================================================================

set -euo pipefail

# 加载 Conda 初始化函数
source /home/zhangzihan/miniconda3/etc/profile.d/conda.sh

#-------------------------------------------------------------------------------
# 配置区 - 请根据实际情况修改
#-------------------------------------------------------------------------------

# Clean Reads (仅用于 OLI 组装，不再用于 mapping)
READS_1_OLI="/data2/zhangzihan2/SRR22936639clean_1.fastq"
READS_2_OLI="/data2/zhangzihan2/SRR22936639clean_2.fastq"

# Original Reads (所有方法的 mapping 统一使用此套 Reads)
READS_1="/data2/zhangzihan2/SRR22936639_1.fastq"
READS_2="/data2/zhangzihan2/SRR22936639_2.fastq"

# 三种方法的组装结果
OLI_CONTIGS="/data2/zhangzihan2/SRR22936639_oli_custom/final.contigs.fa"
DEFAULT_CONTIGS="/data2/zhangzihan2/SRR22936639_ou/final.contigs.fa"
COMPARISON_CONTIGS="/data2/zhangzihan2/SRR22936639_ou1/final.contigs.fa"

# 输出目录
OUTPUT_DIR="/data2/zhangzihan2/SRR22936639_analysis"
CHECKV_DIR="${OUTPUT_DIR}/checkv_results"
MAPPING_DIR="${OUTPUT_DIR}/mapping_results"

# 线程数
THREADS=16

# Conda 环境名称
CHECKV_ENV="checkv"
MAPPING_ENV="mapping_env"

#-------------------------------------------------------------------------------
# 函数定义
#-------------------------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

cleanup() {
    log "清理可能存在的旧文件..."
    if [ -d "$CHECKV_DIR" ]; then
        rm -rf "$CHECKV_DIR"
    fi
    if [ -d "$MAPPING_DIR" ]; then
        rm -rf "$MAPPING_DIR"
    fi
    mkdir -p "$CHECKV_DIR" "$MAPPING_DIR"
}

check_files() {
    log "检查输入文件是否存在..."
    local missing=0
    
    # 检查组装文件
    for file in "$OLI_CONTIGS" "$DEFAULT_CONTIGS" "$COMPARISON_CONTIGS"; do
        if [ ! -f "$file" ]; then
            log "ERROR: 组装文件不存在: $file"
            missing=1
        fi
    done
    
    # 检查 Clean Reads (OLI 组装用)
    if [ ! -f "$READS_1_OLI" ] || [ ! -f "$READS_2_OLI" ]; then
        log "ERROR: OLI Clean Reads 不存在"
        missing=1
    fi
    
    # 检查 Original Reads (所有 mapping 用)
    if [ ! -f "$READS_1" ] || [ ! -f "$READS_2" ]; then
        log "ERROR: Original Reads 不存在"
        missing=1
    fi
    
    if [ $missing -eq 1 ]; then
        log "ERROR: 部分输入文件缺失，程序终止"
        exit 1
    fi
    log "✓ 所有输入文件检查通过"
}

run_checkv() {
    local contigs=$1
    local method=$2
    local output="${CHECKV_DIR}/${method}"
    
    log "Running CheckV for ${method}..."
    
    conda activate "$CHECKV_ENV"
    checkv end_to_end "$contigs" "$output" -t "$THREADS" --quiet
    conda deactivate
    
    # 统计高质量基因组数量 (修复了 grep -c 的换行问题)
    local summary=""
    if [ -f "${CHECKV_DIR}/${method}/quality_summary.tsv" ]; then
        summary="${CHECKV_DIR}/${method}/quality_summary.tsv"
    elif [ -f "${CHECKV_DIR}/${method}/summary.tsv" ]; then
        summary="${CHECKV_DIR}/${method}/summary.tsv"
    elif [ -f "${CHECKV_DIR}/${method}/completeness.tsv" ]; then
        summary="${CHECKV_DIR}/${method}/completeness.tsv"
    fi
    
    if [ -f "$summary" ]; then
        local complete=0 high_quality=0 medium_quality=0
        complete=$(grep -c "Complete" "$summary" 2>/dev/null) || true
        high_quality=$(grep -c "High-quality" "$summary" 2>/dev/null) || true
        medium_quality=$(grep -c "Medium-quality" "$summary" 2>/dev/null) || true
        log "✓ ${method}: Complete=${complete}, High-quality=${high_quality}, Medium-quality=${medium_quality}"
    else
        log "⚠️ WARNING: CheckV summary not found for ${method}"
        if [ -d "$output" ]; then
            log "   输出目录内容："
            ls -la "$output" | head -20
        else
            log "   输出目录不存在: $output"
        fi
    fi
}

run_mapping() {
    # 参数: contigs, method
    local contigs=$1
    local method=$2
    local index="${MAPPING_DIR}/index_${method}"
    local bam="${MAPPING_DIR}/${method}.bam"
    local flagstat="${MAPPING_DIR}/${method}_flagstat.txt"
    
    log "Running Bowtie2 mapping for ${method} (using Original Reads)..."
    
    conda activate "$MAPPING_ENV"
    
    log "  Building index..."
    bowtie2-build "$contigs" "$index" > /dev/null 2>&1
    
    log "  Mapping reads..."
    bowtie2 -x "$index" \
            -1 "$READS_1" \
            -2 "$READS_2" \
            -p "$THREADS" \
            --quiet 2>/dev/null | \
        samtools view -bS - > "$bam"
    
    samtools flagstat "$bam" > "$flagstat"
    
    local mapped=$(grep "mapped (" "$flagstat" | head -1 | awk '{print $1}')
    local total=$(grep "in total" "$flagstat" | awk '{print $1}')
    
    conda deactivate
    
    log "✓ ${method}: Mapped=${mapped}/${total}"
}

generate_summary() {
    log "生成汇总报告..."
    
    local summary_file="${OUTPUT_DIR}/prefilter_summary.txt"
    
    {
        echo "==============================================================================="
        echo "                    真实数据集质量评估汇总报告"
        echo "                    Generated: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "==============================================================================="
        echo ""
        
        # CheckV 结果
        echo "【CheckV 基因组完整性评估】"
        echo "-------------------------------------------------------------------------------"
        printf "%-15s %12s %15s %18s\n" "Method" "Complete" "High-quality" "Medium-quality"
        echo "-------------------------------------------------------------------------------"
        
        for method in "oli" "default" "comparison"; do
            local sumf=""
            if [ -f "${CHECKV_DIR}/${method}/quality_summary.tsv" ]; then
                sumf="${CHECKV_DIR}/${method}/quality_summary.tsv"
            elif [ -f "${CHECKV_DIR}/${method}/summary.tsv" ]; then
                sumf="${CHECKV_DIR}/${method}/summary.tsv"
            elif [ -f "${CHECKV_DIR}/${method}/completeness.tsv" ]; then
                sumf="${CHECKV_DIR}/${method}/completeness.tsv"
            fi
            
            if [ -n "$sumf" ]; then
                local complete=0 high=0 medium=0
                complete=$(grep -c "Complete" "$sumf" 2>/dev/null) || true
                high=$(grep -c "High-quality" "$sumf" 2>/dev/null) || true
                medium=$(grep -c "Medium-quality" "$sumf" 2>/dev/null) || true
                printf "%-15s %12s %15s %18s\n" "$method" "$complete" "$high" "$medium"
            else
                printf "%-15s %12s %15s %18s\n" "$method" "N/A" "N/A" "N/A"
            fi
        done
        
        echo ""
        
        # Mapping 结果
        echo "【Read Mapping 回贴率评估】"
        echo "所有方法均使用同一套 Original Reads 进行 mapping，确保可比性"
        echo "-------------------------------------------------------------------------------"
        printf "%-15s %15s %15s %12s\n" "Method" "Total Reads" "Mapped" "Rate (%)"
        echo "-------------------------------------------------------------------------------"
        
        for method in "oli" "default" "comparison"; do
            local flagstat="${MAPPING_DIR}/${method}_flagstat.txt"
            if [ -f "$flagstat" ]; then
                local total=$(grep "in total" "$flagstat" | awk '{print $1}')
                local mapped=$(grep "mapped (" "$flagstat" | head -1 | awk '{print $1}')
                local rate=$(echo "scale=2; $mapped * 100 / $total" | bc 2>/dev/null || echo "0.00")
                printf "%-15s %15s %15s %11s%%\n" "$method" "$total" "$mapped" "$rate"
            else
                printf "%-15s %15s %15s %12s\n" "$method" "N/A" "N/A" "N/A"
            fi
        done
        
        echo ""
        echo "==============================================================================="
        echo "                           详细文件位置"
        echo "==============================================================================="
        echo "CheckV 结果目录：${CHECKV_DIR}"
        echo "Mapping 结果目录：${MAPPING_DIR}"
        echo ""
        echo "==============================================================================="
        
    } > "$summary_file"
    
    log "✓ 汇总报告已保存至：$summary_file"
    cat "$summary_file"
}

#-------------------------------------------------------------------------------
# 主流程
#-------------------------------------------------------------------------------

main() {
    log "=========================================="
    log "   前置分析脚本启动 (统一 Reads 版本)"
    log "=========================================="
    
    check_files
    cleanup
    mkdir -p "$OUTPUT_DIR" "$CHECKV_DIR" "$MAPPING_DIR"
    
    log ">>> 步骤 1: CheckV 完整性评估..."
    run_checkv "$OLI_CONTIGS" "oli"
    run_checkv "$DEFAULT_CONTIGS" "default"
    run_checkv "$COMPARISON_CONTIGS" "comparison"
    
    log ">>> 步骤 2: 统一使用 Original Reads 进行 Read Mapping..."
    run_mapping "$OLI_CONTIGS" "oli"
    run_mapping "$DEFAULT_CONTIGS" "default"
    run_mapping "$COMPARISON_CONTIGS" "comparison"
    
    log ">>> 步骤 3: 生成最终报告..."
    generate_summary
    
    log "=========================================="
    log "   前置分析完成 ✅"
    log "=========================================="
}

# 运行主函数
main "$@"