#!/usr/bin/env python3
"""
run_pipeline_custom.py: 自定义参数版本的 K-mer 选择流水线
特点：支持通过命令行调整分区参数，不影响原始代码
"""
import sys, os, argparse, logging, glob
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kmer_batch import AdaptiveKmerSelector   

def setup_logging(log_file="run_pipeline_custom.log"):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)-8s | %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w', encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def resolve_spectrum_dir(input_path: str, logger) -> str:
    input_path = os.path.expanduser(input_path)
    if os.path.isdir(input_path):
        logger.info(f"📁 检测到目录: {input_path}")
        return os.path.abspath(input_path)
    if os.path.isfile(input_path):
        logger.info(f"📄 检测到文件: {input_path}")
        return os.path.dirname(os.path.abspath(input_path))
    if any(c in input_path for c in ['*', '?', '[']):
        matches = glob.glob(input_path)
        if not matches:
            logger.error(f"❌ 通配符无匹配: {input_path}")
            return None
        logger.info(f"🔍 匹配到 {len(matches)} 个文件，使用目录: {os.path.dirname(matches[0])}")
        return os.path.dirname(os.path.abspath(matches[0]))
    logger.error(f"❌ 路径不存在: {input_path}")
    return None

def check_spectrum_files(spectrum_dir: str, k_min: int, k_max: int, step: int, logger) -> bool:
    patterns = ["k{k}_spectrum.txt", "spectrum_k{k}.txt", "ntcard_k{k}.txt", "k{k}.txt"]
    found = 0
    for k in range(k_min, k_max + 1, step):
        for pat in patterns:
            if os.path.exists(os.path.join(spectrum_dir, pat.format(k=k))):
                found += 1
                break
    total = (k_max - k_min) // step + 1
    logger.info(f"📊 找到 {found}/{total} 个频谱文件 @ {spectrum_dir}")
    if found == 0:
        logger.error("❌ 未找到任何频谱文件")
        return False
    if found < total * 0.8:
        logger.warning(f"⚠️  仅找到 {found} 个文件，可能影响结果连续性")
    return True

def setup_output_dir(output_path: str, logger) -> str:
    output_dir = os.path.abspath(output_path)
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"📁 输出目录已就绪: {output_dir}")
    return output_dir

def main():
    temp_parser = argparse.ArgumentParser(add_help=False)
    temp_parser.add_argument("--output", "-o", type=str)
    temp_args, _ = temp_parser.parse_known_args()
    
    log_file = "run_pipeline_custom.log"
    if temp_args.output:
        os.makedirs(os.path.abspath(temp_args.output), exist_ok=True)
        log_file = os.path.join(os.path.abspath(temp_args.output), "run_pipeline_custom.log")
    
    logger = setup_logging(log_file=log_file)
    
    parser = argparse.ArgumentParser(
        description="🧬 自适应 K-mer 选择流水线 (自定义参数版)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # 核心参数
    parser.add_argument("--input", "-i", type=str, required=True, help="频谱文件路径")
    parser.add_argument("--output", "-o", type=str, required=True, help="输出目录路径")
    
    # K值范围参数
    parser.add_argument("--k_min", type=int, default=21)
    parser.add_argument("--k_max", type=int, default=141)
    parser.add_argument("--step", type=int, default=2)
    
    # 运行参数
    parser.add_argument("--workers", "-w", type=int, default=4)
    
    # === 自定义分区参数 ===
    parser.add_argument("--min_partitions", type=int, default=6,
                        help="最小分区数，默认: 6")
    parser.add_argument("--max_partitions", type=int, default=12,
                        help="最大分区数，默认: 12")
    parser.add_argument("--min_size", type=int, default=4,
                        help="每个分区最小k值数量，默认: 4")
    parser.add_argument("--sensitivity", type=float, default=0.5,
                        help="自适应扩张灵敏度，默认: 0.5（值越大分区越多）")
    
    parser.add_argument("--result_prefix", type=str, default="kmer_results_custom",
                        help="结果文件名前缀")
    parser.add_argument("--skip_diagnose", action="store_true", help="跳过诊断图表")
    parser.add_argument("--ablation", action="store_true", 
                        help="开启消融模式：跳过分区，直接全局CRITIC选K")
    parser.add_argument("--n_top_k", type=int, default=5,
                        help="消融模式下全局选择的k值数量 (默认: 5)")
    parser.add_argument("--min_gap", type=int, default=10,
                        help="消融模式下防止k值聚集的最小间隔 (默认: 10)")
    
    args = parser.parse_args()

    # 1. 解析输入路径
    spectrum_dir = resolve_spectrum_dir(args.input, logger)
    if not spectrum_dir:
        sys.exit(1)
    
    # 2. 验证频谱文件
    if not check_spectrum_files(spectrum_dir, args.k_min, args.k_max, args.step, logger):
        sys.exit(1)

    # 3. 创建输出目录
    output_dir = setup_output_dir(args.output, logger)

    # 4. 初始化 Selector
    logger.info(f"⚙️  初始化 AdaptiveKmerSelector (k={args.k_min}:{args.k_max}:{args.step})")
    selector = AdaptiveKmerSelector(
        k_range=(args.k_min, args.k_max),
        step=args.step,
        n_workers=args.workers
    )

    # 5. 执行流水线（分步调用，使用自定义参数）
    try:
        logger.info("🚀 开始执行流水线...")
        
        # Stage 1: EM拟合 (保持不变)
        logger.info("[Stage 1] 执行 EM 拟合...")
        ok = selector.fit_em_frequencies(
            spectrum_dir=spectrum_dir,
            init_method='robust',
            n_signal_components=2
        )
        if not ok:
            logger.error("❌ EM 拟合失败")
            sys.exit(1)
        
        # ====================================================================
        # 🧪 分支：消融模式 vs 正常模式 (替换原来的 Stage 1.5 和 Stage 2)
        # ====================================================================
        if args.ablation:
            logger.info(f"🧪 [ABLATION MODE] 激活！跳过分区，直接全局选择 Top {args.n_top_k} 个 k 值...")
            
            # 1. 收集所有收敛的 metrics，计算全局 4 目标指标
            feasible_candidates = []
            for metrics in selector._stage1_metrics:
                if not metrics or not metrics.get('em_converged'):
                    continue
                
                mu = metrics['mu']
                mu_rare = metrics.get('mu_rare', 0.1 * mu)
                sigma = metrics['sigma']
                sigma_rare = metrics.get('sigma_rare', sigma * 0.5)
                
                base_value = (
                    abs(mu - mu_rare) / (sigma + sigma_rare + 1e-6) * 
                    min(mu / max(mu_rare, 1e-6), 3.0) / 
                    max(metrics.get('coverage_heterogeneity', 1.0), 1.0)
                )
                species_resolution_value = base_value if metrics.get('n_signal_components', 1) == 2 else 0.5
                
                candidate = {
                    'k': metrics['k'],
                    'node_retention': metrics['soft_effective_diff'] / metrics.get('distinct_kmers', 1),
                    'error_mass': 1.0 / (1.0 + 0.3*metrics.get('contamination_proxy', 0) + 
                                         0.2 * (1.0 / max(metrics.get('p_err', 0.3), 0.05))),
                    'species_resolution': species_resolution_value,
                    'cutoff_sharpness': np.exp(-metrics['boundary_uncertainty']),
                }
                feasible_candidates.append(candidate)
                
            if not feasible_candidates:
                logger.error("❌ 消融模式：没有找到收敛的候选k值")
                sys.exit(1)
                
            candidates_df = pd.DataFrame(feasible_candidates)
            
            # 2. 直接调用 selector 内部的 CRITIC 权重计算和归一化方法 (Python允许调用_开头的方法)
            global_weights = selector._compute_critic_weights(candidates_df)
            logger.info(f"🌍 Global CRITIC Weights: {global_weights}")
            
            for obj in ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']:
                if obj == 'species_resolution':
                    candidates_df[f'{obj}_log'] = np.log1p(candidates_df[obj])
                    candidates_df[f'{obj}_normalized'] = selector._apply_iqr_normalization(candidates_df[f'{obj}_log'])
                else:
                    candidates_df[f'{obj}_normalized'] = selector._apply_iqr_normalization(candidates_df[obj])
                    
            candidates_df['final_score'] = (
                global_weights['node_retention'] * candidates_df['node_retention_normalized'] +
                global_weights['error_mass']  * candidates_df['error_mass_normalized'] +
                global_weights['species_resolution'] * candidates_df['species_resolution_normalized'] +
                global_weights['cutoff_sharpness'] * candidates_df['cutoff_sharpness_normalized']
            )
            
            # 3. 贪心选择 Top N (带 min_gap 约束，防止选出 31,33,35,37)
            sorted_df = candidates_df.sort_values('final_score', ascending=False).reset_index(drop=True)
            selected_ks = []
            for _, row in sorted_df.iterrows():
                k_val = int(row['k'])
                if not selected_ks:
                    selected_ks.append(k_val)
                else:
                    if all(abs(k_val - sel_k) >= args.min_gap for sel_k in selected_ks):
                        selected_ks.append(k_val)
                if len(selected_ks) >= args.n_top_k:
                    break
                    
            # 如果因为间隔限制选不够，则放宽限制补齐
            if len(selected_ks) < args.n_top_k:
                logger.warning(f"因 min_gap={args.min_gap} 限制仅选出 {len(selected_ks)} 个，放宽限制补齐...")
                for _, row in sorted_df.iterrows():
                    k_val = int(row['k'])
                    if k_val not in selected_ks:
                        selected_ks.append(k_val)
                    if len(selected_ks) >= args.n_top_k:
                        break
                        
            selected_ks.sort()
            selector.selected_kmers = selected_ks  # 同步给 selector 以便后续保存和画图
            logger.info(f"🎯 [ABLATION] 全局选中 {len(selected_ks)} 个 k 值: {selected_ks}")

        else:
            # 正常模式：Stage 1.5 + Stage 2
            logger.info(f"[Stage 1.5] 构建自适应分区 (min_partitions={args.min_partitions}, sensitivity={args.sensitivity})")
            selector.construct_adaptive_partitions(
                min_partitions=args.min_partitions,
                max_partitions=args.max_partitions,
                min_size=args.min_size,
                sensitivity=args.sensitivity
            )
            
            logger.info("[Stage 2] 执行多目标选择...")
            selected_ks = selector.execute_multi_objective_selection(
                method='balanced',
                use_critic=True
            )
            
        # 诊断图表
        diagnose_path = os.path.join(output_dir, f"{args.result_prefix}_diagnosis.png")
        if not args.skip_diagnose:
            selector.diagnose_partitioning(save_path=diagnose_path)

        if not selected_ks:
            logger.error("❌ 未返回有效 k 值")
            sys.exit(1)

        # 保存结果
        logger.info("💾 保存结果...")
        result_prefix = os.path.join(output_dir, args.result_prefix)
        selector.save_results(prefix=result_prefix)
        
        # 列出生成的文件
        logger.info("✅ 流水线执行成功！")
        logger.info(f"🎯 选中 k 值: {selected_ks}")
        logger.info(f"📁 所有输出文件位于: {output_dir}/")
        saved_files = [f for f in os.listdir(output_dir) if f.startswith(args.result_prefix)]
        for f in sorted(saved_files):
            fpath = os.path.join(output_dir, f)
            size = os.path.getsize(fpath) / 1024
            logger.info(f"   ✓ {f} ({size:.1f} KB)")
        
    except KeyboardInterrupt:
        logger.warning("⏹️  用户中断")
        sys.exit(130)
    except Exception as e:
        logger.exception(f"💥 执行失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
