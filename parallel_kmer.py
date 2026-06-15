#!/usr/bin/env python3
"""
debug_single_k.py — 单 k 值参数调试工具

用途：加载单个 k-mer 频谱文件，运行 EM 拟合，打印所有参数和指标
适用于：调试 w_err/w_sig 异常、检查双峰参数、验证 cutoff 合理性

用法：
    python debug_single_k.py spectrum_k21.txt --k 21 --bimodal
    python debug_single_k.py spectrum_k21.txt --k 21 --simple --plot
"""

import sys
import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt

# 导入分析器类
try:
    from analyze_single_k import SparseKmerSpectrumAnalyzer
except ImportError:
    print("❌ 错误: 找不到 analyze_single_k.py")
    print("请确保该文件与 debug_single_k.py 在同一目录，或设置 PYTHONPATH")
    sys.exit(1)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def print_separator(title: str, char: str = "═", width: int = 70):
    """打印带标题的分隔线"""
    print(f"\n{char * width}")
    print(f" {title}")
    print(f"{char * width}")


def print_metric_section(title: str, metrics: dict, keys: list, fmt: str = "{:.4f}"):
    """打印指标分组"""
    print(f"\n📊 {title}")
    print("-" * 40)
    for key in keys:
        if key in metrics:
            val = metrics[key]
            if isinstance(val, float):
                print(f"  {key:25s} : {fmt.format(val)}")
            else:
                print(f"  {key:25s} : {val}")


def plot_spectrum(analyzer: SparseKmerSpectrumAnalyzer, save_path: str = None):
    """绘制频谱 + 后验概率图"""
    freq_vals, counts = analyzer._get_arrays()
    
    # 图1: 原始频谱
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Panel 1: 频谱
    ax1.bar(freq_vals, counts, width=0.8, alpha=0.7, color='steelblue', label='observed')
    ax1.set_xlabel('Frequency')
    ax1.set_ylabel('Count')
    ax1.set_title(f'K={analyzer.k_value} — Raw Spectrum')
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, max(freq_vals) * 1.1)
    
    # Panel 2: 后验概率
    x_vals, post_err = analyzer.compute_posterior()
    post_sig = 1.0 - post_err
    
    ax2.plot(x_vals, post_err, 'r-', label='P(error|x)', linewidth=2)
    ax2.plot(x_vals, post_sig, 'g-', label='P(signal|x)', linewidth=2)
    ax2.axvline(x=analyzer.cutoff, color='orange', linestyle='--', 
                label=f'cutoff={analyzer.cutoff}')
    ax2.set_xlabel('Frequency')
    ax2.set_ylabel('Posterior Probability')
    ax2.set_title('Posterior Probabilities (Error vs Signal)')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, max(freq_vals) * 1.1)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        logger.info(f"📈 诊断图已保存: {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="🔍 Debug single k-mer EM fitting parameters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基础用法（双峰模式）
  python debug_single_k.py spectrum_k21.txt --k 21
  
  # 单峰模式 + 保存诊断图
  python debug_single_k.py spectrum_k21.txt --k 21 --simple --plot --save-fig k21_diag.png
  
  # 自定义参数
  python debug_single_k.py spectrum_k21.txt --k 21 --max-freq 200 --tol 1e-5
        """
    )
    
    parser.add_argument('spectrum_file', help='频谱文件路径 (KmerLight/ntCard/KMC 格式)')
    parser.add_argument('--k', type=int, required=True, help='k-mer 长度值')
    parser.add_argument('--simple', action='store_true', help='使用简单初始化（非 robust）')
    parser.add_argument('--unimodal', action='store_true', help='使用单峰模型（默认双峰）')
    parser.add_argument('--max-freq', type=int, default=None, help='最大频率截断')
    parser.add_argument('--tol', type=float, default=1e-4, help='EM 收敛阈值')
    parser.add_argument('--max-iter', type=int, default=200, help='EM 最大迭代次数')
    parser.add_argument('--plot', action='store_true', help='绘制诊断图')
    parser.add_argument('--save-fig', type=str, default=None, help='保存诊断图路径')
    
    args = parser.parse_args()
    
    # 打印头部
    print_separator(f"🔍 Single K Debug — k={args.k}", "█")
    print(f"📁 频谱文件: {args.spectrum_file}")
    print(f"⚙️  模式: {'单峰 (unimodal)' if args.unimodal else '双峰 (bimodal)'}")
    print(f"🔧 初始化: {'simple' if args.simple else 'robust'}")
    
    # 初始化分析器
    try:
        analyzer = SparseKmerSpectrumAnalyzer(
            spectrum_file=args.spectrum_file,
            max_freq=args.max_freq,
            posterior_threshold=0.5,
            max_em_iter=args.max_iter,
            tol=args.tol
        )
        analyzer.k_value = args.k  # 用于绘图标识
    except Exception as e:
        logger.error(f"❌ 初始化失败: {e}")
        sys.exit(1)
    
    # 打印基础统计
    stats = analyzer.get_spectrum_stats()
    print_separator("📈 Spectrum Statistics")
    print(f"  distinct_kmers          : {stats['distinct_kmers']:,}")
    print(f"  total_occurrences       : {stats['total_kmer_occurrences']:,}")
    print(f"  frequency range         : {stats['min_frequency']} ~ {stats['max_frequency']}")
    print(f"  peak frequency          : {stats['peak_frequency']} (count: {stats['peak_count']:,})")
    
    # 运行完整分析
    print_separator("🔄 Running EM Fitting + Metrics")
    try:
        n_components = 1 if args.unimodal else 2
        method = 'simple' if args.simple else 'robust'
        
        metrics = analyzer.get_metrics(method=method, n_signal_components=n_components)
    except Exception as e:
        logger.error(f"❌ 分析失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # 打印收敛信息
    print_separator("✅ Convergence Status")
    print(f"  converged               : {metrics['em_converged']}")
    print(f"  iterations              : {metrics['em_iterations']}")
    print(f"  n_signal_components     : {metrics['n_signal_components']}")
    
    # 打印核心混合参数
    print_metric_section("🎯 Mixture Model Parameters", metrics, 
                        ['w_err', 'w_sig', 'p_err', 'mu', 'sigma'], 
                        fmt="{:.6f}")
    
    # 打印双峰特有参数
    if metrics['n_signal_components'] == 2:
        print_metric_section("🧬 Rare Species Parameters (Bimodal)", metrics,
                            ['mu_rare', 'sigma_rare', 'w_rare_frac'],
                            fmt="{:.6f}")
    
    # 打印 cutoff 相关
    print_metric_section("✂️  Cutoff & Effective Counts", metrics,
                        ['cutoff', 'error_sum', 'effective_diff', 'soft_effective_diff', 'distinct_kmers'],
                        fmt="{:.0f}" if 'diff' in metrics else "{:.4f}")
    
    # 打印创新指标
    print_metric_section("✨ Innovation Metrics (CRITIC Objectives)", metrics,
                        ['boundary_uncertainty'],
                        fmt="{:.6f}")
    
    # 打印派生指标（用于分区/选K）
    print("\n📐 Derived Metrics for Partitioning/Selection")
    print("-" * 40)
    info = metrics['soft_effective_diff'] / max(metrics['distinct_kmers'], 1)
    snr = metrics['w_sig'] / (metrics['w_err'] + 1e-10)
    unif = 1.0 / (1.0 + metrics['sigma'])
    conf = 1.0 / (1.0 + metrics['boundary_uncertainty'])
    
    print(f"  info (soft ratio)       : {info:.6f}")
    print(f"  snr (w_sig/w_err)       : {snr:.6f}")
    print(f"  unif (1/(1+sigma))      : {unif:.6f}")
    print(f"  conf (1/(1+uncertainty)): {conf:.6f}")
    
    # 诊断图
    if args.plot:
        print_separator("📊 Generating Diagnostic Plot")
        plot_spectrum(analyzer, save_path=args.save_fig)
    
    # 尾部总结
    print_separator("🏁 Summary", "★")
    if metrics['w_err'] > 0.5:
        print(f"⚠️  警告: w_err={metrics['w_err']:.4f} > 0.5，可能表示:")
        print("   • cutoff 设置过低，把稀有物种误判为错误")
        print("   • 频谱质量差或测序错误率高")
        print("   • 建议检查 cutoff 和 mu_rare 参数")
    elif metrics['w_err'] < 0.05:
        print(f"✓ w_err={metrics['w_err']:.4f} 合理（5% 以下）")
    else:
        print(f"✓ w_err={metrics['w_err']:.4f}, w_sig={metrics['w_sig']:.4f} 比例合理")
    
    print(f"\n💡 提示: 如果参数异常，尝试:")
    print("   • --simple 使用简单初始化")
    print("   • --unimodal 强制单峰模型")
    print("   • --max-freq 100 截断高频噪声")
    print("   • 检查频谱文件格式是否正确")
    
    print(f"\n✅ 调试完成！")


if __name__ == "__main__":
    main()