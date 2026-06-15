#!/usr/bin/env python3
"""
AdaptiveKmerSelector: CRITIC-based Multi-Objective Optimization Framework
for K-mer Length Selection in De Bruijn Graph Assembly

Architecture: Two-Stage Fidelity-Preserving Pipeline
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Stage 1 (Coarse, NO cross-k prior):
        → soft_effective_diff  → change-point partitioning  (保真分区)
        → Stage 1 metrics      → CRITIC 4-objective selection (无污染选K)
Four CRITIC objectives (metagenomic-optimized):
    1. node_retention — soft_eff / (distinct_kmers + 0.1*error_sum)
    2. error_mass     — 1 / (1 + 0.5*w_err / p_err)
    3. species_resolution — |μ-μ_rare|/(σ+σ_rare) × peak_ratio (bimodal)
    4. cutoff_sharpness   — exp(-boundary_uncertainty)

Partitioning:
    Default method is 'change_point' — detects derivative discontinuities
    in the soft_effective_diff profile, suitable for monotonic curves.
    Falls back to 'peaks' or fixed partitioning if insufficient structure.

References:
    Diakoulaki et al. (1995). COR, 22(7), 763-770.
    Dempster et al. (1977). JRSS-B, 39(1), 1-22.
    Kirk et al. (2012). Intro to Statistical Quality Control, 7th Ed.
"""

import os
import sys
import json
import logging
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')                          # [FIX] non-interactive backend
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import pearsonr
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from typing import Dict, List, Tuple, Optional

# [FIX] make import flexible (已替换为 analysis)
try:
    from analysis import SparseKmerSpectrumAnalyzer
except ImportError:
    from kmer_analyzer import SparseKmerSpectrumAnalyzer

# ------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------
# 移除这里的 basicConfig，统一由 run_pipeline.py 中的 setup_logging 管理
logger = logging.getLogger(__name__)


# ==================================================================
# Worker function for parallel EM fitting
# ==================================================================
def _fit_kmer_spectrum_worker(k: int, spectrum_file: Optional[str], init_method: str = 'robust', 
                              n_signal_components: int = 2) -> Optional[Dict]:
    """ 
    Parallel worker: fits mixture model for a single k-value.
    [FIX] 强制传入 n_signal_components 控制单/双峰
    [FIX] 移除了 Stage 2 相关的 prev_params 和 smooth_alpha 参数
    """
    try:
        if not spectrum_file or not os.path.exists(spectrum_file):
            logger.debug(f"k={k}: spectrum file not found")
            return None
            
        # 不再传入 prev_params 和 smooth_alpha
        analyzer = SparseKmerSpectrumAnalyzer(
            spectrum_file=spectrum_file,
            max_freq=None,
            posterior_threshold=0.5,
            max_em_iter=200,
            tol=1e-4,
        )
        analyzer.k_value = k
        
        # 强制传入 n_signal_components (宏基因组固定为 2)
        metrics = analyzer.get_metrics(method=init_method, n_signal_components=n_signal_components)
        
        if metrics and metrics.get('soft_effective_diff', 0) <= 0:
            logger.warning(f"k={k}: non-positive soft effective distinctness")
            return None
        return metrics
    except Exception as e:
        logger.exception(f"k={k} failed: {e}")
        return None

# ==================================================================
# Main class
# ==================================================================
class AdaptiveKmerSelector:
    """
    Two-stage adaptive k-mer selector with CRITIC 4-objective optimization.
    """

    def __init__(self, k_range: Tuple[int, int] = (21, 141),
                 step: int = 2, n_workers: int = 8,
                 read_length: Optional[int] = None):
        self.k_values = list(range(k_range[0], k_range[1] + 1, step))
        self.n_workers = n_workers
        self.read_length = read_length

        # Results storage
        self.metrics_list: List[Optional[Dict]] = []
        self.effective_profiles: List[float] = []
        self.adaptive_partitions: List[Tuple[int, int]] = []
        self.selected_kmers: List[int] = []

        # Configuration tracking
        self.partition_config: Dict = {}
        self.selection_config: Dict = {}
        self.critic_weights_history: List[Dict] = []

        # ==============================================================
        # [FIX] Stage 1 results — frozen, never overwritten by Stage 2
        # All decision-making (partitioning + CRITIC) reads from here.
        # ==============================================================
        self._stage1_metrics: List[Optional[Dict]] = []
        self._stage1_profiles: List[float] = []

    # ------------------------------------------------------------------
    # File lookup
    # ------------------------------------------------------------------
    def _locate_spectrum_file(self, directory: str, k: int) -> Optional[str]:
        naming_patterns = [
            f"spectrum_k{k}.txt",
            f"k{k}_spectrum.txt",
            f"ntcard_k{k}.txt",
            f"k{k}.txt"
        ]
        for pattern in naming_patterns:
            filepath = os.path.join(directory, pattern)
            if os.path.exists(filepath):
                return filepath
        return None
    
    def fit_em_frequencies(self, spectrum_dir: str, init_method: str = 'robust', n_signal_components: int = 2) -> bool:
        """Stage 1: parallel EM fitting without any cross-k prior."""
        logger.info(
            f"Parallel EM fitting for {len(self.k_values)} k-values "
            f"(mode={'BIMODAL' if n_signal_components==2 else 'UNIMODAL'})...")
            
        tasks = [(k, self._locate_spectrum_file(spectrum_dir, k)) for k in self.k_values]
        results = {}
        if tasks:
            logger.info(f"Total tasks: {len(tasks)}. K-range covered: from k={tasks[0][0]} to k={tasks[-1][0]}")
        else:
            logger.warning("No valid spectrum tasks generated!")
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # 把 n_signal_components 传给 worker
            futures = {
                executor.submit(_fit_kmer_spectrum_worker, k, sf, init_method, n_signal_components): k 
                for k, sf in tasks
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="EM Fitting"):
                k = futures[future]
                try:
                    result = future.result()
                    results[k] = result
                except Exception as e:
                    logger.error(f"Future k={k} failed: {e}")
                    results[k] = None
                    
        self.metrics_list = [results.get(k) for k in self.k_values]
        self.effective_profiles = [
            m['soft_effective_diff'] if m else np.nan 
            for m in self.metrics_list
        ]
        
        # Freeze Stage 1
        self._stage1_metrics = [m.copy() if m else None for m in self.metrics_list]
        self._stage1_profiles = self.effective_profiles.copy()
        
        converged_count = sum(1 for m in self.metrics_list if m and m.get('em_converged'))
        logger.info(f"✅ EM fitting complete: {converged_count}/{len(self.k_values)} converged")
        return converged_count > 0

    # ==================================================================
    # [FIX] CRITIC weight computation — 4 objectives
    # ==================================================================
    def construct_adaptive_partitions(
        self,
        min_partitions: int = 8,
        max_partitions: int = 15,
        min_size: int = 4,
        sensitivity: float = 0.4
    ) -> List[Tuple[int, int]]:
        """
        自适应分区：Top-Down 斜率优先 + 数据驱动扩张

        阶段 1：强制切到 min_partitions
        阶段 1.5：自适应扩张，若剩余段优先级仍高于阈值则继续切
        阶段 2：收缩合并（超过 max_partitions 时）
        """
        self.partition_config = {
            'method': 'top_down_adaptive',
            'min_partitions': min_partitions,
            'max_partitions': max_partitions,
            'min_size': min_size,
            'sensitivity': sensitivity
        }

        scores = np.array(self._stage1_profiles, dtype=float)
        valid_mask = ~np.isnan(scores)
        if np.sum(valid_mask) < 10:
            n = len(self.k_values)
            step = max(1, n // min_partitions)
            return [(i * step, min((i + 1) * step - 1, n - 1))
                    for i in range(min_partitions)]

        partitions = [(0, len(self.k_values) - 1)]
        cut_priorities_history = []

        def _compute_priorities():
            prios = []
            for s, e in partitions:
                y = scores[s:e + 1]
                x = np.arange(len(y))
                if len(x) < 3:
                    prios.append(0.0)
                    continue
                A = np.column_stack([x, np.ones(len(x))])
                coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
                slope = coef[0]
                y_pred = A @ coef
                mse = np.mean((y - y_pred) ** 2)
                y_range = scores.max() - scores.min() + 1e-10
                slope_priority = abs(slope) / y_range * 50.0
                length_priority = (e - s) * 2.0
                mse_priority = mse * 1000.0
                prios.append(slope_priority + length_priority + mse_priority)
            return prios

        def _find_split(s, e):
            best_split = None
            min_sum_mse = float('inf')
            mse_candidates = []
            for split in range(s + min_size, e - min_size + 1):
                mse_left = self._calculate_segment_mse(s, split, scores)
                mse_right = self._calculate_segment_mse(split + 1, e, scores)
                total = mse_left + mse_right
                mse_candidates.append((split, total))
                if total < min_sum_mse:
                    min_sum_mse = total
                    best_split = split
            if best_split is None or min_sum_mse < 1e-10:
                max_var = -1
                for split, _ in mse_candidates:
                    left_var = np.var(scores[s:split + 1])
                    right_var = np.var(scores[split + 1:e + 1])
                    if left_var + right_var > max_var:
                        max_var = left_var + right_var
                        best_split = split
            if best_split is None:
                best_split = s + int((e - s) * 0.618)
            return best_split

        def _do_cut(priorities):
            worst_idx = priorities.index(max(priorities))
            s, e = partitions[worst_idx]
            if e - s < 2 * min_size:
                lengths = [(i, p[1] - p[0]) for i, p in enumerate(partitions)]
                lengths.sort(key=lambda x: x[1], reverse=True)
                found = False
                for idx, length in lengths:
                    if length >= 2 * min_size:
                        worst_idx = idx
                        s, e = partitions[idx]
                        found = True
                        break
                if not found:
                    return False
            split = _find_split(s, e)
            partitions[worst_idx] = (s, split)
            partitions.insert(worst_idx + 1, (split + 1, e))
            return True

        # 阶段 1：强制切到 min_partitions
        iteration = 0
        while len(partitions) < min_partitions and iteration < 200:
            iteration += 1
            priorities = _compute_priorities()
            cut_priorities_history.append(max(priorities))
            if not _do_cut(priorities):
                break

        # 阶段 1.5：自适应扩张
        if len(partitions) >= min_partitions and len(partitions) < max_partitions:
            if cut_priorities_history:
                threshold = np.median(cut_priorities_history) * sensitivity
            else:
                threshold = 0.0

            while len(partitions) < max_partitions:
                priorities = _compute_priorities()
                worst_priority = max(priorities)
                if worst_priority < threshold:
                    break
                cut_priorities_history.append(worst_priority)
                if not _do_cut(priorities):
                    break

            logger.info(
                f"Adaptive expansion: {min_partitions} -> {len(partitions)} "
                f"(threshold={threshold:.2f}, sensitivity={sensitivity})")

        # 阶段 2：收缩合并
        while len(partitions) > max_partitions:
            min_merge_cost = float('inf')
            merge_idx = -1
            for i in range(len(partitions) - 1):
                s1, e1 = partitions[i]
                s2, e2 = partitions[i + 1]
                cost = self._calculate_segment_mse(s1, e2, scores) - \
                       (self._calculate_segment_mse(s1, e1, scores) +
                        self._calculate_segment_mse(s2, e2, scores))
                if cost < min_merge_cost:
                    min_merge_cost = cost
                    merge_idx = i
            if merge_idx == -1:
                break
            s1, e1 = partitions[merge_idx]
            s2, e2 = partitions[merge_idx + 1]
            partitions[merge_idx] = (s1, e2)
            del partitions[merge_idx + 1]

        result = list(partitions)
        self.adaptive_partitions = result
        logger.info(
            f"Constructed {len(result)} adaptive partitions "
            f"(top_down_adaptive, range=[{min_partitions},{max_partitions}])")
        return self.adaptive_partitions
        
    def _calculate_segment_mse(self, start, end, scores):
        """计算指定区间的线性拟合 MSE（自底向上/Top-Down 分区专用辅助函数）"""
        y = scores[start:end+1]
        x = np.arange(len(y))
        if len(x) < 2:
            return 0.0
        A = np.column_stack([x, np.ones(len(x))])
        coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        return np.mean((y - A @ coef) ** 2)
        
    def _compute_critic_weights(
            self, candidates_df: pd.DataFrame) -> Dict[str, float]:
        """
        CRITIC objective weight determination.

        4 objectives: info, snr, unif, conf
        Weight = contrast_intensity (std) × conflict (1 - |corr|)
        """
        objectives = ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']
        data_matrix = candidates_df[objectives].values

        # Step 1: Min-Max归一化（直接对原始数据操作，符合CRITIC标准）
        normalized = np.zeros_like(data_matrix, dtype=float)
        for j in range(data_matrix.shape[1]):
            col = data_matrix[:, j]  # 直接使用原始数据
            min_val, max_val = col.min(), col.max()
            if max_val > min_val:
                normalized[:, j] = (col - min_val) / (max_val - min_val + 1e-10)
            else:
                normalized[:, j] = 0.5
        # Step 2: Contrast intensity (standard deviation)
        std_devs = np.std(normalized, axis=0)

        # Step 3: Conflict (1 - |correlation|)
        if normalized.shape[0] < 3:
            correlation_matrix = np.eye(len(objectives))
        else:
            correlation_matrix = np.corrcoef(normalized, rowvar=False)
            correlation_matrix = np.nan_to_num(correlation_matrix, nan=0.0)

        conflict_matrix = 1 - np.abs(correlation_matrix)
        conflict_sum = conflict_matrix.sum(axis=1)

        # Step 4: Information content
        information_content = std_devs * conflict_sum

        # Step 5: Normalize to weights
        total_info = information_content.sum()
        if total_info < 1e-10:
            return {obj: 1.0 / len(objectives) for obj in objectives}

        weights = information_content / total_info
        return dict(zip(objectives, weights))

    # ------------------------------------------------------------------
    def _apply_iqr_normalization(self, series: pd.Series) -> pd.Series:
        """IQR-based robust normalization with 1.5×IQR outlier clipping."""
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        if iqr < 1e-10:
            return pd.Series([0.5] * len(series), index=series.index)
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        clipped = series.clip(lower_bound, upper_bound)
        normalized = ((clipped - lower_bound) /
                      (upper_bound - lower_bound + 1e-10))
        return normalized.clip(0, 1)

    # ==================================================================
    # [FIX] Multi-objective selection — ALL from Stage 1 (no prior)
    # ==================================================================
    def execute_multi_objective_selection(
            self,
            method: str = 'balanced',
            weights: Optional[Dict[str, float]] = None,
            use_critic: bool = True) -> List[int]:
        """
        Select optimal k in each partition with constraints:
        - No duplicate k-values
        - Adjacent k-values gap ≤ 28
        - Prioritize CRITIC score

        [FIX] All CRITIC features come from self._stage1_metrics
        (frozen, no cross-k prior contamination).

        Objectives:
            info — PESC soft effective ratio
            snr  — signal-to-noise ratio
            unif — peak sharpness
            conf — boundary confidence
        """
        if weights is None:
            weights = {
                'node_retention': 0.25, 'error_mass': 0.25,
                'species_resolution': 0.25, 'cutoff_sharpness': 0.25
            }

        self.selection_config = {
            'method': method,
            'use_critic': use_critic,
            'manual_weights': weights.copy() if not use_critic else None,
            'objectives': ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']
        }

        self.critic_weights_history = []
        
        # 为每个分区生成候选k值
        partition_candidates = []
        for part_idx, (start_idx, end_idx) in \
                enumerate(self.adaptive_partitions):

            feasible_candidates = []
            for i in range(start_idx, end_idx + 1):
                # [FIX] Read from frozen Stage 1 — no prior contamination
                metrics = self._stage1_metrics[i]
                if not metrics or not metrics.get('em_converged'):
                    continue

                base_value = ((abs(metrics['mu'] - metrics['mu_rare']) / (metrics['sigma'] + metrics['sigma_rare'] + 1e-6) * min(metrics['mu'] / max(metrics.get('mu_rare', 0.1*metrics['mu']), 1e-6), 3.0) / max(metrics.get('coverage_heterogeneity', 1.0), 1.0)))
                species_resolution_value = base_value if metrics.get('n_signal_components', 1) == 2 else 0.5

                candidate = {
                    'k': metrics['k'],
                    # 4 CRITIC objectives (all from Stage 1)
                    # 1. 真实节点保留率（加入重复序列惩罚，安全访问）
                    'node_retention': metrics['soft_effective_diff'] / metrics.get('distinct_kmers', 1) ,
                    # 2. 错误与污染惩罚 (移除 w_err，仅用 contamination_proxy + p_err)
                    'error_mass': 1.0 / (1.0 + 0.3*metrics.get('contamination_proxy', 0) +
                                        0.2 * (1.0 / max(metrics.get('p_err', 0.3), 0.05))),
                    # 3. 物种解离度 (移除无效 conserved_kmer_ratio，加入异质性惩罚)
                    'species_resolution': species_resolution_value,
                    # 4. 截断锐度（保留原逻辑）
                    'cutoff_sharpness': np.exp(-metrics['boundary_uncertainty']),
                    # Auxiliary (not in CRITIC)
                    'effective_diff': metrics['effective_diff'],
                    'soft_effective_diff': metrics['soft_effective_diff'],
                    'sigma': metrics['sigma'],
                    'mu': metrics['mu'],
                    'uncertainty': metrics['boundary_uncertainty']
                }
                feasible_candidates.append(candidate)

            if not feasible_candidates:
                logger.warning(
                    f"Partition {part_idx}: no feasible candidates, using default k")
                # 使用分区中间的k值作为备选
                default_k = self.k_values[(start_idx + end_idx) // 2]
                partition_candidates.append([default_k])
                continue

            candidates_df = pd.DataFrame(feasible_candidates)

            # CRITIC automatic weighting
            if use_critic and len(candidates_df) >= 3:
                current_weights = self._compute_critic_weights(candidates_df)
                self.critic_weights_history.append({
                    'partition': part_idx,
                    'weights': current_weights,
                    'method': 'CRITIC'
                })
                logger.debug(
                    f"Partition {part_idx} CRITIC weights: "
                    f"{current_weights}")
            else:
                current_weights = weights

            # 计算每个候选的得分
            for obj in ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']:
                # 对可能大范围的指标做log压缩，避免异常值主导
                if obj == 'species_resolution':
                    candidates_df[f'{obj}_log'] = np.log1p(candidates_df[obj])
                    candidates_df[f'{obj}_normalized'] = self._apply_iqr_normalization(candidates_df[f'{obj}_log'])
                else:
                    candidates_df[f'{obj}_normalized'] = self._apply_iqr_normalization(candidates_df[obj])

            candidates_df['final_score'] = (
                current_weights['node_retention'] * candidates_df['node_retention_normalized'] +
                current_weights['error_mass']  * candidates_df['error_mass_normalized'] +
                current_weights['species_resolution'] * candidates_df['species_resolution_normalized'] +
                current_weights['cutoff_sharpness'] * candidates_df['cutoff_sharpness_normalized']
            )

            # 按得分排序，获取前3个候选
            sorted_candidates = candidates_df.sort_values('final_score', ascending=False)
            top_candidates = sorted_candidates['k'].tolist()[:3]  # 最优、次优、三级备选
            partition_candidates.append(top_candidates)

        # 按照用户要求的逻辑：先收集所有最优k值，再检查重复和间隔
        # 步骤1：收集每个分区的候选k值（最优、次优、第三优）
        all_candidates = []
        for i, candidates in enumerate(partition_candidates):
            if not candidates:
                logger.warning(f"Partition {i}: no feasible candidates")
                all_candidates.append([])
            else:
                logger.info(f"Partition {i}: candidates={candidates}")
                all_candidates.append(candidates)

        # 步骤2：去重 — 先解决重复k（比较次优替换），再解决间隔
        # 2a. 构建无重复的 optimal_ks
        partition_choices = []
        used_k = set()

        for part_idx, cand_list in enumerate(all_candidates):
            if not cand_list:
                continue

            best_k = cand_list[0]
            if best_k not in used_k:
                partition_choices.append((best_k, part_idx))
                used_k.add(best_k)
            else:
                conflict_part_idx = None
                for _, (prev_k, prev_part) in enumerate(partition_choices):
                    if prev_k == best_k:
                        conflict_part_idx = prev_part
                        break

                if conflict_part_idx is None:
                    partition_choices.append((best_k, part_idx))
                    used_k.add(best_k)
                    continue

                list_a = all_candidates[conflict_part_idx]
                list_b = cand_list

                assigned_b = None
                for candidate in list_b:
                    if candidate not in used_k:
                        assigned_b = candidate
                        break

                if assigned_b is not None:
                    partition_choices.append((assigned_b, part_idx))
                    used_k.add(assigned_b)
                else:
                    assigned_a = None
                    for candidate in list_a:
                        if candidate not in used_k and candidate != best_k:
                            assigned_a = candidate
                            break

                    if assigned_a is not None:
                        partition_choices = [
                            (k, p) for k, p in partition_choices
                            if p != conflict_part_idx
                        ]
                        partition_choices.append((assigned_a, conflict_part_idx))
                        used_k.discard(best_k)
                        used_k.add(assigned_a)
                        partition_choices.append((best_k, part_idx))
                        used_k.add(best_k)
                    else:
                        mid_k = (list_a[0] + list_b[0]) // 2
                        closest_k = min(
                            self.k_values, key=lambda x: abs(x - mid_k))
                        if closest_k not in used_k:
                            partition_choices.append((closest_k, part_idx))
                            used_k.add(closest_k)
                        else:
                            for kv in self.k_values:
                                if kv not in used_k:
                                    partition_choices.append((kv, part_idx))
                                    used_k.add(kv)
                                    break
                        logger.warning(
                            f"Partition {part_idx} vs {conflict_part_idx}: "
                            f"all candidates conflict, using fallback k")

        optimal_ks = [
            k for k, _ in sorted(partition_choices, key=lambda x: x[1])]
        logger.info(f"After dedup: optimal_ks = {optimal_ks}")

        # 2b. 检查相邻最优k值的间隔，添加过渡k值（每次append都去重）
        fixed_ks = [optimal_ks[0]]
        max_level = 3

        for i in range(1, len(optimal_ks)):
            prev_k = fixed_ks[-1]
            curr_k = optimal_ks[i]
            gap = curr_k - prev_k

            if gap <= 28:
                if curr_k not in fixed_ks:
                    fixed_ks.append(curr_k)
            else:
                logger.info(
                    f"Gap {gap} > 28 between {prev_k} and {curr_k}, "
                    f"trying to add transition k")

                transition_added = False

                for level in range(1, max_level):
                    if i > 0 and len(all_candidates[i - 1]) > level:
                        prev_candidate = all_candidates[i - 1][level]
                        if (prev_k < prev_candidate < curr_k
                                and prev_candidate not in fixed_ks):
                            if (abs(prev_k - prev_candidate) <= 28
                                    and abs(prev_candidate - curr_k) <= 28):
                                fixed_ks.append(prev_candidate)
                                if curr_k not in fixed_ks:
                                    fixed_ks.append(curr_k)
                                logger.info(
                                    f"Added transition k={prev_candidate} "
                                    f"from previous partition level {level}")
                                transition_added = True
                                break

                    if len(all_candidates[i]) > level:
                        curr_candidate = all_candidates[i][level]
                        if (prev_k < curr_candidate < curr_k
                                and curr_candidate not in fixed_ks):
                            if (abs(prev_k - curr_candidate) <= 28
                                    and abs(curr_candidate - curr_k) <= 28):
                                fixed_ks.append(curr_candidate)
                                if curr_k not in fixed_ks:
                                    fixed_ks.append(curr_k)
                                logger.info(
                                    f"Added transition k={curr_candidate} "
                                    f"from current partition level {level}")
                                transition_added = True
                                break

                if not transition_added:
                    if (i > 0
                            and len(all_candidates[i - 1]) > 1
                            and len(all_candidates[i]) > 1):
                        prev_candidate = all_candidates[i - 1][1]
                        curr_candidate = all_candidates[i][1]
                        if (prev_k < prev_candidate < curr_candidate < curr_k
                                and prev_candidate not in fixed_ks
                                and curr_candidate not in fixed_ks):
                            if (abs(prev_k - prev_candidate) <= 28
                                    and abs(prev_candidate - curr_candidate) <= 28
                                    and abs(curr_candidate - curr_k) <= 28):
                                fixed_ks.append(prev_candidate)
                                fixed_ks.append(curr_candidate)
                                if curr_k not in fixed_ks:
                                    fixed_ks.append(curr_k)
                                logger.info(
                                    f"Added transition ks="
                                    f"{prev_candidate}, {curr_candidate} "
                                    f"from both partitions")
                                transition_added = True

                if not transition_added:
                    if i > 0 and len(all_candidates[i - 1]) > 2:
                        prev_candidate = all_candidates[i - 1][2]
                        if (prev_k < prev_candidate < curr_k
                                and prev_candidate not in fixed_ks):
                            if (abs(prev_k - prev_candidate) <= 28
                                    and abs(prev_candidate - curr_k) <= 28):
                                fixed_ks.append(prev_candidate)
                                if curr_k not in fixed_ks:
                                    fixed_ks.append(curr_k)
                                logger.info(
                                    f"Added transition k={prev_candidate} "
                                    f"from previous partition level 2")
                                transition_added = True

                    if (not transition_added
                            and len(all_candidates[i]) > 2):
                        curr_candidate = all_candidates[i][2]
                        if (prev_k < curr_candidate < curr_k
                                and curr_candidate not in fixed_ks):
                            if (abs(prev_k - curr_candidate) <= 28
                                    and abs(curr_candidate - curr_k) <= 28):
                                fixed_ks.append(curr_candidate)
                                if curr_k not in fixed_ks:
                                    fixed_ks.append(curr_k)
                                logger.info(
                                    f"Added transition k={curr_candidate} "
                                    f"from current partition level 2")
                                transition_added = True

                if not transition_added:
                    mid_point = (prev_k + curr_k) // 2
                    closest_k = min(
                        self.k_values, key=lambda x: abs(x - mid_point))
                    if closest_k not in fixed_ks:
                        fixed_ks.append(closest_k)
                    if curr_k not in fixed_ks:
                        fixed_ks.append(curr_k)
                    logger.warning(
                        f"All methods failed, using closest k={closest_k} "
                        f"to midpoint {mid_point} as fallback")

        selected_ks = fixed_ks

        self.selected_kmers = selected_ks
        logger.info(
            f"🎯 Selected {len(self.selected_kmers)} optimal k-values: "
            f"{self.selected_kmers}")

        if self.critic_weights_history:
            avg_weights = {
                obj: np.mean([w['weights'][obj]
                              for w in self.critic_weights_history])
                for obj in ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']
            }
            logger.info(f"📊 Average CRITIC weights: {avg_weights}")

        return self.selected_kmers

    # ==================================================================
    # [ABLATION] Global CRITIC Selection (No Partitioning)
    # ==================================================================
    def execute_global_critic_selection(self, n_top_k: int = 5, min_gap: int = 10) -> List[int]:
        """
        Ablation mode: Global CRITIC selection without adaptive partitioning.
        Selects top N k-values based on global 4-objective CRITIC scores.
        
        Parameters
        ----------
        n_top_k : int
            Number of k-values to select globally.
        min_gap : int
            Minimum gap between selected k-values to prevent clustering 
            (e.g., selecting 31, 33, 35, 37, 39).
        """
        logger.info(f"🧪 [ABLATION] Running Global CRITIC Selection (N={n_top_k}, min_gap={min_gap})...")
        
        # 1. 收集所有收敛的 metrics 并构建全局候选 DataFrame
        feasible_candidates = []
        for metrics in self._stage1_metrics:
            if not metrics or not metrics.get('em_converged'):
                continue
            
            # 计算 4 个 CRITIC 指标 (公式与分区模式完全一致，保证控制变量)
            base_value = (
                abs(metrics['mu'] - metrics['mu_rare']) / 
                (metrics['sigma'] + metrics['sigma_rare'] + 1e-6) * 
                min(metrics['mu'] / max(metrics.get('mu_rare', 0.1*metrics['mu']), 1e-6), 3.0) / 
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
            logger.error("No feasible candidates for global selection.")
            return []
            
        candidates_df = pd.DataFrame(feasible_candidates)
        
        # 2. 计算全局 CRITIC 权重
        global_weights = self._compute_critic_weights(candidates_df)
        logger.info(f"🌍 Global CRITIC Weights: {global_weights}")
        
        # 3. 归一化并计算全局综合得分
        for obj in ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']:
            if obj == 'species_resolution':
                candidates_df[f'{obj}_log'] = np.log1p(candidates_df[obj])
                candidates_df[f'{obj}_normalized'] = self._apply_iqr_normalization(candidates_df[f'{obj}_log'])
            else:
                candidates_df[f'{obj}_normalized'] = self._apply_iqr_normalization(candidates_df[obj])
                
        candidates_df['final_score'] = (
            global_weights['node_retention'] * candidates_df['node_retention_normalized'] +
            global_weights['error_mass']  * candidates_df['error_mass_normalized'] +
            global_weights['species_resolution'] * candidates_df['species_resolution_normalized'] +
            global_weights['cutoff_sharpness'] * candidates_df['cutoff_sharpness_normalized']
        )
        
        # 4. 贪心选择 Top N (带 min_gap 约束，防止 k 值过度聚集)
        sorted_df = candidates_df.sort_values('final_score', ascending=False).reset_index(drop=True)
        
        selected_ks = []
        for _, row in sorted_df.iterrows():
            k_val = int(row['k'])
            if not selected_ks:
                selected_ks.append(k_val)
            else:
                # 检查与已选 k 值的最小间隔
                if all(abs(k_val - sel_k) >= min_gap for sel_k in selected_ks):
                    selected_ks.append(k_val)
            if len(selected_ks) >= n_top_k:
                break
                
        # 如果因为 min_gap 限制导致选不够 n_top_k，则放宽限制补齐
        if len(selected_ks) < n_top_k:
            logger.warning(f"Could only select {len(selected_ks)} k-values with min_gap={min_gap}. Relaxing gap constraint...")
            for _, row in sorted_df.iterrows():
                k_val = int(row['k'])
                if k_val not in selected_ks:
                    selected_ks.append(k_val)
                if len(selected_ks) >= n_top_k:
                    break
                    
        selected_ks.sort()
        self.selected_kmers = selected_ks
        logger.info(f"🎯 [ABLATION] Globally selected {len(self.selected_kmers)} k-values: {self.selected_kmers}")
        return self.selected_kmers
    # ------------------------------------------------------------------
    def _select_max_effective(
            self, df: pd.DataFrame, weights: Dict[str, float]) -> int:
        """Fallback: select k with max soft_effective_diff."""
        return int(df.loc[df['soft_effective_diff'].idxmax(), 'k'])

    def _select_balanced_pareto(
            self, df: pd.DataFrame, weights: Dict[str, float]) -> int:
        """Weighted 4-objective Pareto selection with IQR normalization."""
        for obj in ['node_retention', 'error_mass', 'species_resolution', 'cutoff_sharpness']:
            # 对可能大范围的指标做log压缩，避免异常值主导
            if obj == 'species_resolution':
                df[f'{obj}_log'] = np.log1p(df[obj])
                df[f'{obj}_normalized'] = self._apply_iqr_normalization(df[f'{obj}_log'])
            else:
                df[f'{obj}_normalized'] = self._apply_iqr_normalization(df[obj])

        df['final_score'] = (
            weights['node_retention'] * df['node_retention_normalized'] +
            weights['error_mass']  * df['error_mass_normalized'] +
            weights['species_resolution'] * df['species_resolution_normalized'] +
            weights['cutoff_sharpness'] * df['cutoff_sharpness_normalized']
        )
        return int(df.loc[df['final_score'].idxmax(), 'k'])

    # ==================================================================
    # Two-stage EM pipeline
    # ==================================================================
    def fit_two_stage_em(self, spectrum_dir: str, init_method: str = 'robust', n_signal_components: int = 2) -> bool:
        """ 
        Single-Stage Pipeline (Stage 2 removed): 
        Stage 1: Coarse EM (NO cross-k prior) → freeze into self._stage1_metrics
        Stage 1.5: Adaptive partitioning from Stage 1 profiles
        """
        logger.info("🔄 Starting EM fitting and partitioning...")
        
        # === STAGE 1: coarse (no prior) ===
        logger.info("[Stage 1/2] Coarse EM — maximizing fidelity...")
        ok = self.fit_em_frequencies(spectrum_dir, init_method=init_method, n_signal_components=n_signal_components)
        if not ok:
            logger.error("Stage 1 failed")
            return False
            
        # === STAGE 1.5: partition from Stage 1 ===
        logger.info("[Stage 2/2] Adaptive partitioning from Stage 1 profiles...")
        self.construct_adaptive_partitions(min_size=4)
        
        return ok

    # ==================================================================
    # [FIX] Diagnostic visualization
    # ==================================================================
    def diagnose_partitioning(self,
                              save_path: str = 'partitioning_diagnosis.png'):
        """
        Generate 3-panel diagnostic figure:
          1. Stage 1 soft_effective_diff + partition boundaries
          2. First derivative (rate of change)
          3. Selected k-values + their metrics
        """
        if not self._stage1_profiles:
            logger.error("No Stage 1 profiles. Run fit_two_stage_em first.")
            return

        fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

        valid_mask = ~np.isnan(np.array(self._stage1_profiles))
        valid_ks = np.array(self.k_values)[valid_mask]
        valid_profiles = np.array(self._stage1_profiles)[valid_mask]

        # --- Panel 1: raw profile + partitions ---
        axes[0].plot(valid_ks, valid_profiles / 1e6, 'b-o',
                     markersize=3, label='Stage 1 soft_effective_diff')

        colors = plt.cm.Set2(
            np.linspace(0, 1, max(len(self.adaptive_partitions), 1)))
        for pi, (s, e) in enumerate(self.adaptive_partitions):
            ks = self.k_values[s:e + 1]
            axes[0].axvspan(
                ks[0] - 1, ks[-1] + 1, alpha=0.15,
                color=colors[pi % len(colors)],
                label=f'P{pi}: k={ks[0]}-{ks[-1]}')

        axes[0].set_ylabel('soft_effective_diff (M)')
        axes[0].set_title('Stage 1 Profile & Adaptive Partitions')
        axes[0].legend(fontsize=7, loc='upper right')
        axes[0].grid(True, alpha=0.3)

        # --- Panel 2: derivative ---
        dx = np.diff(valid_profiles)
        abs_dx = np.abs(dx)
        dx_ks = valid_ks[:-1]

        axes[1].plot(dx_ks, abs_dx / 1e6, 'r-', alpha=0.5,
                     label='Raw |d/dk|')

        win = min(7, max(3, len(abs_dx) // 6))
        if win % 2 == 0:
            win += 1
        if len(abs_dx) >= win and win >= 5:
            smoothed_dx = savgol_filter(abs_dx, window_length=win,
                                        polyorder=2)
            axes[1].plot(dx_ks, smoothed_dx / 1e6, 'k-', linewidth=2,
                         label=f'Smoothed |d/dk| (win={win})')

        axes[1].set_ylabel('|Δsoft_diff| (M)')
        axes[1].set_title('First Derivative (Rate of Change)')
        axes[1].legend(fontsize=8)
        axes[1].grid(True, alpha=0.3)

        # --- Panel 3: selected k-values ---
        if self.selected_kmers:
            sel_profiles = []
            sel_uncertainties = []
            for sk in self.selected_kmers:
                idx = (self.k_values.index(sk)
                       if sk in self.k_values else -1)
                if idx >= 0 and self._stage1_metrics[idx]:
                    sel_profiles.append(
                        self._stage1_metrics[idx]['soft_effective_diff']
                        / 1e6)
                    sel_uncertainties.append(
                        self._stage1_metrics[idx].get(
                            'boundary_uncertainty', 0))
                else:
                    sel_profiles.append(0)
                    sel_uncertainties.append(0)

            # 使用实际的 k 值作为 x 轴，保持与上两图对齐
            bar_width = 1.5 * (self.k_values[1] - self.k_values[0]) if len(self.k_values) > 1 else 1.5
            
            axes[2].bar(self.selected_kmers, sel_profiles, width=bar_width,
                        color='steelblue', alpha=0.7,
                        label='soft_effective_diff (M)')
            
            # 设置 x 轴刻度，只显示选中的 k 值
            axes[2].set_xticks(self.selected_kmers)
            axes[2].set_xticklabels(
                [f'k={k}' for k in self.selected_kmers],
                rotation=45, fontsize=9)
            
            axes[2].set_ylabel('soft_effective_diff (M)')
            axes[2].set_title('Selected k-values & Metrics')
            axes[2].legend(fontsize=8, loc='upper left')

            ax2 = axes[2].twinx()
            ax2.bar(self.selected_kmers, sel_uncertainties, width=bar_width*0.6,
                    color='orange', alpha=0.5, label='uncertainty')
            ax2.set_ylabel('Boundary Uncertainty')
            ax2.legend(fontsize=8, loc='upper right')
            
            # 统一设置 x 轴范围，使三个图对齐
            if len(valid_ks) > 0:
                axes[2].set_xlim(valid_ks[0] - 2, valid_ks[-1] + 2)

        axes[2].grid(True, axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"📊 Diagnostic plot saved to {save_path}")
        
    def _enforce_megahit_gap_limit(self, selected_ks: List[int], max_gap: int = 28) -> List[int]:
        """
        【智能版】插入过渡k时，选择区间内最优k，而非简单中点
        """
        import math
        
        if len(selected_ks) < 2:
            return selected_ks
        
        fixed_ks = [selected_ks[0]]
        
        for i in range(1, len(selected_ks)):
            prev_k = fixed_ks[-1]
            curr_k = selected_ks[i]
            gap = curr_k - prev_k
            
            if gap > max_gap:
                # 需要插入过渡k
                n_segments = math.ceil(gap / max_gap)
                
                for j in range(1, n_segments):
                    ideal_k = prev_k + j * (gap / n_segments)
                    
                    # 【关键改进】在候选k中找CRITIC综合得分最高的
                    valid_candidates = [
                        k for k in self.k_values 
                        if prev_k < k < curr_k
                    ]
                    
                    if valid_candidates:
                        # 找到这些k对应的metrics
                        best_k = None
                        best_score = -1
                        
                        for candidate_k in valid_candidates:
                            idx = self.k_values.index(candidate_k)
                            metrics = self._stage1_metrics[idx]
                            
                            if metrics and metrics.get('em_converged'):
                                # 计算CRITIC综合得分（与主流程新指标完全一致）
                                score = (
                                    metrics['soft_effective_diff'] / max(metrics.get('distinct_kmers', 1) + 0.1*metrics.get('error_sum', 0), 1) * 0.15 +
                                    1.0 / (1.0 + 0.3*metrics.get('contamination_proxy', 0) + 0.2 * metrics['w_err'] * (1.0 / max(metrics.get('p_err', 0.3), 0.05))) * 0.15 +
                                    ((abs(metrics['mu'] - metrics['mu_rare']) / (metrics['sigma'] + metrics['sigma_rare'] + 1e-6) * min(metrics['mu'] / max(metrics.get('mu_rare', 0.1*metrics['mu']), 1e-6), 3.0)) / max(metrics.get('coverage_heterogeneity', 1.0), 1.0) if metrics.get('n_signal_components', 1) == 2 else 0.5) * 0.5 +
                                    np.exp(-metrics['boundary_uncertainty']) * 0.2
                                )
                                
                                # 优先选择接近ideal_k的
                                distance_penalty = abs(candidate_k - ideal_k) / gap
                                adjusted_score = score * (1 - 0.3 * distance_penalty)
                                
                                if adjusted_score > best_score:
                                    best_score = adjusted_score
                                    best_k = candidate_k
                        
                        if best_k and best_k != fixed_ks[-1]:
                            fixed_ks.append(best_k)
                            logger.info(
                                f"  插入过渡k={best_k} (理想点={ideal_k:.1f}, "
                                f"CRITIC得分={best_score:.3f})"
                            )
            
            # 无论gap是否大于max_gap，都添加当前k值
            fixed_ks.append(curr_k)
        
        return fixed_ks

   
   
     
    # ==================================================================
    # Full pipeline runner
    # ==================================================================
    # def run(self, spectrum_dir: str, init_method: str = 'robust', 
    #         selection_method: str = 'balanced', use_critic: bool = True, 
    #         diagnose: bool = True, diagnose_path: str = 'partitioning_diagnosis.png',
    #         n_signal_components: int = 2) -> List[int]:
    #     """ 
    #     Complete pipeline: fit → partition → select → (optional) diagnose. 
    #     [FIX] 移除了 alpha 参数，新增 n_signal_components 强制双峰
    #     """
    #     # Stage 1 + 1.5 (已移除 Stage 2)
    #     ok = self.fit_two_stage_em(
    #         spectrum_dir, init_method=init_method, n_signal_components=n_signal_components)
            
    #     if not ok:
    #         logger.error("Pipeline failed at EM fitting stage")
    #         return []
            
    #     # CRITIC selection (from Stage 1) - 已集成间隔约束检查
    #     selected = self.execute_multi_objective_selection(
    #         method=selection_method, use_critic=use_critic)
    #     # 同步更新，确保诊断图画的是修正后的
    #     self.selected_kmers = selected           
    #     # Diagnostic
    #     if diagnose:
    #         self.diagnose_partitioning(save_path=diagnose_path)
            
    #     return selected
    def run(self, spectrum_dir: str, init_method: str = 'robust',
            selection_method: str = 'balanced', use_critic: bool = True,
            diagnose: bool = True, diagnose_path: str = 'partitioning_diagnosis.png',
            n_signal_components: int = 2,
            ablation_mode: bool = False, n_top_k: int = 5) -> List[int]:
        """
        Complete pipeline: fit → partition → select → (optional) diagnose.
        
        Parameters
        ----------
        ablation_mode : bool
            If True, skips adaptive partitioning and uses global CRITIC selection.
        n_top_k : int
            Number of k-values to select in ablation_mode.
        """
        # Stage 1 + 1.5 (已移除 Stage 2)
        ok = self.fit_two_stage_em(
            spectrum_dir, init_method=init_method, n_signal_components=n_signal_components)
        if not ok:
            logger.error("Pipeline failed at EM fitting stage")
            return []
            
        # ==================================================================
        # 分支：消融模式 vs 正常模式
        # ==================================================================
        if ablation_mode:
            logger.info("🧪 ABLATION MODE ACTIVATED: Bypassing partitioning...")
            selected = self.execute_global_critic_selection(n_top_k=n_top_k, min_gap=10)
        else:
            # 正常模式：先分区，再局部 CRITIC 选优
            selected = self.execute_multi_objective_selection(
                method=selection_method, use_critic=use_critic)
                
        # 同步更新，确保诊断图画的是修正后的
        self.selected_kmers = selected
        
        # Diagnostic
        if diagnose:
            # 如果是消融模式，诊断图可能不需要画分区，但画出来看看全局分布也好
            self.diagnose_partitioning(save_path=diagnose_path)
            
        return selected
    # ------------------------------------------------------------------
    # Save results to files
    # ------------------------------------------------------------------
    def save_results(self, prefix: str = "kmer_results") -> None:
        """
        Save all analysis results to files for reproducibility.
        
        Parameters
        ----------
        prefix : str
            Output file prefix (default: "kmer_results")
        """
        import json
        import pandas as pd
        
        logger.info(f"💾 Saving results with prefix: {prefix}")
        
        # 1. Save selected k-values
        with open(f"{prefix}_selected_ks.txt", 'w') as f:
            f.write(f"# Selected k-mer lengths (CRITIC-based multi-objective optimization)\n")
            f.write(f"# Generated: {pd.Timestamp.now().isoformat()}\n")
            f.write(f"K_LIST={self.selected_kmers}\n")
        logger.info(f"   ✅ {prefix}_selected_ks.txt")
        
        # 2. Save comprehensive summary JSON
        summary = {
            'metadata': {
                'k_range': [self.k_values[0], self.k_values[-1]] if self.k_values else [],
                'step': self.k_values[1] - self.k_values[0] if len(self.k_values) > 1 else 1,
                'generated_at': pd.Timestamp.now().isoformat(),
                'read_length': self.read_length
            },
            'partitioning': {
                'n_partitions': len(self.adaptive_partitions),
                'config': self.partition_config,
                'partitions_k': [
                    (self.k_values[s], self.k_values[e])
                    for s, e in self.adaptive_partitions
                ] if self.adaptive_partitions else []
            },
            'selection': {
                'config': self.selection_config,
                'critic_weights_history': self.critic_weights_history,
                'selected_ks': self.selected_kmers
            },
            'convergence': {
                'total_k_evaluated': len(self.k_values),
                'converged_count': sum(
                    1 for m in self.metrics_list if m and m.get('em_converged')
                ) if self.metrics_list else 0
            }
        }
        
        with open(f"{prefix}_summary.json", 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        logger.info(f"   ✅ {prefix}_summary.json")
        
        # 3. Export metrics table (CSV)
        rows = []
        for i, (k, metrics) in enumerate(zip(self.k_values, self.metrics_list)):
            if not metrics or not metrics.get('em_converged'):
                continue
            
            # Determine partition membership
            partition_id = None
            for j, (start, end) in enumerate(self.adaptive_partitions):
                if start <= i <= end:
                    partition_id = j
                    break
            
            species_resolution_value = ((abs(metrics.get('mu', 0) - metrics.get('mu_rare', 0)) / (metrics.get('sigma', 1) + metrics.get('sigma_rare', 1) + 1e-6) * min(metrics.get('mu', 1) / max(metrics.get('mu_rare', 0.1*metrics.get('mu', 1)), 1e-6), 3.0) / max(metrics.get('coverage_heterogeneity', 1.0), 1.0)))
            species_resolution_value = species_resolution_value if metrics.get('n_signal_components', 1) == 2 else 0.5
            
            row = {
                'k': k,
                'partition_id': partition_id,
                'effective_diff': metrics.get('effective_diff'),
                'soft_effective_diff': metrics.get('soft_effective_diff'),
                'distinct_kmers': metrics.get('distinct_kmers'),
                'w_sig': metrics.get('w_sig'),
                'w_err': metrics.get('w_err'),
                'p_err': metrics.get('p_err'),
                'sigma': metrics.get('sigma'),
                'sigma_rare': metrics.get('sigma_rare'),
                'mu': metrics.get('mu'),
                'mu_rare': metrics.get('mu_rare'),
                'cutoff': metrics.get('cutoff'),
                'boundary_uncertainty': metrics.get('boundary_uncertainty'),
                'n_signal_components': metrics.get('n_signal_components', 1),
                'contamination_proxy': metrics.get('contamination_proxy', 0),
                'coverage_heterogeneity': metrics.get('coverage_heterogeneity', 1),
                'node_retention': metrics.get('soft_effective_diff', 0) / 
                                 max(metrics.get('distinct_kmers', 1) + 0.1*metrics.get('error_sum', 0), 1),
                'error_mass': 1.0 / (1.0 + 0.3*metrics.get('contamination_proxy', 0) +
                                      0.2 * (1.0 / max(metrics.get('p_err', 0.3), 0.05))),
                'species_resolution': species_resolution_value,
                'cutoff_sharpness': np.exp(-metrics.get('boundary_uncertainty', 0)),
                'is_selected': k in self.selected_kmers
            }
            rows.append(row)
        
        if rows:
            pd.DataFrame(rows).to_csv(f"{prefix}_metrics.csv", index=False)
            logger.info(f"   ✅ {prefix}_metrics.csv ({len(rows)} rows)")
        
        logger.info(f"💾 All results saved with prefix '{prefix}'")