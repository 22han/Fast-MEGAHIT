#!/usr/bin/env python3
"""
SparseKmerSpectrumAnalyzer — Single K-mer Spectrum Mixture Model Fitter

Two-component mixture model:
  - Error component:   Geometric distribution (discrete exponential decay)
  - Genuine component: Normal distribution (coverage peak)

Innovations:
  1. PESC (Posterior-Expected Soft-Counting): continuous soft boundary
     that eliminates step-function artifacts in effective k-mer counting.
  2. CRITIC-cutoff: gradient-based adaptive threshold that replaces
     the fixed posterior_threshold=0.5 with data-driven contrast analysis.
  3. Boundary uncertainty: Shannon entropy of ambiguous k-mers near the
     cutoff, providing a natural confidence metric for downstream CRITIC.

References:
  Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). JRSS-B, 39(1), 1-22.
  Diakoulaki, D., Mavrotas, G., & Papayannakis, L. (1995). COR, 22(7), 763-770.
  Savitzky, A., & Golay, M. J. E. (1964). Anal. Chem., 36(8), 1627-1639.
"""

import sys
import numpy as np
from scipy.stats import geom, norm

class SparseKmerSpectrumAnalyzer:
    """
    Memory-efficient k-mer spectrum analyzer using sparse representation.

    Stores only non-zero frequency counts as {frequency: count} dict,
    suitable for large spectrum files (KmerLight / ntCard / KMC format).
    """

    def __init__(self, spectrum_file=None, spectrum_dict=None,
                 max_freq=None, posterior_threshold=0.5,
                 max_em_iter=200, tol=1e-4):
        """
        Parameters
        ----------
        spectrum_file : str, optional
            Path to KmerLight/ntCard/KMC spectrum file.
        spectrum_dict : dict, optional
            Pre-built {frequency: count} dictionary.
        max_freq : int, optional
            Maximum frequency to consider.
        posterior_threshold : float
            Fallback threshold for hard cutoff (default 0.5).
        max_em_iter : int
            Maximum EM iterations (default 200).
        tol : float
            Convergence tolerance (default 1e-4).
        """
        self.frequencies = {}  # 使用稀疏字典而非稠密数组，只存储有观测到的频率。
        self.distinct_kmers = None 
        self.posterior_threshold = posterior_threshold
        self.max_em_iter = max_em_iter
        self.tol = tol

        # Load data
        if spectrum_file is not None:
            self._parse_file(spectrum_file, max_freq)
        elif spectrum_dict is not None:
            self.frequencies = {int(k): float(v) for k, v in spectrum_dict.items()}
            if max_freq is not None:
                self.frequencies = {k: v for k, v in self.frequencies.items()
                                   if k <= max_freq}
        
        self.max_freq = max(self.frequencies.keys()) if self.frequencies else 0

        # Results placeholders
        self.params = None  #EM拟合结果
        self.cutoff = None   #自适应截断频率
        self.effective_diff = None  #硬计数有效k-mer
        self.error_sum = None
        self.converged = False 

        # [FIX] PESC and uncertainty results
        self.soft_effective_diff = None   #PESC软计数
        self.boundary_uncertainty = None   #边界熵
        self.mu_rare = None
        self.sigma_rare = None
        self.w_rare_frac = 0.3  # Default, will be updated in fit()

    # ------------------------------------------------------------------
    # File parsing
    # ------------------------------------------------------------------
    def _parse_file(self, filepath, max_freq=None):
        """
        Parse spectrum file in streaming fashion (memory efficient).

        Supported formats:
          KMC (2-col):  "freq  count"
          KMC (3-col):  "line  freq  count"  (first col ignored)
          KmerLight:    "f1  12345"
        """
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) < 2:
                        continue

                    if parts[0].startswith('f'):
                        # KmerLight / ntCard
                        try:
                            freq = int(parts[0][1:])
                            count = float(parts[1])
                        except ValueError:
                            continue
                    else:
                        if len(parts) == 2:
                            freq = int(parts[0])
                            count = float(parts[1])
                        elif len(parts) >= 3:
                            freq = int(parts[-2])
                            count = float(parts[-1])
                        else:
                            continue

                    if max_freq is not None and freq > max_freq:
                        continue
                    if count > 0:
                        self.frequencies[freq] = count

            if not self.frequencies:
                raise ValueError(f"No frequency data found in {filepath}")

            self.distinct_kmers = sum(self.frequencies.values())

        except Exception as e:
            raise Exception(f"Error parsing {filepath}: {e}")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _get_arrays(self):
        """Convert sparse dict to sorted numpy arrays."""
        if not self.frequencies:
            raise ValueError("No frequency data available")
        freq_vals = np.array(sorted(self.frequencies.keys()), dtype=np.float64)
        counts = np.array([self.frequencies[f] for f in freq_vals], dtype=np.float64)
        return freq_vals, counts #返回两个频率数组

    def _initialize_params(self):
        """Simple heuristic initialization."""
        freq_vals, counts = self._get_arrays()
        
        total_obs = np.sum(counts)
        if total_obs == 0:
            raise ValueError("No observed k-mers")

        w_err = 0.7

        low_freqs = [f for f in freq_vals if f <= 3]
        if low_freqs:
            low_counts = np.array([self.frequencies[f] for f in low_freqs])
            if np.sum(low_counts) > 0:
                mean_low = np.sum(np.array(low_freqs) * low_counts) / np.sum(low_counts)
                p_err = 1.0 / max(mean_low, 0.5)
            else:
                p_err = 0.5
        else:
            p_err = 0.5
        p_err = np.clip(p_err, 0.1, 0.9)
        peak_idx = np.argmax(counts)
        mu = float(freq_vals[peak_idx])
        sigma = max(1.0, mu * 0.2)
        return w_err, p_err, mu, sigma, total_obs

    def _initialize_params_robust(self, n_signal_components=2):
        """
        修正版：保留 SG 找谷底的逻辑，修正初始参数计算，避免 w_err 过高
        """
        from scipy.signal import savgol_filter, find_peaks
        import numpy as np

        freq_vals, counts = self._get_arrays()
        total_obs = np.sum(counts)
        if total_obs == 0:
            raise ValueError("No observed k-mers")

        # 1. 保留 SG 平滑和找谷底的逻辑（用于找到信号区域的大致起点）
        window = min(11, len(counts) // 4)
        if window % 2 == 0:
            window += 1
        if len(counts) >= window and window >= 5:
            smoothed = savgol_filter(counts, window_length=window, polyorder=3)
        else:
            smoothed = counts.copy()

        prominence = 0.05 * (np.max(smoothed) - np.min(smoothed) + 1e-12)
        valleys, _ = find_peaks(-smoothed, prominence=prominence, distance=3)
        cutoff_idx = valleys[0] if len(valleys) > 0 and valleys[0] < len(freq_vals)-2 else 5
        # cutoff_idx 是谷底位置，代表错误和信号的大致分界

        # 2. 修正 w_err 与 p_err 初始化：基于谷底分界，兼顾统计合理性与 EM 稳定性
        valley_freq = freq_vals[cutoff_idx]  # 谷底对应的实际频率值
        
        # 计算错误权重初值：使用谷底左侧所有频率（f ≤ f_valley）
        err_mask = freq_vals <= valley_freq
        err_counts = counts[err_mask]
        w_err = np.sum(err_counts) / total_obs if np.sum(err_counts) > 0 else 0.2
        # 安全约束：限制上限防止 EM 陷入“错误吞噬信号”的局部最优
        w_err = np.clip(w_err, 0.05, 0.3)
        
        # 计算几何分布衰减参数 p_err：仅用极低频（抗稀有物种污染）
        p_err_mask = freq_vals <= min(5, valley_freq)
        p_err_counts = counts[p_err_mask]
        p_err_freqs = freq_vals[p_err_mask]
        mean_err = np.sum(p_err_freqs * p_err_counts) / np.sum(p_err_counts) if np.sum(p_err_counts) > 0 else 1.5
        p_err = 1.0 / max(mean_err, 0.5)
        p_err = np.clip(p_err, 0.1, 0.9)

        # 3. 修正信号区域参数：用谷底后的频率估算 mu 和 sigma（更接近真实信号）
        sig_mask = freq_vals >= cutoff_idx  # 信号区域从谷底开始
        sig_counts = counts[sig_mask]
        sig_freqs = freq_vals[sig_mask]

        if np.sum(sig_counts) > 0:
            if n_signal_components == 2:
                # 双峰模式：用谷底后的频率找两个峰
                sig_peaks, props = find_peaks(sig_counts, prominence=0.1*np.max(sig_counts), distance=5)
                if len(sig_peaks) >= 2:
                    peak_prominences = props['prominences']
                    top_two_indices = np.argsort(peak_prominences)[-2:]
                    p1_idx, p2_idx = sig_peaks[top_two_indices[0]], sig_peaks[top_two_indices[1]]
                    peak1_idx = min(p1_idx, p2_idx)  # 稀有种（低频峰）
                    peak2_idx = max(p1_idx, p2_idx)  # 优势种（高频峰）
                else:
                    # 备选：用最高峰 + 均值分割
                    peak2_idx = np.argmax(sig_counts)
                    lower_half = sig_counts[:len(sig_counts)//2]
                    peak1_idx = np.argmax(lower_half) if np.max(lower_half) > 0.1*np.max(sig_counts) else max(0, peak2_idx // 2)

                # 优势种参数
                mu = sig_freqs[peak2_idx]
                sigma = max(sig_freqs[peak2_idx] * 0.2, 1.0)
                # 稀有种参数
                self.mu_rare = sig_freqs[peak1_idx] if peak1_idx < len(sig_freqs) else mu * 0.3
                self.sigma_rare = max(self.mu_rare * 0.4, 0.5)
                peak1_height = sig_counts[peak1_idx] if peak1_idx < len(sig_counts) else 0
                peak2_height = sig_counts[peak2_idx]
                self.w_rare_frac_init = np.clip(peak1_height / (peak1_height + peak2_height + 1e-10), 0.05, 0.5)
            else:
                # 单峰模式：用最高峰估算
                peak_local_idx = np.argmax(sig_counts)
                mu = sig_freqs[peak_local_idx]
                median_val = np.median(sig_freqs)
                mad = np.median(np.abs(sig_freqs - median_val))
                sigma = mad * 1.4826
                sigma = max(sigma, 0.5)
        else:
            # 极端情况：没有信号区域，给默认值
            mu = 10.0
            sigma = 2.0

        return w_err, p_err, mu, sigma, total_obs
    # ------------------------------------------------------------------
    # EM fitting
    # ------------------------------------------------------------------
    def fit(self, method='robust', n_signal_components=2):
        """
        Run EM algorithm to estimate mixture parameters.

        Parameters
        ----------
        method : str
            'robust' (default) — Savitzky-Golay smoothing + valley detection
            'simple'           — basic heuristic
        n_signal_components : int
            1 = unimodal signal (single species)
            2 = bimodal signal (dominant + rare species) [default for metagenomes]

        Returns
        -------
        dict : fitted parameters
        """
        # 1. Initialize parameters
        
        if method == 'robust':
            try:
                w_err, p_err, mu, sigma, total_obs = self._initialize_params_robust(
                    n_signal_components=n_signal_components)
            except Exception as e:
                sys.stderr.write(f"Robust init failed ({e}), falling back to simple.\n")
                w_err, p_err, mu, sigma, total_obs = self._initialize_params()
        else:
            w_err, p_err, mu, sigma, total_obs = self._initialize_params()

        freq_vals, counts = self._get_arrays()
        
        # Initialize rare species parameters (for bimodal)
        if n_signal_components == 2:
            mu_rare = getattr(self, 'mu_rare', mu * 0.3)
            sigma_rare = getattr(self, 'sigma_rare', sigma * 0.5)
            w_rare_frac =getattr(self, 'w_rare_frac_init', 0.3)
        else:
            mu_rare = sigma_rare = w_rare_frac = 0.0

        # 2. EM iterations
        for iteration in range(self.max_em_iter):
            # ===== E-step =====
            p_err_vals = geom.pmf(freq_vals, p_err) #几何分布的结果
            
            if n_signal_components == 1:
                # Unimodal signal
                p_sig_vals = norm.pdf(freq_vals, mu, sigma)
                joint_err = w_err * p_err_vals
                joint_sig = (1 - w_err) * p_sig_vals
            else:
                # 🔧 Bimodal signal: dominant + rare
                p_sig_dom = norm.pdf(freq_vals, mu, sigma)  # Dominant species
                p_sig_rare = norm.pdf(freq_vals, mu_rare, sigma_rare)  # Rare species
                
                # Mixture within signal component
                joint_err = w_err * p_err_vals
                joint_sig = (1 - w_err) * ((1 - w_rare_frac) * p_sig_dom + 
                                          w_rare_frac * p_sig_rare)
            
            total_density = joint_err + joint_sig + 1e-12
            gamma = joint_err / total_density  # P(error | x)

            # ===== M-step =====
            # Update error weight
            w_err_new = np.sum(gamma * counts) / total_obs

            # Update error distribution parameter
            sum_gamma = np.sum(gamma * counts)
            sum_gamma_x = np.sum(gamma * counts * freq_vals)
            p_err_new = sum_gamma / sum_gamma_x if sum_gamma_x > 0 else p_err
            p_err_new = np.clip(p_err_new, 0.01, 0.99)

            # Update signal distribution(s)
            weights_sig = (1 - gamma) * counts
            sum_weights = np.sum(weights_sig)
            
            if sum_weights > 0:
                if n_signal_components == 1:
                    # Unimodal update
                    mu_new = np.sum(weights_sig * freq_vals) / sum_weights
                    var_new = np.sum(weights_sig * (freq_vals - mu_new) ** 2) / sum_weights
                    sigma_new = np.sqrt(max(var_new, 0.1))
                else:
                    # 🔧 Bimodal update (proper EM for rare species)
                    
                    # Step 1: 计算每个频率点属于"稀有种"的后验概率（在信号成分内部）
                    # γ_rare[f] = P(rare | signal, freq=f)
                    p_sig_dom = norm.pdf(freq_vals, mu, sigma)
                    p_sig_rare = norm.pdf(freq_vals, mu_rare, sigma_rare)
                    
                    # 避免除零
                    denom = (1 - w_rare_frac) * p_sig_dom + w_rare_frac * p_sig_rare + 1e-12
                    gamma_rare = (w_rare_frac * p_sig_rare) / denom  # P(rare | signal, x)
                    
                    # Step 2: 用 (1-γ)×γ_rare 作为稀有种的权重，更新其参数
                    weights_rare = weights_sig * gamma_rare  # 稀有种的加权计数
                    sum_weights_rare = np.sum(weights_rare)
                    
                    if sum_weights_rare > 10:  # 至少10个有效观测才更新
                        mu_rare_new = np.sum(weights_rare * freq_vals) / sum_weights_rare
                        var_rare_new = np.sum(weights_rare * (freq_vals - mu_rare_new) ** 2) / sum_weights_rare
                        sigma_rare_new = np.sqrt(max(var_rare_new, 0.1))
                        
                        # 更新稀有种权重（在信号成分内的占比）
                        w_rare_frac_new = sum_weights_rare / sum_weights
                        w_rare_frac_new = np.clip(w_rare_frac_new, 0.05, 0.5)  # 限制在5%~50%
                    else:
                        # 稀有种观测不足，保持原参数 + 轻微向主峰收缩（防漂移）
                        mu_rare_new = 0.9 * mu_rare + 0.1 * mu
                        sigma_rare_new = sigma_rare
                        w_rare_frac_new = w_rare_frac
                    
                    # Step 3: 更新主峰参数（用总信号权重 - 稀有种权重）
                    weights_dom = weights_sig * (1 - gamma_rare)
                    sum_weights_dom = np.sum(weights_dom)
                    
                    if sum_weights_dom > 10:
                        mu_new = np.sum(weights_dom * freq_vals) / sum_weights_dom
                        var_new = np.sum(weights_dom * (freq_vals - mu_new) ** 2) / sum_weights_dom
                        sigma_new = np.sqrt(max(var_new, 0.1))
                    else:
                        mu_new = mu
                        sigma_new = sigma
                    
                    # 保存更新后的稀有种参数
                    mu_rare = mu_rare_new
                    sigma_rare = sigma_rare_new
                    w_rare_frac = w_rare_frac_new
            else:
                mu_new = mu
                sigma_new = sigma
     
            # Check convergence
            delta = max(abs(w_err_new - w_err), abs(p_err_new - p_err),
                        abs(mu_new - mu), abs(sigma_new - sigma))
            if delta < self.tol:
                self.converged = True
                w_err, p_err, mu, sigma = w_err_new, p_err_new, mu_new, sigma_new
                if n_signal_components == 2:
                    self.mu_rare = mu_rare
                    self.sigma_rare = sigma_rare
              
                break

            w_err, p_err, mu, sigma = w_err_new, p_err_new, mu_new, sigma_new

        # 3. Store results
        self.params = {
            'w_err': w_err,
            'w_sig': 1 - w_err,
            'p_err': p_err,
            'mu': mu,
            'sigma': sigma,
            'total_obs': total_obs,
            'iterations': iteration + 1,
            'converged': self.converged,
            'n_signal_components': n_signal_components
        }
        
        # Store rare species parameters if bimodal
        if n_signal_components == 2:
            self.params['mu_rare'] = mu_rare
            self.params['sigma_rare'] = sigma_rare
            self.params['w_rare_frac'] = w_rare_frac

        return self.params
    # ------------------------------------------------------------------
    # Posterior computation  计算每个频率下，k‑mer 属于“错误成分”的后验概率

    def compute_posterior(self, max_x=None):
        """
        Compute P(error | x) for x = 1..max_x.
        
        For bimodal models, integrates over both signal components.
        
        Returns
        -------
        tuple : (x_vals, posterior_error_probabilities)
        """
        if self.params is None:
            self.fit()
        if max_x is None:
            max_x = max(int(2 * self.params['mu']), self.max_freq)

        x_vals = np.arange(1, max_x + 1)
        p_err_vals = geom.pmf(x_vals, self.params['p_err']) #每个频率的几何分布结果
        
        # Check if bimodal
        if self.params.get('n_signal_components', 1) == 2:
            # Bimodal signal
            mu_rare = self.params.get('mu_rare', self.params['mu'] * 0.3)
            sigma_rare = self.params.get('sigma_rare', self.params['sigma'] * 0.5)
            w_rare_frac = self.params.get('w_rare_frac', 0.3)
            
            p_sig_dom = norm.pdf(x_vals, self.params['mu'], self.params['sigma'])
            p_sig_rare = norm.pdf(x_vals, mu_rare, sigma_rare)
            p_sig_vals = (1 - w_rare_frac) * p_sig_dom + w_rare_frac * p_sig_rare
        else:
            # Unimodal signal
            p_sig_vals = norm.pdf(x_vals, self.params['mu'], self.params['sigma'])
        
        joint_err = self.params['w_err'] * p_err_vals
        joint_sig = self.params['w_sig'] * p_sig_vals
        posterior_err = joint_err / (joint_err + joint_sig + 1e-12) #P(error∣x) 每个频率的后验概率

        return x_vals, posterior_err
    # ------------------------------------------------------------------
    # [FIX] Innovation 1: CRITIC-inspired adaptive cutoff
    
    def find_cutoff_critic(self):
        """
        CRITIC-inspired adaptive cutoff based on posterior gradient contrast.
        
        [FIX for bimodal] In bimodal mode, ensure cutoff separates error from 
        BOTH rare and dominant species.

        Returns
        -------
        int : adaptive cutoff frequency
        """
        if self.params is None:
            self.fit()
        
        x_vals, post = self.compute_posterior() # post每个频率的后验概率
        
        if len(x_vals) < 3:
            self.cutoff = 1
            return self.cutoff

        # 🔧 检测是否为双峰模式
        is_bimodal = self.params.get('n_signal_components', 1) == 2
        
        if is_bimodal:
            # 双峰模式：保护稀有种，cutoff 不应超过稀有种均值
            mu_rare = self.params.get('mu_rare', self.params['mu'] * 0.3)
            # 关键修复：使用更高的保护阈值，确保稀有物种不被误杀
            rare_protect_limit = int(max(2, mu_rare * 0.9))  # 从 0.7 提高到 0.9
            
            # 策略1：进一步放宽阈值到 0.8（而非 0.7），更严格保护稀有种
            threshold_idx = None
            for i, p in enumerate(post):
                if p < 0.8:  # ← 双峰模式下用 0.8，更保守
                    threshold_idx = i
                    break
            
            # 策略2：梯度法限制在低频区搜索（≤ rare_protect_limit）
            search_end = min(len(post) - 1, rare_protect_limit)
            search_end = max(search_end, 1)
            gradient = np.diff(post[:search_end + 1])
            cliff_idx = int(np.argmin(gradient))
            
            # 保守融合：取靠后的点，但不超过保护阈值
            if threshold_idx is not None:
                final_idx = max(cliff_idx, threshold_idx)
                final_idx = min(final_idx, rare_protect_limit)  # ← 关键保护
            else:
                final_idx = cliff_idx
            
            self.cutoff = int(x_vals[final_idx])
            
        else:
            # 单峰模式：原逻辑保持不变
            gradient = np.diff(post)
            cliff_idx = int(np.argmin(gradient))
            
            threshold_idx = None
            for i, p in enumerate(post):
                if p < self.posterior_threshold:
                    threshold_idx = i
                    break
            
            if threshold_idx is not None:
                self.cutoff = int(x_vals[max(cliff_idx, threshold_idx)])
            else:
                self.cutoff = max(1, int(x_vals[cliff_idx]))

        return self.cutoff
    # ------------------------------------------------------------------
    # [KEPT] Original hard cutoff — fallback only
    # ------------------------------------------------------------------
    # def find_cutoff(self):
    #     """Original hard cutoff using fixed posterior_threshold."""
    #     if self.params is None:
    #         self.fit()
    #     x_vals, post = self.compute_posterior()
    #     cutoff = None
    #     for x, p in zip(x_vals, post):
    #         if p < self.posterior_threshold:
    #             cutoff = int(x)
    #             break
    #     if cutoff is None:
    #         cutoff = int(x_vals[-1])
    #     self.cutoff = cutoff
    #     return cutoff

    # ------------------------------------------------------------------
    # [KEPT] Original hard effective_diff — backward compatible
    # ------------------------------------------------------------------
    def compute_effective_diff(self):
        """
        Hard effective k-mers: count of distinct k-mers with freq >= cutoff.
        [FIX] Ensure distinct_kmers and cutoff are properly initialized.
        """
        # 确保 cutoff 已计算
        if self.cutoff is None:
            self.find_cutoff_critic()
        
        # 确保 distinct_kmers 已设置（关键修复！）
        if self.distinct_kmers is None or self.distinct_kmers <= 0:
            if self.frequencies:
                self.distinct_kmers = sum(self.frequencies.values())
            else:
                logger.warning("No frequency data for distinct_kmers calculation")
                return 0, 0
        
        # 计算低频错误 k-mer 数量 (freq < cutoff)
        error_sum = sum(
            count for freq, count in self.frequencies.items()
            if 1 <= freq < self.cutoff
        )
        
        # 有效 k-mer = 总数 - 错误数
        eff = self.distinct_kmers - error_sum
        
        # 防御性编程：防止负数（可能由于数值误差）
        self.effective_diff = max(0, int(eff))
        self.error_sum = int(error_sum)
        
        return self.effective_diff, self.error_sum

    # ------------------------------------------------------------------
    # [FIX] Innovation 2: Posterior-Expected Soft-Counting (PESC)
    # ------------------------------------------------------------------
    def compute_effective_diff_soft(self):
        """
        [ABLATION: w/o PESC] 
        Fallback to hard counting to demonstrate the necessity of soft boundary.
        Instead of continuous soft boundary, uses binary hard cutoff.
        """
        # 直接复用已经计算好的硬阈值结果 (effective_diff)
        # 因为在 get_metrics() 中，compute_effective_diff() 会先于本函数被调用
        if self.effective_diff is None:
            self.compute_effective_diff()
        
        # 将软计数变量强制赋值为硬计数结果，从而在下游彻底关闭 PESC
        self.soft_effective_diff = float(self.effective_diff)
        return self.soft_effective_diff
    # ------------------------------------------------------------------
    # [FIX] Innovation 3: Boundary uncertainty (Shannon entropy)
    # ------------------------------------------------------------------
    def compute_boundary_uncertainty(self, window=2):
        """
        Shannon entropy of k-mer membership in the boundary region.

        Low  entropy = sharp boundary  = high confidence
        High entropy = ambiguous       = noisy / overlapping components

        This serves as a natural confidence weight (4th CRITIC objective).

        Parameters
        ----------
        window : int
            Half-width of boundary region around cutoff (default 2).

        Returns
        -------
        float : normalized boundary entropy in [0, ln2]
        """
        if self.params is None:
            self.fit()
        if self.cutoff is None:
            self.find_cutoff_critic()

        x_vals, post_err = self.compute_posterior(
            max_x=self.cutoff + window + 1)

        window_freqs = range(
            max(1, self.cutoff - window), self.cutoff + window + 1)
        boundary_entropy = 0.0
        boundary_count = 0

        for f in window_freqs:
            count = self.frequencies.get(f, 0)
            if count > 0:
                idx = f - 1
                p_err = post_err[idx] if 0 <= idx < len(post_err) else 0.0
                p_sig = 1.0 - p_err

                if 0 < p_err < 1:
                    entropy = -(p_err * np.log(p_err) +
                                p_sig * np.log(p_sig))
                else:
                    entropy = 0.0

                boundary_entropy += entropy * count
                boundary_count += count

        self.boundary_uncertainty = (
            boundary_entropy / boundary_count if boundary_count > 0 else 0.0)
        return self.boundary_uncertainty

    # ------------------------------------------------------------------
    # [FIX] Unified metrics interface
    # ------------------------------------------------------------------
    def get_metrics(self, method='robust', n_signal_components=2):
        """
        Run full analysis and return all metrics.

        Always computes:
          - CRITIC adaptive cutoff
          - Hard effective_diff (backward compatible)
          - PESC soft_effective_diff (primary metric for partitioning)
          - Boundary uncertainty (4th CRITIC objective)

        Parameters
        ----------
        method : str
            'robust' or 'simple'
        n_signal_components : int
            1 = unimodal, 2 = bimodal (recommended for metagenomes)

        Returns
        -------
        dict : all analysis results
        """
        if self.params is None:
            # 🔧 关键修复：传递双峰参数
            self.fit(method=method, n_signal_components=n_signal_components)

        # CRITIC adaptive cutoff
        if self.cutoff is None:
            self.find_cutoff_critic()

        # Hard effective_diff (backward compatible)
        if self.effective_diff is None:
            self.compute_effective_diff()

        # PESC soft count
        self.compute_effective_diff_soft()

        # Boundary uncertainty
        self.compute_boundary_uncertainty()

        # 🔍 计算污染代理指标 (freq ∈ [3, μ/3])
        mu_val = self.params.get('mu', 10.0)
        upper_bound = max(3, int(mu_val / 3.0))
        foreign_count = sum(
            count for freq, count in self.frequencies.items()
            if 3 <= freq <= upper_bound
        )
        contamination_proxy = foreign_count / self.distinct_kmers if self.distinct_kmers > 0 else 0.0

        # 📊 覆盖度异质性 (仅双峰模式有效)
        mu_rare_val = self.params.get('mu_rare', 0.0)
        if n_signal_components == 2 and mu_rare_val > 0:
            coverage_heterogeneity = min(mu_val / max(mu_rare_val, 1.0), 5.0)  # 上限5防极端值
        else:
            coverage_heterogeneity = 1.0  # 单峰默认无异质性惩罚

        # 🔧 关键修复：构建包含双峰参数的字典
        metrics = {
            'k': getattr(self, 'k_value', None),
            'distinct_kmers': self.distinct_kmers,
            'w_err': self.params['w_err'],
            'w_sig': self.params['w_sig'],
            'p_err': self.params['p_err'],
            'mu': self.params['mu'],
            'sigma': self.params['sigma'],
            'cutoff': self.cutoff,
            'error_sum': self.error_sum,
            # Traditional (hard)
            'effective_diff': self.effective_diff,
            # Innovation metrics
            'soft_effective_diff': self.soft_effective_diff,
            'boundary_uncertainty': self.boundary_uncertainty,
            # Diagnostics
            'em_converged': self.params['converged'],
            'em_iterations': self.params['iterations'],
            'max_freq_observed': self.max_freq,
            'num_freq_bins': len(self.frequencies),
            # 🔧 新增：双峰模型参数（如果是双峰）
            'n_signal_components': self.params.get('n_signal_components', 1),
            'mu_rare': self.params.get('mu_rare', None),
            'sigma_rare': self.params.get('sigma_rare', None),
            'w_rare_frac': self.params.get('w_rare_frac', None),
            # 🆕 新增：纯频谱统计指标
            'contamination_proxy': contamination_proxy,
            'coverage_heterogeneity': coverage_heterogeneity,
        }
        return metrics

    def get_spectrum_stats(self):
        """Basic statistics about the loaded spectrum."""
        if not self.frequencies:
            return {}
        freqs = list(self.frequencies.keys())
        counts = list(self.frequencies.values())
        return {
            'distinct_kmers': self.distinct_kmers,
            'total_kmer_occurrences': sum(counts),
            'min_frequency': min(freqs),
            'max_frequency': max(freqs),
            'peak_frequency': freqs[np.argmax(counts)],
            'peak_count': max(counts),
            'num_frequency_bins': len(freqs)
        }
