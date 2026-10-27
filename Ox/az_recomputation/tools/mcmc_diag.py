#!/usr/bin/env python3
"""複数チェーン MCMC 診断（Vehtari, Gelman, Simpson, Carpenter & Bürkner 2021, Bayesian Analysis）の純 Python 実装。
R パッケージ posterior（rhat, ess_bulk, ess_tail, ess_basic, mcse_mean）の定義に合わせる。照合は tools/test_mcmc_diag.py（合成データ、R posterior との突き合わせ）。
入力はいずれも「チェーンのリスト（各チェーンは float のリスト）」。
- split_chains: 各チェーンを前半・後半に分ける（奇数長なら中央の 1 ドローを落とす。posterior と同じ）
- z_scale: 全ドロー込みの順位（同順位は平均順位）→ 正規スコア Φ^{-1}((r − 3/8)/(S + 1/4))
- rhat_basic: 古典的 split-R-hat（分割済みチェーンに適用）
- rhat: rank-normalized split R-hat と folded（|x − median|）版の最大値（posterior::rhat）
- ess_basic: Geyer の初期正列・初期単調列による ESS（分割済みチェーンに適用。posterior の .ess と同じ自己共分散の正規化 1/N、改良版 tau）
- ess_bulk: ess_basic(z_scale(split_chains(x)))、ess_tail: 5%・95% 分位の指示関数の ess_basic（分割）の最小値
- mcse_mean: sd / sqrt(ess_basic(split_chains(x)))
- quantile7: R type 7（Ox の quantilec と同じ線形補間、p·(n−1) 位置）
- batch_mcse_within: チェーン内バッチ平均による MCSE（バッチがチェーン境界をまたがない）
"""
import math

def quantile7(v, p, sorted_=False):
    s = v if sorted_ else sorted(v); n = len(s)
    if n == 0: return float("nan")
    i = p * (n - 1); lo = int(math.floor(i)); hi = min(lo + 1, n - 1)
    return s[lo] + (i - lo) * (s[hi] - s[lo])

def ndtri(p):  # Acklam の有理近似（相対誤差 1.15e-9）
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    if p <= 0: return float("-inf")
    if p >= 1: return float("inf")
    if p < 0.02425:
        q = math.sqrt(-2 * math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > 1 - 0.02425:
        q = math.sqrt(-2 * math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p - 0.5; r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)

def split_chains(chains):
    out = []
    for c in chains:
        n = len(c); h = n // 2
        out.append(list(c[:h])); out.append(list(c[n - h:]))
    return out

def z_scale(chains):
    allv = [(v, i, j) for i, ch in enumerate(chains) for j, v in enumerate(ch)]; allv.sort(key=lambda x: x[0]); S = len(allv)
    out = [[0.0] * len(ch) for ch in chains]; k = 0
    while k < S:
        m = k
        while m + 1 < S and allv[m + 1][0] == allv[k][0]: m += 1
        z = ndtri(((k + m) / 2 + 1 - 3 / 8) / (S + 1 / 4))
        for q in range(k, m + 1): out[allv[q][1]][allv[q][2]] = z
        k = m + 1
    return out

def fold(chains):
    med = quantile7([v for c in chains for v in c], 0.5)
    return [[abs(v - med) for v in c] for c in chains]

def _mean(c): return sum(c) / len(c)
def _var(c, m=None):
    m = _mean(c) if m is None else m; n = len(c); return sum((x - m) ** 2 for x in c) / (n - 1) if n > 1 else float("nan")

def rhat_basic(split):  # 分割済みチェーン
    M = len(split); N = min(len(c) for c in split)
    if M < 2 or N < 2: return float("nan")
    split = [c[:N] for c in split]; means = [_mean(c) for c in split]; gm = _mean(means)
    B = N / (M - 1) * sum((m - gm) ** 2 for m in means); W = _mean([_var(c, m) for c, m in zip(split, means)])
    if not (W > 0): return float("nan")
    return math.sqrt(((N - 1) / N * W + B / N) / W)

def rhat(chains):
    sp = split_chains(chains)
    return max(rhat_basic(z_scale(sp)), rhat_basic(z_scale(split_chains(fold(chains)))))

def _acov(c, m, t):  # 自己共分散（1/N 正規化、posterior::autocovariance と同じ）
    N = len(c); return sum((c[i] - m) * (c[i + t] - m) for i in range(N - t)) / N

def ess_basic(split):  # 分割済みチェーン（posterior の .ess）
    M = len(split); N = min(len(c) for c in split)
    if N < 4: return float("nan")
    split = [c[:N] for c in split]; means = [_mean(c) for c in split]
    acov0 = [_acov(c, m, 0) for c, m in zip(split, means)]
    mean_var = _mean(acov0) * N / (N - 1); var_plus = mean_var * (N - 1) / N + (_var(means) if M > 1 else 0.0)
    if not (var_plus > 0): return float("nan")  # 定数列
    def rho(t): return 1 - (mean_var - _mean([_acov(c, m, t) for c, m in zip(split, means)])) / var_plus
    rho_t = {}; t = 0; rho_even = 1.0; rho_t[0] = rho_even; rho_odd = rho(1); rho_t[1] = rho_odd
    while t < N - 5 and not math.isnan(rho_even + rho_odd) and rho_even + rho_odd > 0:
        t += 2; rho_even = rho(t); rho_odd = rho(t + 1)
        if rho_even + rho_odd >= 0: rho_t[t] = rho_even; rho_t[t + 1] = rho_odd
    max_t = t
    if rho_even > 0: rho_t[max_t] = rho_even
    r = [rho_t.get(i, 0.0) for i in range(max_t + 1)]
    t = 0
    while t <= max_t - 4:
        t += 2
        if r[t] + r[t + 1] > r[t - 2] + r[t - 1]: r[t] = (r[t - 2] + r[t - 1]) / 2; r[t + 1] = r[t]
    ess = M * N; tau = -1 + 2 * sum(r[:max_t]) + r[max_t]
    tau = max(tau, 1 / math.log10(ess)); return ess / tau

def ess_bulk(chains): return ess_basic(z_scale(split_chains(chains)))
def ess_tail(chains):
    allv = sorted(v for c in chains for v in c); q05 = quantile7(allv, 0.05, True); q95 = quantile7(allv, 0.95, True)
    e05 = ess_basic(split_chains([[1.0 if v <= q05 else 0.0 for v in c] for c in chains])); e95 = ess_basic(split_chains([[1.0 if v <= q95 else 0.0 for v in c] for c in chains]))
    return float("nan") if (math.isnan(e05) or math.isnan(e95)) else min(e05, e95)  # posterior と同じく片方が定数なら NA
def ess_mean(chains): return ess_basic(split_chains(chains))
def mcse_mean(chains):
    allv = [v for c in chains for v in c]; return math.sqrt(_var(allv) / ess_mean(chains))

def batch_mcse_within(chains, stat, nb_per_chain=10):
    """チェーンごとに nb_per_chain 個の連続バッチを切り、各バッチの stat の標本 SD / sqrt(総バッチ数)。stat は float のリスト → float。"""
    vals = []
    for c in chains:
        b = len(c) // nb_per_chain
        for i in range(nb_per_chain): vals.append(stat(c[i * b:(i + 1) * b]))
    nb = len(vals); m = _mean(vals); return math.sqrt(sum((v - m) ** 2 for v in vals) / (nb - 1) / nb)

def summary(chains):
    """mean, sd, 2.5%, 97.5%, rhat, ess_bulk, ess_tail, mcse_mean（併合ドローの事後要約）"""
    allv = sorted(v for c in chains for v in c); m = _mean(allv)
    return {"mean": m, "sd": math.sqrt(_var(allv, m)), "q025": quantile7(allv, 0.025, True), "q975": quantile7(allv, 0.975, True),
            "rhat": rhat(chains), "ess_bulk": ess_bulk(chains), "ess_tail": ess_tail(chains), "mcse_mean": mcse_mean(chains)}
