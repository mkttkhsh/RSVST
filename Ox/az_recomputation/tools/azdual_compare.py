#!/usr/bin/env python3
"""AZ 系 20 対象の修正再計算（dual 版、§85-86）の比較。§83 の azfix_collect_check.py の比較節を 20 対象に一般化。
2 つの差を分けて報告する（WORKFLOW・§85）:
  (A) dual 内 fixed − orig: 共通ドロー上の予測式修正効果（§83 で posstat 等が bitwise 一致することを確認済みの構造）
  (B) orig − 論文掲載 run: 旧式をそろえた run・計算設定間の差（seed・Ox 版・分割・実行機・burn-in の位置が違う）
  (C) fixed − 論文掲載 run: 掲載値の更新全体（= A + B の窓ごとの和だが、順位・MCS・GW は差分の足し算では説明しない）
  尺度: mean(exp h) は対数比、VaR/ES は相対差、損失は平均の差。QLIKE の尺度から FZ0 の影響の大小は判断しない
窓レベル（aggregated の windows_<model>.csv から）: mean(exp h) の相対差・対数比、VaR/ES 1%・5% の相対差、QLIKE（HL 修正後 proxy 4 種）、FZ0 1%・5%、VaR 違反日数と判定が変わった日。
表レベル（評価の再集計 *_r3corrected_<tag>/ があれば）: 各標本 × proxy の QLIKE・順位・MCS p、FZ0 の順位・MCS p、matched pair（RSV−SV）の向きを、
  *_r3corrected（論文掲載 run）/ _dualorig / _dualfixed で並べる。ψ 感度表（COVID）も同様。
出力: r3_tables/azdual_comparison.md、r3_tables/azdual_comparison_windows.csv
使い方: python3 tools/azdual_compare.py [--tag dualfixed] [--res DIR]   部分集約でも窓レベルは動く（存在する対象だけ）。
"""
import argparse, csv, math, os
ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="dualfixed"); ap.add_argument("--orig-tag", default="dualorig"); ap.add_argument("--res", default=None)
ap.add_argument("--tables-only", action="store_true", help="2026-10-11: 窓レベル（集約の windows_*.csv が必要）を飛ばし表レベルだけ（例: --tag dualfixed_excl --orig-tag dualfixed の除外感度）")
ap.add_argument("--out", default="azdual_comparison", help="出力名（r3_tables/<out>.md と <out>_windows.csv）")
A = ap.parse_args()
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); IJF = os.path.join(ROOT, "20260806IJF_revision")
RES = os.path.abspath(A.res) if A.res else os.path.join(IJF, "azdual_results"); AGG = os.path.join(RES, "aggregated"); OUT = os.path.join(IJF, "r3_tables")
SAMPLES = {"DJIA_main": ("DJIA_prediction", os.path.join(ROOT, "Ox2022", "DJIA_prediction"), os.path.join(ROOT, "R", "prediction_results_2026"), "DJIA", 603),
           "N225_main": ("N225_prediction", os.path.join(ROOT, "Ox2022", "N225_prediction"), os.path.join(ROOT, "R", "prediction_results_2026"), "N225", 590),
           "DJIA_2020": ("DJIA_prediction_2020", os.path.join(IJF, "reproducibility", "Ox", "DJIA_prediction_2020"), os.path.join(IJF, "reproducibility", "R", "prediction_results_djia_2020"), "DJIA", 253),
           "N225_2020": ("N225_prediction_2020", os.path.join(IJF, "reproducibility", "Ox", "N225_prediction_2020"), os.path.join(IJF, "reproducibility", "R", "prediction_results_n225_2020"), "N225", 242)}
MODELS = {"DJIA_main": ["sv_az_sn", "sv_az_st", "rsv_az_sn", "rsv_az_st"], "N225_main": ["sv_az_sn", "sv_az_st", "rsv_az_sn", "rsv_az_st"],
          "DJIA_2020": ["sv_az_sn", "sv_az_st", "rsv_az_sn", "rsv_az_st", "rsv_az_sn_psi", "rsv_az_st_psi"], "N225_2020": ["sv_az_sn", "sv_az_st", "rsv_az_sn", "rsv_az_st", "rsv_az_sn_psi", "rsv_az_st_psi"]}
NAME = {"sv_az_sn": "SV-AZ-SN", "sv_az_st": "SV-AZ-ST", "rsv_az_sn": "RSV-AZ-SN", "rsv_az_st": "RSV-AZ-ST", "rsv_az_sn_psi": "RSV-AZ-SN (psi-free)", "rsv_az_st_psi": "RSV-AZ-ST (psi-free)"}
PROX = ("RV5", "RKth2", "BV", "MedRV")
def rd(p): return list(csv.DictReader(open(p)))
def q(v, p):
    v = sorted(v); k = (len(v) - 1) * p; f = math.floor(k); c = min(f + 1, len(v) - 1); return v[f] + (v[c] - v[f]) * (k - f)
def sd(d):
    m = sum(d) / len(d); return math.sqrt(sum((x - m) ** 2 for x in d) / (len(d) - 1)) if len(d) > 1 else float("nan")
def summ(d, idx, dates):
    a = [abs(x) for x in d]; k = max(range(len(a)), key=lambda i: a[i])
    return (f"mean {sum(d)/len(d):+.5f}, med {q(d,.5):+.5f}, SD {sd(d):.5f}, q01 {q(d,.01):+.5f}, q99 {q(d,.99):+.5f}; abs: med {q(a,.5):.5f}, q95 {q(a,.95):.5f}, max {max(a):.5f} at t={idx[k]} ({dates[idx[k]]})")
def nw_tstat(d, lag):
    n = len(d); m = sum(d) / n; e = [x - m for x in d]; s = sum(x * x for x in e) / n
    for l in range(1, lag + 1):
        w = 1 - l / (lag + 1); s += 2 * w * sum(e[i] * e[i - l] for i in range(l, n)) / n
    return m / math.sqrt(s / n) if s > 0 else float("nan")
def qlike(y, s): return y / s - math.log(y / s) - 1
def fz0(y, v, e, a): return -(1 / (a * e)) * (1 if y <= v else 0) * (v - y) + v / e + math.log(-e) - 1
def _runs(ts):
    out = []; a = ts[0]; p = ts[0]
    for t in ts[1:]:
        if t != p + 1: out.append((a, p)); a = t
        p = t
    out.append((a, p)); return out
def paper_series(pdir, m):  # 論文掲載 run の集約ファイル → 窓別 mean(exp h), var1, var5, es1, es5
    def cols(p):
        with open(p) as f:
            r = csv.reader(f); hdr = next(r); rows = [row for row in r if row]
        off = 1 if hdr[0] in ("", '""') or hdr[0].startswith(",") else 0
        if hdr[0] == "" or off: return [[float(x) for x in row[1:]] for row in rows]
        return [[float(x) for x in row] for row in rows]
    vf, va, es = cols(os.path.join(pdir, f"volf_{m}.csv")), cols(os.path.join(pdir, f"var_{m}.csv")), cols(os.path.join(pdir, f"es_{m}.csv"))
    return [{"mexp": vf[t][2], "var1": va[t][1], "var5": va[t][2], "es1": es[t][1], "es5": es[t][2]} for t in range(len(vf))]
doc = ["# AZ 再計算（dual）の比較: 予測式の効果（fixed − orig、同じドロー）と掲載 run との差（fixed − 論文）", "",
       f"tag = {A.tag} / {A.orig_tag}。QLIKE の proxy は HL 修正後（*_r3corrected の評価 CSV）、FZ0 のリターンは評価 CSV の Return 列。NW のラグ = floor(0.75 n^(1/3)) − 1。", ""]
wrows = []
# ---------- 窓レベル ----------
if not A.tables_only:
    doc += ["## 窓レベル", ""]
    for sample, (sdir, pdir, evdir, idx, n) in SAMPLES.items():
        ev = rd(os.path.join(evdir + "_r3corrected", f"{idx}_forecast_vol.csv")); v1 = rd(os.path.join(evdir + "_r3corrected", f"{idx}_forecast_var1p.csv"))
        dates = [r["date"] for r in ev]; y = [float(r["Return"]) for r in v1]
        lag = max(int(math.floor(0.75 * n ** (1 / 3))) - 1, 0)
        for m in MODELS[sample]:
            pf, po = os.path.join(AGG, A.tag, sdir, f"windows_{m}.csv"), os.path.join(AGG, A.orig_tag, sdir, f"windows_{m}.csv")
            if not os.path.exists(pf): continue
            F = {int(r["t"]): r for r in rd(pf) if r["mean_exp_h"] not in ("", "nan", "NA")}; O = {int(r["t"]): r for r in rd(po)} if os.path.exists(po) else {}
            P = paper_series(pdir, m)
            T = sorted(t for t in F if t in O) if O else []; TP = sorted(F)
            complete = len(F) == n
            cov = "完全（全窓）" if complete else f"部分: {len(F)}/{n} 窓、連続区間 " + ", ".join(f"{a}-{b}" for a, b in _runs(sorted(F)))
            doc += [f"### {sample} {NAME[m]}（{cov}、fixed − orig 比較可能 {len(T)}）", ""]
            if not complete: doc += ["部分データ: 記述統計のみ。不連続な窓を詰めた NW t は表示しない（正式な検定は完全な共通日付系列で行う）。", ""]
            nwt = (lambda d: f"{nw_tstat(d, lag):+.2f}") if complete else (lambda d: "n/a (partial)")
            g = lambda D, t, k: float(D[t][k if k != "mexp" else "mean_exp_h"])
            if T:
                rel = [(g(F, t, "mexp") - g(O, t, "mexp")) / g(O, t, "mexp") for t in T]; lr = [math.log(g(F, t, "mexp") / g(O, t, "mexp")) for t in T]
                doc += ["(A) fixed − orig（同じドロー、予測式の効果）", "", f"- mean(exp h) 相対差: {summ(rel, T, dates)}", f"- log(fixed/orig): {summ(lr, T, dates)}（符号 負 {sum(x < 0 for x in lr)}・正 {sum(x > 0 for x in lr)}）"]
                for k, nm in (("var1", "VaR 1%"), ("var5", "VaR 5%"), ("es1", "ES 1%"), ("es5", "ES 5%")):
                    r = [(g(F, t, k) - g(O, t, k)) / abs(g(O, t, k)) for t in T]; doc.append(f"- {nm} 相対差: {summ(r, T, dates)}")
                doc += ["", "| 損失 | orig | fixed | fixed − orig | NW t |", "|---|---|---|---|---|"]
                for px in PROX:
                    lo = [qlike(float(ev[t][px]), g(O, t, "mexp")) for t in T]; lf = [qlike(float(ev[t][px]), g(F, t, "mexp")) for t in T]; d = [lf[i] - lo[i] for i in range(len(T))]
                    doc.append(f"| QLIKE {px} | {sum(lo)/len(T):.5f} | {sum(lf)/len(T):.5f} | {sum(d)/len(T):+.6f} | {nwt(d)} |")
                for a, kv, ke in ((0.01, "var1", "es1"), (0.05, "var5", "es5")):
                    lo = [fz0(y[t], g(O, t, kv), g(O, t, ke), a) for t in T]; lf = [fz0(y[t], g(F, t, kv), g(F, t, ke), a) for t in T]; d = [lf[i] - lo[i] for i in range(len(T))]
                    vo = sum(y[t] <= g(O, t, kv) for t in T); vf_ = sum(y[t] <= g(F, t, kv) for t in T); chg = sum((y[t] <= g(O, t, kv)) != (y[t] <= g(F, t, kv)) for t in T)
                    doc.append(f"| FZ0 {int(a*100)}% | {sum(lo)/len(T):.5f} | {sum(lf)/len(T):.5f} | {sum(d)/len(T):+.6f} | {nwt(d)} |")
                    doc.append(f"| VaR {int(a*100)}% 違反日数 | {vo} | {vf_} | 判定が変わった日 {chg} | |")
                doc.append("")
            if TP and len(P) == n:
                lp = [math.log(g(F, t, "mexp") / P[t]["mexp"]) for t in TP]
                doc += ["(B) fixed − 論文掲載 run（修正効果 ＋ run 間変動）", "", f"- log(fixed/論文): {summ(lp, TP, dates)}; |log| > 0.5: {sum(abs(x) > 0.5 for x in lp)} 窓, > 1.0: {sum(abs(x) > 1.0 for x in lp)} 窓"]
                if O and T: lo_ = [math.log(g(O, t, "mexp") / P[t]["mexp"]) for t in T]; doc.append(f"- log(orig/論文)（同じ旧式の別 run = run 間変動のみ）: {summ(lo_, T, dates)}; |log| > 0.5: {sum(abs(x) > 0.5 for x in lo_)} 窓, > 1.0: {sum(abs(x) > 1.0 for x in lo_)} 窓")
                doc += ["", "| 損失 | 論文 | fixed | fixed − 論文 | NW t |", "|---|---|---|---|---|"]
                for px in PROX:
                    lpp = [qlike(float(ev[t][px]), P[t]["mexp"]) for t in TP]; lf = [qlike(float(ev[t][px]), g(F, t, "mexp")) for t in TP]; d = [lf[i] - lpp[i] for i in range(len(TP))]
                    doc.append(f"| QLIKE {px} | {sum(lpp)/len(TP):.5f} | {sum(lf)/len(TP):.5f} | {sum(d)/len(TP):+.6f} | {nwt(d)} |")
                for a, kv, ke in ((0.01, "var1", "es1"), (0.05, "var5", "es5")):
                    lpp = [fz0(y[t], P[t][kv], P[t][ke], a) for t in TP]; lf = [fz0(y[t], g(F, t, kv), g(F, t, ke), a) for t in TP]; d = [lf[i] - lpp[i] for i in range(len(TP))]
                    doc.append(f"| FZ0 {int(a*100)}% | {sum(lpp)/len(TP):.5f} | {sum(lf)/len(TP):.5f} | {sum(d)/len(TP):+.6f} | {nwt(d)} |")
                doc.append("")
            for t in TP: wrows.append([sample, NAME[m], t, dates[t], g(F, t, "mexp"), g(O, t, "mexp") if t in O else "", P[t]["mexp"] if len(P) == n else ""] + [g(F, t, k) for k in ("var1", "var5", "es1", "es5")] + [g(O, t, k) if t in O else "" for k in ("var1", "var5", "es1", "es5")])
# ---------- 表レベル ----------
doc += ["## 表レベル（評価の再集計）", "", "各列は *_r3corrected（論文掲載 run）/ _" + A.orig_tag + " / _" + A.tag + "。順位はその表のモデル数の中（小さいほど良い）。MCS p > 0.25 が 75% MCS。", ""]
def read_vol(d, idx, px, sfx):
    p = os.path.join(d, f"{idx}_{px}_forecast_vol_mse_qlike_mcs_gw06test{sfx}.csv")
    if not os.path.exists(p): return None
    rr = list(csv.reader(open(p))); return {r[0]: (float(r[3]), float(r[4])) for r in rr[1:]}
def read_fz(d, idx, lv, sfx, names):
    p = os.path.join(d, f"{idx}_forecast_var_es_{lv}_rate_fz0_mcs_dq_gw06test{sfx}.csv")
    if not os.path.exists(p): return None
    rr = list(csv.reader(open(p))); rows = rr[1:]
    if all(r[0].strip('"').isdigit() for r in rows): assert len(rows) == len(names), (p, len(rows), len(names)); rows = [[names[i]] + r[1:] for i, r in enumerate(rows)]
    return {r[0]: (float(r[2]), float(r[3]), float(r[1])) for r in rows}
for sample, (sdir, pdir, evdir, idx, n) in SAMPLES.items():
    for sfx, label in (("", "主表"), ("_psi_sensitivity", "ψ 感度表")):
        if sfx and "2020" not in sample: continue
        variants = [("論文", evdir + "_r3corrected"), (A.orig_tag, evdir + "_r3corrected_" + A.orig_tag), (A.tag, evdir + "_r3corrected_" + A.tag)]
        tabs = {v: read_vol(d, idx, "RV5", sfx) for v, d in variants}
        if tabs[A.tag] is None: doc.append(f"- {sample} {label}: `{evdir}_r3corrected_{A.tag}` に評価出力なし（未実行）"); continue
        names = list(tabs["論文"].keys()); az = [nm for nm in names if "AZ-S" in nm]
        doc += [f"### {sample} {label}", "", "| 損失 | モデル | QLIKE/FZ0 論文 / " + A.orig_tag + " / " + A.tag + " | 順位 | MCS p | matched pair RSV−SV（論文 / " + A.orig_tag + " / " + A.tag + "）|", "|---|---|---|---|---|---|"]
        for px in PROX:
            T3 = {v: read_vol(d, idx, px, sfx) for v, d in variants}
            for nm in az:
                vals, ranks, mcs, pair = [], [], [], []
                for v, _ in variants:
                    tb = T3[v]
                    if tb is None or nm not in tb: vals.append("-"); ranks.append("-"); mcs.append("-"); pair.append("-"); continue
                    order = sorted(tb, key=lambda k: tb[k][0]); vals.append(f"{tb[nm][0]:.4f}"); ranks.append(str(order.index(nm) + 1)); mcs.append(f"{tb[nm][1]:.2f}")
                    mate = ("RSV-" + nm[3:]) if nm.startswith("SV-") else ("SV-" + nm[4:]) if nm.startswith("RSV-") and "psi" not in nm else None
                    if mate and mate in tb: pair.append(f"{(tb['RSV-' + nm[3:]][0] - tb[nm][0]) if nm.startswith('SV-') else (tb[nm][0] - tb[mate][0]):+.4f}")
                    else: pair.append("-")
                doc.append(f"| QLIKE {px} | {nm} | {' / '.join(vals)} | {' / '.join(ranks)} | {' / '.join(mcs)} | {' / '.join(pair)} |")
        for lv in ("1p", "5p"):
            T3 = {v: read_fz(d, idx, lv, sfx, names) for v, d in variants}
            for nm in az:
                vals, ranks, mcs, pair = [], [], [], []
                for v, _ in variants:
                    tb = T3[v]
                    if tb is None or nm not in tb: vals.append("-"); ranks.append("-"); mcs.append("-"); pair.append("-"); continue
                    order = sorted(tb, key=lambda k: tb[k][0]); vals.append(f"{tb[nm][0]:.4f} (rate {tb[nm][2]*100:.2f}%)"); ranks.append(str(order.index(nm) + 1)); mcs.append(f"{tb[nm][1]:.2f}")
                    mate = ("RSV-" + nm[3:]) if nm.startswith("SV-") else ("SV-" + nm[4:]) if nm.startswith("RSV-") and "psi" not in nm else None
                    if mate and mate in tb: pair.append(f"{(tb['RSV-' + nm[3:]][0] - tb[nm][0]) if nm.startswith('SV-') else (tb[nm][0] - tb[mate][0]):+.4f}")
                    else: pair.append("-")
                doc.append(f"| FZ0 {lv[0]}% | {nm} | {' / '.join(vals)} | {' / '.join(ranks)} | {' / '.join(mcs)} | {' / '.join(pair)} |")
        # 全モデルの順位の変化（AZ 以外も MCS 集合が変わるので）
        for px in ("RV5",):
            tp, tf = read_vol(evdir + "_r3corrected", idx, px, sfx), read_vol(evdir + "_r3corrected_" + A.tag, idx, px, sfx)
            if tp and tf:
                op, of = sorted(tp, key=lambda k: tp[k][0]), sorted(tf, key=lambda k: tf[k][0])
                doc += ["", f"QLIKE {px} の順位（論文 → {A.tag}）: " + ", ".join(f"{nm} {op.index(nm)+1}→{of.index(nm)+1}" for nm in op if op.index(nm) != of.index(nm)) + f"; 75% MCS 論文 {sum(v[1] > 0.25 for v in tp.values())} → {A.tag} {sum(v[1] > 0.25 for v in tf.values())} モデル; 最良 {op[0]} → {of[0]}"]
        doc.append("")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, A.out + "_windows.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "model", "t", "date", "mexp_fixed", "mexp_orig", "mexp_paper", "var1_fixed", "var5_fixed", "es1_fixed", "es5_fixed", "var1_orig", "var5_orig", "es1_orig", "es5_orig"]); w.writerows(wrows)
open(os.path.join(OUT, A.out + ".md"), "w").write("\n".join(doc) + "\n"); print("\n".join(doc))
