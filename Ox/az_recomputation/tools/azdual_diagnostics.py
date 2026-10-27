#!/usr/bin/env python3
"""AZ 再計算（dual 版、§85-86）の新 run に、主標本 MCMC 診断（tools/mcmc_main_diagnostics.py、§84）と同じ規則を適用する。
入力: azdual_results/<target>/（統合済み）: posstath_<m>_<t>.csv（最終行 = h_T: mean, sd, 2.5%, 97.5%, Geweke p, Inef）、h1_y_<m>_<t>.csv（列 1 = h_{T+1} fixed）
      再利用分 DJIA 主標本 SV-AZ-ST は azfix_results/fixed（要確認 10 窓 t = 30, 33, 35, 73, 108, 116, 117, 118, 134, 446 は §85 で確定済み。ここでも同じ規則で再計算して照合）
規則（全モデル同一、事前固定）: F1 h_T Inef > 300、F2 h_T Geweke p < 0.01（記録のみ）、F3 exp(h1) の先頭 3000 ドロー平均 / 残り 12000 平均 が > 2 または < 0.5。
  確認パスの対象 = F1 ∪ F3。F4（他モデル中央値との乖離 > 1.0、入口のみ）は評価 CSV（*_r3corrected_dualfixed、無ければ *_r3corrected）があれば併記。
出力: r3_tables/azdual_diagnostics.md、r3_tables/azdual_diagnostics_windows.csv（全窓）、r3_tables/azdual_confirm_windows.csv（F1∪F3 ＋ 再利用分 10 窓）
使い方: python3 tools/azdual_diagnostics.py [--res DIR] [--ndraw N]（部分統合でも動く。未統合の窓は未評価）
損失・順位・他モデルへの近さで run や seed を選ばない（§85）。
"""
import argparse, csv, math, os, statistics as st
ap = argparse.ArgumentParser(); ap.add_argument("--res", default=None); ap.add_argument("--ndraw", type=int, default=15000); ap.add_argument("--pkg", default=None)
A = ap.parse_args()
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); IJF = os.path.join(ROOT, "20260806IJF_revision")
RES = os.path.abspath(A.res) if A.res else os.path.join(IJF, "azdual_results"); PKG = os.path.abspath(A.pkg) if A.pkg else os.path.join(IJF, "server_upload_azdual")
OUT = os.path.join(IJF, "r3_tables"); NDRAW = A.ndraw; NB = max(1, NDRAW // 5)
REUSE_CONFIRM = [30, 33, 35, 73, 108, 116, 117, 118, 134, 446]
NAME = {"sv_az_sn": "SV-AZ-SN", "sv_az_st": "SV-AZ-ST", "rsv_az_sn": "RSV-AZ-SN", "rsv_az_st": "RSV-AZ-ST", "rsv_az_sn_psi": "RSV-AZ-SN (psi-free)", "rsv_az_st_psi": "RSV-AZ-ST (psi-free)"}
SVC = ["SV-N", "SV-T", "SV-AZ-SN", "SV-AZ-ST", "SV-FS-SN", "SV-FS-ST", "SV-GH-ST"]; RSVC = ["RSV-" + x[3:] for x in SVC]
EVALDIR = {"DJIA_main": (os.path.join(ROOT, "R", "prediction_results_2026"), "DJIA"), "N225_main": (os.path.join(ROOT, "R", "prediction_results_2026"), "N225"),
           "DJIA_2020": (os.path.join(IJF, "reproducibility", "R", "prediction_results_djia_2020"), "DJIA"), "N225_2020": (os.path.join(IJF, "reproducibility", "R", "prediction_results_n225_2020"), "N225")}
def f(x):
    x = x.strip(); return float("inf") if x in ("+.Inf", ".Inf") else float("-inf") if x == "-.Inf" else float("nan") if x in (".NaN", "NaN") else float(x)
def rows(p):
    with open(p) as fh:
        r = csv.reader(fh); next(r); return [[f(x) for x in row[1:]] for row in r if row]
def h1_col(p):
    with open(p) as fh:
        next(fh); return [float(line.split(",", 2)[1]) for line in fh if line.strip()]
man = list(csv.DictReader(open(os.path.join(PKG, "manifest.csv"))))
targets = {}
for r in man: targets.setdefault(r["target"], {"sample": r["sample"], "model": r["model"], "nwin": int(r["nwin"]), "dir": os.path.join(RES, r["target"])})
targets["djia_main_sv_az_st"] = {"sample": "DJIA_main", "model": "sv_az_st", "nwin": 603, "dir": os.path.join(IJF, "azfix_results", "fixed"), "reuse": True}
md = ["# AZ 再計算（dual）の新 run の MCMC 診断（§86）", "", f"RES = {RES}, ndraw = {NDRAW}. 規則はスクリプト冒頭（mcmc_main_diagnostics.py と同一）。F1: h_T Inef > 300、F2: Geweke p < 0.01（記録のみ）、F3: 先頭 {NB} ドロー / 残りの exp(h1) 平均比 > 2 または < 0.5、F4: |log(mean(exp h1) / 同クラス他モデル中央値)| > 1.0（入口、評価 CSV があれば）。確認パス対象 = F1 ∪ F3。", "",
      "| target | 評価窓/全窓 | F1 | F2 | F3 | F1∪F3 | F4 | (F1∪F3)∩F4 | h_T Inef 中央値 / q90 / max | ブロック比 中央値 / min / max | 確認パス対象の窓 |", "|---|---|---|---|---|---|---|---|---|---|---|"]
allw, confirm = [], []
for tg, info in targets.items():
    n, m, d = info["nwin"], info["model"], info["dir"]; name = NAME[m]
    evdir, idx = EVALDIR[info["sample"]]; ev = None
    for sfx in ("_r3corrected_dualfixed", "_r3corrected"):
        p = os.path.join(evdir + sfx, f"{idx}_forecast_vol.csv")
        if os.path.exists(p) and not m.endswith("_psi"): ev = list(csv.DictReader(open(p))); evsrc = sfx; break
    cls = SVC if name.startswith("SV") else RSVC; base = name.replace(" (psi-free)", ""); others = [c for c in cls if c != base]
    W = []
    for t in range(n):
        ph, pp = os.path.join(d, f"posstath_{m}_{t}.csv"), os.path.join(d, f"h1_y_{m}_{t}.csv")
        if not (os.path.exists(ph) and os.path.exists(pp)): continue
        last = rows(ph)[-1]; h1 = h1_col(pp)
        if len(h1) != NDRAW: continue
        e = [math.exp(x) for x in h1]; a = sum(e[:NB]) / NB; b = sum(e[NB:]) / (NDRAW - NB); ratio = a / b if b > 0 else float("nan")
        w = {"target": tg, "model": name, "t": t, "hT_mean": last[0], "hT_sd": last[1], "hT_gew": last[4], "hT_inef": last[5], "block_ratio": ratio, "mexp": sum(e) / NDRAW, "dev": float("nan")}
        if ev is not None and base in ev[0]:
            try: w["dev"] = math.log(w["mexp"] / st.median(float(ev[t][c]) for c in others if c in ev[t]))
            except Exception: pass
        w["F1"] = w["hT_inef"] > 300; w["F2"] = w["hT_gew"] < 0.01; w["F3"] = ratio > 2 or ratio < 0.5; w["F4"] = abs(w["dev"]) > 1.0 if not math.isnan(w["dev"]) else False
        w["confirm"] = w["F1"] or w["F3"]; W.append(w)
    allw += W
    if not W: md.append(f"| {tg} | 0/{n} | | | | | | | | | 未評価 |"); continue
    c13 = [w for w in W if w["confirm"]]; inef = sorted(w["hT_inef"] for w in W); br = sorted(w["block_ratio"] for w in W)
    cw = sorted(w["t"] for w in c13)
    if info.get("reuse"):
        note = f"再利用分。§85 の要確認 10 窓 {REUSE_CONFIRM} と{'一致' if cw == REUSE_CONFIRM else '不一致 → 要確認'}"
        for t in REUSE_CONFIRM:
            if t not in cw: cw.append(t)
        cw = sorted(cw)
    else: note = ""
    for t in cw: confirm.append([tg, name, t, "F1" if any(w["t"] == t and w["F1"] for w in W) else "", "F3" if any(w["t"] == t and w["F3"] for w in W) else "", "reuse-list" if info.get("reuse") and t in REUSE_CONFIRM else ""])
    q = lambda v, p: v[min(len(v) - 1, int(p * (len(v) - 1) + 0.5))]
    md.append(f"| {tg} | {len(W)}/{n} | {sum(w['F1'] for w in W)} | {sum(w['F2'] for w in W)} | {sum(w['F3'] for w in W)} | {len(c13)} | {sum(w['F4'] for w in W)} | {sum(w['F4'] for w in c13)} | {q(inef, .5):.0f} / {q(inef, .9):.0f} / {inef[-1]:.0f} | {q(br, .5):.2f} / {br[0]:.2f} / {br[-1]:.2f} | {cw if len(cw) <= 25 else str(cw[:25]) + '…'} {note} |")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "azdual_diagnostics_windows.csv"), "w", newline="") as fh:
    w = csv.writer(fh); keys = ["target", "model", "t", "hT_mean", "hT_sd", "hT_gew", "hT_inef", "block_ratio", "mexp", "dev", "F1", "F2", "F3", "F4", "confirm"]; w.writerow(keys); [w.writerow([x[k] for k in keys]) for x in allw]
with open(os.path.join(OUT, "azdual_confirm_windows.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["target", "model", "t", "F1", "F3", "note"]); w.writerows(confirm)
md += ["", f"評価窓 {len(allw)}、確認パス対象（F1∪F3 ＋ 再利用分の 10 窓）{len(confirm)} 窓 → `r3_tables/azdual_confirm_windows.csv`。F2 のみの窓は記録だけ（§84）。",
       "確認パス: burn-in 15000・独立 3 チェーン（初期値 3 通り、窓ごとの固定 seed）を dual で回し、§84 補遺 §4 の採用規則（チェーン間 h_T 一致 0.5 SD、Inef < 300、Geweke p > 0.01、ブロック比 0.5〜2）で判定。"]
open(os.path.join(OUT, "azdual_diagnostics.md"), "w").write("\n".join(md) + "\n"); print("\n".join(md))
