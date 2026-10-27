#!/usr/bin/env python3
"""確認パス（§87-89、規則 r3_snippets/azconfirm_rules_20261002.md §7-§8）の判定と併合。
入力: <pkg>/manifest_confirm.csv（予定チェーン集合: target, model, t, chain, seed, driver, stage, nsave, outdir。stage/nsave/outdir が無い旧形式は 1/15000/t<t>_c<k>）、
      <pkg>/Ox/out_confirm/<target>/<outdir>/ の pstate_<m>_<t>.csv（φ σ ρ μ δ ν (ξ ψ σ_u) h_T z0_T λ_T z01 λ1 e_h e_y (e_x)）、h1_y_<m>_<t>.csv（[h1 y1 (x1)] fixed | orig）、
      posstath_<m>_<t>.csv（最終行 = h_T: mean sd 2.5% 97.5% Geweke Inef）、<pkg>/Ox/logs/<driver>.done（完了の印）。
窓ごとの判定対象 = manifest にあるその窓の最大 stage のチェーン集合（段階 2 があれば段階 2 の 3 本だけ。段階 1 と加算しない）。
前提検査（満たさないと判定しない）: chain ID が一意、予定チェーン全部に .done と出力があり行数 = nsave（incomplete）、manifest 外の出力ディレクトリが無い（problem）、
  pstate・h1_y の全行全列が有限（モデル上の固定値 ν=+Inf（SN）・ψ=1（ψ 固定）は明示的に区別し、全ドローで期待値に等しいことを確認）、A の診断値（Inef・Geweke・mean・sd）が有限（invalid）。
判定（標準実装 tools/mcmc_diag.py = Vehtari et al. 2021 / R posterior と一致）:
  A（9/30）: 各チェーン h_T Inef < 300、Geweke p > 0.01、先頭ブロック比 0.5〜2、チェーン間 h_T 平均差 ≤ 0.5 SD。
  B 標準（10/2 13:00 の基準、正式判定）: rank-normalized／folded split R-hat < 1.01、bulk/tail ESS ≥ 400 を h_T・パラメータ・exp(h1)・y1 の全量に、
     MCSE（チェーン内バッチ平均、各チェーン 10 バッチ = 事前固定）: mean(exp h1) ≤ 0.5%、1%・5% の VaR・ES ≤ |値| の 1%。非有限の診断値は不合格。
  B 参考（10/2 13:30 の緩和版、記録のみ）。1% ES・VaR の MCSE はバッチ数 5/10/20/50 の感度も記録（判定は 10 のみ。合格する設定を選ばない）。
  status: accepted（前提検査 OK、チェーン ≥ --min-chains、A ∧ B 標準）、rejected、incomplete（チェーン不足・未完了）、invalid（非有限・形式不正）。
併合（--pool）: accepted の窓だけ、仮出力 → 検証 → 正式名（azconfirm_pool.py）。入力チェーンの md5 と併合ファイルの md5 を confirm_status.csv に結び付け、
  併合が完了・検証できた後に status = accepted / pooled = yes を確定する。accepted でない窓の古い併合出力は消す。
--verify-prefix: 段階 2 の各チェーンの先頭 nsave(段階 1) 行が段階 1 の pstate・h1_y と一致するか（同 seed・同初期値なら bitwise 一致のはず）を報告。
出力: <tables>/azconfirm_check.md、<tables>/azconfirm_check_windows.csv。
使い方: python3 tools/azconfirm_check.py [--pkg DIR] [--ndraw N] [--pool] [--smoke] [--min-chains 3] [--out-dir DIR] [--tables-dir DIR] [--only target:t,...] [--verify-prefix]
"""
import argparse, csv, math, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcmc_diag as D
from azconfirm_pool import pool_window, remove_pooled, read_status, write_status, now, ox_var_es, md5, fmt_md5s
ap = argparse.ArgumentParser()
ap.add_argument("--pkg", default=None); ap.add_argument("--ndraw", type=int, default=15000, help="段階 1 の保存ドロー数（manifest に nsave が無い旧形式の既定）"); ap.add_argument("--pool", action="store_true"); ap.add_argument("--smoke", action="store_true")
ap.add_argument("--min-chains", type=int, default=3); ap.add_argument("--out-dir", default=None); ap.add_argument("--tables-dir", default=None); ap.add_argument("--only", default=None); ap.add_argument("--verify-prefix", action="store_true")
A = ap.parse_args()
IJF = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PKG = os.path.abspath(A.pkg) if A.pkg else os.path.join(IJF, "server_upload_azconfirm")
OUTD = os.path.abspath(A.out_dir) if A.out_dir else os.path.join(IJF, "azdual_results", "confirm"); TAB = os.path.abspath(A.tables_dir) if A.tables_dir else os.path.join(IJF, "r3_tables")
MINCH = A.min_chains; NB_MCSE = 10; NB_SENS = (5, 10, 20, 50)
RHAT_MAX, ESS_MIN, MCSE_MEAN, MCSE_VAR = 1.01, 400, 0.005, 0.01            # B 標準（10/2 13:00）
KEY = ("h_T", "exp_h1", "y1"); RHAT_PAR_REF, ESS_PAR_REF, MCSE_VAR1_REF = 1.05, 100, 0.02  # B 参考（10/2 13:30、記録のみ）
try: CODE = subprocess.run(["git", "-C", IJF, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip() or "n/a"
except Exception: CODE = "n/a"
def f(x):
    x = x.strip(); return float("inf") if x in ("+.Inf", ".Inf", "inf") else float("-inf") if x in ("-.Inf", "-inf") else float("nan") if x in (".NaN", "NaN", "nan", "") else float(x)
def mat(p):
    with open(p) as fh:
        r = csv.reader(fh); next(r); return [[f(x) for x in row[1:]] for row in r if row]
def fmt(x, d=3): return "nan" if (isinstance(x, float) and (math.isnan(x) or math.isinf(x))) else f"{x:.{d}f}"
def fin(x): return isinstance(x, (int, float)) and math.isfinite(x)
# ---- manifest ----
man = list(csv.DictReader(open(os.path.join(PKG, "manifest_confirm.csv"))))
for r in man:
    r["stage"] = int(r.get("stage") or 1); r["nsave"] = int(r.get("nsave") or A.ndraw); r["outdir"] = r.get("outdir") or f"t{r['t']}_c{r['chain']}"; r["chain"] = int(r["chain"]); r["t"] = int(r["t"])
windows = {}
for r in man: windows.setdefault((r["target"], r["model"], r["t"]), []).append(r)
only = set(A.only.split(",")) if A.only else None
md = ["# 確認パスの判定（§87-89）", "", f"pkg = {PKG}, min chains = {MINCH}, code = {CODE}, {now()}",
      "規則: r3_snippets/azconfirm_rules_20261002.md §7-§8（A: 9/30。B 標準 = 10/2 13:00 の基準を標準実装（rank-normalized／folded split R-hat、bulk／tail ESS、Vehtari et al. 2021）で判定。",
      f"B 参考 = 10/2 13:30 の緩和版、正式判定には使わない）。MCSE はチェーン内バッチ平均（各チェーン {NB_MCSE} バッチ、事前固定。感度 {NB_SENS} は記録のみ）。判定対象 = その窓の最大 stage のチェーン集合。", "",
      "| target | t | stage × nsave | chains | status | A: Inef max / Geweke min / block ratio / Δh_T (SD) | B 標準: R-hat max (量) / bulk ESS min (量) / tail ESS min (量) | MCSE rel: mean(exp h) / VaR1 / VaR5 / ES1 / ES5 | 不合格・前提検査 | B 参考 |",
      "|---|---|---|---|---|---|---|---|---|---|"]
rows_out = []; status_path = os.path.join(OUTD, "confirm_status.csv"); status = read_status(status_path) if A.pool else {}; problems = []
def outdir_of(r): return os.path.join(PKG, "Ox", "out_confirm", r["target"], r["outdir"] + ("_smoke" if A.smoke else ""))
# manifest 外の出力ディレクトリ
known = {(r["target"], r["outdir"]) for r in man}
oc = os.path.join(PKG, "Ox", "out_confirm")
if os.path.isdir(oc):
    for tg in sorted(os.listdir(oc)):
        for d in sorted(os.listdir(os.path.join(oc, tg))):
            if d.endswith("_smoke") or d.startswith("_tmp_"): continue
            if (tg, d) not in known: problems.append(f"output directory not in manifest: out_confirm/{tg}/{d}")
for (tg, model, t), rs in sorted(windows.items()):
    if only and f"{tg}:{t}" not in only: continue
    stage = max(r["stage"] for r in rs); plan = [r for r in rs if r["stage"] == stage]; nsave = 300 if A.smoke else plan[0]["nsave"]; NB = max(1, nsave // 5)  # smoke 出力は常に 300 ドロー
    is_rsv = model.startswith("rsv"); is_psi = model.endswith("_psi"); is_st = model.split("_")[2] == "st"
    pnames = ["phi", "sigma", "rho", "mu", "delta"] + (["nu"] if is_st else []) + ((["xi"] + (["psi"] if is_psi else []) + ["sigma_u"]) if is_rsv else [])
    pidx = {"phi": 0, "sigma": 1, "rho": 2, "mu": 3, "delta": 4, "nu": 5, "xi": 6, "psi": 7, "sigma_u": 8}; hT_idx = 9 if is_rsv else 6; ncol = 6 if is_rsv else 4; npcol = 17 if is_rsv else 13
    fixed = {}  # モデル上の固定値（列 → 期待値）。SN: ν = +Inf。ψ 固定（RSV, 非 _psi）: ψ = 1
    if not is_st: fixed[5] = float("inf")
    if is_rsv and not is_psi: fixed[7] = 1.0
    pre = []  # 前提検査の不備
    ids = [r["chain"] for r in plan]
    if len(set(ids)) != len(ids): pre.append(f"duplicate chain id in manifest {ids}")
    if any(len(set(r["nsave"] for r in plan)) != 1 for _ in [0]): pre.append("nsave differs within the planned set")
    PS, H, DIAG, used, in_md5, missing = [], [], [], [], {}, []
    for r in sorted(plan, key=lambda r: r["chain"]):
        k = r["chain"]; d = outdir_of(r); pp = os.path.join(d, f"pstate_{model}_{t}.csv"); hp = os.path.join(d, f"h1_y_{model}_{t}.csv"); sp = os.path.join(d, f"posstath_{model}_{t}.csv")
        done = os.path.exists(os.path.join(PKG, "Ox", "logs", r["driver"][:-3] + ".done"))
        if not (os.path.exists(pp) and os.path.exists(hp) and os.path.exists(sp)): missing.append(f"c{k}:no output"); continue
        if not done and not A.smoke: missing.append(f"c{k}:not .done"); continue
        if any(os.path.getsize(x) == 0 for x in (pp, hp, sp)): missing.append(f"c{k}:empty output file (write failure, e.g. disk full)"); continue  # 2026-10-09: t446_c2_s2 on sakura
        ps = mat(pp); hy = mat(hp); last = mat(sp)[-1]
        if len(ps) != nsave or len(hy) != nsave: pre.append(f"c{k}: rows pstate {len(ps)} / h1_y {len(hy)} != nsave {nsave}"); continue
        if any(len(row) != npcol for row in ps) or any(len(row) != ncol for row in hy): pre.append(f"c{k}: column count (pstate {len(ps[0])}, h1_y {len(hy[0])})"); continue
        bad = sorted({j for row in ps for j in range(npcol) if j not in fixed and not math.isfinite(row[j])}); badh = sorted({j for row in hy for j in range(ncol) if not math.isfinite(row[j])})
        if bad: pre.append(f"c{k}: non-finite pstate columns {bad}")
        if badh: pre.append(f"c{k}: non-finite h1_y columns {badh}")
        for j, v in fixed.items():
            if not all(row[j] == v for row in ps): pre.append(f"c{k}: fixed column {j} not constant at {v}")
        if not all(fin(x) for x in (last[0], last[1], last[4], last[5])) or not last[1] > 0: pre.append(f"c{k}: non-finite A diagnostics (mean {last[0]}, sd {last[1]}, Geweke {last[4]}, Inef {last[5]})")
        if bad or badh: continue
        e = [math.exp(row[0]) for row in hy]
        PS.append(ps); H.append(hy); used.append(k); in_md5[f"c{k}"] = ",".join(md5(p) for p in (pp, hp, sp))
        DIAG.append({"inef": last[5], "gew": last[4], "hT_mean": last[0], "hT_sd": last[1], "block": (sum(e[:NB]) / NB) / (sum(e[NB:]) / (nsave - NB))})
    nch = len(PS)
    # 段階 2 の先頭一致（記録のみ、判定には使わない）
    prefix_note = ""
    if A.verify_prefix and stage >= 2:
        notes = []
        for r in sorted(plan, key=lambda r: r["chain"]):
            r1 = [x for x in rs if x["stage"] == 1 and x["chain"] == r["chain"]]
            if not r1: notes.append(f"c{r['chain']}: no stage-1 row"); continue
            d1, d2 = outdir_of(r1[0]), outdir_of(r)
            res = []
            for kind in ("pstate", "h1_y"):
                p1, p2 = os.path.join(d1, f"{kind}_{model}_{t}.csv"), os.path.join(d2, f"{kind}_{model}_{t}.csv")
                if not (os.path.exists(p1) and os.path.exists(p2)): res.append(f"{kind}: n/a"); continue
                with open(p1) as f1, open(p2) as f2:
                    l1 = f1.readlines(); l2 = f2.readlines()
                n1 = r1[0]["nsave"] + 1; diff = next((i for i in range(min(n1, len(l1), len(l2))) if l1[i] != l2[i]), None)
                res.append(f"{kind}: {'identical first ' + str(n1 - 1) if diff is None and len(l2) >= n1 and len(l1) >= n1 else 'DIFFERS at row ' + str(diff) if diff is not None else 'too short'}")
            notes.append(f"c{r['chain']} " + "; ".join(res))
        prefix_note = " / ".join(notes); print(f"{tg} t={t} prefix check: {prefix_note}")
    base = {"target": tg, "model": model, "t": t, "stage": stage, "nsave": nsave, "chains": nch, "chain_ids": " ".join(map(str, used)), "ndraw_per_chain": nsave, "ndraw_total": nch * nsave, "prefix_check": prefix_note, "fixed_columns": " ".join(f"{j}={v}" for j, v in fixed.items())}
    if pre or missing or nch < MINCH or len(plan) < MINCH:
        st_ = "invalid" if pre else "incomplete"; note = "; ".join(pre + missing) or f"{nch}/{len(plan)} chains"
        md.append(f"| {tg} | {t} | {stage} × {nsave} | {nch}/{len(plan)} | {st_} | | | | {note} | |"); rows_out.append(dict(base, status=st_, verdict=note, fail_items=note))
        if A.pool:
            od = os.path.join(OUTD, tg); gone = remove_pooled(od, model, t)
            if gone: print(f"{tg} t={t}: status {st_}; removed stale pooled output {gone}")
            status[(tg, model, t)] = dict(base, status=st_, verdict=note, okA="", okB_std="", okB_ref="", pooled="no", input_md5=fmt_md5s(in_md5), pooled_md5="", checked_at=now(), code=CODE)
        if st_ == "invalid": problems.append(f"{tg} t={t}: {note}")
        continue
    # ---- A ----
    inef_max = max(dg["inef"] for dg in DIAG); gew_min = min(dg["gew"] for dg in DIAG); br = [dg["block"] for dg in DIAG]
    sdp = math.sqrt(sum(dg["hT_sd"] ** 2 for dg in DIAG) / nch); dh = max(abs(a["hT_mean"] - b["hT_mean"]) for a in DIAG for b in DIAG) / sdp
    failA = [s for s, ok in (("Inef", inef_max < 300), ("Geweke", gew_min > 0.01), ("block", all(0.5 <= x <= 2 for x in br)), ("ΔhT", dh <= 0.5)) if not ok]; okA = not failA
    # ---- B: R-hat / ESS ----
    qty = {nm: [[row[pidx[nm]] for row in ps] for ps in PS] for nm in pnames}; qty["h_T"] = [[row[hT_idx] for row in ps] for ps in PS]
    qty["exp_h1"] = [[math.exp(row[0]) for row in hy] for hy in H]; qty["y1"] = [[row[1] for row in hy] for hy in H]
    rh = {nm: D.rhat(ch) for nm, ch in qty.items()}; be = {nm: D.ess_bulk(ch) for nm, ch in qty.items()}; te = {nm: D.ess_tail(ch) for nm, ch in qty.items()}
    amax = lambda dd: max(dd, key=lambda k: (float("inf") if not fin(dd[k]) else dd[k])); amin = lambda dd: min(dd, key=lambda k: (-1 if not fin(dd[k]) else dd[k]))
    rh_arg, be_arg, te_arg = amax(rh), amin(be), amin(te); rh_max, be_min, te_min = rh[rh_arg], be[be_arg], te[te_arg]
    # ---- B: MCSE（チェーン内バッチ、fixed と orig、バッチ数感度）----
    def mcse_set(ecols, ycols, nb):
        m_e = sum(v for c in ecols for v in c) / sum(len(c) for c in ecols); out = {"mean_exp_h": D.batch_mcse_within(ecols, lambda b: sum(b) / len(b), nb) / m_e, "mean_exp_h_ess": D.mcse_mean(ecols) / m_e}
        y_all = [v for c in ycols for v in c]; var, es = ox_var_es(y_all)
        for i, lab, p in ((1, "1", 0.01), (2, "5", 0.05)):
            out["VaR" + lab] = D.batch_mcse_within(ycols, lambda b: D.quantile7(b, p), nb) / abs(var[i]); out["ES" + lab] = D.batch_mcse_within(ycols, lambda b: ox_var_es(b)[1][i], nb) / abs(es[i])
        return out
    eo = [[math.exp(row[ncol // 2]) for row in hy] for hy in H]; yo = [[row[ncol // 2 + 1] for row in hy] for hy in H]
    mc = mcse_set(qty["exp_h1"], qty["y1"], NB_MCSE); mco = mcse_set(eo, yo, NB_MCSE); sens = {nb: mcse_set(qty["exp_h1"], qty["y1"], nb) for nb in NB_SENS}
    failB = []
    bad_rh = [q for q in qty if not (fin(rh[q]) and rh[q] < RHAT_MAX)]; bad_be = [q for q in qty if not (fin(be[q]) and be[q] >= ESS_MIN)]; bad_te = [q for q in qty if not (fin(te[q]) and te[q] >= ESS_MIN)]
    if bad_rh: failB.append("R-hat(" + ",".join(bad_rh) + ")")
    if bad_be: failB.append("bulkESS(" + ",".join(bad_be) + ")")
    if bad_te: failB.append("tailESS(" + ",".join(bad_te) + ")")
    if not (fin(mc["mean_exp_h"]) and mc["mean_exp_h"] <= MCSE_MEAN): failB.append("MCSE mean(exp h)")
    for k in ("VaR1", "VaR5", "ES1", "ES5"):
        if not (fin(mc[k]) and mc[k] <= MCSE_VAR): failB.append("MCSE " + k)
    okB = not failB
    okB_ref = all(fin(rh[q]) and rh[q] < RHAT_MAX and fin(be[q]) and be[q] >= ESS_MIN and fin(te[q]) and te[q] >= ESS_MIN for q in KEY) and all(fin(rh[q]) and rh[q] < RHAT_PAR_REF and fin(be[q]) and fin(te[q]) and min(be[q], te[q]) >= ESS_PAR_REF for q in pnames) \
        and all(fin(mc[k]) for k in mc) and mc["mean_exp_h"] <= MCSE_MEAN and mc["VaR5"] <= MCSE_VAR and mc["ES5"] <= MCSE_VAR and mc["VaR1"] <= MCSE_VAR1_REF and mc["ES1"] <= MCSE_VAR1_REF
    st_ = "accepted" if (okA and okB) else "rejected"
    verdict = ("採用（併合）" if st_ == "accepted" else "不合格") + ("" if okA else "；A 不合格") + ("" if okB else "；B 標準不合格")
    md.append(f"| {tg} | {t} | {stage} × {nsave} | {nch}/{len(plan)} | {st_} | {inef_max:.0f} / {gew_min:.3f} / {min(br):.2f}-{max(br):.2f} / {dh:.2f} | {fmt(rh_max)} ({rh_arg}) / {fmt(be_min, 0)} ({be_arg}) / {fmt(te_min, 0)} ({te_arg}) | "
              f"{mc['mean_exp_h']*100:.2f}% / {mc['VaR1']*100:.2f}% / {mc['VaR5']*100:.2f}% / {mc['ES1']*100:.2f}% / {mc['ES5']*100:.2f}% | {'; '.join(failA + failB) or '—'} | {'合格' if okB_ref else '不合格'} |")
    row = dict(base, status=st_, verdict=verdict, okA=okA, okB_std=okB, okB_ref=okB_ref, inef_max=inef_max, geweke_min=gew_min, block_min=min(br), block_max=max(br), dhT_sd=dh, rhat_max=rh_max, rhat_arg=rh_arg, bulk_ess_min=be_min, bulk_ess_arg=be_arg,
               tail_ess_min=te_min, tail_ess_arg=te_arg, mcse_mean_rel=mc["mean_exp_h"], mcse_mean_rel_ess=mc["mean_exp_h_ess"], mcse_var1=mc["VaR1"], mcse_var5=mc["VaR5"], mcse_es1=mc["ES1"], mcse_es5=mc["ES5"],
               mcse_mean_rel_orig=mco["mean_exp_h"], mcse_var1_orig=mco["VaR1"], mcse_var5_orig=mco["VaR5"], mcse_es1_orig=mco["ES1"], mcse_es5_orig=mco["ES5"], fail_items="; ".join(failA + failB))
    for nb in NB_SENS: row[f"mcse_es1_nb{nb}"] = sens[nb]["ES1"]; row[f"mcse_var1_nb{nb}"] = sens[nb]["VaR1"]
    for q in qty: row[f"rhat_{q}"] = rh[q]; row[f"bulk_ess_{q}"] = be[q]; row[f"tail_ess_{q}"] = te[q]
    rows_out.append(row)
    if A.pool:
        od = os.path.join(OUTD, tg); pooled = "no"; pmd5 = {}
        if st_ == "accepted":
            try: pmd5 = pool_window(od, model, t, PS, H, ncol, pnames, pidx, hT_idx); pooled = "yes"; print(f"{tg} t={t}: pooled {nch} chains × {nsave} -> {od} ({len(pmd5)} files, verified)")
            except Exception as ex: problems.append(f"{tg} t={t}: pooling failed: {ex}"); st_ = "rejected"; verdict += f"；併合失敗 {ex}"; remove_pooled(od, model, t)
        else:
            gone = remove_pooled(od, model, t)
            if gone: print(f"{tg} t={t}: status {st_}; removed stale pooled output {gone}")
        status[(tg, model, t)] = dict(base, status=st_, verdict=verdict, okA=okA, okB_std=okB, okB_ref=okB_ref, pooled=pooled, input_md5=fmt_md5s(in_md5), pooled_md5=fmt_md5s(pmd5), checked_at=now(), code=CODE)
if A.pool: write_status(status_path, status); print(f"status -> {status_path}")
os.makedirs(TAB, exist_ok=True)
cols = ["target", "model", "t", "stage", "nsave", "chains", "chain_ids", "status", "verdict", "okA", "okB_std", "okB_ref", "inef_max", "geweke_min", "block_min", "block_max", "dhT_sd", "rhat_max", "rhat_arg", "bulk_ess_min", "bulk_ess_arg", "tail_ess_min", "tail_ess_arg",
        "mcse_mean_rel", "mcse_mean_rel_ess", "mcse_var1", "mcse_var5", "mcse_es1", "mcse_es5", "mcse_mean_rel_orig", "mcse_var1_orig", "mcse_var5_orig", "mcse_es1_orig", "mcse_es5_orig", "fail_items", "prefix_check", "fixed_columns"]
extra = sorted({k for r in rows_out for k in r if k not in cols})
with open(os.path.join(TAB, "azconfirm_check_windows.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols + extra); w.writeheader(); w.writerows(rows_out)
md += ["", "status: accepted = 前提検査 OK・チェーン ≥ min chains・A ∧ B 標準（併合は仮出力の検証後に確定）。rejected = 不合格（併合出力は作らない）。incomplete = チェーン不足・未完了。invalid = 非有限値・形式不正・manifest 不整合。",
       "量ごとの R-hat／ESS は azconfirm_check_windows.csv の rhat_*／bulk_ess_*／tail_ess_* 列。orig 側（旧式の予測）の MCSE は *_orig 列、バッチ数感度は mcse_es1_nb*／mcse_var1_nb* 列（判定には使わない）。"]
if problems: md += ["", "PROBLEMS:"] + [f"- {p}" for p in problems]
open(os.path.join(TAB, "azconfirm_check.md"), "w").write("\n".join(md) + "\n"); print("\n".join(md))
sys.exit(1 if problems else 0)
