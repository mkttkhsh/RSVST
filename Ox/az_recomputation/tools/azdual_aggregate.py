#!/usr/bin/env python3
"""AZ 系 20 対象の修正再計算（dual 版、PROGRESS.md §85-86）の窓別出力 → R の merge スクリプトが読む集約ファイル。
入力: azdual_results/<target>/（tools/azdual_collect_check.py で統合済み・検証済みの窓別 12 種）
      再利用分 DJIA 主標本 SV-AZ-ST: azfix_results/fixed（修正式）と azfix_results/orig（旧式、同一 seed・同一 MCMC、§83）
      確認パス（--confirm-dir、既定 azdual_results/confirm）: confirm_status.csv で status == accepted かつ pooled == yes の窓だけ、併合出力（fixed・orig とも同じ採用済みドロー）で差し替える。
        併合ファイルの md5 を状態ファイルの pooled_md5 と照合する。併合出力があるのに accepted でない窓、accepted なのに併合出力が無い窓、行数が ndraw_total と違う窓、md5 不一致は問題。
      確認対象の一覧（--confirm-manifest、既定 server_upload_azconfirm/manifest_confirm.csv）: 確認対象なのに accepted でない窓の値は run の値（単一チェーン）を**暫定**として入れ、
        windows_<model>.csv の confirm_status 列（accepted / rejected / incomplete / invalid / not_run）と PROVISIONAL_<model>.txt に記録する。暫定窓が 1 つでもあれば、その集約は正式確定ではない。
出力: azdual_results/aggregated/<variant>/<sdir>/{volf,var,es}_<model>.csv   variant = dualfixed（修正式、正式評価）| dualorig（旧式、同じドローから）
      sdir = DJIA_prediction | N225_prediction | DJIA_prediction_2020 | N225_prediction_2020（論文掲載 run の集約ファイルと同じ場所名）
      列と行名の書式は論文掲載 run の同名ファイル（Ox2022/<sdir>/、reproducibility/Ox/<sdir>/）のヘッダをそのまま写す（R 側は無変更で読める）:
        volf = mean(h), median(h), mean(exp h), median(exp h) = exp(median h)   var/es = 0.5%, 1%, 5%, 10% 分位（1% = Var2、5% = Var3）
      あわせて windows_<model>.csv（t, date, mean_h, median_h, mean_exp_h, var1, var5, es1, es5, source, ndraw, confirm_status）、PROVISIONAL_<model>.txt（暫定窓があるとき）、
      manifest.csv（対象・窓数・欠損・確認パス採用数・暫定数・入力元・md5）
行の対応: 窓 t の行は常に評価 CSV の t 行目の日付に対応する（行数 = 窓数）。欠けた窓・読めない窓は NaN 行（NA）。窓を飛ばして行をずらすことはしない。
検査: 各窓の volf の mean(h) と h1_y から計算した平均が一致（1e-8）、median(h) との差（1e-6）、非有限値なし、ドロー数（run は ndraw、確認パスは ndraw_total）、全窓そろっていること
      （--allow-partial で欠損窓を NaN 埋めして PARTIAL.txt に記録。PARTIAL.txt の他対象の行は消さない）。
**検査に 1 件でも問題があれば、集約ファイルは一切書かない（既存の出力は変更されない）。** 問題が無いときだけ全部を書く。
使い方: python3 tools/azdual_aggregate.py [--res DIR] [--ndraw N] [--allow-partial] [--targets a,b,...] [--confirm-dir DIR] [--no-confirm] [--confirm-manifest CSV]
点予測は mean(exp h)（§75 の mean 統一）。median 列は書式維持のために出すだけで原稿には使わない。
"""
import argparse, csv, hashlib, math, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from azconfirm_pool import parse_md5s, md5 as file_md5
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None); ap.add_argument("--ndraw", type=int, default=15000); ap.add_argument("--allow-partial", action="store_true")
ap.add_argument("--targets", default=None); ap.add_argument("--pkg", default=None); ap.add_argument("--azfix", default=None)
ap.add_argument("--confirm-dir", default=None, help="確認パスの併合出力（既定 azdual_results/confirm）。confirm_status.csv で accepted かつ pooled の窓だけ差し替える")
ap.add_argument("--no-confirm", action="store_true", help="確認パスの差し替えをしない（run の値のみ。確認対象の窓はすべて暫定）")
ap.add_argument("--confirm-manifest", default=None, help="確認対象の一覧（既定 server_upload_azconfirm/manifest_confirm.csv。無ければ確認対象なしとして扱う）")
A = ap.parse_args()
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); IJF = os.path.join(ROOT, "20260806IJF_revision")
RES = os.path.abspath(A.res) if A.res else os.path.join(IJF, "azdual_results"); PKG = os.path.abspath(A.pkg) if A.pkg else os.path.join(IJF, "server_upload_azdual")
AZFIX = os.path.abspath(A.azfix) if A.azfix else os.path.join(IJF, "azfix_results"); NDRAW = A.ndraw
CONF = None if A.no_confirm else (os.path.abspath(A.confirm_dir) if A.confirm_dir else os.path.join(RES, "confirm"))
CMAN = os.path.abspath(A.confirm_manifest) if A.confirm_manifest else os.path.join(IJF, "server_upload_azconfirm", "manifest_confirm.csv")
SDIR = {"DJIA_main": "DJIA_prediction", "N225_main": "N225_prediction", "DJIA_2020": "DJIA_prediction_2020", "N225_2020": "N225_prediction_2020"}
PAPER = {"DJIA_main": os.path.join(ROOT, "Ox2022", "DJIA_prediction"), "N225_main": os.path.join(ROOT, "Ox2022", "N225_prediction"),
         "DJIA_2020": os.path.join(IJF, "reproducibility", "Ox", "DJIA_prediction_2020"), "N225_2020": os.path.join(IJF, "reproducibility", "Ox", "N225_prediction_2020")}
EVAL = {"DJIA_main": os.path.join(ROOT, "R", "prediction_results_2026_r3corrected", "DJIA_forecast_vol.csv"), "N225_main": os.path.join(ROOT, "R", "prediction_results_2026_r3corrected", "N225_forecast_vol.csv"),
        "DJIA_2020": os.path.join(IJF, "reproducibility", "R", "prediction_results_djia_2020_r3corrected", "DJIA_forecast_vol.csv"), "N225_2020": os.path.join(IJF, "reproducibility", "R", "prediction_results_n225_2020_r3corrected", "N225_forecast_vol.csv")}
def fnum(x):
    x = x.strip(); return float("inf") if x in ("+.Inf", ".Inf") else float("-inf") if x == "-.Inf" else float("nan") if x in (".NaN", "NaN", "NA", "") else float(x)
def rows(p):
    with open(p) as f:
        r = csv.reader(f); next(r); return [[fnum(x) for x in row[1:]] for row in r if row]
def header_style(p):  # 論文掲載 run のファイルのヘッダ行と行名の様式を写す
    with open(p) as f:
        hdr = f.readline().rstrip("\n"); first = f.readline().rstrip("\n")
    rn = None
    if hdr.startswith(","): rn = "t-1" if first.split(",")[0].strip('"').endswith("-1") else "int"
    return hdr, rn
NAN4 = [math.nan] * 4
man = list(csv.DictReader(open(os.path.join(PKG, "manifest.csv"))))
targets = {}
for r in man: targets.setdefault(r["target"], {"sample": r["sample"], "model": r["model"], "nwin": int(r["nwin"])})
targets["djia_main_sv_az_st"] = {"sample": "DJIA_main", "model": "sv_az_st", "nwin": 603, "reuse": True}  # 再利用（azfix_results）
sel = set(A.targets.split(",")) if A.targets else set(targets)
# 確認パスの状態と確認対象
conf_status = {}
if CONF:
    sp = os.path.join(CONF, "confirm_status.csv")
    if os.path.exists(sp):
        for r in csv.DictReader(open(sp)): conf_status[(r["target"], r["model"], int(r["t"]))] = r
planned = set()
if os.path.exists(CMAN):
    for r in csv.DictReader(open(CMAN)): planned.add((r["target"], r["model"], int(r["t"])))
else: print(f"WARNING: confirm manifest not found ({CMAN}); no window is treated as a confirmation target")
problems, manifest, outputs, partial_notes, prov_notes = [], [], [], {}, {}  # outputs: (path, text) を全部ためて最後に書く
for variant in ("dualfixed", "dualorig"):
    for tg, info in targets.items():
        if tg not in sel: continue
        n, model, sample = info["nwin"], info["model"], info["sample"]; is_rsv = model.startswith("rsv"); ncol = 6 if is_rsv else 4
        if info.get("reuse"):
            src = os.path.join(AZFIX, "fixed" if variant == "dualfixed" else "orig"); hcol_run = 0; kinds_run = ("volf", "var", "es")
        else:
            src = os.path.join(RES, tg); hcol_run = 0 if variant == "dualfixed" else ncol // 2
            kinds_run = ("volf", "var", "es") if variant == "dualfixed" else ("volf_orig", "var_orig", "es_orig")
        hcol_conf = 0 if variant == "dualfixed" else ncol // 2; kinds_conf = ("volf", "var", "es") if variant == "dualfixed" else ("volf_orig", "var_orig", "es_orig")
        ev = list(csv.DictReader(open(EVAL[sample])))
        if len(ev) != n: problems.append((variant, tg, f"evaluation CSV rows {len(ev)} != nwin {n}")); continue
        cdir = os.path.join(CONF, tg) if CONF else None
        out = os.path.join(RES, "aggregated", variant, SDIR[sample])
        volf, var, es, win, miss, nconf, prov = [], [], [], [], [], 0, []
        for t in range(n):
            date = ev[t]["date"]; key = (tg, model, t); st = conf_status.get(key); conf_file = cdir and os.path.exists(os.path.join(cdir, f"h1_y_{model}_{t}.csv"))
            use_conf = False
            if cdir:
                accepted = st is not None and st["status"] == "accepted" and st["pooled"] == "yes"
                if accepted and not conf_file: problems.append((variant, tg, t, "confirm status accepted/pooled but pooled h1_y missing"))
                elif conf_file and not accepted: problems.append((variant, tg, t, f"pooled output exists but status is {st['status'] + '/' + st['pooled'] if st else 'absent'}"))
                use_conf = accepted and conf_file
                if use_conf:  # 併合ファイルの md5 を採用記録と照合
                    rec = parse_md5s(st.get("pooled_md5", ""))
                    for fn in (f"h1_y_{model}_{t}.csv",) + tuple(f"{k}_{model}_{t}.csv" for k in kinds_conf):
                        if fn not in rec: problems.append((variant, tg, t, f"pooled_md5 record missing for {fn}")); use_conf = False
                        elif file_md5(os.path.join(cdir, fn)) != rec[fn]: problems.append((variant, tg, t, f"pooled file md5 mismatch: {fn}")); use_conf = False
            cstat = "" if key not in planned else ("accepted" if use_conf else (st["status"] if st else "not_run"))
            if cstat and cstat != "accepted": prov.append((t, cstat))
            def nanrow(src_label, nd):
                volf.append(list(NAN4)); var.append(list(NAN4)); es.append(list(NAN4)); win.append([t, date] + [math.nan] * 7 + [src_label, nd, cstat])
            if use_conf: sdir_t, hcol, kinds, need, label = cdir, hcol_conf, kinds_conf, int(st["ndraw_total"]), "confirm"
            else: sdir_t, hcol, kinds, need, label = src, hcol_run, kinds_run, NDRAW, "run"
            ph = os.path.join(sdir_t, f"h1_y_{model}_{t}.csv")
            if not os.path.exists(ph): miss.append(t); nanrow("missing", 0); continue
            try:
                h = rows(ph); h1 = [r[hcol] for r in h]
                vf = rows(os.path.join(sdir_t, f"{kinds[0]}_{model}_{t}.csv"))[0]; va = rows(os.path.join(sdir_t, f"{kinds[1]}_{model}_{t}.csv"))[0]; e = rows(os.path.join(sdir_t, f"{kinds[2]}_{model}_{t}.csv"))[0]
            except Exception as ex: problems.append((variant, tg, t, f"read error {ex}")); nanrow(label, 0); continue
            nd = len(h1)
            if nd != need or not all(math.isfinite(x) for x in h1): problems.append((variant, tg, t, f"h1_y rows {nd} (expected {need}, {label}) or non-finite")); nanrow(label, nd); continue
            mh = sum(h1) / nd; mexp = sum(math.exp(x) for x in h1) / nd; medh = statistics.median(h1)
            if abs(mh - vf[0]) > 1e-8: problems.append((variant, tg, t, f"volf mean(h) {vf[0]} != {mh}"))
            if abs(medh - vf[1]) > 1e-6: problems.append((variant, tg, t, f"volf median(h) {vf[1]} != {medh}"))
            if len(va) != 4 or len(e) != 4 or not all(math.isfinite(x) for x in va + e): problems.append((variant, tg, t, "var/es not 1x4 finite")); nanrow(label, nd); continue
            volf.append([mh, medh, mexp, math.exp(medh)]); var.append(list(va)); es.append(list(e)); win.append([t, date, mh, medh, mexp, va[1], va[2], e[1], e[2], label, nd, cstat]); nconf += use_conf
        assert len(volf) == len(var) == len(es) == len(win) == n, (tg, variant, len(volf), n)
        if miss:
            if A.allow_partial: partial_notes.setdefault(out, []).append(f"{variant} {tg} {model}: {len(miss)} windows missing {miss[:20]}\n"); print("WARNING partial:", variant, tg, f"{len(miss)} windows missing (e.g. {miss[:6]})")
            else: problems.append((variant, tg, f"{len(miss)} windows missing (e.g. {miss[:6]})")); continue
        for kind, m in (("volf", volf), ("var", var), ("es", es)):
            ref = os.path.join(PAPER[sample], f"{kind}_{model}.csv")
            if not os.path.exists(ref): problems.append((variant, tg, f"paper-run file for the header style not found: {ref}")); continue
            hdr, rn = header_style(ref); fn = os.path.join(out, f"{kind}_{model}.csv"); lines = [hdr + "\n"]
            for t, row in enumerate(m):
                vals = ",".join("NA" if (isinstance(x, float) and math.isnan(x)) else repr(float(x)) for x in row)
                lines.append((f'"{t+1}-1",' if rn == "t-1" else f"{t+1}," if rn == "int" else "") + vals + "\n")
            outputs.append((fn, "".join(lines))); manifest.append([variant, tg, model, SDIR[sample], kind, n, len(miss), nconf, len(prov), "provisional" if prov else ("final" if not miss else "partial"), src, os.path.relpath(fn, RES), hashlib.md5("".join(lines).encode()).hexdigest(), hdr])
        wl = ["t,date,mean_h,median_h,mean_exp_h,var1,var5,es1,es5,source,ndraw,confirm_status\n"] + [",".join("" if (isinstance(x, float) and math.isnan(x)) else str(x) for x in r) + "\n" for r in win]
        outputs.append((os.path.join(out, f"windows_{model}.csv"), "".join(wl)))
        pf = os.path.join(out, f"PROVISIONAL_{model}.txt")
        if prov: prov_notes[pf] = (f"{variant} {tg} {model}: {len(prov)} confirmation-target windows are NOT accepted; their run values (single chain) are provisional, not final:\n" + "".join(f"  t={t} {c}\n" for t, c in prov))
        else: prov_notes[pf] = None
        print(f"{variant} {tg}: {n - len(miss)}/{n} windows ({nconf} from the confirmation pass, {len(prov)} provisional) -> {out}")
print("\nPROBLEMS:", "none" if not problems else len(problems)); [print(" -", p) for p in problems[:40]]
if problems: print("nothing written (existing aggregated files unchanged)"); sys.exit(1)
# ---- 問題なし: ここで初めて書く ----
for fn, text in outputs:
    os.makedirs(os.path.dirname(fn), exist_ok=True); open(fn, "w", newline="").write(text)
for v in ("dualfixed", "dualorig"):  # PARTIAL.txt: 今回の対象の行だけ置き換える（他対象の未完状態は残す）
    for sd in SDIR.values():
        out = os.path.join(RES, "aggregated", v, sd); pt = os.path.join(out, "PARTIAL.txt")
        keep = [l for l in (open(pt).readlines() if os.path.exists(pt) else []) if l.split(" ")[1] not in sel] if os.path.exists(pt) else []
        new = keep + partial_notes.get(out, [])
        if new: os.makedirs(out, exist_ok=True); open(pt, "w").writelines(new)
        elif os.path.exists(pt): os.remove(pt)
for pf, text in prov_notes.items():
    if text: open(pf, "w").write(text)
    elif os.path.exists(pf): os.remove(pf)
mp = os.path.join(RES, "aggregated", "manifest.csv"); os.makedirs(os.path.dirname(mp), exist_ok=True); MH = ["variant", "target", "model", "sdir", "kind", "nwin", "missing", "nconfirm", "nprovisional", "finalization", "source", "file", "md5", "header"]
keep = []
if A.targets and os.path.exists(mp):  # 部分実行: 対象外の既存行は残し、今回の対象の行だけ置き換える
    old_rows = list(csv.reader(open(mp)))
    if old_rows and old_rows[0] == MH: keep = [r for r in old_rows[1:] if r and r[1] not in sel]
with open(mp, "w", newline="") as f:
    w = csv.writer(f); w.writerow(MH); w.writerows(keep); w.writerows(manifest)
nprov = sum(int(r[8]) for r in manifest if r[4] == "volf"); print(f"written: {len(outputs)} files, manifest -> {mp}; provisional windows (volf rows, both variants): {nprov}" + (" -> NOT final" if nprov else ""))
