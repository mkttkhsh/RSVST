#!/usr/bin/env python3
"""確認パスの併合出力と状態ファイル（tools/azconfirm_check.py から呼ぶ。tools/test_*.py も同じ書式で fixture を作る）。規則 r3_snippets/azconfirm_rules_20261002.md §7-§8。
併合出力（azdual_results/confirm/<target>/、dual 版の窓別出力と同じ 7 種 ＋ 併合事後要約 1 種）:
  h1_y_<m>_<t>.csv   チェーン順に連結した予測ドロー（[h1 y1 (x1)] fixed | orig）。行数 = チェーン数 × 保存ドロー数
  pstate_<m>_<t>.csv 同じ順の状態ドロー
  volf / var / es（_orig も）: 併合ドローから Ox の定義で再計算（volf = mean(h1) ~ quantilec(h1) の線形補間中央値、var = quantilec(y1, <.005,.01,.05,.1>)、es = VaR 未満のドローの平均）
  posstat_pooled_<m>_<t>.csv: パラメータと h_T の併合事後要約（mean sd 2.5% 97.5% rhat ess_bulk ess_tail mcse_mean、ヘッダ付き）。
  チェーン個別の posstat / posh / posstath / aprob は複製しない（併合値ではないため。各チェーンの値はパッケージの out_confirm/<target>/<outdir>/ にある）
手順: 仮出力ディレクトリ _tmp_<m>_<t>/ に書く → 読み戻して検証（行数、有限性、volf の mean(h) と h1_y の一致）→ 正式名に移す → 併合ファイルの md5 を返す。
      途中で失敗したら仮出力を消し、正式名のファイルは作らない。
状態ファイル azdual_results/confirm/confirm_status.csv: 窓ごとに 1 行（同じ窓は上書き）。status ∈ {accepted, rejected, incomplete, invalid}、pooled ∈ {yes, no}。
  input_md5 = 入力チェーンのファイル（pstate, h1_y, posstath）の md5 を "c<k>:<md5>,<md5>,<md5>;..." で結び付け、pooled_md5 = 併合ファイルの "<name>:<md5>;..."。
  azdual_aggregate.py --confirm-dir は status == accepted かつ pooled == yes の窓だけを取り込み、併合ファイルの md5 を pooled_md5 と照合する。
"""
import csv, hashlib, math, os, shutil, statistics, time
from mcmc_diag import quantile7, summary

STATUS_FIELDS = ["target", "model", "t", "stage", "nsave", "chains", "chain_ids", "ndraw_per_chain", "ndraw_total", "status", "verdict", "okA", "okB_std", "okB_ref", "pooled", "input_md5", "pooled_md5", "checked_at", "code"]
PROBS = (0.005, 0.01, 0.05, 0.1)
POOL_KINDS = ("h1_y", "pstate", "volf", "var", "es", "volf_orig", "var_orig", "es_orig", "posstat_pooled")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def ox_var_es(y):
    s = sorted(y); var = [quantile7(s, p, True) for p in PROBS]
    es = []
    for v in var:
        tail = [x for x in s if x < v]; es.append(sum(tail) / len(tail) if tail else float("nan"))
    return var, es

def ox_volf(h): return [sum(h) / len(h), quantile7(h, 0.5)]

def savemat(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow([""] + [f"Var{i+1}" for i in range(len(rows[0]))])
        for i, r in enumerate(rows): w.writerow([i + 1] + [repr(float(x)) for x in r])

def _readmat(p):
    with open(p) as fh:
        r = csv.reader(fh); next(r); return [[float(x) if x not in ("inf", "-inf", "nan") else float(x) for x in row[1:]] for row in r if row]

def pool_window(od, model, t, PS, H, ncol, pnames, pidx, hT_idx):
    """PS: チェーン別 pstate 行列、H: チェーン別 h1_y 行列。仮出力 → 検証 → 正式名。戻り値: {ファイル名: md5}。失敗時は例外（仮出力は消す、正式名は作らない）。"""
    stem = f"{model}_{t}"; tmp = os.path.join(od, f"_tmp_{stem}"); shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    try:
        hp = [r for hy in H for r in hy]; psp = [r for ps in PS for r in ps]; ntot = len(hp)
        savemat(os.path.join(tmp, f"h1_y_{stem}.csv"), hp); savemat(os.path.join(tmp, f"pstate_{stem}.csv"), psp)
        expect = {}
        for off, sfx in ((0, ""), (ncol // 2, "_orig")):
            h1 = [r[off] for r in hp]; y1 = [r[off + 1] for r in hp]; var, es = ox_var_es(y1); vf = ox_volf(h1)
            savemat(os.path.join(tmp, f"volf{sfx}_{stem}.csv"), [vf]); savemat(os.path.join(tmp, f"var{sfx}_{stem}.csv"), [var]); savemat(os.path.join(tmp, f"es{sfx}_{stem}.csv"), [es]); expect[sfx] = (vf, var, es)
        with open(os.path.join(tmp, f"posstat_pooled_{stem}.csv"), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["quantity", "mean", "sd", "q025", "q975", "rhat", "ess_bulk", "ess_tail", "mcse_mean"])
            for nm in list(pnames) + ["h_T"]:
                idx = hT_idx if nm == "h_T" else pidx[nm]; s = summary([[r[idx] for r in ps] for ps in PS])
                w.writerow([nm] + [repr(float(s[k])) for k in ("mean", "sd", "q025", "q975", "rhat", "ess_bulk", "ess_tail", "mcse_mean")])
        # ---- 読み戻して検証 ----
        hb = _readmat(os.path.join(tmp, f"h1_y_{stem}.csv")); pb = _readmat(os.path.join(tmp, f"pstate_{stem}.csv"))
        if len(hb) != ntot or len(pb) != ntot: raise RuntimeError(f"pooled rows {len(hb)}/{len(pb)} != {ntot}")
        if len(hb[0]) != ncol: raise RuntimeError("pooled h1_y column count")
        if not all(math.isfinite(x) for r in hb for x in r): raise RuntimeError("non-finite in pooled h1_y")
        for off, sfx in ((0, ""), (ncol // 2, "_orig")):
            vf = _readmat(os.path.join(tmp, f"volf{sfx}_{stem}.csv"))[0]; va = _readmat(os.path.join(tmp, f"var{sfx}_{stem}.csv"))[0]; es = _readmat(os.path.join(tmp, f"es{sfx}_{stem}.csv"))[0]
            h1 = [r[off] for r in hb]
            if abs(vf[0] - sum(h1) / ntot) > 1e-8 or abs(vf[1] - statistics.median(h1)) > 1e-6: raise RuntimeError(f"pooled volf{sfx} does not match pooled h1_y")
            if len(va) != 4 or len(es) != 4 or not all(math.isfinite(x) for x in va + es): raise RuntimeError(f"pooled var/es{sfx} not 1x4 finite")
            if any(abs(a - b) > 1e-12 for a, b in zip(va, expect[sfx][1])) or any(abs(a - b) > 1e-12 for a, b in zip(es, expect[sfx][2])): raise RuntimeError("pooled var/es read-back mismatch")
        # ---- 正式名へ ----
        out = {}
        for k in POOL_KINDS:
            fn = f"{k}_{stem}.csv"; shutil.move(os.path.join(tmp, fn), os.path.join(od, fn)); out[fn] = md5(os.path.join(od, fn))
        shutil.rmtree(tmp, ignore_errors=True); return out
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True); raise

def remove_pooled(od, model, t):
    """その窓の併合出力（仮出力も）を消す。消したファイル名を返す。"""
    stem = f"{model}_{t}"; gone = []
    for k in POOL_KINDS:
        p = os.path.join(od, f"{k}_{stem}.csv")
        if os.path.exists(p): os.remove(p); gone.append(os.path.basename(p))
    tmp = os.path.join(od, f"_tmp_{stem}")
    if os.path.isdir(tmp): shutil.rmtree(tmp); gone.append(f"_tmp_{stem}/")
    return gone

def fmt_md5s(d): return ";".join(f"{k}:{v}" for k, v in sorted(d.items()))
def parse_md5s(s): return dict(kv.split(":", 1) for kv in s.split(";") if kv)

def read_status(path):
    if not os.path.exists(path): return {}
    out = {}
    for r in csv.DictReader(open(path)): out[(r["target"], r["model"], int(r["t"]))] = r
    return out

def write_status(path, rows_by_key):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_FIELDS); w.writeheader()
        for k in sorted(rows_by_key, key=lambda k: (k[0], k[1], k[2])): w.writerow({f: rows_by_key[k].get(f, "") for f in STATUS_FIELDS})

def now(): return time.strftime("%Y-%m-%d %H:%M:%S")
