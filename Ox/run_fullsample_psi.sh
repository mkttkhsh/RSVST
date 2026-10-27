#!/bin/bash
# 2026 revision: launch the 14 single full-sample in-sample psi-free
# estimations (main-sample period, NOT rolling) -- one MCMC run per model,
# not 590 rolling windows. Reports the posterior distribution / 95% CI of
# psi for the sensitivity analysis (Comment 2), to assess whether psi≈1 is
# supported by the data.
#
# psi=1-fixed RSV variants and SV were both considered and dropped
# 2026-07-16: this package's whole purpose is psi's posterior/CI, and
# neither has anything to report there -- SV has no psi at all (no
# realized-measure equation), and a psi=1-fixed model has psi pinned by
# construction, not estimated, so there's no CI to speak of. Both
# exclusions mirror Watanabe-sensei's own reasoning for excluding EGARCH
# (1)(2) -- and then REGARCH-n/REGARCH-t with psi=1 fixed, (3)(5) -- from
# his parallel REGARCH in-sample estimation request (only (4)(6), the
# psi-free REGARCH variants, remain there too).
#
# Model set: RSV-AZ {N,SN,ST,T} + RSV-GH-ST (10) + RSV-FS {SN,ST} (4),
# psi-free, for each of N225 and DJIA = 14.
#
# MCMC settings (2026-07-16, revised up from the rolling-derived
# 15,000/5,000): AZ/GH use 40,000 draws / 10,000 burn-in (total 50,000),
# matching this paper's own RSV simulation-study convention (see
# rsvst202601suppl.tex) -- since this is a single one-off run (not 590
# rolling windows), the extra cost is easily affordable and buys a tighter
# posterior CI, which matters since the whole point is judging whether the
# psi CI excludes 1. FS is untouched (already 90,000/60,000, heavier than
# even the revised AZ/GH setting).
#
# ⚠️ FS MODELS ARE STILL MUCH HEAVIER: 90,000/60,000 vs 40,000/10,000 for
# AZ/GH -- roughly ~2x the cost on top of already being a bigger model,
# likely several hours per model (vs ~2-3h for the revised AZ/GH runs).
# Consider launching the 10 AZ/GH models first and the 4 FS models
# separately/later if the server is busy.
#
# Needs ../Data/{N225,DJIA}_rcc_rv_200906_201909.csv (main-sample RV5 input
# for AZ/GH) and ../Data/{Nikkei,DJIA}_rcc_rv_2009_2019.csv (main-sample
# input for FS -- a different file, FS's own recalculated series, verified
# 2026-07-16 to contain the same underlying return/RV data as the AZ/GH
# file, just different formatting/precision -- see PROGRESS.md) -- already
# present.
#
# Usage: run from inside Ox/ on the server.
#   ./run_fullsample_psi.sh            # launch all 14 models
#   ./run_fullsample_psi.sh MODEL...   # launch only the named model(s), e.g.
#                                       #   ./run_fullsample_psi.sh rsv_az_st_psi_n225 rsv_fs_st_psi_djia
#
# Each job is fully detached (setsid + disown) so it survives logout and is
# immune to Ctrl+Z / kill %N on this script -- do NOT remove that.

set -eu
cd "$(dirname "$0")"

ALL_MODELS=(
  rsv_az_n_psi_n225 rsv_az_sn_psi_n225 rsv_az_st_psi_n225 rsv_az_t_psi_n225 rsv_gh_st_psi_n225
  rsv_az_n_psi_djia rsv_az_sn_psi_djia rsv_az_st_psi_djia rsv_az_t_psi_djia rsv_gh_st_psi_djia
  rsv_fs_st_psi_n225 rsv_fs_sn_psi_n225
  rsv_fs_st_psi_djia rsv_fs_sn_psi_djia
)

MODELS=("$@")
if [ ${#MODELS[@]} -eq 0 ]; then
  MODELS=("${ALL_MODELS[@]}")
fi

mkdir -p N225_fullsample_psi DJIA_fullsample_psi logs
mkdir -p output_emp_N225_fs_st_psi output_emp_N225_fs_sn_psi
mkdir -p output_emp_DJIA_fs_st_psi output_emp_DJIA_fs_sn_psi

for m in "${MODELS[@]}"; do
  script="${m}_fullsample.ox"
  if [ ! -f "$script" ]; then
    echo "skip: $script not found" >&2
    continue
  fi
  echo "launching $script -> logs/${m}_fullsample.log"
  if command -v setsid >/dev/null 2>&1; then
    setsid nohup oxl "$script" < /dev/null > "logs/${m}_fullsample.log" 2>&1 &
  else
    nohup oxl "$script" < /dev/null > "logs/${m}_fullsample.log" 2>&1 &
  fi
  disown
done

echo "launched ${#MODELS[@]} model(s); monitor with: tail -f logs/*_fullsample.log"
