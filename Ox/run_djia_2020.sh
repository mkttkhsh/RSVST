#!/bin/bash
# IJF R2 (2026): launch all 15 DJIA COVID-period (2012-2020) rolling MCMC
# prediction models as background processes, one per model, mirroring the
# copy_oxfile.sh/runpred.sh parallelism already used for the FS models.
#
# Model set mirrors N225's COVID-period coverage: RSV with psi fixed at 1
# (baseline) x5, RSV with psi free (sensitivity) x5, SV (no psi concept --
# it doesn't use the log-RV measurement equation) x5. Added 2026-07-16
# (rsv_az_{n,sn,st,t}_psi) to match N225, which already had all 5 psi-free
# RSV distributions; previously only rsv_gh_st_psi existed for DJIA.
#
# Usage: run from inside Ox/ on the server, after ../Data/ is populated with
#   DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_tick.csv
#   ./run_djia_2020.sh            # launch all 15 models
#   ./run_djia_2020.sh MODEL...   # launch only the named model(s), e.g.
#                                  #   ./run_djia_2020.sh rsv_az_st sv_az_st
#
# Each model runs its own full 253-window rolling loop serially (like the
# main-sample scripts); parallelism is across the 15 models, not within one.
# Needs ~15 free cores; reduce concurrency by passing a subset of models.
#
# Each job is fully detached (setsid + disown) so it survives logout and is
# immune to Ctrl+Z / kill %N on this script -- do NOT remove that.

set -eu
cd "$(dirname "$0")"

ALL_MODELS=(
  rsv_az_n rsv_az_sn rsv_az_st rsv_az_t
  rsv_gh_st
  rsv_az_n_psi rsv_az_sn_psi rsv_az_st_psi rsv_az_t_psi rsv_gh_st_psi
  sv_az_n sv_az_sn sv_az_st sv_az_t
  sv_gh_st
)

MODELS=("$@")
if [ ${#MODELS[@]} -eq 0 ]; then
  MODELS=("${ALL_MODELS[@]}")
fi

mkdir -p DJIA_prediction_2020 logs

for m in "${MODELS[@]}"; do
  script="${m}_djia_prediction_2020.ox"
  if [ ! -f "$script" ]; then
    echo "skip: $script not found" >&2
    continue
  fi
  echo "launching $script -> logs/${m}_djia_2020.log"
  if command -v setsid >/dev/null 2>&1; then
    setsid nohup oxl "$script" < /dev/null > "logs/${m}_djia_2020.log" 2>&1 &
  else
    nohup oxl "$script" < /dev/null > "logs/${m}_djia_2020.log" 2>&1 &
  fi
  disown
done

echo "launched ${#MODELS[@]} model(s); monitor with: tail -f logs/*.log"
