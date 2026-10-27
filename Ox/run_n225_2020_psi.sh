#!/bin/bash
# 2026 revision: launch all 5 N225 COVID-period (2012-2020) psi-free RSV
# models (AZ family: N/SN/ST/T + GH-ST) as background processes, one per
# model. Needs ~5 free cores; each runs its own full 242-window rolling
# loop serially. FS family (SN/ST) is a separate, heavier job — see
# server_upload_n225_2020_fs_{st,sn}_psi/.
#
# Each job is fully detached (setsid + disown) so it survives logout and is
# immune to Ctrl+Z / kill %N on this script -- do NOT remove that.
set -eu
cd "$(dirname "$0")"

MODELS=(
  rsv_az_n_psi_n225_prediction_2020
  rsv_az_sn_psi_n225_prediction_2020
  rsv_az_st_psi_n225_prediction_2020
  rsv_az_t_psi_n225_prediction_2020
  rsv_gh_st_psi_n225_prediction_2020
)

mkdir -p N225_prediction_2020 logs
for m in "${MODELS[@]}"; do
  echo "launching ${m}.ox -> logs/${m}.log"
  if command -v setsid >/dev/null 2>&1; then
    setsid nohup oxl "${m}.ox" < /dev/null > "logs/${m}.log" 2>&1 &
  else
    nohup oxl "${m}.ox" < /dev/null > "logs/${m}.log" 2>&1 &
  fi
  disown
done
echo "launched ${#MODELS[@]} model(s); monitor with: tail -f logs/*.log"
