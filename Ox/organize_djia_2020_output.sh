#!/bin/bash
# IJF R2 (2026): merge_prediction_results_{az,gh}_2020.ox read per-window
# output from DJIA_prediction_2020/files/, but the rolling-forecast scripts
# write directly into DJIA_prediction_2020/ (same convention already used by
# the main-sample and N225 pipelines). Run this once estimation is finished,
# before running the merge_prediction_results_*_2020.ox scripts.
set -eu
cd "$(dirname "$0")"
mkdir -p DJIA_prediction_2020/files
mv DJIA_prediction_2020/*.csv DJIA_prediction_2020/files/
echo "moved $(ls DJIA_prediction_2020/files | wc -l) files into DJIA_prediction_2020/files/"
