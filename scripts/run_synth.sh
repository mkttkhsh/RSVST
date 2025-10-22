#!/usr/bin/env bash
set -euo pipefail
oxl ox/run_estimation.ox --config=config/config.ini.example --mode=synthetic
echo "Done. See results/ for outputs."
