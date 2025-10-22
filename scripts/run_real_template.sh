#!/usr/bin/env bash
set -euo pipefail
# Copy config.ini.example to config.ini and edit as needed.
# Then set your real data directory (licensed; not included) and run:
oxl ox/run_estimation.ox --config=config/config.ini --mode=real --data_dir=/ABS/PATH/TO/REALDATA
