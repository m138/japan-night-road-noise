#!/usr/bin/env bash
# Main analysis. Usage: bash run_all.sh DATA_DIR OUT_DIR [comma-separated prefecture codes]
set -euo pipefail
D=${1:?DATA_DIR}; O=${2:?OUT_DIR}; P=${3:-}
H=$(cd "$(dirname "$0")" && pwd)/src
mkdir -p "$O"
python3 "$H"/validate_national.py "$D" "$H"                    # match monitoring sites to census sections -> $D/national_validation_points.csv
python3 "$H"/cv.py "$D"/national_validation_points.csv "$O"     # cross-validation (four schemes)
python3 "$H"/cv2.py "$D"/national_validation_points.csv "$O"    # comparison of bias-correction schemes
python3 "$H"/noise_national.py "$D" "$H" "$D"/nat $P            # national receptors (carriageway receptors removed) -> $D/nat/sa_XX.parquet
python3 "$H"/make_figs_en.py "$D" "$O"                          # Figs 1-3, Tables 2-3 (English labels; make_figs.py gives Japanese labels)
python3 "$H"/export_areas.py "$D" "$O"                         # municipalities.csv, prefectures.csv, municipalities.gpkg
