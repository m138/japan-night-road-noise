#!/usr/bin/env bash
# Sensitivity analyses of the national exposure (Section 3.10, Supplementary Tables S5 and S6). Run run_all.sh first.
# Usage: bash run_sensitivity.sh DATA_DIR OUT_DIR [comma-separated prefecture codes]
# Needs about 4 x the time and disk of the main run (receptor-level results are kept).
set -euo pipefail
D=${1:?DATA_DIR}; O=${2:?OUT_DIR}; P=${3:-}
H=$(cd "$(dirname "$0")" && pwd)/src; S="$O"/sens; mkdir -p "$S" "$O"/sensitivity
# validation with the ASJ RTN-Model 2023 coefficient (for its own bias correction)
cp "$D"/national_validation_points.csv "$D"/_vp_main.csv
ASJ_RTN=2023 python3 "$H"/validate_national.py "$D" "$H"
mv "$D"/national_validation_points.csv "$D"/national_validation_points_asj2023.csv; mv "$D"/_vp_main.csv "$D"/national_validation_points.csv
V23="$D"/national_validation_points_asj2023.csv
python3 "$H"/sensitivity/sens_national.py "$D" "$H" "$V23" 20261106 500 "$S"/s0 $P   # main seed; keeps 200 m and 500 m energies
for i in 1 2 3; do python3 "$H"/sensitivity/sens_national.py "$D" "$H" "$V23" $i 200 "$S"/s$i $P; done   # other random seeds 1, 2, 3
python3 "$H"/sensitivity/agg_sens.py "$S" "$O"/sensitivity v4         # Table S5 (shares, mean Lnight, rho, largest rank change)
python3 "$H"/sensitivity/rank_stats.py "$S" "$O"/sensitivity/rank_stats.csv   # Table S5 (median rank change, top-100 kept)
python3 "$H"/sensitivity/dens_check.py "$S" "$D"/estat_all/shp "$O"/sensitivity   # Table S6 (density-adjusted agreement)
python3 "$H"/sensitivity/profile_split.py "$D"/census27.csv > "$O"/sensitivity/profile_split.txt   # out-of-sample check of the nighttime profile
