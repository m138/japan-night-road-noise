#!/usr/bin/env bash
# Detailed models with buildings in Osaka (27), Hyogo (28) and Wakayama (30) (Sections 3.8, 4.4; Supplementary Tables S2, S2b, S6, Fig. S1).
# Run run_all.sh first. Usage: bash run_detailed.sh DATA_DIR OUT_DIR OSM_PBF
# OSM_PBF: an OpenStreetMap extract covering the three prefectures (for example the Geofabrik Kansai extract). Data (c) OpenStreetMap contributors, ODbL 1.0.
set -euo pipefail
D=${1:?DATA_DIR}; O=${2:?OUT_DIR}; PBF=${3:?OSM_PBF}
H=$(cd "$(dirname "$0")" && pwd)/src; X="$O"/detailed; mkdir -p "$X" "$O"/three_prefectures
for pc in 27 28 30; do
  python3 "$H"/detailed/osm_bld_pref.py "$PBF" "$D"/estat_all/shp $pc "$X"/bld_$pc.parquet
  python3 "$H"/detailed/detailed.py "$D" "$H" $pc "$X"/bld_$pc.parquet "$X"
done
python3 "$H"/detailed/cmp_3pref.py "$D"/nat "$X" "$O"/three_prefectures
python3 "$H"/detailed/make_figS1.py "$O"/three_prefectures/municipal_3pref_v4.csv "$O"
python3 "$H"/sensitivity/dens_check_3pref.py "$D"/estat_all/shp "$O"/three_prefectures/municipal_3pref_v4.csv "$O"/three_prefectures/density_check_3pref.csv
# Osaka with local roads and the comparison with the official area-wide assessment (Supplementary Text S1, Table S3).
# Needs two extra files in DATA_DIR (see README): kasyo27.csv and osaka_souon_2019_2023.csv.
if [ -f "$D"/kasyo27.csv ] && [ -f "$D"/osaka_souon_2019_2023.csv ]; then
  python3 "$H"/detailed/osaka_validate_spatial.py "$D"
  python3 "$H"/detailed/osaka_local_eval.py "$D" "$H" "$X"/bld_27.parquet "$X" "$PBF" "$D"/osaka_validation_points_spatial.csv "$O"/osaka_local
else
  echo "skip Osaka local roads: kasyo27.csv or osaka_souon_2019_2023.csv not found in $D"
fi
