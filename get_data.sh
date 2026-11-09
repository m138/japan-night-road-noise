#!/usr/bin/env bash
# 入力データをすべて公開元から取得する。使い方: bash get_data.sh DATA_DIR [都道府県コード ...]
# 都道府県コードを省略すると全国（01〜47）を取得する。
set -euo pipefail
D=${1:?DATA_DIR}; shift || true
PREFS=${*:-$(seq -w 1 47)}
mkdir -p "$D"/geo "$D"/estat_all/zip "$D"/estat_all/shp "$D"/tn
# 1) 道路交通センサス R3 の区間線形（MIT。元データは国土交通省、公共データ利用規約1.0）
curl -fL -o "$D"/geo/traffic_census_2021_converted.parquet \
  https://github.com/shiwaku/mlit-road-traffic-census-converter/releases/download/data-v1/traffic_census_2021_converted.parquet
# 2) 時間帯別交通量（大阪府、24時間観測区間）。夜間の時間配分に使う
curl -fL https://www.mlit.go.jp/road/census/r3/data/csv/zkntrf27.csv | iconv -f cp932 -t utf-8 > "$D"/census27.csv
# 3) 令和2年国勢調査 小地域境界（人口つき、JGD2011、shape）。e-Stat
for p in $PREFS; do
  curl -fL -o "$D"/estat_all/zip/$p.zip \
    "https://www.e-stat.go.jp/gis/statmap-search/data?dlserveyId=A002005212020&code=$p&coordSys=1&format=shape&downloadType=5&datum=2011"
  unzip -o -q "$D"/estat_all/zip/$p.zip -d "$D"/estat_all/shp
done
# 4) 自動車騒音の常時監視（2019〜2023年度）。環境GIS＋（国立環境研究所）
for y in 2019 2020 2021 2022 2023; do
  curl -fL --retry 5 --retry-delay 3 --retry-all-errors -o "$D"/tn/tn$y.zip https://tenbou.nies.go.jp/gisplus/tn/tn$y.zip
  unzip -o -q "$D"/tn/tn$y.zip -d "$D"/tn/tn$y
done
echo "done: $D"
