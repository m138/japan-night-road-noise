# 日本全国の夜間道路交通騒音（Lnight）推計モデル

公開データだけから、日本の全住民の夜間（22〜6時）道路交通騒音暴露を推計し、自動車騒音の常時監視の実測値で検証するコードです。論文「Nationwide nighttime road traffic noise exposure in Japan estimated from open data」（Transportation Research Part D に投稿予定）の再現用です。英語版は README.md です。

> 作成：松下大輔（Daisuke Matsushita, ORCID 0000-0002-6765-9366）。公開前の草稿です。DOIと論文の引用情報はまだ決まっていません（[ ] の箇所）。

## できること
- 道路交通センサス（令和3年度）の約9.6万区間を音源として、ASJ RTN-Model 2018の式で夜間の音響パワーレベルを計算します。
- 令和2年国勢調査の小地域に、約10人に1点の受音点（高さ4 m、乱数シード20261106）を置いて、Lnightを計算します。車道上に落ちた受音点は除き、その小地域の人口は残りの受音点に割り振ります。
- 常時監視の夜間LAeq（2019〜2023年度、12,648件）と照合し、偏りを補正します。一般道は地域ブロック別、自動車専用道路は道路種別に補正します。
- 4つの方法で交差検証します（都道府県押し出し、地域ブロック押し出し、年度分割、地点単位10分割）。
- 小地域・市区町村・都道府県ごとに、45〜65 dB以上の人口割合と人口加重平均Lnightを集計します。
- 感度分析：計算範囲500 m、平均Lnightの下限30/35/45 dB、別の乱数シード（1、2、3）、ASJ RTN-Model 2023の係数、車道上の受音点を残す場合、人口密度を調整した市区町村順位の一致。
- 大阪府・兵庫県・和歌山県でOpenStreetMapの建物を使う詳細版。大阪府では生活道路を加えた版と、国の面的評価との比較。

## 使い方
まず大阪府だけで試せます（約3分）。結果は examples/osaka/ と比べてください。
```
conda env create -f environment.yml && conda activate jnrn   # または pip install -r requirements.txt（Python 3.13）
bash get_data.sh data 27          # 大阪府の入力データ（約140 MB）
bash run_all.sh data out 27       # 照合 → 交差検証 → 暴露計算 → 図表・地域別集計
```
全国で実行する場合：
```
bash get_data.sh data                    # 全国の入力データ（約1 GB）
bash run_all.sh data out                 # 主解析
bash run_sensitivity.sh data out         # 感度分析（補足表S5・S6）
bash run_detailed.sh data out kansai.osm.pbf   # 3府県の詳細版（補足表S2・S3、図S1）
```
- 実測照合と交差検証は常に全国の実測値を使います。大阪府だけの試行でも表2と図1は論文と同じになります。
- 出力：`data/nat/sa_XX.parquet`（小地域）、`out/municipalities.csv`・`out/prefectures.csv`（45〜65 dB以上の割合と平均Lnight）、`out/municipalities.gpkg`（同じ値の市区町村ポリゴン、GIS用）、表2・表3と図1〜3。
- 推計値だけが必要な場合は、モデルを動かさずにZenodo（[DOI]）からダウンロードしてください。
- run_detailed.sh には3府県を含むOpenStreetMapの抽出ファイルが必要です。大阪府の生活道路版には、get_data.sh で取得しない2つのファイル（kasyo27.csv：道路交通センサスの大阪府の箇所別基本表、osaka_souon_2019_2023.csv：常時監視の大阪府の記録）も必要で、ない場合はその部分を飛ばします。

## ファイル
| ファイル | 内容 |
|---|---|
| src/noise_lib.py | 夜間時間配分、音源（ASJ RTN-Model 2018。環境変数 ASJ_RTN=2023 で2023年版の係数）、線音源の積分 |
| src/validate_national.py | 常時監視地点とセンサス区間の対応づけと照合 |
| src/cv.py, src/cv2.py | 交差検証、補正方法の比較 |
| src/noise_national.py | 全国の受音点計算（車道上の受音点を除く） |
| src/make_figs_en.py, src/make_figs.py | 図1〜3、表2・表3（英語／日本語ラベル） |
| src/export_areas.py | 市区町村・都道府県の集計（CSV、GeoPackage） |
| src/sensitivity/ | 感度分析（sens_national.py、agg_sens.py、rank_stats.py、dens_check.py、dens_check_3pref.py、profile_split.py） |
| src/detailed/ | 3府県の詳細版と大阪府の生活道路版（osm_bld_pref.py、detailed.py、cmp_3pref.py、make_figS1.py、osaka_validate_spatial.py、osaka_local_eval.py） |

## モデルの主な仮定と限界
- 半自由空間の有限長線音源で、建物・地表面・防音壁による減衰は入れていません。偏りは実測との照合で補正しています。
- 音源は密粒舗装の式（2車種分類）です。速度は指定最高速度とし、自動車専用道路は定常、一般道は非定常の式を使います。
- 夜間の時間配分は大阪府の24時間観測区間の値を全国に使っています。
- 全国版はセンサス対象外の生活道路を含みません。センサス道路から離れた地域では暴露の下限とみてください。
- 補正後の精度（都道府県押し出し）は、一般道でRMSE 4.48 dB、自動車専用道路で5.89 dBです。個々の住宅の予測より、地域の比較に向きます。

## ASJ RTN-Model 2018について
音源の係数は次の論文の表2.3によります。このコードは日本音響学会の公式実装ではありません。
Sakamoto S. Road traffic noise prediction model "ASJ RTN-Model 2018": Report of the Research Committee on Road Traffic Noise. Acoust Sci Tech. 2020;41(3):529–589. doi:10.1250/ast.41.529

## 入力データの出典と利用条件
DATA_SOURCES.md を見てください。入力データ自体はこのリポジトリに含めず、get_data.sh で公開元から取得します。

## ライセンス
- コード：MIT License（LICENSE）
- 推計結果のデータ（Zenodoで公開予定）：CC BY 4.0。OpenStreetMapから派生した結果はODbLの求める出典表示をつけて公開します。

## 引用
[論文の書誌情報とZenodoのDOI]
