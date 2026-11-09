# 入力データの出典と利用条件

| データ | 公開元 | 利用条件（2026-11-08時点で確認した範囲） | 公表時の表示 |
|---|---|---|---|
| 道路交通センサス R3 区間線形（traffic_census_2021_converted.parquet） | GitHub: shiwaku/mlit-road-traffic-census-converter（release data-v1） | 変換コードとデータはMIT。元データは国土交通省で、公共データ利用規約（第1.0版）に基づく | 出典：国土交通省「令和3年度 全国道路・街路交通情勢調査」を加工 |
| 時間帯別交通量（大阪府 zkntrf27.csv） | 国土交通省 道路局 | 政府標準利用規約（第2.0版、CC BY 4.0互換） | 出典：国土交通省「令和3年度 全国道路・街路交通情勢調査」を加工 |
| 令和2年国勢調査 小地域境界 | e-Stat（総務省統計局） | 政府標準利用規約（第2.0版、CC BY 4.0互換） | 出典：政府統計の総合窓口（e-Stat）「令和2年国勢調査 小地域（町丁・字等別）境界データ」を加工 |
| 自動車騒音の常時監視（souonYYYY.csv） | 一次出典は環境省（水・大気環境局 自動車環境対策課）。環境GIS＋（国立環境研究所）から取得 | 「自由に利用可。公表時はクレジット表示」（環境GIS＋のダウンロードページの利用規約） | 出典：環境省「自動車騒音の常時監視結果」（国立環境研究所 環境GIS＋より取得）を加工 |

- 3府県（大阪府・兵庫県・和歌山県）の詳細版と大阪府の生活道路版はOpenStreetMap（© OpenStreetMap contributors、ODbL 1.0）を使います。建物・道路のデータそのものは配布せず、そこから得た結果はODbLの求める出典表示をつけて公開します。全国版には使いません。
- Detailed models (Osaka, Hyogo, Wakayama) use OpenStreetMap data (© OpenStreetMap contributors, ODbL 1.0). The building and road data are not redistributed; results derived from them are shared with the attribution that the ODbL requires. The national model does not use them.
- 利用条件は2026-11-08に確認した。いずれも出典を表示すれば加工・再配布できる。ただし、このリポジトリは入力データを含めず、get_data.sh で公開元から取得する。
