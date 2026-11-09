# Nationwide nighttime road traffic noise (Lnight) in Japan from open data

Code for the paper "Nationwide nighttime road traffic noise exposure in Japan estimated from open data" (Daisuke Matsushita, ORCID 0000-0002-6765-9366; submitted to *Transportation Research Part D*). It estimates the nighttime (22:00–06:00) road traffic noise exposure of every resident of Japan from public data only, and checks the estimates against the national continuous monitoring of road traffic noise. A Japanese version of this file is in README_ja.md.

> Status: pre-publication. The citation and the Zenodo DOI of the estimates will be added when the paper is published ([ ] below).

## What the code does
- Uses about 96,000 sections of the 2021 Road Traffic Census as sources and calculates their nighttime sound power with ASJ RTN-Model 2018.
- Places receptors at a height of 4 m in the 2020 census small areas, at about one receptor per 10 residents (random seed 20261106). Receptors that fall on a carriageway are removed, and the population of the small area is divided among the remaining receptors.
- Matches the estimates with the measured nighttime LAeq of the continuous monitoring (fiscal years 2019–2023; 12,648 records at 8,761 sites) and corrects the bias by regional block (ordinary roads) and by road type (motorways).
- Cross-validates the correction in four ways: leave-one-prefecture-out, leave-one-block-out, split by fiscal year and 10-fold by monitoring site.
- Summarizes, for small areas, municipalities and prefectures, the shares of residents at 45–65 dB or higher and the population-weighted mean Lnight.
- Sensitivity analyses: 500 m calculation range, floors of 30/35/45 dB for the mean Lnight, other random seeds (1, 2, 3), the ASJ RTN-Model 2023 coefficient, receptors on the carriageway kept, and the agreement of municipal rankings after adjustment for population density.
- Detailed models with OpenStreetMap buildings in Osaka, Hyogo and Wakayama, and, for Osaka, a version with local roads and a comparison with the official area-wide assessment.

## Quick start (Osaka only, about 3 minutes)
```
conda env create -f environment.yml && conda activate jnrn   # or: pip install -r requirements.txt (Python 3.13)
bash get_data.sh data 27        # download the inputs for Osaka (about 140 MB)
bash run_all.sh data out 27     # validation -> cross-validation -> exposure -> figures, tables and area summaries
```
Compare `out/` with `examples/osaka/`. Validation and cross-validation always use the national monitoring data, so Table 2 and Fig. 1 are those of the paper even in this test run.

## Full run
```
bash get_data.sh data                           # all inputs (about 1 GB)
bash run_all.sh data out                        # main analysis for the whole country
bash run_sensitivity.sh data out                # sensitivity analyses (Supplementary Tables S5 and S6)
bash run_detailed.sh data out kansai.osm.pbf    # detailed models in three prefectures (Supplementary Tables S2, S3, Fig. S1)
```
- The national model takes from a few seconds to a few tens of seconds per prefecture (4 cores).
- run_detailed.sh needs an OpenStreetMap extract that covers the three prefectures. The Osaka local-road part also needs two files that get_data.sh does not download: `kasyo27.csv` (the section table of the 2021 Road Traffic Census for Osaka, from the census web site of the Ministry of Land, Infrastructure, Transport and Tourism, converted to UTF-8) and `osaka_souon_2019_2023.csv` (the Osaka records of the continuous monitoring, fiscal years 2019–2023, taken from `data/tn`). Without them this part is skipped.

## Outputs of run_all.sh
| File | Content |
|---|---|
| data/nat/sa_XX.parquet | Small areas: population weight `w`, population at 45–65 dB or higher `w45`…`w65`, `wL` (population × Lnight with a 40 dB floor) |
| out/municipalities.csv, out/prefectures.csv | Shares of residents at 45–65 dB or higher (%) and mean Lnight (dB) |
| out/municipalities.gpkg | The same municipal values with polygons (JGD2011), for GIS |
| out/table2_accuracy.csv, out/table3_exposure_by_pref.csv | Tables 2 and 3 of the paper |
| out/fig1–3 (png, eps) | Figures 1–3 of the paper |

If you only need the estimates, download them from Zenodo ([DOI]) instead of running the model.

## Files
| File | Content | Paper |
|---|---|---|
| src/noise_lib.py | Nighttime profile, source model (ASJ RTN-Model 2018; `ASJ_RTN=2023` switches to the 2023 coefficient), line-source integration | 3.2–3.4 |
| src/validate_national.py | Matching of monitoring sites with census sections | 3.6, 4.1 |
| src/cv.py, src/cv2.py | Cross-validation and comparison of correction schemes | 3.7, 4.2 |
| src/noise_national.py | National receptors and exposure | 3.5, 4.3 |
| src/make_figs_en.py (make_figs.py) | Figs 1–3, Tables 2 and 3 (English / Japanese labels) | 4.1–4.3 |
| src/export_areas.py | Municipal and prefectural summaries (CSV, GeoPackage) | 4.3 |
| src/sensitivity/ | sens_national.py, agg_sens.py, rank_stats.py, dens_check.py, dens_check_3pref.py, profile_split.py | 3.3, 3.10, 4.7, Tables S5, S6 |
| src/detailed/ | osm_bld_pref.py, detailed.py, cmp_3pref.py, make_figS1.py, osaka_validate_spatial.py, osaka_local_eval.py | 3.8, 4.4, 4.5, Text S1, Tables S2, S3, Fig. S1 |

## Main assumptions and limitations
- Finite line sources over a hard half-space; no attenuation by buildings, ground or barriers in the national model. The bias is corrected with the monitoring data.
- Dense-graded asphalt and two vehicle classes. Speed is the posted limit; motorways use the steady-flow equation, ordinary roads the non-steady-flow equation.
- The nighttime share of traffic comes from the 24-hour counts in Osaka Prefecture and is used for the whole country.
- Local roads not covered by the census are not included in the national model, so the national shares are a lower bound for areas far from census roads.
- After correction, the leave-one-prefecture-out RMSE is 4.48 dB for ordinary roads and 5.89 dB for motorways. The estimates are meant for comparing areas, not for predicting the level at a single dwelling.

## ASJ RTN-Model
The source coefficients follow Table 2.3 of Sakamoto (2020). This code is not an official implementation of the Acoustical Society of Japan.
- Sakamoto, S., 2020. Road traffic noise prediction model "ASJ RTN-Model 2018": report of the Research Committee on Road Traffic Noise. Acoust. Sci. Technol. 41, 529–589. https://doi.org/10.1250/ast.41.529
- Okada, Y., Anai, K., Itiki, T., Yamauchi, K., Yasuda, Y., Yokota, T., 2025. Road traffic noise prediction model "ASJ RTN-Model 2023": report of the Research Committee on Road Traffic Noise. Acoust. Sci. Technol. 46, 365–429. https://doi.org/10.1250/ast.e24.125

## Input data
See DATA_SOURCES.md. The inputs are not included in this repository; get_data.sh downloads them from their sources.

## License
- Code: MIT License (LICENSE)
- Estimates (to be deposited in Zenodo): CC BY 4.0. Results derived from OpenStreetMap are shared under the attribution that the ODbL requires.

## Citation
[Citation of the paper and the Zenodo DOI]
