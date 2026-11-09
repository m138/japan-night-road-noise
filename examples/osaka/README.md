# Expected output of the Osaka test run

Produced by

```
bash get_data.sh data 27      # about 1 minute, about 140 MB
bash run_all.sh data out 27   # about 1.5 minutes on 4 cores
```

on 2026-11-09 (Python 3.13). Compare your `out/` with these files.

- table2_accuracy.csv, fig1_validation_scatter.png: validation and cross-validation always use the national monitoring data, so they equal Table 2 and Fig. 1 of the paper.
- table3_exposure_by_pref.csv, prefectures.csv: Osaka only. Shares at 55, 60 and 65 dB or higher of 23.7%, 9.2% and 1.6% and a mean Lnight of 47.6 dB, as in Table 3 of the paper.
- municipalities.csv, fig2_map_municipal_p60.png: the 72 municipalities of Osaka (the map shows only Osaka; geopandas warns that the other prefectures are empty, which is expected).
