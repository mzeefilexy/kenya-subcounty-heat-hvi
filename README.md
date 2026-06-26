# Kenya heat-exposure screening reproducibility package

This release supports the manuscript *Where severe heat and priority-population burden coincide in Kenya: a biometeorological subcounty screening analysis, 2015-2025*.

The archived Zenodo record for this study is: https://zenodo.org/records/19202174

## What is included

- Final revised manuscript and STROBE checklist.
- Reproducible analysis scripts, without personal paths or credentials.
- Publication figures, main tables, supplementary results, reference-reconciliation register, software requirements, and Zenodo metadata.

## What is deliberately not included

No Copernicus/ERA5, WorldPop, Relative Wealth Index, or administrative-boundary source files are redistributed. These are third-party inputs and must be obtained under their original terms. No API keys, browser data, or personal files are included.

## Data inputs required

Set `HVI_DATA_DIR` to a folder containing `kenya_adm2_boundaries.geojson`, `kenya_subcounty_population_and_events_2015-2030.csv`, and `UTCI_files/UTCI_daily_stats_*_degC.nc`. Download ERA5 hourly 2-m temperature and dew-point temperature through the Copernicus Climate Data Store for the Heat Index comparator. The scripts expect the ERA5 files under `data/era5_sensitivity/`.

## Reproduction order

1. Create an environment: `pip install -r requirements.txt`.
2. Export `HVI_DATA_DIR=/path/to/authorised/input-data`.
3. Run `python scripts/run_resubmission_analyses.py` for UTCI trends, FDR, Moran's I, denominator scenarios, tables, and map.
4. Run `python scripts/download_era5_heatindex_2025.py`, then `python scripts/calculate_heat_index_2025.py` for the 2025 Heat Index comparator.
5. Run `python scripts/audit_references_crossref.py` and `python scripts/reconcile_references.py` for the reference register.

## Interpretation

Outputs are ecological, subcounty-level screening priorities. They are not causal effects, individual risk predictions, or observed pregnancy counts. The Heat Index analysis is a 2025 cross-metric robustness test and does not replace the 1991-2025 UTCI trend analysis.

## Licences

Code is MIT licensed. Documentation, derived tables, figures, and supplementary outputs are CC BY 4.0. Third-party source data retain their original licences.
