"""Focused 2025 Heat Index sensitivity analysis for the hazard-gated priorities."""
from __future__ import annotations

import csv
import json
import math
import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "pydeps"))
import numpy as np
from netCDF4 import Dataset, num2date
from shapely.geometry import Point, shape

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "era5_sensitivity"
OUT = ROOT / "results"
BOUNDARY = Path(os.environ["HVI_DATA_DIR"]).expanduser() / "kenya_adm2_boundaries.geojson"
BASE = OUT / "subcounty_priority_hazard_gated_all_groups.csv"


def heat_index_f(temp_f, rh):
    """NWS Heat Index algorithm; humidity is paired with each hourly temperature."""
    simple = 0.5 * (temp_f + 61.0 + (temp_f - 68.0) * 1.2 + rh * 0.094)
    hi = np.where((simple + temp_f) / 2 < 80, simple, 0.0)
    use = (simple + temp_f) / 2 >= 80
    t, r = temp_f[use], rh[use]
    hi[use] = (-42.379 + 2.04901523*t + 10.14333127*r - .22475541*t*r
               - .00683783*t*t - .05481717*r*r + .00122874*t*t*r
               + .00085282*t*r*r - .00000199*t*t*r*r)
    low = use & (rh < 13) & (temp_f > 80) & (temp_f < 112)
    hi[low] -= ((13-rh[low])/4) * np.sqrt((17-np.abs(temp_f[low]-95))/17)
    high = use & (rh > 85) & (temp_f > 80) & (temp_f < 87)
    hi[high] += ((rh[high]-85)/10) * ((87-temp_f[high])/5)
    return hi


def average_desc_ranks(values):
    order = sorted(range(len(values)), key=lambda i: -values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def main():
    daily = []
    lat = lon = None
    for month in range(1, 13):
        path = DATA / ("era5_kenya_hourly_2025_01_extracted.nc" if month == 1 else f"era5_kenya_hourly_2025_{month:02d}.nc")
        with Dataset(path) as ds:
            temp_c = np.asarray(ds["t2m"][:], dtype=float) - 273.15
            dew_c = np.asarray(ds["d2m"][:], dtype=float) - 273.15
            lat, lon = np.asarray(ds["latitude"][:]), np.asarray(ds["longitude"][:])
            dates = num2date(ds["valid_time"][:], ds["valid_time"].units, only_use_cftime_datetimes=False)
        rh = np.clip(100 * np.exp((17.625*dew_c)/(243.04+dew_c) - (17.625*temp_c)/(243.04+temp_c)), 0, 100)
        hi_f = heat_index_f(temp_c * 1.8 + 32, rh)
        days = np.array([date.strftime("%Y-%m-%d") for date in dates])
        daily.extend(np.nanmax(hi_f[days == day], axis=0) for day in np.unique(days))
    daily = np.stack(daily)

    geo = json.loads(BOUNDARY.read_text(encoding="utf-8"))
    points = [Point(x, y) for y in lat for x in lon]
    hi = {}
    for feature in geo["features"]:
        name, geom = feature["properties"]["shapeName"], shape(feature["geometry"])
        mask = np.array([geom.covers(point) for point in points])
        method = "grid-centres within polygon"
        if not mask.any():
            centre = geom.representative_point()
            mask[np.argmin([centre.distance(point) for point in points])] = True
            method = "nearest grid centre (small polygon)"
        values_f = daily.reshape(365, -1)[:, mask].mean(axis=1)
        hi[name] = {
            "aggregation_method": method,
            "mean_daily_max_heat_index_c": float(((values_f-32)/1.8).mean()),
            "p95_daily_max_heat_index_c": float(np.quantile((values_f-32)/1.8, .95)),
            "annual_days_dailymax_heat_index_ge_32c": int(np.sum(values_f >= 89.6)),
        }

    with BASE.open(encoding="utf-8", newline="") as fh:
        base = [row for row in csv.DictReader(fh) if row["period"] == "2025"]
    rows, summary = [], []
    for group in sorted({row["group"] for row in base}):
        subset = [row for row in base if row["group"] == group]
        days = [hi[row["adm2_name"]]["annual_days_dailymax_heat_index_ge_32c"] for row in subset]
        lo, high = min(days), max(days)
        new_scores = [0 if day == 0 else (day-lo)*100/(high-lo) for day in days]
        new_priority = [0 if day == 0 else .55*score + .30*float(row["S_share"]) + .15*float(row["PovScore"])
                        for row, day, score in zip(subset, days, new_scores)]
        old_priority = [float(row["priority_hazard_gated"]) for row in subset]
        old_ranks, new_ranks = average_desc_ranks(old_priority), average_desc_ranks(new_priority)
        rho = float(np.corrcoef(old_ranks, new_ranks)[0, 1])
        old_top = {subset[i]["adm2_name"] for i in sorted(range(len(subset)), key=lambda i: (-old_priority[i], subset[i]["adm2_name"]))[:20]}
        new_top = {subset[i]["adm2_name"] for i in sorted(range(len(subset)), key=lambda i: (-new_priority[i], subset[i]["adm2_name"]))[:20]}
        summary.append({"group": group, "year": 2025, "heat_index_exposed_subcounties": sum(day > 0 for day in days),
                        "rank_correlation_with_utci_priority": rho, "top20_overlap_with_utci_priority": len(old_top & new_top)})
        for row, day, score, priority, rank in zip(subset, days, new_scores, new_priority, new_ranks):
            metric = hi[row["adm2_name"]]
            rows.append({"adm2_name": row["adm2_name"], "group": group, "year": 2025, **metric,
                         "heat_index_hazard_score": score, "priority_heat_index_gated": priority,
                         "rank_heat_index_gated": rank, "priority_utci_gated": float(row["priority_hazard_gated"]),
                         "rank_utci_gated": float(row["rank_hazard_gated"])})
    for filename, fieldnames, records in [
        ("heat_index_2025_subcounty_metric_sensitivity.csv", list(rows[0]), rows),
        ("heat_index_2025_priority_robustness_summary.csv", list(summary[0]), summary),
    ]:
        with (OUT / filename).open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames); writer.writeheader(); writer.writerows(records)
    print(summary)


if __name__ == "__main__":
    main()
