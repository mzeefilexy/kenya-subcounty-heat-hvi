"""Reproducible analyses available from the supplied Kenya HVI source files.

Creates: 35-year UTCI trend/FDR results, global Moran's I diagnostics,
pregnancy-denominator scenarios, and hazard-gated ranking tables/figure.
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pydeps"))

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date
from shapely.geometry import Point, shape
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPolygon


ROOT = Path(__file__).resolve().parents[1]
DROPBOX = Path(os.environ["HVI_DATA_DIR"]).expanduser()
PKG = ROOT / ".review_package" / "PLOS_Climate_HVI_resubmission_package_2026-06-21" / "02_source_data_and_archives"
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)


def norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def mk_test(values: np.ndarray) -> tuple[float, float, float]:
    """Tie-corrected Mann-Kendall Z, two-sided p, and Sen slope."""
    x = np.asarray(values, dtype=float)
    n = len(x)
    s = sum(np.sign(x[j] - x[i]) for i in range(n - 1) for j in range(i + 1, n))
    _, counts = np.unique(x, return_counts=True)
    var_s = (n * (n - 1) * (2 * n + 5) - sum(c * (c - 1) * (2 * c + 5) for c in counts)) / 18.0
    z = 0.0 if var_s == 0 else ((s - 1) / math.sqrt(var_s) if s > 0 else ((s + 1) / math.sqrt(var_s) if s < 0 else 0.0))
    p = 2 * (1 - norm_cdf(abs(z)))
    slopes = [(x[j] - x[i]) / (j - i) for i in range(n - 1) for j in range(i + 1, n)]
    return z, p, float(np.median(slopes))


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def load_annual_h38() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # H38/H40/H46 are exceedances of the daily *maximum* UTCI.  Do not use
    # the daily-mean archive: its values cannot reproduce high-heat days.
    files = sorted((DROPBOX / "UTCI_files").glob("UTCI_daily_stats_*_degC.nc"))
    year_parts, value_parts = [], []
    lat = lon = None
    for nc in files:
        with Dataset(nc) as ds:
            lat = np.asarray(ds.variables["lat"][:])
            lon = np.asarray(ds.variables["lon"][:])
            dates = num2date(ds.variables["time"][:], ds.variables["time"].units, only_use_cftime_datetimes=False)
            year_parts.append(np.array([d.year for d in dates]))
            value_parts.append(np.asarray(ds.variables["utci_max_c"][:], dtype=float))
    years = np.concatenate(year_parts)
    utci = np.concatenate(value_parts, axis=0)
    unique_years = np.unique(years)
    h38 = np.stack([(utci[years == y] > 38).sum(axis=0) for y in unique_years])
    h40 = np.stack([(utci[years == y] > 40).sum(axis=0) for y in unique_years])
    h46 = np.stack([(utci[years == y] > 46).sum(axis=0) for y in unique_years])
    return unique_years, lat, lon, h38, h40, h46


def load_polygons():
    geo = json.loads((DROPBOX / "kenya_adm2_boundaries.geojson").read_text(encoding="utf-8"))
    records = []
    for f in geo["features"]:
        records.append((f["properties"]["shapeName"], shape(f["geometry"])))
    return records


def polygon_annual(years, lat, lon, h38_grid, h40_grid, h46_grid, polygons):
    pts = [Point(x, y) for y in lat for x in lon]
    rows = []
    for name, geom in polygons:
        mask = np.array([geom.covers(pt) for pt in pts])
        # Small/urban ADM2 units can be smaller than an ERA5-HEAT grid cell.
        # Use the nearest grid centre only for those units and retain the method
        # flag so this resolution limitation is visible in the output.
        method = "grid-centres within polygon"
        if not mask.any():
            centre = geom.representative_point()
            nearest = int(np.argmin([centre.distance(pt) for pt in pts]))
            mask[nearest] = True
            method = "nearest grid centre (small polygon)"
        for yi, year in enumerate(years):
            rows.append({
                "adm2_name": name, "year": int(year),
                "aggregation_method": method,
                "H38_polygon": float(h38_grid[yi].ravel()[mask].mean()),
                "H40_polygon": float(h40_grid[yi].ravel()[mask].mean()),
                "H46_polygon": float(h46_grid[yi].ravel()[mask].mean()),
            })
    return pd.DataFrame(rows)


def queen_neighbors(polygons):
    names = [n for n, _ in polygons]
    geoms = [g for _, g in polygons]
    nbrs = [[] for _ in geoms]
    for i in range(len(geoms)):
        for j in range(i + 1, len(geoms)):
            if geoms[i].bounds[2] < geoms[j].bounds[0] or geoms[j].bounds[2] < geoms[i].bounds[0] or geoms[i].bounds[3] < geoms[j].bounds[1] or geoms[j].bounds[3] < geoms[i].bounds[1]:
                continue
            if geoms[i].touches(geoms[j]):
                nbrs[i].append(j)
                nbrs[j].append(i)
    return names, nbrs


def morans_i(values, neighbors, permutations=999, seed=20260621):
    x = np.asarray(values, dtype=float)
    z = x - x.mean()
    n = len(x)
    s0 = sum(len(a) for a in neighbors)
    def calc(v):
        zz = v - v.mean()
        den = (zz * zz).sum()
        num = sum(zz[i] * zz[j] / len(neighbors[i]) for i in range(n) for j in neighbors[i] if neighbors[i])
        return (n / n) * num / den  # row-standardised weights have S0=n
    observed = calc(x)
    rng = np.random.default_rng(seed)
    sims = np.array([calc(rng.permutation(x)) for _ in range(permutations)])
    p = (1 + np.sum(np.abs(sims - sims.mean()) >= abs(observed - sims.mean()))) / (permutations + 1)
    return observed, p, float(sims.mean()), float(sims.std(ddof=1)), s0


def rank_desc(s):
    return s.rank(method="min", ascending=False).astype(int)


def top_overlap(a, b, n=20):
    return len(set(a.nsmallest(n, "rank").adm2_name) & set(b.nsmallest(n, "rank").adm2_name))


def pearson(x, y):
    return float(np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1])


def pregnancy_scenarios(panel):
    base = panel[panel.group.eq("pregnancies")].copy()
    pop = pd.read_csv(DROPBOX / "kenya_subcounty_population_and_events_2015-2030.csv")
    pop = pop.rename(columns={"Est_annual_pregs": "preg_base", "total_infants": "infants_proxy"})
    d = base.merge(pop[["adm2_name", "year", "preg_base", "infants_proxy"]], on=["adm2_name", "year"], how="left")
    d["preg_infant_share"] = d.groupby("year").apply(lambda g: g["preg_base"].sum() * g["infants_proxy"] / g["infants_proxy"].sum(), include_groups=False).reset_index(level=0, drop=True)
    d["preg_blended"] = 0.5 * d["preg_base"] + 0.5 * d["preg_infant_share"]
    records, long = [], []
    for year, g in d.groupby("year"):
        for scenario, col in [("WRA-share base", "preg_base"), ("Infant-share", "preg_infant_share"), ("50:50 blended", "preg_blended")]:
            x = g[["adm2_name", "H38", "Haz", "Score_Pov", col]].copy()
            x["Score_N_scenario"] = (x[col] - x[col].min()) / (x[col].max() - x[col].min()) * 100
            x["priority"] = (0.55*x["Haz"] + 0.30*x["Score_N_scenario"] + 0.15*x["Score_Pov"]).where(x.H38 > 0, 0)
            x["rank"] = rank_desc(x.priority)
            x["year"] = year; x["scenario"] = scenario
            long.append(x[["adm2_name", "year", "scenario", "H38", col, "priority", "rank"]].rename(columns={col:"pregnancy_denominator"}))
    long = pd.concat(long, ignore_index=True)
    for year, g in long.groupby("year"):
        b = g[g.scenario.eq("WRA-share base")]
        for scen in ["Infant-share", "50:50 blended"]:
            s = g[g.scenario.eq(scen)]
            m = b.merge(s, on="adm2_name", suffixes=("_base","_scenario"))
            exposed = m[m.H38_base > 0]
            records.append({"year":year,"scenario":scen,"rank_correlation_exposed":pearson(exposed.rank_base, exposed.rank_scenario),"top20_overlap_exposed":top_overlap(b[b.H38>0],s[s.H38>0]),"national_total_ratio":s.pregnancy_denominator.sum()/b.pregnancy_denominator.sum()})
    return long, pd.DataFrame(records)


def make_map(polygons, gated):
    geom = dict(polygons)
    groups = [("pregnancies", "Pregnant women"), ("u5", "Children under 5"), ("age60plus", "Adults 60+")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), constrained_layout=True)
    for ax, (code, label) in zip(axes, groups):
        d = gated[(gated.group.eq(code)) & (gated.period.astype(str).eq("2025"))].copy()
        values = d.set_index("adm2_name")["priority_hazard_gated"].to_dict()
        patches, colors = [], []
        for name, g in polygons:
            geoms = [g] if g.geom_type == "Polygon" else list(g.geoms)
            for gg in geoms:
                patches.append(MplPolygon(np.asarray(gg.exterior.coords), closed=True))
                colors.append(values.get(name, np.nan))
        pc = PatchCollection(patches, cmap="YlOrRd", edgecolor="#666666", linewidth=0.15)
        pc.set_array(np.asarray(colors)); pc.set_clim(0, max(gated.priority_hazard_gated))
        ax.add_collection(pc); ax.autoscale_view(); ax.set_aspect("equal"); ax.set_axis_off(); ax.set_title(label, fontsize=12, weight="bold")
        for _, row in d.nsmallest(20,"rank_hazard_gated").iterrows():
            c = geom[row.adm2_name].representative_point()
            ax.text(c.x,c.y,str(int(row.rank_hazard_gated)),fontsize=4.5,ha="center",va="center")
    cb=fig.colorbar(pc,ax=axes,shrink=.7,pad=.01); cb.set_label("Hazard-gated screening-priority score")
    fig.suptitle("Hazard-gated subcounty screening priorities, Kenya 2025",fontsize=15,weight="bold")
    fig.savefig(OUT / "Figure_hazard_gated_screening_priorities_2025.png",dpi=350,bbox_inches="tight")
    plt.close(fig)


def main():
    years, lat, lon, h38_grid, h40_grid, h46_grid = load_annual_h38()
    polygons = load_polygons()
    annual = polygon_annual(years, lat, lon, h38_grid, h40_grid, h46_grid, polygons)
    annual.to_csv(OUT / "utci_polygon_annual_threshold_sensitivity_1991_2025.csv", index=False)

    # 35-year MK trends and multiplicity correction.
    trend_rows=[]
    for name,g in annual.groupby("adm2_name"):
        z,p,slope=mk_test(g.sort_values("year").H38_polygon.values)
        trend_rows.append({"adm2_name":name,"mk_z":z,"p_unadjusted":p,"sen_slope_days_per_year":slope})
    trends=pd.DataFrame(trend_rows); trends["p_bh_fdr"]=bh_fdr(trends.p_unadjusted.values); trends["significant_bh_fdr_0_05"]=trends.p_bh_fdr<.05
    trends.to_csv(OUT / "subcounty_H38_trends_1991_2025_bh_fdr.csv",index=False)

    # Spatial autocorrelation for 2025 polygon H38 and the revised priority scores.
    names,nbrs=queen_neighbors(polygons)
    h2025=annual[annual.year.eq(2025)].set_index("adm2_name").reindex(names).H38_polygon.values
    gated=pd.read_csv(OUT / "subcounty_priority_hazard_gated_all_groups.csv")
    spatial=[]
    for label,values in [("H38 polygon mean",h2025)]:
        I,p,simmean,simstd,edges=morans_i(values,nbrs); spatial.append({"outcome":label,"year":2025,"morans_I":I,"permutation_p_two_sided":p,"permutations":999,"queen_links_directed":edges,"null_mean":simmean,"null_sd":simstd})
    for group in ["pregnancies","u5","age60plus"]:
        vals=gated[(gated.group.eq(group))&(gated.period.astype(str).eq("2025"))].set_index("adm2_name").reindex(names).priority_hazard_gated.values
        I,p,simmean,simstd,edges=morans_i(vals,nbrs); spatial.append({"outcome":"Hazard-gated priority: "+group,"year":2025,"morans_I":I,"permutation_p_two_sided":p,"permutations":999,"queen_links_directed":edges,"null_mean":simmean,"null_sd":simstd})
    pd.DataFrame(spatial).to_csv(OUT / "spatial_autocorrelation_moransI_2025.csv",index=False)

    # Threshold sensitivity: same UTCI metric, explicitly not a different heat metric.
    threshold = annual[annual.year.between(2015,2025)].groupby("adm2_name")[["H38_polygon","H40_polygon","H46_polygon"]].mean().reset_index()
    threshold.to_csv(OUT / "utci_threshold_sensitivity_mean_2015_2025.csv",index=False)

    pl,ps = pregnancy_scenarios(pd.read_csv(PKG / "subcounty_HVI_2015_2025_all_groups.csv"))
    pl.to_csv(OUT / "pregnancy_denominator_scenarios_all_years.csv",index=False)
    ps.to_csv(OUT / "pregnancy_denominator_scenario_sensitivity.csv",index=False)

    primary=gated[(gated.period.astype(str)=="2025")].copy()
    primary.sort_values(["group","rank_hazard_gated","adm2_name"]).to_csv(OUT / "Table_hazard_gated_priority_rankings_2025.csv",index=False)
    primary.groupby("group",group_keys=False).apply(lambda x:x.nsmallest(20,"rank_hazard_gated"),include_groups=False).to_csv(OUT / "Table_top20_hazard_gated_priority_2025.csv",index=False)
    make_map(polygons,gated)
    print("Completed analysis outputs in",OUT)
    print("FDR significant trends:",int(trends.significant_bh_fdr_0_05.sum()))
    print(pd.DataFrame(spatial).to_string(index=False))
    print(ps[ps.year.eq(2025)].to_string(index=False))


if __name__ == "__main__":
    main()
