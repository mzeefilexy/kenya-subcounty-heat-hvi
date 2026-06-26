"""Rebuild the HVI outputs for a hazard-gated screening framework.

This script preserves the submitted components and creates two transparent
screening scores. It does not claim validation against health outcomes.
"""

from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / ".review_package" / "PLOS_Climate_HVI_resubmission_package_2026-06-21" / "02_source_data_and_archives" / "adm2_HVI_rankings_all_groups.csv"
OUT = ROOT / "revision_workspace" / "outputs"


def rank_descending(series: pd.Series) -> pd.Series:
    """Dense-free rank: tied values retain the same minimum rank."""
    return series.rank(method="min", ascending=False).astype(int)


def top_overlap(a: pd.DataFrame, b: pd.DataFrame, n: int = 20) -> int:
    aa = set(a.nsmallest(n, "rank_current")["adm2_name"])
    bb = set(b.nsmallest(n, "rank_revised")["adm2_name"])
    return len(aa & bb)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SOURCE)

    required = {"adm2_name", "group", "period", "H38_bar", "Haz", "Mag", "PovScore", "HVI", "Rank"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected columns: {sorted(missing)}")

    df = df.copy()
    df["heat_exposed"] = df["H38_bar"] > 0
    # Primary score: preserve submitted weights only when severe heat exists.
    df["priority_hazard_gated"] = df["HVI"].where(df["heat_exposed"], 0.0)
    # Equal-domain sensitivity: heat hazard, population magnitude, and
    # socioeconomic disadvantage contribute equally, but heat remains necessary.
    df["priority_equal_weight_gated"] = (
        (df["Haz"] + df["Mag"] + df["PovScore"]) / 3
    ).where(df["heat_exposed"], 0.0)

    df["rank_hazard_gated"] = (
        df.groupby(["group", "period"])["priority_hazard_gated"]
        .transform(rank_descending)
    )
    df["rank_equal_weight_gated"] = (
        df.groupby(["group", "period"])["priority_equal_weight_gated"]
        .transform(rank_descending)
    )
    # A positive H38 condition is explicit in the output rather than implicit.
    df["screening_class"] = "No severe heat exposure (not ranked for heat priority)"
    df.loc[df["heat_exposed"], "screening_class"] = "Severe-heat screening priority"

    output_columns = [
        "adm2_name", "group", "period", "H38_bar", "H46_bar", "S_share", "RWI_bar",
        "Haz", "Mag", "PovScore", "HVI", "Rank", "heat_exposed", "screening_class",
        "priority_hazard_gated", "rank_hazard_gated",
        "priority_equal_weight_gated", "rank_equal_weight_gated",
    ]
    df[output_columns].to_csv(OUT / "subcounty_priority_hazard_gated_all_groups.csv", index=False)

    # Sensitivity table for each group-period combination.
    records = []
    for (group, period), subset in df.groupby(["group", "period"], sort=True):
        s = subset.copy()
        exposed = s[s["heat_exposed"]].copy()
        # Compare ranks among the heat-exposed units only. This avoids assigning
        # substantive meaning to the ranks of zero-hazard locations.
        exposed["rank_current"] = rank_descending(exposed["HVI"])
        exposed["rank_revised"] = rank_descending(exposed["priority_equal_weight_gated"])
        records.append({
            "group": group,
            "period": period,
            "n_subcounties": len(s),
            "n_heat_exposed": len(exposed),
            "n_zero_hazard": int((~s["heat_exposed"]).sum()),
            "best_zero_hazard_submitted_rank": int(s.loc[~s["heat_exposed"], "Rank"].min()) if (~s["heat_exposed"]).any() else pd.NA,
            "spearman_submitted_vs_equal_weight_among_exposed": exposed["HVI"].rank().corr(exposed["priority_equal_weight_gated"].rank(), method="pearson"),
            "top20_overlap_submitted_vs_equal_weight_among_exposed": top_overlap(exposed, exposed),
        })
    summary = pd.DataFrame(records)
    summary.to_csv(OUT / "priority_index_sensitivity_summary.csv", index=False)

    # The manuscript's primary period is 2025; extract compact reporting tables.
    primary = df[df["period"].astype(str) == "2025"].copy()
    primary.sort_values(["group", "rank_hazard_gated", "adm2_name"]).to_csv(
        OUT / "priority_hazard_gated_2025_rankings.csv", index=False
    )

    print(f"Wrote {len(df):,} rows to {OUT}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
