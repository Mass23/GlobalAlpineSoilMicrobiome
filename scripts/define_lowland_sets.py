#!/usr/bin/env python3
"""
Lowland replicates: random draws of lowland samples that are EVENLY SPACED in latitude / longitude, repeated n_rep times.

Run     python3 define_lowland_sets.py <sample_table.tsv> <other_data.csv> <lowland_sets.csv>
                [--seed 23] [--ratio 1.0] [--cell-deg 5] [--n-rep 20] [--thin-m 100]
Reads   sample_table.tsv   tab-separated (no quoting); needs columns Sample_acc and Alpine ("Alpine" / "Lowland")
        other_data.csv     one row per point; needs site, lat, lon, elevation_m (written by download_others.py)
Writes  lowland_sets.csv   one row per SAMPLE (site): group, lat, lon, elevation_m, one yes/no column per replicate
                           (rep_01 ... rep_<n_rep>), n_rep_selected and subset_any
        lowland_sets_summary.csv   what every replicate looks like (cells covered, latitude, elevation), the overlap between
                                   replicates and the size of their union, for checking

One replicate = every alpine sample + n lowland samples, n = round(ratio x number of alpine samples), drawn like this:
  * the lowland samples are spread as evenly as possible over a grid of cell-deg x cell-deg latitude / longitude cells: one
    sample per cell in turn (cells with few samples give all they have), so densely sampled regions (Europe, North America)
    stop dominating. WHICH sample of a cell is taken, and the order of the cells in the last, partial round, is random;
  * no two drawn lowland samples are closer than thin-m metres (and none share coordinates): a sample that lies within thin-m
    of one already drawn is skipped for the next one of its cell (thin-m 0 = only identical coordinates are skipped);
  * every replicate has its own seed-controlled random stream, so the replicates differ (a region with few samples gives the
    same samples in every replicate; a densely sampled one gives different ones) and the whole table is reproducible.
All alpine samples are "yes" in every replicate. subset_any = alpine, or lowland in at least one replicate (this is the set of
samples that the soil download fetches: the union of the replicates, see the summary for its size).
n_rep_selected = in how many replicates the sample was drawn (for an alpine sample: n_rep).
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R_EARTH = 6371008.8                                   # metres


def xyz(lat, lon):
    """3D coordinates on the globe (metres): exact distances everywhere, also near the poles and the dateline"""
    la = np.radians(np.asarray(lat, float)); lo = np.radians(np.asarray(lon, float))
    return (R_EARTH * np.cos(la) * np.cos(lo), R_EARTH * np.cos(la) * np.sin(lo), R_EARTH * np.sin(la))


def draw_even(df: pd.DataFrame, n: int, rng, cell_deg: float, thin_m: float) -> set:
    """n sites as evenly as possible over lat/lon cells (round-robin, one per cell and round), none within thin_m of another"""
    n = min(n, len(df))
    if n <= 0:
        return set()
    cell = (pd.Series(np.floor(df["lat"].to_numpy() / cell_deg).astype(int).astype(str), index=df.index) + "_" +
            pd.Series(np.floor(df["lon"].to_numpy() / cell_deg).astype(int).astype(str), index=df.index))
    d = df.assign(cell=cell).iloc[rng.permutation(len(df))].copy()      # random order inside every cell
    d["round"] = d.groupby("cell").cumcount()                           # 0 = the first site of its cell, 1 = the second, ...
    cells = d["cell"].unique()
    pos = dict(zip(rng.permutation(cells), range(len(cells))))          # random cell order in a partial round
    d["cell_pos"] = d["cell"].map(pos)
    d = d.sort_values(["round", "cell_pos"], kind="stable")
    if thin_m is None or thin_m < 0:
        return set(d["site"].iloc[:n])
    # walk through the candidates in that order; skip a site that lies within thin_m of one already accepted
    X, Y, Z = xyz(d["lat"], d["lon"])
    r2 = float(max(thin_m, 0.0)) ** 2
    cs = max(float(thin_m), 1.0)                                         # hash cells of the radius' size: only 27 neighbours to check
    seen: dict = {}
    keep = []
    for i in range(len(d)):
        cx, cy, cz = int(X[i] // cs), int(Y[i] // cs), int(Z[i] // cs)
        too_close = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (px, py, pz) in seen.get((cx + dx, cy + dy, cz + dz), ()):
                        if (X[i] - px) ** 2 + (Y[i] - py) ** 2 + (Z[i] - pz) ** 2 <= r2:
                            too_close = True
                            break
                    if too_close:
                        break
                if too_close:
                    break
            if too_close:
                break
        if too_close:
            continue
        seen.setdefault((cx, cy, cz), []).append((X[i], Y[i], Z[i]))
        keep.append(i)
        if len(keep) >= n:
            break
    return set(d["site"].iloc[keep])


def read_sites(sample_tsv: Path, other_csv: Path) -> pd.DataFrame:
    tab = pd.read_csv(sample_tsv, sep="\t", quoting=csv.QUOTE_NONE, low_memory=False, usecols=["Sample_acc", "Alpine"])
    group = tab.groupby("Sample_acc")["Alpine"].apply(lambda g: "Alpine" if (g == "Alpine").any() else "Lowland")
    sites = group.rename("group").rename_axis("site").reset_index()
    oth = pd.read_csv(other_csv, usecols=lambda c: c in {"site", "lat", "lon", "elevation_m"})
    missing = {"site", "lat", "lon", "elevation_m"} - set(oth.columns)
    if missing:
        raise ValueError(f"{other_csv} lacks column(s) {sorted(missing)}")
    per_site = oth.groupby("site").agg(lat=("lat", "first"), lon=("lon", "first"),
                                       elevation_m=("elevation_m", "median")).reset_index()
    return sites.merge(per_site, on="site", how="left")


def summarise(name: str, d: pd.DataFrame, cell_deg: float) -> dict:
    cells = (np.floor(d["lat"] / cell_deg).astype("Int64").astype(str) + "_" + np.floor(d["lon"] / cell_deg).astype("Int64").astype(str))
    e = d["elevation_m"].dropna()
    q = e.quantile([0.05, 0.5, 0.95]) if len(e) else pd.Series([np.nan] * 3)
    return {"set": name, "n": len(d), "n_cells": cells.nunique(),
            "share_in_top_cell": round(cells.value_counts(normalize=True).iloc[0], 3) if len(d) else np.nan,
            "lat_min": d["lat"].min(), "lat_max": d["lat"].max(), "lon_min": d["lon"].min(), "lon_max": d["lon"].max(),
            "elev_q05": q.iloc[0], "elev_median": q.iloc[1], "elev_q95": q.iloc[2]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_table"); ap.add_argument("other_data"); ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--ratio", type=float, default=1.0, help="lowland samples per alpine sample")
    ap.add_argument("--cell-deg", type=float, default=5.0, help="lat/lon cell size of the even draws (degrees)")
    ap.add_argument("--n-rep", type=int, default=20, help="number of replicates")
    ap.add_argument("--thin-m", type=float, default=100.0, help="no two drawn lowland samples closer than this (metres); 0 = only identical coordinates")
    a = ap.parse_args()

    sites = read_sites(Path(a.sample_table), Path(a.other_data))
    alpine = sites[sites["group"] == "Alpine"]
    lowland = sites[(sites["group"] == "Lowland") & sites["lat"].notna() & sites["lon"].notna()].sort_values("site").reset_index(drop=True)
    if alpine.empty:
        raise ValueError("no alpine sample found")
    n = int(round(a.ratio * len(alpine)))
    print(f"{len(alpine)} alpine samples; {len(lowland)} lowland candidates; {a.n_rep} replicates of n = {n} lowland samples "
          f"(even over {a.cell_deg:g}-degree cells, at most one per {a.thin_m:g} m)")
    lost = int(((sites["group"] == "Lowland") & (sites["lat"].isna() | sites["lon"].isna())).sum())
    if lost:
        print(f"  note: {lost} lowland sample(s) have no coordinates in {a.other_data} and are ignored")
    if len(lowland) < n:
        print(f"  note: only {len(lowland)} lowland candidates: each replicate takes all of them that pass the distance rule")

    out = sites.copy()
    draws = []
    for r in range(1, a.n_rep + 1):
        chosen = draw_even(lowland, n, np.random.default_rng([a.seed, 2, r]), a.cell_deg, a.thin_m)
        draws.append(chosen)
        out[f"rep_{r:02d}"] = np.where((out["group"] == "Alpine") | out["site"].isin(chosen), "yes", "no")
    rep_cols = [f"rep_{r:02d}" for r in range(1, a.n_rep + 1)]
    out["n_rep_selected"] = (out[rep_cols] == "yes").sum(axis=1)
    out["subset_any"] = np.where((out["group"] == "Alpine") | (out["n_rep_selected"] > 0), "yes", "no")
    out = out[["site", "group", "lat", "lon", "elevation_m"] + rep_cols + ["n_rep_selected", "subset_any"]]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False, na_rep="NA")

    rows = [summarise("alpine (reference)", alpine, a.cell_deg)]
    for r, chosen in enumerate(draws, 1):
        rows.append(summarise(f"lowland rep_{r:02d}", lowland[lowland["site"].isin(chosen)], a.cell_deg))
    union = set().union(*draws) if draws else set()
    rows.append(summarise("lowland union of the replicates", lowland[lowland["site"].isin(union)], a.cell_deg))
    rows.append(summarise("all lowland candidates", lowland, a.cell_deg))
    summary = pd.DataFrame(rows)
    summary.to_csv(Path(a.out).with_name("lowland_sets_summary.csv"), index=False)

    jac = [len(draws[i] & draws[j]) / max(len(draws[i] | draws[j]), 1) for i in range(len(draws)) for j in range(i + 1, len(draws))]
    in_all = len(set.intersection(*draws)) if draws else 0
    pd.set_option("display.width", 200)
    print(summary.iloc[[0, 1, 2, len(rows) - 2, len(rows) - 1]].round(1).to_string(index=False))
    print(f"replicates: mean overlap between two replicates (Jaccard) {np.mean(jac):.2f}; union {len(union)} lowland samples; "
          f"{in_all} lowland samples are in every replicate")
    print(f"samples with soil data to fetch (alpine + union of the replicates): {int((out['subset_any'] == 'yes').sum())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

