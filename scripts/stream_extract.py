#!/usr/bin/env python3
"""Streaming extractor: download and process one variable (or one set of bands) at a time.

This keeps peak disk usage small by removing downloaded tiles after extracting their values.

Usage: python scripts/stream_extract.py config/extract_config.yaml results/sampled.parquet
"""
import sys
import os
import yaml
import pandas as pd
import subprocess
from urllib.parse import urlparse

# Ensure package import works when running scripts from repo root without installing the package
_repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
src_path = os.path.join(_repo_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from globalalpine import extract_pipeline as ep


def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def is_url(p):
    return isinstance(p, str) and (p.startswith('http://') or p.startswith('https://'))


def download_to(path_url, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if os.path.exists(out_path):
        return out_path
    print(f"Downloading {path_url} -> {out_path}")
    subprocess.check_call(["curl", "-fSL", "--retry", "3", "-o", out_path, path_url])
    return out_path


def resolve_local(p):
    # if url, map to data/tiles/<basename>, else return as-is
    if is_url(p):
        fn = os.path.basename(p.split('?')[0])
        return os.path.join('data', 'tiles', fn)
    return p


def process_soil_var(df, varname, band_map_cfg, depth_col, chunk_size=5000):
    # Build local band map (download if needed)
    local_map = {}
    downloaded = []
    for band_label, p in band_map_cfg.items():
        local = resolve_local(p)
        if is_url(p):
            download_to(p, local)
            downloaded.append(local)
        local_map[band_label] = local
    # call extraction for this single var
    out = ep.extract_soilgrids_for_points(df, coords_cols=("lon","lat"), depth_col=depth_col, variables_to_band_paths={varname: local_map}, chunk_size=chunk_size)
    # keep only id, coords, and the two columns for var
    keep_cols = ["lon","lat"]
    depth_colname = f"soilgrid_depth_{varname}"
    nodepth_colname = f"soilgrid_nodepth_{varname}"
    part = out[keep_cols + [depth_colname, nodepth_colname]]
    return part, downloaded


def process_singleband_vars(df, var_map_cfg, prefix, chunk_size=5000):
    # var_map_cfg: {name: path}
    parts = {}
    downloaded_all = []
    for varname, p in (var_map_cfg or {}).items():
        local = resolve_local(p)
        if is_url(p):
            download_to(p, local)
            downloaded_all.append(local)
        # sample single band
        submap = {varname: local}
        out = ep.sample_singleband_rasters(df, coords_cols=("lon","lat"), var_map=submap, chunk_size=chunk_size)
        colname = f"{prefix}_{varname}"
        part = out[["lon","lat", colname]]
        parts[varname] = part
    return parts, downloaded_all


def process_dem(df, dem_cfg, chunk_size=5000, dem_window=3):
    downloaded = []
    if dem_cfg is None:
        return None, []
    dem_local = resolve_local(dem_cfg)
    if is_url(dem_cfg):
        download_to(dem_cfg, dem_local)
        downloaded.append(dem_local)
    out = ep.sample_dem_for_df(df, dem_local, coords_cols=("lon","lat"), window=dem_window)
    part = out[["lon","lat","dem_elevation","dem_slope_deg","dem_aspect_deg","dem_tpi","dem_hillshade"]]
    return part, downloaded


def main(argv):
    if len(argv) < 3:
        print("Usage: stream_extract.py <config.yaml> <output_parquet>")
        sys.exit(2)
    cfg = load_config(argv[1])
    out_path = argv[2]
    pts_path = cfg.get('points_csv')
    if pts_path is None:
        raise ValueError('points_csv must be set')
    df = pd.read_csv(pts_path)
    # ensure lon/lat exist
    if not set(['lon','lat']).issubset(df.columns):
        raise ValueError('points CSV must contain lon and lat columns')

    depth_col = cfg.get('depth_col', 'depth_cm')
    # If the configured depth column is missing from the points CSV, disable per-point depth
    if depth_col and depth_col not in df.columns:
        print(f"Warning: depth column '{depth_col}' not found in points CSV; per-point depth disabled.")
        depth_col = None
    chunk_size = cfg.get('chunk_size', 5000)
    dem_window = cfg.get('dem_window', 3)

    # ensure results directories exist up-front to avoid Snakemake latency/parent-dir races
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    parts_dir = os.path.join(out_dir, 'parts') if out_dir else 'results/parts'
    os.makedirs(parts_dir, exist_ok=True)

    part_files = []

    # Soil variables: process one variable (all its bands) at a time
    soil_cfg = cfg.get('soil_vars') or {}
    for varname, band_map in soil_cfg.items():
        print(f"Processing soil variable {varname}")
        part_df, downloaded = process_soil_var(df, varname, band_map, depth_col, chunk_size=chunk_size)
        part_file = os.path.join(parts_dir, f"soil_{varname}.parquet")
        part_df.to_parquet(part_file, index=False)
        part_files.append(part_file)
        # remove downloaded files to free space
        for f in downloaded:
            try:
                os.remove(f)
            except Exception:
                pass

    # CHELSA
    chelsa_cfg = cfg.get('chelsa') or {}
    chelsa_parts, chelsa_downloaded = process_singleband_vars(df, chelsa_cfg, prefix='chelsa', chunk_size=chunk_size)
    for varname, part_df in chelsa_parts.items():
        part_file = os.path.join(parts_dir, f"chelsa_{varname}.parquet")
        part_df.to_parquet(part_file, index=False)
        part_files.append(part_file)
    for f in chelsa_downloaded:
        try:
            os.remove(f)
        except Exception:
            pass

    # DEM
    dem_cfg = cfg.get('dem')
    dem_part, dem_downloaded = process_dem(df, dem_cfg, chunk_size=chunk_size, dem_window=dem_window)
    if dem_part is not None:
        part_file = os.path.join(parts_dir, "dem.parquet")
        dem_part.to_parquet(part_file, index=False)
        part_files.append(part_file)
    for f in dem_downloaded:
        try:
            os.remove(f)
        except Exception:
            pass

    # Landcover
    lc_cfg = cfg.get('landcover') or {}
    lc_parts, lc_downloaded = process_singleband_vars(df, lc_cfg, prefix='lc', chunk_size=chunk_size)
    for varname, part_df in lc_parts.items():
        part_file = os.path.join(parts_dir, f"lc_{varname}.parquet")
        part_df.to_parquet(part_file, index=False)
        part_files.append(part_file)
    for f in lc_downloaded:
        try:
            os.remove(f)
        except Exception:
            pass

    # NDVI
    ndvi_cfg = cfg.get('ndvi') or {}
    ndvi_parts, ndvi_downloaded = process_singleband_vars(df, ndvi_cfg, prefix='ndvi', chunk_size=chunk_size)
    for varname, part_df in ndvi_parts.items():
        part_file = os.path.join(parts_dir, f"ndvi_{varname}.parquet")
        part_df.to_parquet(part_file, index=False)
        part_files.append(part_file)
    for f in ndvi_downloaded:
        try:
            os.remove(f)
        except Exception:
            pass

    # Snow
    snow_cfg = cfg.get('snow') or {}
    snow_parts, snow_downloaded = process_singleband_vars(df, snow_cfg, prefix='snow', chunk_size=chunk_size)
    for varname, part_df in snow_parts.items():
        part_file = os.path.join(parts_dir, f"snow_{varname}.parquet")
        part_df.to_parquet(part_file, index=False)
        part_files.append(part_file)
    for f in snow_downloaded:
        try:
            os.remove(f)
        except Exception:
            pass

    # Merge parts horizontally on lon/lat
    print(f"Merging {len(part_files)} part files into {out_path}")
    merged = df[['lon','lat']].copy()
    for pf in part_files:
        part = pd.read_parquet(pf)
        # avoid duplicate lon/lat columns when merging
        cols = [c for c in part.columns if c not in ['lon','lat']]
        merged = merged.join(part[cols])
    # Save merged
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    merged.to_parquet(out_path, index=False)
    print('Done')


if __name__ == '__main__':
    main(sys.argv)
