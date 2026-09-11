#!/usr/bin/env python3
"""Runner script used by Snakemake to execute the extraction pipeline.

Usage:
  python scripts/run_extract.py <config_yaml> <output_path>

The config yaml should follow config/extract_config.yaml template in the repo.
"""
import sys
import os
import yaml
import pandas as pd

from globalalpine import extract_pipeline as ep


def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def main(argv):
    if len(argv) < 3:
        print("Usage: run_extract.py <config_yaml> <output_path>")
        sys.exit(2)
    cfg_path = argv[1]
    out_path = argv[2]

    cfg = load_config(cfg_path)
    pts_path = cfg.get('points_csv')
    if pts_path is None:
        raise ValueError('points_csv not set in config')
    df = pd.read_csv(pts_path)

    soil_map = cfg.get('soil_vars')
    chelsa_map = cfg.get('chelsa')
    dem_path = cfg.get('dem')
    landcover_map = cfg.get('landcover')
    ndvi_map = cfg.get('ndvi')
    snow_map = cfg.get('snow')
    depth_col = cfg.get('depth_col', 'depth_cm')

    print(f"Loaded {len(df)} points from {pts_path}")

    out_df = ep.extract_all(
        df,
        soil_vars_map=soil_map,
        chelsa_map=chelsa_map,
        dem_path=dem_path,
        landcover_map=landcover_map,
        ndvi_map=ndvi_map,
        snow_map=snow_map,
        coords_cols=("lon", "lat"),
        depth_col=depth_col,
        chunk_size=cfg.get('chunk_size', 5000),
        dem_window=cfg.get('dem_window', 3),
    )

    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # save as parquet for efficiency
    out_df.to_parquet(out_path, index=False)
    print(f"Wrote output to {out_path}")


if __name__ == '__main__':
    main(sys.argv)
