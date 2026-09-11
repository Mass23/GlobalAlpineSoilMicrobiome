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

# Ensure package import works when running scripts from repo root without installing the package
_repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
src_path = os.path.join(_repo_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

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

    # Resolve any http(s) paths in the config to local data/tiles/<basename> so Snakemake can
    # download them and mark them temporary. If a config entry is already a local path, keep it.
    def resolve_path(p):
        import os
        if isinstance(p, str) and (p.startswith('http://') or p.startswith('https://')):
            fn = os.path.basename(p.split('?')[0])
            return os.path.join('data', 'tiles', fn)
        return p

    def resolve_mapping(m):
        if m is None:
            return None
        if isinstance(m, dict):
            out = {}
            for k, v in m.items():
                out[k] = resolve_path(v)
            return out
        return m

    # soil_vars is nested: var -> band_label -> path
    soil_map_cfg = cfg.get('soil_vars') or {}
    soil_map = {}
    for var, band_map in soil_map_cfg.items():
        soil_map[var] = {}
        for band_label, p in band_map.items():
            soil_map[var][band_label] = resolve_path(p)

    chelsa_map_cfg = cfg.get('chelsa')
    chelsa_map = resolve_mapping(chelsa_map_cfg)
    dem_path_cfg = cfg.get('dem')
    dem_path = resolve_path(dem_path_cfg) if dem_path_cfg else None
    landcover_map = resolve_mapping(cfg.get('landcover'))
    ndvi_map = resolve_mapping(cfg.get('ndvi'))
    snow_map = resolve_mapping(cfg.get('snow'))
    depth_col = cfg.get('depth_col', 'depth_cm')

    print(f"Loaded {len(df)} points from {pts_path}")
    print("Resolved soil/chelsa/dem paths; any http URLs should be available under data/tiles/")

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
