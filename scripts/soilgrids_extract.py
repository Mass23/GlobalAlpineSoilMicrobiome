#!/usr/bin/env python3
"""Sample SoilGrids VRTs at points in data/points.csv and write results/results_soilgrids.csv.
Best-effort implementation using rasterio to open remote VRTs listed in config/extract_config.yaml.

Notes:
- Expects config/extract_config.yaml to contain a mapping `soil_vars` like in the repo.
- Points CSV must contain columns: site, lon, lat, optional depth_cm.
- This is a best-effort parser; remote VRT access depends on GDAL network support in the conda env.
"""
import os
import sys
import yaml
import pandas as pd
import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT

CFG = 'config/extract_config.yaml'
PTS = 'data/points.csv'
OUT = 'results/soilgrids_sampled.csv'

if not os.path.exists(CFG):
    raise SystemExit('Missing config/extract_config.yaml')
if not os.path.exists(PTS):
    raise SystemExit('Missing data/points.csv')

cfg = yaml.safe_load(open(CFG))
soil_vars = cfg.get('soil_vars', {})
if not soil_vars:
    raise SystemExit('No soil_vars found in config/extract_config.yaml')

pts = pd.read_csv(PTS)
if 'site' not in pts.columns:
    raise SystemExit("points.csv must contain a 'site' column")
if not {'lon','lat'}.issubset(pts.columns):
    raise SystemExit("points.csv must contain 'lon' and 'lat' columns")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

# helper: sample a single raster (url) at lon/lat list, return values array
def sample_raster(url, lon, lat):
    try:
        with rasterio.Env():
            with rasterio.open(url) as src:
                # ensure same crs, use src.sample on lon/lat in src CRS
                # src.sample expects [(x,y), ...] in src CRS (likely EPSG:4326 here)
                coords = list(zip(lon, lat))
                vals = [v[0] if v is not None else np.nan for v in src.sample(coords)]
                return np.array(vals, dtype=float)
    except Exception as e:
        # return nans if cannot sample
        print(f'Warning: failed to open/sample {url}: {e}', file=sys.stderr)
        return np.full(len(lon), np.nan)

# thicknesses match config ordering in extract_config.yaml
thicknesses = [5,10,15,30,40,100]

results = pd.DataFrame({'site': pts['site']})
lon = pts['lon'].values
lat = pts['lat'].values

for var, bands in soil_vars.items():
    # expect bands ordering keys: "0-5","5-15","15-30","30-60","60-100","100-200"
    band_keys = ['0-5','5-15','15-30','30-60','60-100','100-200']
    urls = [bands.get(k) for k in band_keys]
    # sample each band
    band_vals = np.vstack([sample_raster(u, lon, lat) if u else np.full(len(lon), np.nan) for u in urls])
    # compute nodepth (0-30 cm) thickness-weighted average of first 3 bands
    w = np.array(thicknesses[:3]) / 30.0
    nodepth = np.nansum(band_vals[:3,:] * w[:,None], axis=0) / np.nansum(np.where(np.isnan(band_vals[:3,:]),0,w) ,axis=0)
    # where all three are nan, result should be nan
    nodepth[np.all(np.isnan(band_vals[:3,:]),axis=0)] = np.nan
    results[f'soilgrid_nodepth_{var}'] = nodepth
    # if depth provided, do linear interpolation using band midpoints
    if 'depth_cm' in pts.columns:
        mids = np.array([(0+5)/2,(5+15)/2,(15+30)/2,(30+60)/2,(60+100)/2,(100+200)/2])
        depth_vals = []
        for i, d in enumerate(pts['depth_cm'].values):
            if np.isnan(d):
                depth_vals.append(np.nan)
                continue
            # use available (non-nan) bands up to the band containing depth
            valid = ~np.isnan(band_vals[:,i])
            if not np.any(valid):
                depth_vals.append(np.nan); continue
            # if d <=30, prefer the first three bands
            if d <= 30:
                use_idx = np.where(valid[:3])[0]
                if len(use_idx)==0:
                    depth_vals.append(np.nan); continue
                x = mids[use_idx]
                y = band_vals[use_idx, i]
            else:
                # include bands up to containing band
                cum_th = np.cumsum(thicknesses)
                band_idx = np.where(cum_th >= d)[0]
                if len(band_idx)==0:
                    band_idx = len(thicknesses)-1
                max_idx = band_idx[0]
                use_idx = np.where(valid[:max_idx+1])[0]
                if len(use_idx)==0:
                    depth_vals.append(np.nan); continue
                x = mids[use_idx]
                y = band_vals[use_idx, i]
            # linear approx
            try:
                val = np.interp(d, x, y)
            except Exception:
                val = np.nan
            depth_vals.append(float(val))
        results[f'soilgrid_depth_{var}'] = depth_vals
    else:
        results[f'soilgrid_depth_{var}'] = [np.nan]*len(pts)

results.to_csv(OUT, index=False)
# write parquet as well for downstream merging
try:
    results.to_parquet('results/soilgrids_sampled.parquet', index=False)
    print('Wrote', OUT, 'and results/soilgrids_sampled.parquet')
except Exception as e:
    print('Warning: failed to write parquet:', e, file=sys.stderr)
    print('Wrote', OUT)
