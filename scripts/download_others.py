#!/usr/bin/env python3
"""Placeholder other-data extractor. Produces NA columns for elevation, slope, aspect, ndvi.
Deterministic: reads data/points.csv (site,lon,lat) and writes CSV+Parquet with fixed columns.
"""
import os
import sys
import pandas as pd
import numpy as np

PTS = 'data/points.csv'
OUT_CSV = 'results/other_data.csv'
OUT_PARQ = 'results/other_data.parquet'

if not os.path.exists(PTS):
    raise SystemExit('Missing data/points.csv')
pts = pd.read_csv(PTS)
if 'site' not in pts.columns:
    raise SystemExit("points.csv must contain 'site' column")

out = pd.DataFrame({'site': pts['site']})
# deterministic placeholders (NA) for required columns
for col in ['elevation','slope','aspect','ndvi']:
    out[col] = np.nan

out.to_csv(OUT_CSV, index=False)
try:
    out.to_parquet(OUT_PARQ, index=False)
except Exception:
    pass
print('Wrote', OUT_CSV)
