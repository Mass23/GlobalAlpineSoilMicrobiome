# GlobalAlpineSoilMicrobiome

Python-first tooling for downloading and sampling geospatial environmental
predictors for alpine soil microbiome workflows.

## Scope

- **Python**: all online data gathering, GeoTIFF handling, raster sampling,
  and covariate extraction from GPS coordinates or point tables.
- **R**: downstream microbiome analysis, visualisation, and statistics
  (reads covariate tables exported by this package).

## Supported data sources

| Source | Coverage | Access |
|---|---|---|
| [SoilGrids v2](https://soilgrids.org) | Global, 250 m | ISRIC S3 COGs |
| [CHELSA v2.1](https://chelsa-climate.org) | Global, ~1 km | WSL cloud COGs |
| [Copernicus DEM GLO-30/90](https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model) | Global, 30/90 m | AWS S3 COGs |

## Installation

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.10.  Runtime dependencies: `numpy`, `pandas`, `requests`,
`rasterio`.

## Quick start

### Single point

```python
from globalalpine import (
    resolve_soilgrids_url,
    resolve_chelsa_url,
    resolve_copernicus_dem_url,
    sample_geotiff_at_point,
    get_terrain_from_dem,
    fetch_copernicus_dem_at_point,
)

lon, lat = 9.8, 46.8  # somewhere in the Swiss Alps

# --- SoilGrids: soil pH at 0-5 cm ---
url = resolve_soilgrids_url("phh2o", depth="0-5cm", quantile="mean")
ph = sample_geotiff_at_point(url, lon, lat)
print(ph)  # {1: 6.2}

# --- CHELSA: mean annual temperature (bio1) ---
url = resolve_chelsa_url("bio1")
bio1 = sample_geotiff_at_point(url, lon, lat)
print(bio1)

# --- Copernicus DEM elevation + terrain ---
# Option A: stream directly
elev = fetch_copernicus_dem_at_point(lon, lat, resolution=30)
print(f"elevation: {elev} m")

# Option B: download tile first, then derive slope & aspect
tile_url = resolve_copernicus_dem_url(lon, lat, resolution=30)
# (download to a local cache, then:)
# terrain = get_terrain_from_dem("/path/to/tile.tif", lon, lat)
```

### Batch points from a DataFrame

```python
import pandas as pd
from globalalpine import fetch_soilgrids_at_points, fetch_chelsa_at_points

points = pd.DataFrame({
    "site": ["A", "B", "C"],
    "lon":  [9.8, 8.1, 11.2],
    "lat":  [46.8, 47.3, 46.1],
})

# Add SoilGrids pH column
points = fetch_soilgrids_at_points("phh2o", points, depth="0-5cm", quantile="mean")

# Add CHELSA bio1 column
points = fetch_chelsa_at_points("bio1", points)

print(points)
```

### Export for R

```python
from globalalpine import export_covariates

export_covariates(points, "covariates.csv")          # CSV
export_covariates(points, "covariates.parquet")       # Parquet
```

In R:

```r
library(readr)
covariates <- read_csv("covariates.csv")
```

### Local caching

Pass `cache_dir` to any `fetch_*` function to download GeoTIFFs once and
reuse them:

```python
from globalalpine import fetch_soilgrids_at_point

value = fetch_soilgrids_at_point("soc", lon=9.8, lat=46.8, cache_dir="./cache")
```

## Project layout

```
src/globalalpine/
├── __init__.py          public API re-exports
├── coordinates.py       GPSPoint, BoundingBox, validate_point, validate_bbox
├── raster_extract.py    sample_geotiff_at_point / _at_points
├── soilgrids.py         resolve_soilgrids_url, fetch_soilgrids_at_point/points
├── chelsa.py            resolve_chelsa_url, fetch_chelsa_at_point/points
├── copernicus.py        resolve_copernicus_dem_url, fetch_copernicus_dem_at_point/points
├── terrain.py           get_terrain_from_dem (elevation, slope, aspect)
├── downloads.py         download_file streaming helper
├── export.py            export_covariates (CSV / Parquet)
└── cli.py               globalalpine --lon … --lat …
tests/
├── test_coordinates.py
├── test_raster_extract.py
├── test_terrain.py
└── test_url_resolvers.py
```

## Running tests

```bash
pytest
```

## Development notes

Dataset access logic is separated from raster extraction logic so that
URL templates, credentials, or local file paths can be swapped independently.
All raster I/O goes through `rasterio`; remote COGs are streamed via GDAL's
`/vsicurl/` virtual filesystem without requiring a full download.
