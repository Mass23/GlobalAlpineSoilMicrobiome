# GlobalAlpineSoilMicrobiome

Python-first tooling for downloading and sampling geospatial environmental predictors for alpine soil microbiome workflows.

## Scope

- **Python**: all online data gathering, GeoTIFF handling, raster sampling, and covariate extraction from GPS coordinates or point tables.
- **R**: downstream microbiome analysis, visualization, and statistics.

## Supported data sources

- SoilGrids
- CHELSA
- Copernicus DEM-derived terrain products

## Planned capabilities

- Download or reference online raster files
- Sample rasters at points or within bounding boxes
- Derive elevation, slope, and aspect from Copernicus DEM
- Export tidy covariate tables for R

## Project layout

- `src/` for the Python package
- `tests/` for lightweight checks
- `README.md` for usage and examples

## Development notes

The repository is intentionally structured so dataset access logic can be separated from raster extraction logic. This makes it easier to adapt to changing URLs, credentials, or local file workflows.
