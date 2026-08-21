"""GlobalAlpine – Python utilities for alpine soil microbiome covariates.

Public API
----------
sample_geotiff_at_point(path_or_url, lon, lat)
    Return the raster value(s) at a single coordinate.

sample_geotiff_at_points(path_or_url, points_df)
    Return a DataFrame of raster values for every row in *points_df*.

resolve_soilgrids_url(property_name, depth, quantile)
    Build the WCS download URL for a SoilGrids layer.

resolve_chelsa_url(variable, period)
    Build the CHELSA download URL for a climate variable.

resolve_copernicus_dem_url(lon, lat, resolution)
    Build the Copernicus DEM tile URL that covers a given point.

get_terrain_from_dem(dem_path, lon, lat)
    Derive elevation, slope, and aspect at a point from a local DEM GeoTIFF.
"""

from .coordinates import GPSPoint, BoundingBox, validate_point, validate_bbox
from .raster_extract import sample_geotiff_at_point, sample_geotiff_at_points
from .soilgrids import resolve_soilgrids_url, fetch_soilgrids_at_point, fetch_soilgrids_at_points
from .chelsa import resolve_chelsa_url, fetch_chelsa_at_point, fetch_chelsa_at_points
from .copernicus import resolve_copernicus_dem_url, fetch_copernicus_dem_at_point, fetch_copernicus_dem_at_points
from .terrain import get_terrain_from_dem
from .export import export_covariates

__all__ = [
    "GPSPoint",
    "BoundingBox",
    "validate_point",
    "validate_bbox",
    "sample_geotiff_at_point",
    "sample_geotiff_at_points",
    "resolve_soilgrids_url",
    "fetch_soilgrids_at_point",
    "fetch_soilgrids_at_points",
    "resolve_chelsa_url",
    "fetch_chelsa_at_point",
    "fetch_chelsa_at_points",
    "resolve_copernicus_dem_url",
    "fetch_copernicus_dem_at_point",
    "fetch_copernicus_dem_at_points",
    "get_terrain_from_dem",
    "export_covariates",
]
