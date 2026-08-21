"""Copernicus DEM raster access helpers.

The Copernicus DEM (GLO-30 and GLO-90) is distributed as 1°×1° COG tiles on
an AWS S3 bucket.  Tile filenames are derived from the bottom-left corner of
each 1° cell.

Functions
---------
resolve_copernicus_dem_url(lon, lat, resolution)
    Build the tile URL for the 1°×1° tile containing *lon*/*lat*.

fetch_copernicus_dem_at_point(lon, lat, resolution, cache_dir)
    Return the elevation value at a single point.

fetch_copernicus_dem_at_points(points_df, resolution, cache_dir)
    Return a DataFrame with elevation appended.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .coordinates import validate_point
from .raster_extract import sample_geotiff_at_point, sample_geotiff_at_points
from .downloads import download_file

# ---------------------------------------------------------------------------
# URL templates
# ---------------------------------------------------------------------------

_GLO30_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"
_GLO90_BASE = "https://copernicus-dem-90m.s3.amazonaws.com"


def resolve_copernicus_dem_url(
    lon: float,
    lat: float,
    resolution: int = 30,
) -> str:
    """Return the HTTPS URL for the Copernicus DEM tile that covers *lon*/*lat*.

    The Copernicus GLO-30 (30 m) and GLO-90 (90 m) tiles are named after the
    bottom-left corner of each 1°×1° cell, using the convention::

        Copernicus_DSM_COG_{res}_{lat_dir}{lat_abs:02d}_00_{lon_dir}{lon_abs:03d}_00_DEM.tif

    Parameters
    ----------
    lon, lat:
        WGS-84 longitude / latitude in decimal degrees.
    resolution:
        DEM resolution in metres: ``30`` (GLO-30) or ``90`` (GLO-90).

    Returns
    -------
    str
        Full HTTPS URL to the COG tile.

    Raises
    ------
    ValueError
        If *resolution* is not ``30`` or ``90``, or if the coordinates are
        out of range.

    Examples
    --------
    >>> resolve_copernicus_dem_url(9.8, 46.8)
    'https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N46_00_E009_00_DEM/Copernicus_DSM_COG_10_N46_00_E009_00_DEM.tif'
    """
    validate_point(lon, lat)
    if resolution not in (30, 90):
        raise ValueError(f"resolution must be 30 or 90, got {resolution}")

    base = _GLO30_BASE if resolution == 30 else _GLO90_BASE
    res_tag = 10 if resolution == 30 else 30

    # Tile origin is the *floor* of the coordinate
    tile_lat = int(math.floor(lat))
    tile_lon = int(math.floor(lon))

    lat_dir = "N" if tile_lat >= 0 else "S"
    lon_dir = "E" if tile_lon >= 0 else "W"
    lat_abs = abs(tile_lat)
    lon_abs = abs(tile_lon)

    stem = (
        f"Copernicus_DSM_COG_{res_tag}_{lat_dir}{lat_abs:02d}_00"
        f"_{lon_dir}{lon_abs:03d}_00_DEM"
    )
    return f"{base}/{stem}/{stem}.tif"


def fetch_copernicus_dem_at_point(
    lon: float,
    lat: float,
    resolution: int = 30,
    cache_dir: Optional[Union[str, Path]] = None,
) -> float:
    """Return the Copernicus DEM elevation at a single point.

    Parameters
    ----------
    lon, lat:
        WGS-84 longitude / latitude.
    resolution:
        DEM resolution in metres: ``30`` or ``90``.
    cache_dir:
        Optional local directory for caching the downloaded tile.

    Returns
    -------
    float
        Elevation in metres (``float("nan")`` if no-data).
    """
    validate_point(lon, lat)
    url = resolve_copernicus_dem_url(lon, lat, resolution)
    source = _resolve_source(url, cache_dir)
    band_values = sample_geotiff_at_point(source, lon, lat)
    return band_values.get(1, float("nan"))


def fetch_copernicus_dem_at_points(
    points_df: pd.DataFrame,
    resolution: int = 30,
    cache_dir: Optional[Union[str, Path]] = None,
    lon_col: str = "lon",
    lat_col: str = "lat",
) -> pd.DataFrame:
    """Return *points_df* with an elevation column appended.

    Each point may fall in a different 1°×1° tile.  Tiles are grouped so each
    unique tile is opened only once.

    Parameters
    ----------
    points_df:
        DataFrame with longitude / latitude columns.
    resolution:
        DEM resolution in metres: ``30`` or ``90``.
    cache_dir:
        Optional local directory for caching downloaded tiles.
    lon_col, lat_col:
        Names of the longitude / latitude columns.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with column ``"elevation_m"`` appended.
    """
    if lon_col not in points_df.columns:
        raise KeyError(f"Column '{lon_col}' not found in points_df")
    if lat_col not in points_df.columns:
        raise KeyError(f"Column '{lat_col}' not found in points_df")

    result = points_df.copy()
    elevations: list[float] = [float("nan")] * len(result)

    # Group by tile URL to minimise file-open overhead
    urls = [
        resolve_copernicus_dem_url(float(row[lon_col]), float(row[lat_col]), resolution)
        for _, row in points_df.iterrows()
    ]
    result["_tile_url"] = urls

    for url, group in result.groupby("_tile_url"):
        source = _resolve_source(str(url), cache_dir)
        sampled = sample_geotiff_at_points(
            source, group[[lon_col, lat_col]], lon_col=lon_col, lat_col=lat_col
        )
        for idx, val in zip(group.index, sampled["band_1"]):
            elevations[result.index.get_loc(idx)] = float(val)

    result = result.drop(columns=["_tile_url"])
    result["elevation_m"] = elevations
    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_source(url: str, cache_dir: Optional[Union[str, Path]]) -> Union[str, Path]:
    if cache_dir is None:
        return url
    dest = Path(cache_dir) / Path(url).name
    if not dest.exists():
        download_file(url, dest)
    return dest
