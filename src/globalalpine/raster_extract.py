"""Core raster-extraction helpers using *rasterio*.

Functions
---------
sample_geotiff_at_point(path_or_url, lon, lat)
    Return a dict mapping band index → value for one coordinate.

sample_geotiff_at_points(path_or_url, points_df)
    Return a DataFrame with one column per band added to *points_df*.

Both functions accept a local file path **or** a ``/vsicurl/``-prefixed URL so
that rasterio can stream from a remote HTTPS source without a full download.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import rasterio

from .coordinates import validate_point

# rasterio can open remote HTTPS files via GDAL's /vsicurl/ virtual filesystem.
# Callers can pass a plain https:// URL and this module will add the prefix.
_VSICURL_PREFIX = "/vsicurl/"


def _open_path(path_or_url: Union[str, Path]) -> str:
    """Return the path/URL string suitable for *rasterio.open*."""
    s = str(path_or_url)
    if s.startswith("http://") or s.startswith("https://"):
        return _VSICURL_PREFIX + s
    return s


def sample_geotiff_at_point(
    path_or_url: Union[str, Path],
    lon: float,
    lat: float,
) -> dict[int, float]:
    """Return raster values at a single point.

    Parameters
    ----------
    path_or_url:
        Local file path or HTTPS URL to a Cloud-Optimised GeoTIFF (or any
        rasterio-readable raster).
    lon, lat:
        WGS-84 longitude and latitude in decimal degrees.

    Returns
    -------
    dict
        ``{band_index: value}`` for every band in the raster.  Missing-data
        values (``nodata``) are returned as ``float("nan")``.

    Raises
    ------
    ValueError
        If the coordinates are outside the WGS-84 valid range.
    rasterio.errors.RasterioIOError
        If the file cannot be opened.
    """
    validate_point(lon, lat)
    src_path = _open_path(path_or_url)

    with rasterio.open(src_path) as src:
        # sample() returns an iterator of arrays, one per (lon, lat) pair
        values = next(src.sample([(lon, lat)]))
        nodata = src.nodata
        result: dict[int, float] = {}
        for band_idx, val in enumerate(values, start=1):
            fval = float(val)
            if nodata is not None and np.isclose(fval, float(nodata)):
                fval = float("nan")
            result[band_idx] = fval
    return result


def sample_geotiff_at_points(
    path_or_url: Union[str, Path],
    points_df: pd.DataFrame,
    lon_col: str = "lon",
    lat_col: str = "lat",
    band_prefix: str = "band_",
) -> pd.DataFrame:
    """Return a copy of *points_df* with raster band columns appended.

    Parameters
    ----------
    path_or_url:
        Local file path or HTTPS URL to the raster.
    points_df:
        DataFrame that must contain *lon_col* and *lat_col* columns.
    lon_col, lat_col:
        Names of the longitude / latitude columns in *points_df*.
    band_prefix:
        Prefix for the new column names (e.g. ``"band_"`` → ``"band_1"``,
        ``"band_2"``…).

    Returns
    -------
    pd.DataFrame
        Input DataFrame with band columns appended.

    Raises
    ------
    KeyError
        If *lon_col* or *lat_col* are absent from *points_df*.
    ValueError
        If any coordinate row is outside the WGS-84 valid range.
    """
    if lon_col not in points_df.columns:
        raise KeyError(f"Column '{lon_col}' not found in points_df")
    if lat_col not in points_df.columns:
        raise KeyError(f"Column '{lat_col}' not found in points_df")

    for _, row in points_df.iterrows():
        validate_point(float(row[lon_col]), float(row[lat_col]))

    src_path = _open_path(path_or_url)

    with rasterio.open(src_path) as src:
        coords = list(zip(points_df[lon_col].astype(float), points_df[lat_col].astype(float)))
        samples = list(src.sample(coords))
        nodata = src.nodata
        n_bands = src.count

    band_cols: dict[str, list[float]] = {f"{band_prefix}{b}": [] for b in range(1, n_bands + 1)}
    for arr in samples:
        for b_idx, val in enumerate(arr, start=1):
            fval = float(val)
            if nodata is not None and np.isclose(fval, float(nodata)):
                fval = float("nan")
            band_cols[f"{band_prefix}{b_idx}"].append(fval)

    result = points_df.copy()
    for col, vals in band_cols.items():
        result[col] = vals
    return result
