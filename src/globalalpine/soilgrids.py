"""SoilGrids raster access helpers.

SoilGrids v2 distributes Cloud-Optimised GeoTIFFs (COGs) through the ISRIC
S3 bucket.  The URL template is::

    https://files.isric.org/soilgrids/latest/data/{property}/{property}_{depth}_{quantile}_250m_s0..0cm_1950..2017_v0.2.tif

where valid values are documented at https://www.isric.org/explore/soilgrids.

Functions
---------
resolve_soilgrids_url(property_name, depth, quantile)
    Return the HTTPS URL string for a SoilGrids COG layer.

fetch_soilgrids_at_point(property_name, lon, lat, depth, quantile, cache_dir)
    Return the raster value at a single point, downloading if needed.

fetch_soilgrids_at_points(property_name, points_df, depth, quantile, cache_dir)
    Return a DataFrame with the SoilGrids band column appended.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .coordinates import validate_point
from .raster_extract import sample_geotiff_at_point, sample_geotiff_at_points
from .downloads import download_file

# ---------------------------------------------------------------------------
# URL templates
# ---------------------------------------------------------------------------

_BASE = "https://files.isric.org/soilgrids/latest/data"

# Supported property names (non-exhaustive)
PROPERTIES = [
    "bdod",   # Bulk density
    "cec",    # Cation exchange capacity
    "cfvo",   # Coarse fragments
    "clay",   # Clay content
    "nitrogen",
    "ocd",    # Organic carbon density
    "ocs",    # Organic carbon stocks
    "phh2o",  # pH (water)
    "sand",   # Sand content
    "silt",   # Silt content
    "soc",    # Soil organic carbon
    "wv0010", # Volumetric water content (pF1)
    "wv0033", # Volumetric water content (pF2)
    "wv1500", # Volumetric water content (pF4.2)
]

# Standard depth intervals
DEPTHS = ["0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"]

# Standard quantiles / summary statistics
QUANTILES = ["Q0.05", "Q0.5", "Q0.95", "mean", "uncertainty"]


def resolve_soilgrids_url(
    property_name: str,
    depth: str = "0-5cm",
    quantile: str = "mean",
) -> str:
    """Return the HTTPS URL for a SoilGrids COG layer.

    Parameters
    ----------
    property_name:
        SoilGrids property identifier (e.g. ``"phh2o"``, ``"soc"``).
    depth:
        Depth interval string, e.g. ``"0-5cm"``.
    quantile:
        Quantile or summary statistic, e.g. ``"mean"``, ``"Q0.5"``.

    Returns
    -------
    str
        Full HTTPS URL to the Cloud-Optimised GeoTIFF.

    Examples
    --------
    >>> resolve_soilgrids_url("phh2o", depth="0-5cm", quantile="mean")
    'https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_0-5cm_mean_250m.tif'
    """
    filename = f"{property_name}_{depth}_{quantile}_250m.tif"
    return f"{_BASE}/{property_name}/{filename}"


def fetch_soilgrids_at_point(
    property_name: str,
    lon: float,
    lat: float,
    depth: str = "0-5cm",
    quantile: str = "mean",
    cache_dir: Optional[Union[str, Path]] = None,
) -> float:
    """Return the SoilGrids value at a single point.

    If *cache_dir* is provided the GeoTIFF is downloaded once and reused.
    Without *cache_dir* the raster is streamed via ``/vsicurl/`` (requires
    GDAL network access).

    Parameters
    ----------
    property_name:
        SoilGrids property identifier (e.g. ``"phh2o"``, ``"soc"``).
    lon, lat:
        WGS-84 longitude / latitude in decimal degrees.
    depth:
        Depth interval string.
    quantile:
        Quantile or summary statistic.
    cache_dir:
        Optional local directory for caching the downloaded GeoTIFF.

    Returns
    -------
    float
        Raster value at the point (``float("nan")`` if no-data).
    """
    validate_point(lon, lat)
    url = resolve_soilgrids_url(property_name, depth, quantile)
    source = _resolve_source(url, cache_dir)
    band_values = sample_geotiff_at_point(source, lon, lat)
    return band_values.get(1, float("nan"))


def fetch_soilgrids_at_points(
    property_name: str,
    points_df: pd.DataFrame,
    depth: str = "0-5cm",
    quantile: str = "mean",
    cache_dir: Optional[Union[str, Path]] = None,
    lon_col: str = "lon",
    lat_col: str = "lat",
) -> pd.DataFrame:
    """Return *points_df* with a SoilGrids value column appended.

    Parameters
    ----------
    property_name:
        SoilGrids property identifier.
    points_df:
        DataFrame with longitude / latitude columns.
    depth:
        Depth interval string.
    quantile:
        Quantile or summary statistic.
    cache_dir:
        Optional local directory for caching the downloaded GeoTIFF.
    lon_col, lat_col:
        Names of the longitude / latitude columns.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with column ``"{property_name}_{depth}_{quantile}"``
        appended.
    """
    url = resolve_soilgrids_url(property_name, depth, quantile)
    source = _resolve_source(url, cache_dir)
    col_name = f"{property_name}_{depth}_{quantile}"
    result = sample_geotiff_at_points(
        source, points_df, lon_col=lon_col, lat_col=lat_col, band_prefix="band_"
    )
    result = result.rename(columns={"band_1": col_name})
    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_source(url: str, cache_dir: Optional[Union[str, Path]]) -> Union[str, Path]:
    """Return local path if *cache_dir* given, else return the remote URL."""
    if cache_dir is None:
        return url
    dest = Path(cache_dir) / Path(url).name
    if not dest.exists():
        download_file(url, dest)
    return dest
