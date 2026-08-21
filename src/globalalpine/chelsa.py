"""CHELSA climate raster access helpers.

CHELSA v2.1 distributes Cloud-Optimised GeoTIFFs on the WSL server.

Variable naming follows the CHELSA v2.1 documentation:
https://chelsa-climate.org/

Functions
---------
resolve_chelsa_url(variable, period)
    Return the HTTPS URL for a CHELSA COG layer.

fetch_chelsa_at_point(variable, lon, lat, period, cache_dir)
    Return the raster value at a single point.

fetch_chelsa_at_points(variable, points_df, period, cache_dir)
    Return a DataFrame with the CHELSA band column appended.
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

_BASE = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies"

# Monthly period labels used in the URL path
_PERIODS = {
    "1981-2010": "1981-2010",
}

# Non-exhaustive list of supported bioclimate variables
BIOCLIM_VARS = [f"bio{i}" for i in range(1, 20)]

# Monthly variables: prefix "pr" (precip), "tas" (mean temp), "tasmax", "tasmin"
MONTHLY_VARS = (
    [f"pr_{m:02d}" for m in range(1, 13)]
    + [f"tas_{m:02d}" for m in range(1, 13)]
    + [f"tasmax_{m:02d}" for m in range(1, 13)]
    + [f"tasmin_{m:02d}" for m in range(1, 13)]
)


def resolve_chelsa_url(
    variable: str,
    period: str = "1981-2010",
) -> str:
    """Return the HTTPS URL for a CHELSA v2.1 COG layer.

    Parameters
    ----------
    variable:
        CHELSA variable name, e.g. ``"bio1"`` (mean annual temperature) or
        ``"pr_01"`` (January precipitation).
    period:
        Climatological period string, e.g. ``"1981-2010"``.

    Returns
    -------
    str
        Full HTTPS URL to the Cloud-Optimised GeoTIFF.

    Examples
    --------
    >>> resolve_chelsa_url("bio1")
    'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif'
    """
    if variable.startswith("bio"):
        subfolder = "bio"
        filename = f"CHELSA_{variable}_{period}_V.2.1.tif"
    else:
        # monthly variable like "pr_01"
        parts = variable.split("_")
        subfolder = parts[0]
        filename = f"CHELSA_{variable}_{period}_V.2.1.tif"
    return f"{_BASE}/{period}/{subfolder}/{filename}"


def fetch_chelsa_at_point(
    variable: str,
    lon: float,
    lat: float,
    period: str = "1981-2010",
    cache_dir: Optional[Union[str, Path]] = None,
) -> float:
    """Return the CHELSA value at a single point.

    Parameters
    ----------
    variable:
        CHELSA variable name.
    lon, lat:
        WGS-84 longitude / latitude in decimal degrees.
    period:
        Climatological period string.
    cache_dir:
        Optional local directory for caching the downloaded GeoTIFF.

    Returns
    -------
    float
        Raster value at the point (``float("nan")`` if no-data).
    """
    validate_point(lon, lat)
    url = resolve_chelsa_url(variable, period)
    source = _resolve_source(url, cache_dir)
    band_values = sample_geotiff_at_point(source, lon, lat)
    return band_values.get(1, float("nan"))


def fetch_chelsa_at_points(
    variable: str,
    points_df: pd.DataFrame,
    period: str = "1981-2010",
    cache_dir: Optional[Union[str, Path]] = None,
    lon_col: str = "lon",
    lat_col: str = "lat",
) -> pd.DataFrame:
    """Return *points_df* with a CHELSA value column appended.

    Parameters
    ----------
    variable:
        CHELSA variable name.
    points_df:
        DataFrame with longitude / latitude columns.
    period:
        Climatological period string.
    cache_dir:
        Optional local directory for caching.
    lon_col, lat_col:
        Names of the longitude / latitude columns.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with column ``"{variable}_{period}"`` appended.
    """
    url = resolve_chelsa_url(variable, period)
    source = _resolve_source(url, cache_dir)
    col_name = f"{variable}_{period}"
    result = sample_geotiff_at_points(
        source, points_df, lon_col=lon_col, lat_col=lat_col, band_prefix="band_"
    )
    result = result.rename(columns={"band_1": col_name})
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
