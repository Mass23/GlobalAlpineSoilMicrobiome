"""Terrain metrics derived from a Copernicus DEM GeoTIFF.

Functions
---------
get_terrain_from_dem(dem_path, lon, lat)
    Derive elevation, slope (°), and aspect (°) at a point from a local DEM.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Union

import numpy as np
import rasterio
from rasterio.windows import Window

from .coordinates import validate_point


def _slope_aspect_from_window(window: np.ndarray, cellsize_m: float) -> tuple[float, float]:
    """Compute slope (°) and aspect (°) from a 3×3 DEM window.

    Uses the standard Horn (1981) finite-difference algorithm.

    Parameters
    ----------
    window:
        3×3 numpy array of elevation values in metres.
    cellsize_m:
        Approximate cell size in metres.

    Returns
    -------
    tuple[float, float]
        ``(slope_deg, aspect_deg)`` where aspect is measured clockwise from
        north (0–360°).
    """
    if window.shape != (3, 3):
        raise ValueError("window must be 3×3")

    # East–west gradient
    dzdx = (
        (window[0, 2] + 2 * window[1, 2] + window[2, 2])
        - (window[0, 0] + 2 * window[1, 0] + window[2, 0])
    ) / (8.0 * cellsize_m)

    # North–south gradient (positive northward)
    dzdy = (
        (window[0, 0] + 2 * window[0, 1] + window[0, 2])
        - (window[2, 0] + 2 * window[2, 1] + window[2, 2])
    ) / (8.0 * cellsize_m)

    slope_deg = math.degrees(math.atan(math.sqrt(dzdx**2 + dzdy**2)))

    # Aspect: clockwise from north
    aspect_rad = math.atan2(dzdx, dzdy)
    aspect_deg = math.degrees(aspect_rad)
    if aspect_deg < 0:
        aspect_deg += 360.0

    return slope_deg, aspect_deg


def get_terrain_from_dem(
    dem_path: Union[str, Path],
    lon: float,
    lat: float,
) -> dict[str, float]:
    """Derive elevation, slope, and aspect at a point from a DEM GeoTIFF.

    Reads a 3×3 neighbourhood around the target pixel to compute slope and
    aspect using the Horn (1981) algorithm.

    Parameters
    ----------
    dem_path:
        Local path to a GeoTIFF DEM (e.g. a Copernicus GLO-30 tile).
    lon, lat:
        WGS-84 longitude / latitude in decimal degrees.

    Returns
    -------
    dict
        ``{"elevation_m": float, "slope_deg": float, "aspect_deg": float}``

    Raises
    ------
    ValueError
        If the coordinates are out of valid range or the pixel is too close to
        the raster edge to extract a 3×3 window.
    """
    validate_point(lon, lat)

    with rasterio.open(str(dem_path)) as src:
        row_f, col_f = src.index(lon, lat)
        row, col = int(row_f), int(col_f)

        if row < 1 or col < 1 or row >= src.height - 1 or col >= src.width - 1:
            raise ValueError(
                f"Point ({lon}, {lat}) is too close to the raster edge to compute terrain metrics."
            )

        # Read a 3×3 window around the central pixel
        window = Window(col - 1, row - 1, 3, 3)
        data = src.read(1, window=window).astype(float)

        nodata = src.nodata
        if nodata is not None:
            data[np.isclose(data, float(nodata))] = float("nan")

        # Approximate cell size in metres from the transform
        res = src.res  # (row_res, col_res) in CRS units
        # For geographic CRS, convert degrees to metres at this latitude
        try:
            crs = src.crs
            is_geographic = crs.is_geographic if crs is not None else True
        except Exception:
            is_geographic = True

        if is_geographic:
            # 1° ≈ 111_320 m; use the smaller of the two resolutions
            cellsize_m = float(min(res)) * 111_320.0
        else:
            cellsize_m = float(min(res))

    elevation = float(data[1, 1])
    slope_deg, aspect_deg = _slope_aspect_from_window(data, cellsize_m)

    return {
        "elevation_m": elevation,
        "slope_deg": slope_deg,
        "aspect_deg": aspect_deg,
    }
