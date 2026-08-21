"""Tests for terrain derivation helpers."""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_bounds

from globalalpine.terrain import _slope_aspect_from_window, get_terrain_from_dem


def _make_flat_dem(tmp_path: Path, elevation: float = 1000.0) -> Path:
    """Write a uniform-elevation DEM tile (5×5) for testing."""
    path = tmp_path / "dem.tif"
    transform = from_bounds(-1, -1, 1, 1, 5, 5)
    data = np.full((1, 5, 5), elevation, dtype=np.float32)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=5,
        width=5,
        count=1,
        dtype=np.float32,
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999.0,
    ) as dst:
        dst.write(data)
    return path


def test_slope_aspect_flat_terrain():
    window = np.full((3, 3), 1000.0)
    slope, aspect = _slope_aspect_from_window(window, cellsize_m=30.0)
    assert math.isclose(slope, 0.0, abs_tol=1e-9)


def test_slope_aspect_wrong_shape():
    with pytest.raises(ValueError, match="3×3"):
        _slope_aspect_from_window(np.zeros((4, 4)), cellsize_m=30.0)


def test_get_terrain_from_dem_flat(tmp_path):
    dem = _make_flat_dem(tmp_path, elevation=500.0)
    result = get_terrain_from_dem(dem, lon=0.0, lat=0.0)
    assert "elevation_m" in result
    assert "slope_deg" in result
    assert "aspect_deg" in result
    assert math.isclose(result["elevation_m"], 500.0, rel_tol=1e-5)
    assert math.isclose(result["slope_deg"], 0.0, abs_tol=1e-6)


def test_get_terrain_from_dem_edge_raises(tmp_path):
    """Point at the raster edge should raise ValueError."""
    dem = _make_flat_dem(tmp_path)
    # Use a coordinate outside or at the boundary
    with pytest.raises(ValueError):
        get_terrain_from_dem(dem, lon=-0.99, lat=-0.99)
