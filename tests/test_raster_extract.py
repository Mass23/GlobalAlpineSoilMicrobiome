"""Tests for raster-extraction helpers using synthetic GeoTIFFs."""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_bounds

from globalalpine.raster_extract import sample_geotiff_at_point, sample_geotiff_at_points
import pandas as pd


def _make_test_tif(tmp_path: Path, value: float = 42.0) -> Path:
    """Write a tiny 10×10 single-band GeoTIFF centred near (0°, 0°)."""
    path = tmp_path / "test.tif"
    transform = from_bounds(-1, -1, 1, 1, 10, 10)
    data = np.full((1, 10, 10), value, dtype=np.float32)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=10,
        width=10,
        count=1,
        dtype=np.float32,
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999.0,
    ) as dst:
        dst.write(data)
    return path


def test_sample_at_point_returns_value(tmp_path):
    tif = _make_test_tif(tmp_path, value=42.0)
    result = sample_geotiff_at_point(tif, lon=0.0, lat=0.0)
    assert 1 in result
    assert math.isclose(result[1], 42.0, rel_tol=1e-5)


def test_sample_at_point_nodata_returns_nan(tmp_path):
    tif = _make_test_tif(tmp_path, value=-9999.0)
    result = sample_geotiff_at_point(tif, lon=0.0, lat=0.0)
    assert math.isnan(result[1])


def test_sample_at_point_invalid_coords(tmp_path):
    tif = _make_test_tif(tmp_path)
    with pytest.raises(ValueError):
        sample_geotiff_at_point(tif, lon=200.0, lat=0.0)


def test_sample_at_points_returns_dataframe(tmp_path):
    tif = _make_test_tif(tmp_path, value=7.5)
    df = pd.DataFrame({"lon": [0.0, 0.1], "lat": [0.0, 0.1]})
    result = sample_geotiff_at_points(tif, df)
    assert "band_1" in result.columns
    assert len(result) == 2
    assert math.isclose(result["band_1"].iloc[0], 7.5, rel_tol=1e-5)


def test_sample_at_points_missing_column(tmp_path):
    tif = _make_test_tif(tmp_path)
    df = pd.DataFrame({"x": [0.0], "y": [0.0]})
    with pytest.raises(KeyError):
        sample_geotiff_at_points(tif, df)
