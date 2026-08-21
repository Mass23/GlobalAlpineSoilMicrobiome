"""Tests for URL-resolution helpers (no network access)."""

from globalalpine.soilgrids import resolve_soilgrids_url
from globalalpine.chelsa import resolve_chelsa_url
from globalalpine.copernicus import resolve_copernicus_dem_url


def test_soilgrids_url():
    url = resolve_soilgrids_url("phh2o", depth="0-5cm", quantile="mean")
    assert "phh2o" in url
    assert "0-5cm" in url
    assert "mean" in url
    assert url.startswith("https://")


def test_chelsa_url_bio():
    url = resolve_chelsa_url("bio1", period="1981-2010")
    assert "bio1" in url
    assert "1981-2010" in url
    assert url.startswith("https://")


def test_chelsa_url_monthly():
    url = resolve_chelsa_url("pr_01", period="1981-2010")
    assert "pr_01" in url
    assert url.startswith("https://")


def test_copernicus_url_north_east():
    url = resolve_copernicus_dem_url(lon=9.8, lat=46.8, resolution=30)
    assert "N46" in url
    assert "E009" in url
    assert url.startswith("https://")
    assert url.endswith(".tif")


def test_copernicus_url_south_west():
    url = resolve_copernicus_dem_url(lon=-73.5, lat=-15.2, resolution=30)
    assert "S16" in url
    assert "W074" in url


def test_copernicus_url_invalid_resolution():
    import pytest
    with pytest.raises(ValueError, match="resolution"):
        resolve_copernicus_dem_url(9.8, 46.8, resolution=60)
