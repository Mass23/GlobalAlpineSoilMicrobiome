"""Tests for coordinate validation helpers."""

import pytest

from globalalpine.coordinates import validate_point, validate_bbox


def test_valid_point():
    validate_point(9.8, 46.8)


def test_valid_point_boundary():
    validate_point(-180.0, -90.0)
    validate_point(180.0, 90.0)


def test_invalid_longitude_high():
    with pytest.raises(ValueError, match="longitude"):
        validate_point(181.0, 0.0)


def test_invalid_longitude_low():
    with pytest.raises(ValueError, match="longitude"):
        validate_point(-181.0, 0.0)


def test_invalid_latitude_high():
    with pytest.raises(ValueError, match="latitude"):
        validate_point(0.0, 91.0)


def test_invalid_latitude_low():
    with pytest.raises(ValueError, match="latitude"):
        validate_point(0.0, -91.0)


def test_valid_bbox():
    validate_bbox(9.0, 46.0, 10.0, 47.0)


def test_bbox_min_lon_equals_max_lon():
    with pytest.raises(ValueError, match="min_lon"):
        validate_bbox(9.0, 46.0, 9.0, 47.0)


def test_bbox_min_lat_equals_max_lat():
    with pytest.raises(ValueError, match="min_lat"):
        validate_bbox(9.0, 46.0, 10.0, 46.0)


def test_bbox_invalid_coordinate():
    with pytest.raises(ValueError):
        validate_bbox(-200.0, 46.0, 10.0, 47.0)
