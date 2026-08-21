"""Coordinate validation helpers."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GPSPoint:
    """A single geographic point with longitude and latitude (WGS-84)."""

    lon: float
    lat: float


@dataclass(frozen=True)
class BoundingBox:
    """An axis-aligned bounding box in WGS-84 degrees."""

    min_lon: float
    min_lat: float
    max_lon: float
    max_lat: float


def validate_point(lon: float, lat: float) -> None:
    """Raise *ValueError* if *lon* or *lat* are out of the WGS-84 range."""
    if not (-180.0 <= lon <= 180.0):
        raise ValueError(f"longitude must be in [-180, 180], got {lon}")
    if not (-90.0 <= lat <= 90.0):
        raise ValueError(f"latitude must be in [-90, 90], got {lat}")


def validate_bbox(min_lon: float, min_lat: float, max_lon: float, max_lat: float) -> None:
    """Raise *ValueError* if the bounding box is invalid."""
    validate_point(min_lon, min_lat)
    validate_point(max_lon, max_lat)
    if min_lon >= max_lon:
        raise ValueError(f"min_lon ({min_lon}) must be less than max_lon ({max_lon})")
    if min_lat >= max_lat:
        raise ValueError(f"min_lat ({min_lat}) must be less than max_lat ({max_lat})")
