"""Simple CLI entry point."""

from __future__ import annotations

import argparse
import sys

from .coordinates import validate_point


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="globalalpine",
        description="GlobalAlpineSoilMicrobiome – validate a GPS coordinate.",
    )
    parser.add_argument("--lon", type=float, required=True, help="Longitude (WGS-84)")
    parser.add_argument("--lat", type=float, required=True, help="Latitude (WGS-84)")
    args = parser.parse_args(argv)

    try:
        validate_point(args.lon, args.lat)
        print(f"Valid point: lon={args.lon}, lat={args.lat}")
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
