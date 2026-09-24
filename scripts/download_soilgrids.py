#!/usr/bin/env python3
"""
SoilGrids 2.0 soil properties at sampling points.

All paths below are under PROJECT_DIR (/home/renku/work/GlobalAlpineSoilMicrobiome).
Reads   data/points.csv              columns lat, lon, depth_cm (may be empty); others kept
Writes  data/soil/soilgrids_soil.csv points + two sets of columns per property:
          <prop>_30cm_avg   thickness-weighted mean of 0-5, 5-15, 15-30 cm (5:10:15)
          <prop>_depth_val  value of the SoilGrids layer containing depth_cm
                            (layer top <= depth < layer bottom, 200 cm included);
                            empty when depth_cm is empty or outside 0-200 cm
        data/tmp/soilgrids_cache.jsonl  fetched pixels; lets an interrupted run resume

Properties (mean prediction, converted to conventional units as in ISRIC's table):
  bdod kg/dm3, cec cmol(c)/kg, cfvo vol%, clay %, sand %, silt %, nitrogen g/kg, phh2o pH,
  soc g/kg, ocd kg/m3, wv0010 / wv0033 / wv1500 vol% (water content at 10 / 33 / 1500 kPa)
  ocs kg/m2: organic carbon STOCK, only published for 0-30 cm -> ocs_30cm_avg only
SoilGrids layers: 0-5, 5-15, 15-30, 30-60, 60-100, 100-200 cm, 250 m resolution.

Data access goes through the `soilgrids` package's WCS services (maps.isric.org). For each
unique coordinate a 3x3-pixel tile (~250 m pixels) centred on the point is requested in
EPSG:4326 and the centre pixel is read. The WCS connection is opened once per property;
the package's get_coverage_data() would re-download the service description on every call
and open results with rioxarray/pyproj, which needs PROJ's database.
"""
# no `from __future__` import: keeps the file usable with Snakemake's `script:` preamble too
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.io import MemoryFile
from soilgrids import SoilGrids

# --------------------------------------------------------------------------- settings
PROJECT_DIR = Path("/home/renku/work/GlobalAlpineSoilMicrobiome")
POINTS_CSV = PROJECT_DIR / "data/points.csv"
OUT_CSV = PROJECT_DIR / "data/soil/soilgrids_soil.csv"
TMP_DIR = PROJECT_DIR / "data/tmp"
CACHE_JSONL = TMP_DIR / "soilgrids_cache.jsonl"

# property -> factor dividing the stored integers (ISRIC conversion table)
PROPERTIES = {
    "bdod": 100, "cec": 10, "cfvo": 10, "clay": 10, "sand": 10, "silt": 10,
    "nitrogen": 100, "phh2o": 10, "soc": 10, "ocd": 10,
    "wv0010": 10, "wv0033": 10, "wv1500": 10,
}
OCS_FACTOR = 10                                     # ocs: t/ha -> kg/m2, 0-30 cm only
LAYERS = [(0, 5), (5, 15), (15, 30), (30, 60), (60, 100), (100, 200)]
TOP30 = LAYERS[:3]
STAT = "mean"

HALF_WIDTH_DEG = 0.003375                           # 3 x 3 tile of ~0.00225 deg (~250 m) pixels
N_WORKERS = 4                                       # parallel requests; keep modest for ISRIC
TRIES = 3

# the package's catalogue lacks the water-retention layers; same WCS pattern at ISRIC
for wv, kpa in [("wv0010", 10), ("wv0033", 33), ("wv1500", 1500)]:
    SoilGrids.MAP_SERVICES.setdefault(wv, {
        "name": f"Volumetric water content at {kpa} kPa",
        "link": f"https://maps.isric.org/mapserv?map=/map/{wv}.map",
        "units": "10^-3 cm3/cm3",
    })

log = logging.getLogger("soilgrids")


def coverage_id(prop: str, top: int, bottom: int) -> str:
    return f"{prop}_{top}-{bottom}cm_{STAT}"


def coord_key(lon: float, lat: float) -> str:
    return f"{lon:.6f},{lat:.6f}"


# --------------------------------------------------------------------------- WCS access
class Services:
    """One WCS connection per property, opened lazily through the soilgrids package."""

    def __init__(self):
        self._wcs, self._crs = {}, {}

    def get(self, prop):
        if prop not in self._wcs:
            wcs, coverages = SoilGrids._get_service_and_coverage_list(prop)
            crs = None
            for cov in coverages:                   # EPSG:4326 code as the service spells it
                codes = [c.getcodeurn() for c in wcs.contents[cov].supportedCRS]
                crs = next((c for c in codes if c.endswith(":4326")), None)
                break
            if crs is None:
                raise RuntimeError(f"{prop}: WCS offers no EPSG:4326 output")
            self._wcs[prop], self._crs[prop] = (wcs, set(coverages)), crs
        return self._wcs[prop], self._crs[prop]


def read_centre(tif_bytes: bytes):
    """centre pixel of a small GeoTIFF tile; None for nodata"""
    with MemoryFile(tif_bytes) as mem, mem.open() as src:
        arr = src.read(1, masked=True)
        v = arr[arr.shape[0] // 2, arr.shape[1] // 2]
        if np.ma.is_masked(v) or int(v) == -32768:
            return None
        return int(v)


def fetch(services: Services, prop: str, cov: str, lon: float, lat: float):
    (wcs, _), crs = services.get(prop)
    bbox = (lon - HALF_WIDTH_DEG, lat - HALF_WIDTH_DEG, lon + HALF_WIDTH_DEG, lat + HALF_WIDTH_DEG)
    for attempt in range(1, TRIES + 1):
        try:
            resp = wcs.getCoverage(identifier=cov, crs=crs, bbox=bbox, width=3, height=3,
                                   response_crs=crs, format="GEOTIFF_INT16", timeout=120)
            body = resp.read()
            if "tiff" not in (resp.info().get("Content-Type", "") or "").lower():
                raise RuntimeError(body[:300].decode("utf-8", "replace"))
            return read_centre(body)
        except Exception:
            if attempt == TRIES:
                raise
            time.sleep(10 * attempt)


# --------------------------------------------------------------------------- cache
class Cache:
    def __init__(self, path: Path):
        self.path, self.data = path, {}
        if path.exists():
            with open(path) as fh:
                for line in fh:
                    try:
                        r = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    self.data[(r["cov"], r["ck"])] = r["raw"]

    def add(self, cov, ck, raw):
        self.data[(cov, ck)] = raw
        with open(self.path, "a") as fh:
            fh.write(json.dumps({"cov": cov, "ck": ck, "raw": raw}) + "\n")


# --------------------------------------------------------------------------- points
def read_points(path: Path):
    df = pd.read_csv(path)
    missing = {"lat", "lon", "depth_cm"} - set(df.columns)
    if missing:
        raise ValueError(f"{path} lacks column(s) {sorted(missing)}; required: lat, lon, depth_cm")
    lat = pd.to_numeric(df["lat"], errors="coerce")
    lon = pd.to_numeric(df["lon"], errors="coerce")
    bad = lat.isna() | lon.isna() | ~lat.between(-90, 90) | ~lon.between(-180, 180)
    if bad.any():
        raise ValueError(f"invalid lat/lon in rows {df.index[bad].tolist()[:10]} (decimal degrees, WGS84)")
    depth = pd.to_numeric(df["depth_cm"], errors="coerce")
    return df, lat.to_numpy(float), lon.to_numpy(float), depth.to_numpy(float)


def layer_for(depth: float):
    if not np.isfinite(depth) or depth < 0 or depth > LAYERS[-1][1]:
        return None
    for top, bottom in LAYERS:
        if top <= depth < bottom:
            return (top, bottom)
    return LAYERS[-1]                               # depth == 200


# --------------------------------------------------------------------------- main
def main() -> int:
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, force=True,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("rasterio").setLevel(logging.ERROR)   # PROJ chatter; values need no CRS
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    df, lat, lon, depth = read_points(POINTS_CSV)
    coords = sorted({coord_key(lon[i], lat[i]) for i in range(len(df))})
    log.info("%d points, %d unique coordinates, %d with a depth", len(df), len(coords),
             int(np.isfinite(depth).sum()))

    jobs = [(p, coverage_id(p, t, b)) for p in PROPERTIES for t, b in LAYERS]
    jobs.append(("ocs", "ocs_0-30cm_mean"))
    cache = Cache(CACHE_JSONL)
    todo = [(p, cov, ck) for p, cov in jobs for ck in coords if (cov, ck) not in cache.data]
    log.info("%d requests needed (%d cached)", len(todo), len(jobs) * len(coords) - len(todo))

    services, failed, unavailable = Services(), {}, set()
    with rasterio.Env(), ThreadPoolExecutor(N_WORKERS) as pool:
        futures = {}
        for p, cov, ck in todo:
            try:
                (_, covs), _ = services.get(p)          # opens the service once, in main thread
            except Exception as e:
                if p not in unavailable:
                    log.error("service %s unavailable: %s", p, e)
                unavailable.add(p)
                continue
            if cov not in covs:
                unavailable.add(cov)
                continue
            lo, la = map(float, ck.split(","))
            futures[pool.submit(fetch, services, p, cov, lo, la)] = (cov, ck)
        for i, fut in enumerate(as_completed(futures), 1):
            cov, ck = futures[fut]
            try:
                cache.add(cov, ck, fut.result())
            except Exception as e:
                failed[(cov, ck)] = str(e)
                log.error("FAILED %s at %s: %s", cov, ck, str(e)[:200])
            if i % 100 == 0 or i == len(futures):
                log.info("  %d/%d requests done", i, len(futures))
    if unavailable:
        log.warning("not offered by the WCS (left empty): %s", ", ".join(sorted(unavailable)))

    # ---- assemble -----------------------------------------------------------
    def value(prop, top, bottom, ck, factor):
        raw = cache.data.get((coverage_id(prop, top, bottom), ck))
        return np.nan if raw is None else raw / factor

    avg = {f"{p}_30cm_avg": [] for p in PROPERTIES}
    avg["ocs_30cm_avg"] = []
    at_depth = {f"{p}_depth_val": [] for p in PROPERTIES}
    for i in range(len(df)):
        ck = coord_key(lon[i], lat[i])
        layer = layer_for(depth[i])
        for p, f in PROPERTIES.items():
            v = [value(p, t, b, ck, f) for t, b in TOP30]
            w = [b - t for t, b in TOP30]
            avg[f"{p}_30cm_avg"].append(np.nan if np.isnan(v).any() else float(np.dot(v, w) / sum(w)))
            at_depth[f"{p}_depth_val"].append(value(p, *layer, ck, f) if layer else np.nan)
        raw = cache.data.get(("ocs_0-30cm_mean", ck))
        avg["ocs_30cm_avg"].append(np.nan if raw is None else raw / OCS_FACTOR)

    out = pd.concat([df, pd.DataFrame(avg).round(4), pd.DataFrame(at_depth).round(4)], axis=1)
    out.to_csv(OUT_CSV, index=False)
    log.info("wrote %s", OUT_CSV)

    n_bad_depth = int((np.isfinite(depth) & np.array([layer_for(d) is None for d in depth])).sum())
    if n_bad_depth:
        log.warning("%d depth_cm value(s) outside 0-200 cm -> _depth_val left empty", n_bad_depth)
    if failed:
        log.error("%d request(s) failed; rerun to retry them (cached values are reused)", len(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
