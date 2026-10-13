#!/usr/bin/env python3
"""
SoilGrids 2.0 soil properties at sampling points.

Run     python3 download_soilgrids.py [points.csv soilgrids_data.csv [properties]]
        Snakemake passes the two paths; the tmp/ cache is written next to the output. `properties` is an
        optional comma-separated list (e.g. phh2o,soc,nitrogen,clay,sand,silt,cec,bdod); empty = all.
        SOILGRIDS_WORKERS (environment variable) sets the number of parallel requests (default 4).

Speed: one WCS request per coordinate x property x layer is the whole cost (about 2.6 requests/s with 4
workers). Only the layers the output needs are requested: the top three (0-5, 5-15, 15-30 cm) for every
coordinate, plus a deeper layer only where a sample's depth_cm falls in it. Choosing fewer properties is
the other big lever (each property = 3 requests per coordinate); the log shows the speed and an ETA.
Reads   points.csv                   columns lat, lon, depth_cm (may be empty); others kept
Writes  soilgrids_data.csv           points + two sets of columns per property:
          <prop>_30cm_avg   thickness-weighted mean of 0-5, 5-15, 15-30 cm (5:10:15)
          <prop>_depth_val  value of the SoilGrids layer containing depth_cm
                            (layer top <= depth < layer bottom, 200 cm included);
                            empty when depth_cm is empty or outside 0-200 cm
        tmp/soilgrids_cache.jsonl       fetched pixels; lets an interrupted run resume

Failures are of two kinds, handled differently:
  * the server ANSWERED with an OGC service exception (e.g. a point outside the map: SoilGrids has no soil
    on the ice sheets): asking again gives the same answer, so it is NOT retried; the value is NA, the
    exception text is kept in the cache, and it does not make the run fail unless more than
    MAX_REFUSED_FRACTION of the coordinates are refused (that would point at a bug, not at the map edge);
  * network / server trouble (timeouts, 5xx, a non-XML page): retried TRIES times with a pause, then
    reported as FAILED (exit code 1; a rerun retries only those).

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
import os
import re
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.io import MemoryFile
from soilgrids import SoilGrids

# --------------------------------------------------------------------------- settings
# Snakemake: download_soilgrids.py <points.csv> <soilgrids_data.csv>; the defaults are for running by hand
PROJECT_DIR = Path("/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome")
POINTS_CSV = Path(sys.argv[1]) if len(sys.argv) >= 3 else PROJECT_DIR / "data/points.csv"
OUT_CSV = Path(sys.argv[2]) if len(sys.argv) >= 3 else PROJECT_DIR / "data/soilgrids_data.csv"
TMP_DIR = OUT_CSV.parent / "tmp"
CACHE_JSONL = TMP_DIR / "soilgrids_cache.jsonl"

# property -> factor dividing the stored integers (ISRIC conversion table)
PROPERTIES = {
    "bdod": 100, "cec": 10, "cfvo": 10, "clay": 10, "sand": 10, "silt": 10,
    "nitrogen": 100, "phh2o": 10, "soc": 10, "ocd": 10,
    "wv0010": 10, "wv0033": 10, "wv1500": 10,
}
if len(sys.argv) >= 4 and sys.argv[3].strip():        # optional 3rd argument: which properties to fetch
    _wanted = [p.strip() for p in sys.argv[3].split(",") if p.strip()]
    _unknown = [p for p in _wanted if p not in PROPERTIES]
    if _unknown:
        raise SystemExit(f"unknown SoilGrids property/properties {_unknown}; known: {sorted(PROPERTIES)}")
    PROPERTIES = {p: PROPERTIES[p] for p in _wanted}
OCS_FACTOR = 10                                     # ocs: t/ha -> kg/m2, 0-30 cm only
LAYERS = [(0, 5), (5, 15), (15, 30), (30, 60), (60, 100), (100, 200)]
TOP30 = LAYERS[:3]
STAT = "mean"

HALF_WIDTH_DEG = 0.003375                           # 3 x 3 tile of ~0.00225 deg (~250 m) pixels
N_WORKERS = int(os.environ.get("SOILGRIDS_WORKERS", 4))   # parallel requests; keep modest for ISRIC
TRIES = 3                                           # for network / server trouble only (see docstring)
MAX_REFUSED_FRACTION = 0.05                         # more refused coordinates than this = suspect, exit 1
RETRY_REFUSED = False                               # True: ask again for coordinates refused in an earlier run

# the package's catalogue lacks the water-retention layers; same WCS pattern at ISRIC
for wv, kpa in [("wv0010", 10), ("wv0033", 33), ("wv1500", 1500)]:
    SoilGrids.MAP_SERVICES.setdefault(wv, {
        "name": f"Volumetric water content at {kpa} kPa",
        "link": f"https://maps.isric.org/mapserv?map=/map/{wv}.map",
        "units": "10^-3 cm3/cm3",
    })

log = logging.getLogger("soilgrids")


class ServiceRefused(Exception):
    """the WCS answered with an OGC service exception: the same request will get the same answer"""


def service_exception(body: bytes):
    """text of an OGC ServiceException / ExceptionReport, or None if the body is not one"""
    text = body[:20000].decode("utf-8", "replace")
    if "ExceptionReport" not in text:
        return None
    parts = re.findall(r"<(?:ows:)?(?:ServiceException|ExceptionText)\b([^>]*)>(.*?)</(?:ows:)?(?:ServiceException|ExceptionText)>",
                       text, flags=re.S)
    msg = " | ".join(" ".join((f"[{re.sub(chr(34), '', a).strip()}]" if a.strip() else "", m.strip())).strip()
                     for a, m in parts)
    return re.sub(r"\s+", " ", msg or text)[:500]


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
                refusal = service_exception(body)
                if refusal is not None:
                    raise ServiceRefused(refusal)
                raise RuntimeError(body[:300].decode("utf-8", "replace"))
            return read_centre(body)
        except ServiceRefused:
            raise                                   # the server answered: asking again changes nothing
        except Exception:
            if attempt == TRIES:
                raise
            time.sleep(10 * attempt)


# --------------------------------------------------------------------------- cache
class Cache:
    def __init__(self, path: Path):
        self.path, self.data, self.refused = path, {}, {}
        if path.exists():
            with open(path) as fh:
                for line in fh:
                    try:
                        r = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    key = (r["cov"], r["ck"])
                    if r.get("refused") is not None:        # a recorded refusal (value NA)
                        self.refused[key] = r["refused"]
                        if RETRY_REFUSED:
                            continue                        # not cached: it is asked again
                    self.data[key] = r["raw"]

    def add(self, cov, ck, raw):
        self.data[(cov, ck)] = raw
        with open(self.path, "a") as fh:
            fh.write(json.dumps({"cov": cov, "ck": ck, "raw": raw}) + "\n")

    def add_refused(self, cov, ck, message):
        self.data[(cov, ck)] = None
        self.refused[(cov, ck)] = message
        with open(self.path, "a") as fh:
            fh.write(json.dumps({"cov": cov, "ck": ck, "raw": None, "refused": message}) + "\n")


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

    # layers needed per coordinate: the top 30 cm (for <prop>_30cm_avg) and, only where a sample has a depth_cm,
    # the layer containing it (for <prop>_depth_val). The other layers would never reach the output.
    layers_needed = {}
    for i in range(len(df)):
        ck = coord_key(lon[i], lat[i])
        layers = layers_needed.setdefault(ck, set(TOP30))
        layer = layer_for(depth[i])
        if layer:
            layers.add(layer)
    cache = Cache(CACHE_JSONL)
    todo, needed = [], 0
    for ck, layers in layers_needed.items():
        wanted = [(p, coverage_id(p, t, b)) for p in PROPERTIES for t, b in sorted(layers)]
        wanted.append(("ocs", "ocs_0-30cm_mean"))
        needed += len(wanted)
        todo += [(p, cov, ck) for p, cov in wanted if (cov, ck) not in cache.data]
    log.info("%d requests needed (%d cached) for %d properties: %s", len(todo), needed - len(todo),
             len(PROPERTIES), ",".join(PROPERTIES))

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
        t0 = time.monotonic()
        for i, fut in enumerate(as_completed(futures), 1):
            cov, ck = futures[fut]
            try:
                cache.add(cov, ck, fut.result())
            except ServiceRefused as e:
                cache.add_refused(cov, ck, str(e))          # NA, not retried, not a failure by itself
            except Exception as e:
                failed[(cov, ck)] = str(e)
                log.error("FAILED %s at %s: %s", cov, ck, str(e)[:200])
            if i % 100 == 0 or i == len(futures):
                rate = i / max(time.monotonic() - t0, 1e-9)
                log.info("  %d/%d requests done (%.1f/s, about %.1f h left)", i, len(futures), rate,
                         (len(futures) - i) / rate / 3600)
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
    out.to_csv(OUT_CSV, index=False, na_rep="NA")
    log.info("wrote %s", OUT_CSV)

    too_many_refused = False
    if cache.refused:
        refused_coords = sorted({ck for (_, ck) in cache.refused})
        lats = [float(ck.split(",")[1]) for ck in refused_coords]
        log.warning("%d of %d coordinate(s) were refused by the WCS for at least one layer (NA in the output); "
                    "latitudes %.1f to %.1f", len(refused_coords), len(coords), min(lats), max(lats))
        for msg, n in Counter(cache.refused.values()).most_common(3):
            log.warning("  %d x  %s", n, msg)
        if len(refused_coords) > MAX_REFUSED_FRACTION * len(coords):
            log.error("more than %.0f%% of the coordinates were refused: that points at a bug or a service "
                      "problem, not at the edge of the map", 100 * MAX_REFUSED_FRACTION)
            too_many_refused = True

    n_bad_depth = int((np.isfinite(depth) & np.array([layer_for(d) is None for d in depth])).sum())
    if n_bad_depth:
        log.warning("%d depth_cm value(s) outside 0-200 cm -> _depth_val left empty", n_bad_depth)
    if failed:
        log.error("%d request(s) failed; rerun to retry them (cached values are reused)", len(failed))
        return 1
    return 1 if too_many_refused else 0


if __name__ == "__main__":
    sys.exit(main())