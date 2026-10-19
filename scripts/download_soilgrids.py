#!/usr/bin/env python3
"""
SoilGrids 2.0 soil properties at sampling points, read straight from ISRIC's global map files.

Run     python3 download_soilgrids.py [points.csv soilgrids_data.csv [properties [lowland_sets.csv]]]
        Snakemake passes the paths; the tmp/ cache is written next to the output. `properties` is an
        optional comma-separated list (e.g. phh2o,soc,nitrogen,clay,sand,silt,cec,bdod); empty = all.
        `lowland_sets.csv` (from define_lowland_sets.py) restricts the download, see Subset.
        SOILGRIDS_WORKERS (environment variable) = number of layers read in parallel (default 8).

Reads   points.csv                   columns lat, lon, depth_cm (may be empty); others kept
Writes  soilgrids_data.csv           points + two sets of columns per property:
          <prop>_30cm_avg   thickness-weighted mean of 0-5, 5-15, 15-30 cm (5:10:15)
          <prop>_depth_val  value of the SoilGrids layer containing depth_cm
                            (layer top <= depth < layer bottom, 200 cm included);
                            empty when depth_cm is empty or outside 0-200 cm
        tmp/soilgrids_tiles_cache.jsonl   values read so far; lets an interrupted run resume

How the data are read (the same idea as the CHELSA script): ISRIC publishes every layer as a global VRT in the
Interrupted Goode Homolosine projection (https://files.isric.org/soilgrids/latest/data/<prop>/<prop>_<depth>_mean.vrt),
which GDAL reads remotely and only fetches the parts that are asked for. Each layer is opened ONCE and the pixels of
all the coordinates that need it are read: the coordinates are converted to Homolosine, grouped by 512 x 512-pixel
block, and each group is read with one window read. No per-point request, no WCS, no `soilgrids` package.
The value is the native 250 m pixel that contains the point (no resampling). Points outside the map (the
ice sheets, open water) or on a no-data pixel get NA; this is counted and logged, and never retried.

Subset: with a lowland_sets.csv (one row per sample, column subset_any), only the points of the samples marked
"yes" (every alpine sample + the lowland samples of any of the lowland replicates) are read; the other rows stay NA.
Without that file every row is read (ONLY_SUBSET = False also reads everything).

Only the layers the output needs are read: the top three (0-5, 5-15, 15-30 cm) for every coordinate, plus a
deeper layer only where a sample's depth_cm falls in it.

Failures: network trouble is retried (GDAL's own HTTP retries, then TRIES attempts per window). What still fails
is reported (exit code 1) and a rerun reads only what is missing. A layer file that does not exist on the server
(HTTP 404) is skipped with a warning: its columns stay NA and it does not make the run fail.

Properties (mean prediction, converted to conventional units as in ISRIC's table):
  bdod kg/dm3, cec cmol(c)/kg, cfvo vol%, clay %, sand %, silt %, nitrogen g/kg, phh2o pH,
  soc g/kg, ocd kg/m3, wv0010 / wv0033 / wv1500 vol% (water content at 10 / 33 / 1500 kPa)
  ocs kg/m2: organic carbon STOCK, only published for 0-30 cm -> ocs_30cm_avg only
SoilGrids layers: 0-5, 5-15, 15-30, 30-60, 60-100, 100-200 cm, 250 m resolution.
"""
# no `from __future__` import: keeps the file usable with Snakemake's `script:` preamble too
import json
import logging
import os
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import rasterio.warp
from rasterio.windows import Window

# --------------------------------------------------------------------------- settings
# Snakemake: download_soilgrids.py <points.csv> <soilgrids_data.csv>; the defaults are for running by hand
PROJECT_DIR = Path("/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome")
POINTS_CSV = Path(sys.argv[1]) if len(sys.argv) >= 3 else PROJECT_DIR / "data/points.csv"
OUT_CSV = Path(sys.argv[2]) if len(sys.argv) >= 3 else PROJECT_DIR / "data/soilgrids_data.csv"
TMP_DIR = OUT_CSV.parent / "tmp"
CACHE_JSONL = TMP_DIR / "soilgrids_tiles_cache.jsonl"

# property -> factor dividing the stored integers (ISRIC conversion table)
PROPERTIES = {
    "bdod": 100, "cec": 10, "cfvo": 10, "clay": 10, "sand": 10, "silt": 10,
    "nitrogen": 100, "phh2o": 10, "soc": 10, "ocd": 10,
    "wv0010": 10, "wv0033": 10, "wv1500": 10,
}
if len(sys.argv) >= 4 and sys.argv[3].strip():        # optional 3rd argument: which properties to read
    _wanted = [p.strip() for p in sys.argv[3].split(",") if p.strip()]
    _unknown = [p for p in _wanted if p not in PROPERTIES]
    if _unknown:
        raise SystemExit(f"unknown SoilGrids property/properties {_unknown}; known: {sorted(PROPERTIES)}")
    PROPERTIES = {p: PROPERTIES[p] for p in _wanted}
OCS_FACTOR = 10                                     # ocs: t/ha -> kg/m2, 0-30 cm only
LAYERS = [(0, 5), (5, 15), (15, 30), (30, 60), (60, 100), (100, 200)]
TOP30 = LAYERS[:3]
STAT = "mean"
NODATA = -32768                                     # SoilGrids Int16 no-data

VRT_BASE = "https://files.isric.org/soilgrids/latest/data"
LONLAT = "+proj=longlat +datum=WGS84 +no_defs"
IGH = "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"   # ISRIC's Homolosine; no PROJ database needed

N_WORKERS = int(os.environ.get("SOILGRIDS_WORKERS", 8))   # layers read in parallel
TRIES = 3                                           # attempts per window / per file open, after GDAL's own retries
BLOCK = 512                                         # pixels: coordinates are grouped by block of this size
MAX_PIXELS_PER_POINT = 4096                         # a group window larger than this per point is read point by point
GDAL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",       # do not list the remote directory
    CPL_VSIL_CURL_CACHE_SIZE="200000000",
    GDAL_HTTP_MAX_RETRY="5",
    GDAL_HTTP_RETRY_DELAY="3",
    GDAL_HTTP_TIMEOUT="120",
    GDAL_CACHEMAX=512,                              # MB; rasterio wants an int here
    VSI_CACHE="TRUE",
    VSI_CACHE_SIZE="134217728",
)
SUBSETS_CSV = Path(sys.argv[4]) if len(sys.argv) >= 5 and sys.argv[4].strip() else None   # lowland_sets.csv
SUBSET_COLUMN = "subset_any"                        # column of that file: "yes" = read this sample
ONLY_SUBSET = True                                  # False: read every row even if the file is given
PREFLIGHT_POINT = (8.5, 46.5)                       # lon, lat of a point with soil data (Switzerland)

log = logging.getLogger("soilgrids")


def coverage_id(prop: str, top: int, bottom: int) -> str:
    return f"{prop}_{top}-{bottom}cm_{STAT}"


def coord_key(lon: float, lat: float) -> str:
    return f"{lon:.6f},{lat:.6f}"


# --------------------------------------------------------------------------- reading the layers
def layer_url(cov: str) -> str:
    """GDAL path of a layer's global VRT, e.g. phh2o_0-5cm_mean -> .../phh2o/phh2o_0-5cm_mean.vrt"""
    return f"/vsicurl/{VRT_BASE}/{cov.split('_')[0]}/{cov}.vrt"


def open_layer(url: str):
    return rasterio.open(url)


def http_status(url: str):
    """HTTP status of a HEAD request on a /vsicurl/ path (None if it is not an http(s) address or unreachable)"""
    if not url.startswith("/vsicurl/http"):
        return None
    try:
        req = urllib.request.Request(url[len("/vsicurl/"):], method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status
    except urllib.error.HTTPError as e:
        return e.code
    except Exception:
        return None


class LayerMissing(Exception):
    """the layer file does not exist on the server (HTTP 404)"""


def open_with_retry(url: str):
    last = None
    for attempt in range(1, TRIES + 1):
        try:
            return open_layer(url)
        except Exception as e:                       # rasterio.errors.RasterioIOError and friends
            last = e
            if http_status(url) == 404:
                raise LayerMissing(url) from e
            if attempt < TRIES:
                time.sleep(5 * attempt)
    raise RuntimeError(f"cannot open {url}: {last}") from last


def to_igh(lons, lats):
    xs, ys = rasterio.warp.transform(LONLAT, IGH, [float(v) for v in lons], [float(v) for v in lats])
    return np.asarray(xs, float), np.asarray(ys, float)


def pixel_index(src, xs, ys):
    """(row, col) of the pixel that contains each (x, y); the raster is assumed north-up (VRT of tiles)"""
    t = src.transform
    cols = np.floor((xs - t.c) / t.a)
    rows = np.floor((ys - t.f) / t.e)
    return rows, cols


def read_group(src, rows, cols):
    """raw integers at the given pixels; one window read when the points are close, point by point otherwise"""
    r0, r1, c0, c1 = int(rows.min()), int(rows.max()), int(cols.min()), int(cols.max())
    h, w = r1 - r0 + 1, c1 - c0 + 1
    out = []
    if h * w <= MAX_PIXELS_PER_POINT * len(rows):
        arr = src.read(1, window=Window(c0, r0, w, h))
        vals = arr[(rows - r0).astype(int), (cols - c0).astype(int)]
    else:
        vals = np.array([src.read(1, window=Window(int(c), int(r), 1, 1))[0, 0] for r, c in zip(rows, cols)])
    return [None if int(v) == NODATA else int(v) for v in vals]


def read_layer(cov: str, cks: list, xy: dict, cache, progress):
    """read one layer for all the coordinates in `cks`; returns (n_value, n_nodata, n_outside, failed_cks)"""
    n_value = n_nodata = n_outside = n_logged = 0
    failed = []
    with rasterio.Env(**GDAL_ENV):
        src = open_with_retry(layer_url(cov))
        try:
            xs = np.array([xy[ck][0] for ck in cks]); ys = np.array([xy[ck][1] for ck in cks])
            rows, cols = pixel_index(src, xs, ys)
            inside = (rows >= 0) & (rows < src.height) & (cols >= 0) & (cols < src.width) & np.isfinite(rows) & np.isfinite(cols)
            for i in np.flatnonzero(~inside):                    # outside the map: NA, deterministic
                cache.add(cov, cks[i], None); n_outside += 1
            progress(int((~inside).sum()))
            idx = np.flatnonzero(inside)
            groups = {}
            for i in idx:
                groups.setdefault((int(rows[i]) // BLOCK, int(cols[i]) // BLOCK), []).append(i)
            for members in groups.values():
                members = np.array(members)
                for attempt in range(1, TRIES + 1):
                    try:
                        vals = read_group(src, rows[members], cols[members])
                        break
                    except Exception as e:
                        if attempt == TRIES:
                            failed += [cks[i] for i in members]
                            n_logged += 1
                            if n_logged <= 3:                    # the first few, then only the total at the end
                                log.error("FAILED %s (%d points near %s): %s", cov, len(members), cks[members[0]], str(e)[:200])
                            vals = None
                        else:
                            time.sleep(5 * attempt)
                if vals is not None:
                    for i, v in zip(members, vals):
                        cache.add(cov, cks[i], v)
                        n_value += v is not None; n_nodata += v is None
                    progress(len(members))
        finally:
            src.close()
    if failed:
        log.error("%s: %d of %d values could not be read", cov, len(failed), len(cks))
    return n_value, n_nodata, n_outside, failed


# --------------------------------------------------------------------------- cache
class Cache:
    def __init__(self, path: Path):
        self.path, self.data, self.lock = path, {}, threading.Lock()
        if path.exists():
            with open(path) as fh:
                for line in fh:
                    try:
                        r = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    self.data[(r["cov"], r["ck"])] = r["raw"]

    def add(self, cov, ck, raw):
        with self.lock:
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
    logging.getLogger("rasterio").setLevel(logging.ERROR)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    df, lat, lon, depth = read_points(POINTS_CSV)
    if ONLY_SUBSET and SUBSETS_CSV is not None:          # only the flagged samples are read; the rest stay NA
        if "site" not in df.columns:
            raise ValueError(f"{POINTS_CSV} needs a `site` column to apply {SUBSETS_CSV}")
        sets = pd.read_csv(SUBSETS_CSV)
        yes_sites = set(sets.loc[sets[SUBSET_COLUMN].astype(str).str.strip().str.lower() == "yes", "site"].astype(str))
        use = df["site"].astype(str).isin(yes_sites).to_numpy()
    else:
        use = np.ones(len(df), dtype=bool)
    coords = sorted({coord_key(lon[i], lat[i]) for i in np.flatnonzero(use)})
    log.info("%d points, %d of them read (%s), %d unique coordinates, %d with a depth", len(df), int(use.sum()),
             SUBSET_COLUMN if use.sum() < len(df) else "all rows", len(coords), int(np.isfinite(depth[use]).sum()))

    # layers needed per coordinate: the top 30 cm (for <prop>_30cm_avg) and, only where a sample has a depth_cm,
    # the layer containing it (for <prop>_depth_val). The other layers would never reach the output.
    layers_needed = {}
    for i in np.flatnonzero(use):
        ck = coord_key(lon[i], lat[i])
        layers = layers_needed.setdefault(ck, set(TOP30))
        layer = layer_for(depth[i])
        if layer:
            layers.add(layer)
    cache = Cache(CACHE_JSONL)
    jobs, needed = {}, 0
    for ck, layers in layers_needed.items():
        wanted = [coverage_id(p, t, b) for p in PROPERTIES for t, b in sorted(layers)] + ["ocs_0-30cm_mean"]
        needed += len(wanted)
        for cov in wanted:
            if (cov, ck) not in cache.data:
                jobs.setdefault(cov, []).append(ck)
    total = sum(len(v) for v in jobs.values())
    log.info("%d values to read in %d layers (%d cached) for %d properties: %s", total, len(jobs), needed - total,
             len(PROPERTIES), ",".join(PROPERTIES))

    failed, missing_layers, stats = [], [], {}
    if jobs:
        # coordinates -> Homolosine, once for all layers
        todo_cks = sorted({ck for v in jobs.values() for ck in v})
        lo, la = np.array([list(map(float, ck.split(","))) for ck in todo_cks]).T
        x, y = to_igh(lo, la)
        xy = {ck: (xi, yi) for ck, xi, yi in zip(todo_cks, x, y)}

        # preflight: open one layer and read a point with soil data, so a problem shows up at once and clearly
        first = sorted(jobs)[0]
        try:
            with rasterio.Env(**GDAL_ENV):
                src = open_with_retry(layer_url(first))
                try:
                    log.info("layer %s: %d x %d px, %s, nodata %s, CRS: %s", first, src.width, src.height, src.dtypes[0],
                             src.nodata, (src.crs.to_string()[:70] if src.crs else "none"))
                    px, py = to_igh([PREFLIGHT_POINT[0]], [PREFLIGHT_POINT[1]])
                    r, c = pixel_index(src, px, py)
                    probe = read_group(src, r, c)[0] if (0 <= r[0] < src.height and 0 <= c[0] < src.width) else None
                finally:
                    src.close()
        except LayerMissing:
            raise RuntimeError(f"layer file {first} does not exist on the server: check VRT_BASE ({VRT_BASE})")
        except Exception as e:
            raise RuntimeError(f"SoilGrids preflight failed on {first}: {e}") from e
        if probe is None:
            raise RuntimeError(f"preflight: no soil value at lon {PREFLIGHT_POINT[0]}, lat {PREFLIGHT_POINT[1]} in {first}: "
                               f"the projection or the file layout is not what this script assumes")
        log.info("preflight OK: raw value %d at lon %.1f, lat %.1f", probe, *PREFLIGHT_POINT)

        done, t0, last_log, lock = 0, time.monotonic(), time.monotonic(), threading.Lock()

        def progress(n):
            nonlocal done, last_log
            with lock:
                done += n
                if time.monotonic() - last_log > 30:
                    last_log = time.monotonic()
                    rate = done / max(last_log - t0, 1e-9)
                    log.info("  %d/%d values (%.0f/s, about %.1f h left)", done, total, rate, (total - done) / rate / 3600)

        with ThreadPoolExecutor(N_WORKERS) as pool:
            futures = {pool.submit(read_layer, cov, cks, xy, cache, progress): cov for cov, cks in jobs.items()}
            for fut in as_completed(futures):
                cov = futures[fut]
                try:
                    n_value, n_nodata, n_outside, bad = fut.result()
                    stats[cov] = (n_value, n_nodata, n_outside)
                    failed += [(cov, ck) for ck in bad]
                except LayerMissing:
                    missing_layers.append(cov)
                    log.warning("layer %s does not exist on the server (404): left NA", cov)
                except Exception as e:
                    failed += [(cov, ck) for ck in jobs[cov]]
                    log.error("FAILED layer %s: %s", cov, str(e)[:200])
        log.info("read %d values in %.1f min", done, (time.monotonic() - t0) / 60)
        if stats:
            n_v = sum(s[0] for s in stats.values()); n_n = sum(s[1] for s in stats.values()); n_o = sum(s[2] for s in stats.values())
            log.info("values: %d with data, %d on a no-data pixel (water, ice, bare rock...), %d outside the map", n_v, n_n, n_o)

    # ---- assemble -----------------------------------------------------------
    def value(prop, top, bottom, ck, factor):
        raw = cache.data.get((coverage_id(prop, top, bottom), ck))
        return np.nan if raw is None else raw / factor

    avg = {f"{p}_30cm_avg": [] for p in PROPERTIES}
    avg["ocs_30cm_avg"] = []
    at_depth = {f"{p}_depth_val": [] for p in PROPERTIES}
    for i in range(len(df)):
        if not use[i]:                                  # not in the subset: NA everywhere
            for p in PROPERTIES:
                avg[f"{p}_30cm_avg"].append(np.nan)
                at_depth[f"{p}_depth_val"].append(np.nan)
            avg["ocs_30cm_avg"].append(np.nan)
            continue
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

    first_prop = next(iter(PROPERTIES))
    ref = f"{first_prop}_30cm_avg"
    n_na = int(out.loc[use, ref].isna().sum())
    if use.sum():
        log.info("%s is NA for %d of %d points (%.1f%%): no-data pixels and points outside the map",
                 ref, n_na, int(use.sum()), 100 * n_na / use.sum())
    n_bad_depth = int((np.isfinite(depth) & np.array([layer_for(d) is None for d in depth])).sum())
    if n_bad_depth:
        log.warning("%d depth_cm value(s) outside 0-200 cm -> _depth_val left empty", n_bad_depth)
    if missing_layers:
        log.warning("layers not on the server: %s", ", ".join(sorted(missing_layers)))
    if failed:
        log.error("%d value(s) failed; rerun to read them (cached values are reused)", len(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())