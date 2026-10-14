#!/usr/bin/env python3
"""
Terrain (Copernicus DEM GLO-30) at sampling points.

Run     python3 download_others.py [points.csv other_data.csv]   (Snakemake passes the two paths; the
        tmp/ cache is written next to the output)
Reads   points.csv      columns lat, lon; all other columns are kept in the output (the date is not used)
Writes  other_data.csv  points + the terrain columns below
        tmp/other_data_cache.jsonl   fetched pixels; lets an interrupted run resume

Copernicus DEM GLO-30 (2019-2021 TanDEM-X-derived, WGS84 lon/lat tiles, no login). DEM_SOURCE = "aws"
(default) reads its public COG tiles from the AWS Open Data bucket by tile name: no search, no token, no
account. DEM_SOURCE = "planetary" reads the same product through Microsoft Planetary Computer's STAC
catalogue (which may answer 403 to anonymous bursts: then set PC_SDK_SUBSCRIPTION_KEY, see make_session).

Terrain variables (from a small window of the DEM centred on each point):
  elevation_m          centre-pixel elevation
  slope_deg            Horn's (1981) method on the 3x3 neighbourhood
  aspect_deg           compass bearing of steepest descent (0=N, 90=E); empty on ~flat ground
  aspect_northness      cos(aspect) -- for regressions instead of the circular aspect_deg
  aspect_eastness       sin(aspect)
  tpi_m                 Topographic Position Index: elevation - mean of a ~200 m ring around it
  roughness_m           terrain ruggedness index: RMS elevation difference to the 3x3 neighbours
A point within ~1 pixel of a DEM tile edge gets a smaller, boundless-padded window. A point without
terrain values is classified and logged (and listed in dem_gaps.csv next to the output):
  no_tile        the DEM has no tile there (open water, some tiny islands)       -> NA, not retried
  no_data_pixel  the tile exists but the pixel at the point is no-data           -> NA
  failed         network / server trouble (retried; exit code 1, a rerun retries only those)
Terrain does not depend on the sampling date, so every point gets values whatever its date looks like.

Other DEM-derived variables worth adding later if useful: a wetness index (needs a full drainage
network from a mosaicked regional DEM, not a per-point window) and snow-cover duration (already in
the CHELSA extraction as `scd`).
"""
# no `from __future__` import: keeps the file usable with Snakemake's `script:` preamble too
import json
import logging
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.windows import Window
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# --------------------------------------------------------------------------- settings
POINTS_CSV = Path("/home/renku/work/GlobalAlpineSoilMicrobiome/data/points.csv")
OUT_CSV = Path("/home/renku/work/GlobalAlpineSoilMicrobiome/data/other_data.csv")
if len(sys.argv) >= 3:                     # Snakemake: download_others.py <points.csv> <other_data.csv>
    POINTS_CSV, OUT_CSV = Path(sys.argv[1]), Path(sys.argv[2])
TMP_DIR = OUT_CSV.parent / "tmp"
CACHE_JSONL = TMP_DIR / "other_data_cache.jsonl"

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"   # only for DEM_SOURCE = "planetary"

DEM_COLLECTION = "cop-dem-glo-30"
DEM_ASSET = "data"
DEM_TPI_RADIUS = 3                              # 7x7 window (~ +-105 m at 30 m pixels)

VSICURL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    GDAL_HTTP_MAX_RETRY="5",
    GDAL_HTTP_RETRY_DELAY="3",
    GDAL_HTTP_TIMEOUT="120",
)
N_WORKERS = 3                                   # polite: bursts are what gets anonymous clients refused
TRIES = 3
DEM_SOURCE = "aws"                              # "aws" (public tiles, no search) | "planetary" (STAC search)
DEM_AWS_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"
MAX_CONSECUTIVE_FAILURES = 25                   # stop instead of hammering a service that keeps refusing
DEG_M = 111_320.0                               # metres per degree of latitude (WGS84, near enough)

log = logging.getLogger("other_data")


def coord_key(lon: float, lat: float) -> str:
    return f"{lon:.6f},{lat:.6f}"


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers["User-Agent"] = "GlobalAlpineSoilMicrobiome-download_others/1.0 (research use)"
    key = os.environ.get("PC_SDK_SUBSCRIPTION_KEY")      # optional: less restricted rate limiting on Planetary Computer
    if key:
        s.headers["Ocp-Apim-Subscription-Key"] = key
    retry = Retry(total=6, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504))
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    return s


# --------------------------------------------------------------------------- STAC + signing (DEM_SOURCE = "planetary")
class Planetary:
    def __init__(self, session: requests.Session):
        self.s = session

    def search(self, collection: str, lon: float, lat: float, eps: float):
        bbox = [lon - eps, lat - eps, lon + eps, lat + eps]
        body = {"collections": [collection], "bbox": bbox, "limit": 50}
        for attempt in range(1, TRIES + 1):
            try:
                r = self.s.post(STAC_URL, json=body, timeout=60)
                r.raise_for_status()
                return r.json().get("features", [])
            except requests.RequestException:
                if attempt == TRIES:
                    raise
                time.sleep(5 * attempt)

    def sign(self, collection: str, href: str) -> str:
        import planetary_computer              # only needed for this source, so the package is optional
        return planetary_computer.sign(href)


# --------------------------------------------------------------------------- DEM
def read_window(url: str, lon: float, lat: float, half: int):
    """(half*2+1)-square window centred on (lon, lat), boundless (nodata outside the file)"""
    with rasterio.Env(**VSICURL_ENV), rasterio.open(f"/vsicurl/{url}") as src:
        row, col = src.index(lon, lat)
        win = Window(col - half, row - half, half * 2 + 1, half * 2 + 1)
        arr = src.read(1, window=win, boundless=True, fill_value=src.nodata, masked=True)
        px_x, px_y = abs(src.transform.a), abs(src.transform.e)
    return arr, px_x, px_y


def dem_tile_url(lon: float, lat: float) -> str:
    """public Copernicus GLO-30 COG tile (1 x 1 degree, named by its south-west corner)"""
    lat0, lon0 = math.floor(lat), math.floor(lon)
    name = (f"Copernicus_DSM_COG_10_{'N' if lat0 >= 0 else 'S'}{abs(lat0):02d}_00_"
            f"{'E' if lon0 >= 0 else 'W'}{abs(lon0):03d}_00_DEM")
    return f"{DEM_AWS_BASE}/{name}/{name}.tif"


NO_TILE = "no_tile"                             # fetch_dem_window result: the DEM has no tile at this point
_TILE_EXISTS: dict = {}                         # tile url -> bool (a 1 x 1 degree tile is asked about once)


def tile_exists(pc: Planetary, url: str) -> bool:
    """True / False for a tile that is / is not in the bucket; anything else (5xx, timeouts) raises"""
    if url in _TILE_EXISTS:
        return _TILE_EXISTS[url]
    r = pc.s.head(url, timeout=60, allow_redirects=True)
    if r.status_code == 200:
        ok = True
    elif r.status_code in (403, 404):           # S3 answers 403 for a key that does not exist
        ok = False
    else:
        raise RuntimeError(f"HTTP {r.status_code} for {url}")
    _TILE_EXISTS[url] = ok
    return ok


def fetch_dem_window(pc: Planetary, lon: float, lat: float):
    """a record {grid, px_x_deg, px_y_deg}, or NO_TILE; read errors are NOT swallowed (they are real failures)"""
    if DEM_SOURCE == "aws":
        url = dem_tile_url(lon, lat)
        if not tile_exists(pc, url):
            return NO_TILE
    else:
        feats = pc.search(DEM_COLLECTION, lon, lat, eps=0.01)
        if not feats:
            return NO_TILE
        url = pc.sign(DEM_COLLECTION, feats[0]["assets"][DEM_ASSET]["href"])
    for attempt in range(1, TRIES + 1):
        try:
            arr, px_x_deg, px_y_deg = read_window(url, lon, lat, DEM_TPI_RADIUS)
            break
        except Exception:                        # GDAL has already retried the HTTP request itself
            if attempt == TRIES:
                raise
            time.sleep(5 * attempt)
    return {
        "grid": [[None if np.ma.is_masked(v) else float(v) for v in row] for row in arr],
        "px_x_deg": px_x_deg, "px_y_deg": px_y_deg,
    }


def terrain_from_window(rec: dict, lat: float):
    g = np.array([[np.nan if v is None else v for v in row] for row in rec["grid"]])
    r = DEM_TPI_RADIUS
    if not np.isfinite(g[r, r]):
        return None
    elevation = float(g[r, r])
    dx = rec["px_x_deg"] * DEG_M * math.cos(math.radians(lat))
    dy = rec["px_y_deg"] * DEG_M

    c = g[r - 1:r + 2, r - 1:r + 2]              # 3x3 around the centre
    slope = aspect_deg = north = east = None
    if np.isfinite(c).all():
        dzdx = ((c[0, 2] + 2 * c[1, 2] + c[2, 2]) - (c[0, 0] + 2 * c[1, 0] + c[2, 0])) / (8 * dx)
        dzdy = ((c[2, 0] + 2 * c[2, 1] + c[2, 2]) - (c[0, 0] + 2 * c[0, 1] + c[0, 2])) / (8 * dy)
        slope = math.degrees(math.atan(math.hypot(dzdx, dzdy)))
        if slope > 1e-6:
            aspect_deg = (90.0 - math.degrees(math.atan2(dzdy, -dzdx))) % 360.0
            north = math.cos(math.radians(aspect_deg))
            east = math.sin(math.radians(aspect_deg))

    ring = g.copy(); ring[r, r] = np.nan
    valid = np.isfinite(ring)
    tpi = elevation - float(np.nanmean(ring)) if valid.any() else None
    rough = float(np.sqrt(np.nanmean((ring[valid] - elevation) ** 2))) if valid.any() else None

    return dict(elevation_m=elevation, slope_deg=slope, aspect_deg=aspect_deg,
                aspect_northness=north, aspect_eastness=east, tpi_m=tpi, roughness_m=rough)


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
                    self.data[(r["kind"], r["key"])] = r["value"]

    def add(self, kind, key, value):
        self.data[(kind, key)] = value
        with open(self.path, "a") as fh:
            fh.write(json.dumps({"kind": kind, "key": key, "value": value}) + "\n")


# --------------------------------------------------------------------------- points
def read_points(path: Path):
    df = pd.read_csv(path)
    missing = {"lat", "lon"} - set(df.columns)
    if missing:
        raise ValueError(f"{path} lacks column(s) {sorted(missing)}; required: lat, lon")
    lat = pd.to_numeric(df["lat"], errors="coerce")
    lon = pd.to_numeric(df["lon"], errors="coerce")
    bad = lat.isna() | lon.isna() | ~lat.between(-90, 90) | ~lon.between(-180, 180)
    if bad.any():
        raise ValueError(f"invalid lat/lon in rows {df.index[bad].tolist()[:10]} (decimal degrees, WGS84)")
    return df, lat.to_numpy(float), lon.to_numpy(float)


# --------------------------------------------------------------------------- main
def main() -> int:
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, force=True,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("rasterio").setLevel(logging.ERROR)   # PROJ chatter; reads here need no CRS
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    df, lat, lon = read_points(POINTS_CSV)
    dem_keys = sorted({coord_key(lon[i], lat[i]) for i in range(len(df))})
    log.info("%d points, %d unique coordinates", len(df), len(dem_keys))

    cache = Cache(CACHE_JSONL)
    # an earlier version cached read errors as "None" (= no DEM): those are fetched again
    dem_todo = [k for k in dem_keys if cache.data.get(("dem", k)) is None and ("dem_gap", k) not in cache.data]
    log.info("DEM: %d cached, %d to fetch", len(dem_keys) - len(dem_todo), len(dem_todo))

    session = make_session()
    pc = Planetary(session)
    failed = {}
    if dem_todo:                                   # one known land tile (the Alps): fail early, with the reason
        try:
            ok = fetch_dem_window(pc, 8.5, 46.5)
        except Exception as e:
            raise RuntimeError(f"DEM preflight failed (DEM_SOURCE = {DEM_SOURCE!r}): {e}") from e
        if ok == NO_TILE:
            raise RuntimeError(f"DEM preflight: no DEM tile at 8.5E 46.5N (DEM_SOURCE = {DEM_SOURCE!r}); "
                               f"check the network / the tile URL {dem_tile_url(8.5, 46.5)}")
    with ThreadPoolExecutor(N_WORKERS) as pool:
        futs = {}
        for ck in dem_todo:
            lo, la = map(float, ck.split(","))
            futs[pool.submit(fetch_dem_window, pc, lo, la)] = ("dem", ck)
        consecutive = 0
        for i, fut in enumerate(as_completed(futs), 1):
            kind, key = futs[fut]
            try:
                res = fut.result()
                if res == NO_TILE:
                    cache.add("dem_gap", key, NO_TILE)      # deterministic: not retried
                else:
                    cache.add(kind, key, res)
                consecutive = 0
            except Exception as e:
                failed[(kind, key)] = str(e)
                consecutive += 1
                log.error("FAILED %s %s: %s", kind, key, str(e)[:200])
                if consecutive >= MAX_CONSECUTIVE_FAILURES:
                    for f in futs:
                        f.cancel()
                    log.error("%d requests in a row failed: stopping (cached values are kept; rerun later). "
                              "If the DEM source keeps refusing: try DEM_SOURCE = 'aws', or set "
                              "PC_SDK_SUBSCRIPTION_KEY for 'planetary'.", consecutive)
                    break
            if i % 20 == 0 or i == len(futs):
                log.info("  %d/%d requests done", i, len(futs))

    # ---- assemble -----------------------------------------------------------
    terrain_cols = ["elevation_m", "slope_deg", "aspect_deg", "aspect_northness", "aspect_eastness",
                    "tpi_m", "roughness_m"]
    out_rows, reasons = [], []
    for i in range(len(df)):
        ck = coord_key(lon[i], lat[i])
        rec = cache.data.get(("dem", ck))
        terr = terrain_from_window(rec, lat[i]) if rec else None
        if terr is not None:
            out_rows.append(terr); reasons.append("")
        else:
            out_rows.append({c: np.nan for c in terrain_cols})
            reasons.append(("no_data_pixel" if rec else cache.data.get(("dem_gap", ck), "failed")))
    out = pd.concat([df, pd.DataFrame(out_rows)], axis=1)
    for c in terrain_cols:
        out[c] = out[c].astype(float).round(4)
    out.to_csv(OUT_CSV, index=False, na_rep="NA")
    log.info("wrote %s", OUT_CSV)

    gaps = pd.DataFrame({"lon": lon, "lat": lat, "reason": reasons})
    gaps = gaps[gaps["reason"] != ""]
    if len(gaps):
        per_coord = gaps.groupby(["lon", "lat", "reason"]).size().rename("n_points").reset_index()
        per_coord.to_csv(OUT_CSV.with_name("dem_gaps.csv"), index=False)
        for reason, g in per_coord.groupby("reason"):
            ex = ", ".join(f"({r.lon:.3f}E {r.lat:.3f}N)" for r in g.head(6).itertuples())
            log.warning("%d point(s) / %d coordinate(s) without terrain values: %s   e.g. %s",
                        int(g.n_points.sum()), len(g), reason, ex)
        log.info("all of them are listed in %s", OUT_CSV.with_name("dem_gaps.csv"))
    if failed:
        log.error("%d request(s) failed; rerun to retry them (cached values are reused)", len(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())