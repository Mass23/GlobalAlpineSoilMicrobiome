#!/usr/bin/env python3
"""
Terrain (Copernicus DEM GLO-30) and vegetation index (MODIS NDVI) at sampling points.

All paths below are relative to the Snakefile's folder (matches chelsa_climate.py /
soilgrids_soil.py, called "other_data" in the Snakefile).
Reads   data/points.csv                              columns lat, lon, sample_date (yyyy-mm-dd)
Writes  ../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/other_data.csv
        data/tmp/other_data_cache.jsonl               fetched pixels; lets an interrupted run resume

MODIS doesn't publish a DEM; the nearest global 30 m product is Copernicus DEM GLO-30 (2019-2021
TanDEM-X-derived, WGS84 lon/lat tiles, no login). It's read through Microsoft Planetary Computer's
public STAC catalogue and Azure blob storage (free per-request SAS tokens, no account needed).

Terrain variables (from a small window of the DEM centred on each point):
  elevation_m          centre-pixel elevation
  slope_deg            Horn's (1981) method on the 3x3 neighbourhood
  aspect_deg           compass bearing of steepest descent (0=N, 90=E); empty on ~flat ground
  aspect_northness      cos(aspect) -- for regressions instead of the circular aspect_deg
  aspect_eastness       sin(aspect)
  tpi_m                 Topographic Position Index: elevation - mean of a ~200 m ring around it
  roughness_m           terrain ruggedness index: RMS elevation difference to the 3x3 neighbours
A point within ~1 pixel of a DEM tile edge gets a smaller, boundless-padded window; too few valid
pixels around a point (open water, a DEM void) leaves these columns empty and is noted in the log.

Vegetation index: MODIS MOD13Q1 v061, 250 m, 16-day composites (search: modis-13Q1-061), read the
same way. For each point the composite closest in time to sample_date is used, searching outward
(+-96 days) until one with a good/marginal pixel-reliability flag is found; otherwise the nearest
available value is used and flagged.
  ndvi                  NDVI of the composite used (-1 to 1)
  ndvi_date             date of that composite (yyyy-mm-dd)
  ndvi_days_offset       composite date minus sample_date, in days (0 = exact 16-day period)
  ndvi_reliability       0 good, 1 marginal, 2 snow/ice, 3 cloudy (MODIS pixel_reliability band)

Other DEM-derived variables worth adding later if useful: a wetness index (needs a full drainage
network from a mosaicked regional DEM, not a per-point window) and snow-cover duration (already in
the CHELSA extraction as `scd`).
"""
# no `from __future__` import: keeps the file usable with Snakemake's `script:` preamble too
import json
import logging
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.windows import Window
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import planetary_computer

# --------------------------------------------------------------------------- settings
POINTS_CSV = Path("/home/renku/work/GlobalAlpineSoilMicrobiome/data/points.csv")
OUT_CSV = Path("/home/renku/work/GlobalAlpineSoilMicrobiome/data/others/other_data.csv")
TMP_DIR = Path("data/tmp")
CACHE_JSONL = TMP_DIR / "other_data_cache.jsonl"

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token"

DEM_COLLECTION = "cop-dem-glo-30"
DEM_ASSET = "data"
DEM_TPI_RADIUS = 3                              # 7x7 window (~ +-105 m at 30 m pixels)

NDVI_COLLECTION = "modis-13Q1-061"
NDVI_ASSET = "250m_16_days_NDVI"
NDVI_QA_ASSET = "250m_16_days_pixel_reliability"
NDVI_SCALE = 0.0001
NDVI_FILL = -3000
NDVI_SEARCH_DAYS = 96                           # +-6 sixteen-day composites
NDVI_GOOD_RELIABILITY = {0, 1}                  # 0 good, 1 marginal (see docstring)
# MODIS sinusoidal grid: a plain spherical projection with a closed-form formula, so no CRS/PROJ
# lookup is needed (unlike a general reprojection) -- see the MODIS Land products user guide.
MODIS_SPHERE_R = 6371007.181

VSICURL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    GDAL_HTTP_MAX_RETRY="5",
    GDAL_HTTP_RETRY_DELAY="3",
    GDAL_HTTP_TIMEOUT="120",
)
N_WORKERS = 6
TRIES = 3
DEG_M = 111_320.0                               # metres per degree of latitude (WGS84, near enough)

log = logging.getLogger("other_data")


def coord_key(lon: float, lat: float) -> str:
    return f"{lon:.6f},{lat:.6f}"


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=6, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504))
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    return s


# --------------------------------------------------------------------------- STAC + signing
class Planetary:
    def __init__(self, session: requests.Session):
        self.s = session
        self._tokens: dict = {}                 # collection -> (token, expiry)

    def search(self, collection: str, lon: float, lat: float, eps: float, datetime_range=None):
        bbox = [lon - eps, lat - eps, lon + eps, lat + eps]
        body = {"collections": [collection], "bbox": bbox, "limit": 50}
        if datetime_range:
            body["datetime"] = datetime_range
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
        return planetary_computer.sign(href)
    #def sign(self, collection: str, href: str) -> str:
    #    token, expiry = self._tokens.get(collection, (None, None))
    #    if token is None or datetime.now(timezone.utc) > expiry - timedelta(minutes=2):
    #        r = self.s.get(f"{SAS_URL}/{collection}", timeout=30)
    #        r.raise_for_status()
    #        body = r.json()
    #        token = body["token"]
    #        expiry = datetime.fromisoformat(body["msft:expiry"].replace("Z", "+00:00"))
    #        self._tokens[collection] = (token, expiry)
    #    return href + ("&" if "?" in href else "?") + token


# --------------------------------------------------------------------------- DEM
def read_window(url: str, lon: float, lat: float, half: int):
    """(half*2+1)-square window centred on (lon, lat), boundless (nodata outside the file)"""
    with rasterio.Env(**VSICURL_ENV), rasterio.open(f"/vsicurl/{url}") as src:
        row, col = src.index(lon, lat)
        win = Window(col - half, row - half, half * 2 + 1, half * 2 + 1)
        arr = src.read(1, window=win, boundless=True, fill_value=src.nodata, masked=True)
        px_x, px_y = abs(src.transform.a), abs(src.transform.e)
    return arr, px_x, px_y


def fetch_dem_window(pc: Planetary, lon: float, lat: float):
    feats = pc.search(DEM_COLLECTION, lon, lat, eps=0.01)
    if not feats:
        return None
    href = pc.sign(DEM_COLLECTION, feats[0]["assets"][DEM_ASSET]["href"])
    arr, px_x_deg, px_y_deg = read_window(href, lon, lat, DEM_TPI_RADIUS)
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


# --------------------------------------------------------------------------- NDVI
def sinusoidal_xy(lon: float, lat: float):
    lon_r, lat_r = math.radians(lon), math.radians(lat)
    return MODIS_SPHERE_R * lon_r * math.cos(lat_r), MODIS_SPHERE_R * lat_r


def read_pixel_sinu(url: str, x: float, y: float):
    with rasterio.Env(**VSICURL_ENV), rasterio.open(f"/vsicurl/{url}") as src:
        row, col = src.index(x, y)               # affine-only; the sinusoidal CRS itself is unused
        if not (0 <= row < src.height and 0 <= col < src.width):
            return None
        v = src.read(1, window=Window(col, row, 1, 1), masked=True)[0, 0]
        return None if np.ma.is_masked(v) else int(v)


def fetch_ndvi(pc: Planetary, lon: float, lat: float, sample_date):
    start = (sample_date - timedelta(days=NDVI_SEARCH_DAYS)).strftime("%Y-%m-%dT00:00:00Z")
    end = (sample_date + timedelta(days=NDVI_SEARCH_DAYS)).strftime("%Y-%m-%dT23:59:59Z")
    feats = pc.search(NDVI_COLLECTION, lon, lat, eps=0.002, datetime_range=f"{start}/{end}")
    if not feats:
        return None
    x, y = sinusoidal_xy(lon, lat)

    #def item_date(f):
    #    return datetime.fromisoformat(f["properties"]["datetime"].replace("Z", "+00:00")).date()
    def item_date(f):
        props = f["properties"]

        dt = props.get("datetime")
        if dt is None:
            dt = props.get("start_datetime") or props.get("end_datetime")

        if dt is None:
            raise ValueError(f"No usable datetime in STAC item {f.get('id')}")

        return datetime.fromisoformat(dt.replace("Z", "+00:00")).date()
        
    feats.sort(key=lambda f: abs((item_date(f) - sample_date).days))
    best = None
    for f in feats:
        href = pc.sign(NDVI_COLLECTION, f["assets"][NDVI_ASSET]["href"])
        raw = read_pixel_sinu(href, x, y)
        if raw is None or raw == NDVI_FILL:
            continue
        qa_href = pc.sign(NDVI_COLLECTION, f["assets"][NDVI_QA_ASSET]["href"])
        rel = read_pixel_sinu(qa_href, x, y)
        rec = {"ndvi": raw * NDVI_SCALE, "ndvi_date": item_date(f).isoformat(),
               "ndvi_days_offset": (item_date(f) - sample_date).days,
               "ndvi_reliability": rel}
        if best is None:
            best = rec                            # nearest-in-time fallback if nothing is "good"
        if rel in NDVI_GOOD_RELIABILITY:
            return rec
    return best


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
    missing = {"lat", "lon", "sample_date"} - set(df.columns)
    if missing:
        raise ValueError(f"{path} lacks column(s) {sorted(missing)}; required: lat, lon, sample_date")
    lat = pd.to_numeric(df["lat"], errors="coerce")
    lon = pd.to_numeric(df["lon"], errors="coerce")
    bad = lat.isna() | lon.isna() | ~lat.between(-90, 90) | ~lon.between(-180, 180)
    if bad.any():
        raise ValueError(f"invalid lat/lon in rows {df.index[bad].tolist()[:10]} (decimal degrees, WGS84)")
    dates = pd.to_datetime(df["sample_date"].astype(str).str.strip(), format="%Y-%m-%d", errors="coerce")
    if dates.isna().any():
        bad_d = df.loc[dates.isna(), "sample_date"].astype(str).tolist()
        raise ValueError(f"{len(bad_d)} date(s) not in yyyy-mm-dd: {bad_d[:10]}")
    return df, lat.to_numpy(float), lon.to_numpy(float), dates


# --------------------------------------------------------------------------- main
def main() -> int:
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, force=True,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("rasterio").setLevel(logging.ERROR)   # PROJ chatter; reads here need no CRS
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    df, lat, lon, dates = read_points(POINTS_CSV)
    dem_keys = sorted({coord_key(lon[i], lat[i]) for i in range(len(df))})
    ndvi_keys = sorted({(coord_key(lon[i], lat[i]), dates.iloc[i].date().isoformat()) for i in range(len(df))})
    log.info("%d points, %d unique coordinates, %d unique coordinate x date pairs", len(df), len(dem_keys),
             len(ndvi_keys))

    cache = Cache(CACHE_JSONL)
    dem_todo = [k for k in dem_keys if ("dem", k) not in cache.data]
    ndvi_todo = [k for k in ndvi_keys if ("ndvi", "|".join(k)) not in cache.data]
    log.info("DEM: %d cached, %d to fetch. NDVI: %d cached, %d to fetch",
             len(dem_keys) - len(dem_todo), len(dem_todo), len(ndvi_keys) - len(ndvi_todo), len(ndvi_todo))

    session = make_session()
    pc = Planetary(session)
    failed = {}
    with ThreadPoolExecutor(N_WORKERS) as pool:
        futs = {}
        for ck in dem_todo:
            lo, la = map(float, ck.split(","))
            futs[pool.submit(fetch_dem_window, pc, lo, la)] = ("dem", ck)
        for ck, d in ndvi_todo:
            lo, la = map(float, ck.split(","))
            futs[pool.submit(fetch_ndvi, pc, lo, la, datetime.fromisoformat(d).date())] = ("ndvi", f"{ck}|{d}")
        for i, fut in enumerate(as_completed(futs), 1):
            kind, key = futs[fut]
            try:
                cache.add(kind, key, fut.result())
            except Exception as e:
                failed[(kind, key)] = str(e)
                log.error("FAILED %s %s: %s", kind, key, str(e)[:200])
            if i % 20 == 0 or i == len(futs):
                log.info("  %d/%d requests done", i, len(futs))

    # ---- assemble -----------------------------------------------------------
    terrain_cols = ["elevation_m", "slope_deg", "aspect_deg", "aspect_northness", "aspect_eastness",
                    "tpi_m", "roughness_m"]
    ndvi_cols = ["ndvi", "ndvi_date", "ndvi_days_offset", "ndvi_reliability"]
    out_rows = []
    no_dem, no_ndvi = 0, 0
    for i in range(len(df)):
        ck = coord_key(lon[i], lat[i])
        d = dates.iloc[i].date().isoformat()
        row = {}
        rec = cache.data.get(("dem", ck))
        terr = terrain_from_window(rec, lat[i]) if rec else None
        if terr is None:
            no_dem += 1
            row.update({c: np.nan for c in terrain_cols})
        else:
            row.update(terr)
        nrec = cache.data.get(("ndvi", f"{ck}|{d}"))
        if nrec is None:
            no_ndvi += 1
            row.update({c: np.nan for c in ndvi_cols})
        else:
            row.update(nrec)
        out_rows.append(row)

    out = pd.concat([df, pd.DataFrame(out_rows)], axis=1)
    for c in ["elevation_m", "slope_deg", "aspect_deg", "aspect_northness", "aspect_eastness",
              "tpi_m", "roughness_m", "ndvi"]:
        out[c] = out[c].astype(float).round(4)
    out.to_csv(OUT_CSV, index=False)
    log.info("wrote %s", OUT_CSV)

    if no_dem:
        log.warning("%d point(s) with no usable DEM pixel (tile gap, void, or ocean)", no_dem)
    if no_ndvi:
        log.warning("%d point(s) with no MODIS composite within +-%d days", no_ndvi, NDVI_SEARCH_DAYS)
    if failed:
        log.error("%d request(s) failed; rerun to retry them (cached values are reused)", len(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
