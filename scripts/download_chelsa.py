#!/usr/bin/env python3
"""
Extract CHELSA V2.1 climate values at dated GPS points.

Input   CSV with lat, lon, site, date (dd-mm-yyyy)          e.g. data/points.csv
Output  same rows + one column per variable                  e.g. data/climate/chelsa_climate.csv
        + <output>_provenance.csv (long format: which file(s) fed each value)
        + chelsa_download_plan.csv next to the output (the unique files to fetch)

Which layer is used, per point and variable
  monthly variables  clt cmi hurs pet pr rsds sfcWind spei12 spi12 tas tasmax tasmin vpd
      1. observed CHELSA monthly time series for the sampling month/year, if published
      2. otherwise the SSP monthly climatology (default ssp370) of the 30-yr period
         containing the sampling year (2011-2040 / 2041-2070 / 2071-2100), mean over GCMs
  annual variables   bio01..bio19 fcf fgd scd
      1. 1981-2010 climatology if the sampling year is <= 2010
      2. otherwise the SSP climatology of the period containing the year, mean over GCMs
  Cases CHELSA does not publish (e.g. SSP monthly vpd) are left empty and flagged
  in the provenance table.

Workflow
  1. list the CHELSA bucket (per variable) instead of hard-coding file names
  2. resolve every point x variable to the GeoTIFF(s) it needs -> unique download plan
  3. for each unique GeoTIFF: download once -> sample all its points -> delete
     sampled values go to a .jsonl cache, so an interrupted run resumes where it stopped

Usage
  python scripts/chelsa_climate.py --points data/points.csv \
         --output data/climate/chelsa_climate.csv [--dry-run] [--mode stream]
  or from Snakemake via `script:` (see chelsa.smk)
"""
# no `from __future__` import: Snakemake's `script:` prepends a preamble, which would break it
import argparse
import json
import logging
import re
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.warp import transform as warp_transform
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# --------------------------------------------------------------------------- constants
BIO_VARS = [f"bio{i:02d}" for i in range(1, 20)] + ["fcf", "fgd", "scd"]
MONTHLY_VARS = ["clt", "cmi", "hurs", "pet", "pr", "rsds", "sfcWind",
                "spei12", "spi12", "tas", "tasmax", "tasmin", "vpd"]
ALL_VARS = BIO_VARS + MONTHLY_VARS

GCMS = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
BASELINE = (1981, 2010)
FUTURE_PERIODS = [(2011, 2040), (2041, 2070), (2071, 2100)]

# other spellings a variable can have in folder / file names on the server
ALIASES = {"pet": ["pet_penman", "pet"]}

DEFAULTS = dict(
    base_url="https://os.unil.cloud.switch.ch/chelsa02",  # CHELSA V2.1 object store (S3 API)
    root="chelsa/global",
    ssp="ssp370",
    gcms=GCMS,
    mode="download",            # "download" = full GeoTIFF, "stream" = read pixels via /vsicurl/
    date_format="%d-%m-%Y",
    lat_col=None, lon_col=None, date_col=None,
    tmpdir=None,                # default: <output dir>/tmp_chelsa
    cache=None,                 # default: <output dir>/.chelsa_cache.jsonl
    no_scale=False,
    dry_run=False,
    log=None,
)

VSICURL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
    GDAL_HTTP_MAX_RETRY="5",
    GDAL_HTTP_RETRY_DELAY="3",
    GDAL_HTTP_TIMEOUT="120",
)

log = logging.getLogger("chelsa")


# --------------------------------------------------------------------------- helpers
def norm(text: str) -> str:
    """lower-case, every non-alphanumeric run -> '_', padded with '_' for token matching"""
    return "_" + re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_") + "_"


def variants(var: str) -> list[str]:
    if var.startswith("bio"):
        n = int(var[3:])
        return [f"bio{n:02d}", f"bio{n}"]
    return ALIASES.get(var, [var])


def var_pattern(var: str) -> str:
    return "(?:" + "|".join(re.escape(norm(v).strip("_")) for v in variants(var)) + ")"


def future_period(year: int):
    for a, b in FUTURE_PERIODS:
        if a <= year <= b:
            return (a, b)
    return None


def coord_key(lon: float, lat: float) -> str:
    return f"{lon:.6f},{lat:.6f}"


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=6, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=("GET", "HEAD"))
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    s.headers["User-Agent"] = "chelsa-point-extractor/1.0"
    return s


# --------------------------------------------------------------------------- bucket listing
class Bucket:
    """Anonymous S3 ListObjects (v1) client with per-prefix caching."""

    def __init__(self, session: requests.Session, base_url: str):
        self.s = session
        self.base = base_url.rstrip("/")
        self.sizes: dict[str, int] = {}
        self._cache: dict = {}

    def url(self, key: str) -> str:
        return f"{self.base}/{key}"

    def list(self, prefix: str, delimiter: str | None = None):
        ck = (prefix, delimiter)
        if ck in self._cache:
            return self._cache[ck]
        files, dirs, marker = {}, [], None
        while True:
            params = {"prefix": prefix}
            if delimiter:
                params["delimiter"] = delimiter
            if marker:
                params["marker"] = marker
            r = self.s.get(self.base, params=params, timeout=60)
            r.raise_for_status()
            root = ET.fromstring(r.content)
            for e in root.iter():                      # drop the S3 XML namespace
                e.tag = e.tag.split("}", 1)[-1]
            page = []
            for c in root.findall("Contents"):
                key = c.findtext("Key")
                files[key] = int(c.findtext("Size") or 0)
                page.append(key)
            for p in root.findall("CommonPrefixes"):
                d = p.findtext("Prefix")
                dirs.append(d)
                page.append(d)
            if root.findtext("IsTruncated") != "true" or not page:
                break
            marker = root.findtext("NextMarker") or max(page)
        self.sizes.update(files)
        self._cache[ck] = (files, dirs)
        return files, dirs


class Catalog:
    """Finds CHELSA GeoTIFFs by listing folders, tolerant to renamed/moved files."""

    def __init__(self, bucket: Bucket, root: str):
        self.b = bucket
        self.root = root.strip("/") + "/"
        _, dirs = bucket.list(self.root, "/")
        self.products = {d[len(self.root):].strip("/").lower(): d for d in dirs}
        log.info("CHELSA products under %s: %s", bucket.url(self.root), ", ".join(sorted(self.products)))
        if "monthly" not in self.products:
            raise RuntimeError(
                f"No 'monthly' folder under {bucket.url(self.root)} (found: {sorted(self.products)}). "
                "The CHELSA storage layout may have changed; browse https://envicloud.wsl.ch "
                "and adjust --base-url / --root.")
        self._folders: dict = {}

    def folders(self, product: str) -> dict:
        if product not in self._folders:
            pre = self.products.get(product)
            if pre is None:
                self._folders[product] = {}
            else:
                _, dirs = self.b.list(pre, "/")
                self._folders[product] = {d[len(pre):].strip("/").lower(): d for d in dirs}
        return self._folders[product]

    def tifs(self, products, var: str) -> dict:
        out = {}
        for product in products:
            folders = self.folders(product)
            for v in variants(var):
                pre = folders.get(v.lower())
                if pre:
                    files, _ = self.b.list(pre)
                    out.update({k: s for k, s in files.items() if k.lower().endswith(".tif")})
        return out


# --------------------------------------------------------------------------- resolution
@dataclass
class Source:
    kind: str                                   # timeseries | baseline_1981-2010 | ssp370_2011-2040 | missing
    keys: list = field(default_factory=list)
    note: str = ""


class Resolver:
    def __init__(self, catalog: Catalog, ssp: str, gcms: list[str]):
        self.c = catalog
        self.ssp = ssp.lower()
        self.gcms = list(gcms)
        self._memo: dict = {}
        self._warned: set = set()

    def _pick(self, hits, what):
        if not hits:
            return None
        hits = sorted(hits, key=lambda k: ("v_2_1" not in norm(k), len(k), k))
        if len(hits) > 1 and what not in self._warned:
            log.warning("%d candidate files for %s -> using %s", len(hits), what, hits[0])
            self._warned.add(what)
        return hits[0]

    def _ssp_keys(self, products, var, period, month=None):
        a, b = period
        files = self.c.tifs(products, var)
        vrx = re.compile(f"_{var_pattern(var)}_")
        keys, missing = [], []
        for g in self.gcms:
            gtok = norm(g)

            def ok(k):
                n, base = norm(k), norm(Path(k).name)
                if f"_{self.ssp}_" not in n or gtok not in n or f"_{a}_{b}_" not in n:
                    return False
                if not vrx.search(base):
                    return False
                return month is None or f"_{month:02d}_" in base

            k = self._pick([k for k in files if ok(k)], f"{var} {g} {self.ssp} {a}-{b} m{month}")
            if k:
                keys.append(k)
            else:
                missing.append(g)
        note = f"no file for GCM(s): {', '.join(missing)}" if keys and missing else ""
        return keys, note

    def monthly(self, var: str, year: int, month: int) -> Source:
        mk = ("m", var, year, month)
        if mk in self._memo:
            return self._memo[mk]
        files = self.c.tifs(["monthly"], var)
        rx = re.compile(rf"_{var_pattern(var)}_{month:02d}_{year}_(?!\d{{4}}_)")
        hit = self._pick([k for k in files if rx.search(norm(Path(k).name))], f"{var} {month:02d}-{year}")
        if hit:
            src = Source("timeseries", [hit])
        else:
            per = future_period(year)
            if per is None:
                src = Source("missing", note="month not in the time series, and no SSP period covers this year")
            else:
                keys, note = self._ssp_keys(["climatologies"], var, per, month)
                src = (Source(f"{self.ssp}_{per[0]}-{per[1]}", keys, note) if keys else
                       Source("missing", note=f"month not in the time series; no {self.ssp} monthly layer for {var}"))
        self._memo[mk] = src
        return src

    def annual(self, var: str, year: int) -> Source:
        mk = ("a", var, year)
        if mk in self._memo:
            return self._memo[mk]
        products = ["bioclim", "climatologies"]
        if year <= BASELINE[1]:
            a, b = BASELINE
            vrx = re.compile(f"_{var_pattern(var)}_")
            files = self.c.tifs(products, var)
            hits = [k for k in files
                    if f"_{a}_{b}_" in norm(k) and "_ssp" not in norm(k) and vrx.search(norm(Path(k).name))]
            hit = self._pick(hits, f"{var} {a}-{b}")
            note = "sampling year before 1981, 1981-2010 used" if year < a else ""
            src = Source(f"baseline_{a}-{b}", [hit], note) if hit else Source("missing", note=f"no {a}-{b} layer")
        else:
            per = future_period(year)
            if per is None:
                src = Source("missing", note="year beyond 2100")
            else:
                keys, note = self._ssp_keys(products, var, per)
                src = (Source(f"{self.ssp}_{per[0]}-{per[1]}", keys, note) if keys else
                       Source("missing", note=f"no {self.ssp} {per[0]}-{per[1]} layer for {var}"))
        self._memo[mk] = src
        return src


# --------------------------------------------------------------------------- cache of sampled values
class ValueCache:
    """key -> {coord_key: value}; appended as JSON lines so a crashed run can resume."""

    def __init__(self, path: Path | None):
        self.path = path
        self.data: dict[str, dict] = {}
        self.meta: dict[str, dict] = {}
        if path and path.exists():
            with open(path) as fh:
                for line in fh:
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:        # half-written last line after a crash
                        continue
                    self.data.setdefault(rec["key"], {}).update(rec["values"])
                    self.meta[rec["key"]] = {"scale": rec.get("scale"), "offset": rec.get("offset")}

    def missing(self, key: str, cks) -> list[str]:
        have = self.data.get(key, {})
        return sorted(c for c in cks if c not in have)

    def add(self, key: str, values: dict, meta: dict):
        self.data.setdefault(key, {}).update(values)
        self.meta[key] = meta
        if self.path:
            with open(self.path, "a") as fh:
                fh.write(json.dumps({"key": key, "values": values, **meta}) + "\n")

    def get(self, key: str, ck: str):
        return self.data.get(key, {}).get(ck)


# --------------------------------------------------------------------------- IO: download + sample
def download(session: requests.Session, url: str, dest: Path, expected: int | None, tries: int = 4) -> Path:
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, tries + 1):
        try:
            with session.get(url, stream=True, timeout=(30, 600)) as r:
                r.raise_for_status()
                with open(part, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 23):
                        fh.write(chunk)
            got = part.stat().st_size
            if expected and got != expected:
                raise OSError(f"incomplete download: {got} of {expected} bytes")
            part.replace(dest)
            return dest
        except (requests.RequestException, OSError) as e:
            part.unlink(missing_ok=True)
            if attempt == tries:
                raise
            wait = 15 * attempt
            log.warning("  attempt %d/%d failed (%s); retrying in %ds", attempt, tries, e, wait)
            time.sleep(wait)


def sample(path: str, cks: list[str], apply_scale: bool = True):
    lons = [float(c.split(",")[0]) for c in cks]
    lats = [float(c.split(",")[1]) for c in cks]
    with rasterio.open(path) as src:
        xs, ys = lons, lats
        if src.crs is not None and src.crs.to_epsg() != 4326:
            xs, ys = warp_transform("EPSG:4326", src.crs, lons, lats)
        scale = float(src.scales[0]) if (apply_scale and src.scales) else 1.0
        offset = float(src.offsets[0]) if (apply_scale and src.offsets) else 0.0
        vals = []
        for v in src.sample(zip(xs, ys), indexes=1, masked=True):
            if np.ma.getmaskarray(v)[0] or not np.isfinite(float(v[0])):
                vals.append(None)
            else:
                vals.append(float(v[0]) * scale + offset)
    return dict(zip(cks, vals)), {"scale": scale, "offset": offset}


# --------------------------------------------------------------------------- points
def read_points(path, lat_col, lon_col, date_col, date_format):
    df = pd.read_csv(path)
    lower = {c.strip().lower(): c for c in df.columns}

    def find(given, candidates, what):
        if given:
            if given not in df.columns:
                raise KeyError(f"column '{given}' not in {path} (columns: {list(df.columns)})")
            return given
        for c in candidates:
            if c in lower:
                return lower[c]
        raise KeyError(f"no {what} column in {path} (columns: {list(df.columns)}); set it explicitly")

    lat_c = find(lat_col, ["lat", "latitude", "y"], "latitude")
    lon_c = find(lon_col, ["lon", "long", "longitude", "lng", "x"], "longitude")
    date_c = find(date_col, ["date", "sampling_date"], "date")

    dates = pd.to_datetime(df[date_c].astype(str).str.strip(), format=date_format, errors="coerce")
    if dates.isna().any():
        bad = df.loc[dates.isna(), date_c].astype(str).tolist()
        raise ValueError(f"{len(bad)} date(s) not in format {date_format}: {bad[:10]}")
    lat = pd.to_numeric(df[lat_c], errors="coerce")
    lon = pd.to_numeric(df[lon_c], errors="coerce")
    bad = lat.isna() | lon.isna() | ~lat.between(-90, 90) | ~lon.between(-180, 180)
    if bad.any():
        raise ValueError(f"invalid coordinates in rows {df.index[bad].tolist()[:10]} "
                         f"(expected decimal degrees, WGS84)")
    site_c = lower.get("site")
    return df, lat.to_numpy(float), lon.to_numpy(float), dates, site_c, date_c


# --------------------------------------------------------------------------- main
def run(a) -> int:
    out_csv = Path(a.output)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    prov_csv = Path(a.provenance) if getattr(a, "provenance", None) else \
        out_csv.with_name(out_csv.stem + "_provenance.csv")
    plan_csv = out_csv.with_name("chelsa_download_plan.csv")
    cache_path = Path(a.cache) if a.cache else out_csv.parent / ".chelsa_cache.jsonl"
    tmp_root = Path(a.tmpdir) if a.tmpdir else out_csv.parent / "tmp_chelsa"

    df, lat, lon, dates, site_c, date_c = read_points(a.points, a.lat_col, a.lon_col, a.date_col, a.date_format)
    log.info("%d points read from %s (%s to %s)", len(df), a.points,
             dates.min().strftime("%d-%m-%Y"), dates.max().strftime("%d-%m-%Y"))

    # ---- 1+2: list the bucket and resolve every point x variable -------------
    session = make_session()
    bucket = Bucket(session, a.base_url)
    resolver = Resolver(Catalog(bucket, a.root), a.ssp, a.gcms)

    plans, need, need_vars = [], {}, {}
    for i in range(len(df)):
        y, m = dates.iloc[i].year, dates.iloc[i].month
        ck = coord_key(lon[i], lat[i])
        p = {v: resolver.annual(v, y) for v in BIO_VARS}
        p.update({v: resolver.monthly(v, y, m) for v in MONTHLY_VARS})
        for var, src in p.items():
            for k in src.keys:
                need.setdefault(k, set()).add(ck)
                need_vars.setdefault(k, set()).add(var)
        plans.append(p)

    kinds = Counter(src.kind.split("_")[0] for p in plans for src in p.values())
    total_gb = sum(bucket.sizes.get(k, 0) for k in need) / 1e9
    log.info("values to fill: %s", dict(kinds))
    log.info("download plan: %d unique GeoTIFFs, %.1f GB in total", len(need), total_gb)
    pd.DataFrame(
        [{"key": k, "url": bucket.url(k), "size_mb": round(bucket.sizes.get(k, 0) / 1e6, 1),
          "n_coords": len(c), "variables": " ".join(sorted(need_vars[k]))} for k, c in sorted(need.items())]
    ).to_csv(plan_csv, index=False)
    log.info("plan written to %s", plan_csv)
    if a.dry_run:
        return 0

    # ---- 3: download each file once, sample, delete ------------------------
    cache = ValueCache(cache_path)
    todo = [(k, cache.missing(k, c)) for k, c in sorted(need.items())]
    todo = [(k, c) for k, c in todo if c]
    log.info("%d files already in cache, %d to fetch (%.1f GB)", len(need) - len(todo), len(todo),
             sum(bucket.sizes.get(k, 0) for k, _ in todo) / 1e9)

    failed = {}
    tmp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=tmp_root) as tmp:
        for i, (key, cks) in enumerate(todo, 1):
            url = bucket.url(key)
            log.info("[%d/%d] %s (%.0f MB, %d points)", i, len(todo), Path(key).name,
                     bucket.sizes.get(key, 0) / 1e6, len(cks))
            try:
                if a.mode == "stream":
                    with rasterio.Env(**VSICURL_ENV):
                        vals, meta = sample(f"/vsicurl/{url}", cks, not a.no_scale)
                else:
                    dest = Path(tmp) / Path(key).name
                    try:
                        download(session, url, dest, bucket.sizes.get(key))
                        vals, meta = sample(str(dest), cks, not a.no_scale)
                    finally:
                        dest.unlink(missing_ok=True)
            except Exception as e:                      # keep going, report at the end
                log.error("  FAILED %s: %s", key, e)
                failed[key] = str(e)
                continue
            cache.add(key, vals, meta)

    # ---- assemble -----------------------------------------------------------
    out = df.copy()
    clash = [v for v in ALL_VARS if v in out.columns]
    if clash:
        log.warning("overwriting existing columns in points file: %s", clash)
    cols = {v: np.full(len(df), np.nan) for v in ALL_VARS}
    prov = []
    for i, p in enumerate(plans):
        ck = coord_key(lon[i], lat[i])
        for var, src in p.items():
            vals = [cache.get(k, ck) for k in src.keys]
            ok = [v for v in vals if v is not None]
            if ok:
                cols[var][i] = float(np.mean(ok))
            note = src.note
            if any(k in failed for k in src.keys):
                note = (note + "; " if note else "") + "download failed"
            prov.append({
                "row": i,
                "site": df[site_c].iloc[i] if site_c else None,
                "date": df[date_c].iloc[i],
                "lat": lat[i], "lon": lon[i],
                "variable": var,
                "source": src.kind,
                "n_files": len(src.keys),
                "n_values": len(ok),
                "value": cols[var][i],
                "scale_offset": ";".join(f"{cache.meta.get(k, {}).get('scale')}/{cache.meta.get(k, {}).get('offset')}"
                                         for k in src.keys),
                "files": ";".join(Path(k).name for k in src.keys),
                "note": note,
            })
    for v in ALL_VARS:
        out[v] = np.round(cols[v], 4)
    out.to_csv(out_csv, index=False)
    pd.DataFrame(prov).to_csv(prov_csv, index=False)
    log.info("wrote %s and %s", out_csv, prov_csv)

    empty = [v for v in ALL_VARS if out[v].isna().all()]
    if empty:
        log.warning("variables with no value at any point: %s (see provenance)", " ".join(empty))
    if failed:
        log.error("%d file(s) failed; rerun to retry them (cached values are reused)", len(failed))
        return 1
    return 0


def setup_logging(logfile=None):
    handlers = [logging.StreamHandler(sys.stderr)]
    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(logfile, mode="w"))
    logging.basicConfig(level=logging.INFO, handlers=handlers, force=True,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")


def args_from_snakemake(smk) -> argparse.Namespace:
    p = smk.params
    a = dict(DEFAULTS)
    for name in DEFAULTS:
        if hasattr(p, name):
            a[name] = getattr(p, name)
    a["points"] = str(getattr(smk.input, "points", smk.input[0]))
    a["output"] = str(getattr(smk.output, "csv", smk.output[0]))
    a["provenance"] = str(getattr(smk.output, "provenance", "")) or None
    a["log"] = str(smk.log[0]) if len(smk.log) else None
    if isinstance(a["gcms"], str):
        a["gcms"] = a["gcms"].split(",")
    return argparse.Namespace(**a)


def args_from_cli(argv=None) -> argparse.Namespace:
    d = DEFAULTS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", default="data/points.csv")
    ap.add_argument("--output", default="data/climate/chelsa_climate.csv")
    ap.add_argument("--provenance", default=None, help="default: <output>_provenance.csv")
    ap.add_argument("--ssp", default=d["ssp"])
    ap.add_argument("--gcms", default=",".join(d["gcms"]), help="comma-separated; values are averaged")
    ap.add_argument("--mode", choices=["download", "stream"], default=d["mode"])
    ap.add_argument("--base-url", default=d["base_url"])
    ap.add_argument("--root", default=d["root"])
    ap.add_argument("--date-format", default=d["date_format"])
    ap.add_argument("--lat-col"), ap.add_argument("--lon-col"), ap.add_argument("--date-col")
    ap.add_argument("--tmpdir"), ap.add_argument("--cache"), ap.add_argument("--log")
    ap.add_argument("--no-scale", action="store_true", help="return raw stored integers")
    ap.add_argument("--dry-run", action="store_true", help="only build and write the download plan")
    a = ap.parse_args(argv)
    a.gcms = [g.strip() for g in a.gcms.split(",") if g.strip()]
    return a


if __name__ == "__main__":
    args = args_from_snakemake(snakemake) if "snakemake" in globals() else args_from_cli()  # noqa: F821
    setup_logging(args.log)
    sys.exit(run(args))