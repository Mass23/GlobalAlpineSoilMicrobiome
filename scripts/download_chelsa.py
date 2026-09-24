#!/usr/bin/env python3
"""
CHELSA V2.1 climate at dated sampling points.

All paths below are under PROJECT_DIR (/home/renku/work/GlobalAlpineSoilMicrobiome).
Reads   data/points.csv                               columns lat, lon, sample_date (yyyy-mm-dd); others kept
Writes  data/climate/chelsa_climate.csv               points + bio01..bio19 fcf fgd scd clt ... vpd
        data/climate/chelsa_climate_provenance.csv    per point x variable: source, anchors, files
        data/tmp/chelsa_plan.csv                      unique CHELSA files the run reads
        data/tmp/chelsa_cache.jsonl                   sampled values; lets an interrupted run resume

Which layer is used, per point and variable
  monthly variables  clt cmi hurs pet pr rsds sfcWind spei12 spi12 tas tasmax tasmin vpd
      observed CHELSA month of the sampling date if published, otherwise projection (below)
  annual variables   bio01..bio19 fcf fgd scd
      1981-2010 climatology if the sampling year is <= 2010, otherwise projection (below)

Projection = SSP370, median of the 5 CHELSA GCMs, linearly interpolated in time:
  each 30-yr period mean sits at its period centre (2011-2040 -> 2025.5, + month centre);
  the sampling time is interpolated between the two anchors that bracket it
    before the 2011-2040 centre : past anchor -> 2011-2040
        monthly: mean of the same calendar month over the last 5 published years
        annual : 1981-2010 climatology at 1995.5 (no observed annual series exists)
    between two period centres  : e.g. 2011-2040 -> 2041-2070
    after the 2071-2100 centre  : held at the 2071-2100 value
  times are decimal years: month m of year y = y + (m - 0.5) / 12, annual values = y + 0.5
The monthly past anchor uses the 5 years ending at the last COMPLETE observed year (all 12 months
published), so a partly published final year is ignored.
Monthly variables without SSP monthly layers (hurs clt cmi pet rsds sfcWind vpd) get their period
value by the delta method from CHELSA's annual-mean layers (e.g. hurs_mean), per GCM, then median:
    1981-2010 climatology of that month + (GCM 2011-2040 mean - 1981-2010 mean)   hurs, clt, cmi
    1981-2010 climatology of that month x (GCM 2011-2040 mean / 1981-2010 mean)   pet, rsds, sfcWind, vpd
and are interpolated exactly like the others. spei12/spi12 have no projection at all: they fall
back to the mean of the recent observed years (source 'recent_obs').
tas/tasmax/tasmin are stored in Kelvin in the monthly files and are converted to degC.
Anything still unavailable is left empty; the provenance note says why.

GeoTIFFs are read remotely through GDAL /vsicurl/ (only the needed pixels), each file opened
once for all points that need it. Values are in physical units (each file's scale/offset applied).
All paths hang off PROJECT_DIR (absolute), so the script can be run from anywhere.
"""
# no `from __future__` import: keeps the file usable with Snakemake's `script:` preamble too
import json
import logging
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# --------------------------------------------------------------------------- settings
PROJECT_DIR = Path("/home/renku/work/GlobalAlpineSoilMicrobiome")
POINTS_CSV = PROJECT_DIR / "data/points.csv"
OUT_CSV = PROJECT_DIR / "data/climate/chelsa_climate.csv"
PROV_CSV = PROJECT_DIR / "data/climate/chelsa_climate_provenance.csv"
TMP_DIR = PROJECT_DIR / "data/tmp"
PLAN_CSV = TMP_DIR / "chelsa_plan.csv"
CACHE_JSONL = TMP_DIR / "chelsa_cache.jsonl"

# CHELSA V2.1 object store (public S3 bucket run by WSL) and the folder of the global
# products inside it; browsable at https://envicloud.wsl.ch
CHELSA_BUCKET = "https://os.unil.cloud.switch.ch/chelsa02"
CHELSA_PREFIX = "chelsa/global/"

SSP = "ssp370"
GCMS = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
N_RECENT = 5                       # observed years averaged for the monthly past anchor
DATE_FORMAT = "%Y-%m-%d"

BIO_VARS = [f"bio{i:02d}" for i in range(1, 20)] + ["fcf", "fgd", "scd"]
MONTHLY_VARS = ["clt", "cmi", "hurs", "pet", "pr", "rsds", "sfcWind",
                "spei12", "spi12", "tas", "tasmax", "tasmin", "vpd"]
ALL_VARS = BIO_VARS + MONTHLY_VARS

BASELINE = (1981, 2010)
FUTURE_PERIODS = [(2011, 2040), (2041, 2070), (2071, 2100)]
ALIASES = {"pet": ["pet_penman", "pet"]}      # spellings used in CHELSA folder / file names

KELVIN_VARS = {"tas", "tasmax", "tasmin"}      # monthly files store K (offset 0): converted to degC
# monthly variables without SSP monthly layers: projected via the change of CHELSA's annual-mean
# layer (delta method); additive for bounded/signed variables, multiplicative for positive ones
MEAN_LAYER = {"hurs": "hurs_mean", "clt": "clt_mean", "cmi": "cmi_mean", "pet": "pet_penman_mean",
              "rsds": "rsds_mean", "sfcWind": "sfcWind_mean", "vpd": "vpd_mean"}
DELTA_OP = {"hurs": "add", "clt": "add", "cmi": "add",
            "pet": "mul", "rsds": "mul", "sfcWind": "mul", "vpd": "mul"}

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


def dec_time(year: float, month=None) -> float:
    """decimal-year centre of a month (or of the year if month is None)"""
    return year + ((month - 0.5) / 12 if month else 0.5)


def period_centre(period, month=None) -> float:
    return dec_time((period[0] + period[1]) / 2, month)


def future_period(year: int):
    for a, b in FUTURE_PERIODS:
        if a <= year <= b:
            return (a, b)
    return None


def coord_key(lon: float, lat: float) -> str:
    return f"{lon:.6f},{lat:.6f}"


# --------------------------------------------------------------------------- bucket listing
class Bucket:
    """Anonymous S3 ListObjects (v1) client with per-prefix caching."""

    def __init__(self, base_url: str):
        self.base = base_url.rstrip("/")
        self.s = requests.Session()
        retry = Retry(total=6, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504))
        self.s.mount("https://", HTTPAdapter(max_retries=retry))
        self.s.mount("http://", HTTPAdapter(max_retries=retry))
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

    def __init__(self, bucket: Bucket):
        self.b = bucket
        _, dirs = bucket.list(CHELSA_PREFIX, "/")
        self.products = {d[len(CHELSA_PREFIX):].strip("/").lower(): d for d in dirs}
        log.info("CHELSA products under %s: %s", bucket.url(CHELSA_PREFIX), ", ".join(sorted(self.products)))
        if "monthly" not in self.products:
            raise RuntimeError(
                f"No 'monthly' folder under {bucket.url(CHELSA_PREFIX)} (found: {sorted(self.products)}). "
                "The CHELSA storage layout may have changed; check https://envicloud.wsl.ch "
                "and update CHELSA_BUCKET / CHELSA_PREFIX.")
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
    """How one point x variable is filled.

    anchors: [(time, [keys], how)]
      how = "mean" | "median"            aggregate of the anchor's files
      how = {"op": "add"|"mul", "clim": key, "base": key}
            delta anchor: keys are per-GCM annual-mean projections; each GCM gives
            clim (+|x) (gcm - base) resp. (gcm / base); the median over GCMs is used
    A single anchor with time None is used as is; otherwise the anchor values are
    linearly interpolated at time t.
    """
    kind: str               # timeseries | baseline_1981-2010 | interp_ssp370 | interp_delta_ssp370 | recent_obs | missing
    anchors: list = field(default_factory=list)
    t: float | None = None
    note: str = ""

    @property
    def keys(self) -> list:
        out = []
        for _, ks, how in self.anchors:
            out += ks
            if isinstance(how, dict):
                out += [how["clim"], how["base"]]
        return out


def tokens_ok(tokens):
    """basename matcher: every token must appear as a whole token in the normalised file name"""
    return lambda base_norm: all(f"_{t}_" in base_norm for t in tokens)


class Resolver:
    def __init__(self, catalog: Catalog):
        self.c = catalog
        self._memo: dict = {}
        self._obs_index: dict = {}
        self._warned: set = set()

    def _pick(self, hits, what):
        if not hits:
            return None
        hits = sorted(hits, key=lambda k: ("v_2_1" not in norm(k), len(k), k))
        if len(hits) > 1 and what not in self._warned:
            log.warning("%d candidate files for %s -> using %s", len(hits), what, hits[0])
            self._warned.add(what)
        return hits[0]

    # ---- file finders
    def _observed(self, var) -> dict:
        """{(month, year): key} for the published monthly time series of `var`"""
        if var not in self._obs_index:
            rx = re.compile(rf"_{var_pattern(var)}_(\d{{2}})_(\d{{4}})_(?!\d{{4}}_)")
            found: dict = {}
            for k in self.c.tifs(["monthly"], var):
                m = rx.search(norm(Path(k).name))
                if m:
                    found.setdefault((int(m.group(1)), int(m.group(2))), []).append(k)
            self._obs_index[var] = {my: self._pick(ks, f"{var} {my[0]:02d}-{my[1]}") for my, ks in found.items()}
        return self._obs_index[var]

    def _last_complete_year(self, obs):
        full = [y for y in {y for (_, y) in obs} if all((m, y) in obs for m in range(1, 13))]
        return max(full) if full else None

    def _baseline_key(self, var, name_ok=None):
        a, b = BASELINE
        name_ok = name_ok or (lambda n: re.search(f"_{var_pattern(var)}_", n))
        files = self.c.tifs(["bioclim", "climatologies"], var)
        hits = [k for k in files
                if f"_{a}_{b}_" in norm(k) and "_ssp" not in norm(k) and name_ok(norm(Path(k).name))]
        return self._pick(hits, f"{var} {a}-{b}")

    def _clim_month_key(self, var, month):
        """1981-2010 monthly climatology of a monthly variable"""
        a, b = BASELINE
        rx = re.compile(rf"_{var_pattern(var)}_{month:02d}_{a}_{b}_")
        hits = [k for k in self.c.tifs(["climatologies"], var)
                if "_ssp" not in norm(k) and rx.search(norm(Path(k).name))]
        return self._pick(hits, f"{var} {month:02d} {a}-{b}")

    def _ssp_keys(self, products, var, period, month=None, name_ok=None):
        a, b = period
        files = self.c.tifs(products, var)
        name_ok = name_ok or (lambda n: re.search(f"_{var_pattern(var)}_", n))
        keys, missing = [], []
        for g in GCMS:
            gtok = norm(g)

            def ok(k):
                n, base = norm(k), norm(Path(k).name)
                if f"_{SSP}_" not in n or gtok not in n or f"_{a}_{b}_" not in n:
                    return False
                if not name_ok(base):
                    return False
                return month is None or f"_{month:02d}_" in base

            k = self._pick([k for k in files if ok(k)], f"{var} {g} {SSP} {a}-{b} m{month}")
            if k:
                keys.append(k)
            else:
                missing.append(g)
        note = f"{a}-{b}: no file for GCM(s) {', '.join(missing)}" if keys and missing else ""
        return keys, note

    # ---- future anchors (one per 30-yr period)
    def _ssp_anchor(self, products, var, month):
        def build(period):
            keys, note = self._ssp_keys(products, var, period, month)
            return ((period_centre(period, month), keys, "median") if keys else None), note
        return build

    def _delta_anchor(self, var, month):
        """monthly projection = 1981-2010 month climatology (+|x) change of the GCM annual mean"""
        layer = MEAN_LAYER[var]
        name_ok = tokens_ok(norm(layer).strip("_").split("_"))
        clim = self._clim_month_key(var, month)
        base = self._baseline_key(layer, name_ok)

        def build(period):
            if not clim or not base:
                what = "monthly 1981-2010 climatology" if not clim else f"{layer} 1981-2010"
                return None, f"delta: no {what} layer"
            keys, note = self._ssp_keys(["bioclim", "climatologies"], layer, period, None, name_ok)
            if not keys:
                return None, f"delta: no {layer} {period[0]}-{period[1]} layer"
            how = {"op": DELTA_OP[var], "clim": clim, "base": base}
            return (period_centre(period, month), keys, how), note
        return build

    def _projection(self, year, month, build, past_anchor, kind):
        """bracket the sampling time with period-centre anchors (see module docstring).
        build(period) -> (anchor or None, note); past_anchor() -> (time, [keys], label) or None"""
        if future_period(year) is None:
            return Source("missing", note="year beyond 2100")
        t = dec_time(year, month)
        centres = [period_centre(p, month) for p in FUTURE_PERIODS]
        nxt = next((i for i, c in enumerate(centres) if c >= t), None)
        if nxt is None:
            last = FUTURE_PERIODS[-1]
            wanted, notes = [len(FUTURE_PERIODS) - 1], [f"after the {last[0]}-{last[1]} centre: held constant"]
        else:
            wanted, notes = ([nxt - 1, nxt] if nxt > 0 else [nxt]), []
        anchors = []
        for i in wanted:
            anchor, note = build(FUTURE_PERIODS[i])
            if anchor:
                anchors.append(anchor)
            elif not note:
                note = f"no {FUTURE_PERIODS[i][0]}-{FUTURE_PERIODS[i][1]} layer"
            if note:
                notes.append(note)
        if not anchors:
            return Source("missing", note="; ".join(notes))
        if nxt == 0:
            past = past_anchor()
            if past:
                anchors.insert(0, (past[0], past[1], "mean"))
                notes.insert(0, f"past anchor: {past[2]}")
            else:
                notes.insert(0, "no past anchor found: projection value used as is")
        return Source(kind, anchors, t, "; ".join(n for n in notes if n))

    # ---- public
    def monthly(self, var: str, year: int, month: int) -> Source:
        mk = ("m", var, year, month)
        if mk in self._memo:
            return self._memo[mk]
        obs = self._observed(var)
        last_full = self._last_complete_year(obs) if obs else None

        def past():
            """same calendar month, N_RECENT years ending at the last complete (12-month) year"""
            end = min(y for y in (last_full, year - 1) if y is not None)
            yrs = sorted(y for (m, y) in obs if m == month and y <= end)[-N_RECENT:]
            if not yrs:
                return None
            return (dec_time(float(np.mean(yrs)), month), [obs[(month, y)] for y in yrs],
                    f"observed {month:02d}/{yrs[0]}-{yrs[-1]} ({len(yrs)} yr; last complete year {last_full})")

        if (month, year) in obs:
            src = Source("timeseries", [(None, [obs[(month, year)]], "mean")])
        elif not obs:
            src = Source("missing", note=f"no monthly files found for {var} in the bucket")
        elif future_period(year) is None:
            src = Source("missing", note="month not in the time series, and no SSP period covers this year")
        else:
            src = self._projection(year, month, self._ssp_anchor(["climatologies"], var, month),
                                   past, f"interp_{SSP}")
            if src.kind == "missing" and var in MEAN_LAYER:
                src = self._projection(year, month, self._delta_anchor(var, month),
                                       past, f"interp_delta_{SSP}")
            if src.kind == "missing":
                p = past()
                if p:
                    src = Source("recent_obs", [(None, p[1], "mean")],
                                 note=f"no {SSP} projection for {var}: mean of {p[2]}")
                else:
                    src.note = f"month not in the time series; {src.note}"
        self._memo[mk] = src
        return src

    def annual(self, var: str, year: int) -> Source:
        mk = ("a", var, year)
        if mk in self._memo:
            return self._memo[mk]
        if year <= BASELINE[1]:
            a, b = BASELINE
            hit = self._baseline_key(var)
            note = "sampling year before 1981, 1981-2010 used" if year < a else ""
            src = (Source(f"baseline_{a}-{b}", [(None, [hit], "mean")], note=note) if hit else
                   Source("missing", note=f"no {a}-{b} layer"))
        else:
            def past():
                hit = self._baseline_key(var)
                return (period_centre(BASELINE), [hit], f"{BASELINE[0]}-{BASELINE[1]} climatology") if hit else None
            src = self._projection(year, None, self._ssp_anchor(["bioclim", "climatologies"], var, None),
                                   past, f"interp_{SSP}")
        self._memo[mk] = src
        return src


def evaluate(src: Source, var: str, ck: str, cache) -> tuple:
    """value for one point; returns (value, n_values, anchor description, converted_to_celsius)"""
    kelvin = False

    def val(k):
        nonlocal kelvin
        v = cache.get(k, ck)
        if v is not None and var in KELVIN_VARS and abs(cache.meta.get(k, {}).get("offset") or 0.0) < 1e-9:
            kelvin = True                  # file stores Kelvin (offset 0): report degC like the bio layers
            v -= 273.15
        return v

    pts, n = [], 0
    for t_anchor, keys, how in src.anchors:
        fut = [v for v in (val(k) for k in keys) if v is not None]
        if isinstance(how, dict):
            c, b = val(how["clim"]), val(how["base"])
            if c is None or b is None:
                fut = []
            elif how["op"] == "add":
                fut = [c + (f - b) for f in fut]
            else:
                fut = [c * f / b for f in fut if b > 0]
        n += len(fut)
        if fut:
            agg = np.mean if how == "mean" else np.median
            pts.append((t_anchor, float(agg(fut))))
    if not pts:
        return np.nan, 0, "", kelvin
    if src.t is None:
        return pts[0][1], n, "", kelvin
    times, vals = zip(*sorted(pts))
    value = float(np.interp(src.t, times, vals))       # clamps outside the anchor range
    desc = " -> ".join(f"{tt:.2f}:{vv:.4g}" for tt, vv in zip(times, vals))
    return value, n, desc, kelvin


# --------------------------------------------------------------------------- cache of sampled values
class ValueCache:
    """key -> {coord_key: value}; appended as JSON lines so a crashed run can resume."""

    def __init__(self, path: Path):
        self.path = path
        self.data: dict[str, dict] = {}
        self.meta: dict[str, dict] = {}
        if path.exists():
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
        with open(self.path, "a") as fh:
            fh.write(json.dumps({"key": key, "values": values, **meta}) + "\n")

    def get(self, key: str, ck: str):
        return self.data.get(key, {}).get(ck)


# --------------------------------------------------------------------------- remote sampling
def sample(url: str, cks: list[str]):
    """read the pixel under each point from a remote GeoTIFF, scale/offset applied"""
    lons = [float(c.split(",")[0]) for c in cks]
    lats = [float(c.split(",")[1]) for c in cks]
    with rasterio.open(f"/vsicurl/{url}") as src:
        # CHELSA global layers are WGS84 lon/lat grids. Checked via the extent rather than
        # the CRS, because CRS lookups need PROJ's proj.db, which some containers lack.
        b = src.bounds
        if not (-181 <= b.left < b.right <= 181 and -91 <= b.bottom < b.top <= 91):
            raise ValueError(f"extent {tuple(round(x, 2) for x in b)} is not a lon/lat grid")
        scale = float(src.scales[0]) if src.scales else 1.0
        offset = float(src.offsets[0]) if src.offsets else 0.0
        vals = []
        for v in src.sample(zip(lons, lats), indexes=1, masked=True):
            if np.ma.getmaskarray(v)[0] or not np.isfinite(float(v[0])):
                vals.append(None)
            else:
                vals.append(float(v[0]) * scale + offset)
    return dict(zip(cks, vals)), {"scale": scale, "offset": offset}


# --------------------------------------------------------------------------- points
def read_points(path: Path):
    df = pd.read_csv(path)
    missing = {"lat", "lon", "sample_date"} - set(df.columns)
    if missing:
        raise ValueError(f"{path} lacks column(s) {sorted(missing)}; required: lat, lon, sample_date")
    dates = pd.to_datetime(df["sample_date"].astype(str).str.strip(), format=DATE_FORMAT, errors="coerce")
    if dates.isna().any():
        bad = df.loc[dates.isna(), "sample_date"].astype(str).tolist()
        raise ValueError(f"{len(bad)} date(s) not in yyyy-mm-dd: {bad[:10]}")
    lat = pd.to_numeric(df["lat"], errors="coerce")
    lon = pd.to_numeric(df["lon"], errors="coerce")
    bad = lat.isna() | lon.isna() | ~lat.between(-90, 90) | ~lon.between(-180, 180)
    if bad.any():
        raise ValueError(f"invalid lat/lon in rows {df.index[bad].tolist()[:10]} (decimal degrees, WGS84)")
    return df, lat.to_numpy(float), lon.to_numpy(float), dates


# --------------------------------------------------------------------------- main
def main() -> int:
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, force=True,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    # GDAL/PROJ chatter (e.g. "Cannot find proj.db") is harmless here: sampling only uses the
    # pixel grid. Real read errors still raise and are logged as FAILED below.
    logging.getLogger("rasterio").setLevel(logging.ERROR)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    PROV_CSV.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    df, lat, lon, dates = read_points(POINTS_CSV)
    log.info("%d points read from %s (%s to %s)", len(df), POINTS_CSV,
             dates.min().strftime("%Y-%m-%d"), dates.max().strftime("%Y-%m-%d"))

    # ---- 1+2: list the bucket and resolve every point x variable -------------
    bucket = Bucket(CHELSA_BUCKET)
    resolver = Resolver(Catalog(bucket))
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

    log.info("values to fill: %s", dict(Counter(src.kind.split("_")[0] for p in plans for src in p.values())))
    log.info("plan: %d unique GeoTIFFs", len(need))
    pd.DataFrame(
        [{"key": k, "url": bucket.url(k), "size_mb": round(bucket.sizes.get(k, 0) / 1e6, 1),
          "n_coords": len(c), "variables": " ".join(sorted(need_vars[k]))} for k, c in sorted(need.items())]
    ).to_csv(PLAN_CSV, index=False)

    # ---- 3: open each file once (remote), sample all its points --------------
    cache = ValueCache(CACHE_JSONL)
    todo = [(k, cache.missing(k, c)) for k, c in sorted(need.items())]
    todo = [(k, c) for k, c in todo if c]
    log.info("%d files already in cache, %d to read", len(need) - len(todo), len(todo))
    failed = {}
    with rasterio.Env(**VSICURL_ENV):
        for i, (key, cks) in enumerate(todo, 1):
            log.info("[%d/%d] %s (%d points)", i, len(todo), Path(key).name, len(cks))
            try:
                vals, meta = sample(bucket.url(key), cks)
            except Exception as e:                      # keep going, report at the end
                log.error("  FAILED %s: %s", key, e)
                failed[key] = str(e)
                continue
            cache.add(key, vals, meta)

    # ---- assemble -----------------------------------------------------------
    out = df.copy()
    cols = {v: np.full(len(df), np.nan) for v in ALL_VARS}
    prov = []
    for i, p in enumerate(plans):
        ck = coord_key(lon[i], lat[i])
        for var, src in p.items():
            cols[var][i], n_ok, anchors, kelvin = evaluate(src, var, ck, cache)
            note = src.note
            if kelvin:
                note = (note + "; " if note else "") + "K -> degC"
            if any(k in failed for k in src.keys):
                note = (note + "; " if note else "") + "read failed"
            prov.append({
                "row": i,
                "site": df["site"].iloc[i] if "site" in df.columns else None,
                "sample_date": df["sample_date"].iloc[i],
                "lat": lat[i], "lon": lon[i],
                "variable": var,
                "source": src.kind,
                "value": cols[var][i],
                "t_sample": round(src.t, 3) if src.t is not None else None,
                "anchors": anchors,           # time:value of each interpolation anchor
                "n_files": len(src.keys),
                "n_values": n_ok,
                "scale_offset": ";".join(f"{cache.meta.get(k, {}).get('scale')}/{cache.meta.get(k, {}).get('offset')}"
                                         for k in src.keys),
                "files": ";".join(Path(k).name for k in src.keys),
                "note": note,
            })
    for v in ALL_VARS:
        out[v] = np.round(cols[v], 4)
    out.to_csv(OUT_CSV, index=False)
    pd.DataFrame(prov).to_csv(PROV_CSV, index=False)
    log.info("wrote %s and %s", OUT_CSV, PROV_CSV)

    empty = [v for v in ALL_VARS if out[v].isna().all()]
    if empty:
        log.warning("variables with no value at any point: %s (see provenance)", " ".join(empty))
    if failed:
        log.error("%d file(s) failed; rerun to retry them (cached values are reused)", len(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())