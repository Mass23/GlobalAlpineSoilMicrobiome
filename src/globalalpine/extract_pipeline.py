"""
Unified extraction utilities for SoilGrids, CHELSA, Copernicus DEM, landcover, MODIS NDVI, and snow-related rasters.

Design goals:
- SoilGrids: per-variable per-depth-band COGs -> depth-interpolated value when point has depth, else 0-30 cm thickness-weighted average.
- CHELSA: sample supplied rasters (BIO1..BIO19, monthly, freezing vars, precipitation, temperature)
- DEM: sample elevation and compute local slope, aspect, TPI and hillshade from small window around pixel
- Landcover/NDVI/Snow: simple per-point sampling from raster COGs

Requires: rasterio, numpy, pandas, tqdm
"""
from typing import Dict, List, Tuple, Optional
import math
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from tqdm import tqdm

# Typical SoilGrids depth bands and bounds (top, bottom) in cm
DEFAULT_BANDS = [
    ("0-5", 0.0, 5.0),
    ("5-15", 5.0, 15.0),
    ("15-30", 15.0, 30.0),
    ("30-60", 30.0, 60.0),
    ("60-100", 60.0, 100.0),
    ("100-200", 100.0, 200.0),
]


def sample_band_rasters(coords: List[Tuple[float, float]], band_paths: List[str]) -> np.ndarray:
    """
    Sample each raster in band_paths at coords and return (n_points, n_bands) array.
    band_paths order must match depth order.
    Supports local files and remote COGs.
    """
    n = len(coords)
    m = len(band_paths)
    out = np.full((n, m), np.nan, dtype=float)
    for j, path in enumerate(band_paths):
        with rasterio.open(path) as src:
            vals = list(src.sample(coords))
            # sample returns array-like per-band; take first band
            arr = np.array([v[0] if hasattr(v, '__iter__') else v for v in vals], dtype=float)
            arr[~np.isfinite(arr)] = np.nan
            out[:, j] = arr
    return out


def compute_depth_interp_single(band_values: np.ndarray, mids: np.ndarray, depth_cm: float) -> float:
    """
    Linear interpolation of band_values at depth_cm along band midpoints.
    If depth outside range, clamp to nearest band value.
    """
    if np.isnan(depth_cm):
        return np.nan
    if np.all(np.isnan(band_values)):
        return np.nan
    valid = ~np.isnan(band_values)
    if valid.sum() == 0:
        return np.nan
    vals = band_values[valid]
    ms = mids[valid]
    if ms.size == 1:
        return float(vals[0])
    if depth_cm <= ms.min():
        return float(vals[ms.argmin()])
    if depth_cm >= ms.max():
        return float(vals[ms.argmax()])
    idx = np.searchsorted(ms, depth_cm)
    i0 = idx - 1
    i1 = idx
    m0, m1 = ms[i0], ms[i1]
    v0, v1 = vals[i0], vals[i1]
    t = (depth_cm - m0) / (m1 - m0) if (m1 - m0) != 0 else 0.0
    return float(v0 * (1 - t) + v1 * t)


def thickness_weighted_avg_0_30(band_values: np.ndarray, tops: np.ndarray, bottoms: np.ndarray, target_top: float = 0.0, target_bottom: float = 30.0) -> float:
    """
    Compute thickness-weighted average of band_values over [target_top, target_bottom).
    Treat band_values as per-band average; overlap is used as weight.
    """
    ov_top = np.maximum(tops, target_top)
    ov_bottom = np.minimum(bottoms, target_bottom)
    overlap = np.maximum(0.0, ov_bottom - ov_top)
    valid = (~np.isnan(band_values)) & (overlap > 0)
    if not np.any(valid):
        return np.nan
    weights = overlap[valid]
    vals = band_values[valid]
    return float(np.sum(vals * weights) / np.sum(weights))


def extract_soilgrids_for_points(
    df: pd.DataFrame,
    coords_cols: Tuple[str, str] = ("lon", "lat"),
    depth_col: Optional[str] = None,
    variables_to_band_paths: Dict[str, Dict[str, str]] = None,
    bands_spec: List[Tuple[str, float, float]] = DEFAULT_BANDS,
    chunk_size: int = 5000,
) -> pd.DataFrame:
    """
    df: points with lon/lat and optional depth (cm) column.
    variables_to_band_paths: mapping var_name -> mapping band_label -> local COG path
    Returns original df plus columns:
      soilgrid_depth_<var>   (interpolated at df[depth_col] or NaN if no depth)
      soilgrid_nodepth_<var> (0-30 cm weighted average)
    """
    if variables_to_band_paths is None:
        raise ValueError("variables_to_band_paths required")

    coords = list(df[[coords_cols[0], coords_cols[1]]].itertuples(index=False, name=None))
    n = len(coords)
    depths = df[depth_col].to_numpy(dtype=float) if depth_col else np.array([np.nan] * n)

    band_labels = [b[0] for b in bands_spec]
    tops = np.array([b[1] for b in bands_spec], dtype=float)
    bottoms = np.array([b[2] for b in bands_spec], dtype=float)
    mids = (tops + bottoms) / 2.0

    result = df.copy()

    for var, band_map in variables_to_band_paths.items():
        band_paths = []
        for label in band_labels:
            if label not in band_map:
                raise ValueError(f"Variable {var} missing band {label} path")
            band_paths.append(band_map[label])

        depth_values = np.full(n, np.nan, dtype=float)
        nodepth_values = np.full(n, np.nan, dtype=float)

        for start in tqdm(range(0, n, chunk_size), desc=f"sampling {var}"):
            end = min(n, start + chunk_size)
            coords_chunk = coords[start:end]
            sampled = sample_band_rasters(coords_chunk, band_paths)

            for i in range(sampled.shape[0]):
                vals = sampled[i, :]
                if depth_col:
                    d = depths[start + i]
                    depth_values[start + i] = compute_depth_interp_single(vals, mids, d)
                nodepth_values[start + i] = thickness_weighted_avg_0_30(vals, tops, bottoms)

        result[f"soilgrid_depth_{var}"] = depth_values
        result[f"soilgrid_nodepth_{var}"] = nodepth_values

    return result


# ---- CHELSA / simple single-band rasters -------------------------------------------------

def sample_singleband_rasters(df: pd.DataFrame, coords_cols: Tuple[str, str], var_map: Dict[str, str], chunk_size: int = 5000) -> pd.DataFrame:
    """
    var_map: mapping var_name -> raster path (COG)
    Produces columns named chelsa_<var> for each entry.
    """
    coords = list(df[[coords_cols[0], coords_cols[1]]].itertuples(index=False, name=None))
    n = len(coords)
    result = df.copy()
    for var, path in var_map.items():
        vals = np.full(n, np.nan, dtype=float)
        for start in tqdm(range(0, n, chunk_size), desc=f"sampling {var}"):
            end = min(n, start + chunk_size)
            chunk = coords[start:end]
            with rasterio.open(path) as src:
                svals = list(src.sample(chunk))
                arr = np.array([v[0] if hasattr(v, '__iter__') else v for v in svals], dtype=float)
                arr[~np.isfinite(arr)] = np.nan
                vals[start:end] = arr
        result[f"chelsa_{var}"] = vals
    return result


# ---- DEM derived metrics ----------------------------------------------------------------

def sample_dem_and_derivs(dem_path: str, coords: List[Tuple[float, float]], window: int = 3) -> Dict[str, np.ndarray]:
    """
    For each coord, returns elevation, slope_deg, aspect_deg, tpi, hillshade.
    window must be odd (e.g., 3,5,7). Reads a small window around the pixel and computes gradients.
    """
    assert window % 2 == 1
    half = window // 2
    n = len(coords)
    elevs = np.full(n, np.nan, dtype=float)
    slopes = np.full(n, np.nan, dtype=float)
    aspects = np.full(n, np.nan, dtype=float)
    tpis = np.full(n, np.nan, dtype=float)
    hills = np.full(n, np.nan, dtype=float)

    with rasterio.open(dem_path) as src:
        band1 = 1
        for i, (lon, lat) in enumerate(tqdm(coords, desc="sampling DEM")):
            try:
                col, row = src.index(lon, lat)
            except Exception:
                # outside bounds
                continue
            # window bounds
            col_off = col - half
            row_off = row - half
            win = Window(col_off, row_off, window, window)
            try:
                arr = src.read(band1, window=win, boundless=True, fill_value=np.nan)
            except Exception:
                arr = None
            if arr is None or np.all(np.isnan(arr)):
                continue
            # central elevation
            r0 = half
            c0 = half
            cen = float(arr[r0, c0])
            elevs[i] = cen
            # compute TPI: center - mean(neighbors)
            nbrs = np.delete(arr.flatten(), window * r0 + c0)
            tpi = cen - np.nanmean(nbrs) if np.any(np.isfinite(nbrs)) else np.nan
            tpis[i] = float(tpi) if np.isfinite(tpi) else np.nan
            # compute gradients using numpy.gradient; need pixel sizes
            transform = src.transform
            xres = transform.a
            yres = abs(transform.e)
            # Replace nan with nearest valid? Use nan-aware gradient by filling small gaps
            if np.isnan(arr).any():
                # if too many nans, skip
                if np.sum(np.isfinite(arr)) < (window * window) // 2:
                    slopes[i] = np.nan
                    aspects[i] = np.nan
                    hills[i] = np.nan
                    continue
                # fill nans with local mean
                arr = np.where(np.isfinite(arr), arr, np.nanmean(arr))
            # gradient: arr.shape -> (rows (y), cols (x))
            # dz/dy, dz/dx
            gy, gx = np.gradient(arr, yres, xres)
            # take center gradients
            dzdx = gx[r0, c0]
            dzdy = gy[r0, c0]
            slope_rad = math.atan(math.hypot(dzdx, dzdy))
            slope_deg = math.degrees(slope_rad)
            slopes[i] = float(slope_deg)
            # aspect: degrees clockwise from north
            aspect_rad = math.atan2(dzdy, -dzdx)
            aspect_deg = math.degrees(aspect_rad)
            if aspect_deg < 0:
                aspect_deg = 90.0 - aspect_deg
            aspects[i] = float(aspect_deg % 360.0)
            # hillshade: sun azimuth 315, altitude 45
            az = math.radians(315.0)
            alt = math.radians(45.0)
            hs = (math.cos(alt) * math.cos(slope_rad) + math.sin(alt) * math.sin(slope_rad) * math.cos(az - math.radians(aspects[i])))
            hills[i] = float(max(0.0, hs))

    return {
        "elevation": elevs,
        "slope_deg": slopes,
        "aspect_deg": aspects,
        "tpi": tpis,
        "hillshade": hills,
    }


def sample_dem_for_df(df: pd.DataFrame, dem_path: str, coords_cols: Tuple[str, str] = ("lon", "lat"), window: int = 3) -> pd.DataFrame:
    coords = list(df[[coords_cols[0], coords_cols[1]]].itertuples(index=False, name=None))
    out = sample_dem_and_derivs(dem_path, coords, window=window)
    res = df.copy()
    res["dem_elevation"] = out["elevation"]
    res["dem_slope_deg"] = out["slope_deg"]
    res["dem_aspect_deg"] = out["aspect_deg"]
    res["dem_tpi"] = out["tpi"]
    res["dem_hillshade"] = out["hillshade"]
    return res


# ---- Generic sampler for single-band rasters (landcover, NDVI, snow cover) ----------------

def sample_single_rasters_to_df(df: pd.DataFrame, var_map: Dict[str, str], coords_cols: Tuple[str, str] = ("lon", "lat"), prefix: str = "misc", chunk_size: int = 5000) -> pd.DataFrame:
    coords = list(df[[coords_cols[0], coords_cols[1]]].itertuples(index=False, name=None))
    n = len(coords)
    res = df.copy()
    for var, path in var_map.items():
        vals = np.full(n, np.nan, dtype=float)
        for start in tqdm(range(0, n, chunk_size), desc=f"sampling {var}"):
            end = min(n, start + chunk_size)
            chunk = coords[start:end]
            with rasterio.open(path) as src:
                svals = list(src.sample(chunk))
                arr = np.array([v[0] if hasattr(v, '__iter__') else v for v in svals], dtype=float)
                arr[~np.isfinite(arr)] = np.nan
                vals[start:end] = arr
        col = f"{prefix}_{var}"
        res[col] = vals
    return res


# ---- Top-level orchestrator ---------------------------------------------------------------

def extract_all(
    df: pd.DataFrame,
    soil_vars_map: Dict[str, Dict[str, str]],
    chelsa_map: Dict[str, str],
    dem_path: Optional[str],
    landcover_map: Optional[Dict[str, str]] = None,
    ndvi_map: Optional[Dict[str, str]] = None,
    snow_map: Optional[Dict[str, str]] = None,
    coords_cols: Tuple[str, str] = ("lon", "lat"),
    depth_col: Optional[str] = "depth_cm",
    chunk_size: int = 5000,
    dem_window: int = 3,
) -> pd.DataFrame:
    """
    Run all extractions and return augmented DataFrame.

    - soil_vars_map: variable -> band_label -> COG path (bands must match DEFAULT_BANDS or supplied bands)
    - chelsa_map: variable -> COG path (e.g. BIO1..BIO19, freezing vars, monthly)
    - dem_path: path to DEM COG (Copernicus) or None
    - landcover_map: variable -> path (e.g. LC class)
    - ndvi_map: variable -> path (e.g. MODIS long-term mean)
    - snow_map: variable -> path (e.g. snow frequency)
    """
    out = df.copy()
    # SoilGrids
    if soil_vars_map:
        out = extract_soilgrids_for_points(out, coords_cols=coords_cols, depth_col=depth_col, variables_to_band_paths=soil_vars_map, chunk_size=chunk_size)
    # CHELSA
    if chelsa_map:
        out = sample_singleband_rasters(out, coords_cols, chelsa_map, chunk_size=chunk_size)
    # DEM
    if dem_path:
        out = sample_dem_for_df(out, dem_path, coords_cols=coords_cols, window=dem_window)
    # Landcover
    if landcover_map:
        out = sample_single_rasters_to_df(out, landcover_map, coords_cols=coords_cols, prefix="lc", chunk_size=chunk_size)
    # NDVI
    if ndvi_map:
        out = sample_single_rasters_to_df(out, ndvi_map, coords_cols=coords_cols, prefix="ndvi", chunk_size=chunk_size)
    # Snow
    if snow_map:
        out = sample_single_rasters_to_df(out, snow_map, coords_cols=coords_cols, prefix="snow", chunk_size=chunk_size)

    return out


# Example usage note (not executed):
# df = pd.read_csv('points.csv')  # must have lon, lat, optional depth_cm
# out = extract_all(df, soil_vars_map=..., chelsa_map=..., dem_path='copernicus_dem.tif', landcover_map={'lc': 'lc.tif'}, ndvi_map={'mean': 'ndvi_mean.tif'}, snow_map={'freq': 'snow_freq.tif'})
# out.to_parquet('sampled.parquet', index=False)
