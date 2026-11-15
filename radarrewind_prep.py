#!/usr/bin/env python3
"""
RadarRewind data prep: NISAR GCOV  ->  small, aligned, map-ready GeoTIFFs.
SI notice: this code was written with help from an AI assistant (Claude by Anthropic) and reviewed/tested by the team.

What it does
  1. Logs in to NASA Earthdata (free account) and finds NISAR L2 GCOV radar
     passes over your area and dates.
  2. Downloads them (or uses .h5 files you already have).
  3. Cuts out your area, converts radar backscatter to decibels (dB),
     and puts every date on the SAME lat/lon grid so dates can be compared
     pixel by pixel.
  4. Writes one GeoTIFF per date:  out/YYYYMMDD_<pol>.tif
     -> drag these files into the RadarRewind web app.
  5. Optional --analyze: writes a change map (PNG + GeoTIFF) and stats.json
     comparing the first and last date.

Install
  pip install earthaccess h5py numpy rasterio pyproj pillow

Examples
  # Nile Delta farmland, Jan-Jun 2026
  python radarrewind_prep.py --bbox 30.95 30.95 31.25 31.15 \
      --start 2026-01-01 --end 2026-06-30 --max 8 --analyze

  # Use .h5 files you already downloaded
  python radarrewind_prep.py --bbox 30.95 30.95 31.25 31.15 --files raw/*.h5
"""
import argparse
import glob
import json
import math
import os
import re
import sys
from collections import defaultdict

import h5py
import numpy as np
import rasterio
from rasterio.transform import Affine, from_origin
from rasterio.warp import Resampling, reproject
from pyproj import Transformer

DEFAULT_SHORT_NAME = "NISAR_L2_GCOV_PROVISIONAL_V1"   # calibrated release
FALLBACK_SHORT_NAME = "NISAR_L2_GCOV_BETA_V1"         # pre-calibration release


# ----------------------------------------------------------------- search/download
def search_and_download(bbox, start, end, short_name, max_granules, raw_dir):
    try:
        import earthaccess
    except ImportError:
        sys.exit("Install earthaccess first:  pip install earthaccess")

    earthaccess.login()  # asks for your Earthdata username/password once
    names = [short_name] if short_name else [DEFAULT_SHORT_NAME, FALLBACK_SHORT_NAME]
    results = []
    for name in names:
        results = earthaccess.search_data(
            short_name=name,
            bounding_box=tuple(bbox),
            temporal=(start, end),
            count=max_granules,
        )
        print(f"[search] {name}: {len(results)} granule(s)")
        if results:
            break
    if not results:
        sys.exit("No NISAR GCOV passes found. Try a bigger area or wider dates.")
    os.makedirs(raw_dir, exist_ok=True)
    files = earthaccess.download(results, local_path=raw_dir)
    return [str(f) for f in files if str(f).endswith(".h5")]


# ----------------------------------------------------------------- reading GCOV
def find_grid_group(h5):
    """Return the HDF5 group that holds xCoordinates/yCoordinates + covariance terms."""
    hits = []

    def visit(name, obj):
        if isinstance(obj, h5py.Group) and "xCoordinates" in obj and "yCoordinates" in obj:
            hits.append(name)

    h5.visititems(visit)
    if not hits:
        raise ValueError("No xCoordinates/yCoordinates grid found - is this a GCOV file?")
    hits.sort(key=lambda n: ("frequencyA" not in n, n))  # prefer frequencyA
    return h5[hits[0]]


def read_epsg(grid):
    if "projection" in grid:
        proj = grid["projection"]
        for key in ("epsg_code", "epsg"):
            if key in proj.attrs:
                return int(np.array(proj.attrs[key]).ravel()[0])
        try:
            return int(np.array(proj[()]).ravel()[0])
        except Exception:
            pass
    raise ValueError("Could not read the EPSG code of the GCOV grid")


def pick_pol(grid, pol):
    wanted = [pol.upper() * 2] if pol else ["HHHH", "VVVV", "HVHV", "VHVH"]
    for term in wanted:
        if term in grid:
            return term
    raise ValueError(f"No covariance term {wanted} in file (has: {list(grid.keys())})")


def acquisition_date(path, h5):
    m = re.search(r"_(\d{8})T(\d{6})_", os.path.basename(path))
    if m:
        return m.group(1)
    for key in ("zeroDopplerStartTime",):
        for name in (f"science/LSAR/identification/{key}",):
            if name in h5:
                txt = h5[name][()]
                txt = txt.decode() if isinstance(txt, bytes) else str(txt)
                return txt[:10].replace("-", "")
    return "unknown"


def read_window(path, bbox, pol):
    """Cut the bbox out of one GCOV file. Returns linear backscatter + its geo-transform."""
    w, s, e, n = bbox
    with h5py.File(path, "r") as h5:
        grid = find_grid_group(h5)
        term = pick_pol(grid, pol)
        epsg = read_epsg(grid)
        x = grid["xCoordinates"][()]
        y = grid["yCoordinates"][()]
        to_proj = Transformer.from_crs(4326, epsg, always_xy=True)
        lons = np.concatenate([np.linspace(w, e, 9), np.full(9, e), np.linspace(e, w, 9), np.full(9, w)])
        lats = np.concatenate([np.full(9, s), np.linspace(s, n, 9), np.full(9, n), np.linspace(n, s, 9)])
        px, py = to_proj.transform(lons, lats)
        dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
        pad = 2
        cols = np.where((x >= px.min() - abs(dx) * pad) & (x <= px.max() + abs(dx) * pad))[0]
        rows = np.where((y >= py.min() - abs(dy) * pad) & (y <= py.max() + abs(dy) * pad))[0]
        if cols.size == 0 or rows.size == 0:
            return None
        c0, c1, r0, r1 = cols[0], cols[-1] + 1, rows[0], rows[-1] + 1
        data = grid[term][r0:r1, c0:c1].astype("float32")
        date = acquisition_date(path, h5)
    data[~np.isfinite(data) | (data <= 0)] = np.nan
    transform = Affine(dx, 0, x[c0] - dx / 2, 0, dy, y[r0] - dy / 2)  # coords are pixel centres
    return dict(data=data, transform=transform, epsg=epsg, date=date, term=term)


# ----------------------------------------------------------------- regrid + write
def to_common_grid(win, bbox, res):
    w, s, e, n = bbox
    width, height = int(math.ceil((e - w) / res)), int(math.ceil((n - s) / res))
    dst = np.full((height, width), np.nan, dtype="float32")
    reproject(
        source=win["data"], destination=dst,
        src_transform=win["transform"], src_crs=f"EPSG:{win['epsg']}", src_nodata=np.nan,
        dst_transform=from_origin(w, n, res, res), dst_crs="EPSG:4326", dst_nodata=np.nan,
        resampling=Resampling.average,  # average in linear power = less speckle
    )
    return dst


def box_smooth(lin, k):
    """NaN-aware k x k mean filter in linear power: reduces radar speckle ("salt and pepper")."""
    if k <= 1:
        return lin
    valid = np.isfinite(lin)
    v = np.where(valid, lin, 0).astype("float64")
    c = valid.astype("float64")
    r = k // 2

    def boxsum(a):
        a = np.pad(a, r, mode="constant")
        s = a.cumsum(0).cumsum(1)
        s = np.pad(s, ((1, 0), (1, 0)))
        return s[k:, k:] - s[:-k, k:] - s[k:, :-k] + s[:-k, :-k]

    num, den = boxsum(v), boxsum(c)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (num / den).astype("float32")
    out[~valid] = np.nan
    return out


def to_db(linear):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (10 * np.log10(linear)).astype("float32")


def write_tif(path, arr, bbox, res, dtype="float32", nodata=np.nan, desc=None):
    w, s, e, n = bbox
    profile = dict(driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1, dtype=dtype,
                   crs="EPSG:4326", transform=from_origin(w, n, res, res), nodata=nodata,
                   compress="deflate")
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)
        if desc:
            dst.update_tags(DESCRIPTION=desc)


# ----------------------------------------------------------------- change analysis
CLASSES = {
    1: ("New open water", (30, 90, 220)),
    2: ("Water gone", (250, 200, 40)),
    3: ("Darker (flooding, bare soil, cleared/burned)", (235, 90, 50)),
    4: ("Brighter (crop growth, construction, wet soil)", (40, 200, 220)),
}


def classify(before_db, after_db, change_db=3.0, water_db=-18.0):
    diff = after_db - before_db
    valid = np.isfinite(diff)
    cls = np.zeros(diff.shape, dtype="uint8")
    was_w, is_w = before_db < water_db, after_db < water_db
    cls[valid & (diff <= -change_db)] = 3
    cls[valid & (diff >= change_db)] = 4
    cls[valid & is_w & ~was_w] = 1
    cls[valid & was_w & ~is_w] = 2
    return cls, valid


def pixel_area_km2(bbox, res):
    lat = (bbox[1] + bbox[3]) / 2
    return (res * 111.32 * math.cos(math.radians(lat))) * (res * 110.57)


def analyze(stack, bbox, res, out_dir, change_db, water_db):
    dates = sorted(stack)
    d0, d1 = dates[0], dates[-1]
    cls, valid = classify(stack[d0], stack[d1], change_db, water_db)
    a = pixel_area_km2(bbox, res)
    stats = {"before": d0, "after": d1, "change_threshold_db": change_db,
             "water_threshold_db": water_db, "valid_area_km2": round(float(valid.sum() * a), 3),
             "classes": {}}
    for k, (name, _) in CLASSES.items():
        stats["classes"][name] = round(float((cls == k).sum() * a), 3)
    write_tif(os.path.join(out_dir, f"change_{d0}_{d1}.tif"), cls, bbox, res, dtype="uint8", nodata=0)
    try:
        from PIL import Image
        rgba = np.zeros(cls.shape + (4,), dtype="uint8")
        for k, (_, rgb) in CLASSES.items():
            rgba[cls == k] = (*rgb, 255)
        Image.fromarray(rgba).save(os.path.join(out_dir, f"change_{d0}_{d1}.png"))
    except ImportError:
        pass
    with open(os.path.join(out_dir, "stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2))


# ----------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bbox", nargs=4, type=float, required=True, metavar=("W", "S", "E", "N"))
    p.add_argument("--start", default="2025-10-01")
    p.add_argument("--end", default="2026-11-15")
    p.add_argument("--short-name", default=None, help=f"default: {DEFAULT_SHORT_NAME}, then {FALLBACK_SHORT_NAME}")
    p.add_argument("--pol", default=None, help="HH, VV, HV ... (default: first available co-pol)")
    p.add_argument("--res", type=float, default=0.0002, help="output pixel size in degrees (0.0002 ~ 20 m)")
    p.add_argument("--max", type=int, default=8, help="max granules to download")
    p.add_argument("--smooth", type=int, default=5, help="speckle filter window in pixels (odd; 1 = off)")
    p.add_argument("--files", nargs="*", help="use these local GCOV .h5 files instead of downloading")
    p.add_argument("--out", default="out")
    p.add_argument("--analyze", action="store_true", help="also write change map + stats.json (first vs last date)")
    p.add_argument("--change-db", type=float, default=3.0)
    p.add_argument("--water-db", type=float, default=-18.0)
    a = p.parse_args()

    bbox = a.bbox
    if not (bbox[0] < bbox[2] and bbox[1] < bbox[3]):
        sys.exit("--bbox must be W S E N (west < east, south < north)")
    cells = ((bbox[2] - bbox[0]) / a.res) * ((bbox[3] - bbox[1]) / a.res)
    if cells > 12e6:
        sys.exit("Area too big for the browser at this --res. Use a smaller bbox or a larger --res.")
    os.makedirs(a.out, exist_ok=True)

    files = []
    for f in a.files or []:
        files.extend(glob.glob(f))
    if not files:
        files = search_and_download(bbox, a.start, a.end, a.short_name, a.max, os.path.join(a.out, "raw"))

    per_date = defaultdict(list)
    term_used = None
    for f in sorted(files):
        try:
            win = read_window(f, bbox, a.pol)
        except Exception as ex:  # keep going if one file is odd
            print(f"[skip] {os.path.basename(f)}: {ex}")
            continue
        if win is None:
            print(f"[skip] {os.path.basename(f)}: does not cover the bbox")
            continue
        term_used = win["term"]
        per_date[win["date"]].append(to_common_grid(win, bbox, a.res))
        print(f"[read] {os.path.basename(f)} -> {win['date']} {win['term']}")

    if not per_date:
        sys.exit("Nothing covered your area.")

    stack = {}
    pol = (term_used or "HHHH")[:2]
    for date, layers in sorted(per_date.items()):
        with np.errstate(invalid="ignore"):
            lin = np.nanmean(np.stack(layers), axis=0)  # mosaic frames from the same date
        db = to_db(box_smooth(lin, a.smooth))
        stack[date] = db
        path = os.path.join(a.out, f"{date}_{pol}.tif")
        write_tif(path, db, bbox, a.res, desc=f"NISAR GCOV gamma0 {pol} (dB) {date}")
        cover = float(np.isfinite(db).mean() * 100)
        print(f"[write] {path}  ({cover:.0f}% of area covered)")

    if a.analyze and len(stack) >= 2:
        analyze(stack, bbox, a.res, a.out, a.change_db, a.water_db)
    print(f"\nDone. Drag the .tif files in '{a.out}/' into the RadarRewind app.")


if __name__ == "__main__":
    main()
