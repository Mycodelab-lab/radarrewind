# RadarRewind: see what changed, through clouds

**NASA Space Apps Challenge 2026 · Challenge: Dancing with the SARs**

RadarRewind is an interactive web app that uses **NASA-ISRO NISAR** L-band radar to track and visualize surface change anywhere on Earth: floods, farming cycles, land reclamation, construction and clearing.

**Live app:** `https://<your-github-username>.github.io/radarrewind/`

## What it does
1. **Find passes:** searches NASA's live Common Metadata Repository (CMR) for real NISAR acquisitions over the map view. It draws each footprint, shows a timeline of passes and the typical revisit time, and links to every file.
2. **See change:** load one radar map per date and pick a "before" and an "after". Every pixel is labelled *new open water*, *water gone*, *darker* or *brighter*, with areas in km². You can press ▶ Play to animate the dates, and click any spot to see its radar history with a plain-language explanation.
3. **Prep script** (`radarrewind_prep.py`): downloads NISAR L2 GCOV files with `earthaccess`, cuts out your area, reduces speckle, converts to dB, puts every date on one lat/lon grid, and writes GeoTIFFs plus an optional change map and `stats.json`.

## Run it
```bash
pip install earthaccess h5py numpy rasterio pyproj pillow
python radarrewind_prep.py --bbox 30.95 30.95 31.25 31.15 --start 2026-01-01 --end 2026-06-30 --max 8 --analyze
```
Open `index.html` (or the live link) → tab **2 · See change** → drop the `.tif` files from `out/`.

## Method
- Input: NISAR L2 GCOV γ⁰ backscatter (HH, or VV if HH is missing), linear power.
- Cut to the area of interest, then reprojected with an *average* of linear power onto a shared EPSG:4326 grid (default 0.0002° ≈ 20 m).
- 5×5 NaN-aware box filter for speckle, then converted to dB.
- Change classes (thresholds adjustable in the app):
  - new water: before ≥ T_w and after < T_w (default T_w = −18 dB)
  - water gone: the reverse
  - darker: after − before ≤ −T_c, brighter: ≥ +T_c (default T_c = 3 dB)
- Area = pixel count × pixel area at the area's mean latitude.

## Tested
- The prep script was run on synthetic GCOV-structured HDF5 files with a known 120 km² flood. It recovered 118.7 km² of new water inside the area of interest. The web app computes the same number from the same files.
- The **simulated sample** in the app is made-up data for learning the controls and is labelled as such. It is not a NISAR measurement.

## Limits
- Radar "dark" also happens on smooth surfaces such as tarmac and calm sand, and wind on water can make water look brighter. That's why the thresholds are adjustable and results should be read alongside the satellite basemap.
- GCOV Beta products are pre-calibration. Use the Provisional collection where it's available.

## Use of Super Intelligence (SI)
Code and text in this repository were written with help from an AI assistant (Claude by Anthropic). The team reviewed, tested and validated all of it. No AI-generated images, audio or video are used.

## Data & tools
NISAR L2 GCOV (ASF DAAC, NASA Earthdata), NASA CMR search API, earthaccess, h5py, rasterio, pyproj, Leaflet, geotiff.js, Esri World Imagery and OpenStreetMap basemaps.
