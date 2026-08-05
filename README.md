# NMME Zarr Archive

Local Zarr archives of monthly [NMME](https://www.cpc.ncep.noaa.gov/products/NMME/)
(North American Multi-Model Ensemble) forecast output, sourced from the
[IRIDL OPeNDAP server](https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/).
Currently archives SST and 2-m air temperature (tref) for seven models spanning
the full NMME hindcast + real-time forecast period.

Storing the data locally in Zarr avoids repeated IRIDL requests and enables fast
offline analysis of forecast skill, ensemble spread, and climatologies.

## Models

| Model | Group name |
|-------|------------|
| CanSIPS-IC4 / CanESM5 | `CanSIPS-IC4-CanESM5` |
| CanSIPS-IC4 / GEM5.2-NEMO | `CanSIPS-IC4-GEM52NEMO` |
| GFDL-SPEAR | `GFDL-SPEAR` |
| COLA-RSMAS / CCSM4 | `COLA-RSMAS-CCSM4` |
| COLA-RSMAS / CESM1 | `COLA-RSMAS-CESM1` |
| NASA GEOSS2S | `NASA-GEOSS2S` |
| NCEP CFSv2 | `NCEP-CFSv2` |

## Data Layout

Both stores live under `data/` (not committed to git — build locally).

**`data/nmme_sst.zarr`** — sea surface temperature
- Variable: `sst(S, M, L, Y, X)` float32 °C
- One Zarr group per model (names above)
- `S` = start month ("months since 1960-01-01", 360-day calendar)
- `M` = ensemble member; `L` = lead (months); `Y`/`X` = lat/lon
- Coordinate `T(S, L)` = valid time; grid: global 1° (360 × 181)

**`data/nmme_tref.zarr`** — 2-m air temperature
- Variable: `tref(S, M, L, Y, X)` float32 **Kelvin** (native IRIDL units)
- Same groups, coordinates, and grid as the SST store

## Environment

```bash
# Create the environment (pangeo-local with zarr 3.x)
mamba env create -f environment.yml   # if provided, else install manually

# Run any script
mamba run -n pangeo-local python code/<script>.py
```

Key packages: `xarray`, `zarr` (v3), `numpy`, `netCDF4`, `numcodecs`, `matplotlib`.
`pangeo-local` currently has `zarr 3.2.1` + `xarray 2026.4.0` — this pairing matters, see gotcha below.

**Gotcha:** `xarray`'s zarr backend must be new enough for whatever `zarr` version
is installed. Older `xarray` reading/writing a zarr-format-3 store under
`zarr>=3.2` fails with `AttributeError: 'Float32' object has no attribute
'value'` — zarr 3.2 changed how array dtype metadata is represented
internally, and older xarray releases assumed the previous representation.
This is a local environment/dependency mismatch, not a store or IRIDL
problem. It bit the `pangeo-2025` env (`xarray 2025.3.1` + `zarr 3.2.1`,
from its `zarr>=3` unpinned dependency), which is why this project now uses
`pangeo-local` instead.

## Usage

### Build archives from scratch

```bash
# SST — all 7 models (~hours)
mamba run -n pangeo-local python code/build_archive.py

# tref — preflight one model first to verify URLs and dimensions
mamba run -n pangeo-local python code/build_archive.py \
    --var tref --models NASA-GEOSS2S --block-size 1
mamba run -n pangeo-local python code/build_archive.py --var tref

# Resume a single model (safe to re-run; skips completed starts)
mamba run -n pangeo-local python code/build_archive.py --models NASA-GEOSS2S
```

### Monthly update

```bash
mamba run -n pangeo-local python code/update_archive.py           # both stores (sst + tref)
mamba run -n pangeo-local python code/update_archive.py --var tref  # tref only

# Update a single model
mamba run -n pangeo-local python code/update_archive.py --models NCEP-CFSv2

# Re-fetch last N starts (e.g. if members were incomplete at build time)
mamba run -n pangeo-local python code/update_archive.py --recheck-n 5

# Bypass a stale IRIDL Squid cache entry (see Data Quality Notes)
mamba run -n pangeo-local python code/update_archive.py --var sst --models NASA-GEOSS2S --bust-cache

# Backfill starts that were empty upstream at fetch time, once they appear
# (--recheck-n overwrites the last N starts in place; count back far enough
#  to cover the gap)
mamba run -n pangeo-local python code/update_archive.py --var sst --models GFDL-SPEAR --recheck-n 14
```

If `Remote S` in the log is stuck below the start count shown on the model's
IRIDL page, the update is reading a stale cache entry and will report
`New S: 0` with nothing to append. Re-run the affected model with
`--bust-cache`.

### Sanity check

```bash
mamba run -n pangeo-local python code/sanity_check.py --recent-only  # fast
mamba run -n pangeo-local python code/sanity_check.py                 # full
```

Figures are written to `plots/sanity/`.

## Opening the Stores

```python
import xarray as xr

# Single model (decode_times=False avoids calendar='360' issue)
ds_sst  = xr.open_zarr('data/nmme_sst.zarr',  group='NASA-GEOSS2S', decode_times=False)
ds_tref = xr.open_zarr('data/nmme_tref.zarr', group='NASA-GEOSS2S', decode_times=False)

# All models via DataTree
dt_sst  = xr.open_datatree('data/nmme_sst.zarr',  engine='zarr', decode_times=False)
dt_tref = xr.open_datatree('data/nmme_tref.zarr', engine='zarr', decode_times=False)
```

See `code/sanity_check.py` for worked examples of efficiently computing indices
(Niño 3.4, global mean temperature) from the stores.

## Code Structure

| File | Description |
|------|-------------|
| `code/nmme_models.py` | Model registry; `VARIABLES` table; `resolve_model(group, var)` |
| `code/iridl_io.py` | OPeNDAP helpers: `fetch_data_block`, retry, DDS dim probing |
| `code/build_archive.py` | Initial archive creation; resumes via `_filled` sentinel |
| `code/update_archive.py` | Monthly incremental update; QAs touched starts for corrupt members (WARNING-only) and re-consolidates store metadata |
| `code/sanity_check.py` | Visual sanity check and worked example of reading the stores |

Behavioral specs for each script live in `specs/`.

## Data Quality Notes

- `build_archive.py` logs a WARNING for any forecast start with all-zero or constant data. Check warnings before running analysis.
- IRIDL serves data through a Squid proxy that caches on exact URL text, producing two distinct symptoms:
  - **Zeros / constant fields** in a rebuild — increment `_TREF_HIND_BUST` / `_TREF_FCST_BUST` (or `_SST_HIND_BUST` / `_SST_FCST_BUST`) in `code/nmme_models.py` and re-run `build_archive.py`.
  - **A stale, short `S` axis** during a monthly update — `Remote S` sits below the start count on the IRIDL page, `New S: 0`, and the newest forecast start is silently never appended. Re-run with `--bust-cache`; no code constant change needed. Confirmed on NASA-GEOSS2S sst 2026-08-05: the cached `.FORECAST/.MONTHLY/.sst` URL reported `S = 114` (through 1 Jul 2026) while a cache-busted request reported `S = 115` (through 1 Aug 2026).
  - Note the static `_*_BUST` constants bust the cache only once — Squid then caches that exact busted URL too. `--bust-cache` uses a fresh per-run timestamp token and is the reliable option.
- NCEP-CFSv2: a real data gap around 2010–2011 is visible in tref.
- GFDL-SPEAR: some recent forecast starts may be missing at build time. GFDL posts SPEAR only periodically, so IRIDL can expose an `S` value whose data is not yet there — the missing data is genuinely missing upstream, not a fetch or cache failure. (Distinguish from the stale-`S`-axis symptom above: there, a cache-busted `.dds` request reports *more* starts than the cached one; here, cached and busted agree and the start is simply empty.)
- **Empty starts are written as NaN placeholders, and backfill is manual.** `S` is an append dimension, so a start skipped at write time can never be inserted in order later — the placeholder has to be written when the `S` value first appears. There is deliberately no automatic backfill pass: when the data shows up upstream, re-fetch it by hand with `--models <MODEL> --recheck-n <N>`, choosing `N` large enough to reach back over the gap. `recheck_tail()` overwrites the last `N` starts in place unconditionally, so this fills the NaNs without touching the `S` axis. Known open gaps, GFDL-SPEAR sst as of 2026-08-05 — 2025-06, 2025-08, 2025-10, 2025-11, 2025-12, 2026-02, 2026-03, 2026-05, 2026-06 (all-NaN across every member and lead; 418 of 427 local starts are valid). Note the intermittent pattern: SPEAR has posted roughly every second or third month since mid-2025, and none of these has filled in over the following year, so treat them as likely permanent rather than pending. A backfill reaching the earliest of them needs `--recheck-n 14` (2025-06 is 13 starts back from 2026-07).
- Downstream, `nmme_enso`'s `config.load_nino34_ssta` prints any interior all-NaN starts in its load banner, so a gap is visible at analysis time without reading update logs. Analysis scripts report these gaps; they never mask them.
