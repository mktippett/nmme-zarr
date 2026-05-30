# CLAUDE.md

Project type: research

## Project Overview

Local Zarr archives of monthly NMME model output (currently SST and 2-m air
temperature), sourced from the IRIDL OPeNDAP server. Enables offline analysis
of forecast skill, ensemble spread, and climatologies without repeated IRIDL
requests. Build/update scripts accept `--var {sst,tref}` to select the archive.

## Environment

- Python via mamba (environment: `pangeo-2025`)
- Run scripts with: `mamba run -n pangeo-2025 python code/<script>.py`
- Key packages: xarray, zarr (v3), numpy, netCDF4, numcodecs, matplotlib

## Data

**SST store:** `data/nmme_sst.zarr`
- One group per model: CanSIPS-IC4-CanESM5, CanSIPS-IC4-GEM52NEMO, GFDL-SPEAR,
  COLA-RSMAS-CCSM4, COLA-RSMAS-CESM1, NASA-GEOSS2S, NCEP-CFSv2
- Variable `sst(S, M, L, Y, X)` float32 °C; coord `T(S, L)` = valid time
- S in "months since 1960-01-01" (360-day calendar); L in months
- Grid: global 1° (360 × 181)

**tref store:** `data/nmme_tref.zarr`
- Same seven groups and coordinate conventions as the SST store
- Variable `tref(S, M, L, Y, X)` float32 **Kelvin** (native IRIDL units, no conversion)
- GFDL-SPEAR tref URL uses plain `.tref` (confirmed; SST is the exception with `_regridded`)

**Observations:** ERSSTv5 Niño 3.4 read live from IRIDL (not archived locally).

## Common Commands

```bash
# Build SST archive (all 7 models, ~hours)
mamba run -n pangeo-2025 python code/build_archive.py

# Build tref archive — preflight first to verify URLs and dims
mamba run -n pangeo-2025 python code/build_archive.py \
    --var tref --models NASA-GEOSS2S --block-size 1  # preflight
mamba run -n pangeo-2025 python code/build_archive.py --var tref

# Build or resume a single model
mamba run -n pangeo-2025 python code/build_archive.py --models NASA-GEOSS2S

# Monthly update (append new starts, re-fetch last 2)
mamba run -n pangeo-2025 python code/update_archive.py          # sst
mamba run -n pangeo-2025 python code/update_archive.py --var tref

# Update a single model
mamba run -n pangeo-2025 python code/update_archive.py --models NCEP-CFSv2

# Re-fetch last N starts (e.g. if members were incomplete at build time)
mamba run -n pangeo-2025 python code/update_archive.py --recheck-n 5

# Sanity check (visual verification + worked example of index computation)
mamba run -n pangeo-2025 python code/sanity_check.py --recent-only  # fast: recent forecast only
mamba run -n pangeo-2025 python code/sanity_check.py                 # full: timeseries + recent
```

> **Analysis moved**: ENSO/t2m analysis (`skill_n34.py`, `global_temperature_enso_nmme.py`, etc.)
> has been relocated to `~/claude/enso-t2m`. This project is now data-archive only.

## Scripts

| Script | Description |
|--------|-------------|
| `code/nmme_models.py` | Model registry; `VARIABLES` table; `resolve_model(group, var)` |
| `code/iridl_io.py` | OPeNDAP helpers: `fetch_data_block`, retry, DDS dim probing |
| `code/build_archive.py` | Initial archive creation (`--var sst\|tref`); resumes via `_filled` sentinel |
| `code/update_archive.py` | Monthly incremental update (`--var sst\|tref`) |
| `code/sanity_check.py` | Visual sanity check; also a worked example of reading stores and computing indices |

## Spec Files

Behavioral specifications live in `specs/` (one per script). Update the
corresponding spec the same session as any code change and add a row to its
Synchronization Log.

| Spec | Covers |
|------|--------|
| `specs/build_archive.md` | `build_archive.py` — initial archive creation |
| `specs/update_archive.md` | `update_archive.py` — monthly update |
| `specs/nmme_tref.md` | tref archive layout, URLs, build commands |
| `specs/sanity_check.md` | `sanity_check.py` — index helpers, figure layout, coordinate notes |

## Opening the Stores

```python
import xarray as xr

# Single model, single variable (decode_times=False avoids calendar='360' issue)
ds_sst  = xr.open_zarr('data/nmme_sst.zarr',  group='NASA-GEOSS2S', decode_times=False)
ds_tref = xr.open_zarr('data/nmme_tref.zarr', group='NASA-GEOSS2S', decode_times=False)

# All models via DataTree
dt_sst  = xr.open_datatree('data/nmme_sst.zarr',  engine='zarr', decode_times=False)
dt_tref = xr.open_datatree('data/nmme_tref.zarr', engine='zarr', decode_times=False)
```

## Data Quality

The zarr stores are the source of truth. When downstream analysis reveals bad
values (e.g., implausible constants, fills), **fix the source, not the analysis**.

- Do not add masks or guards in analysis scripts to work around known store defects.
  (`tref.where(tref != 0)` in an analysis script is hiding a problem, not solving it.)
- Instead: identify the root cause in `build_archive.py` / `update_archive.py`,
  add QA to catch it on future builds, and rebuild the affected groups.
- `build_archive.py` logs WARNING for any start with all-zero or constant data.
  Check those warnings before running analysis.

Known issues — **resolved 2026-05-28**:
Previous tref store had all-zero starts due to IRIDL Squid cache returning zeros.
Store was rebuilt clean on 2026-05-27/28 (all 7 models, no warnings).

Expected NaN (real data gaps, not defects):
- NCEP-CFSv2: gap around 2010–2011 visible in tref.mean('M') plot; real at build time
- GFDL-SPEAR: some recent forecast starts missing; real at build time

**If zeros reappear in a future rebuild**, increment `_TREF_HIND_BUST` / `_TREF_FCST_BUST`
in `code/nmme_models.py` before rebuilding to force IRIDL to bypass the Squid cache.

## zarr v3 Notes

The `pangeo-2025` environment uses zarr 3.x. Key differences from zarr 2:
- Compressor: `zarr.codecs.ZstdCodec(level=5)`, encoding key `compressors=[...]`
- Append: `mode='a'` with `append_dim=` (not `mode='r+'`)
- `_filled` is `uint8`; cast with `.astype(bool)` before logical operations
