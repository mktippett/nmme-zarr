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
```

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
| `code/update_archive.py` | Monthly incremental update |
| `code/sanity_check.py` | Visual sanity check and worked example of reading the stores |

Behavioral specs for each script live in `specs/`.

## Data Quality Notes

- `build_archive.py` logs a WARNING for any forecast start with all-zero or constant data. Check warnings before running analysis.
- If zeros appear in a tref rebuild, increment `_TREF_HIND_BUST` / `_TREF_FCST_BUST` in `code/nmme_models.py` to bypass the IRIDL Squid cache.
- NCEP-CFSv2: a real data gap around 2010–2011 is visible in tref.
- GFDL-SPEAR: some recent forecast starts may be missing at build time.
