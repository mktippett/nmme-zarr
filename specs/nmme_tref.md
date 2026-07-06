# Spec: nmme_tref zarr archive

## 1. Purpose

Parallel archive to `data/nmme_sst.zarr` for 2-m reference air temperature
(`tref`), using the same seven NMME models and identical build/update
machinery. Enables offline analysis of surface temperature forecast skill
without repeated IRIDL requests.

## 2. Store layout

**Path:** `data/nmme_tref.zarr`  
**Groups:** same seven model groups as the SST store (see `nmme_models.py`).

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `tref` | (S, M, L, Y, X) | K | 2-m reference air temperature (native Kelvin) |
| `T` | (S, L) | months since 1960-01-01 | Valid time = S + L |
| `_filled` | (S,) | uint8 | 1 = block written, 0 = not yet fetched |

Coordinates S, M, L, Y, X: same conventions as SST store (360-day calendar,
1° grid, 181 × 360).

## 3. IRIDL URL pattern

Base URL: `{NMME_BASE}/.<ModelPath>/.{HINDCAST|FORECAST}/.MONTHLY/.tref/N/pop/dods`

The `/N/pop/` segment (where N is an integer) changes the Squid cache key
without changing the returned data, bypassing stale cached responses. Constants
`_TREF_HIND_BUST` and `_TREF_FCST_BUST` in `nmme_models.py` control these
values (currently `"1"` and `"5"` respectively; increment both before any
rebuild if IRIDL cache corruption is suspected).

NCEP-CFSv2 has an extra path segment:
`.NCEP-CFSv2/.HINDCAST/.PENTAD_SAMPLES/.MONTHLY/.tref/N/pop/dods`

COLA models (combined=True) use the forecast URL for both streams, with no
`.HINDCAST` segment.

**GFDL-SPEAR note:** SST uses `.sst_regridded`; tref uses plain `.tref`
(confirmed; `sst_regridded` is the exception, not the pattern).

## 4. Build / update

```bash
# Preflight: verify one model, one block
mamba run -n pangeo-local python code/build_archive.py \
    --var tref --models NASA-GEOSS2S --block-size 1

# Full build (all 7 models)
mamba run -n pangeo-local python code/build_archive.py --var tref

# Monthly update
mamba run -n pangeo-local python code/update_archive.py --var tref
```

## 5. Known data-quality issues and rebuild history

### IRIDL Squid cache zeros (original build, pre-2026-05-27)
The original tref store was built without cache-busting URLs and with zarr
`fill_value=0.0` (the default). IRIDL's Squid proxy returned all-zero
responses for some model/year combinations; these were stored as 0.0 K rather
than NaN.

| Model | Affected starts | Pattern |
|-------|-----------------|---------|
| CanSIPS-IC4-CanESM5 | 48 starts (Jan 2007–Dec 2010) | all-zero |
| CanSIPS-IC4-GEM52NEMO | 60 starts (scattered 2007–2025) | all-zero |
| GFDL-SPEAR | 84 starts (1999–2025 partial) | all-zero |
| COLA×2, NASA, NCEP | 0 | clean |

SST store was unaffected by *this* zero-fill incident (different URL paths;
different cache state at the time). A separate Squid staleness symptom was
later found on NCEP-CFSv2's sst forecast URL (stale/short `S` axis hiding the
newest start, not zero-fill) — see `specs/build_archive.md` §5/§7, 2026-07-06.

**Downstream symptom:** zeros in the tref climatology (clim ~260 K instead of
~299 K), producing ~40 K anomaly offsets and flat lag correlations for the
affected models.

**Fix applied 2026-05-27:**
1. Added `fill_value=np.float32("nan")` to zarr encoding in `build_archive.py`.
2. Added `_qa_slab()` QA in `build_archive.py` (WARNING on all-zero or constant starts).
3. Updated all tref URLs in `nmme_models.py` to use `/N/pop/dods` cache-busting.
4. Deleted and rebuilt `data/nmme_tref.zarr` from scratch.

### Real data gaps (not defects)
- **NCEP-CFSv2**: NaN gap around 2010–2011 visible as a white stripe in the
  tref mean plot. Real at build time; NaN is the correct representation.
- **GFDL-SPEAR**: some recent forecast starts (2024–2025) missing; real at
  build time.

## 6. Differences from SST archive

| Property | SST | tref |
|----------|-----|------|
| Store path | `data/nmme_sst.zarr` | `data/nmme_tref.zarr` |
| Units | Celsius_scale | K |
| `k_offset` | True (CanSIPS only) | False (all models) |
| GFDL var_src | `sst_regridded` | `tref` |
| standard_name | `sea_surface_temperature` | `air_temperature` |
| URL cache-bust | NCEP-CFSv2 only, since 2026-07-06 (`/N/pop/dods`) | all 7 models (see Section 3) |

## 7. Synchronization Log

| Date | Change | Spec updated |
|------|--------|--------------|
| 2026-04-22 | Initial tref spec created alongside `--var tref` support in build/update scripts | 2026-04-22 |
| 2026-05-27 | Documented Squid cache-zero incident; recorded affected models; added Sections 3 (cache-busting URL), 5 (data-quality issues); updated GFDL var_src note; removed "(verify)" qualifier | 2026-05-27 |
| 2026-07-06 | Corrected §5 "SST unaffected" note and §6 cache-bust row: NCEP-CFSv2's sst URLs now use `/N/pop/` too (separate stale-S-axis incident; full details in `specs/build_archive.md`) | 2026-07-06 |

## 8. Verification snippet

Run after a full build to confirm no zero-fill and sensible Kelvin values:

```python
import xarray as xr, zarr, numpy as np, sys
sys.path.insert(0, "code")
from nmme_models import MODELS

store = "data/nmme_tref.zarr"
root = zarr.open_group(store, mode="r")

for m in MODELS:
    group = m["group"]
    if group not in root:
        print(f"MISSING: {group}"); continue
    ds = xr.open_zarr(store, group=group, decode_times=False)
    assert ds.tref.dims == ("S", "M", "L", "Y", "X"), f"{group}: wrong dims"
    assert ds.tref.attrs["units"] == "K", f"{group}: units not K"

    # _filled: all written
    assert root[group]["_filled"][:].astype(bool).all(), f"{group}: blocks missing"

    # Kelvin range: a hindcast start, tropical band, first lead
    sf = ds.S.values.astype(float)
    s0 = int(np.where(sf >= (1991-1960)*12)[0][0])
    t = ds.tref.isel(S=s0, L=0, M=0).sel(Y=slice(-30, 30)).compute()
    assert float(t.min()) > 240, f"{group}: tref too cold ({float(t.min()):.1f} K)"
    assert float(t.max()) < 340, f"{group}: tref too warm ({float(t.max()):.1f} K)"

    # No all-zero starts: check tropical mean over all starts
    trop = ds.tref.isel(M=0, L=0).sel(Y=slice(-22.5, 22.5)).mean(["Y","X"]).compute()
    n_zero = int((trop < 1.0).sum())   # < 1 K catches genuine zeros
    if n_zero > 0:
        print(f"WARN {group}: {n_zero} starts with trop-mean < 1 K (possible zeros)")
    else:
        print(f"OK   {group}: clean  shape={ds.tref.shape}")
    ds.close()

print("tref verification done.")
```
