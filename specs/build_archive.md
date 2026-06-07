# Spec: build_archive.py

## 1. Purpose

Create a local Zarr archive of monthly forecasts from seven active NMME
models for a chosen variable (`sst` or `tref`), enabling efficient offline
analysis without repeated IRIDL OPeNDAP requests. Supports resume after
interruption.

## 2. Inputs

| Source | Details |
|--------|---------|
| IRIDL OPeNDAP DDS endpoints | Public (no auth). Used to determine S counts and dim values before downloading data. Accessed via `iridl_io.dds_dims` and `xr.open_dataset`. |
| IRIDL OPeNDAP dods endpoints | One per model stream (hindcast / forecast or combined). Accessed via `iridl_io.fetch_data_block` in blocks of `--block-size` S indices (default 12). |
| `nmme_models.py` | Group name, URLs, `var_src`, `k_offset`, member counts, chunk shapes for each of the 7 models. `VARIABLES` maps each variable name to its store path and zarr attrs. `resolve_model(group, var)` returns the flat dict consumed by the script. |

All grids: 360×181 (1°), S in "months since 1960-01-01" (360-day calendar),
L in months, M integer member index.

## 3. Outputs

| File/Object | Contents | Format |
|-------------|---------|--------|
| `data/nmme_<var>.zarr/<group>/<var>` | Float32 variable data, dims (S, M, L, Y, X) | Zarr array, Zstd level 5 |
| `data/nmme_<var>.zarr/<group>/T` | Float32 target time = S + L, dims (S, L) | Zarr array |
| `data/nmme_<var>.zarr/<group>/_filled` | uint8 sentinel, 1 where block written, dim (S,) | Zarr array |
| Coordinate arrays S, M, L, Y, X | As raw IRIDL values | Zarr arrays |
| Group `.zattrs` | model, URLs, created, last_updated | JSON |

Seven groups total (see `nmme_models.py` for group names).
Actual on-disk size after compression (2026-06): ~76 GB (`nmme_sst.zarr`)
+ ~100 GB (`nmme_tref.zarr`) ≈ 176 GB total.

## 4. Algorithm

1. **Parse `--var`** (default `sst`): look up `VARIABLES[var]` for default store
   path and zarr attrs. Call `resolve_model(group, var)` for each model to get
   a flat dict containing `h_url`, `f_url`, `var_src`, `k_offset`, `var_name`,
   `var_attrs`.
2. **For each model** (selected via `--models` or all 7):
   a. Query the S coordinate from each stream (H and F) via a small
      `xr.open_dataset` call (metadata only).
   b. Take the sorted union `s_all = union(s_hind, s_fcst)`.  
      `M_max = max(M_hind, M_fcst)`.
   c. **Initialize group** (if not present) by writing an empty Dataset
      (all NaN, false `_filled`) with the correct shape and chunk encoding
      via `to_zarr(mode="w")`.  Variable name and attrs come from `m["var_name"]`
      and `m["var_attrs"]`.  T = S + L is computed and stored at this point.
   d. **Fill blocks**: iterate `blk_start` in `range(0, len(s_all), block_size)`.
      - Skip if `all(_filled[blk_start:blk_end])`.
      - Partition S values in the block into hindcast vs. forecast using
        `np.isin`.
      - Call `fetch_data_block(url, var_src, si0, si1, M_stream)` for each
        non-empty partition.  `si0/si1` are the indices within the stream's
        own S axis, found via `np.searchsorted`.
      - Assemble a slab of shape `(n_s, M_max, L, Y, X)`, placing hindcast
        data into `slab[:, :M_hind, ...]` and forecast data into
        `slab[:, :M_fcst, ...]`.  The remaining M slots stay NaN.
      - Subtract 273.15 if `k_offset=True`.
      - Run `_qa_slab(slab, s_slice_vals, group_name)`: emit a WARNING for
        any start whose finite values are all-zero or spatially constant
        (fill-value signatures from IRIDL Squid cache misses).
      - Write slab and set `_filled[blk_start:blk_end] = True` via direct
        zarr array indexing.
3. Update the group attribute `last_updated` timestamp.

`fetch_data_block` (in `iridl_io.py`) opens the URL, does `.isel(S=slice)`,
transposes to `(S, M, L, Y, X)`, and loads. It retries up to 5 times with
exponential backoff starting at 2 s.

## 5. Constants & Scientific Rationale

| Name | Value | Why |
|------|-------|-----|
| `DEFAULT_BLOCK` | 12 | One year of starts per request. Worst-case raw size ≈ 78 MB (CanSIPS, M=20), safely under any IRIDL per-request cap. Aligns with seasonal groupbys. |
| Zstd level | 5 | Good compression of smooth SST fields with modest CPU overhead. zarr v3: `zarr.codecs.ZstdCodec(level=5)`, encoding key `compressors=[...]`. |
| `fill_value` | `np.float32("nan")` | Unwritten chunks must be NaN, not 0. zarr's default float32 fill_value is 0.0, which is physically plausible for anomaly fields and silently masquerades as data. Explicitly setting NaN means any un-fetched S position is correctly missing. |
| `DEFAULT_TIMEOUT` | 600 s | OPeNDAP transfers of ~78 MB on academic networks can be slow. |
| `DEFAULT_RETRIES` | 5 | IRIDL occasionally drops connections mid-transfer. |
| `DEFAULT_BACKOFF` | 2.0 s (doubles) | Gives server time to recover without hammering. |
| `_TREF_HIND_BUST` / `_TREF_FCST_BUST` | `"1"` / `"5"` (in `nmme_models.py`) | IRIDL uses a Squid caching proxy. Adding `/N/pop/` before `/dods` changes the cache key without changing the data, bypassing stale cache entries that return zeros. Increment both before any tref rebuild if cache corruption is suspected; use distinct values for hind and fcst so both keys are unique. |
| Chunk S | 12 | Aligns with calendar year; monthly appends fill predictably. |
| Chunk M | 4–5 | Keeps per-chunk size 52–78 MB (float32). |
| Chunk L, Y, X | full | One full map per (S, M) tile; enables efficient map reads. |

## 6. Edge Cases & Error Handling

- **Resume**: `_filled` sentinel is checked before each block. A run
  interrupted mid-block leaves that block's data in the zarr but `_filled`
  still False, so the block is re-fetched cleanly on resume (safe because
  we overwrite the entire block). `_filled` is stored as `uint8`; when
  reading, cast with `.astype(bool)` — zarr v3 returns integers and `~int(1) = -2`
  is truthy, which breaks `np.where(~filled)` without the cast.
- **Member count mismatch** (NASA H=4, F=10; SPEAR H=15, F=30): the slab is
  allocated to `M_max` and only the relevant sub-range is populated.
  NaN padding makes the per-S member availability transparent.
- **Lead count mismatch**: taken from `m["L_max"]` in `nmme_models.py`;
  NASA has 9, others 10 or 12.  Verify with DDS before assuming.
- **COLA models combined**: `h_url=None`, `combined=True`; no hindcast/
  forecast partitioning needed.
- **GFDL-SPEAR variable rename**: `var_src="sst_regridded"` is fetched and
  stored under `var_name="sst"` in the zarr. For tref, GFDL uses plain `.tref`
  (confirmed; `sst_regridded` is the exception, not the pattern).
- **K → °C**: CanSIPS-IC4 SST is delivered in Kelvin (`k_offset=True`).
  For tref, all models use `k_offset=False` and data is stored as native K.
- **Non-contiguous S in a block**: possible if IRIDL has gaps. `np.searchsorted`
  finds the position correctly; `fetch_data_block` fetches the contiguous
  range `[si0, si1)` which may include gaps. Gaps remain NaN in the slab.
- **IRIDL Squid cache returning zeros**: on repeated requests IRIDL may serve a
  stale cached response that contains all-zero data (verified for tref in 2026).
  `_qa_slab` detects this at write time and logs a WARNING. The fix is to
  increment `_TREF_HIND_BUST` / `_TREF_FCST_BUST` in `nmme_models.py` and
  rebuild. Do **not** add downstream masks; fix the source.
- **fill_value=0 in pre-2026-05 stores**: stores built before the `fill_value=nan`
  fix have zarr `fill_value=0.0` in their encoding metadata. Unwritten chunks in
  those stores silently return 0 instead of NaN. The only safe remedy is a full
  rebuild; there is no in-place patch for zarr fill_value.

## 7. Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|-------------|
| 2026-04-06 | Initial implementation | 2026-04-06 |
| 2026-04-06 | zarr v3 fixes: `ZstdCodec`, `compressors=[...]`, `dtype='uint8'` for `_filled`, calendar patched to `'360_day'` at write time | 2026-04-06 |
| 2026-04-22 | Parameterised for multiple variables: `--var {sst,tref}` CLI flag; `VARIABLES` table and `resolve_model` in `nmme_models.py`; `var_name`/`var_attrs` keys in model dict replace hardcoded `"sst"` literals; `fetch_sst_block` renamed `fetch_data_block` | 2026-04-22 |
| 2026-05-27 | Added `fill_value=np.float32("nan")` to zarr encoding in `init_group` — prevents unwritten chunks from silently returning 0.0 instead of NaN. Added `_qa_slab()`: logs WARNING for any start with all-zero or spatially constant finite values (IRIDL Squid cache-miss signature). Added `_TREF_HIND_BUST`/`_TREF_FCST_BUST` cache-busting URL constants to `nmme_models.py`; all 7 tref URLs updated to include `/N/pop/dods` pattern. tref store rebuilt from scratch with these fixes. | 2026-05-27 |

---

## Verification snippet

```python
# Run after build for one model (e.g. NASA-GEOSS2S)
import xarray as xr, zarr, numpy as np

store = "data/nmme_sst.zarr"
ds = xr.open_zarr(store, group="NASA-GEOSS2S")

# Shape checks
assert ds.sst.dims == ("S", "M", "L", "Y", "X")
assert ds.sst.shape[1] == 10   # M_max
assert ds.sst.shape[2] == 9    # L
assert ds.sst.shape[3] == 181  # Y
assert ds.sst.shape[4] == 360  # X

# _filled: all True
root = zarr.open_group(store, mode="r")
assert root["NASA-GEOSS2S"]["_filled"][:].all(), "Some blocks not filled"

# T = S + L (sample)
T = ds["T"].values
S = ds["S"].values
L = ds["L"].values
np.testing.assert_allclose(T[0, :], S[0] + L, rtol=1e-5)
np.testing.assert_allclose(T[-1, :], S[-1] + L, rtol=1e-5)

# SST range sanity (-2 to 35 °C)
assert float(ds.sst.min()) > -5, "SST too cold"
assert float(ds.sst.max()) < 40, "SST too warm"
print("All checks passed.")
```
