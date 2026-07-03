# Spec: update_archive.py

## 1. Purpose

Incrementally update an existing NMME zarr archive (built by
`build_archive.py`) to include new monthly forecasts issued since the last
update, and to refresh the most recent two starts in case IRIDL has revised
them (late-arriving members, corrections). By default updates **both** stores
(`sst` and `tref`); pass `--var sst` or `--var tref` to restrict to one.
Each variable has its own independent store.

## 2. Inputs

| Source | Details |
|--------|---------|
| IRIDL OPeNDAP dods endpoints | Same URLs as `build_archive.py`. Used to retrieve current S coordinates and new data. |
| `data/nmme_sst.zarr`, `data/nmme_tref.zarr` | Existing zarr stores (both updated by default). Must already have all model groups. |
| `nmme_models.py` | Model config via `resolve_model(group, var)` (same as build). |

## 3. Outputs

| Change | Description |
|--------|-------------|
| Appended S slices | New start times appended to each model group via `append_dim="S"`. |
| Overwritten tail | Last 2 S indices overwritten in-place (zarr region write). |
| `last_updated` attr | Updated per-group. |

## 4. Algorithm

For each model:

1. **Retrieve remote S** from IRIDL (same `get_remote_s` calls as build).
2. **Retrieve local S** from `zarr.open_group(...)[group]["S"][:]`.
3. **New starts** = `setdiff1d(remote_S, local_S)`.
4. If new starts exist:
   - For each new S value, determine which stream (H or F) it belongs to
     via `np.isin`.
   - Fetch one start at a time (small request, simpler logic than block mode).
   - Assemble NaN-padded slab `(1, M_max, L, 181, 360)`.
   - Subtract 273.15 if `k_offset`.
   - Collect all new starts into a Dataset (variable named `m["var_name"]`,
     attrs from `m["var_attrs"]`) and call
     `ds_new.to_zarr(store, group=group, mode="a", append_dim="S")`.
5. **Recheck tail**: for the last `RECHECK_TAIL=2` starts in the local store
   (after any append), re-fetch and overwrite in-place via direct zarr array
   indexing `group[m["var_name"]][local_pos:local_pos+1] = slab`.
6. Update `last_updated` attribute.

## 5. Constants & Scientific Rationale

| Name | Value | Why |
|------|-------|-----|
| `RECHECK_TAIL` | 2 (default) | IRIDL sometimes updates the most recent 1–2 months as new members arrive. Override with `--recheck-n N` when more starts had partial members at build time. |
| Timeout, retries, backoff | Same as build | Same rationale. |

## 6. Edge Cases & Error Handling

- **No new starts**: only the tail recheck runs. The script exits normally.
- **New model added to IRIDL** (model retired/added): if a group is missing
  from the local store, `update_archive.py` will fail cleanly with an error.
  Run `build_archive.py --models <new_group>` first.
- **Remote S shorter than local S**: rare (IRIDL retraction). The script
  only appends; it does not delete local starts. Log a warning.
- **Append creates partial chunk**: zarr handles trailing partial chunks
  correctly; no special treatment needed.
- **T coordinate on append**: the appended Dataset includes `T[S, L] = S + L`
  so the T array stays consistent without rewriting old data.
- **xarray/zarr version skew**: `ds_new.to_zarr(..., mode="a", append_dim="S")`
  opens existing store variables (`open_store_variable`) as part of the
  append, which decodes `_FillValue` against the array's zarr dtype metadata.
  Under `zarr>=3.2`, that metadata is a `ZDType` object (e.g. `Float32(...)`)
  with no `.value` attribute; xarray releases older than the zarr 3.2 dtype
  refactor raise `AttributeError: 'Float32' object has no attribute 'value'`
  on every read/write of the store, not just appends. Requires a matched
  xarray/zarr pair (verified: `xarray 2026.4.0` + `zarr 3.2.1` in
  `pangeo-local`). Not an IRIDL issue.

## 7. Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|-------------|
| 2026-04-06 | Initial implementation | 2026-04-06 |
| 2026-04-06 | zarr v3 fixes: `ZstdCodec`, `compressors=[...]`, `mode='a'` for `append_dim`; `--recheck-n N` CLI argument added | 2026-04-06 |
| 2026-04-22 | Parameterised for multiple variables: `--var {sst,tref}` CLI flag; `resolve_model` replaces direct `MODEL_BY_GROUP` lookup; `m["var_name"]`/`m["var_attrs"]` replace hardcoded `"sst"` literals; `fetch_sst_block` renamed `fetch_data_block` | 2026-04-22 |
| 2026-05-05 | Removed `encoding` argument entirely from `append_new_starts`; xarray raises `ValueError` if encoding is supplied for any variable that already exists in the store during append | 2026-05-05 |
| 2026-05-05 | Added `_patch_calendar()` helper; xarray also reads existing S/T metadata during append validation and cftime rejects `calendar='360'` (requires `'360_day'`); patch fixes the store in-place before each append | 2026-05-05 |
| 2026-05-30 | `--var` changed to `nargs="+"` defaulting to all variables; bare invocation now updates both `sst` and `tref` stores; `--store` now errors if multiple vars are selected | 2026-05-30 |
| 2026-07-03 | Switched documented environment from `pangeo-2025` to `pangeo-local` (usage docstring, README) after `pangeo-2025`'s unpinned `zarr>=3` drifted to 3.2.1 while its `xarray` (2025.3.1) did not, breaking all reads/writes to zarr-format-3 stores; documented as an edge case above | 2026-07-03 |

---

## Verification snippet

```python
# Run immediately after update_archive.py on an already-complete store.
# Should report zero new starts and identical tail after recheck.
import zarr, numpy as np

store = "data/nmme_sst.zarr"
root = zarr.open_group(store, mode="r")

for group_name in root.group_keys():
    g = root[group_name]
    local_s = g["S"][:]
    print(f"{group_name}: {len(local_s)} starts, "
          f"last S = {local_s[-1]:.1f}, "
          f"last_updated = {g.attrs.get('last_updated', 'n/a')}")

# Manual check: open one model and confirm SST is not all-NaN for recent starts
import xarray as xr
ds = xr.open_zarr(store, group="NCEP-CFSv2")
recent = ds.sst.isel(S=slice(-3, None), M=0, L=0)
assert not recent.isnull().all(), "Recent SST is all NaN — update may have failed"
print("Tail check passed.")
```
