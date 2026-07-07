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
| QA warnings | `qa_member_consistency` WARNINGs for members inconsistent with their ensemble on the touched starts (log only — data never modified). |
| `last_updated` attr | Updated per-group. |
| Consolidated metadata | Root `zarr.json` re-consolidated after each model's update so consolidated-metadata readers see the run's writes. |

## 4. Algorithm

For each model:

0. **Optional cache-bust** (`--bust-cache`): if set, generate one token per
   model per run (`str(int(time.time()))`) and rewrite `h_url`/`f_url` via
   `bust_url()` (in `iridl_io.py`), which inserts a fresh `/<token>/pop/`
   segment before the trailing `/dods`. This token is stacked on top of any
   static bust segment already baked into the URL by `nmme_models.py`
   (`_TREF_HIND_BUST` etc.) and is used consistently for the remote-S probe
   and all subsequent data fetches in that model's run.
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
6. **Member-consistency QA** on every start touched this run (new appends ∪
   rechecked tail): compute the cos-lat-weighted Niño-3.4 box average
   (X∈[190,240], Y∈[−5,5], same box as `sanity_check.n34_average`) per
   (member, lead) from the *local store after writing*, and flag any member
   whose |deviation from the ensemble median| exceeds
   `QA_Z_THRESH * max(1.4826 * MAD, QA_SCALE_FLOOR)` at any lead. WARNING
   log only — flagged data is never masked or modified (fix-at-source
   discipline: verify against IRIDL, report upstream, repair later with
   `--recheck-n`).
7. Update `last_updated` attribute and re-consolidate store metadata
   (`zarr.consolidate_metadata`), so consolidated-metadata readers — in
   particular fingerprint-based downstream caches such as nmme_enso's
   `load_nino34_ssta` — see this run's writes.

## 5. Constants & Scientific Rationale

| Name | Value | Why |
|------|-------|-----|
| `RECHECK_TAIL` | 2 (default) | IRIDL sometimes updates the most recent 1–2 months as new members arrive. Override with `--recheck-n N` when more starts had partial members at build time. |
| Timeout, retries, backoff | Same as build | Same rationale. |
| Cache-bust token | `str(int(time.time()))` per model per run, only when `--bust-cache` set | Needs to be unique per run (unlike the static `_TREF_*_BUST`/`_SST_*_BUST` constants in `nmme_models.py`, which bust the cache once and then become stale cache keys themselves once Squid has seen that exact URL). Wall-clock time is a convenient, always-fresh, human-debuggable token; the numeric value has no semantic meaning to IRIDL. |
| `QA_Z_THRESH` | 6.0 | Robust-z threshold for the member-consistency QA. Calibrated against the two known corrupt CESM1 members (z ≈ 13–56, flagged) vs. the largest legitimate deviations found in a full-store lead-0 scan — lagged-ensemble init-date spread (CFSv2/GEOSS2S, up to ~1.4 °C at lead 0 with correspondingly larger MAD) and long-lead ensemble divergence (both pass at z < 6). |
| `QA_SCALE_FLOOR` | 0.15 °C | Floor on the MAD-based spread scale. At lead 0 non-lagged ensembles are near-identical (median lead-0 member range 0.12–0.38 °C across models), so an unfloored MAD would flag numerical noise; the floor sets the minimum flaggable deviation to `6 × 0.15 = 0.9 °C`, comfortably above real lead-0 agreement and below the ≥ 2.8 °C corrupt cases. |
| QA region | Niño-3.4 box (X∈[190,240], Y∈[−5,5]) | Cheap scalar per (member, lead), sensitive to the tropical-Pacific field this archive exists to serve; same box as `sanity_check.n34_average`. Applied to whichever variable is being updated (the check compares differences, so the °C/K offset distinction doesn't matter). |

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
- **Stale Squid cache (`Remote S` stuck below IRIDL's page count)**: IRIDL's
  Squid proxy caches on exact URL text, including any static `/pop/`
  bust segment baked into `nmme_models.py`. Once that exact busted URL has
  been requested, it too becomes a stale cache key — the symptom recurs on
  any model/variable, not just the ones previously documented in the
  project CLAUDE.md (NCEP-CFSv2 sst, tref zero-fills). Confirmed on
  2026-07-07 for NASA-GEOSS2S tref: the plain run reported `Remote S: 545`
  (matching local, so "nothing new") while IRIDL's own data-selection page
  showed a forecast start through 1 Jul 2026. Re-running with `--bust-cache`
  immediately picked up the missing start (`Remote S: 546`). Diagnostic: compare
  the script's `Remote S` count against the IRIDL variable's data-selection
  HTML page (not just `dds_dims()`, which is served through the same cache).
  Fix: rerun the affected model(s) with `--bust-cache`; no code constant
  needs to be hand-incremented.
- **Why the QA is a robust z at every lead, not a fixed lead-0 threshold**:
  the two motivating corrupt members (see Synchronization Log 2026-07-07)
  had opposite signatures — one was +2.8 °C wrong *at* lead 0, the other was
  correct at lead 0 and ~3 °C wrong only at leads 2–8. A fixed lead-0
  threshold misses the second entirely, and any fixed all-lead threshold
  either false-alarms on legitimate long-lead ensemble divergence (±1.5–2.5 °C
  from the median is normal at lead ~10 during ENSO events) or on
  lagged-ensemble models (NCEP-CFSv2 is 4 members × 6 init dates; NASA-GEOSS2S
  similar — their members legitimately disagree by >1 °C at nominal lead 0
  during rapid SST transitions, in coherent blocks by init date). The
  MAD-based scale absorbs both: wherever legitimate spread is large, the
  tolerance is large.
- **QA flags are advisory, never corrective**: per the project's
  fix-at-source discipline, a flag means "verify this member against IRIDL
  and report upstream if confirmed" — the store keeps exactly what IRIDL
  delivered so the defect stays visible and reportable.
- **Consolidated-metadata staleness (fixed 2026-07-07)**: xarray's `to_zarr`
  re-consolidates root metadata on append, but a recheck-only run (raw zarr
  region writes + the attrs stamp) previously did not. Readers that open the
  store through consolidated metadata — including nmme_enso's
  `_store_fingerprint`, which keys its expensive-load cache on `last_updated`
  — saw a stale stamp, so a recheck that repaired data *without* appending
  new starts would not invalidate downstream caches (silent stale-data
  hazard). `update_model` now ends with `zarr.consolidate_metadata()`
  unconditionally. (zarr emits a `ZarrUserWarning` that consolidated
  metadata is not part of the v3 spec — expected, harmless, and pre-existing:
  the stores have carried consolidated metadata since they were built with
  xarray.)
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
| 2026-07-07 | Added `--bust-cache` CLI flag and `bust_url()` helper (`iridl_io.py`); generates a fresh timestamp token per model per run to bypass stale Squid cache entries, generalizing the previously model-specific static-constant workaround | 2026-07-07 |
| 2026-07-07 | Added `qa_member_consistency()` (robust-z member-consistency QA on the Niño-3.4 box, run on every start touched per model; `QA_Z_THRESH=6.0`, `QA_SCALE_FLOOR=0.15`; WARNING-only). Motivated by two corrupt COLA-RSMAS-CESM1 members delivered by IRIDL and confirmed unchanged under a cache-busted re-fetch (2026-03 M=5: +2.79 °C at leads 0.5–3.5; 2026-07 M=1: −3.04 °C at leads 1.5–8.5); both flagged by the new QA, reported upstream by MKT. | 2026-07-07 |
| 2026-07-07 | `update_model` now ends with `zarr.consolidate_metadata()`; previously a recheck-only run never refreshed root consolidated metadata, so fingerprint-based downstream caches (nmme_enso `load_nino34_ssta`) saw a stale `last_updated` and would not invalidate after an in-place data repair. Verified end-to-end: nmme_enso's `_store_fingerprint` sees the new stamp immediately after a recheck-only run. | 2026-07-07 |

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
