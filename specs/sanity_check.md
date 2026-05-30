# Spec: sanity_check.py

## 1. Purpose

Provide a fast visual sanity check that the NMME SST and tref zarr stores are
correctly populated and physically plausible.  Also serves as a worked example
of how to open the stores and compute scalar indices.

## 2. Inputs

| Source | Path | Variables | Notes |
|--------|------|-----------|-------|
| SST store | `data/nmme_sst.zarr/<group>` | `sst(S, M, L, Y, X)` float32 °C | 7 model groups |
| tref store | `data/nmme_tref.zarr/<group>` | `tref(S, M, L, Y, X)` float32 K | Same 7 groups |

Opened with `decode_times=False` (raw float S axis, no calendar decoding needed).

## 3. Outputs

| File | Contents |
|------|---------|
| `plots/sanity/store_timeseries.pdf` | Figure 1: N34 SST and global-mean tref time series |
| `plots/sanity/store_recent.pdf` | Figure 2: tropical-mean SST and global-mean tref vs lead for most-recent start |

Default output directory is `plots/sanity/`; override with `--out-dir`.

## 4. Algorithm

### Build pattern — `_build_timeseries` and `_build_recent`

Both figures use the same **merge pattern**:

```python
for group in groups:
    ds = xr.open_zarr(store, group=group, decode_times=False)
    da = avg_fn(ds[var]).astype(np.float64).compute()   # spatial reduction in loop
    da = da.drop_vars([...]).assign_coords(model=group).expand_dims('model')
    ds_list.append(da.to_dataset(name=var))
merged = xr.merge(ds_list)   # union S-axis; NaN where model has no data → shows gaps
```

Spatial reduction inside the loop keeps peak memory to one model at a time.
`xr.merge` (not `xr.concat`) is correct because models have different S ranges;
the union join NaN-pads missing starts, which produces the visible white gaps in
Figure 1.  The same NaN structure across M reveals per-model member counts.

### Figure layout — `_two_section_figure`

Both figures use a single-page PDF with two sections stacked vertically via
`fig.subfigures(2, 1)`.  Each section is a `nrows × col_wrap` grid of axes filled
with `da.plot(x=..., y=..., ax=ax, add_colorbar=False)`.  One shared colorbar
(2nd–98th percentile across all models) is added per section via `sf.colorbar()`.

`da.plot(ax=ax)` (not `col='model'`) is used because xarray's `col=` FacetGrid
creates and owns its own figure and cannot be embedded in a combined layout.

### Figure 1 — (S, L) pcolormesh (ensemble mean over M)

1. `_build_timeseries(store_sst, 'sst', n34_average, groups)` → `(model, S, M, L)`.
2. `_build_timeseries(store_tref, 'tref', global_mean, groups)` → `(model, S, M, L)`.
3. Top section: `ds_sst.sst.mean('M').sel(model=g).plot(x='S', y='L', ax=ax)` per group.
4. Bottom section: same for `ds_tref.tref.mean('M')`.

### Figure 2 — (M, L) pcolormesh (most-recent start, per model)

**Improvement over post-merge `isel(S=-1)`:** `isel(S=-1)` is applied
*per model before the merge*, so each panel shows that model's own most-recent
start.  Post-merge `isel(S=-1)` would pick the union's latest S and NaN-fill
models whose last start is earlier — hiding real data.

1. `_build_recent(store_sst, 'sst', tropical_mean, groups)` → `(model, M, L)` + per-model years.
2. `_build_recent(store_tref, 'tref', global_mean, groups)` → same.
3. Top section: `ds_sst.sst.sel(model=g).plot(x='L', y='M', ax=ax)` per group.
4. Bottom section: same for `ds_tref.tref`.
5. Suptitle shows the range of per-model most-recent start years.

### Index helper functions

All helpers accept `(..., Y, X)` DataArrays; return `(...)` with Y and X reduced.
`sortby("Y")` is called before any `sel(Y=slice(...))` for coord-order safety.

| Function | Region / weighting |
|----------|--------------------|
| `n34_average(x)` | Y∈[−5, 5], X∈[190, 240]; cos-lat weighted |
| `global_mean(x)` | All Y and X; cos-lat weighted |
| `tropical_mean(x, tlat=22.5)` | Y∈[−tlat, tlat]; cos-lat weighted; tlat is a parameter |

## 5. Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| `TLAT` | 22.5° | Tropical band for tropical_mean; consistent with enso-t2m analysis scripts |
| `COL_WRAP` | 3 | 3-column FacetGrid; override with `--col-wrap` |
| `isel(S=-1)` | positional last start | Per-model; safe on raw-float S axis |

### Coordinate name notes (discrepancies from reference pseudo-code)

The reference used `latitude`/`longitude` dim names and `col='model'` (not a store
dimension).  The store uses `Y` (latitude) and `X` (longitude).  `col='model'` is
made available by the merge pattern: `assign_coords(model=group).expand_dims('model')`
+ `xr.merge` before plotting.

The reference `ds.sst.mean('M').plot(x='S')` would fail because after spatial
reduction the array is `(S, M, L)` — still 3D.  The pattern reduces spatially
inside the loop, then `mean('M')` at plot time gives `(model, S, L)`, which
xarray renders as a pcolormesh with `x='S', y='L'` per panel.  Adding `y=` was
the only change needed from the reference code.

## 6. Edge Cases & Error Handling

- **Store not found**: `sys.exit(1)` with an error message before any I/O.
- **`--models` group not in store**: logged as WARNING; remaining valid groups are processed.
- **No valid groups**: `sys.exit(1)`.
- **Unused panels** (`n` not a multiple of `col_wrap`): hidden with `ax.set_visible(False)`.
- **`--recent-only`**: skips Figure 1; only Figure 2 is produced.

## 7. Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-05-30 | Initial implementation | Yes |
| 2026-05-30 | Rewrite: switch from manual line/pcolormesh loops to merge pattern + xarray FacetGrid. `isel(S=-1)` moved before merge in `_build_recent` to show each model's own most-recent start. `--col-wrap` CLI flag added. | Yes |
| 2026-05-30 | Second rewrite: replaced PdfPages two-page approach and FacetGrid `col=` with `fig.subfigures(2,1)` + `da.plot(ax=ax)` per panel. Both figures are now single-page with SST top / tref bottom sections. One colorbar per section via `sf.colorbar()`. `_two_section_figure` helper centralises the layout. | Yes |
