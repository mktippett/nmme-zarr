# Spec: extract_seasonal_anomaly.py

## 1. Purpose

Answer one-off data requests from the NMME stores. For one forecast start and
season, write each member's seasonal-mean anomaly for every model to netCDF,
together with each model's climatology and the multi-model mean, and plot
ensemble-mean anomaly maps. The first request was the Oct 2026 start, DJF
tref, 1991–2020 climatology, global and North America.

## 2. Inputs

| Source | Path | Variables | Filters |
|--------|------|-----------|---------|
| tref or sst store | `data/nmme_{var}.zarr/<group>` | `{var}(S, M, L, Y, X)` | S = the forecast start plus the same calendar-month starts in the climatology years; L = the season's leads |

All groups in the store are used. A model that lacks any of the season's leads
is skipped with a message. Opened via `xr.open_datatree(..., decode_times=False)`.

## 3. Outputs

| File | Contents | Format |
|------|----------|--------|
| `data/requests/nmme_{var}_anom_{season}_{start}_start_clim{clim}{suffix}.nc` | `{var}_anom(model, member, lat, lon)`, `{var}_clim(model, lat, lon)`, `{var}_anom_mme(lat, lon)`, `n_members(model)`, `n_clim_samples(model)`; global attrs record start, season, leads, clim period, and region bounds for a regional file | netCDF4, float32, zlib level 4 + shuffle |
| `plots/requests/<same stem>.png` | 4×2 maps: ensemble-mean anomaly per model (title shows member count) plus the multi-model mean | PNG, 150 dpi |

One netCDF and one PNG are written per entry in `--regions`; `suffix` is empty
for `global` and `_namerica` for `namerica`. `data/` is gitignored, so the
outputs are not under version control.

## 4. Algorithm

1. Season leads: for each season month m and start month s0,
   `L = ((m − s0) mod 12) + 0.5`. DJF from an Oct start gives L = 2.5, 3.5, 4.5.
   A season that wraps past the start month is rejected.
2. Start values: `S = (year − 1960)·12 + (month − 1)`. Every climatology start
   and the forecast start must exist in the store, or the script raises an error.
3. For each model, load `{var}.sel(S=[clim starts + forecast start], L=leads)`.
   Seasonal mean = mean over L with `skipna=False`, so a member missing any
   lead gets NaN.
4. Climatology `C(Y, X)` = mean of the seasonal means over all climatology
   starts and all members (NaN members skipped). The climatology depends on
   model, start month and lead.
5. Anomaly = forecast-start seasonal mean of each member − C.
6. Stack models along `model` with an outer join on M, so members a model lacks
   are NaN (M renamed `member`, Y→`lat`, X→`lon`). `n_members` = members with
   any finite value; `n_clim_samples` = member-starts with any finite value in
   the climatology.
7. Multi-model mean = mean over `member` (skipna), then an equal-weight mean
   over `model`.
8. For each region: subset lat/lon (none for `global`), write the netCDF, and plot.
9. Plot: contourf of the ensemble mean on Robinson (global, centred at 180°) or
   Lambert conformal (central lon −100°, central lat 45°, with country and state
   borders). Regional maps are drawn from the **global** field and clipped with
   `set_extent`, so the lat/lon box, which is curved on the projection, fills
   the frame; the corners therefore show data just outside the netCDF box.

## 5. Constants & Scientific Rationale

| Name | Value | Why |
|------|-------|-----|
| `--clim` default | 1991–2020 | Requested by the user (current WMO normal period) |
| `REGIONS["namerica"]` | lat 10–75°N, lon 190–310°E (170°W–50°W) | Chosen by Claude to cover Alaska to Newfoundland and Mexico/Central America; user did not specify bounds |
| Plot levels | −4 to 4 K, step 0.5, extend both | Chosen by Claude for the Oct 2026 request (Niño-3.4 ensemble means reach ~4.5 K); not tuned to other requests |
| netCDF compression | zlib level 4, shuffle | Chosen by Claude; reduces the 57 MB global file to ~18 MB |

## 6. Edge Cases & Error Handling

- Member counts differ between hindcast and forecast (GFDL-SPEAR: 15 hindcast,
  30 forecast; NASA-GEOSS2S: 4–10 hindcast members by year). The climatology
  averages all available hindcast members, so it is unbiased but noisier where
  members are few.
- NCEP-CFSv2 has 24 members except for November starts (28), so members 25–28
  are NaN by design (see CLAUDE.md Data Quality).
- Forecast members delivered late by IRIDL are absent until
  `update_archive.py --recheck-n` re-fetches the start; rerun this script afterward.
- `--var sst`: land points are NaN in the store and remain NaN.
- GFDL-SPEAR tref is NaN on the ±90° rows for every start in the store
  (±89° are populated), so its climatology, anomalies and the multi-model mean
  at the poles are NaN or come from the other 6 models. This is passed through
  as stored, not masked (see CLAUDE.md Data Quality).

**Open item (2026-10-07):** NASA-GEOSS2S's Oct 2026 DJF tref anomalies are much
larger and noisier than the other models' (widespread >4 K over the Southern
Ocean and Antarctica, sharp cold patches near Greenland, the Weddell Sea and the
US East Coast). Not investigated. Decision needed: compare NASA's recent
Oct-start DJF anomalies with the other models to test for a hindcast/forecast
system mismatch, and if one exists, decide whether to note it in deliverables
or use a shorter, recent climatology for that model.

## 7. Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-10-07 | Initial implementation: per-member seasonal anomalies, model dimension, climatology, multi-model mean, netCDF with compression | Yes |
| 2026-10-07 | Added the 4×2 ensemble-mean map figure | Yes |
| 2026-10-07 | Replaced `--out` with `--regions` (`global`, `namerica` presets); regional maps drawn from the global field and clipped | Yes |

## Verification Snippet

```python
import numpy as np, xarray as xr
ds = xr.open_dataset("data/requests/nmme_tref_anom_DJF_2026-10_start_clim1991-2020.nc")
assert ds.sizes["model"] == 7
# finite-member count matches n_members for each model
nfin = ds.tref_anom.notnull().any(["lat", "lon"]).sum("member")
assert (nfin == ds.n_members).all()
# MME equals the member-then-model mean
mme = ds.tref_anom.mean("member").mean("model")
assert np.allclose(mme, ds.tref_anom_mme, equal_nan=True, atol=1e-5)
# climatology complete everywhere except GFDL-SPEAR's NaN pole rows
interior = ds.tref_clim.sel(lat=slice(-89, 89))
assert interior.notnull().all()
print("OK", ds.n_members.to_series().to_dict())
```
