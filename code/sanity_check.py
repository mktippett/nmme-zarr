#!/usr/bin/env python
"""
sanity_check.py — Visual sanity check for the NMME SST and tref zarr stores.

Also serves as a worked example of how to open the stores and compute scalar indices.

Figure 1  plots/sanity/store_timeseries.pdf  (skip with --recent-only)
  Top:    Niño 3.4 SST ensemble-mean as (S, L) pcolormesh, one panel per model
  Bottom: Global-mean tref ensemble-mean as (S, L) pcolormesh, one panel per model

Figure 2  plots/sanity/store_recent.pdf
  Top:    Tropical-mean SST as (M, L) pcolormesh, most-recent start per model
  Bottom: Global-mean tref as (M, L) pcolormesh, most-recent start per model

Both figures are single-page PDFs with two sections stacked vertically, one
colorbar per section.  Layout uses fig.subfigures(2, 1) + da.plot(ax=ax) per
panel rather than xarray's col= FacetGrid (which owns its own figure and cannot
be embedded in a combined layout).

Coordinate notes (discrepancies from reference pseudo-code):
  * Store dims are Y / X, not latitude / longitude.
  * 'model' is not a zarr dimension; added via expand_dims + merge.
  * tref is in Kelvin (native IRIDL, ≈ 287 K).
  * tlat was undefined in reference pseudo-code; fixed as parameter (22.5°).

Run:
  mamba run -n pangeo-2025 python code/sanity_check.py --recent-only
  mamba run -n pangeo-2025 python code/sanity_check.py
  mamba run -n pangeo-2025 python code/sanity_check.py --models NASA-GEOSS2S GFDL-SPEAR
"""

import argparse
import logging
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent))

log = logging.getLogger(__name__)

_DEFAULT_STORE_SST  = Path(__file__).parent.parent / "data" / "nmme_sst.zarr"
_DEFAULT_STORE_TREF = Path(__file__).parent.parent / "data" / "nmme_tref.zarr"
_DEFAULT_OUT_DIR    = Path(__file__).parent.parent / "plots" / "sanity"

COL_WRAP = 3
TLAT     = 22.5   # tropical latitude bound for tropical_mean (°)


# ---------------------------------------------------------------------------
# Index helpers — store dim names are Y and X (not latitude/longitude)
# ---------------------------------------------------------------------------

def n34_average(x):
    """Niño 3.4 box mean: X∈[190,240], Y∈[−5,5], cos-lat weighted."""
    x = x.sortby("Y")
    w = np.cos(np.deg2rad(x.Y))
    y = x.sel(X=slice(190, 240)).sel(Y=slice(-5, 5)).weighted(w).mean(["X", "Y"])
    y.attrs = x.attrs.copy()
    return y


def global_mean(x):
    """Cos-lat-weighted global mean over all (Y, X)."""
    x = x.sortby("Y")
    w = np.cos(np.deg2rad(x.Y))
    y = x.weighted(w).mean(["X", "Y"])
    y.attrs = x.attrs.copy()
    return y


def tropical_mean(x, tlat=TLAT):
    """Cos-lat-weighted mean over Y∈[−tlat, tlat]."""
    x = x.sortby("Y")
    w = np.cos(np.deg2rad(x.Y))
    y = x.sel(Y=slice(-tlat, tlat)).weighted(w).mean(["X", "Y"])
    y.attrs = x.attrs.copy()
    return y


# ---------------------------------------------------------------------------
# Build merged datasets with a 'model' dimension
# ---------------------------------------------------------------------------

_DROP = ["X", "Y", "_filled", "T"]


def _build_timeseries(store, var, avg_fn, groups):
    """Apply avg_fn to each model's full (S, M, L) field and merge across models.

    S-axis is the union of all model S-axes; NaN-padding shows gaps in the plot.
    """
    ds_list = []
    for group in groups:
        log.info("[%s] loading %s …", group, var)
        ds = xr.open_zarr(str(store), group=group, decode_times=False)
        da = avg_fn(ds[var]).astype(np.float64).compute()   # (S, M, L)
        ds.close()
        da = da.drop_vars([v for v in _DROP if v in da.coords])
        da = da.assign_coords(model=group).expand_dims("model")
        ds_list.append(da.to_dataset(name=var))
    merged = xr.merge(ds_list)
    merged["S"].attrs.update({"units": "months since 1960-01-01", "calendar": "360"})
    return merged


def _build_recent(store, var, avg_fn, groups):
    """Apply avg_fn to each model's most-recent start (isel(S=-1) before merge).

    Selecting S=-1 per model before merging ensures every panel shows real data.
    Post-merge isel(S=-1) would pick the union's latest S and NaN-fill earlier models.
    Returns (Dataset with dims (model, M, L), dict {group: approx_decimal_year}).
    """
    ds_list = []
    s_years = {}
    for group in groups:
        log.info("[%s] loading recent %s …", group, var)
        ds = xr.open_zarr(str(store), group=group, decode_times=False)
        s_years[group] = 1960.0 + float(ds.S.values[-1]) / 12.0
        da = avg_fn(ds[var].isel(S=-1)).astype(np.float64).compute()  # (M, L)
        ds.close()
        da = da.drop_vars([v for v in _DROP + ["S"] if v in da.coords])
        da = da.assign_coords(model=group).expand_dims("model")
        ds_list.append(da.to_dataset(name=var))
    return xr.merge(ds_list), s_years


# ---------------------------------------------------------------------------
# Core layout: one figure, two sections, ax= per panel
# ---------------------------------------------------------------------------

def _two_section_figure(groups, top_das, bot_das,
                        top_x, top_y, bot_x, bot_y,
                        top_title, bot_title,
                        top_cbar, bot_cbar,
                        suptitle, col_wrap):
    """Build a single figure with top and bottom sections, each a grid of pcolormesh panels.

    Uses subfigures(2,1) so each section gets its own suptitle and colorbar.
    da.plot(ax=ax) (not col=) gives layout control without xarray owning the figure.
    """
    n = len(groups)
    nrows = math.ceil(n / col_wrap)

    fig = plt.figure(figsize=(4.5 * col_wrap + 1, 3.5 * nrows * 2),
                     layout="constrained")
    sf_top, sf_bot = fig.subfigures(2, 1, hspace=0.08)

    axes_t = sf_top.subplots(nrows, col_wrap, squeeze=False)
    axes_b = sf_bot.subplots(nrows, col_wrap, squeeze=False)

    # Shared colorscale per section (2nd–98th percentile across all models)
    def _bounds(das):
        vals = np.concatenate([da.values[np.isfinite(da.values)].ravel() for da in das])
        return float(np.percentile(vals, 2)), float(np.percentile(vals, 98))

    vmin_t, vmax_t = _bounds(top_das)
    vmin_b, vmax_b = _bounds(bot_das)

    last_t = last_b = None
    vis_t, vis_b = [], []

    for i, (group, da_t, da_b) in enumerate(zip(groups, top_das, bot_das)):
        row, col = divmod(i, col_wrap)

        ax = axes_t[row, col]
        im = da_t.plot(x=top_x, y=top_y, ax=ax,
                       vmin=vmin_t, vmax=vmax_t,
                       add_colorbar=False, rasterized=True)
        ax.set_title(f"model = {group}", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.set_xlabel(ax.get_xlabel(), fontsize=7)
        ax.set_ylabel(ax.get_ylabel(), fontsize=7)
        last_t = im
        vis_t.append(ax)

        ax = axes_b[row, col]
        im = da_b.plot(x=bot_x, y=bot_y, ax=ax,
                       vmin=vmin_b, vmax=vmax_b,
                       add_colorbar=False, rasterized=True)
        ax.set_title(f"model = {group}", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.set_xlabel(ax.get_xlabel(), fontsize=7)
        ax.set_ylabel(ax.get_ylabel(), fontsize=7)
        last_b = im
        vis_b.append(ax)

    # Hide unused panels
    for j in range(n, nrows * col_wrap):
        row, col = divmod(j, col_wrap)
        axes_t[row, col].set_visible(False)
        axes_b[row, col].set_visible(False)

    sf_top.colorbar(last_t, ax=vis_t, label=top_cbar, shrink=0.6)
    sf_bot.colorbar(last_b, ax=vis_b, label=bot_cbar, shrink=0.6)
    sf_top.suptitle(top_title, fontsize=11, fontweight="bold")
    sf_bot.suptitle(bot_title, fontsize=11, fontweight="bold")
    fig.suptitle(suptitle, fontsize=12, fontweight="bold")
    fig.set_facecolor("white")
    return fig


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_timeseries(groups, store_sst, store_tref, out_path, col_wrap=COL_WRAP):
    """Figure 1: (S, L) pcolormesh — N34 SST (top) and global tref (bottom)."""
    ds_sst  = _build_timeseries(store_sst,  "sst",  n34_average, groups)
    ds_tref = _build_timeseries(store_tref, "tref", global_mean,  groups)

    top_das = [ds_sst.sst.mean("M").sel(model=g)   for g in groups]
    bot_das = [ds_tref.tref.mean("M").sel(model=g) for g in groups]

    fig = _two_section_figure(
        groups, top_das, bot_das,
        "S", "L", "S", "L",
        "Nino 3.4 totals  (SST °C, ens. mean over M)",
        "Global mean tref totals  (K, ens. mean over M)",
        "SST (°C)", "tref (K)",
        "NMME store sanity check — all starts",
        col_wrap,
    )
    _save(fig, out_path)


def plot_recent(groups, store_sst, store_tref, out_path, col_wrap=COL_WRAP):
    """Figure 2: (M, L) pcolormesh — tropical SST (top) and global tref (bottom)."""
    ds_sst,  yrs_sst  = _build_recent(store_sst,  "sst",  tropical_mean, groups)
    ds_tref, yrs_tref = _build_recent(store_tref, "tref", global_mean,   groups)

    top_das = [ds_sst.sst.sel(model=g)    for g in groups]
    bot_das = [ds_tref.tref.sel(model=g)  for g in groups]

    def _yr_range(d):
        lo, hi = min(d.values()), max(d.values())
        return f"{lo:.2f}" if lo == hi else f"{lo:.2f} – {hi:.2f}"

    fig = _two_section_figure(
        groups, top_das, bot_das,
        "L", "M", "L", "M",
        f"Tropical mean SST — most recent start per model  (≈ {_yr_range(yrs_sst)})",
        f"Global mean tref — most recent start per model  (≈ {_yr_range(yrs_tref)})",
        "SST (°C)", "tref (K)",
        "NMME store sanity check — most recent start",
        col_wrap,
    )
    _save(fig, out_path)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _save(fig, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    log.info("Saved %s", out_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Visual sanity check for the NMME zarr stores")
    parser.add_argument("--recent-only", action="store_true",
                        help="Skip Figure 1 (timeseries); produce only Figure 2 (recent)")
    parser.add_argument("--store-sst",  default=str(_DEFAULT_STORE_SST), metavar="PATH")
    parser.add_argument("--store-tref", default=str(_DEFAULT_STORE_TREF), metavar="PATH")
    parser.add_argument("--out-dir",    type=Path, default=_DEFAULT_OUT_DIR, metavar="DIR")
    parser.add_argument("--models",     nargs="*", default=None, metavar="GROUP")
    parser.add_argument("--col-wrap",   type=int,  default=COL_WRAP, metavar="N")
    parser.add_argument("--log-level",  default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")

    store_sst, store_tref = Path(args.store_sst), Path(args.store_tref)
    for store in (store_sst, store_tref):
        if not store.exists():
            log.error("Store not found: %s", store)
            sys.exit(1)

    import zarr
    all_groups = sorted(zarr.open_group(str(store_sst), mode="r").group_keys())
    if args.models:
        groups = [g for g in all_groups if g in set(args.models)]
        missing = set(args.models) - set(groups)
        if missing:
            log.warning("Groups not found in store: %s", sorted(missing))
    else:
        groups = all_groups

    if not groups:
        log.error("No model groups to process.")
        sys.exit(1)

    log.info("Models: %s", groups)

    if not args.recent_only:
        plot_timeseries(groups, store_sst, store_tref,
                        args.out_dir / "store_timeseries.pdf", args.col_wrap)

    plot_recent(groups, store_sst, store_tref,
                args.out_dir / "store_recent.pdf", args.col_wrap)

    log.info("Done.")


if __name__ == "__main__":
    main()
