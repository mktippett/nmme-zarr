#!/usr/bin/env python
"""
update_archive.py — monthly update of an NMME zarr archive.

For each model:
  1. Query IRIDL DDS to get the current S count for hindcast + forecast.
  2. Determine which S values are missing from the local store.
  3. Append missing starts to the zarr group (append_dim="S").
  4. Re-fetch and overwrite the last 2 starts already in the store
     to absorb late-arriving members or IRIDL corrections.

Usage
-----
  mamba run -n pangeo-local python code/update_archive.py            # both stores, all models
  mamba run -n pangeo-local python code/update_archive.py --var sst  # sst only
  mamba run -n pangeo-local python code/update_archive.py --var sst tref  # explicit both
  mamba run -n pangeo-local python code/update_archive.py --models NCEP-CFSv2
  mamba run -n pangeo-local python code/update_archive.py --recheck-n 5
"""

import argparse
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr
import zarr
import netCDF4

sys.path.insert(0, str(Path(__file__).parent))
from nmme_models import MODELS, MODEL_BY_GROUP, resolve_model, VARIABLES
from iridl_io import fetch_data_block, bust_url

log = logging.getLogger(__name__)

RECHECK_TAIL = 2   # number of latest starts to re-fetch even if present


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_remote_s(url: str) -> np.ndarray:
    """Return remote S coordinate as float32 array."""
    netCDF4.default_timeout = 60
    ds = xr.open_dataset(url, engine="netcdf4", decode_times=False)
    s = ds["S"].values.copy().astype(np.float32)
    ds.close()
    return s


def get_local_s(store_path: Path, group_name: str) -> np.ndarray:
    """Return the S coordinate stored in the local zarr group."""
    root = zarr.open_group(str(store_path), mode="r")
    return root[group_name]["S"][:].astype(np.float32)


def _patch_calendar(store_path: Path, group_name: str) -> None:
    """Fix calendar='360' → '360_day' in S and T metadata (cftime rejects '360')."""
    root = zarr.open_group(str(store_path), mode="r+")
    grp = root[group_name]
    for arr_name in ("S", "T"):
        if arr_name not in grp:
            continue
        attrs = dict(grp[arr_name].attrs)
        if attrs.get("calendar") == "360":
            attrs["calendar"] = "360_day"
            grp[arr_name].attrs.update(attrs)


def append_new_starts(store_path: Path, m: dict,
                      new_s: np.ndarray, M_max: int,
                      L_max: int, s_all_remote: np.ndarray,
                      s_hind: np.ndarray, s_fcst: np.ndarray,
                      Y_vals: np.ndarray, X_vals: np.ndarray,
                      L_vals: np.ndarray, S_attrs: dict) -> None:
    """Fetch new S values and append them to the zarr group."""
    n_new = len(new_s)
    Y, X = len(Y_vals), len(X_vals)
    slab = np.full((n_new, M_max, L_max, Y, X), np.nan, dtype=np.float32)
    T_new = new_s[:, None] + L_vals[None, :]

    for i, sv in enumerate(new_s):
        in_hind = np.isin([sv], s_hind)[0]
        in_fcst = np.isin([sv], s_fcst)[0]

        def fetch_one(url, stream_s, M_stream, target_row):
            pos = int(np.searchsorted(stream_s, sv))
            data, _, _, _, _, _ = fetch_data_block(url, m["var_src"], pos, pos + 1, M_stream)
            slab[target_row, :M_stream, :, :, :] = data[0]

        if m["combined"]:
            pos = int(np.searchsorted(s_all_remote, sv))
            data, _, _, _, _, _ = fetch_data_block(
                m["f_url"], m["var_src"], pos, pos + 1, m["M_fcst"])
            slab[i, :m["M_fcst"], :, :, :] = data[0]
        else:
            if in_hind:
                fetch_one(m["h_url"], s_hind, m["M_hind"], i)
            if in_fcst:
                fetch_one(m["f_url"], s_fcst, m["M_fcst"], i)

    if m["k_offset"]:
        slab -= 273.15

    S_attrs_fixed = dict(S_attrs)
    if S_attrs_fixed.get("calendar") == "360":
        S_attrs_fixed["calendar"] = "360_day"

    ds_new = xr.Dataset(
        {
            m["var_name"]: xr.DataArray(slab, dims=["S", "M", "L", "Y", "X"],
                                        attrs=m["var_attrs"]),
            "_filled": xr.DataArray(np.ones(n_new, dtype=np.uint8), dims=["S"]),
            "T": xr.DataArray(T_new.astype(np.float32), dims=["S", "L"],
                              attrs=dict(units=S_attrs_fixed.get("units",
                                                                  "months since 1960-01-01"),
                                         calendar="360_day")),
        },
        coords={
            "S": xr.DataArray(new_s, dims=["S"], attrs=S_attrs_fixed),
            "L": xr.DataArray(L_vals, dims=["L"]),
            "M": xr.DataArray(np.arange(1, M_max + 1, dtype=np.int32), dims=["M"]),
            "Y": xr.DataArray(Y_vals, dims=["Y"]),
            "X": xr.DataArray(X_vals, dims=["X"]),
        },
    )

    # xarray reads existing S metadata during append validation; cftime rejects
    # calendar='360' (requires '360_day').  Patch the store in-place if needed.
    _patch_calendar(store_path, m["group"])

    log.info("[%s] Appending %d new S values ...", m["group"], n_new)
    # No encoding on append: all variables already exist in the store and
    # xarray raises ValueError if encoding is provided for existing variables.
    ds_new.to_zarr(str(store_path), group=m["group"],
                   mode="a", append_dim="S")


def recheck_tail(store_path: Path, m: dict,
                 local_s: np.ndarray, s_hind: np.ndarray,
                 s_fcst: np.ndarray, s_all_remote: np.ndarray, M_max: int,
                 n: int = RECHECK_TAIL) -> None:
    """Re-fetch and overwrite the last n starts in the local store."""
    tail = local_s[-n:]
    group = zarr.open_group(str(store_path), mode="r+")[m["group"]]

    for sv in tail:
        local_pos = int(np.searchsorted(local_s, sv))
        slab = np.full((1, M_max, m["L_max"], 181, 360), np.nan, dtype=np.float32)
        in_hind = np.isin([sv], s_hind)[0]
        in_fcst = np.isin([sv], s_fcst)[0]

        if m["combined"]:
            pos = int(np.searchsorted(s_all_remote, sv))
            data, _, _, _, _, _ = fetch_data_block(
                m["f_url"], m["var_src"], pos, pos + 1, m["M_fcst"])
            slab[0, :m["M_fcst"], :, :, :] = data[0]
        else:
            if in_hind:
                pos = int(np.searchsorted(s_hind, sv))
                data, _, _, _, _, _ = fetch_data_block(
                    m["h_url"], m["var_src"], pos, pos + 1, m["M_hind"])
                slab[0, :m["M_hind"], :, :, :] = data[0]
            if in_fcst:
                pos = int(np.searchsorted(s_fcst, sv))
                data, _, _, _, _, _ = fetch_data_block(
                    m["f_url"], m["var_src"], pos, pos + 1, m["M_fcst"])
                slab[0, :m["M_fcst"], :, :, :] = data[0]

        if m["k_offset"]:
            slab -= 273.15

        group[m["var_name"]][local_pos:local_pos + 1, :, :, :, :] = slab
        log.info("[%s] Rechecked tail S index %d (value %.1f).",
                 m["group"], local_pos, sv)


# ---------------------------------------------------------------------------
# Per-model update
# ---------------------------------------------------------------------------

def update_model(store_path: Path, m: dict, recheck_n: int = RECHECK_TAIL,
                 bust_cache: bool = False) -> None:
    log.info("=== %s ===", m["group"])

    if bust_cache:
        token = str(int(time.time()))
        m = dict(m)
        for key in ("h_url", "f_url"):
            if m.get(key):
                m[key] = bust_url(m[key], token)
        log.info("[%s] Cache-busting this run with token %s", m["group"], token)

    # --- remote S ---
    if m["combined"]:
        s_all_remote = get_remote_s(m["f_url"])
        s_hind = s_all_remote
        s_fcst = np.array([], dtype=np.float32)
    else:
        s_hind = get_remote_s(m["h_url"])
        s_fcst = get_remote_s(m["f_url"])
        s_all_remote = np.union1d(s_hind, s_fcst)

    # --- local S ---
    local_s = get_local_s(store_path, m["group"])
    new_s = np.setdiff1d(s_all_remote, local_s)

    log.info("[%s] Remote S: %d  Local S: %d  New S: %d",
             m["group"], len(s_all_remote), len(local_s), len(new_s))

    M_max = max(m["M_hind"], m["M_fcst"])

    if len(new_s) > 0:
        ref_url = m["f_url"] if m["combined"] else m["h_url"]
        netCDF4.default_timeout = 60
        ds_ref = xr.open_dataset(ref_url, engine="netcdf4", decode_times=False)
        Y_vals = ds_ref["Y"].values.astype(np.float32)
        X_vals = ds_ref["X"].values.astype(np.float32)
        L_vals = ds_ref["L"].values.astype(np.float32)
        S_attrs = dict(ds_ref["S"].attrs)
        ds_ref.close()
        append_new_starts(store_path, m, new_s, M_max, m["L_max"],
                          s_all_remote, s_hind, s_fcst,
                          Y_vals, X_vals, L_vals, S_attrs)
        # Re-read local_s to include the newly appended
        local_s = get_local_s(store_path, m["group"])

    recheck_tail(store_path, m, local_s, s_hind, s_fcst, s_all_remote, M_max,
                 n=recheck_n)

    # Update timestamp
    root = zarr.open_group(str(store_path), mode="r+")
    root[m["group"]].attrs["last_updated"] = datetime.now(timezone.utc).isoformat()
    log.info("[%s] Update complete.", m["group"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Monthly update of NMME zarr archive")
    parser.add_argument("--var", nargs="+", default=list(VARIABLES),
                        choices=list(VARIABLES),
                        help="Variable(s) to update (default: all)")
    parser.add_argument("--models", nargs="*", default=None)
    parser.add_argument("--store", default=None,
                        help="Path to zarr store; only valid when a single --var is given")
    parser.add_argument("--recheck-n", type=int, default=RECHECK_TAIL,
                        help="Number of tail starts to re-fetch (default: %(default)s)")
    parser.add_argument("--bust-cache", action="store_true",
                        help="Append a fresh per-run cache-busting token to IRIDL URLs, "
                             "bypassing stale Squid cache entries (use if Remote S looks "
                             "stuck below the count shown on the IRIDL data page)")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    if args.store and len(args.var) > 1:
        log.error("--store can only be used when a single --var is specified.")
        sys.exit(1)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    selected = args.models if args.models else [m["group"] for m in MODELS]
    for group_name in selected:
        if group_name not in MODEL_BY_GROUP:
            log.error("Unknown model group: %s", group_name)
            sys.exit(1)

    for var in args.var:
        store_path = Path(args.store) if args.store else VARIABLES[var]["default_store"]
        if not store_path.exists():
            log.error("Store not found: %s  Run build_archive.py first.", store_path)
            sys.exit(1)
        log.info("===== Updating %s store =====", var)
        for group_name in selected:
            update_model(store_path, resolve_model(group_name, var),
                         recheck_n=args.recheck_n, bust_cache=args.bust_cache)


if __name__ == "__main__":
    main()
