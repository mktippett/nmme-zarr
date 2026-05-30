#!/usr/bin/env python
"""
build_archive.py — create (or resume) an NMME zarr archive.

Usage
-----
  mamba run -n pangeo-2025 python code/build_archive.py
  mamba run -n pangeo-2025 python code/build_archive.py --var tref
  mamba run -n pangeo-2025 python code/build_archive.py --models NASA-GEOSS2S
  mamba run -n pangeo-2025 python code/build_archive.py --models NASA-GEOSS2S --block-size 6

  --var        variable to archive: sst or tref (default: sst)
  --models     subset of group names from nmme_models.py (default: all 7)
  --store      path to output zarr store (default: data/nmme_<var>.zarr)
  --block-size number of S indices per OPeNDAP request (default: 12)

The store is created once with the full shape inferred from the IRIDL DDS
endpoints, then filled block-by-block.  Interrupted runs resume automatically
by checking the _filled sentinel.
"""

import argparse
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr
import zarr

# Allow running from either code/ or the project root
sys.path.insert(0, str(Path(__file__).parent))
from nmme_models import MODELS, MODEL_BY_GROUP, resolve_model, VARIABLES
from iridl_io import dds_dims, fetch_data_block

log = logging.getLogger(__name__)

DEFAULT_BLOCK = 12   # S indices per request (12 = 1 year, keeps each request ≤ ~78 MB)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_s_values(url: str) -> np.ndarray:
    """Return the S coordinate values from the IRIDL DDS (raw float32)."""
    import urllib.request, re
    base = url.rstrip("/")
    if base.endswith("/dods"):
        base = base[:-5]
    dds_url = base + "/dods.dds"
    das_url = base + "/dods.das"

    # DDS gives us the count; DAS gives us the units.
    # But we need the actual values, so we open a tiny slice.
    import netCDF4
    netCDF4.default_timeout = 60
    ds = xr.open_dataset(url, engine="netcdf4", decode_times=False)
    s = ds["S"].values.copy()
    ds.close()
    return s.astype(np.float32)


def build_s_union(m: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Query IRIDL to get the combined S axis for a model.

    Returns
    -------
    s_all : sorted union of S values across hindcast + forecast
    s_hind: S values from hindcast stream (or s_all if combined)
    s_fcst: S values from forecast stream (or empty if combined)
    M_max : max ensemble members across both streams
    L_max : lead count from hindcast/forecast (take from m dict)
    """
    log.info("[%s] Querying S coordinates ...", m["group"])
    if m["combined"]:
        s_all = get_s_values(m["f_url"])
        s_hind = s_all
        s_fcst = np.array([], dtype=np.float32)
    else:
        s_hind = get_s_values(m["h_url"])
        s_fcst = get_s_values(m["f_url"])
        s_all = np.union1d(s_hind, s_fcst)
    M_max = max(m["M_hind"], m["M_fcst"])
    return s_all, s_hind, s_fcst, M_max, m["L_max"]


def init_group(store_path: Path, m: dict, s_all: np.ndarray,
               M_max: int, L_max: int,
               Y_vals: np.ndarray, X_vals: np.ndarray,
               L_vals: np.ndarray, S_attrs: dict) -> None:
    """Create an empty (NaN-filled) zarr group for one model."""
    import zarr.codecs
    compressor = zarr.codecs.ZstdCodec(level=5)

    n_S = len(s_all)
    chunks = m["chunks"]
    chunk_tuple = (
        chunks["S"], chunks["M"], chunks["L"], chunks["Y"], chunks["X"]
    )

    # Patch calendar: IRIDL uses '360' but xarray needs '360_day'
    S_attrs_fixed = dict(S_attrs)
    if S_attrs_fixed.get("calendar") == "360":
        S_attrs_fixed["calendar"] = "360_day"

    # Compute T = S + L as a 2-D array
    T_vals = s_all[:, None] + L_vals[None, :]   # (S, L)
    T_attrs = dict(
        units=S_attrs_fixed.get("units", "months since 1960-01-01"),
        calendar="360_day",
        long_name="Target (valid) time",
        standard_name="forecast_reference_time",
    )

    ds_empty = xr.Dataset(
        {
            m["var_name"]: xr.DataArray(
                data=np.full((n_S, M_max, L_max, len(Y_vals), len(X_vals)),
                             np.nan, dtype=np.float32),
                dims=["S", "M", "L", "Y", "X"],
                attrs=m["var_attrs"],
            ),
            "_filled": xr.DataArray(
                data=np.zeros(n_S, dtype=bool),
                dims=["S"],
                attrs=dict(long_name="Block has been written"),
            ),
            "T": xr.DataArray(
                data=T_vals.astype(np.float32),
                dims=["S", "L"],
                attrs=T_attrs,
            ),
        },
        coords={
            "S": xr.DataArray(s_all, dims=["S"], attrs=S_attrs_fixed),
            "L": xr.DataArray(L_vals, dims=["L"],
                              attrs=dict(units="months", long_name="Lead")),
            "M": xr.DataArray(np.arange(1, M_max + 1, dtype=np.int32),
                              dims=["M"],
                              attrs=dict(long_name="Ensemble Member")),
            "Y": xr.DataArray(Y_vals, dims=["Y"],
                              attrs=dict(units="degree_north",
                                         standard_name="latitude")),
            "X": xr.DataArray(X_vals, dims=["X"],
                              attrs=dict(units="degree_east",
                                         standard_name="longitude")),
        },
        attrs=dict(
            model=m["label"],
            nmme_hindcast_url=m.get("h_url", "combined"),
            nmme_forecast_url=m["f_url"],
            created=datetime.now(timezone.utc).isoformat(),
            last_updated="",
        ),
    )

    log.info("[%s] Initialising zarr group (shape %s) ...", m["group"],
             ds_empty[m["var_name"]].shape)
    encoding = {
        m["var_name"]: dict(compressors=[compressor], chunks=chunk_tuple, dtype="float32",
                            fill_value=np.float32("nan")),
        "_filled": dict(dtype="uint8"),
        "T": dict(compressors=[compressor], dtype="float32"),
    }
    ds_empty.to_zarr(
        str(store_path), group=m["group"], mode="w",
        encoding=encoding,
    )


def fill_group(store_path: Path, m: dict, s_all: np.ndarray,
               s_hind: np.ndarray, s_fcst: np.ndarray,
               M_max: int, block_size: int) -> None:
    """Fill a zarr group block-by-block, skipping already-filled blocks."""
    group = zarr.open_group(str(store_path), mode="r+")[m["group"]]
    filled_flags = group["_filled"][:]

    n_S = len(s_all)
    for blk_start in range(0, n_S, block_size):
        blk_end = min(blk_start + block_size, n_S)

        # Skip if entire block is already filled
        if np.all(filled_flags[blk_start:blk_end]):
            log.info("[%s] S[%d:%d] already filled, skipping.",
                     m["group"], blk_start, blk_end)
            continue

        s_slice_vals = s_all[blk_start:blk_end]

        # Determine source URL and indices for each S in this block
        _fill_block(store_path, m, s_all, s_hind, s_fcst,
                    M_max, blk_start, blk_end, s_slice_vals, filled_flags)

        # Re-read flags (fill_block writes them)
        filled_flags = zarr.open_group(str(store_path), mode="r")[m["group"]]["_filled"][:]

    # Update last_updated attribute
    root = zarr.open_group(str(store_path), mode="r+")
    root[m["group"]].attrs["last_updated"] = datetime.now(timezone.utc).isoformat()


def _qa_slab(slab, s_slice_vals, group_name):
    """Warn if any start in the slab has all-zero or spatially constant data.

    Checks only finite (non-NaN) values; NaN-only starts are skipped because
    they represent S positions outside both hindcast and forecast streams.
    """
    for i, s_val in enumerate(s_slice_vals):
        finite = slab[i][np.isfinite(slab[i])]
        if finite.size == 0:
            continue
        year = int(s_val) // 12 + 1960
        if np.all(finite == 0):
            log.warning("[%s] S=%.0f (year %d): all-zero — possible missing/fill data from IRIDL",
                        group_name, float(s_val), year)
        elif finite.std() == 0:
            log.warning("[%s] S=%.0f (year %d): constant=%.4g — possible missing/fill data from IRIDL",
                        group_name, float(s_val), year, float(finite[0]))


def _fill_block(store_path, m, s_all, s_hind, s_fcst,
                M_max, blk_start, blk_end, s_slice_vals, filled_flags):
    """Fetch and write one S block, handling H/F boundary and member padding."""
    # Partition the S values in this block into hindcast and forecast
    in_hind = np.isin(s_slice_vals, s_hind)
    in_fcst = np.isin(s_slice_vals, s_fcst)

    # Allocate output slab for this block
    n_s = blk_end - blk_start
    L_max = m["L_max"]
    Y = 181
    X = 360
    slab = np.full((n_s, M_max, L_max, Y, X), np.nan, dtype=np.float32)

    def write_stream(url, mask, M_stream):
        """Fetch a subset of S from one stream (H or F) and place into slab."""
        if not np.any(mask):
            return
        # global indices in s_all for this stream partition
        global_idx = np.where(mask)[0]   # relative to block
        # We need the index of each S value within the stream's own S array
        stream_s = s_hind if url == m.get("h_url") else s_fcst
        if m["combined"]:
            stream_s = s_all
        pos_in_stream = np.searchsorted(stream_s, s_slice_vals[mask])
        # Fetch in one call (indices are a contiguous sub-range if in order)
        si0 = int(pos_in_stream[0])
        si1 = int(pos_in_stream[-1]) + 1
        log.info("[%s] Fetching %s S[%d:%d] ...", m["group"],
                 "hindcast" if url == m.get("h_url") else "forecast", si0, si1)
        data, _, _, _, _, _ = fetch_data_block(url, m["var_src"], si0, si1,
                                              M_stream)
        # data shape: (si1-si0, M_stream, L_max, Y, X)
        # Place into slab at the correct M and S positions
        for k, rel_idx in enumerate(global_idx):
            src_s = k if (si1 - si0 == len(global_idx)) else np.searchsorted(
                np.arange(si0, si1), pos_in_stream[k])
            slab[rel_idx, :M_stream, :, :, :] = data[src_s]

    if m["combined"]:
        write_stream(m["f_url"], np.ones(n_s, dtype=bool), m["M_fcst"])
    else:
        write_stream(m["h_url"], in_hind, m["M_hind"])
        write_stream(m["f_url"], in_fcst, m["M_fcst"])

    # Apply K → °C conversion
    if m["k_offset"]:
        slab -= 273.15

    _qa_slab(slab, s_slice_vals, m["group"])

    # Write to zarr
    group = zarr.open_group(str(store_path), mode="r+")[m["group"]]
    group[m["var_name"]][blk_start:blk_end, :, :, :, :] = slab
    group["_filled"][blk_start:blk_end] = True
    log.info("[%s] S[%d:%d] written.", m["group"], blk_start, blk_end)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_model(store_path: Path, m: dict, block_size: int) -> None:
    log.info("=== %s ===", m["group"])
    s_all, s_hind, s_fcst, M_max, L_max = build_s_union(m)
    log.info("[%s] S count: %d  M_max: %d  L: %d",
             m["group"], len(s_all), M_max, L_max)

    # Get Y, X, L coordinate values from hindcast (or forecast if combined)
    ref_url = m["f_url"] if m["combined"] else m["h_url"]
    log.info("[%s] Fetching dim metadata from %s ...", m["group"], ref_url)
    import netCDF4
    netCDF4.default_timeout = 60
    ds_ref = xr.open_dataset(ref_url, engine="netcdf4", decode_times=False)
    Y_vals = ds_ref["Y"].values.astype(np.float32)
    X_vals = ds_ref["X"].values.astype(np.float32)
    L_vals = ds_ref["L"].values.astype(np.float32)
    S_attrs = dict(ds_ref["S"].attrs)
    ds_ref.close()

    # Check if group already initialised
    root = zarr.open_group(str(store_path), mode="a")
    if m["group"] not in root:
        init_group(store_path, m, s_all, M_max, L_max, Y_vals, X_vals, L_vals, S_attrs)
    else:
        log.info("[%s] Group already exists, resuming fill.", m["group"])

    fill_group(store_path, m, s_all, s_hind, s_fcst, M_max, block_size)
    log.info("[%s] Done.", m["group"])


def main():
    parser = argparse.ArgumentParser(description="Build NMME zarr archive")
    parser.add_argument("--var", default="sst", choices=list(VARIABLES),
                        help="Variable to archive (default: sst)")
    parser.add_argument("--models", nargs="*", default=None,
                        help="Group names to process (default: all)")
    parser.add_argument("--store", default=None,
                        help="Path to output zarr store (default: data/nmme_<var>.zarr)")
    parser.add_argument("--block-size", type=int, default=DEFAULT_BLOCK,
                        help="Number of S indices per OPeNDAP request")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    store_path = Path(args.store) if args.store else VARIABLES[args.var]["default_store"]
    store_path.parent.mkdir(parents=True, exist_ok=True)

    selected = args.models if args.models else [m["group"] for m in MODELS]
    for group_name in selected:
        if group_name not in MODEL_BY_GROUP:
            log.error("Unknown model group: %s", group_name)
            sys.exit(1)
        build_model(store_path, resolve_model(group_name, args.var), args.block_size)


if __name__ == "__main__":
    main()
