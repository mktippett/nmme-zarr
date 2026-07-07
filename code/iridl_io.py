"""
iridl_io.py — low-level helpers for reading IRIDL OPeNDAP endpoints.

Public API:
  dds_dims(url)             → dict of dim_name → size
  open_dods(url, ...)       → xr.Dataset (decode_times=False), with retries
  fetch_data_block(url, var_src, s_indices, m_max) → np.ndarray (S, M, L, Y, X)
"""

import re
import time
import logging
import urllib.request

import netCDF4
import numpy as np
import xarray as xr

log = logging.getLogger(__name__)

# Default network settings
DEFAULT_TIMEOUT = 600   # seconds per request
DEFAULT_RETRIES = 5
DEFAULT_BACKOFF = 2.0   # seconds; doubles each retry

# IRIDL seems to choke on requests larger than ~500 MB;
# requesting 12 starts × 28 members × 10 leads × 181 × 360 × 4 B ≈ 52 MB is safe.
# build_archive drives block size; this module doesn't enforce it.


def bust_url(url: str, token: str) -> str:
    """Insert a fresh /<token>/pop/ segment before the trailing /dods.

    IRIDL's Squid cache keys on the full URL text. Pushing an arbitrary
    number onto the ingrid stack and popping it is a numeric no-op but
    changes the cache key, forcing a fresh fetch past a stale cache entry.
    Safe to stack on top of any /pop/ segments already baked into the URL
    (e.g. the static _TREF_HIND_BUST-style constants in nmme_models.py).
    """
    if url.endswith("/dods"):
        return f"{url[:-len('/dods')]}/{token}/pop/dods"
    return f"{url}/{token}/pop"


def dds_dims(url: str) -> dict[str, int]:
    """Return {dim_name: size} by parsing the OPeNDAP DDS text endpoint.

    Works without authentication (public endpoint).
    URL can be the base dods URL (with or without /dods suffix).
    """
    dds_url = url.rstrip("/")
    if not dds_url.endswith(".dds"):
        dds_url = dds_url.replace("/dods", "") + "/dods.dds"
        if "/dods.dds" not in dds_url:
            dds_url = dds_url + "/dods.dds"
    # Simpler: just append .dds to whatever we were given
    base = url.rstrip("/")
    if base.endswith("/dods"):
        base = base[:-5]
    dds_url = base + "/dods.dds"

    log.debug("DDS: %s", dds_url)
    with urllib.request.urlopen(dds_url, timeout=30) as resp:
        text = resp.read().decode("utf-8")

    dims = {}
    # Match lines like: Float32 S[S = 414];
    for m in re.finditer(r"Float32 (\w+)\[(\w+) = (\d+)\]", text):
        var, dim, size = m.group(1), m.group(2), int(m.group(3))
        if var == dim:   # coordinate declarations look like this
            dims[dim] = size
    return dims


def open_dods(url: str, retries: int = DEFAULT_RETRIES,
              backoff: float = DEFAULT_BACKOFF,
              timeout: int = DEFAULT_TIMEOUT) -> xr.Dataset:
    """Open an IRIDL OPeNDAP URL as an xr.Dataset (decode_times=False).

    Retries on any exception with exponential backoff.
    """
    netCDF4.default_timeout = timeout
    wait = backoff
    last_exc = None
    for attempt in range(retries):
        try:
            ds = xr.open_dataset(url, engine="netcdf4", decode_times=False)
            return ds
        except Exception as exc:
            last_exc = exc
            log.warning("Attempt %d/%d failed for %s: %s", attempt + 1, retries, url, exc)
            if attempt < retries - 1:
                time.sleep(wait)
                wait *= 2
    raise RuntimeError(f"All {retries} attempts failed for {url}") from last_exc


def _build_subset_url(base_url: str, s_indices: list[int],
                      m_max: int, L: int, Y: int, X: int) -> str:
    """Build an OPeNDAP constraint expression URL for a list of S indices.

    IRIDL OPeNDAP allows [start:stride:stop] slice notation per dim.
    For a non-contiguous S list we must issue one request per contiguous run
    (handled by the caller).  Here we expect s_indices to be a contiguous range.
    """
    s0, s1 = s_indices[0], s_indices[-1]
    # Constraint: [S0:S1][0:M-1][0:L-1][0:Y-1][0:X-1]
    # Dim order in IRIDL NMME datasets is always [S][M or L][...] but varies.
    # We request the full variable without a CE and use xarray subsetting,
    # because the dim ordering varies by model. Slicing S via the URL
    # does require knowing the dim position.  Use the "S/(range)" IRIDL
    # filtering instead — but that's only available via the DL ingrid language,
    # not vanilla OPeNDAP CE. So we just open the full URL and isel in Python.
    # See fetch_data_block below.
    raise NotImplementedError("use fetch_data_block directly")


def fetch_data_block(
    url: str,
    var_src: str,
    s_start: int,
    s_stop: int,   # exclusive, python-slice convention
    M_max: int,
    retries: int = DEFAULT_RETRIES,
    backoff: float = DEFAULT_BACKOFF,
    timeout: int = DEFAULT_TIMEOUT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Download a block of start times from an IRIDL OPeNDAP endpoint.

    Parameters
    ----------
    url      : full /dods URL for the model stream (hindcast or forecast)
    var_src  : IRIDL variable name (e.g. 'sst' or 'sst_regridded')
    s_start  : first S index (0-based integer index into the S coordinate)
    s_stop   : one-past-last S index
    M_max    : total members after padding (to check shape)

    Returns
    -------
    data : float32 ndarray shaped (S, M, L, Y, X) — M axis may be < M_max
           if this stream has fewer members; padding is the caller's job.
    S_vals, M_vals, L_vals, Y_vals, X_vals : 1-D coordinate arrays (raw,
           not decoded) for the returned block.
    """
    netCDF4.default_timeout = timeout
    wait = backoff
    last_exc = None

    for attempt in range(retries):
        try:
            ds = xr.open_dataset(url, engine="netcdf4", decode_times=False)
            blk = ds[var_src].isel(S=slice(s_start, s_stop)).load()
            # Drop degenerate extra dims (e.g. Z=2m level in tref datasets)
            extra = [d for d in blk.dims if d not in ("S", "M", "L", "Y", "X")]
            if extra:
                blk = blk.squeeze(extra)
            # Normalize dim order to (S, M, L, Y, X)
            blk = blk.transpose("S", "M", "L", "Y", "X")
            data = blk.values.astype(np.float32)
            S_vals = ds["S"].values[s_start:s_stop]
            M_vals = ds["M"].values
            L_vals = ds["L"].values
            Y_vals = ds["Y"].values
            X_vals = ds["X"].values
            ds.close()
            return data, S_vals, M_vals, L_vals, Y_vals, X_vals
        except Exception as exc:
            last_exc = exc
            log.warning(
                "Attempt %d/%d failed for %s [S=%d:%d]: %s",
                attempt + 1, retries, url, s_start, s_stop, exc,
            )
            if attempt < retries - 1:
                time.sleep(wait)
                wait *= 2

    raise RuntimeError(
        f"All {retries} attempts failed for {url} S={s_start}:{s_stop}"
    ) from last_exc
