"""
nmme_models.py — single source of truth for the 7 NMME model entries.

Each entry describes:
  group     : zarr group name (safe ASCII, no dots)
  label     : human-readable model name
  combined  : True when hindcast and forecast share a single URL
  M_hind    : ensemble member count in hindcast stream
  M_fcst    : ensemble member count in forecast stream (M used = max of both)
  L_max     : number of lead months
  chunks    : zarr chunk shape dict {S, M, L, Y, X}
  vars      : per-variable sub-dict keyed by variable name ("sst", "tref", ...)
              Each sub-dict has: h_url, f_url, var_src, k_offset

VARIABLES maps each supported variable to its default store path and zarr attrs.

Use resolve_model(group, var) to get a flat dict ready for build_archive / update_archive.
"""

from pathlib import Path

NMME_BASE = "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME"

# IRIDL uses a Squid caching proxy. Appending an integer + /pop/ to a URL
# changes the cache key without changing the returned data, forcing a fresh
# fetch past stale cache entries. Increment these constants before any rebuild
# if you suspect the previous build cached bad data.
_TREF_HIND_BUST = "1"   # appended as /<n>/pop/ before /dods in tref hindcast URLs
_TREF_FCST_BUST = "5"   # different value for forecast so both keys are unique
DATA_DIR = Path(__file__).parent.parent / "data"

VARIABLES = {
    "sst": dict(
        default_store=DATA_DIR / "nmme_sst.zarr",
        attrs=dict(
            units="Celsius_scale",
            long_name="Sea Surface Temperature",
            standard_name="sea_surface_temperature",
        ),
    ),
    "tref": dict(
        default_store=DATA_DIR / "nmme_tref.zarr",
        attrs=dict(
            units="K",
            long_name="Reference (2-m) Air Temperature",
            standard_name="air_temperature",
        ),
    ),
}

MODELS = [
    dict(
        group="CanSIPS-IC4-CanESM5",
        label="CanSIPS-IC4 / CanESM5",
        combined=False,
        M_hind=20, M_fcst=20, L_max=12,
        chunks=dict(S=12, M=4, L=12, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=f"{NMME_BASE}/.CanSIPS-IC4/.CanESM5/.HINDCAST/.MONTHLY/.sst/dods",
                f_url=f"{NMME_BASE}/.CanSIPS-IC4/.CanESM5/.FORECAST/.MONTHLY/.sst/dods",
                var_src="sst", k_offset=True,
            ),
            tref=dict(
                h_url=f"{NMME_BASE}/.CanSIPS-IC4/.CanESM5/.HINDCAST/.MONTHLY/.tref/{_TREF_HIND_BUST}/pop/dods",
                f_url=f"{NMME_BASE}/.CanSIPS-IC4/.CanESM5/.FORECAST/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
    dict(
        group="CanSIPS-IC4-GEM52NEMO",
        label="CanSIPS-IC4 / GEM5.2-NEMO",
        combined=False,
        M_hind=20, M_fcst=20, L_max=12,
        chunks=dict(S=12, M=4, L=12, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=f"{NMME_BASE}/.CanSIPS-IC4/.GEM5.2-NEMO/.HINDCAST/.MONTHLY/.sst/dods",
                f_url=f"{NMME_BASE}/.CanSIPS-IC4/.GEM5.2-NEMO/.FORECAST/.MONTHLY/.sst/dods",
                var_src="sst", k_offset=True,
            ),
            tref=dict(
                h_url=f"{NMME_BASE}/.CanSIPS-IC4/.GEM5.2-NEMO/.HINDCAST/.MONTHLY/.tref/{_TREF_HIND_BUST}/pop/dods",
                f_url=f"{NMME_BASE}/.CanSIPS-IC4/.GEM5.2-NEMO/.FORECAST/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
    dict(
        group="GFDL-SPEAR",
        label="GFDL-SPEAR",
        combined=False,
        M_hind=15, M_fcst=30, L_max=12,
        chunks=dict(S=12, M=5, L=12, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=f"{NMME_BASE}/.GFDL-SPEAR/.HINDCAST/.MONTHLY/.sst_regridded/dods",
                f_url=f"{NMME_BASE}/.GFDL-SPEAR/.FORECAST/.MONTHLY/.sst_regridded/dods",
                var_src="sst_regridded", k_offset=False,
            ),
            tref=dict(
                h_url=f"{NMME_BASE}/.GFDL-SPEAR/.HINDCAST/.MONTHLY/.tref/{_TREF_HIND_BUST}/pop/dods",
                f_url=f"{NMME_BASE}/.GFDL-SPEAR/.FORECAST/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
    dict(
        group="COLA-RSMAS-CCSM4",
        label="COLA-RSMAS-CCSM4",
        combined=True,
        M_hind=10, M_fcst=10, L_max=12,
        chunks=dict(S=12, M=5, L=12, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=None,
                f_url=f"{NMME_BASE}/.COLA-RSMAS-CCSM4/.MONTHLY/.sst/dods",
                var_src="sst", k_offset=False,
            ),
            tref=dict(
                h_url=None,
                f_url=f"{NMME_BASE}/.COLA-RSMAS-CCSM4/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
    dict(
        group="COLA-RSMAS-CESM1",
        label="COLA-RSMAS-CESM1",
        combined=True,
        M_hind=10, M_fcst=10, L_max=12,
        chunks=dict(S=12, M=5, L=12, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=None,
                f_url=f"{NMME_BASE}/.COLA-RSMAS-CESM1/.MONTHLY/.sst/dods",
                var_src="sst", k_offset=False,
            ),
            tref=dict(
                h_url=None,
                f_url=f"{NMME_BASE}/.COLA-RSMAS-CESM1/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
    dict(
        group="NASA-GEOSS2S",
        label="NASA-GEOSS2S",
        combined=False,
        M_hind=4, M_fcst=10, L_max=9,
        chunks=dict(S=12, M=5, L=9, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=f"{NMME_BASE}/.NASA-GEOSS2S/.HINDCAST/.MONTHLY/.sst/dods",
                f_url=f"{NMME_BASE}/.NASA-GEOSS2S/.FORECAST/.MONTHLY/.sst/dods",
                var_src="sst", k_offset=False,
            ),
            tref=dict(
                h_url=f"{NMME_BASE}/.NASA-GEOSS2S/.HINDCAST/.MONTHLY/.tref/{_TREF_HIND_BUST}/pop/dods",
                f_url=f"{NMME_BASE}/.NASA-GEOSS2S/.FORECAST/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
    dict(
        group="NCEP-CFSv2",
        label="NCEP-CFSv2",
        combined=False,
        M_hind=28, M_fcst=28, L_max=10,
        chunks=dict(S=12, M=4, L=10, Y=181, X=360),
        vars=dict(
            sst=dict(
                h_url=f"{NMME_BASE}/.NCEP-CFSv2/.HINDCAST/.PENTAD_SAMPLES/.MONTHLY/.sst/dods",
                f_url=f"{NMME_BASE}/.NCEP-CFSv2/.FORECAST/.PENTAD_SAMPLES/.MONTHLY/.sst/dods",
                var_src="sst", k_offset=False,
            ),
            tref=dict(
                h_url=f"{NMME_BASE}/.NCEP-CFSv2/.HINDCAST/.PENTAD_SAMPLES/.MONTHLY/.tref/{_TREF_HIND_BUST}/pop/dods",
                f_url=f"{NMME_BASE}/.NCEP-CFSv2/.FORECAST/.PENTAD_SAMPLES/.MONTHLY/.tref/{_TREF_FCST_BUST}/pop/dods",
                var_src="tref", k_offset=False,
            ),
        ),
    ),
]

# Lookup by group name
MODEL_BY_GROUP = {m["group"]: m for m in MODELS}


def resolve_model(group: str, var: str) -> dict:
    """Flatten the registry for one (group, var) pair.

    Returns a dict with the same top-level keys the build/update scripts consume:
    group, label, combined, M_hind, M_fcst, L_max, chunks,
    h_url, f_url, var_src, k_offset, var_name, var_attrs.
    """
    m = MODEL_BY_GROUP[group]
    flat = {k: v for k, v in m.items() if k != "vars"}
    flat.update(m["vars"][var])
    flat["var_name"] = var
    flat["var_attrs"] = VARIABLES[var]["attrs"]
    return flat
