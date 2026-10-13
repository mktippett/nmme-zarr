#!/usr/bin/env python
"""
cfsv2_completeness.py — check NCEP-CFSv2 starts for missing members and repair.

CFSv2 is a lagged ensemble: runs start every 5th day at 00/06/12/18Z, so each
calendar month gets 24 members, except November which gets 28 (73 start dates
a year; see Saha et al. CFSv2 calibration-data document). The earliest lagged
members may lack the last lead (9.5), so a member counts as complete when it
has leads 0.5–8.5.

Short starts have two known causes, with different fixes:
  1. IRIDL's web cache holds a copy taken before the last members arrived.
     A cache-busted re-fetch (update_archive's --bust-cache mechanism) fixes it.
  2. IRIDL's downstream copy of the start was built before the original files
     were complete. Only a server-side refresh of that start fixes it; the
     command is kept in the gitignored code/iridl_private.py.

To keep the two separable, repair is staged per still-short start:
  stage A: cache-busted re-fetch → recheck
  stage B (--refresh, only if still short after A): server-side refresh →
           cache-busted re-fetch with a new token → recheck. IRIDL currently
           accepts the refresh only from a browser, so in practice stage B
           prints the URLs to open; rerunning --repair afterwards is stage B's
           re-fetch.
A start fixed in A was cache-only; one fixed only in B points at the
downstream copy. (Not proof — IRIDL could have refreshed on its own in
between — but the A-then-B order rules out the cache explanation.)

Usage
-----
  mamba run -n pangeo-local python code/cfsv2_completeness.py                 # report only
  mamba run -n pangeo-local python code/cfsv2_completeness.py --repair        # + stage A
  mamba run -n pangeo-local python code/cfsv2_completeness.py --repair --refresh  # + stage B
  mamba run -n pangeo-local python code/cfsv2_completeness.py --var tref --since 2021-01
"""

import argparse
import logging
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import xarray as xr
import zarr

sys.path.insert(0, str(Path(__file__).parent))
from nmme_models import NMME_BASE, VARIABLES, resolve_model
from iridl_io import bust_url
from update_archive import (get_remote_s, get_local_s, refetch_starts,
                            qa_member_consistency, _s_label)

try:
    from iridl_private import REFRESH_CMD
except ImportError:
    REFRESH_CMD = None

log = logging.getLogger(__name__)

GROUP = "NCEP-CFSv2"
N_MEMBERS = 24            # members per start month …
N_MEMBERS_NOV = 28        # … except November
REQUIRED_LEADS = slice(0.5, 8.5)   # 9.5 legitimately absent for early members
PROBE_X, PROBE_Y = 200.0, 0.0      # equatorial Pacific: ocean for sst, valid for tref
REFRESH_TIMEOUT = 1800             # s; the refresh returns silently when done
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


# ---------------------------------------------------------------------------
# Completeness check
# ---------------------------------------------------------------------------

def check_store(store_path: Path, var: str, since: float | None = None,
                s_values: np.ndarray | None = None) -> dict[float, str]:
    """Return {S value: description} for starts whose members are incomplete.

    Presence is tested at one probe grid point: IRIDL gaps are whole fields
    (a member or a member-lead), never partial maps.
    """
    ds = xr.open_zarr(str(store_path), group=GROUP, decode_times=False)
    da = ds[var].sel(X=PROBE_X, Y=PROBE_Y)
    if s_values is not None:
        da = da.sel(S=s_values)
    elif since is not None:
        da = da.sel(S=slice(since, None))
    present = da.notnull().load()                                   # (S, M, L)

    month = (present.S.astype(int) % 12) + 1
    n_exp = xr.where(month == 11, N_MEMBERS_NOV, N_MEMBERS)
    expected = present.M <= n_exp                                   # (S, M)
    complete = present.sel(L=REQUIRED_LEADS).all("L")              # (S, M)
    any_data = present.any("L")

    missing = expected & ~any_data           # expected member, no data at all
    partial = expected & any_data & ~complete
    bad_s = (missing | partial).any("M")

    short = {}
    for sv in present.S.where(bad_s, drop=True).values:
        miss_m = present.M.where(missing.sel(S=sv), drop=True).values.astype(int)
        parts = []
        if miss_m.size:
            parts.append(f"missing members {miss_m.tolist()}")
        for mv in present.M.where(partial.sel(S=sv), drop=True).values:
            gaps = present.L.where(~present.sel(S=sv, M=mv)
                                   .sel(L=REQUIRED_LEADS), drop=True).values
            parts.append(f"M={int(mv)} lacks leads {gaps.tolist()}")
        n_ok = int((expected & complete).sel(S=sv).sum())
        short[float(sv)] = (f"{n_ok}/{int(n_exp.sel(S=sv))} complete; "
                            + "; ".join(parts))
    return short


def report(short: dict[float, str], var: str, header: str) -> None:
    if not short:
        log.info("[%s] %s: all starts complete.", var, header)
        return
    log.warning("[%s] %s: %d incomplete start(s):", var, header, len(short))
    for sv, desc in short.items():
        log.warning("    %s (S=%d): %s", _s_label(sv), sv, desc)


# ---------------------------------------------------------------------------
# Repair
# ---------------------------------------------------------------------------

def _stream_base(var: str, forecast: bool) -> str:
    stream = "FORECAST" if forecast else "HINDCAST"
    return f"{NMME_BASE}/.NCEP-CFSv2/.{stream}/.PENTAD_SAMPLES/.MONTHLY/.{var}"


def refresh_url(var: str, sv: float, forecast: bool) -> str:
    y, mo = 1960 + int(sv) // 12, int(sv) % 12
    sel = f"S/(0000%201%20{MONTHS[mo]}%20{y})/VALUE"
    return f"{_stream_base(var, forecast)}/{sel}/{REFRESH_CMD}"


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


def server_refresh(url: str) -> str:
    """Issue the server-side refresh; return 'ok', 'login' or 'failed'.

    As of 2026-10 IRIDL redirects every non-browser client to /auth/login
    for this command ('login'); the caller then lists the URL for opening
    in a browser. Kept as a request so it starts working if that changes.
    """
    opener = urllib.request.build_opener(_NoRedirect)
    t0 = time.time()
    try:
        with opener.open(url, timeout=REFRESH_TIMEOUT) as resp:
            status = resp.status
    except urllib.error.HTTPError as exc:
        status, loc = exc.code, exc.headers.get("Location", "")
        if status in (301, 302, 303, 307) and "/auth/login" in loc:
            return "login"
        log.warning("    refresh failed: HTTP %d (%s)", status, url)
        return "failed"
    except Exception as exc:
        log.warning("    refresh failed: %s (%s)", exc, url)
        return "failed"
    log.info("    refresh HTTP %d after %.0f s", status, time.time() - t0)
    return "ok" if status == 200 else "failed"


def busted_refetch(store_path: Path, var: str, s_values: list[float]) -> None:
    """Re-fetch the given starts through a fresh cache-busting token."""
    m = dict(resolve_model(GROUP, var))
    token = str(int(time.time()))
    for key in ("h_url", "f_url"):
        m[key] = bust_url(m[key], token)
    s_hind, s_fcst = get_remote_s(m["h_url"]), get_remote_s(m["f_url"])
    local_s = get_local_s(store_path, GROUP)
    log.info("[%s] Re-fetching %d start(s) with cache-bust token %s",
             var, len(s_values), token)
    refetch_starts(store_path, m, np.array(s_values, dtype=np.float32), local_s,
                   s_hind, s_fcst, np.union1d(s_hind, s_fcst),
                   max(m["M_hind"], m["M_fcst"]))


def repair(store_path: Path, var: str, short: dict[float, str],
           do_refresh: bool, browser_urls: list[str]) -> dict[float, str]:
    """Staged repair; refresh URLs IRIDL would only take from a browser are
    appended to browser_urls."""
    s_list = list(short)
    touched = set(s_list)

    # Stage A: cache-busted re-fetch alone
    busted_refetch(store_path, var, s_list)
    after_a = check_store(store_path, var, s_values=np.array(s_list))
    for sv in s_list:
        if sv not in after_a:
            log.info("[%s] %s fixed by cache-busted re-fetch (stage A).",
                     var, _s_label(sv))
    report(after_a, var, "after stage A")
    remaining = after_a

    # Stage B: server-side refresh, then a re-fetch past the web cache
    if remaining and do_refresh:
        if REFRESH_CMD is None:
            log.warning("[%s] No code/iridl_private.py — skipping server-side "
                        "refresh. Ask IRIDL to refresh the starts above.", var)
        else:
            s_fcst = get_remote_s(resolve_model(GROUP, var)["f_url"])
            accepted = []
            for sv in remaining:
                url = refresh_url(var, sv, forecast=bool(np.isin(sv, s_fcst)))
                log.info("[%s] Server-side refresh of %s", var, _s_label(sv))
                outcome = server_refresh(url)
                if outcome == "ok":
                    accepted.append(sv)
                elif outcome == "login":
                    browser_urls.append(url)
            if accepted:
                busted_refetch(store_path, var, accepted)
                after_b = check_store(store_path, var,
                                      s_values=np.array(list(remaining)))
                for sv in accepted:
                    if sv not in after_b:
                        log.info("[%s] %s fixed after server-side refresh "
                                 "(stage B; was still short after stage A).",
                                 var, _s_label(sv))
                report(after_b, var, "after stage B")
                remaining = after_b

    m = resolve_model(GROUP, var)
    qa_member_consistency(store_path, m, np.array(sorted(touched), dtype=np.float32))
    zarr.consolidate_metadata(str(store_path))
    return remaining


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--var", nargs="+", default=list(VARIABLES),
                        choices=list(VARIABLES))
    parser.add_argument("--since", default=None,
                        help="Only check starts from YYYY-MM on (default: all)")
    parser.add_argument("--repair", action="store_true",
                        help="Re-fetch incomplete starts with cache busting (stage A)")
    parser.add_argument("--refresh", action="store_true",
                        help="With --repair: server-side refresh of starts still "
                             "incomplete after stage A, then re-fetch (stage B)")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    since = None
    if args.since:
        y, mo = map(int, args.since.split("-"))
        since = float((y - 1960) * 12 + mo - 1)

    browser_urls = []
    for var in args.var:
        store_path = VARIABLES[var]["default_store"]
        short = check_store(store_path, var, since=since)
        report(short, var, "store check")
        if short and args.repair:
            repair(store_path, var, short, do_refresh=args.refresh,
                   browser_urls=browser_urls)

    if browser_urls:
        log.warning("IRIDL redirected the server-side refresh to login. Open "
                    "each URL below in a browser (it returns silently when "
                    "done), then rerun with --repair to re-fetch and recheck:")
        print("\n".join(browser_urls))


if __name__ == "__main__":
    main()
