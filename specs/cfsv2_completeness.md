# Spec: cfsv2_completeness.py

## 1. Purpose

Check that every NCEP-CFSv2 start in the sst and tref stores holds its full
lagged ensemble, and repair starts that IRIDL delivered incomplete, keeping
the two known causes (web-cache staleness vs. a stale IRIDL downstream copy)
separately attributable.

## 2. Inputs

| File / source | Relevant content | Filters |
|---|---|---|
| `data/nmme_{sst,tref}.zarr`, group `NCEP-CFSv2` | `{var}(S, M, L, Y, X)` | one probe point X=200, Y=0; optional `--since YYYY-MM` |
| IRIDL OPeNDAP `NCEP-CFSv2/.{HINDCAST,FORECAST}/.PENTAD_SAMPLES/.MONTHLY/.{var}` | full fields for the flagged starts | cache-busted URL (fresh `/<token>/pop/`) |
| `code/iridl_private.py` (gitignored, optional) | `REFRESH_CMD`: private IRIDL command that rebuilds IRIDL's downstream copy of one start | absent → stage B reports only |

## 3. Outputs

| Output | Contents | Format |
|---|---|---|
| Log (stderr) | per incomplete start: `n_ok/n_expected complete`, missing members, members lacking required leads; per stage results; QA flags | text |
| Zarr stores (with `--repair`) | flagged starts overwritten in place with a cache-busted re-fetch; metadata re-consolidated | zarr v3 |
| stdout (with `--refresh`) | refresh URLs IRIDL redirected to login, one per line, to open in a browser | text |

## 4. Algorithm

1. **Check.** For each start, `present(S,M,L) = notnull` at the probe point.
   Expected members `n_exp = 28` if the start month is November, else 24;
   member `M` is expected iff `M ≤ n_exp`. A member is *complete* iff it has
   data at every lead in 0.5–8.5 (9.5 optional). A start is incomplete if any
   expected member has no data at all (*missing*) or has data but is not
   complete (*partial*). Members `> n_exp` are not checked.
2. **Stage A** (`--repair`): re-fetch all incomplete starts through a fresh
   cache-bust token (`bust_url`, same mechanism as `update_archive.py
   --bust-cache`) via `update_archive.refetch_starts`, then re-check. Starts
   now complete are logged "fixed by cache-busted re-fetch (stage A)".
3. **Stage B** (`--repair --refresh`, only for starts still incomplete): build
   `…/.{stream}/.PENTAD_SAMPLES/.MONTHLY/.{var}/S/(0000 1 Mon YYYY)/VALUE/<REFRESH_CMD>`
   (stream = FORECAST if the start is in the remote forecast S axis, else
   HINDCAST) and request it without following redirects. HTTP 200 → accepted;
   redirect to `/auth/login` → URL appended to the browser list. Accepted
   starts are re-fetched with a new token and re-checked; those now complete
   are logged "fixed after server-side refresh (stage B)". When IRIDL
   redirects (the current behavior), the browser URLs are printed at the end;
   opening them and re-running `--repair` completes stage B by hand.
4. Run `update_archive.qa_member_consistency` on all touched starts, then
   `zarr.consolidate_metadata`.

Attribution logic: a start that survives stage A and is fixed only after the
refresh cannot be explained by web-cache staleness. This is not proof of
causation (IRIDL could refresh on its own in between), only the ordering
that rules out the cache explanation.

## 5. Constants & Scientific Rationale

| Name | Value | Why |
|---|---|---|
| `N_MEMBERS`, `N_MEMBERS_NOV` | 24, 28 | CFSv2 runs start every 5th day at 00/06/12/18Z; the monthly lagged ensemble is 24 members, November 28 (Saha, van den Dool & Johansson, *CFSv2 Retrospective Forecast Calibration Climatologies*, 2011: "an ensemble size of 24 forecasts for each month, except November which has 28"). Confirmed in the store: sst/tref Novembers carry 28 after repair. |
| `REQUIRED_LEADS` | 0.5–8.5 | In forecast starts members 1–12 have 9 leads and 13–24 have 10; the missing lead 9.5 for the earliest lagged members is by design (MKT, 2026-10-13), so 9.5 is not required. |
| `PROBE_X, PROBE_Y` | 200, 0 | Equatorial Pacific: ocean, so valid for sst as well as tref. One point suffices because IRIDL gaps observed so far are whole fields (a member or a member-lead). Chunks are full-globe, so the probe does not reduce I/O; it keeps the result small. |
| `REFRESH_TIMEOUT` | 1800 s | Implementation choice (the refresh returns silently only when done); no measured run time yet. |

## 6. Edge Cases & Error Handling

- **Refresh is browser-only (2026-10-13).** Python and curl requests to the
  refresh URL are redirected to `/auth/login` regardless of user agent or
  trailing slash, from both the Claude session and the user's own shell.
  The browser succeeds without a deliberate login (presumably a session
  cookie). The script still tries the request so it works if IRIDL changes.
- **Pasting long shell one-liners breaks** (terminal wraps split the command);
  hence the URLs are printed one per line rather than as a loop.
- **Single member-lead holes** that survive a cache-busted re-fetch are most
  likely gaps in IRIDL's original files; leave them as NaN (fix-at-source
  discipline: do not mask downstream). See Open item.
- Gitignored private file absent → stage B logs a warning and only reports.
- The command name is never written into committed files.

### Open item (2026-10-13)

Nine starts hold every member but lack one or two member-leads, unchanged
under repeated cache-busted re-fetches: sst 2018-06 (M16 L7.5), 2020-12
(M1 L5.5), 2026-05 (M12 L8.5); tref 1984-11 (M11 L0.5), 1993-03 (M17,M18
L3.5), 2019-08 (M17 L1.5), 2020-12 (M1 L5.5), 2021-03 (M4 L0.5), 2025-01
(M23 L2.5). None has had a refresh that actually ran (the user's browser
refresh of tref 2021-03 still showed the hole, but timing vs. the store is
not established). **Decision:** run `--repair --refresh`, open the printed
URLs in a browser, re-run `--repair`; any hole still present gets listed in
CLAUDE.md "Expected NaN" as an upstream gap.

## 7. Synchronization Log

| Date | Code change | Spec updated |
|---|---|---|
| 2026-10-13 | Initial implementation. Motivated by 39 tref starts (2022-08, 2023-01→2024-12, 2025-08→2026-09) and sst 2026-10 holding 16–20 of 24 members from stale IRIDL caches, fixed by cache-busted re-fetch; sst 2026-08 needed a browser refresh. | 2026-10-13 |

## Verification snippet

```python
# Run after a --repair pass: no expected member may be entirely missing.
import sys; sys.path.insert(0, "code")
from cfsv2_completeness import check_store
from nmme_models import VARIABLES
for var in VARIABLES:
    short = check_store(VARIABLES[var]["default_store"], var)
    missing = {s: d for s, d in short.items() if "missing members" in d}
    assert not missing, f"{var}: starts with missing members: {missing}"
    print(f"{var}: {len(short)} start(s) with member-lead holes only")
```
