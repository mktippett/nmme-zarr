"""Extract seasonal-mean NMME forecast anomalies for one start month to netCDF.

For each model in the store: average the requested season's leads for every
member, subtract the model's lead-dependent climatology (mean over all
hindcast members and all same-calendar-month starts in the climatology
period), and stack models along a `model` dimension (member padded with NaN).

Writes one netCDF (data/requests/) and one 4x2 map figure (plots/requests/)
per requested region.

Example (Oct 2026 start, DJF, 1991-2020 climatology, global + North America):
    mamba run -n pangeo-local python code/extract_seasonal_anomaly.py \
        --var tref --start 2026-10 --season DJF --regions global namerica
"""

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
MONTHS = "JFMAMJJASONDJFMAMJJASOND"

# lat/lon bounds (degrees_north, degrees_east on the store's 0-359 grid);
# None = full globe. "suffix" is appended to output filenames.
REGIONS = {
    "global": {"lat": None, "lon": None, "suffix": "", "label": ""},
    "namerica": {"lat": (10, 75), "lon": (190, 310), "suffix": "_namerica",
                 "label": "N. America "},
}


def season_months(season: str) -> list[int]:
    """'DJF' -> [12, 1, 2] (1-based calendar months)."""
    i = MONTHS.find(season.upper())
    if i < 0:
        raise ValueError(f"unrecognized season {season!r}")
    return [(i + k) % 12 + 1 for k in range(len(season))]


def s_value(year: int, month: int) -> float:
    """Start coordinate S in months since 1960-01."""
    return float((year - 1960) * 12 + (month - 1))


def model_anomaly(ds, var, start_year, start_month, leads, clim_years):
    """Seasonal-mean anomaly (member, Y, X) and climatology (Y, X) for one model."""
    clim_S = [s_value(y, start_month) for y in range(clim_years[0], clim_years[1] + 1)]
    fcst_S = s_value(start_year, start_month)
    missing = sorted(set(clim_S + [fcst_S]) - set(ds.S.values.tolist()))
    if missing:
        raise ValueError(f"starts missing from store: S={missing}")

    da = ds[var].sel(S=clim_S + [fcst_S], L=leads).load()
    # skipna=False: a member missing any lead of the season gets no seasonal mean
    seas = da.mean("L", skipna=False)
    clim = seas.sel(S=clim_S).mean(["S", "M"])
    anom = seas.sel(S=fcst_S, drop=True) - clim
    n_clim = int(seas.sel(S=clim_S).notnull().any(["Y", "X"]).sum())
    return anom, clim, n_clim


def subset_region(ds, region):
    r = REGIONS[region]
    if r["lat"] is None:
        return ds
    return ds.sel(lat=slice(*r["lat"]), lon=slice(*r["lon"])).assign_attrs(
        region=region, region_lat=str(r["lat"]), region_lon=str(r["lon"]))


def plot_anomalies(ds_out, var, region, out_path):
    """4x2 maps: ensemble-mean anomaly for each model plus the multi-model mean.

    Takes the global ds_out and clips the view to the region, so a projected
    (non-rectangular) lat/lon box is filled to the frame edge.
    """
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    ens = ds_out[f"{var}_anom"].mean("member")
    panels = [(ens.sel(model=m), f"{m} ({int(ds_out.n_members.sel(model=m))})")
              for m in ens.model.values]
    panels.append((ds_out[f"{var}_anom_mme"], "Multi-model mean"))

    r = REGIONS[region]
    if r["lat"] is None:
        proj, figsize = ccrs.Robinson(central_longitude=180), (11, 11.5)
    else:
        proj, figsize = ccrs.LambertConformal(central_longitude=-100, central_latitude=45), (9, 13)

    levels = np.arange(-4, 4.25, 0.5)
    fig, axes = plt.subplots(4, 2, figsize=figsize, subplot_kw={"projection": proj},
                             layout="constrained")
    for ax, (da, title) in zip(axes.flat, panels):
        cf = da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels,
                              cmap="RdBu_r", extend="both", add_colorbar=False)
        ax.coastlines(linewidth=0.5)
        if r["lat"] is None:
            ax.set_global()
        else:
            import cartopy.feature as cfeature
            ax.add_feature(cfeature.BORDERS, linewidth=0.4)
            ax.add_feature(cfeature.STATES, linewidth=0.2, edgecolor="0.4")
            lon0, lon1 = (x - 360 for x in r["lon"])
            ax.set_extent([lon0, lon1, *r["lat"]], crs=ccrs.PlateCarree())
        ax.set_title(title, fontsize=10)
    fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.6, aspect=40,
                 label=f"{var} anomaly ({ds_out[f'{var}_anom'].attrs['units']})")
    fig.suptitle(f"NMME {r['label']}{var} {ds_out.attrs['season']} ensemble-mean anomaly, "
                 f"{ds_out.attrs['forecast_start']} start "
                 f"(clim {ds_out.attrs['clim_period']}; members in parentheses)",
                 fontsize=11)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--var", default="tref", choices=["tref", "sst"])
    p.add_argument("--start", required=True, help="forecast start, YYYY-MM")
    p.add_argument("--season", required=True, help="consecutive-month season, e.g. DJF")
    p.add_argument("--clim", default="1991-2020", help="climatology start years, YYYY-YYYY")
    p.add_argument("--regions", nargs="+", default=["global"], choices=list(REGIONS))
    args = p.parse_args()

    start_year, start_month = map(int, args.start.split("-"))
    clim_years = tuple(map(int, args.clim.split("-")))
    months = season_months(args.season)
    leads = [((m - start_month) % 12) + 0.5 for m in months]
    if np.diff(leads).min() != 1:
        raise ValueError(f"season {args.season} wraps past start month {start_month}")

    store = ROOT / "data" / f"nmme_{args.var}.zarr"
    dt = xr.open_datatree(store, engine="zarr", decode_times=False)

    anoms, clims, n_mem, n_clim = [], [], [], []
    for name in dt.children:
        ds = dt[name].to_dataset()
        if not set(leads) <= set(ds.L.values.tolist()):
            print(f"{name}: leads {leads} not all available, skipping")
            continue
        anom, clim, nc = model_anomaly(ds, args.var, start_year, start_month, leads, clim_years)
        nm = int(anom.notnull().any(["Y", "X"]).sum())
        print(f"{name}: {nm} forecast members, {nc} hindcast member-starts in climatology")
        anoms.append(anom.expand_dims(model=[name]))
        clims.append(clim.expand_dims(model=[name]))
        n_mem.append(nm)
        n_clim.append(nc)

    # outer join on M pads members a model lacks with NaN
    anom = xr.concat(anoms, "model", join="outer").rename(M="member", Y="lat", X="lon")
    clim = xr.concat(clims, "model").rename(Y="lat", X="lon")
    models = anom.model.values
    units = dt[models[0]][args.var].attrs.get("units", "")
    long_name = dt[models[0]][args.var].attrs.get("long_name", args.var)

    ds_out = xr.Dataset(
        {
            f"{args.var}_anom": anom.astype("float32").assign_attrs(
                units=units, long_name=f"{args.season} mean {long_name} anomaly"),
            f"{args.var}_clim": clim.astype("float32").assign_attrs(
                units=units,
                long_name=f"{args.season} mean {long_name} climatology, "
                          f"{clim_years[0]}-{clim_years[1]} hindcasts"),
            f"{args.var}_anom_mme": anom.mean("member").mean("model").astype("float32").assign_attrs(
                units=units,
                long_name=f"{args.season} mean {long_name} anomaly, multi-model mean "
                          "(equal weight per model ensemble mean)"),
            "n_members": xr.DataArray(n_mem, coords={"model": models}, dims="model",
                                      attrs={"long_name": "forecast members present"}),
            "n_clim_samples": xr.DataArray(n_clim, coords={"model": models}, dims="model",
                                           attrs={"long_name": "member-starts in climatology"}),
        },
        attrs={
            "title": f"NMME {args.var} {args.season} anomalies, {args.start} start",
            "source": "NMME via IRIDL; local archive " + store.name,
            "forecast_start": args.start,
            "season": args.season,
            "clim_period": args.clim,
            "leads": str(leads),
            "climatology": (f"per model and grid point: mean over all members and "
                            f"start years {args.clim} with start month {start_month:02d}"),
            "anomaly": "each forecast member's seasonal mean minus its model's climatology",
        },
    )
    ds_out.lat.attrs.update(units="degrees_north")
    ds_out.lon.attrs.update(units="degrees_east")

    stem = f"nmme_{args.var}_anom_{args.season}_{args.start}_start_clim{args.clim}"
    for region in args.regions:
        ds_reg = subset_region(ds_out, region)
        name = stem + REGIONS[region]["suffix"]
        out = ROOT / "data" / "requests" / f"{name}.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        enc = {v: {"zlib": True, "complevel": 4} for v in ds_reg.data_vars}
        ds_reg.to_netcdf(out, encoding=enc)
        print(f"wrote {out}")
        plot_anomalies(ds_out, args.var, region, ROOT / "plots" / "requests" / f"{name}.png")


if __name__ == "__main__":
    main()
