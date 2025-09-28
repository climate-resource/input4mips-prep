"""
Re-write emissions files received from IIASA (Jarmo :))
"""

from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np
from input4mips_prep import (
    copy_attributes,
    copy_dimension,
    copy_dimensions,
    copy_variable,
)
from input4mips_validation.io import (
    generate_creation_timestamp,
    generate_tracking_id,
)


def add_correct_time_bounds(
    ds: netCDF4.Dataset,
    ds_fixed: netCDF4.Dataset,
    bounds_time: str = "time_bnds",
) -> netCDF4.Dataset:
    """
    Add the correct time bounds

    Parameters
    ----------
    ds
        Dataset from which to get the original data

    ds_fixed
        Fixed dataset from which to get the correct time bounds

    bounds_time
        Name of the variable to use for time bounds

    Returns
    -------
    :
        `ds_fixed` with time bounds applied

        Note that the operation is in place,
        so returning this is just a convenience.
    """
    if "bound" not in ds.dimensions:
        raise AssertionError

    ds_fixed.createVariable(bounds_time, ds["time"].datatype, ("time", "bound"))

    time_int = ds["time"][:]
    time_date = netCDF4.num2date(
        time_int,
        units=ds["time"].getncattr("units"),
        calendar=ds["time"].getncattr("calendar"),
    )
    bounds_time_fixed_l = []
    # Loop, crazy slow, whatever
    n_months_in_year = 12
    for v_int, v_date in zip(time_int, time_date):
        start_of_bound_date = type(v_date)(v_date.year, v_date.month, 1)
        if v_date.month == n_months_in_year:
            end_of_bound_date = type(v_date)(v_date.year + 1, 1, 1)

        else:
            end_of_bound_date = type(v_date)(v_date.year, v_date.month + 1, 1)

        bounds_time_fixed_l.append(
            [
                netCDF4.date2num(
                    v,
                    units=ds["time"].getncattr("units"),
                    calendar=ds["time"].getncattr("calendar"),
                )
                for v in [start_of_bound_date, end_of_bound_date]
            ]
        )

    bounds_time_fixed = np.array(bounds_time_fixed_l)
    ds_fixed.variables[bounds_time][:] = bounds_time_fixed

    return ds_fixed


def rewrite_anthro_file(
    fp: Path,
    out_dir: Path,
    verbose: bool = True,
):
    ds = netCDF4.Dataset(fp)

    fixed_file = out_dir / fp.name.replace("_gn", "-0-3-0_gn")
    ds_fixed = netCDF4.Dataset(fixed_file, "w")

    if verbose:
        print(f"Rewriting {fp} to {fixed_file}")

    copy_attributes(ds, ds_fixed)
    ds_fixed.setncattr("source_version", "0.3.0")
    ds_fixed.setncattr("contact", "kikstra@iiasa.ac.at")
    ds_fixed.setncattr("creation_date", generate_creation_timestamp())
    ds_fixed.setncattr("tracking_id", generate_tracking_id())
    # Urgh out of date cf-checker
    ds_fixed.setncattr("Conventions", "CF-1.8")

    for dimension in ds.dimensions:
        if dimension != "sector":
            copy_dimension(ds, ds_fixed, dimension)

    in_sectors = list(ds["sector"][:])
    assert in_sectors == [  # noqa: S101
        "Agriculture",
        "Energy",
        "Industrial",
        "Transportation",
        "Residential, Commercial, Other",
        "Solvents Production and Application",
        "Waste",
        "International Shipping",
        "Other Capture and Removal",
    ], "Have to shuffle order too probably to avoid confusing people"

    ds_fixed.createDimension("sector", len(ds["sector"][:]))
    ds_fixed.createVariable("sector", int, ("sector",))
    ds_fixed["sector"][:] = np.arange(len(ds["sector"][:]))
    ds_fixed["sector"].setncattr("long_name", "sector")
    ds_fixed["sector"].setncattr("bounds", "sector_bnds")
    ds_fixed["sector"].setncattr(
        "ids", "; ".join(f"{i}: {sector}" for i, sector in enumerate(in_sectors))
    )

    ds_fixed.createVariable("sector_bnds", float, ("sector", "bound"))
    ds_fixed["sector_bnds"][:] = np.vstack(
        [
            np.arange(len(ds["sector"][:])) - 0.5,
            np.arange(len(ds["sector"][:])) + 0.5,
        ]
    ).T

    bounds_vars_except_time = [
        v for v in ds.variables if "bnds" in v and "time" not in v
    ]
    bounds_time = "time_bnds"
    other_vars = [
        v for v in ds.variables if v not in [*bounds_vars_except_time, bounds_time]
    ]

    for variable in bounds_vars_except_time:
        copy_variable(
            ds,
            ds_fixed,
            variable,
            # Value taken from variable instead rather than being repeated
            copy_fill_value=False,
            copy_attrs=False,
        )

    ds_fixed = add_correct_time_bounds(ds, ds_fixed)

    for variable in other_vars:
        if variable == "sector":
            continue

        if verbose:
            print(f"Copying {variable}")

        copy_variable(ds, ds_fixed, variable)

    v_name = fp.name.split("_")[0]
    v_name = ds.getncattr("variable_id")
    ds_fixed[v_name].setncattr(
        "long_name", f"{v_name.split('_')[0]} anthropogenic emissions"
    )
    # TODO: check this
    ds_fixed[v_name].setncattr("units", "kg m-2 s-1")

    ds.close()
    ds_fixed.close()


def rewrite_air_anthro_file(
    fp: Path,
    out_dir: Path,
    verbose: bool = True,
):
    ds = netCDF4.Dataset(fp)

    fixed_file = out_dir / fp.name.replace("_gn", "-0-3-0_gn")
    ds_fixed = netCDF4.Dataset(fixed_file, "w")

    if verbose:
        print(f"Rewriting {fp} to {fixed_file}")

    copy_attributes(ds, ds_fixed)
    ds_fixed.setncattr("source_version", "0.3.0")
    ds_fixed.setncattr("contact", "kikstra@iiasa.ac.at")
    ds_fixed.setncattr("creation_date", generate_creation_timestamp())
    ds_fixed.setncattr("tracking_id", generate_tracking_id())

    copy_dimensions(ds, ds_fixed)

    bounds_vars_except_time = [
        v for v in ds.variables if "bnds" in v and "time" not in v
    ]
    bounds_time = "time_bnds"
    other_vars = [
        v for v in ds.variables if v not in [*bounds_vars_except_time, bounds_time]
    ]

    for variable in bounds_vars_except_time:
        copy_variable(
            ds,
            ds_fixed,
            variable,
            # Value taken from variable instead rather than being repeated
            copy_fill_value=False,
            copy_attrs=False,
        )

    ds_fixed = add_correct_time_bounds(ds, ds_fixed)

    # Add level bounds
    if "bound" not in ds.dimensions:
        raise AssertionError

    level_diff = np.diff(ds["level"][:])
    # Check levels all same size,
    # if not, logic below needs updating
    np.testing.assert_allclose(level_diff, level_diff[0], atol=1e-10)
    level_lower = np.round(ds["level"] - level_diff[0] / 2, 8)
    level_lower[level_lower == 0.0] = np.abs(level_lower[level_lower == 0.0])
    level_upper = level_lower + level_diff[0]
    level_bnds_vals = np.vstack([level_lower, level_upper]).T

    ds_fixed.createVariable("level_bnds", ds["level"].datatype, ("level", "bound"))
    ds_fixed.variables["level_bnds"][:] = level_bnds_vals

    for variable in other_vars:
        if verbose:
            print(f"Copying {variable}")

        copy_variable(ds, ds_fixed, variable)

    ds.close()
    ds_fixed.close()


def main():
    """
    Re-write the files
    """
    SOURCE_DIR = Path("../input4MIPs_CVs/iiasa-20250928")
    TMP_DIR = Path("../input4MIPs_CVs/iiasa-20250928-rewrite")
    OUT_DIR = Path("../input4MIPs_CVs/iiasa-20250928-drs")

    TMP_DIR.mkdir(exist_ok=True, parents=True)

    # for f in SOURCE_DIR.glob("*AIR-anthro*.nc"):
    #     rewrite_air_anthro_file(f, TMP_DIR)

    for f in SOURCE_DIR.glob("*em-anthro*.nc"):
        rewrite_anthro_file(f, TMP_DIR)


if __name__ == "__main__":
    main()
