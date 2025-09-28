"""
Re-write emissions files received from IIASA (Jarmo :))
"""

from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np
from input4mips_prep import copy_attributes, copy_dimensions, copy_variable
from input4mips_validation.io import (
    generate_creation_timestamp,
    generate_tracking_id,
)


def rewrite_anthro_file(fp: Path, out_dir: Path):
    ds = netCDF4.Dataset(fp)
    breakpoint()

    ds.close()


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

    # Fix the time bounds
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
    for v_int, v_date in zip(time_int, time_date):
        start_of_bound_date = type(v_date)(v_date.year, v_date.month, 1)
        if v_date.month == 12:
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
    for f in SOURCE_DIR.glob("*AIR-anthro*.nc"):
        rewrite_air_anthro_file(f, TMP_DIR)

    # for f in SOURCE_DIR.glob("*em-anthro*.nc"):
    #     rewrite_anthro_file(f, TMP_DIR)


if __name__ == "__main__":
    main()
