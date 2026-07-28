"""
Re-write vl-cf emissions files
"""

from __future__ import annotations

import shutil
from pathlib import Path

import netCDF4
import numpy as np


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

    # ds_fixed.createVariable(bounds_time, ds["time"].datatype, ("time", "bound"))

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


def rewrite_file(
    fp: Path,
    out_dir: Path,
    verbose: bool = True,
):
    """
    Re-write anthropogenic emissions file

    Parameters
    ----------
    fp
        File path to rewrite

    out_dir
        Directory in which to write the rewritten file

    verbose
        Should lots of information be printed?
    """
    fixed_file = out_dir / fp.name
    if verbose:
        print(f"Rewriting {fp} to {fixed_file}")

    shutil.copy2(fp, fixed_file)
    ds = netCDF4.Dataset(fixed_file, "a")
    ds_original = netCDF4.Dataset(fixed_file, "r")

    ds = add_correct_time_bounds(ds_original, ds)

    ds["time"].setncattr("bounds", "time_bnds")
    ds["time_bnds"].delncattr("units")
    ds["time_bnds"].delncattr("calendar")

    ds.close()
    ds_original.close()


def main():
    """
    Re-write the files
    """
    SOURCE_DIR = Path("../input4MIPs_CVs/vl-cf-emms")
    TMP_DIR = Path("../input4MIPs_CVs/vl-cf-emms-clean")
    # OUT_DIR = Path("../input4MIPs_CVs/vl-cf-emms-rw")

    # CV_SOURCE = "../input4MIPs_CVs"

    TMP_DIR.mkdir(exist_ok=True, parents=True)

    for f in SOURCE_DIR.glob("*.nc"):
        rewrite_file(f, TMP_DIR)

    # OUT_DIR.mkdir(exist_ok=True, parents=True)
    # for tmp_file in TMP_DIR.glob("*.nc"):
    #     out_file = rewrite_file_in_drs(
    #         tmp_file, output_root=OUT_DIR, cv_source=CV_SOURCE
    #     )
    #
    #     validation_res = get_validate_file_result(
    #         out_file,
    #         cv_source=CV_SOURCE,
    #         # xr_variable_processor=xr_variable_processor,
    #         # frequency_metadata_keys=frequency_metadata_keys,
    #         # bounds_info=bounds_info,
    #         # Allow issues with units for sector
    #         allow_cf_checker_warnings=True,
    #     )
    #     validation_res.raise_if_errors()


if __name__ == "__main__":
    main()
