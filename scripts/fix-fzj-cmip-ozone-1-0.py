"""
Fix raw FZJ-CMIP-ozone-1-0 data
"""

from __future__ import annotations

import shutil
from pathlib import Path

import iris
import netCDF4
from input4mips_prep import (
    copy_attributes,
    copy_dimensions,
    copy_variable,
    reverse_variable_dimensions,
    rewrite_time_from_bounded_to_climatology,
)
from input4mips_validation.cvs import load_cvs
from input4mips_validation.inference.from_data import (
    infer_time_start_time_end_for_filename,
)
from input4mips_validation.validation.file import get_validate_file_result
from input4mips_validation.xarray_helpers.iris import ds_from_iris_cubes


def rewrite_clim_qbo_files(
    clim_qbo_files: list[Path],
    out_dir: Path,
    correct_contact: str,
    verbose: bool = True,
) -> list[Path]:
    """
    Rewrite the climatology files including the QBO signal
    """
    res = []
    for f in clim_qbo_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        ds_fixed = netCDF4.Dataset(fixed_file, "w")

        if verbose:
            print(f"Rewriting {f} to {fixed_file}")

        copy_attributes(ds, ds_fixed)
        ds_fixed.setncattr("contact", correct_contact)

        copy_dimensions(ds, ds_fixed)

        bounds_vars = [v for v in ds.variables if "bnds" in v]
        other_vars = [v for v in ds.variables if v not in bounds_vars]

        for variable in bounds_vars:
            reverse_variable_dimensions(ds, ds_fixed, variable)

        for variable in other_vars:
            copy_variable(ds, ds_fixed, variable)

        ds.close()
        ds_fixed.close()
        res.append(fixed_file)

    return res


def rewrite_clim_files(
    clim_files: list[Path],
    out_dir: Path,
    correct_contact: str,
    verbose: bool = True,
) -> list[Path]:
    """
    Rewrite the climatology files
    """
    res = []
    for f in clim_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        ds_fixed = netCDF4.Dataset(fixed_file, "w")

        if verbose:
            print(f"Rewriting {f} to {fixed_file}")

        copy_attributes(ds, ds_fixed)
        ds_fixed.setncattr("contact", correct_contact)

        copy_dimensions(ds, ds_fixed)

        bounds_vars = [v for v in ds.variables if "bnds" in v]
        other_vars = [v for v in ds.variables if v not in bounds_vars]

        for variable in [v for v in bounds_vars if "time" not in v]:
            reverse_variable_dimensions(ds, ds_fixed, variable)

        for variable in other_vars:
            copy_variable(ds, ds_fixed, variable)

        rewrite_time_from_bounded_to_climatology(
            ds_fixed, original_time_bounds=ds["time_bnds"]
        )

        ds.close()
        ds_fixed.close()
        res.append(fixed_file)

    return res


def rewrite_historical_files(
    historical_files: list[Path],
    out_dir: Path,
    correct_contact: str,
    verbose: bool = True,
) -> list[Path]:
    """
    Rewrite the historical files
    """
    res = []
    for f in historical_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        ds_fixed = netCDF4.Dataset(fixed_file, "w")

        if verbose:
            print(f"Rewriting {f} to {fixed_file}")

        copy_attributes(ds, ds_fixed)
        ds_fixed.setncattr("contact", correct_contact)

        copy_dimensions(ds, ds_fixed)
        bounds_vars = [v for v in ds.variables if "bnds" in v]
        other_vars = [v for v in ds.variables if v not in bounds_vars]

        for variable in bounds_vars:
            reverse_variable_dimensions(ds, ds_fixed, variable)

        for variable in other_vars:
            copy_variable(ds, ds_fixed, variable)

        ds.close()
        ds_fixed.close()
        res.append(fixed_file)

    return res


def rewrite_file_in_drs(
    file: Path, output_root: Path, cv_source: str, verbose: bool = True
) -> Path:
    """
    Re-write files in the DRS - ensures that filenames are correct too

    Also does some other checks implicitly
    """
    cvs = load_cvs(cv_source=cv_source)

    ds = ds_from_iris_cubes(
        iris.load(file),
        # xr_variable_processor=xr_variable_processor,
        raw_file=file,
        time_dimension="time",
    )

    time_start, time_end = infer_time_start_time_end_for_filename(
        ds=ds,
        frequency_metadata_key="frequency",
        no_time_axis_frequency="fx",
        time_dimension="time",
    )

    full_file_path = cvs.DRS.get_file_path(
        root_data_dir=output_root,
        available_attributes=ds.attrs,
        time_start=time_start,
        time_end=time_end,
    )

    full_file_path.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy(file, full_file_path)

    if verbose:
        print(f"File written according to the DRS in {full_file_path}")

    return full_file_path


def main() -> None:
    """
    Apply the fixes
    """
    REPO_ROOT = Path(__file__).parents[1]

    IN_DIR = REPO_ROOT / "../input4MIPs_CVs/fzj-ozone"
    TMP_DIR = REPO_ROOT / "../input4MIPs_CVs/fzj-ozone-tmp"
    OUT_FULL_TREE_DIR = REPO_ROOT / "../input4MIPs_CVs/fzj-ozone-rw-full-tree"
    OUT_DIR = REPO_ROOT / "../input4MIPs_CVs/fzj-ozone-rw"

    CV_SOURCE = "gh:main"
    SOURCE_ID = "FZJ-CMIP-ozone-1-0"

    cvs = load_cvs(cv_source=CV_SOURCE)

    files_to_rewrite = tuple(IN_DIR.rglob("*.nc"))
    clim_qbo_files = [v for v in files_to_rewrite if "qbo" in str(v)]
    clim_files = [
        v for v in files_to_rewrite if "qbo" not in str(v) and "clim" in str(v)
    ]
    historical_files = [
        f for f in files_to_rewrite if f not in [*clim_files, *clim_qbo_files]
    ]

    TMP_DIR.mkdir(exist_ok=True, parents=True)

    rewritten_clim_qbo_files = rewrite_clim_qbo_files(
        clim_qbo_files,
        out_dir=TMP_DIR,
        correct_contact=cvs.source_id_entries[SOURCE_ID].values.contact,
    )

    rewritten_clim_files = rewrite_clim_files(
        clim_files,
        out_dir=TMP_DIR,
        correct_contact=cvs.source_id_entries[SOURCE_ID].values.contact,
    )

    rewritten_historical_files = rewrite_historical_files(
        historical_files,
        out_dir=TMP_DIR,
        correct_contact=cvs.source_id_entries[SOURCE_ID].values.contact,
    )

    OUT_DIR.mkdir(exist_ok=True, parents=True)
    for tmp_file in [
        *rewritten_clim_qbo_files,
        *rewritten_clim_files,
        *rewritten_historical_files,
    ]:
        renamed_file = rewrite_file_in_drs(
            tmp_file, output_root=OUT_FULL_TREE_DIR, cv_source=CV_SOURCE
        )
        out_file = OUT_DIR / Path(renamed_file).name
        shutil.copy(renamed_file, out_file)

        validation_res = get_validate_file_result(
            out_file,
            cv_source=CV_SOURCE,
            # xr_variable_processor=xr_variable_processor,
            # frequency_metadata_keys=frequency_metadata_keys,
            # bounds_info=bounds_info,
            # There is an issue with the data and the bounds.
            # The data's lat co-ordinate is
            # `[-90, -87.5, ... 87.5, 90]`.
            # As a result, it isn't clear what the cell bounds are.
            # They can't be e.g.
            # `[[-90, -87.5], [-87.5, -85], ..., [87.5, 90], [90, 92.5]]`
            # because there is no latitude greater than 90.
            # As a result, this data fails
            # (rightly so, it's unclear what this data represents).
            # allow_cf_checker_warnings=False,
            allow_cf_checker_warnings=True,
        )
        validation_res.raise_if_errors()


if __name__ == "__main__":
    main()
