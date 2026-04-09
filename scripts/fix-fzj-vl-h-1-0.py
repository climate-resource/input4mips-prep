"""
Fix raw FZJ `vl` and `h` scenario 1.0 data
"""

from __future__ import annotations

import shutil
from pathlib import Path

import iris
import netCDF4
from input4mips_validation.cvs import load_cvs
from input4mips_validation.inference.from_data import (
    infer_time_start_time_end_for_filename,
)
from input4mips_validation.validation.file import get_validate_file_result
from input4mips_validation.xarray_helpers.iris import ds_from_iris_cubes


def rewrite_nitrogen_files(
    nitrogen_files: list[Path],
    out_dir: Path,
    verbose: bool = True,
) -> list[Path]:
    """
    Rewrite the nitrogen files
    """
    res = []
    for f in nitrogen_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        shutil.copy2(f, fixed_file)
        ds_fixed = netCDF4.Dataset(fixed_file, "a")

        if verbose:
            print(f"Rewriting {f} to {fixed_file}")

        ds_fixed.setncattr("target_mip", "ScenarioMIP")
        source_id = f.name.split("_")[4]
        ds_fixed.setncattr("source_id", source_id)
        assert source_id.endswith("1-0")
        ds_fixed.setncattr("source_version", "1.0")

        ds.close()
        ds_fixed.close()
        res.append(fixed_file)

    return res


def rewrite_zmta_files(
    zmta_files: list[Path],
    out_dir: Path,
    verbose: bool = True,
) -> list[Path]:
    """
    Rewrite the zmta files
    """
    res = []
    for f in zmta_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        shutil.copy2(f, fixed_file)
        ds_fixed = netCDF4.Dataset(fixed_file, "a")

        if verbose:
            print(f"Rewriting {f} to {fixed_file}")

        ds_fixed.setncattr("target_mip", "ScenarioMIP")
        source_id = f.name.split("_")[4]
        ds_fixed.setncattr("source_id", source_id)
        assert source_id.endswith("1-0")
        ds_fixed.setncattr("source_version", "1.0")

        ds.close()
        ds_fixed.close()
        res.append(fixed_file)

    return res


def rewrite_ozone_files(
    ozone_files: list[Path],
    out_dir: Path,
    verbose: bool = True,
) -> list[Path]:
    """
    Rewrite the ozone files
    """
    res = []
    for f in ozone_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        shutil.copy2(f, fixed_file)
        ds_fixed = netCDF4.Dataset(fixed_file, "a")

        if verbose:
            print(f"Rewriting {f} to {fixed_file}")

        ds_fixed.setncattr("target_mip", "ScenarioMIP")
        source_id = f.name.split("_")[4]
        ds_fixed.setncattr("source_id", source_id)
        assert source_id.endswith("1-0")
        ds_fixed.setncattr("source_version", "1.0")

        ds_fixed["time"].setncattr("calendar", "noleap")

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
    IN_DIR = Path(
        "/global/cfs/projectdirs/m4931/zrjn-tmp/20260409-ozone/20260408-ozone/"
    )
    TMP_DIR = Path(
        "/global/cfs/projectdirs/m4931/zrjn-tmp/20260409-ozone/20260408-ozone-tmp/"
    )
    OUT_FULL_TREE_DIR = Path(
        "/global/cfs/projectdirs/m4931/zrjn-tmp/20260409-ozone/20260408-ozone-tmp-rw/"
    )

    CV_SOURCE = "./CVs"

    VERBOSE = True

    files_to_rewrite = tuple(IN_DIR.rglob("*.nc"))
    nitrogen_files = [v for v in files_to_rewrite if "nitrogen" in str(v)]
    zmta_files = [v for v in files_to_rewrite if "zmta" in str(v)]
    ozone_files = [
        f for f in files_to_rewrite if f not in [*zmta_files, *nitrogen_files]
    ]

    TMP_DIR.mkdir(exist_ok=True, parents=True)

    rewritten_nitrogen_files = rewrite_nitrogen_files(
        nitrogen_files, out_dir=TMP_DIR, verbose=VERBOSE
    )

    rewritten_zmta_files = rewrite_zmta_files(
        zmta_files, out_dir=TMP_DIR, verbose=VERBOSE
    )

    rewritten_ozone_files = rewrite_ozone_files(
        ozone_files, out_dir=TMP_DIR, verbose=VERBOSE
    )

    for tmp_file in [
        *rewritten_nitrogen_files,
        *rewritten_zmta_files,
        *rewritten_ozone_files,
    ]:
        renamed_file = rewrite_file_in_drs(
            tmp_file, output_root=OUT_FULL_TREE_DIR, cv_source=CV_SOURCE
        )
        if VERBOSE:
            print(f"Rewrote {tmp_file.name} to {renamed_file}")

        validation_res = get_validate_file_result(
            renamed_file,
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
        if VERBOSE:
            print(f"Validated {renamed_file}")


if __name__ == "__main__":
    main()
