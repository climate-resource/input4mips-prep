"""
Fix raw FZJ-CMIP-ozone-1-0 data
"""

from __future__ import annotations

from pathlib import Path

import netCDF4
from input4mips_prep import copy_attributes, copy_dimensions, copy_variables
from input4mips_validation.cvs import load_cvs
from input4mips_validation.validation.file import get_validate_file_result


def rewrite_historical_files(
    historical_files: list[Path],
    out_dir: Path,
    correct_contact: str,
) -> list[Path]:
    """
    Rewrite the historical files
    """
    res = []
    for f in historical_files:
        ds = netCDF4.Dataset(f, "r")

        fixed_file = out_dir / f.name
        ds_fixed = netCDF4.Dataset(fixed_file, "w")

        copy_attributes(ds, ds_fixed)
        ds_fixed.setncattr("contact", correct_contact)

        copy_dimensions(ds, ds_fixed)
        # breakpoint()
        copy_variables(ds, ds_fixed)

        ds.close()
        ds_fixed.close()
        res.append(fixed_file)

    return res


def main() -> None:
    """
    Apply the fixes
    """
    REPO_ROOT = Path(__file__).parents[1]

    IN_DIR = REPO_ROOT / "../input4MIPs_CVs/fzj-ozone"
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

    OUT_DIR.mkdir(exist_ok=True, parents=True)

    rewritten_historical_files = rewrite_historical_files(
        historical_files,
        out_dir=OUT_DIR,
        correct_contact=cvs.source_id_entries[SOURCE_ID].values.contact,
    )

    for rewritten_file in [*rewritten_historical_files]:
        validation_res = get_validate_file_result(
            rewritten_file,
            cv_source=CV_SOURCE,
            # xr_variable_processor=xr_variable_processor,
            # frequency_metadata_keys=frequency_metadata_keys,
            # bounds_info=bounds_info,
            allow_cf_checker_warnings=False,
        )
        validation_res.raise_if_errors()


if __name__ == "__main__":
    main()
