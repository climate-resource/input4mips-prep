"""
Time co-ordinate handling
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import netCDF4


def rewrite_time_from_bounded_to_climatology(
    ds: netCDF4.Dataset,
    original_time_bounds: netCDF4.Variable,
    frequency: str | None = None,
) -> netCDF4.Dataset:
    """
    Re-write the time dimension from being bounded to being a climatology

    Parameters
    ----------
    ds
        Dataset on which to perform the translation

    original_time_bounds
        Original time bounds

    frequency
        Frequency to set

        If not supplied, we use the existing frequency plus "C"

    Returns
    -------
    :
        `ds` with updated attributes

        Note that the operation occurs in place,
        so the input `ds` object is also affected
        (returning `ds` is done for convenience)
    """
    if frequency is None:
        frequency = f"{ds.getncattr('frequency')}C"

    ds.setncattr("frequency", frequency)

    ds["time"].setncattr("climatology", "climatology_bounds")
    ds["time"].delncattr("bounds")

    ds.createDimension("nv", 2)
    ds.createVariable("climatology_bounds", ds["time"].datatype, ("time", "nv"))
    ds["climatology_bounds"][:] = original_time_bounds[:].T

    return ds
