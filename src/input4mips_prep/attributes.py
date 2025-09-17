"""
Attribute manipulation
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import netCDF4


def copy_attributes(
    source: netCDF4.Dataset, target: netCDF4.Dataset
) -> netCDF4.Dataset:
    """
    Copy attributes from source to target

    Parameters
    ----------
    source
        Source dataset from which to copy attributes

    target
        Target dataset to which to copy attributes

    Returns
    -------
    :
        `target` with updated attributes

        Note that the operation occurs in place,
        so the input `target` object is also affected
        (returning `target` is done for convenience)
    """
    for name in source.ncattrs():
        target.setncattr(name, source.getncattr(name))
