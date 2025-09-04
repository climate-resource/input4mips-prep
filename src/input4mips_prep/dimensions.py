"""
Dimension manipulation
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import netCDF4


def copy_dimensions(
    source: netCDF4.Dataset, target: netCDF4.Dataset
) -> netCDF4.Dataset:
    """
    Copy dimensions from source to target

    Parameters
    ----------
    source
        Source dataset from which to copy dimensions

    target
        Target dataset to which to copy dimensions

    Returns
    -------
    :
        `target` with updated dimensions

        Note that the operation occurs in place,
        so the input `target` object is also affected
        (returning `target` is done for convenience)
    """
    for name, dimension in source.dimensions.items():
        target.createDimension(
            name, (len(dimension) if not dimension.isunlimited else None)
        )
