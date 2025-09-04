"""
Variable manipulation
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import netCDF4


def copy_variables(source: netCDF4.Dataset, target: netCDF4.Dataset) -> netCDF4.Dataset:
    """
    Copy variables from source to target

    Parameters
    ----------
    source
        Source dataset from which to copy variables

    target
        Target dataset to which to copy variables

    Returns
    -------
    :
        `target` with updated variables

        Note that the operation occurs in place,
        so the input `target` object is also affected
        (returning `target` is done for convenience)
    """
    for name, variable in source.variables.items():
        dimensions = variable.dimensions

        create_kwargs = {}
        variable_attrs_to_copy = variable.ncattrs()
        if "_FillValue" in variable_attrs_to_copy:
            fill_value = variable.getncattr("_FillValue")
            create_kwargs["fill_value"] = fill_value
            variable_attrs_to_copy.pop(variable_attrs_to_copy.index("_FillValue"))

        target.createVariable(name, variable.datatype, dimensions)
        # Copy the data
        # unclear if this only works with 1D data
        target.variables[name][:] = variable[:]

        for attr in variable_attrs_to_copy:
            target[name].setncattr(attr, variable.getncattr(attr))
