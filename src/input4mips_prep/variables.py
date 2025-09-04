"""
Variable manipulation
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import netCDF4


def copy_variable(
    source: netCDF4.Dataset, target: netCDF4.Dataset, variable: str
) -> netCDF4.Dataset:
    """
    Copy a variable from source to target

    Parameters
    ----------
    source
        Source dataset from which to copy the variable

    target
        Target dataset to which to copy the variable

    variable
        Variable to copy

    Returns
    -------
    :
        `target` with variable now included

        Note that the operation occurs in place,
        so the input `target` object is also affected
        (returning `target` is done for convenience)
    """
    variable_nc = source.variables[variable]

    dimensions = variable_nc.dimensions

    create_kwargs = {}
    variable_attrs_to_copy = variable_nc.ncattrs()
    if "_FillValue" in variable_attrs_to_copy:
        fill_value = variable_nc.getncattr("_FillValue")
        create_kwargs["fill_value"] = fill_value
        variable_attrs_to_copy.pop(variable_attrs_to_copy.index("_FillValue"))

    target.createVariable(variable, variable_nc.datatype, dimensions)
    # Copy the data
    # unclear if this only works with 1D data
    target.variables[variable][:] = variable_nc[:]

    for attr in variable_attrs_to_copy:
        target[variable].setncattr(attr, variable_nc.getncattr(attr))

    return target


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
        target = copy_variable(source, target, variable)

    return target


def reverse_variable_dimensions(
    source: netCDF4.Dataset, target: netCDF4.Dataset, variable: str
) -> netCDF4.Dataset:
    """
    Copy a variable from source to target, reversing the dimensions in the process

    Parameters
    ----------
    source
        Source dataset from which to copy the variable

    target
        Target dataset to which to write the transposed variable

    variable
        Variable to copy

    Returns
    -------
    :
        `target` with variable now included

        Note that the operation occurs in place,
        so the input `target` object is also affected
        (returning `target` is done for convenience)
    """
    variable_nc = source.variables[variable]

    dimensions = variable_nc.dimensions

    create_kwargs = {}
    variable_attrs_to_copy = variable_nc.ncattrs()
    if "_FillValue" in variable_attrs_to_copy:
        fill_value = variable_nc.getncattr("_FillValue")
        create_kwargs["fill_value"] = fill_value
        variable_attrs_to_copy.pop(variable_attrs_to_copy.index("_FillValue"))

    dimensions_new = dimensions[::-1]
    target.createVariable(variable, variable_nc.datatype, dimensions_new)
    # Tranpose the data
    target.variables[variable][:] = variable_nc[:].T

    for attr in variable_attrs_to_copy:
        target[variable].setncattr(attr, variable_nc.getncattr(attr))

    return target
