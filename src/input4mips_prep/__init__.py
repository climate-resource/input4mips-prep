"""
Supporting functions for scripts used to help preparing input4MIPs data
"""

from __future__ import annotations

from input4mips_prep.attributes import copy_attributes
from input4mips_prep.dimensions import copy_dimension, copy_dimensions
from input4mips_prep.time import rewrite_time_from_bounded_to_climatology
from input4mips_prep.variables import (
    copy_variable,
    copy_variable_iiasa,
    copy_variables,
    reverse_variable_dimensions,
)

__all__ = [
    "copy_attributes",
    "copy_dimension",
    "copy_dimensions",
    "copy_variable",
    "copy_variable_iiasa",
    "copy_variables",
    "reverse_variable_dimensions",
    "rewrite_time_from_bounded_to_climatology",
]
