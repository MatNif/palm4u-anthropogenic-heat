"""Tools for testing"""
from pprint import pprint
from typing import List, Optional

from deepdiff import DeepDiff
from netCDF4 import Dataset  # type: ignore
from numpy import ma


def ncdf_equal(
    ncdf_ref: str,
    ncdf_com: str,
    check_metadata: bool = True,
    check_fields: bool = True,
    metadata_significant_digits: Optional[int] = None,
    metadata_exclude_regex_paths: Optional[List[str]] = None,
) -> bool:
    """Test if two netCDF files are equal in terms of values and metadata."""
    nc_data_ref = Dataset(ncdf_ref, "r")
    nc_data_com = Dataset(ncdf_com, "r")

    diff = {}
    if check_metadata:
        # metadata diff, a Dataset does not yet load the values of the fields
        # exclude "_grpid" and "_varid" because they are expected to differ and
        #  do not have an impact on the content of a netCDF file
        exclude = ["_grpid", "_varid"]
        if metadata_exclude_regex_paths is not None:
            exclude.extend(metadata_exclude_regex_paths)
        metadata_diff = DeepDiff(
            nc_data_ref,
            nc_data_com,
            exclude_regex_paths=exclude,
            significant_digits=metadata_significant_digits,
            number_format_notation="e",
        )
        diff.update(metadata_diff.to_dict())

    if check_fields:
        # value diff
        value_diff = {}
        keys_ref = nc_data_ref.variables.keys()
        keys_com = nc_data_com.variables.keys()
        for variable in keys_ref & keys_com:
            field_ref = nc_data_ref.variables[variable][:]
            field_com = nc_data_com.variables[variable][:]
            if not ma.allclose(field_ref, field_com):
                value_diff.update({variable: "field values differ"})
        diff.update(value_diff)

    nc_data_ref.close()
    nc_data_com.close()

    if diff:
        pprint(diff)
        return False

    return True
