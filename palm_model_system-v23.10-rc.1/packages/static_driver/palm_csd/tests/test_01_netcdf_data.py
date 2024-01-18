"""Run netcdf_data tests"""
import os

import pytest
from netCDF4 import Dataset  # type: ignore
from numpy import ma

from palm_csd.netcdf_data import NCDFDimension, NCDFVariable, remove_existing_file
from tests.tools import ncdf_equal


def test_01_variables():
    """Test dimensions and variables"""
    # with defined values
    x_dim = NCDFDimension(
        name="x",
        values=ma.MaskedArray([1.0, 3.0, 5.0]),
        datatype="f4",
        standard_name="projection_x_coordinate",
        long_name="x",
        units="m",
    )
    assert x_dim.size == 3
    assert len(x_dim) == 3

    # without defined values
    y_dim = NCDFDimension(
        "y",
        datatype="f4",
        standard_name="projection_y_coordinate",
        long_name="y",
        units="m",
    )
    with pytest.raises(ValueError):
        print(y_dim.size)
    with pytest.raises(ValueError):
        len(y_dim)

    # write to file
    to_file_dimension = "tests/01_netcdf_data/dimension.nc"
    remove_existing_file(to_file_dimension)
    nc_data = Dataset(to_file_dimension, "a", format="NETCDF4")
    x_dim.to_dataset(nc_data)
    # fail because no values
    with pytest.raises(ValueError):
        y_dim.to_dataset(nc_data)

    # assign values
    y_dim.values = ma.MaskedArray([2.0, 4.0, 6.0, 8.0])
    assert y_dim.size == 4
    assert len(y_dim) == 4
    y_dim.to_dataset(nc_data)

    nc_data.close()

    # create variable
    buildings_var = NCDFVariable(
        name="buildings_2d",
        dimensions=(y_dim, x_dim),
        datatype="f4",
        fillvalue=-9999.0,
        long_name="buildings",
        units="m",
        res_orig=2,
        lod=1,
        coordinates="E_UTM N_UTM lon lat",
        grid_mapping="crs",
    )
    to_file_variable = "tests/01_netcdf_data/variable.nc"
    remove_existing_file(to_file_variable)

    # no defined values
    with pytest.raises(ValueError):
        buildings_var.to_file(nc_data)

    buildings_var.values = ma.arange(1, 13, 1).reshape([4, 3])
    # no defined filename
    with pytest.raises(ValueError):
        buildings_var.to_file(nc_data)

    buildings_var.to_file(filename=to_file_variable)

    buildings_id = NCDFVariable(
        name="buildings_id",
        dimensions=(y_dim, x_dim),
        values=ma.MaskedArray(range(10, 22)).reshape([4, 3]),
        datatype="i",
        fillvalue=9999,
        long_name="buildings id",
        units="",
        coordinates="E_UTM N_UTM lon lat",
        grid_mapping="crs",
        filename=to_file_variable,
    )

    buildings_id.to_file()

    # compare output files
    assert ncdf_equal(to_file_dimension, "tests/01_netcdf_data/output/dimension.nc")
    assert ncdf_equal(to_file_variable, "tests/01_netcdf_data/output/variable.nc")

    os.remove(to_file_dimension)
    os.remove(to_file_variable)
