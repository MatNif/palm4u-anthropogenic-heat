"""Module with objects that could be read or written to netCDF"""
import os
from dataclasses import dataclass
from typing import ClassVar, List, Optional, Tuple

from netCDF4 import Dataset  # type: ignore
from numpy import ma
import numpy.typing as npt


@dataclass
class NCDFDimension:
    """A dimension that could written to netCDF. A corresponding dimension variable is included.
    Its values can be stored in the attribute values or supplied when written.
    """

    name: str
    datatype: str

    long_name: Optional[str] = None
    units: Optional[str] = None
    standard_name: Optional[str] = None
    values: Optional[npt.NDArray] = None

    _metadata: ClassVar[List[str]] = ["long_name", "units", "standard_name"]

    @property
    def size(self) -> int:
        """Get the size of the dimension. It is derived from the values attribute."""
        if self.values is None:
            raise ValueError("Values of dimension " + self.name + " not defined.")
        return self.values.size

    def __len__(self) -> int:
        return self.size

    def to_dataset(self, nc_data: Dataset, values: Optional[ma.MaskedArray] = None) -> None:
        """Add dimension and its dimension variable to a netCDF Dataset if it does not already
        include it. The values are either supplied in the function call or taken from the values
        attribute.
        """
        if self.name not in nc_data.dimensions:
            if values is not None:
                to_write = values
            elif self.values is not None:
                to_write = self.values
            else:
                raise ValueError("Values of dimension " + self.name + " not defined.")

            print("Writing dimension " + self.name + " to file...")

            nc_data.createDimension(self.name, len(to_write))

            nc_var = nc_data.createVariable(self.name, self.datatype, self.name)
            nc_var[:] = to_write

            for attr in self._metadata:
                attr_value = getattr(self, attr)
                if attr_value is not None:
                    nc_var.setncattr(attr, attr_value)


@dataclass
class NCDFVariable:
    """A variable that could written to netCDF. It includes its metadata and dimensions.
    Its values can be stored in the attribute values or supplied when written. A default filename
    can be also stored.
    """

    name: str
    dimensions: Tuple[NCDFDimension, ...]
    datatype: str
    fillvalue: float
    long_name: str
    units: str

    values: Optional[npt.NDArray] = None

    coordinates: Optional[str] = None
    grid_mapping: Optional[str] = None
    lod: Optional[int] = None
    res_orig: Optional[float] = None
    standard_name: Optional[str] = None

    filename: Optional[str] = None

    _metadata: ClassVar[List[str]] = [
        "long_name",
        "units",
        "standard_name",
        "res_orig",
        "lod",
        "coordinates",
        "grid_mapping",
    ]

    def to_file(
        self, values: Optional[npt.ArrayLike] = None, filename: Optional[str] = None
    ) -> None:
        """Write variable to a netCDF file, which is openend and closed. The file is either
        specified in the function call or taken from the default filename. If the variable was not
        yet added the file, it is added with its dimensions, otherwise its values are overwritten.
        The values are either supplied in the function call or taken from the values attribute.
        """

        if values is not None:
            to_write = values
        elif self.values is not None:
            to_write = self.values
        else:
            raise ValueError("Values of variable " + self.name + " not defined.")

        if filename is not None:
            to_file = filename
        elif self.filename is not None:
            to_file = self.filename
        else:
            raise ValueError("Output filename for variable " + self.name + " not defined.")

        try:
            nc_data = Dataset(to_file, "a", format="NETCDF4")
        except FileNotFoundError:
            print("Error. Could not open file: " + to_file + ". Aborting...")
            raise

        print("Writing array " + self.name + " to file...")
        if self.name not in nc_data.variables:
            for nc_dim in self.dimensions:
                nc_dim.to_dataset(nc_data)

            nc_var = nc_data.createVariable(
                self.name,
                self.datatype,
                (o.name for o in self.dimensions),
                fill_value=self.fillvalue,
            )

            for attr in self._metadata:
                attr_value = getattr(self, attr)
                if attr_value is not None:
                    nc_var.setncattr(attr, attr_value)

        else:
            nc_var = nc_data.variables[self.name]

        if len(self.dimensions) == 1:
            nc_var[:] = to_write
        elif len(self.dimensions) == 2:
            nc_var[:, :] = to_write
        elif len(self.dimensions) == 3:
            nc_var[:, :, :] = to_write
        else:
            raise NotImplementedError

        nc_data.close()

    def from_file(self, filename: Optional[str] = None) -> ma.MaskedArray:
        """Return the variable from a netCDF file. The file is either specified in the
        function call or taken from the default filename.
        """

        if filename is not None:
            from_file = filename
        elif self.filename is not None:
            from_file = self.filename
        else:
            raise ValueError("Input filename for variable " + self.name + " not defined.")

        try:
            nc_data = Dataset(from_file, "r", format="NETCDF4")
        except FileNotFoundError:
            print("Error. Could not open file: " + from_file + ". Aborting...")
            raise

        tmp_array = nc_data.variables[self.name][:]
        nc_data.close()

        return tmp_array


@dataclass
class NCDFCoordinateReferenceSystem:
    """A coordinate reference system that can be written to a netCDF file."""

    long_name: str
    grid_mapping_name: str
    semi_major_axis: float
    inverse_flattening: float
    longitude_of_prime_meridian: float
    longitude_of_central_meridian: float
    scale_factor_at_central_meridian: float
    latitude_of_projection_origin: float
    false_easting: float
    false_northing: float
    spatial_ref: str
    units: str
    epsg_code: str

    filename: Optional[str] = None

    def to_file(self, filename: Optional[str] = None) -> None:
        """Write CRS to a netCDF file, which is openend and closed. The file is either specified
        in the function call or taken from the default filename.
        """

        if filename is not None:
            to_file = filename
        elif self.filename is not None:
            to_file = self.filename
        else:
            raise ValueError("Output filename for CRS not defined.")

        try:
            nc_data = Dataset(to_file, "a", format="NETCDF4")
        except FileNotFoundError:
            print("Error. Could not open file: " + to_file + ". Aborting...")
            raise

        print("Writing crs to file...")

        nc_var = nc_data.createVariable("crs", "i")

        nc_var.long_name = self.long_name
        nc_var.grid_mapping_name = self.grid_mapping_name
        nc_var.semi_major_axis = self.semi_major_axis
        nc_var.inverse_flattening = self.inverse_flattening
        nc_var.longitude_of_prime_meridian = self.longitude_of_prime_meridian
        nc_var.longitude_of_central_meridian = self.longitude_of_central_meridian
        nc_var.scale_factor_at_central_meridian = self.scale_factor_at_central_meridian
        nc_var.latitude_of_projection_origin = self.latitude_of_projection_origin
        nc_var.false_easting = self.false_easting
        nc_var.false_northing = self.false_northing
        nc_var.spatial_ref = self.spatial_ref
        nc_var.units = self.units
        nc_var.epsg_code = self.epsg_code

        nc_data.close()


def remove_existing_file(filename) -> None:
    """Remove a file if it exists."""
    try:
        os.remove(filename)
    except FileNotFoundError:
        pass
