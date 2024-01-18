from typing import Optional, Tuple

from netCDF4 import Dataset, Variable  # type: ignore
from numpy import ma

from palm_csd.csd_config import CSDConfig, CSDConfigAttributes, CSDConfigDomain, CSDConfigInput
from palm_csd.geo_converter import GeoConverter
from palm_csd.netcdf_data import (
    NCDFCoordinateReferenceSystem,
    NCDFDimension,
    NCDFVariable,
    remove_existing_file,
)


class CSDDomain:
    """A domain that stores all its configurations, output dimensions and variables."""

    name: str

    config: CSDConfigDomain
    input_config: CSDConfigInput
    attributes: CSDConfigAttributes

    # TODO: Python 3.11: Use Self to indicate same type as class
    parent: Optional["CSDDomain"]

    geo_converter: Optional[GeoConverter]

    rotation_angle: float
    origin_x: Optional[float]
    origin_y: Optional[float]
    origin_lon: Optional[float]
    origin_lat: Optional[float]
    origin_z: Optional[float]

    x: NCDFDimension
    y: NCDFDimension
    z: NCDFDimension
    zlad: NCDFDimension

    nsurface_fraction: NCDFDimension
    nbuilding_pars: NCDFDimension
    nvegetation_pars: NCDFDimension
    nwater_pars: NCDFDimension

    lat: NCDFVariable
    lon: NCDFVariable

    x_global: NCDFVariable
    y_global: NCDFVariable

    E_UTM: NCDFVariable
    N_UTM: NCDFVariable

    zt: NCDFVariable

    buildings_2d: NCDFVariable
    building_id: NCDFVariable
    building_type: NCDFVariable
    buildings_3d: NCDFVariable

    surface_fraction: NCDFVariable

    vegetation_type: NCDFVariable
    pavement_type: NCDFVariable
    water_type: NCDFVariable
    soil_type: NCDFVariable
    street_type: NCDFVariable
    street_crossing: NCDFVariable

    building_pars: NCDFVariable
    vegetation_pars: NCDFVariable
    water_pars: NCDFVariable

    lad: NCDFVariable
    bad: NCDFVariable
    tree_id: NCDFVariable
    tree_type: NCDFVariable

    def __init__(
        self,
        name: str,
        config: CSDConfig,
        parent: Optional["CSDDomain"] = None,
    ) -> None:
        self.name = name

        # configurations
        self.config = config.domain_dict[name]
        self.input_config = config.input_of_domain(self.config)

        self.attributes = config.attributes

        self.parent = parent

        # converter for geo data
        if config.settings.epsg is None:
            self.geo_converter = None
        else:
            # Find root parent
            if self.parent is not None:
                parent_geoconverter = self.parent.geo_converter
                if parent_geoconverter is None:
                    raise ValueError("Parent domain has no geo converter")
                parent_tmp = self.parent
                while parent_tmp.parent is not None:
                    parent_tmp = parent_tmp.parent
                root_parent_geoconverter = parent_tmp.geo_converter
                if root_parent_geoconverter is None:
                    raise ValueError("Root parent domain has no geo converter")
            else:
                parent_geoconverter = None
                root_parent_geoconverter = None
            print(f"Setting up of coordinate calculation for domain {self.name}...")
            self.geo_converter = GeoConverter(
                self.config, config.settings, parent_geoconverter, root_parent_geoconverter
            )

        self.rotation_angle = config.settings.rotation_angle

        self.check_consistency()

        # dimensions
        self.x = NCDFDimension(
            name="x",
            datatype="f4",
            standard_name="projection_x_coordinate",
            long_name="x",
            units="m",
        )
        self.y = NCDFDimension(
            name="y",
            datatype="f4",
            standard_name="projection_y_coordinate",
            long_name="y",
            units="m",
        )
        self.z = NCDFDimension(name="z", datatype="f4", long_name="z", units="m")

        self.nsurface_fraction = NCDFDimension(name="nsurface_fraction", datatype="i")
        self.nbuilding_pars = NCDFDimension(name="nbuilding_pars", datatype="i")
        self.nvegetation_pars = NCDFDimension(name="nvegetation_pars", datatype="i")
        self.nwater_pars = NCDFDimension(name="nwater_pars", datatype="i")

        self.zlad = NCDFDimension(name="zlad", datatype="f4")

        dimensions_yx = (self.y, self.x)
        dimensions_zladyx = (self.zlad, self.y, self.x)

        # variables
        self.lat = self._variable_float(
            name="lat",
            dimensions=dimensions_yx,
            long_name="latitude",
            standard_name="latitude",
            units="degrees_north",
            full=False,
        )
        self.lon = self._variable_float(
            name="lon",
            dimensions=dimensions_yx,
            long_name="longitude",
            standard_name="longitude",
            units="degrees_east",
            full=False,
        )

        self.x_global = self._variable_float(
            name="x_UTM",
            dimensions=(self.x,),
            long_name="easting",
            standard_name="projection_x_coordinate",
            units="m",
            full=False,
        )

        self.y_global = self._variable_float(
            name="y_UTM",
            dimensions=(self.y,),
            long_name="northing",
            standard_name="projection_y_coordinate",
            units="m",
            full=False,
        )

        self.E_UTM = self._variable_float(
            name="E_UTM",
            dimensions=dimensions_yx,
            long_name="easting",
            standard_name="projection_x_coordinate",
            units="m",
            full=False,
        )

        self.N_UTM = self._variable_float(
            name="N_UTM",
            dimensions=dimensions_yx,
            long_name="northing",
            standard_name="projection_y_coordinate",
            units="m",
            full=False,
        )

        self.zt = self._variable_float(
            name="zt",
            dimensions=dimensions_yx,
            long_name="orography",
            units="m",
        )

        self.buildings_2d = self._variable_float(
            name="buildings_2d",
            dimensions=dimensions_yx,
            long_name="buildings",
            units="m",
            lod=1,
        )

        self.building_id = self._variable_int(
            name="building_id",
            dimensions=dimensions_yx,
            long_name="building id",
            units="",
        )

        self.building_type = self._variable_byte(
            name="building_type",
            dimensions=dimensions_yx,
            long_name="building type",
            units="",
        )

        self.buildings_3d = self._variable_byte(
            name="buildings_3d",
            dimensions=(self.z, self.y, self.x),
            long_name="buildings 3d",
            units="",
            lod=2,
        )

        self.surface_fraction = self._variable_float(
            name="surface_fraction",
            dimensions=(self.nsurface_fraction, self.y, self.x),
            long_name="surface fraction",
            units="1",
        )

        self.vegetation_type = self._variable_byte(
            name="vegetation_type",
            dimensions=dimensions_yx,
            long_name="vegetation type",
            units="",
        )

        self.pavement_type = self._variable_byte(
            name="pavement_type",
            dimensions=dimensions_yx,
            long_name="pavement type",
            units="",
        )

        self.water_type = self._variable_byte(
            name="water_type",
            dimensions=dimensions_yx,
            long_name="water type",
            units="",
        )

        self.soil_type = self._variable_byte(
            name="soil_type",
            dimensions=dimensions_yx,
            long_name="soil type",
            units="",
        )

        self.street_type = self._variable_byte(
            name="street_type",
            dimensions=dimensions_yx,
            long_name="street type",
            units="",
        )

        self.street_crossing = self._variable_byte(
            name="street_crossing",
            dimensions=dimensions_yx,
            long_name="street crossings",
            units="",
        )

        self.building_pars = self._variable_float(
            name="building_pars",
            dimensions=(self.nbuilding_pars, self.y, self.x),
            long_name="building_pars",
            units="",
        )

        self.vegetation_pars = self._variable_float(
            name="vegetation_pars",
            dimensions=(self.nvegetation_pars, self.y, self.x),
            long_name="vegetation_pars",
            units="",
        )

        self.water_pars = self._variable_float(
            name="water_pars",
            dimensions=(self.nwater_pars, self.y, self.x),
            long_name="water_pars",
            units="",
        )

        self.lad = self._variable_float(
            name="lad",
            dimensions=dimensions_zladyx,
            long_name="leaf area density",
            units="m2 m-3",
        )

        self.bad = self._variable_float(
            name="bad",
            dimensions=dimensions_zladyx,
            long_name="basal area density",
            units="m2 m-3",
        )

        self.tree_id = self._variable_int(
            name="tree_id",
            dimensions=dimensions_zladyx,
            long_name="tree id",
            units="",
        )

        self.tree_type = self._variable_byte(
            name="tree_type",
            dimensions=dimensions_zladyx,
            long_name="tree type",
            units="",
        )

    def check_consistency(self) -> None:
        # Geographical coordinate input

        if (
            self.input_config.file_x_UTM is None
            or self.input_config.file_y_UTM is None
            or self.input_config.file_lon is None
            or self.input_config.file_lat is None
        ):
            if self.geo_converter is None:
                raise ValueError(
                    f"{self.name}: Not all coordinate inputs provided "
                    "but target reference system not defined."
                )

    def _variable_float(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        full: bool = True,
    ) -> NCDFVariable:
        """Helper functions that returns a variables with some predefined atrributes
        for float values.
        """

        default_values = {
            "datatype": "f4",
            "fillvalue": -9999.0,
            "filename": self.config.filename,
        }
        if full:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def _variable_int(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        full: bool = True,
    ) -> NCDFVariable:
        """Helper functions that returns a variables with some predefined atrributes
        for int values.
        """

        default_values = {
            "datatype": "i",
            "fillvalue": -9999,
            "filename": self.config.filename,
        }
        if full:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def _variable_byte(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        full: bool = True,
    ) -> NCDFVariable:
        """Helper functions that returns a variables with some predefined atrributes
        for byte values.
        """

        default_values = {
            "datatype": "b",
            "fillvalue": -127,
            "filename": self.config.filename,
        }
        if full:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def remove_existing_output(self) -> None:
        """Remove configured output file if it exists."""
        remove_existing_file(self.config.filename)

    def write_global_attributes(self) -> None:
        """Write global attributes to the netcdf filename, None attributes are not added."""

        print("Writing global attributes to file...")

        nc_data = Dataset(self.config.filename, "a", format="NETCDF4")

        nc_data.setncattr("Conventions", "CF-1.7")

        all_attributes = vars(self.attributes)
        for attribute in all_attributes:
            if all_attributes[attribute] is not None:
                nc_data.setncattr(attribute, all_attributes[attribute])

        # add additional attributes
        for attribute in [
            "rotation_angle",
            "origin_x",
            "origin_y",
            "origin_lon",
            "origin_lat",
            "origin_z",
        ]:
            if getattr(self, attribute) is not None:
                nc_data.setncattr(attribute, getattr(self, attribute))
            else:
                raise Exception(f"Attribute {attribute} not set.")

        nc_data.close()

    def write_crs_to_file(self) -> None:
        """Write crs in CF convention to the netcdf filename. Values are taken from
        geo_converter's dst_crs.
        """

        print("Writing crs to file...")

        if self.geo_converter is None:
            raise ValueError("geoconverter must not be None.")
        crs_dict = self.geo_converter.dst_crs_to_cf()

        try:
            nc_data = Dataset(self.config.filename, "a", format="NETCDF4")
        except FileNotFoundError:
            print(f"Error. Could not open file: {self.config.filename}. Aborting...")
            raise

        nc_var = nc_data.createVariable("crs", "i")

        # Add long_name
        nc_var.setncattr("long_name", "coordinate reference system")

        # Add crs information
        for key, value in crs_dict.items():
            nc_var.setncattr(key, value)

        nc_data.close()

    def read_from_file_3d(
        self,
        filename: Optional[str],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
        y0: Optional[int] = None,
        y1: Optional[int] = None,
        z0: Optional[int] = None,
        z1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 3d variable from a netCDF file, which is openend and closed. If the filename
        is None, the values of the returned array are all masked. The default boundary
        coordinates are taken from the containing domain. If complete, the full variable is read.
        """

        if x0 is None:
            x0 = self.config.x0
        if x1 is None:
            x1 = self.config.x1
        if y0 is None:
            y0 = self.config.y0
        if y1 is None:
            y1 = self.config.y1
        if z0 is None:
            if not complete:
                raise NotImplementedError
        if z1 is None:
            if not complete:
                raise NotImplementedError

        if filename is not None:
            try:
                nc_data = Dataset(filename, "r", format="NETCDF4")
            except FileNotFoundError:
                print("Error: " + filename + ". No such file. Aborting...")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 3)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:, :, :]
            else:
                tmp_array = variable[z0 : (z1 + 1), y0 : (y1 + 1), x0 : (x1 + 1)]  # type: ignore
            nc_data.close()
        else:
            if complete:
                raise ValueError("filename needs to given when complete==True")
            tmp_array = ma.masked_all((z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1))  # type: ignore

        return tmp_array

    def read_from_file_2d(
        self,
        filename: Optional[str],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
        y0: Optional[int] = None,
        y1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 2d variable from a netCDF file, which is openend and closed. If the filename
        is None, the values of the returned array are all masked. The default boundary
        coordinates are taken from the containing domain. If complete, the full variable is read.
        """

        if x0 is None:
            x0 = self.config.x0
        if x1 is None:
            x1 = self.config.x1
        if y0 is None:
            y0 = self.config.y0
        if y1 is None:
            y1 = self.config.y1

        if filename is not None:
            try:
                nc_data = Dataset(filename, "r", format="NETCDF4")
            except FileNotFoundError:
                print("Error: " + filename + ". No such file. Aborting...")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 2)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:, :]
            else:
                tmp_array = variable[y0 : (y1 + 1), x0 : (x1 + 1)]
            nc_data.close()
        else:
            if complete:
                raise ValueError("filename needs to given when complete==True")
            tmp_array = ma.masked_all((y1 - y0 + 1, x1 - x0 + 1))

        return tmp_array

    def read_from_file_1d(
        self,
        filename: Optional[str],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 1d variable from a netCDF file, which is openend and closed. If the filename
        is None, the values of the returned array are all masked. The default boundary
        coordinates are taken from the containing domain. If complete, the full variable is read.
        """

        if x0 is None:
            x0 = self.config.x0
        if x1 is None:
            x1 = self.config.x1

        if filename is not None:
            try:
                nc_data = Dataset(filename, "r", format="NETCDF4")
            except FileNotFoundError:
                print("Error: " + filename + ". No such file. Aborting...")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 2)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:]
            else:
                tmp_array = variable[x0 : (x1 + 1)]
            nc_data.close()
        else:
            if complete:
                raise ValueError("filename needs to given when complete==True")
            tmp_array = ma.masked_all(x1 - x0 + 1)

        return tmp_array

    def read_from_file_crs(
        self, filename: Optional[str] = None, varname: Optional[str] = None
    ) -> NCDFCoordinateReferenceSystem:
        """Return coordinate reference system from a netCDF file, which is opened and closed."""

        if filename is not None:
            from_file = filename
        elif self.input_config.file_x_UTM is not None:
            from_file = self.input_config.file_x_UTM
        else:
            raise ValueError("filename or input_config.file_x_UTM needs to be not None")

        try:
            nc_data = Dataset(from_file, "r", format="NETCDF4")
        except FileNotFoundError:
            print("Error: " + from_file + ". No such file. Aborting...")
            raise

        if varname is None:
            variable = _find_variable_name(nc_data, 2)
        else:
            variable = nc_data.variables[varname]
        crs_from_file = nc_data.variables[variable.grid_mapping]

        # Get EPSG code from crs
        try:
            epsg_code = crs_from_file.epsg_code
        except AttributeError:
            epsg_code = "unknown"
            if crs_from_file.spatial_ref.find("ETRS89", 0, 100) and crs_from_file.spatial_ref.find(
                "UTM", 0, 100
            ):
                if crs_from_file.spatial_ref.find("28N", 0, 100) != -1:
                    epsg_code = "EPSG:25828"
                elif crs_from_file.spatial_ref.find("29N", 0, 100) != -1:
                    epsg_code = "EPSG:25829"
                elif crs_from_file.spatial_ref.find("30N", 0, 100) != -1:
                    epsg_code = "EPSG:25830"
                elif crs_from_file.spatial_ref.find("31N", 0, 100) != -1:
                    epsg_code = "EPSG:25831"
                elif crs_from_file.spatial_ref.find("32N", 0, 100) != -1:
                    epsg_code = "EPSG:25832"
                elif crs_from_file.spatial_ref.find("33N", 0, 100) != -1:
                    epsg_code = "EPSG:25833"
                elif crs_from_file.spatial_ref.find("34N", 0, 100) != -1:
                    epsg_code = "EPSG:25834"
                elif crs_from_file.spatial_ref.find("35N", 0, 100) != -1:
                    epsg_code = "EPSG:25835"
                elif crs_from_file.spatial_ref.find("36N", 0, 100) != -1:
                    epsg_code = "EPSG:25836"
                elif crs_from_file.spatial_ref.find("37N", 0, 100) != -1:
                    epsg_code = "EPSG:25837"

        crs_var = NCDFCoordinateReferenceSystem(
            long_name="coordinate reference system",
            grid_mapping_name=crs_from_file.grid_mapping_name,
            semi_major_axis=crs_from_file.semi_major_axis,
            inverse_flattening=crs_from_file.inverse_flattening,
            longitude_of_prime_meridian=crs_from_file.longitude_of_prime_meridian,
            longitude_of_central_meridian=crs_from_file.longitude_of_central_meridian,
            scale_factor_at_central_meridian=crs_from_file.scale_factor_at_central_meridian,
            latitude_of_projection_origin=crs_from_file.latitude_of_projection_origin,
            false_easting=crs_from_file.false_easting,
            false_northing=crs_from_file.false_northing,
            spatial_ref=crs_from_file.spatial_ref,
            units="m",
            epsg_code=epsg_code,
            filename=self.config.filename,
        )

        nc_data.close()

        return crs_var


def _find_variable_name(nc_data: Dataset, ndim: int) -> Variable:
    """Find the Variable of the input Dataset with the given number of dimensions.
    Exclude dimension variables.
    """

    dimension_names = list(nc_data.dimensions.keys())

    variables_correct_dim = []
    for name, variable in nc_data.variables.items():
        if len(variable.dimensions) == ndim and name not in dimension_names:
            variables_correct_dim.append(name)

    nfound = len(variables_correct_dim)
    if nfound == 0:
        raise ValueError(f"No suitable variable with {ndim} dimension found.")
    elif nfound > 1:
        raise ValueError(f"Found {nfound} suitable variables when expecting 1.")
    return nc_data.variables[variables_correct_dim[0]]
