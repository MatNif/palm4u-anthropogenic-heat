#!/usr/bin/env python3
# -------------------------------------------------------------------- #
# LICENSE INFORMATION
# 
# @TODO x/y coordinates of point souces in metre or grid coordinates?
# -------------------------------------------------------------------- #

import datetime
import math

from netCDF4 import Dataset
import numpy as np


class AHDriver:
    """This is an example script to generate anthropogenic heat emissions 
    drivers for PALM.

    You can use it as a starting point for creating your setup specific
    driver.
    """

    def __init__(self):
        """Open the driver as netCDF4 file. Here, you have to
        give the full path to the static driver that shall be created.
        Existing file with same name is deleted.
        """
        print('Opening file...')
        self.nc_file = Dataset('example_ah_file.nc', 'w', format='NETCDF4')

    def write_global_attributes(self):
        """Write global attributes to driver."""
        print("Writing global attributes...")

        # Optional global attributes
        # --------------------------
        self.nc_file.title = 'Example PALM ah driver'
        self.nc_file.author = 'Tobias Gronemeier, Mathias Niffeler'
        self.nc_file.institution = 'iMA Richter Roeckle GmbH Co KG, ' \
            'Singapore ETH Centre'
        self.nc_file.comment = 'Testing setup'
        self.nc_file.creation_date = datetime.datetime.utcnow().strftime('%y-%m-%d %H:%M:%S %z')
        self.nc_file.history = ''
        self.nc_file.keywords = 'example, PALM-4U'
        self.nc_file.license = ''
        self.nc_file.palm_version = ''
        self.nc_file.references = ''
        self.nc_file.source = ''
        self.nc_file.version = '1'

        # Mandatory global attributes
        # ---------------------------
        self.nc_file.Conventions = 'CF-1.7'
        self.nc_file.origin_lat = 1.303687    # (overwrite initialization_parameters)
        self.nc_file.origin_lon = 103.773933  # Used to initialize Coriolis parameter
        self.nc_file.origin_time = '2024-01-01 10:00:00 +08'
        self.nc_file.origin_x = 363594.0
        self.nc_file.origin_y = 144130.0
        self.nc_file.origin_z = 0.0
        self.nc_file.rotation_angle = 0.0

    def define_dimensions(self):
        """Set dimensions on which variables are defined."""
        print("Writing dimensions...")

        # Specify size
        # ------------
        self.nbuildings = 2
        self.nstreets = 3
        self.npoints = 5
        self.nt = 24

        # Define coordinates
        # ------------------
        self.nc_file.createDimension('building_id', self.nbuildings)
        self.building_id = self.nc_file.createVariable('building_id', 'i4', ('building_id',))
        self.building_id.long_name = 'id of buildings that emit anthropogenic heat'
        self.building_id.units = '1'
        self.building_id[:] = np.arange(1, (self.nbuildings)+1, 1)

        self.nc_file.createDimension('street_id', self.nstreets)
        self.street_id = self.nc_file.createVariable('street_id', 'i4', ('street_id',))
        self.street_id.long_name = 'id of streets that emit anthropogenic heat'
        self.street_id.units = '1'
        self.street_id[:] = np.arange(0, (self.nstreets), 1)

        self.nc_file.createDimension('point_id', self.npoints)
        self.point_id = self.nc_file.createVariable('point_id', 'i4', ('point_id',))
        self.point_id.long_name = 'id of points that emit anthropogenic heat'
        self.point_id.units = '1'
        self.point_id[:] = np.arange(0, (self.npoints), 1)

        self.nc_file.createDimension('time', self.nt)
        self.time = self.nc_file.createVariable('time', 'f4', ('time',))
        self.time.long_name = 'time'
        self.time.standard_name = 'time'
        self.time.units = 'seconds'
        self.time.axis = 'T'
        self.time[:] = np.arange(0, (self.nt), 1) * 3600.0


    def define_variables(self):
        """Define variables for the ah driver.

        Variables give anthropogenic heat emissions per building, street or 
        point depending on time and respective id.
        """

        print("Writing variables...")

        # Define variables
        # -----------------
        nc_building_ah = self.nc_file.createVariable(
            'building_ah', 'i4', ('time', 'building_id'), fill_value=-9999)
        nc_building_ah.long_name = "anthropogenic heat emission from building roof"
        nc_building_ah.units = "kWh"
        nc_building_ah[:, :] = nc_building_ah._FillValue

        nc_street_ah = self.nc_file.createVariable(
            'street_ah', 'i4', ('time', 'street_id'), fill_value=-9999)
        nc_street_ah.long_name = "anthropogenic heat emission from street segments"
        nc_street_ah.units = "kWh"
        nc_street_ah[:, :] = nc_street_ah._FillValue

        nc_point_ah = self.nc_file.createVariable(
            'point_ah', 'i4', ('time', 'point_id'), fill_value=-9999)
        nc_point_ah.long_name = "anthropogenic heat emission from a specific ground tile"
        nc_point_ah.units = "kWh"
        nc_point_ah[:, :] = nc_point_ah._FillValue

        nc_point_x = self.nc_file.createVariable(
            'point_x', 'f8', ('point_id',), fill_value=-9999)
        nc_point_x.long_name = "x coordinate from point source (UTM coordinates)"
        nc_point_x.units = "m"
        nc_point_x[:] = nc_point_x._FillValue

        nc_point_y = self.nc_file.createVariable(
            'point_y', 'f8', ('point_id',), fill_value=-9999)
        nc_point_y.long_name = "y coordinate from point source (UTM coordinates)"
        nc_point_y.units = "m"
        nc_point_y[:] = nc_point_y._FillValue

        # AH emissions
        # --------------
        nc_building_ah[:, 0] = 1000
        nc_building_ah[:, 1] = 2000

        nc_street_ah[:, 0] = 100
        nc_street_ah[:, 1] = 200
        nc_street_ah[:, 2] = 300

        nc_point_ah[:, 0] = 10
        nc_point_ah[:, 1] = 20
        nc_point_ah[:, 2] = 30
        nc_point_ah[:, 3] = 40
        nc_point_ah[:, 4] = 50

        # Define coordinates of emission points at the grid centres and assume 2 m grid spacing.
        # (coordinates are, however, saved as UTM coordinates)
        dx = dy = 2.0

        nc_point_x[:] = (np.array([4, 9, 25, 22, 26]) + 0.5) * dx + self.nc_file.origin_x
        nc_point_y[:] = (np.array([4, 9,  2, 13, 25]) + 0.5) * dy + self.nc_file.origin_y


    def finalize(self):
        """Close file."""
        print("Closing file...")

        self.nc_file.close()


if __name__ == '__main__':
    driver = AHDriver()
    driver.write_global_attributes()
    driver.define_dimensions()
    driver.define_variables()
    driver.finalize()
