#!/usr/bin/env python3
# -------------------------------------------------------------------- #
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2021 Leibniz Universitaet Hannover
# -------------------------------------------------------------------- #

import datetime
import math

from netCDF4 import Dataset
import numpy as np


class StaticDriver:
    """This script generates a static driver for the anthropogenic heat 
    emission test case.
    """

    def __init__(self):
        """Open the static driver as netCDF4 file.
        """
        print('Opening file...')
        self.nc_file = Dataset('ah_emission_test_static', 'w', format='NETCDF4')

    def write_global_attributes(self):
        """Write global attributes to static driver."""
        print("Writing global attributes...")

        # Optional global attributes
        # --------------------------
        self.nc_file.title = 'PALM static driver for AHE module test case'
        self.nc_file.author = 'Tobias Gronemeier, Mathias Niffeler'
        self.nc_file.institution = 'iMA Richter Roeckle GmbH Co KG, ' \
            'Singapore ETH Centre'
        self.nc_file.comment = 'Generic city setup; coordinates use EPSG:3414'
        self.nc_file.creation_date = datetime.datetime.utcnow().strftime('%y-%m-%d %H:%M:%S %z')
        self.nc_file.history = ''
        self.nc_file.keywords = 'test case, PALM-4U'
        self.nc_file.license = ''
        self.nc_file.palm_version = ''
        self.nc_file.references = ''
        self.nc_file.source = ''
        self.nc_file.version = '1'

        # Mandatory global attributes
        # ---------------------------
        self.nc_file.Conventions = 'CF-1.7'
        self.nc_file.origin_lat = 1.3
        self.nc_file.origin_lon = 103.8
        self.nc_file.origin_time = '2019-03-06 12:00:00 +08'
        self.nc_file.origin_x = 363594.0
        self.nc_file.origin_y = 144130.0
        self.nc_file.origin_z = 0.0
        self.nc_file.rotation_angle = 0.0

    def define_dimensions(self):
        """Set dimensions on which variables are defined."""
        print("Writing dimensions...")

        # Specify general grid parameters
        # These values must equal to those set in the initialization_parameters
        self.nx = 31
        self.ny = 31
        self.nz = 60
        dx = 2
        dy = 2
        dz = 2

        # Create soil grid (only relevant if land surface module is used)
        dz_soil = np.array((0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86))
        zsoil_fullLayers = np.zeros_like(dz_soil)
        zsoil_fullLayers = np.around(
            [np.sum(dz_soil[:zs]) for zs in np.arange(1, len(dz_soil)+1)], 2)
        zsoil_array = zsoil_fullLayers - dz_soil/2.

        # Coordinates
        # -----------
        self.nc_file.createDimension('x', self.nx+1)
        self.x = self.nc_file.createVariable('x', 'f4', ('x',))
        self.x.long_name = 'distance to origin in x-direction'
        self.x.units = 'm'
        self.x.axis = 'X'
        self.x[:] = np.arange(0, (self.nx+1)*dx, dx) + 0.5 * dx

        self.nc_file.createDimension('y', self.ny+1)
        self.y = self.nc_file.createVariable('y', 'f4', ('y',))
        self.y.long_name = 'distance to origin in y-direction'
        self.y.units = 'm'
        self.y.axis = 'Y'
        self.y[:] = np.arange(0, (self.ny+1)*dy, dy) + 0.5 * dy

        z_array = np.append(0, np.arange(dz/2, (self.nz)*dz, dz))
        self.nc_file.createDimension('z', self.nz+1)
        self.z = self.nc_file.createVariable('z', 'f4', ('z',))
        self.z.long_name = 'height above origin'
        self.z.units = 'm'
        self.z.axis = 'Z'
        self.z.positive = 'up'
        self.z[:] = z_array

        zlad_array = self.z[:6]
        self.nc_file.createDimension('zlad', len(zlad_array))
        self.zlad = self.nc_file.createVariable('zlad', 'f4', ('zlad',))
        self.zlad.long_name = 'height above ground'
        self.zlad.units = 'm'
        self.zlad.axis = 'Z'
        self.zlad.positive = 'up'
        self.zlad[:] = zlad_array

        self.nc_file.createDimension('nsurface_fraction', 3)
        self.nsurface_fraction = self.nc_file.createVariable(
            'nsurface_fraction', 'i4', ('nsurface_fraction',))
        self.nsurface_fraction[:] = np.arange(3)


    def define_variables(self):
        """Define variables for the static driver.

        building_2d_array = np.ones((self.ny+1, self.nx+1)) * -9999.0
        south_wall, north_wall, left_wall, right_wall = 20, 25, 20, 25
        building_2d_array[
            south_wall:north_wall,
            left_wall:right_wall
            ] = 50
        nc_buildings_2d = self.nc_file.createVariable(
            'buildings_2d', 'f4', ('y', 'x'), fill_value=-9999.0)
        nc_buildings_2d.lod = 1
        nc_buildings_2d[:, :] = building_2d_array
        """
        print("Writing variables...")

        # Topography set-up
        # -----------------
        nc_building_id = self.nc_file.createVariable(
            'building_id', 'i4', ('y', 'x'), fill_value=-9999)
        nc_building_id.long_name = "building id number"
        nc_building_id.units = "1"
        nc_building_id[:, :] = nc_building_id._FillValue

        nc_buildings_2d = self.nc_file.createVariable(
            'buildings_2d', 'f4', ('y', 'x'), fill_value=-9999.0)
        nc_buildings_2d.long_name = "building height"
        nc_buildings_2d.units = "m"
        nc_buildings_2d.lod = np.int32(1)
        nc_buildings_2d[:, :] = nc_buildings_2d._FillValue

        nc_zt = self.nc_file.createVariable(
            'zt', 'f4', ('y', 'x'), fill_value=-9999.0)
        nc_zt.long_name = 'terrain height'
        nc_zt.units = 'm'
        nc_zt[:, :] = nc_zt._FillValue

        # Main surface clasification
        # --------------------------
        nc_building_type = self.nc_file.createVariable(
            'building_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_building_type.long_name = "building type classification"
        nc_building_type.units = "1"
        nc_building_type[:, :] = nc_building_type._FillValue

        nc_pavement_type = self.nc_file.createVariable(
            'pavement_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_pavement_type.long_name = "pavement type classification"
        nc_pavement_type.units = "1"
        nc_pavement_type[:, :] = nc_pavement_type._FillValue

        nc_soil_type = self.nc_file.createVariable(
            'soil_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_soil_type.long_name = "soil type classification"
        nc_soil_type.units = "1"
        nc_soil_type.lod = np.int32(1)
        nc_soil_type[:, :] = nc_soil_type._FillValue

        # nc_street_type = self.nc_file.createVariable(
        #     'street_type', 'i1', ('y', 'x'),
        #     fill_value=-127)
        # nc_street_type.long_name = "street type classification"
        # nc_street_type.units = "1"
        # nc_street_type[:, :] = nc_street_type._FillValue

        nc_surface_fraction = self.nc_file.createVariable(
            'surface_fraction', 'f4', ('nsurface_fraction', 'y', 'x'), fill_value=-9999.0)
        nc_surface_fraction.long_name = "surface fraction"
        nc_surface_fraction.units = "1"
        nc_surface_fraction[0, :, :] = nc_surface_fraction._FillValue  # vegetation fraction
        nc_surface_fraction[1, :, :] = nc_surface_fraction._FillValue  # pavement fraction
        nc_surface_fraction[2, :, :] = nc_surface_fraction._FillValue  # water fraction

        nc_vegetation_type = self.nc_file.createVariable(
            'vegetation_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_vegetation_type.long_name = "vegetation type classification"
        nc_vegetation_type.units = "1"
        nc_vegetation_type[:, :] = nc_vegetation_type._FillValue

        nc_water_type = self.nc_file.createVariable(
            'water_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_water_type.long_name = "water type classification"
        nc_water_type.units = "1"
        nc_water_type[:, :] = nc_water_type._FillValue

        # Vegetation parameters
        # ---------------------
        nc_lad = self.nc_file.createVariable(
            'lad', 'f4', ('zlad', 'y', 'x'), fill_value=-9999.0)
        nc_lad.long_name = "leaf area density"
        nc_lad.units = "m2 m-3"
        nc_lad[:, :, :] = nc_lad._FillValue


        # Set topography
        # --------------
        nc_zt[:, :] = 0.0

        # Set buildings
        # --------------
        nc_building_id[:8, 12:16] = 1
        nc_building_id[4:8, 16:20] = 1
        nc_building_id[8:18, 20:28] = 2

        nc_buildings_2d[:, :] = np.where(
            nc_building_id[:, :] > 0,
            nc_building_id[:, :] * 10,
            nc_buildings_2d[:, :])

        nc_building_type[:, :] = np.where(
            nc_building_id[:, :] > nc_building_id._FillValue,
            nc_building_id[:, :],
            nc_building_type[:, :])

        # Set surface types
        # -----------------
        nc_pavement_type[8:11, :20] = 1
        nc_pavement_type[11:, 16:20] = 2

        for i in range(0, 7):
            try:
                nc_pavement_type[25+i, 20+i] = 3
            except:
                continue
            try:
                nc_pavement_type[26+i, 20+i] = 3
            except:
                continue
            try:
                nc_pavement_type[27+i, 20+i] = 3
            except:
                continue

        nc_vegetation_type[:, :] = np.where(
            nc_pavement_type[:, :] > nc_pavement_type._FillValue,
            nc_vegetation_type[:, :],
            3)

        nc_soil_type[:, :] = np.where(
            nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
            1,
            nc_soil_type[:, :])
        nc_soil_type[:, :] = np.where(
            nc_pavement_type[:, :] > nc_pavement_type._FillValue,
            1,
            nc_soil_type[:, :])

        nc_surface_fraction[:, :, :] = 0
        nc_surface_fraction[0, :, :] = np.where(
            nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
            1,
            nc_surface_fraction[0, :, :])
        nc_surface_fraction[1, :, :] = np.where(
            nc_pavement_type[:, :] > nc_pavement_type._FillValue,
            1,
            nc_surface_fraction[1, :, :])
        nc_surface_fraction[2, :, :] = np.where(
            nc_water_type[:, :] > nc_water_type._FillValue,
            1,
            nc_surface_fraction[2, :, :])

        # Set vegetation
        # --------------
        lad_profile = [0.0, 0.01070122, 0.1070122, 0.3130108, 0.3879193, 0.1712195]
        for k in range(len(self.zlad)):
            nc_lad[k, 23:25, 7:9] = lad_profile[k]


    def finalize(self):
        """Close file."""
        print("Closing file...")

        self.nc_file.close()


if __name__ == '__main__':
    driver = StaticDriver()
    driver.write_global_attributes()
    driver.define_dimensions()
    driver.define_variables()
    driver.finalize()
