#!/usr/bin/env python3
# ------------------------------------------------------------------------------ #
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2021  Leibniz Universitaet Hannover
# ------------------------------------------------------------------------------ #
#
# Description:
# ------------
# Processing tool for creating PIDS conform static drivers from rastered netCDF
# input
#
# @Author Bjoern Maronga (maronga@muk.uni-hannover.de)
#
# @TODO Make input files optional
# @TODO Allow for ASCII input of terrain height and building height
# @TODO Modularize reading config file
# @TODO Convert to object-oriented treatment (buidings, trees)
# @TODO Automatically shift child domains so that their origin lies intersects
#       a edge note of the parent grid
# ------------------------------------------------------------------------------ #

import math
from typing import Dict, List

import numpy as np
import numpy.ma as ma

from palm_csd.canopy_generator import process_patch
from palm_csd.csd_config import CSDConfig, CSDConfigSettings
from palm_csd.csd_domain import CSDDomain
from palm_csd.tools import (
    blend_array_2d,
    height_to_z_grid,
    check_arrays_2,
    check_consistency_3,
    check_consistency_4,
    interpolate_2d,
    ma_isin,
    make_3d_from_2d,
    make_3d_from_bridges_2d,
)
from palm_csd.vegetation import DomainTree

############################################################


def create_driver(input_configuration_file: str) -> None:
    """Main routine for creating the static driver."""

    defaultvalues = {
        "lat": float(-9999.0),
        "lon": float(-9999.0),
        "E_UTM": float(-9999.0),
        "N_UTM": float(-9999.0),
        "zt": float(0.0),
        "buildings_2d": float(0.0),
        "buildings_3d": 0,
        "bridges_2d": float(0.0),
        "building_id": int(0),
        "bridges_id": int(0),
        "building_type": 1,
        "nsurface_fraction": int(-9999),
        "vegetation_type": 3,
        "vegetation_height": float(-9999.0),
        "pavement_type": 1,
        "water_type": 1,
        "street_type": 1,
        "street_crossings": 0,
        "soil_type": 1,
        "surface_fraction": float(0.0),
        "buildings_pars": float(-9999.0),
        "tree_data": float(-9999.0),
        "tree_type": np.byte(-127),
        "patch_type": 0,
        "vegetation_pars": float(-9999.0),
        "water_pars": float(-9999.0),
        "water_temperature": 283.0,
    }

    #  vegetation_types that are considered high vegetation
    vt_high_vegetation = [
        4,  # evergreen needleleaf trees
        5,  # deciduous needleleaf trees
        6,  # evergreen broadleaf trees
        7,  # deciduous broadleaf trees
        17,  # mixed forest/woodland
        18,  # interrupted forest
    ]

    # Read configuration file and set parameters accordingly
    config = CSDConfig(input_configuration_file)

    def add_domain_and_parents(name: str, domains: Dict[str, CSDDomain]) -> None:
        if name not in domains:
            parent = None
            parent_name = config.domain_dict[name].domain_parent
            if parent_name is not None:
                add_domain_and_parents(parent_name, domains)
                parent = domains[parent_name]
            domains[name] = CSDDomain(name, config, parent)

    # Create domains
    domains: Dict[str, CSDDomain] = {}
    for name in config.domain_dict:
        add_domain_and_parents(name, domains)

    # Initialize domain tree's defaults
    DomainTree.populate_defaults()

    # Minium of terrain height. This value will be substracted for all data sets
    zt_min = minimum_terrain_height_global(domains)
    print("Shift terrain heights by -" + str(zt_min))

    # Loop over domains, domains are independent of each other except zt_min calculated above
    for domain in domains.values():
        domain.remove_existing_output()

        process_coordinates(domain, zt_min)

        process_buildings_bridges(domain, config.settings, defaultvalues)

        process_types(domain, defaultvalues, vt_high_vegetation)
        process_street_type_crossing(domain, defaultvalues)

        if domain.config.vegetation_on_roofs:
            process_vegetation_roof(domain, config.settings)
        process_lai(domain)
        if domain.config.street_trees:
            process_single_trees(domain, config.settings, vt_high_vegetation)
        if domain.config.generate_vegetation_patches:
            process_vegetation_patches(domain, config.settings, defaultvalues, vt_high_vegetation)
        process_lai_final(domain)

        if (
            domain.config.water_temperature_per_water_type is not None
            or domain.input_config.file_water_temperature is not None
        ):
            process_water_temperature(domain, defaultvalues)

        consistency_check_update_surface_fraction(domain)
        domain.write_global_attributes()


def minimum_terrain_height_global(domains: Dict[str, CSDDomain]) -> float:
    """Calculate minimum terrain height of given domains."""

    zt_min = math.inf
    for domain in domains.values():
        # Read terrain height
        zt = domain.read_from_file_2d(domain.input_config.file_zt)
        zt_min = min(zt_min, min(zt.flatten()))

    return zt_min


def process_coordinates(domain: CSDDomain, zt_min: float) -> None:
    """Process coordinates and terrain height of a given domain.
    z_min is subtracted from the terrain height. Store the results.
    """

    # Use origin_x and origin_y to calculate UTM and lon/lat coordinates
    if domain.geo_converter is not None:
        domain.origin_x = domain.geo_converter.origin_x
        domain.origin_y = domain.geo_converter.origin_y
        domain.origin_lon = domain.geo_converter.origin_lon
        domain.origin_lat = domain.geo_converter.origin_lat

        # Global x and y coordinates (cell centre) relative to root parent domain
        # Used only for zt interpolation below
        (
            domain.x_global.values,
            domain.y_global.values,
        ) = domain.geo_converter.global_palm_coordinates()

        # Coordinates
        e_UTM, n_UTM, lon, lat = domain.geo_converter.geographic_coordinates()

        # Write CRS
        domain.write_crs_to_file()

    else:
        # Get coordinates near origin
        x_UTM_origin = domain.read_from_file_2d(
            domain.input_config.file_x_UTM,
            x0=domain.config.x0,
            x1=domain.config.x0 + 1,
            y0=domain.config.y0,
            y1=domain.config.y0 + 1,
        )
        y_UTM_origin = domain.read_from_file_2d(
            domain.input_config.file_y_UTM,
            x0=domain.config.x0,
            x1=domain.config.x0 + 1,
            y0=domain.config.y0,
            y1=domain.config.y0 + 1,
        )
        lat_origin = domain.read_from_file_2d(
            domain.input_config.file_lat,
            x0=domain.config.x0,
            x1=domain.config.x0 + 1,
            y0=domain.config.y0,
            y1=domain.config.y0 + 1,
        )
        lon_origin = domain.read_from_file_2d(
            domain.input_config.file_lon,
            x0=domain.config.x0,
            x1=domain.config.x0 + 1,
            y0=domain.config.y0,
            y1=domain.config.y0 + 1,
        )

        # Calculate position of origin. Added as global attributes later
        domain.origin_x = float(x_UTM_origin[0, 0]) - 0.5 * (
            float(x_UTM_origin[0, 1]) - float(x_UTM_origin[0, 0])
        )
        domain.origin_y = float(y_UTM_origin[0, 0]) - 0.5 * (
            float(y_UTM_origin[1, 0]) - float(y_UTM_origin[0, 0])
        )
        domain.origin_lon = float(lon_origin[0, 0]) - 0.5 * (
            float(lon_origin[0, 1]) - float(lon_origin[0, 0])
        )
        domain.origin_lat = float(lat_origin[0, 0]) - 0.5 * (
            float(lat_origin[1, 0]) - float(lat_origin[0, 0])
        )

        # Read x and y values
        domain.x_global.values = domain.read_from_file_1d(domain.input_config.file_x_UTM, "x")
        domain.y_global.values = domain.read_from_file_1d(
            domain.input_config.file_y_UTM, "y", x0=domain.config.y0, x1=domain.config.y1
        )

        # Read and write lon, lat and UTM coordinates
        lat = domain.read_from_file_2d(domain.input_config.file_lat)
        lon = domain.read_from_file_2d(domain.input_config.file_lon)

        e_UTM = domain.read_from_file_2d(domain.input_config.file_x_UTM)
        n_UTM = domain.read_from_file_2d(domain.input_config.file_y_UTM)

        # Write CRS
        crs = domain.read_from_file_crs()
        crs.to_file()

    # Shift x and y coordinates for x and y local cell centre coordinates of domain
    # Used as output dimensions
    domain.x.values = (
        domain.x_global.values
        - min(domain.x_global.values.flatten())
        + domain.config.pixel_size / 2.0
    )
    domain.y.values = (
        domain.y_global.values
        - min(domain.y_global.values.flatten())
        + domain.config.pixel_size / 2.0
    )

    domain.lat.to_file(lat)
    domain.lon.to_file(lon)

    domain.E_UTM.to_file(e_UTM)
    domain.N_UTM.to_file(n_UTM)

    # Read, process and write terrain height (zt)
    domain.zt.values = domain.read_from_file_2d(domain.input_config.file_zt)
    domain.zt.values = domain.zt.values - zt_min
    domain.origin_z = float(zt_min)

    # If necessary, interpolate parent domain terrain height on child domain grid and blend
    # the two
    if domain.config.interpolate_terrain:
        if domain.parent is None:
            raise Exception("Interpolation of terrain height requires a parent domain")
        if domain.parent.x_global.values is None or domain.parent.y_global.values is None:
            raise Exception(f"x_UTM or y_UTM of parent {domain.parent.name} not calculated")
        if domain.parent.zt.values is None:
            raise Exception(f"zt of parent {domain.parent.name} not calculated")

        tmp_x0 = np.searchsorted(domain.parent.x_global.values, domain.x_global.values[0]) - 1
        tmp_y0 = np.searchsorted(domain.parent.y_global.values, domain.y_global.values[0]) - 1
        tmp_x1 = np.searchsorted(domain.parent.x_global.values, domain.x_global.values[-1]) + 1
        tmp_y1 = np.searchsorted(domain.parent.y_global.values, domain.y_global.values[-1]) + 1

        if tmp_x0 < 0:
            raise Exception(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the left border"
            )
        if tmp_y0 < 0:
            raise Exception(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the bottom border"
            )
        if tmp_x1 > domain.parent.x_global.values.shape[0]:
            raise Exception(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the right border"
            )
        if tmp_y1 > domain.parent.y_global.values.shape[0]:
            raise Exception(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the top border"
            )

        tmp_x = domain.parent.x_global.values[tmp_x0:tmp_x1]
        tmp_y = domain.parent.y_global.values[tmp_y0:tmp_y1]

        zt_parent = domain.parent.zt.values[tmp_y0:tmp_y1, tmp_x0:tmp_x1]

        # Interpolate array and bring to PALM grid of child domain
        zt_ip = interpolate_2d(
            zt_parent, tmp_x, tmp_y, domain.x_global.values, domain.y_global.values
        )
        zt_ip = height_to_z_grid(zt_ip, domain.parent.config.dz)

        # Shift the child terrain height according to the parent mean terrain height
        print(
            f"Shifting terrain height of domain {domain.name}: -{np.mean(domain.zt.values)} "
            + f"+{np.mean(zt_ip)}"
        )
        domain.zt.values = domain.zt.values - np.mean(domain.zt.values) + np.mean(zt_ip)

        # Blend over the parent and child terrain height within a radius of 50 px (or less if
        # domain is smaller than 50 px)
        domain.zt.values = blend_array_2d(
            domain.zt.values, zt_ip, min(50, min(domain.zt.values.shape) * 0.5)
        )

    # If necessary, bring terrain height to PALM's vertical grid. This is either forced by
    # the user or implicitly by using interpolation for a child domain
    if domain.config.use_palm_z_axis:
        domain.zt.values = height_to_z_grid(domain.zt.values, domain.config.dz)

    domain.zt.to_file()  # zt values stored directly on NCDFVariable


def process_buildings_bridges(
    domain: CSDDomain, settings: CSDConfigSettings, defaultvalues: Dict[str, float]
) -> None:
    """Process buildings and bridges of a given domain. Store the results."""

    buildings_2d = domain.read_from_file_2d(domain.input_config.file_buildings_2d)
    building_id = domain.read_from_file_2d(domain.input_config.file_building_id)

    building_type = domain.read_from_file_2d(domain.input_config.file_building_type)
    ma.masked_greater_equal(building_type, 254, copy=False)
    building_type = ma.where(building_type < 1, defaultvalues["building_type"], building_type)

    check = check_arrays_2(buildings_2d, building_id)
    # make masks equal
    if not check:
        buildings_2d.mask = ma.mask_or(building_id.mask, buildings_2d.mask)
        # copy mask from building_2d to building_id
        building_id.mask = buildings_2d.mask.copy()
        print("Data check #1 " + str(check_arrays_2(buildings_2d, building_id)))

    check = check_arrays_2(buildings_2d, building_type)
    if not check:
        building_type.mask = ma.mask_or(buildings_2d.mask, building_type.mask)
        building_type = ma.where(
            building_type.mask & ~buildings_2d.mask, defaultvalues["building_type"], building_type
        )
        print("Data check #2 " + str(check_arrays_2(buildings_2d, building_type)))

    domain.buildings_2d.to_file(buildings_2d)
    domain.building_id.to_file(building_id)
    domain.building_type.to_file(building_type)

    # Create 3d buildings if necessary. In that course, read bridge objects and add them to
    # building layer
    if domain.config.buildings_3d:
        bridges_2d = domain.read_from_file_2d(domain.input_config.file_bridges_2d)
        bridges_id = domain.read_from_file_2d(domain.input_config.file_bridges_id)

        ma.masked_equal(bridges_2d, 0.0, copy=False)
        building_id = ma.where(bridges_2d.mask, building_id, bridges_id)

        if ma.any(~buildings_2d.mask):
            buildings_3d, z = make_3d_from_2d(buildings_2d, domain.x, domain.y, domain.config.dz)
            if ma.any(~bridges_2d.mask):
                buildings_3d = make_3d_from_bridges_2d(
                    buildings_3d,
                    bridges_2d,
                    domain.x,
                    domain.y,
                    domain.config.dz,
                    settings.bridge_width,
                )
            else:
                print("Skipping creation of 3D bridges (no bridges in domain)")

            domain.z.values = z
            domain.building_id.to_file(building_id)
            domain.buildings_3d.to_file(buildings_3d)

        else:
            print("Skipping creation of 3D buildings (no buildings in domain)")


def process_types(
    domain: CSDDomain, defaultvalues: Dict[str, float], vt_high_vegetation: List[int]
) -> None:
    """Read vegetation type, water_type, pavement_type and soil_type and make fields consistent.
    Store the results.
    """

    building_type = domain.building_type.from_file()

    vegetation_type = domain.read_from_file_2d(domain.input_config.file_vegetation_type)
    ma.masked_equal(vegetation_type, 255, copy=False)
    vegetation_type = ma.where(
        vegetation_type < 1, defaultvalues["vegetation_type"], vegetation_type
    )

    pavement_type = domain.read_from_file_2d(domain.input_config.file_pavement_type)
    ma.masked_equal(pavement_type, 255, copy=False)
    pavement_type = ma.where(pavement_type < 1, defaultvalues["pavement_type"], pavement_type)

    water_type = domain.read_from_file_2d(domain.input_config.file_water_type)
    ma.masked_equal(water_type, 255, copy=False)
    water_type = ma.where(water_type < 1, defaultvalues["water_type"], water_type)

    soil_type = domain.read_from_file_2d(domain.input_config.file_soil_type)
    ma.masked_equal(soil_type, 255, copy=False)
    soil_type = ma.where(soil_type < 1, defaultvalues["soil_type"], soil_type)
    # use default values for masked values, made consistent below
    soil_type = ma.where(soil_type.mask, defaultvalues["soil_type"], soil_type)

    # Make arrays consistent
    # #1 Set vegetation type to masked for pixel where a pavement type is set
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~pavement_type.mask)

    # #2 Set vegetation type to masked for pixel where a building type is set
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~building_type.mask)

    # #3 Set vegetation type to masked for pixel where a water type is set
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~water_type.mask)

    # #4 Remove pavement for pixels with buildings
    pavement_type.mask = ma.mask_or(pavement_type.mask, ~building_type.mask)

    # #5 Remove pavement for pixels with water.
    pavement_type.mask = ma.mask_or(pavement_type.mask, ~water_type.mask)

    # #6 Remove water for pixels with buildings
    water_type.mask = ma.mask_or(water_type.mask, ~building_type.mask)

    # Correct vegetation_type when a vegetation height is available and is indicative of low
    # vegetation
    vegetation_height = domain.read_from_file_2d(domain.input_config.file_vegetation_height)

    # correct vegetation_type depending on vegetation_height
    # ma.where gives ma.masked when its first argument is ma.masked. We don't want this for
    # vegetation_height here so do an extra check and use .data
    vegetation_type = ma.where(
        (~vegetation_height.mask)
        & (vegetation_height.data == 0.0)
        & ma_isin(vegetation_type, vt_high_vegetation),
        3,
        vegetation_type,
    )
    ma.masked_where(
        (vegetation_height == 0.0) & ma_isin(vegetation_type, vt_high_vegetation),
        vegetation_height,
        copy=False,
    )

    # Check for consistency and fill empty fields with default vegetation type
    consistency_array, test = check_consistency_4(
        vegetation_type, building_type, pavement_type, water_type
    )

    if test:
        vegetation_type = ma.where(
            consistency_array == 0, defaultvalues["vegetation_type"], vegetation_type
        )
        consistency_array, test = check_consistency_4(
            vegetation_type, building_type, pavement_type, water_type
        )

    # #7 Todo: remove soil_type for pixels with vegetation_type and pavement_type
    soil_type = ma.where(vegetation_type.mask & pavement_type.mask, ma.masked, soil_type)

    # Check for consistency and fill empty fields with default vegetation type
    consistency_array, test = check_consistency_3(vegetation_type, pavement_type, soil_type)

    # Create surface_fraction array
    domain.nsurface_fraction.values = ma.arange(0, 3)
    surface_fraction = ma.ones((domain.nsurface_fraction.size, domain.y.size, domain.x.size))

    surface_fraction[0, :, :] = ma.where(vegetation_type.mask, 0.0, 1.0)
    surface_fraction[1, :, :] = ma.where(pavement_type.mask, 0.0, 1.0)
    surface_fraction[2, :, :] = ma.where(water_type.mask, 0.0, 1.0)

    domain.surface_fraction.to_file(surface_fraction)
    domain.vegetation_type.to_file(vegetation_type)
    domain.pavement_type.to_file(pavement_type)
    domain.water_type.to_file(water_type)
    domain.soil_type.to_file(soil_type)

    # pixels with bridges get building_type = 7 = bridge. This does not change the _type
    # setting for the under-bridge area
    # NOTE: when bridges are present the consistency check will fail at the moment
    if domain.config.buildings_3d:
        if np.any(~building_type.mask):
            bridges_2d = domain.read_from_file_2d(domain.input_config.file_bridges_2d)
            ma.masked_equal(bridges_2d, 0.0, copy=False)
            building_type = ma.where(bridges_2d.mask, building_type, 7)
            domain.building_type.to_file(building_type)


def process_street_type_crossing(domain: CSDDomain, defaultvalues: Dict[str, float]) -> None:
    """Process street type and street crossings of a given domain. Store the results."""

    street_type = domain.read_from_file_2d(domain.input_config.file_street_type)
    ma.masked_equal(street_type, 255, copy=False)
    street_type = ma.where(street_type < 1, defaultvalues["street_type"], street_type)

    pavement_type = domain.pavement_type.from_file()
    street_type.mask = ma.mask_or(pavement_type.mask, street_type.mask)

    domain.street_type.to_file(street_type)

    street_crossings = domain.read_from_file_2d(domain.input_config.file_street_crossings)
    ma.masked_equal(street_crossings, 255, copy=False)
    street_crossings = ma.where(
        street_crossings < 1, defaultvalues["street_crossings"], street_crossings
    )

    domain.street_crossing.to_file(street_crossings)


def process_vegetation_roof(domain: CSDDomain, settings: CSDConfigSettings) -> None:
    """Process vegetation on roofs of a given domain. Store the results."""

    green_roofs = domain.read_from_file_2d(domain.input_config.file_vegetation_on_roofs)
    buildings_2d = domain.buildings_2d.from_file()

    domain.nbuilding_pars.values = ma.arange(0, 150)
    building_pars = ma.ones((domain.nbuilding_pars.size, domain.y.size, domain.x.size))
    building_pars[:, :, :] = ma.masked

    # assign green fraction on roofs
    building_pars[3, :, :] = ma.where((~buildings_2d.mask) & (green_roofs != 0.0), 1, ma.masked)

    # set wall fraction to 0 where green fraction defined
    building_pars[89, :, :] = ma.where(
        ~(building_pars.mask[3, :, :]),
        0.0,
        ma.masked,
    )

    # set window fraction to 0 where green fraction defined
    building_pars[102, :, :] = ma.where(
        ~(building_pars.mask[3, :, :]),
        0.0,
        ma.masked,
    )

    # assign leaf area index for vegetation on roofs
    building_pars[4, :, :] = ma.where(
        (~(building_pars.mask[3, :, :])) & (green_roofs >= 0.5),
        settings.lai_roof_intensive,
        ma.masked,
    )
    building_pars[4, :, :] = ma.where(
        (~(building_pars.mask[3, :, :])) & (green_roofs < 0.5),
        settings.lai_roof_extensive,
        building_pars[4, :, :],
    )

    domain.building_pars.to_file(building_pars)


def process_lai(domain: CSDDomain) -> None:
    """Process LAI of a given domain. Store the results."""

    lai = domain.read_from_file_2d(domain.input_config.file_lai)

    vegetation_type = domain.vegetation_type.from_file()

    lai.mask = ma.mask_or(vegetation_type.mask, lai.mask)

    domain.nvegetation_pars.values = ma.arange(0, 12)
    vegetation_pars = ma.ones((domain.nvegetation_pars.size, domain.y.size, domain.x.size))
    vegetation_pars[:, :, :] = ma.masked

    vegetation_pars[1, :, :] = lai

    # Write out first version of LAI. Will later be overwritten.
    domain.vegetation_pars.to_file(vegetation_pars)


def process_single_trees(
    domain: CSDDomain, settings: CSDConfigSettings, vt_high_vegetation
) -> None:
    """Read tree data and create LAD and BAD arrays using the canopy generator for a given domain.
    Store the results.
    """

    lai = domain.read_from_file_2d(domain.input_config.file_lai)

    # Read all tree parameters from file. They are defined at the centre of the tree.
    # Data correction and modification is done in generate_tree below
    tree_height_centre = domain.read_from_file_2d(domain.input_config.file_tree_height)
    tree_crown_diameter_centre = domain.read_from_file_2d(
        domain.input_config.file_tree_crown_diameter
    )
    tree_trunk_diameter_centre = domain.read_from_file_2d(
        domain.input_config.file_tree_trunk_diameter
    )
    tree_type_centre = domain.read_from_file_2d(domain.input_config.file_tree_type)

    # Check, if trunk diameter is too large (e.g. being in cm as in an older version)
    # TODO: develop general input data check and remove this here
    max_tree_trunk_diameter = ma.max(tree_trunk_diameter_centre)
    if max_tree_trunk_diameter > 14.0:
        raise ValueError(
            f"Input trunk diameter is too large with a maximum value of {max_tree_trunk_diameter}. "
            + "Note that the input unit has changed from cm to m recently."
        )

    # Read patch height for general height of resolve vegetation, not used for trees
    patch_height = domain.read_from_file_2d(domain.input_config.file_patch_height)

    # For vegetation pixel with missing height information (-1), set default patch height
    patch_height = ma.where(patch_height == -1.0, settings.patch_height_default, patch_height)

    # Compare patch height array with vegetation type and correct accordingly
    vegetation_type = domain.vegetation_type.from_file()

    # For zero-height patches, set vegetation_type to short grass and remove these pixels
    # from the patch height array
    # ma.where gives ma.masked when its first argument is ma.masked. We don't want this for
    # patch_height here so do an extra check and use .data
    vegetation_type = ma.where(
        (~patch_height.mask)
        & (patch_height.data == 0.0)
        & ma_isin(vegetation_type, vt_high_vegetation),
        3,
        vegetation_type,
    )
    ma.masked_where(
        (patch_height == 0.0) & ma_isin(vegetation_type, vt_high_vegetation),
        patch_height,
        copy=False,
    )

    max_tree_height = ma.max(tree_height_centre)
    max_patch_height = ma.max(patch_height)
    max_canopy_height = ma.max(ma.append(max_tree_height, max_patch_height))

    if max_canopy_height is ma.masked:
        print("No street trees generated in domain " + domain.config.name)
        return

    # Create array for vegetation canopy heights
    zlad = ma.arange(
        0,
        math.floor(max_canopy_height / domain.config.dz) * domain.config.dz + 2 * domain.config.dz,
        domain.config.dz,
    )
    zlad[1:] = zlad[1:] - 0.5 * domain.config.dz
    domain.zlad.values = zlad

    # Create common arrays for LAD and BAD as well as arrays for tree IDs and types
    dimensions = (len(domain.zlad), len(domain.y), len(domain.x))
    lad = ma.masked_all(dimensions)
    bad = ma.masked_all(dimensions)
    tree_id = ma.masked_all(dimensions)
    tree_type = ma.masked_all(dimensions)

    # Centre of a tree?
    tree_pixels = np.where(
        ~tree_height_centre.mask
        | ~tree_type_centre.mask
        | ~tree_crown_diameter_centre.mask
        | ~tree_trunk_diameter_centre.mask,
        True,
        False,
    )
    number_of_trees = np.sum(tree_pixels)

    if number_of_trees == 0:
        return

    # Create a DomainTree for each single tree
    print("Start generating " + str(number_of_trees) + " trees...")
    trees: List[DomainTree] = []
    # counter for tree IDs and adjusted trees
    DomainTree.reset_counter()
    for i in range(0, len(domain.x)):
        for j in range(0, len(domain.y)):
            if tree_pixels[j, i]:
                tree = DomainTree.generate_tree(
                    i=i,
                    j=j,
                    type=tree_type_centre[j, i],
                    shape=ma.masked,
                    height=tree_height_centre[j, i],
                    lai=lai[j, i],
                    crown_diameter=tree_crown_diameter_centre[j, i],
                    trunk_diameter=tree_trunk_diameter_centre[j, i],
                    config=domain.config,
                    settings=settings,
                )
                if tree is not None:
                    trees.append(tree)

    DomainTree.check_counter(domain.config)

    for tree in trees:
        tree.generate_store_3d_fields(lad, bad, tree_id, tree_type, domain.config)

    # Remove LAD volumes that are inside buildings
    if not domain.config.overhanging_trees:
        buildings_2d = domain.buildings_2d.from_file()
        for k in range(0, len(zlad)):
            ma.masked_where(~buildings_2d.mask, lad[k, :, :], copy=False)
            ma.masked_where(~buildings_2d.mask, bad[k, :, :], copy=False)
            ma.masked_where(~buildings_2d.mask, tree_id[k, :, :], copy=False)

    # Write results to file
    domain.lad.to_file(lad)
    domain.bad.to_file(bad)
    domain.tree_id.to_file(tree_id)
    domain.tree_type.to_file(tree_type)


def process_vegetation_patches(
    domain: CSDDomain,
    settings: CSDConfigSettings,
    defaultvalues: Dict[str, float],
    vt_high_vegetation: List[int],
) -> None:
    """Create vegetation patches for locations with high vegetation type of a given domain.
    Store the results.
    """

    patch_height = domain.read_from_file_2d(domain.input_config.file_patch_height)

    # For vegetation pixel with missing height information (-1), set default patch height
    patch_height = ma.where(patch_height == -1.0, settings.patch_height_default, patch_height)

    max_patch_height = ma.max(patch_height)
    if max_patch_height is not ma.masked:
        # Call canopy generator for single trees only if there is any tree height
        # available in the domain. This does not guarantee that there are street trees
        # that can be processed. This is checked in the canopy generator.
        # Load vegetation_type and lad array (at level z = 0) for re-processing
        vegetation_type = domain.vegetation_type.from_file()
        lad = domain.lad.from_file()
        tree_id = domain.tree_id.from_file()
        tree_types = domain.tree_type.from_file()
        zlad = domain.read_from_file_1d(domain.config.filename, "zlad", complete=True)
        vegetation_pars = domain.vegetation_pars.from_file()
        lai = domain.read_from_file_2d(domain.input_config.file_lai)
        patch_type_2d = domain.read_from_file_2d(domain.input_config.file_patch_type)

        # patch_type_2d: use high vegetation vegetation_type if patch_type is missing
        # Note: vegetation_type corrected above
        patch_type_2d = ma.where(
            patch_type_2d.mask & ma_isin(vegetation_type, vt_high_vegetation),
            -vegetation_type,
            patch_type_2d,
        )
        # patch_type_2d: use default value for the rest of the pixels
        patch_type_2d = ma.where(patch_type_2d.mask, defaultvalues["patch_type"], patch_type_2d)

        # Determine all pixels that do not already have an LAD but which are high
        # vegetation to a dummy value of 1.0 and remove all other pixels
        lai_high = ma.where(
            (lad.mask[0, :, :])
            & (
                ma_isin(vegetation_type, vt_high_vegetation)
                & (patch_height.mask | (patch_height.data >= domain.config.dz))
            ),
            1.0,
            ma.masked,
        )

        # Treat all pixels where short grass is defined, but where a patch_height >= dz
        # is found, as high vegetation (often the case in backyards)
        lai_high = ma.where(
            lai_high.mask
            & ~patch_height.mask
            & (patch_height.data >= domain.config.dz)
            & ~vegetation_type.mask
            & (vegetation_type.data == 3),
            1.0,
            lai_high,
        )

        # If overhanging trees are allowed, assign pixels with patch_height > dz that are
        # not included in vegetation_type to high
        if domain.config.overhanging_trees:
            lai_high = ma.where(
                lai_high.mask & ~patch_height.mask & (patch_height.data >= domain.config.dz),
                1.0,
                lai_high,
            )

        # Now, assign either the default LAI for high vegetation or keep 1.0 from the
        # lai_high array.
        lai_high = ma.where(
            ~lai_high.mask & lai.mask, settings.lai_high_vegetation_default, lai_high
        )

        # If LAI values are available in the LAI array, write them on the lai_high array
        lai_high = ma.where(~lai_high.mask & ~lai.mask, lai, lai_high)

        # Define a patch height wherever it is missing, but where a high vegetation LAI
        # was set
        patch_height = ma.where(
            ~lai_high.mask & patch_height.mask, settings.patch_height_default, patch_height
        )

        # Remove pixels where street trees were already set
        ma.masked_where(~lad.mask[0, :, :], patch_height, copy=False)

        # Remove patch heights that have no lai_high value
        ma.masked_where(lai_high.mask, patch_height, copy=False)

        # For missing LAI values, set either the high vegetation default or the low
        # vegetation default
        lai_high = ma.where(
            lai_high.mask & ~patch_height.mask & (patch_height.data > 2.0),
            settings.lai_high_vegetation_default,
            lai_high,
        )
        lai_high = ma.where(
            lai_high.mask & ~patch_height.mask & (patch_height.data <= 2.0),
            settings.lai_low_vegetation_default,
            lai_high,
        )

        if ma.max(patch_height) >= (2.0 * domain.config.dz):
            print("    start calculating LAD (this might take some time)")

            lad_patch, patch_id, patch_types, patch_nz, status = process_patch(
                domain.config.dz,
                patch_height,
                patch_type_2d,
                vegetation_type,
                max(zlad),
                lai_high,
                settings.lai_alpha,
                settings.lai_beta,
            )

            # Set negative ids for patches
            patch_id = ma.where(patch_id.mask, patch_id, -patch_id)

            # 2D loop in order to avoid memory problems with large arrays
            for iii in range(0, domain.config.nx + 1):
                for jj in range(0, domain.config.ny + 1):
                    tree_id[0 : patch_nz + 1, jj, iii] = ma.where(
                        lad.mask[0 : patch_nz + 1, jj, iii],
                        patch_id[0 : patch_nz + 1, jj, iii],
                        tree_id[0 : patch_nz + 1, jj, iii],
                    )
                    tree_types[0 : patch_nz + 1, jj, iii] = ma.where(
                        lad.mask[0 : patch_nz + 1, jj, iii],
                        patch_types[0 : patch_nz + 1, jj, iii],
                        tree_types[0 : patch_nz + 1, jj, iii],
                    )
                    lad[0 : patch_nz + 1, jj, iii] = ma.where(
                        lad.mask[0 : patch_nz + 1, jj, iii],
                        lad_patch[0 : patch_nz + 1, jj, iii],
                        lad[0 : patch_nz + 1, jj, iii],
                    )

        # Remove high vegetation wherever it is replaced by a leaf area density. This
        # should effectively remove all high vegetation pixels
        vegetation_type = ma.where(
            ~lad.mask[0, :, :] & ~vegetation_type.mask,
            settings.vegetation_type_below_trees,
            vegetation_type,
        )

        # If desired, remove all high vegetation.
        # TODO: check if this is still necessary
        if not domain.config.allow_high_vegetation:
            vegetation_type = ma.where(
                ~vegetation_type.mask & ma_isin(vegetation_type, vt_high_vegetation),
                3,
                vegetation_type,
            )

        # Set default low LAI for pixels with an LAD (short grass below trees)
        lai_low = ma.where(lad.mask[0, :, :], lai, settings.lai_low_vegetation_default)

        # Fill low vegetation pixels without LAI set or with LAI = 0 with default value
        lai_low = ma.where(
            (lai_low.mask | (lai_low == 0.0)) & ~vegetation_type.mask,
            settings.lai_low_vegetation_default,
            lai_low,
        )

        # Remove lai for pixels that have no vegetation_type
        lai_low = ma.where(~vegetation_type.mask & (vegetation_type != 1), lai_low, ma.masked)

        # Overwrite lai in vegetation_parameters
        vegetation_pars[1, :, :] = ma.copy(lai_low)
        domain.vegetation_pars.to_file(vegetation_pars)

        # Overwrite lad and id arrays
        domain.lad.to_file(lad)
        domain.tree_id.to_file(tree_id)
        domain.tree_type.to_file(tree_types)

        domain.vegetation_type.to_file(vegetation_type)

    else:
        print("No tree patches found in domain " + domain.config.name)


def process_lai_final(domain: CSDDomain) -> None:
    """Remove LAI where a bare soil was set of a given domain."""

    vegetation_type = domain.vegetation_type.from_file()
    vegetation_pars = domain.vegetation_pars.from_file()
    lai = vegetation_pars[1, :, :]

    # Remove lai for pixels that have no vegetation_type
    ma.masked_where(vegetation_type.mask | (vegetation_type == 1), lai, copy=False)

    # Overwrite lai in vegetation_parameters
    vegetation_pars[1, :, :] = ma.copy(lai)
    domain.vegetation_pars.to_file(vegetation_pars)


def process_water_temperature(domain: CSDDomain, defaultvalues: Dict[str, float]) -> None:
    """Process water temperatures of a given domain. Store the results."""

    # Read water type from output file and create water_pars
    water_type = domain.water_type.from_file()

    domain.nwater_pars.values = ma.arange(0, 7)
    water_pars = ma.masked_all((domain.nwater_pars.size, domain.y.size, domain.x.size))

    # Assign water temperature
    # First, set default value for all water surfaces
    water_pars[0, :, :] = ma.where(~water_type.mask, defaultvalues["water_temperature"], ma.masked)

    # Set specific water temperature per type as assigned in config
    if domain.config.water_temperature_per_water_type is not None:
        for (
            water_type_index,
            water_temperature,
        ) in domain.config.water_temperature_per_water_type.items():
            if water_temperature is not None:
                water_pars[0, :, :] = ma.where(
                    water_type == water_type_index, water_temperature, water_pars[0, :, :]
                )

    # Set water temperature based on input file
    if domain.input_config.file_water_temperature is not None:
        water_temperature_from_file = domain.read_from_file_2d(
            domain.input_config.file_water_temperature
        )
        water_temperature_from_file.mask = ma.mask_or(
            water_temperature_from_file.mask, water_type.mask
        )
        water_temperature_from_file = ma.where(
            (water_temperature_from_file < 265.0) & (water_temperature_from_file > 373.15),
            ma.masked,
            water_temperature_from_file,
        )
        water_pars[0, :, :] = ma.where(
            ~water_temperature_from_file.mask, water_temperature_from_file, water_pars[0, :, :]
        )

    domain.water_pars.to_file(water_pars)


def consistency_check_update_surface_fraction(domain: CSDDomain) -> None:
    """Do consistency check and update surface fractions for a given domain."""

    vegetation_type = domain.vegetation_type.from_file()
    pavement_type = domain.pavement_type.from_file()
    building_type = domain.building_type.from_file()
    water_type = domain.water_type.from_file()
    soil_type = domain.soil_type.from_file()

    # Check for consistency and fill empty fields with default vegetation type
    consistency_array, test = check_consistency_4(
        vegetation_type, building_type, pavement_type, water_type
    )

    # Check for consistency and fill empty fields with default vegetation type
    consistency_array, test = check_consistency_3(vegetation_type, pavement_type, soil_type)

    surface_fraction = domain.surface_fraction.from_file()
    surface_fraction[0, :, :] = ma.where(vegetation_type.mask, 0.0, 1.0)
    surface_fraction[1, :, :] = ma.where(pavement_type.mask, 0.0, 1.0)
    surface_fraction[2, :, :] = ma.where(water_type.mask, 0.0, 1.0)
    domain.surface_fraction.to_file(surface_fraction)
