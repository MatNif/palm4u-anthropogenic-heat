"""Run the Berlin test cases"""
import itertools
import os
import shutil
from typing import Dict, Generator, List, Optional, Tuple

import pytest
import rasterio
import rasterio.crs
import yaml
from netCDF4 import Dataset  # type: ignore

from palm_csd.create_driver import create_driver
from palm_csd.csd_config import _reset_all_config_counters
from palm_csd.geo_converter import GeoConverter
from palm_csd.netcdf_data import remove_existing_file
from tests.tools import ncdf_equal

test_folder = "tests/99_full_application/"
test_folder_ref = test_folder + "output/"

epsg = 25833
# coordinates of the root domain
origin_x_root = 386891.5
origin_y_root = 5818569.0
origin_lon_root = 13.333505900572572
origin_lat_root = 52.50549956009953
# coordinates of the child domain
# aligned coordinates fit to the root domain
# invalid coordinates are slightly shifted and do not align
# when using the GeoConverter, the invalid coordinates should be corrected to the aligned ones
origin_x_nest_aligned = 389456.5
origin_y_nest_aligned = 5819499.0
origin_lon_nest_aligned = 13.3709716985343
origin_lat_nest_aligned = 52.5143830934953
origin_x_nest_invalid = 389462.5
origin_y_nest_invalid = 5819496.0
origin_lon_nest_invalid = 13.371061076093232
origin_lat_nest_invalid = 52.51435735065383
lower_left_x_nest = 2565.0
lower_left_y_nest = 930.0


@pytest.fixture(autouse=True)
def config_counters():
    """Reset config class counters to allow independent test runs."""
    yield
    _reset_all_config_counters()


def modify_configuration(
    config_in: str,
    config_out: str,
    to_delete: Optional[List[List[str]]],
    to_set: Optional[List[List]],
) -> None:
    """Generate a configuration YAML with `to_delete` deleted and `to_set` set in the
    original `config_in`. The new configuration is saved in `config_out`.
    """

    with open(config_in, "r", encoding="utf-8") as file:
        complete_dict = yaml.safe_load(file)

    def nested_del(dic: Dict, keys: List):
        for key in keys[:-1]:
            dic = dic[key]
        del dic[keys[-1]]

    def nested_set(dic: Dict, keys: List, value):
        for key in keys[:-1]:
            dic = dic.setdefault(key, {})
        dic[keys[-1]] = value

    if to_delete is not None:
        for td in to_delete:
            nested_del(complete_dict, td)
    if to_set is not None:
        for ts in to_set:
            nested_set(complete_dict, ts[0], ts[1])

    with open(config_out, "w", encoding="utf-8") as file:
        yaml.dump(complete_dict, file)


@pytest.fixture
def configuration_wrong_range_tree_trunk_diameter() -> Generator[Tuple[str, str, str], None, None]:
    """Generate a configuration file with the tree trunk diameter wrong range input file."""

    config_in = test_folder + "berlin_tiergarten.yml"
    config_out = test_folder + "berlin_tiergarten_wrong_range_tree_trunk_diameter.yml"
    file_out = "berlin_tiergarten_wrong_range_tree_trunk_diameter"
    file_ref = test_folder_ref + "berlin_tiergarten"

    to_set = [
        [
            ["input_01", "file_tree_trunk_diameter"],
            "Berlin_trees_trunk_clean_15m_DLR_wrong_range.nc",
        ],
        [
            ["input_02", "file_tree_trunk_diameter"],
            "Berlin_trees_trunk_clean_3m_DLR_wrong_range.nc",
        ],
        [["output", "file_out"], file_out],
    ]

    modify_configuration(config_in, config_out, None, to_set)
    yield config_out, test_folder + file_out, file_ref
    os.remove(config_out)


@pytest.fixture
def configuration_no_coordinates_input(
    request,
) -> Generator[Tuple[str, str, str], None, None]:
    """Generate a configuration file with defined origin_x, origin_y, epsg_code
    but not coordinate files.
    """

    run = request.param

    config_in = test_folder + "berlin_tiergarten.yml"
    config_out = test_folder + f"berlin_tiergarten_no_coordinates_{run}.yml"
    file_out = f"berlin_tiergarten_no_coordinates_{run}"
    file_ref = test_folder_ref + "berlin_tiergarten"

    # remove coordinate inputs from config
    to_delete = []
    for file_input in ["input_01", "input_02"]:
        for file in ["file_x_UTM", "file_y_UTM", "file_lon", "file_lat"]:
            to_delete.append([file_input, file])

    # set new values depending on the run parameter
    to_set = [[["output", "file_out"], file_out]]
    # to_set.extend([[["settings", "rotation_angle"], rotation_angle]])
    if run == "full":
        to_set.extend(
            [
                [["settings", "epsg"], epsg],
                [["domain_root", "origin_x"], origin_x_root],
                [["domain_root", "origin_y"], origin_y_root],
                [["domain_root", "origin_lon"], origin_lon_root],
                [["domain_root", "origin_lat"], origin_lat_root],
                [["domain_N02", "origin_x"], origin_x_nest_aligned],
                [["domain_N02", "origin_y"], origin_y_nest_aligned],
                [["domain_N02", "origin_lon"], origin_lon_nest_aligned],
                [["domain_N02", "origin_lat"], origin_lat_nest_aligned],
            ]
        )
    elif run == "origin_xy":
        to_set.extend(
            [
                [["settings", "epsg"], epsg],
                [["domain_root", "origin_x"], origin_x_root],
                [["domain_root", "origin_y"], origin_y_root],
                [["domain_N02", "origin_x"], origin_x_nest_aligned],
                [["domain_N02", "origin_y"], origin_y_nest_aligned],
            ]
        )
    elif run == "origin_lonlat":
        to_set.extend(
            [
                [["settings", "epsg"], epsg],
                [["domain_root", "origin_lon"], origin_lon_root],
                [["domain_root", "origin_lat"], origin_lat_root],
                [["domain_N02", "origin_lon"], origin_lon_nest_aligned],
                [["domain_N02", "origin_lat"], origin_lat_nest_aligned],
            ]
        )
    elif run == "no_epsg":
        to_set.extend(
            [
                [["domain_root", "origin_x"], origin_x_root],
                [["domain_root", "origin_y"], origin_y_root],
                [["domain_N02", "origin_x"], origin_x_nest_aligned],
                [["domain_N02", "origin_y"], origin_y_nest_aligned],
            ]
        )
    elif run == "no_origin":
        to_set.extend(
            [
                [["settings", "epsg"], epsg],
                [["domain_root", "origin_x"], origin_x_root],
                [["domain_root", "origin_lat"], origin_lat_root],
                [["domain_N02", "origin_y"], origin_y_nest_aligned],
                [["domain_N02", "origin_lon"], origin_lon_nest_aligned],
            ]
        )
    else:
        raise ValueError("Unknown parameter")

    modify_configuration(config_in, config_out, to_delete, to_set)

    yield config_out, test_folder + file_out, file_ref

    # remove the created file
    os.remove(config_out)


@pytest.fixture
def configuration_rotation(
    request,
) -> Generator[Tuple[str, str, str], None, None]:
    """Generate a configuration file with rotation and a result file from the default
    and the diff.
    """

    run = request.param[0]
    rotation_angle = request.param[1]

    config_in = test_folder + "berlin_tiergarten.yml"
    config_out = test_folder + f"berlin_tiergarten_no_coordinates_{rotation_angle}_{run}.yml"
    file_out = f"berlin_tiergarten_no_coordinates_{rotation_angle}_{run}"
    # reference file with be generated from the `orig` non-rotated case and a diff
    file_ref_orig = test_folder_ref + "berlin_tiergarten"
    file_ref_diff = test_folder_ref + f"diff_berlin_tiergarten_no_coordinates_{rotation_angle}"
    file_ref = test_folder_ref + f"berlin_tiergarten_no_coordinates_{rotation_angle}_{run}"

    # need to move the nested domain to be included in the parent domain
    # with the following, the nested domain is at the same relative position within the
    # parent domain as in the non-rotated case
    # all results except coordinates should be equal to non-rotated case
    diff_origin_x_invalid_rot, diff_origin_y_invalid_rot = GeoConverter._rotate(
        origin_x_nest_invalid - origin_x_root, origin_y_nest_invalid - origin_y_root, rotation_angle
    )
    origin_x_nest_rot_invalid = origin_x_root + diff_origin_x_invalid_rot
    origin_y_nest_rot_invalid = origin_y_root + diff_origin_y_invalid_rot

    # repeat for corrected coordinates
    diff_origin_x_aligned_rot, diff_origin_y_aligned_rot = GeoConverter._rotate(
        origin_x_nest_aligned - origin_x_root, origin_y_nest_aligned - origin_y_root, rotation_angle
    )
    origin_x_nest_rot_aligned = origin_x_root + diff_origin_x_aligned_rot
    origin_y_nest_rot_aligned = origin_y_root + diff_origin_y_aligned_rot

    lon, lat = GeoConverter._transform_points(
        rasterio.crs.CRS.from_epsg(epsg),
        GeoConverter.crs_wgs84,
        [origin_x_nest_rot_aligned],
        [origin_y_nest_rot_aligned],
    )
    origin_lon_nest_rot_aligned = lon[0]
    origin_lat_nest_rot_aligned = lat[0]

    # remove coordinate inputs from config
    to_delete = []
    for file_input in ["input_01", "input_02"]:
        for file in ["file_x_UTM", "file_y_UTM", "file_lon", "file_lat"]:
            to_delete.append([file_input, file])

    # set new values, use invalid origin_?_nest here to test the correction
    to_set = [[["output", "file_out"], file_out]]
    to_set.extend(
        [
            [["settings", "rotation_angle"], rotation_angle],
            [["settings", "epsg"], epsg],
            [["domain_root", "origin_x"], origin_x_root],
            [["domain_root", "origin_y"], origin_y_root],
        ]
    )
    if run == "lower_left":
        to_set.extend(
            [
                [["domain_N02", "lower_left_x"], lower_left_x_nest],
                [["domain_N02", "lower_left_y"], lower_left_y_nest],
            ]
        )
    elif run == "origin_xy":
        to_set.extend(
            [
                [["domain_N02", "origin_x"], origin_x_nest_rot_invalid],
                [["domain_N02", "origin_y"], origin_y_nest_rot_invalid],
            ]
        )
    else:
        raise ValueError("Unknown parameter")

    modify_configuration(config_in, config_out, to_delete, to_set)

    # generate result file from non-rotated case and diff file
    # these diff files are generated after the runs with the following command:
    # for i in berlin_tiergarten_no_coordinates_*; do
    #     ncks -O -L9 -v E_UTM,N_UTM,lat,lon $i output/diff_$i
    # done
    for nest in ["_root", "_N02"]:
        shutil.copyfile(file_ref_orig + nest, file_ref + nest)
        ds = Dataset(file_ref + nest, "a")
        ds_diff = Dataset(file_ref_diff + nest, "r")
        for variable in ["E_UTM", "N_UTM", "lat", "lon"]:
            ds[variable][:] = ds_diff[variable][:]
        ds.rotation_angle = ds_diff.rotation_angle
        if nest == "_N02":
            ds.origin_x = origin_x_nest_rot_aligned
            ds.origin_y = origin_y_nest_rot_aligned
            ds.origin_lon = origin_lon_nest_rot_aligned
            ds.origin_lat = origin_lat_nest_rot_aligned
        ds.close()
        ds_diff.close()

    yield config_out, test_folder + file_out, file_ref

    # remove configuation and result files
    os.remove(config_out)
    for nest in ["_root", "_N02"]:
        os.remove(file_ref + nest)


def test_complete_run():
    """Run the Berlin test case and compare with correct output"""

    create_driver(test_folder + "berlin_tiergarten.yml")

    assert ncdf_equal(
        test_folder_ref + "berlin_tiergarten_root",
        test_folder + "berlin_tiergarten_root",
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        test_folder_ref + "berlin_tiergarten_N02",
        test_folder + "/berlin_tiergarten_N02",
    ), "Nest driver does not comply with reference"

    os.remove(test_folder + "berlin_tiergarten_root")
    os.remove(test_folder + "berlin_tiergarten_N02")


@pytest.mark.parametrize(
    "configuration_no_coordinates_input",
    ["full", "origin_xy", "origin_lonlat"],
    indirect=True,
)
def test_no_coordinates_successful(configuration_no_coordinates_input):
    """Run the Berlin test case without coordinate input but with
    enough information to calculate it."""

    create_driver(configuration_no_coordinates_input[0])

    output_root = configuration_no_coordinates_input[1] + "_root"
    output_nest = configuration_no_coordinates_input[1] + "_N02"

    output_root_ref = configuration_no_coordinates_input[2] + "_root"
    output_nest_ref = configuration_no_coordinates_input[2] + "_N02"

    # exclude crs from comparison because it is generate with pyproj and not read from file
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Nest driver does not comply with reference"

    # check crs manually
    with Dataset(output_root) as nc_data:
        crs_root = nc_data.variables["crs"].__dict__
    with Dataset(output_nest) as nc_data:
        crs_nest = nc_data.variables["crs"].__dict__

    crs_ref = {
        "long_name": "coordinate reference system",
        "crs_wkt": 'PROJCRS["ETRS89 / UTM zone 33N",BASEGEOGCRS["ETRS89",'
        'DATUM["European Terrestrial Reference System 1989",'
        'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
        'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4258]],'
        'CONVERSION["UTM zone 33N",METHOD["Transverse Mercator",ID["EPSG",9807]],'
        'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],'
        'ID["EPSG",8801]],PARAMETER["Longitude of natural origin",15,'
        'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
        'PARAMETER["Scale factor at natural origin",0.9996,SCALEUNIT["unity",1],'
        'ID["EPSG",8805]],PARAMETER["False easting",500000,LENGTHUNIT["metre",1],'
        'ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],'
        'ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],'
        'AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",25833]]',
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.314140356,
        "inverse_flattening": 298.257222101,
        "reference_ellipsoid_name": "GRS 1980",
        "longitude_of_prime_meridian": 0.0,
        "prime_meridian_name": "Greenwich",
        "geographic_crs_name": "ETRS89",
        "horizontal_datum_name": "European Terrestrial Reference System 1989",
        "projected_crs_name": "ETRS89 / UTM zone 33N",
        "grid_mapping_name": "transverse_mercator",
        "latitude_of_projection_origin": 0.0,
        "longitude_of_central_meridian": 15.0,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "epsg_code": "EPSG:25833",
        "units": "m",
    }

    assert crs_root == crs_ref, "Root crs does not comply with reference"
    assert crs_nest == crs_ref, "Nest crs does not comply with reference"

    os.remove(output_root)
    os.remove(output_nest)


@pytest.mark.parametrize(
    "configuration_no_coordinates_input",
    ["no_epsg", "no_origin"],
    indirect=True,
)
def test_no_coordinates_failing(configuration_no_coordinates_input):
    """Run the Berlin test case without coordinate input and
    not enough information to calculate it. This should raise an error.
    """

    with pytest.raises(ValueError):
        create_driver(configuration_no_coordinates_input[0])


run = ["lower_left", "origin_xy"]
angles = [30, 165, 200, 320]
combinations = list(itertools.product(run, angles))
names = [x + " " + str(y) for x, y in combinations]


@pytest.mark.parametrize(
    "configuration_rotation",
    combinations,
    ids=names,
    indirect=True,
)
def test_rotation(configuration_rotation):
    """Run the Berlin test case with rotation."""

    create_driver(configuration_rotation[0])

    output_root = configuration_rotation[1] + "_root"
    output_nest = configuration_rotation[1] + "_N02"

    output_root_ref = configuration_rotation[2] + "_root"
    output_nest_ref = configuration_rotation[2] + "_N02"

    # exclude crs from comparison because it is generate with pyproj and not read from file
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Nest driver does not comply with reference"

    # check crs manually
    with Dataset(output_root) as nc_data:
        crs_root = nc_data.variables["crs"].__dict__
    with Dataset(output_nest) as nc_data:
        crs_nest = nc_data.variables["crs"].__dict__

    crs_ref = {
        "long_name": "coordinate reference system",
        "crs_wkt": 'PROJCRS["ETRS89 / UTM zone 33N",BASEGEOGCRS["ETRS89",'
        'DATUM["European Terrestrial Reference System 1989",'
        'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
        'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4258]],'
        'CONVERSION["UTM zone 33N",METHOD["Transverse Mercator",ID["EPSG",9807]],'
        'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],'
        'ID["EPSG",8801]],PARAMETER["Longitude of natural origin",15,'
        'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
        'PARAMETER["Scale factor at natural origin",0.9996,SCALEUNIT["unity",1],'
        'ID["EPSG",8805]],PARAMETER["False easting",500000,LENGTHUNIT["metre",1],'
        'ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],'
        'ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],'
        'AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",25833]]',
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.314140356,
        "inverse_flattening": 298.257222101,
        "reference_ellipsoid_name": "GRS 1980",
        "longitude_of_prime_meridian": 0.0,
        "prime_meridian_name": "Greenwich",
        "geographic_crs_name": "ETRS89",
        "horizontal_datum_name": "European Terrestrial Reference System 1989",
        "projected_crs_name": "ETRS89 / UTM zone 33N",
        "grid_mapping_name": "transverse_mercator",
        "latitude_of_projection_origin": 0.0,
        "longitude_of_central_meridian": 15.0,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "epsg_code": "EPSG:25833",
        "units": "m",
    }

    assert crs_root == crs_ref, "Root crs does not comply with reference"
    assert crs_nest == crs_ref, "Nest crs does not comply with reference"

    os.remove(output_root)
    os.remove(output_nest)


def test_wrong_range_tree_trunk_diameter(configuration_wrong_range_tree_trunk_diameter):
    """Run the Berlin test case with a wrong range of the tree trunk diameter.
    This should raise an error.
    """
    with pytest.raises(ValueError):
        create_driver(configuration_wrong_range_tree_trunk_diameter[0])

    remove_existing_file(configuration_wrong_range_tree_trunk_diameter[1] + "_root")
    remove_existing_file(configuration_wrong_range_tree_trunk_diameter[1] + "_N02")
