"""Module to handle the palm_csd configuration

Several classes CSDConfig* are defined to handle a specific section of the palm_csd configuration.
The reading of the input yaml file and storing of the several objects is done by CSDConfig.
"""

import os
from dataclasses import dataclass, field
from math import floor

# remove Dict here and replace by dict below once Python >=3.9 could be used
from typing import Dict, Optional
import yaml


class CSDConfigElement:
    """Basic class for configuration elements. It defines a counter to check how often the class
    has been used and applies certain modifications of the read in options."""

    _type: str
    counter: int = 0
    _unique: bool = True

    def __post_init__(self):
        if type(self)._unique and type(self).counter == 1:
            raise ValueError(f"More than 1 configuration section of type {type(self)._type} found")

        # increase number of processed configs
        type(self).counter += 1

        # append path value to every file entry
        if hasattr(self, "path"):
            if self.path is not None:  # type: ignore
                for key in vars(self).keys():
                    if key.startswith("file") and getattr(self, key) is not None:
                        setattr(self, key, self.path + "/" + getattr(self, key))  # type: ignore

    @classmethod
    def _reset_counter(cls):
        """Reset the counter for the class"""
        cls.counter = 0


@dataclass
class CSDConfigAttributes(CSDConfigElement):
    """Class for global attributes in the static driver"""

    # pylint: disable=too-many-instance-attributes
    # configuration elements as attributes by design

    _type = "attributes"

    author: Optional[str] = None
    contact_person: Optional[str] = None
    institution: Optional[str] = None
    acronym: Optional[str] = None
    campaign: Optional[str] = None
    location: Optional[str] = None
    site: Optional[str] = None
    comment: Optional[str] = None
    data_content: Optional[str] = None
    dependencies: Optional[str] = None
    keywords: Optional[str] = None
    references: Optional[str] = None
    source: Optional[str] = None
    palm_version: Optional[str] = None
    origin_time: Optional[str] = None

    # calculated later, domain dependent, not read from the input yaml
    origin_x: Optional[float] = field(default=None, init=False)
    origin_y: Optional[float] = field(default=None, init=False)
    origin_z: Optional[float] = field(default=None, init=False)
    origin_lat: Optional[float] = field(default=None, init=False)
    origin_lon: Optional[float] = field(default=None, init=False)


@dataclass
class CSDConfigSettings(CSDConfigElement):
    """Class for settings in the static driver"""

    # pylint: disable=too-many-instance-attributes
    # configuration elements as attributes by design

    _type = "settings"

    bridge_width: float = 3.0
    lai_roof_extensive: float = 0.8
    lai_roof_intensive: float = 2.5
    lai_tree_lower_threshold: float = 0.0
    lai_low_vegetation_default: float = 1.0
    lai_high_vegetation_default: float = 6.0
    lai_alpha: float = 5.0
    lai_beta: float = 3.0
    patch_height_default: float = 10.0
    rotation_angle: float = 0.0
    season: str = "summer"
    vegetation_type_below_trees: int = 3
    debug_mode: bool = False
    epsg: Optional[int] = None


@dataclass
class CSDConfigOutput(CSDConfigElement):
    """Class for output configuration of the static driver"""

    _type = "output"

    file_out: str
    path: Optional[str] = None
    version: Optional[int] = None


@dataclass
class CSDConfigInput(CSDConfigElement):
    """Class for input data configuration for the static driver"""

    # pylint: disable=too-many-instance-attributes
    # configuration elements as attributes by design

    _type = "input"
    _unique = False

    file_zt: str
    file_buildings_2d: str
    file_bridges_2d: str
    file_building_id: str
    file_bridges_id: str
    file_building_type: str
    file_vegetation_type: str
    file_vegetation_height: str
    file_pavement_type: str
    file_water_type: str
    file_street_type: str
    file_street_crossings: str
    file_vegetation_on_roofs: str
    file_tree_type: str
    file_patch_height: str
    pixel_size: float
    file_x_UTM: Optional[str] = None
    file_y_UTM: Optional[str] = None
    file_lat: Optional[str] = None
    file_lon: Optional[str] = None
    file_lai: Optional[str] = None
    file_patch_type: Optional[str] = None
    file_soil_type: Optional[str] = None
    file_tree_crown_diameter: Optional[str] = None
    file_tree_trunk_diameter: Optional[str] = None
    file_tree_height: Optional[str] = None
    file_water_temperature: Optional[str] = None
    path: Optional[str] = None


@dataclass
class CSDConfigDomain(CSDConfigElement):
    """Class for domain configuration of the static driver"""

    # pylint: disable=too-many-instance-attributes
    # configuration elements as attributes by design

    _type = "domain"
    _unique = False

    name: str

    pixel_size: float
    input_lower_left_x: float
    input_lower_left_y: float
    nx: int
    ny: int
    dz: float
    allow_high_vegetation: bool = False
    buildings_3d: bool = False
    domain_parent: Optional[str] = None
    generate_vegetation_patches: bool = True
    interpolate_terrain: bool = False
    lower_left_x: Optional[float] = None
    lower_left_y: Optional[float] = None
    origin_lon: Optional[float] = None
    origin_lat: Optional[float] = None
    origin_x: Optional[float] = None
    origin_y: Optional[float] = None
    overhanging_trees: bool = True
    remove_low_lai_tree: bool = False
    street_trees: bool = True
    use_palm_z_axis: bool = False
    vegetation_on_roofs: bool = True
    water_temperature_per_water_type: Optional[Dict[int, float]] = None

    # output filename generated from output section and domain name later and not read from
    # input yaml file
    filename: str = field(init=False)

    def __post_init__(self):
        super().__post_init__()

        self.x0 = int(floor(self.input_lower_left_x / self.pixel_size))
        self.y0 = int(floor(self.input_lower_left_y / self.pixel_size))
        self.x1 = self.x0 + self.nx
        self.y1 = self.y0 + self.ny

        if self.interpolate_terrain and not self.use_palm_z_axis:
            self.use_palm_z_axis = True
            print("+++ Overwrite user setting for use_palm_z_axis")


class CSDConfig:
    """Class to collect all palm_csd configuration sections"""

    attributes: CSDConfigAttributes
    settings: CSDConfigSettings
    output: CSDConfigOutput
    input_dict: Dict[str, CSDConfigInput]
    domain_dict: Dict[str, CSDConfigDomain]

    def __init__(self, configuration_file: str):
        # Check if configuration files exists and quit otherwise
        if not os.path.isfile(configuration_file):
            print("Error. No configuration file " + configuration_file + " found.")
            raise FileNotFoundError
        with open(configuration_file, "r", encoding="utf-8") as file:
            complete_dict = yaml.safe_load(file)

        self.domain_dict = {}
        self.input_dict = {}
        for key in complete_dict.keys():
            key_splitted = key.split("_")
            if len(key_splitted) > 2:
                raise ValueError(key + "includes too many separators")
            # use match/case for the following once Python >=3.10 could be used
            if key_splitted[0] == "attributes":
                self.attributes = CSDConfigAttributes(**complete_dict[key])
            elif key_splitted[0] == "settings":
                self.settings = CSDConfigSettings(**complete_dict[key])
            elif key_splitted[0] == "output":
                self.output = CSDConfigOutput(**complete_dict[key])
            elif key_splitted[0] == "input":
                self.input_dict[key_splitted[1]] = CSDConfigInput(**complete_dict[key])
            elif key_splitted[0] == "domain":
                self.domain_dict[key_splitted[1]] = CSDConfigDomain(
                    **complete_dict[key], name=key_splitted[1]
                )
            else:
                raise ValueError("Unknown configuration section " + key_splitted[0])

        # set output file names for each domain
        for domain_config in self.domain_dict.values():
            domain_config.filename = self.output.file_out + "_" + domain_config.name

    def input_of_domain(self, domain_config: CSDConfigDomain) -> CSDConfigInput:
        """Return the fitting input configuration for a domain according to the pixel_size"""
        for input_config in self.input_dict.values():
            if input_config.pixel_size == domain_config.pixel_size:
                return input_config
        raise ValueError("No fitting input configuration section for domain " + domain_config.name)

    def input_of_parent_domain(self, domain_config: CSDConfigDomain) -> Optional[CSDConfigInput]:
        """Return the input configuration for the parent of a domain according to the pixel_size"""
        domain_parent = domain_config.domain_parent
        if domain_parent is None:
            return None
        return self.input_of_domain(self.domain_dict[domain_parent])

    def domain_of_parent_domain(self, domain_config: CSDConfigDomain) -> Optional[CSDConfigDomain]:
        """Return the domain configuration for the parent of a domain according to the pixel_size"""
        domain_parent = domain_config.domain_parent
        if domain_parent is None:
            return None
        return self.domain_dict[domain_parent]


def _reset_all_config_counters() -> None:
    """Reset counters of all config classes"""
    for cls in [
        CSDConfigAttributes,
        CSDConfigSettings,
        CSDConfigOutput,
        CSDConfigInput,
        CSDConfigDomain,
    ]:
        cls._reset_counter()
