from dataclasses import dataclass
from importlib.resources import open_text
from math import pi
from typing import ClassVar, List, Union

import numpy as np
import numpy.ma as ma
import numpy.ma.core as ma_core
import scipy.integrate as integrate

from palm_csd.csd_config import CSDConfigDomain, CSDConfigSettings

# TODO Python 3.11 supports Self in typing
# from typing import Self


@dataclass
class Tree:
    """A tree object with size parameters."""

    # defines the general shape of the tree and can be one of the following types:
    # 1 sphere or ellipsoid
    # 2 cylinder
    # 3 cone
    # 4 inverted cone
    # 5 paraboloid (rounded cone)
    # 6 inverted paraboloid (inverted rounded cone)
    shape: int
    # ratio of maximum crown height to maximum crown diameter
    crown_ratio: float
    # crown diameter (m)
    crown_diameter: float
    # total height of the tree including trunk (m)
    height: float
    # height where the leaf area density is maximum relative to total tree height
    lad_max_height: float
    # ratio of basal area in the crown area to the leaf area
    bad_scale: float
    # trunk diameter at breast height (1.4 m) (m)
    trunk_diameter: float


@dataclass
class ReferenceTree(Tree):
    """A tree object with additional reference parameters."""

    # name of the tree type
    species: str
    # default leaf area index fully leafed
    lai_summer: float
    # default winter-time leaf area index
    lai_winter: float


@dataclass
class DomainTree(Tree):
    """A tree object with all the necessary parameters to calculate the LAD and BAD fields."""

    # actual leaf area index
    lai: float
    # x coordinate of the tree
    i: int
    # y coordinate of the tree
    j: int
    # id of the tree
    id: int
    # type of the tree
    type: int

    # Table of default values
    defaults: ClassVar[List[ReferenceTree]] = []

    # Parameters of the LAD and BAD generator
    # Extinction coefficients are experimental
    sphere_extinction: ClassVar[float] = 0.6
    cone_extinction: ClassVar[float] = 0.2
    # based on Lalic and Mihailovic (2004)
    ml_n_low: ClassVar[float] = 0.5
    ml_n_high: ClassVar[float] = 6.0

    # Counters
    shallow_tree_count: ClassVar[int] = 0
    low_lai_count: ClassVar[int] = 0
    mod_count: ClassVar[int] = 0
    id_count: ClassVar[int] = 0

    @classmethod
    def generate_tree(
        cls,
        i,
        j,
        type: Union[int, ma_core.MaskedConstant],
        shape: Union[int, ma_core.MaskedConstant],
        height: Union[float, ma_core.MaskedConstant],
        lai: Union[float, ma_core.MaskedConstant],
        crown_diameter: Union[float, ma_core.MaskedConstant],
        trunk_diameter: Union[float, ma_core.MaskedConstant],
        config: CSDConfigDomain,
        settings: CSDConfigSettings,
    ):  # -> Optional[Self]:  # TODO add with Python 3.11
        """Generate a tree object. Input values are checked and default values are used
        if necessary. Depending on the set-up, trees with too low tree_height or tree_lai
        are removed.
        """

        # Increase the tree ID counter
        cls.id_count += 1

        # Check for missing data in the input and set default values if needed
        if type is ma.masked or type == -1:
            type_checked = 0
        else:
            type_checked = int(type)

        if shape is ma.masked:
            shape_checked = cls.defaults[type_checked].shape
        else:
            shape_checked = int(shape)
        if shape_checked > 6 or shape_checked < 0:
            raise ValueError(f"Tree shape must be between 0 and 6 instead of {shape_checked}.")

        if height is ma.masked or height <= 0.0:
            height_checked = cls.defaults[type_checked].height
        else:
            height_checked = float(height)

        if lai is ma.masked:
            if settings.season == "summer":
                lai_checked = cls.defaults[type_checked].lai_summer
            elif settings.season == "winter":
                lai_checked = cls.defaults[type_checked].lai_winter
            else:
                raise ValueError(
                   f"Season must either be 'summer' or 'winter' instead of {settings.season}."
                )
        else:
            lai_checked = float(lai)

        if crown_diameter is ma.masked or crown_diameter <= 0.0:
            crown_diameter_checked = cls.defaults[type_checked].crown_diameter
        else:
            crown_diameter_checked = float(crown_diameter)

        if trunk_diameter is ma.masked or trunk_diameter <= 0.0:
            trunk_diameter_checked = cls.defaults[type_checked].trunk_diameter
        else:
            trunk_diameter_checked = float(trunk_diameter)

        # Very small trees are ignored
        if height_checked <= (0.5 * config.dz):
            cls.shallow_tree_count += 1
            print("    Shallow tree found. Action: ignore.")
            return None

        # Check tree_lai
        # Tree LAI lower than threshold?
        if lai_checked < settings.lai_tree_lower_threshold:
            # Deal with low lai tree
            cls.mod_count += 1
            if config.remove_low_lai_tree:
                # Skip this tree
                print(f"Removed tree with LAI = {lai_checked:0.3f} at ({i}, {j}).")
                return None
            else:
                # Use type specific default
                if settings.season == "summer":
                    lai_checked = cls.defaults[type_checked].lai_summer
                elif settings.season == "winter":
                    lai_checked = cls.defaults[type_checked].lai_winter
                else:
                    raise ValueError(
                        f"Season must either be 'summer' or 'winter' instead of {settings.season}."
                    )
                print(f"Adjusted tree to LAI = {lai_checked:0.3f} at ({i}, {j}).")

        # Warn about a tree with lower LAI than we would expect in winter
        if lai_checked < cls.defaults[type_checked].lai_winter:
            cls.low_lai_count += 1
            print(
                f"Found tree with LAI = {lai_checked:0.3f} (tree-type specific default winter LAI "
                + f"of {cls.defaults[type_checked].lai_winter:0.2}) at ({i}, {j})."
            )

        # Assign values that are not defined as user input from lookup table
        crown_ratio_checked = cls.defaults[type_checked].crown_ratio
        lad_max_height_checked = cls.defaults[type_checked].lad_max_height
        bad_scale_checked = cls.defaults[type_checked].bad_scale

        return cls(
            type=type_checked,
            shape=shape_checked,
            crown_ratio=crown_ratio_checked,
            crown_diameter=crown_diameter_checked,
            height=height_checked,
            lad_max_height=lad_max_height_checked,
            bad_scale=bad_scale_checked,
            trunk_diameter=trunk_diameter_checked,
            lai=lai_checked,
            i=i,
            j=j,
            id=cls.id_count,
        )

    @classmethod
    def reset_counter(cls) -> None:
        """Reset the tree counters."""
        cls.shallow_tree_count = 0
        cls.low_lai_count = 0
        cls.mod_count = 0
        cls.id_count = 0

    @classmethod
    def check_counter(cls, config: CSDConfigDomain) -> None:
        """Print a summary of the tree and tree LAI adjustments."""
        if cls.shallow_tree_count > 0:
            print(f"Removed {cls.shallow_tree_count} shallow trees with height < 1/2 dz.")
        if cls.mod_count > 0:
            if config.remove_low_lai_tree:
                print(f"Removed {cls.mod_count} trees due to low LAI.")
            else:
                print(f"Adjusted LAI of {cls.mod_count} trees.")
        if cls.low_lai_count > 0:
            print(
                f"Warning: Found {cls.low_lai_count} trees with LAI lower then the "
                + "tree-type specific default winter LAI. "
                + "Consider adjusting lai_tree_lower_threshold and remove_low_lai_tree.",
            )

    def generate_store_3d_fields(
        self,
        lad_global: ma.MaskedArray,
        bad_global: ma.MaskedArray,
        id_global: ma.MaskedArray,
        type_global: ma.MaskedArray,
        config: CSDConfigDomain,
    ) -> None:
        """Generate the LAD and BAD profile. Store them in the global arrays given as input.
        Also store tree id and type with the same extent as the LAD and BAD arrays in the
        respective input.
        """

        # Calculate crown height and height of the crown center
        crown_height = self.crown_ratio * self.crown_diameter
        if crown_height > self.height:
            crown_height = self.height

        crown_center = self.height - crown_height * 0.5

        # Calculate height of maximum LAD
        z_lad_max = self.lad_max_height * self.height

        # Calculate the maximum LAD after Lalic and Mihailovic (2004)
        lad_max_part_1 = integrate.quad(
            lambda z: ((self.height - z_lad_max) / (self.height - z)) ** self.ml_n_high
            * np.exp(self.ml_n_high * (1.0 - (self.height - z_lad_max) / (self.height - z))),
            0.0,
            z_lad_max,
        )
        lad_max_part_2 = integrate.quad(
            lambda z: ((self.height - z_lad_max) / (self.height - z)) ** self.ml_n_low
            * np.exp(self.ml_n_low * (1.0 - (self.height - z_lad_max) / (self.height - z))),
            z_lad_max,
            self.height,
        )

        lad_max = self.lai / (lad_max_part_1[0] + lad_max_part_2[0])

        # Define position of tree and its output domain
        nx = int(self.crown_diameter / config.pixel_size) + 2
        nz = int(self.height / config.dz) + 2

        # Add one grid point if diameter is an odd value
        if (self.crown_diameter % 2.0) != 0.0:
            nx = nx + 1

        # Create local domain of the tree's LAD
        x = np.arange(0, nx * config.pixel_size, config.pixel_size)
        x[:] = x[:] - 0.5 * config.pixel_size
        y = x

        z = np.arange(0, nz * config.dz, config.dz)
        z[1:] = z[1:] - 0.5 * config.dz

        # Define center of the tree position inside the local LAD domain
        location_x = x[int(nx / 2)]
        location_y = y[int(nx / 2)]

        # Calculate LAD profile after Lalic and Mihailovic (2004).
        # Will be later used for normalization
        lad_profile = np.arange(0, nz, 1.0)
        lad_profile[:] = 0.0

        for k in range(1, nz - 1):
            if (z[k] > 0.0) & (z[k] < z_lad_max):
                n = self.ml_n_high
            else:
                n = self.ml_n_low

            lad_profile[k] = (
                lad_max
                * ((self.height - z_lad_max) / (self.height - z[k])) ** n
                * np.exp(n * (1.0 - (self.height - z_lad_max) / (self.height - z[k])))
            )

        # Create lad array and populate according to the specific tree shape. 
        # NOTE This is still experimental
        lad_local = ma.masked_all((nz, nx, nx))
        bad_local = ma.copy(lad_local)

        # Branch for spheres and ellipsoids. 
        # A symmetric LAD sphere is created assuming an LAD extinction towards the center of the 
        # tree, representing the effect of sunlight extinction and thus less leaf mass inside the
        # tree crown. 
        # NOTE Extinction coefficients are experimental.
        if self.shape == 1:
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(0, nz):
                        r_test = np.sqrt(
                            (x[i] - location_x) ** 2 / (self.crown_diameter * 0.5) ** 2
                            + (y[j] - location_y) ** 2 / (self.crown_diameter * 0.5) ** 2
                            + (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** (2)
                        )
                        if r_test <= 1.0:
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.sphere_extinction * (1.0 - r_test)
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for cylinder shapes
        elif self.shape == 2:
            k_min = int((crown_center - crown_height * 0.5) / config.dz)
            k_max = int((crown_center + crown_height * 0.5) / config.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        r_test = np.sqrt(
                            (x[i] - location_x) ** 2 / (self.crown_diameter * 0.5) ** 2
                            + (y[j] - location_y) ** 2 / (self.crown_diameter * 0.5) ** 2
                        )
                        if r_test <= 1.0:
                            r_test3 = np.sqrt(
                                (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** 2
                            )
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.sphere_extinction * (1.0 - max(r_test, r_test3))
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for cone shapes
        elif self.shape == 3:
            k_min = int((crown_center - crown_height * 0.5) / config.dz)
            k_max = int((crown_center + crown_height * 0.5) / config.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k - k_min
                        r_test = (
                            (x[i] - location_x) ** 2
                            + (y[j] - location_y) ** 2
                            - ((self.crown_diameter * 0.5) ** 2 / crown_height**2)
                            * (z[k_rel] - crown_height) ** 2
                        )
                        if r_test <= 0.0:
                            r_test2 = np.sqrt(
                                (x[i] - location_x) ** 2 / (self.crown_diameter * 0.5) ** 2
                                + (y[j] - location_y) ** 2 / (self.crown_diameter * 0.5) ** 2
                            )
                            r_test3 = np.sqrt(
                                (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** 2
                            )
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.cone_extinction
                                * (1.0 - max((r_test + 1.0), r_test2, r_test3))
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for inverted cone shapes.
        # TODO: what is r_test2 and r_test3 used for? Debugging needed!
        elif self.shape == 4:
            k_min = int((crown_center - crown_height * 0.5) / config.dz)
            k_max = int((crown_center + crown_height * 0.5) / config.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k_max - k
                        r_test = (
                            (x[i] - location_x) ** 2
                            + (y[j] - location_y) ** 2
                            - ((self.crown_diameter * 0.5) ** 2 / crown_height**2)
                            * (z[k_rel] - crown_height) ** 2
                        )
                        if r_test <= 0.0:
                            r_test2 = np.sqrt(
                                (x[i] - location_x) ** 2 / (self.crown_diameter * 0.5) ** 2
                                + (y[j] - location_y) ** 2 / (self.crown_diameter * 0.5) ** 2
                            )
                            r_test3 = np.sqrt(
                                (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** 2
                            )
                            lad_local[k, j, i] = lad_max * np.exp(-self.cone_extinction * (-r_test))
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for paraboloid shapes
        elif self.shape == 5:
            k_min = int((crown_center - crown_height * 0.5) / config.dz)
            k_max = int((crown_center + crown_height * 0.5) / config.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k - k_min
                        r_test = (
                            (x[i] - location_x) ** 2 + (y[j] - location_y) ** (2)
                        ) * crown_height / (self.crown_diameter * 0.5) ** 2 - z[k_rel]
                        if r_test <= 0.0:
                            lad_local[k, j, i] = lad_max * np.exp(-self.cone_extinction * (-r_test))
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for inverted paraboloid shapes
        elif self.shape == 6:
            k_min = int((crown_center - crown_height * 0.5) / config.dz)
            k_max = int((crown_center + crown_height * 0.5) / config.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k_max - k
                        r_test = (
                            (x[i] - location_x) ** 2 + (y[j] - location_y) ** (2)
                        ) * crown_height / (self.crown_diameter * 0.5) ** 2 - z[k_rel]
                        if r_test <= 0.0:
                            lad_local[k, j, i] = lad_max * np.exp(-self.cone_extinction * (-r_test))
                        else:
                            lad_local[k, j, i] = ma.masked

        else:
            raise ValueError("Unknown tree shape.")

        # Leave if no LAD was generated
        if ma.all(lad_local.mask):
            return

        # Indicate a defined LAD in a column by setting lowest value to 0
        for i in range(0, nx):
            for j in range(0, nx):
                if ma.any(~lad_local.mask[:, j, i]):
                    lad_local[0, j, i] = 0.0

        # Normalize the LAD profile so that the vertically integrated Lalic and Mihailovic (2004) is
        # reproduced by the LAD array. Deactivated for now.
        # for i in range(0,nx):
        # for j in range(0,nx):
        # lad_clean = np.where(lad_loc[:,j,i] == fillvalues["tree_data"],0.0,lad_loc[:,j,i])
        # lai_from_int = integrate.simps (lad_clean, z)
        # print(lai_from_int)
        # for k in range(0,nz):
        # if ( np.any(lad_loc[k,j,i] > 0.0) ):
        # lad_loc[k,j,i] = np.where(
        #     (lad_loc[k,j,i] != fillvalues["tree_data"]),
        #     lad_loc[k,j,i] / lai_from_int * lad_profile[k],
        #     lad_loc[k,j,i]
        #     )

        # Create BAD array and populate. 
        # TODO: revise as low LAD inside the foliage does not result in low BAD values.
        bad_local = (1.0 - (lad_local / (ma.max(lad_local) + 0.01))) * 0.1

        # Overwrite grid cells that are occupied by the tree trunk
        radius = self.trunk_diameter * 0.5
        for i in range(0, nx):
            for j in range(0, nx):
                for k in range(0, nz):
                    if z[k] <= crown_center:
                        r_test = np.sqrt((x[i] - location_x) ** 2 + (y[j] - location_y) ** 2)
                        if r_test == 0.0:
                            if self.trunk_diameter <= config.pixel_size:
                                bad_local[k, j, i] = radius**2 * pi
                            else:
                                # WORKAROUND: divide remaining circle area over the 8 surrounding
                                # valid_pixels
                                bad_local[k, j - 1 : j + 2, i - 1 : i + 2] = radius**2 * pi / 8.0
                                # for the central pixel fill the pixel
                                bad_local[k, j, i] = config.pixel_size**2
                        # elif ( r_test <= radius ):
                        # TODO: calculate circle segment of grid points cut by the grid

        # Calculate the position of the local 3d tree array within the full
        # domain in order to achieve correct mapping and cutting off at the edges
        # of the full domain
        lad_loc_nx = int(len(x) / 2)
        lad_loc_ny = int(len(y) / 2)
        lad_loc_nz = int(len(z))

        odd_x = int(len(x) % 2)
        odd_y = int(len(y) % 2)

        ind_l_x = max(0, (self.i - lad_loc_nx))
        ind_l_y = max(0, (self.j - lad_loc_ny))
        ind_r_x = min(lad_global.shape[2] - 1, self.i + lad_loc_nx - 1 + odd_x)
        ind_r_y = min(lad_global.shape[1] - 1, self.j + lad_loc_ny - 1 + odd_y)

        out_l_x = ind_l_x - (self.i - lad_loc_nx)
        out_l_y = ind_l_y - (self.j - lad_loc_ny)
        out_r_x = len(x) - 1 + ind_r_x - (self.i + lad_loc_nx - 1 + odd_x)
        out_r_y = len(y) - 1 + ind_r_y - (self.j + lad_loc_ny - 1 + odd_y)

        lad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            ~lad_local.mask[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            lad_local[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            lad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
        bad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            ~bad_local.mask[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            bad_local[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            bad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
        id_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            ~lad_local.mask[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            self.id,
            id_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
        type_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            ~lad_local.mask[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            self.type,
            type_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )

    @classmethod
    def populate_defaults(cls) -> None:
        """Read default tree species data from file"""

        # Read csv from palm_csd.data. Use files instead of open_text for Python >=3.9
        with open_text("palm_csd.data", "tree_defaults.csv") as tree_csv:
            tree_data = np.genfromtxt(
                tree_csv,
                delimiter=",",
                dtype=None,
                names=True,
                skip_header=10,
                encoding="utf-8",
            )
        for tree in tree_data:
            cls.defaults.append(
                ReferenceTree(
                    species=tree["species"],
                    shape=tree["shape"],
                    crown_ratio=tree["crown_ratio"],
                    crown_diameter=tree["crown_diameter"],
                    height=tree["height"],
                    lai_summer=tree["lai_summer"],
                    lai_winter=tree["lai_winter"],
                    lad_max_height=tree["lad_max_height"],
                    bad_scale=tree["bad_scale"],
                    trunk_diameter=tree["trunk_diameter"],
                )
            )
