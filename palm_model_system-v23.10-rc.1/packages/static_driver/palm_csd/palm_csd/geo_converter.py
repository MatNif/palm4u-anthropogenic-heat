from typing import Any, ClassVar, Dict, List, Optional, Tuple

import affine
import numpy as np
import numpy.typing as npt
import pyproj
import rasterio
import rasterio.crs
import rasterio.transform
import rasterio.warp

from palm_csd.csd_config import CSDConfigDomain, CSDConfigSettings


class GeoConverter:
    """Class to convert between different coordinate systems."""

    dst_crs: rasterio.CRS
    dst_transform: rasterio.Affine
    crs_wgs84: ClassVar[rasterio.CRS] = rasterio.crs.CRS.from_epsg(4326)

    pixel_size: float
    rotation_angle: float

    dst_width: int
    dst_height: int

    lower_left_x: float
    lower_left_y: float

    origin_x: float
    origin_y: float
    origin_lon: float
    origin_lat: float

    def __init__(
        self,
        domain_config: CSDConfigDomain,
        settings: CSDConfigSettings,
        parent: Optional["GeoConverter"] = None,
        root_parent: Optional["GeoConverter"] = None,
    ):
        if settings.epsg is None:
            raise ValueError("No EPSG code given.")
        # destination crs
        self.dst_crs = rasterio.crs.CRS.from_epsg(settings.epsg)

        self.pixel_size = domain_config.pixel_size
        self.rotation_angle = settings.rotation_angle

        # domain width and height
        self.dst_width = domain_config.nx + 1
        self.dst_height = domain_config.ny + 1

        if parent is None:
            self.lower_left_x = 0.0
            self.lower_left_y = 0.0
            if domain_config.origin_x is not None and domain_config.origin_y is not None:
                self.origin_x = domain_config.origin_x
                self.origin_y = domain_config.origin_y
        else:
            if root_parent is None:
                raise ValueError("Root parent not given.")

            # Check compatibility of parent and child
            if parent.dst_crs != self.dst_crs:
                raise ValueError("Parent and child have different CRS.")

            if parent.rotation_angle != self.rotation_angle:
                raise ValueError("Parent and child have different rotation angle.")

            if not parent.pixel_size % self.pixel_size == 0.0:
                raise ValueError(
                    "Parent pixel size is not an integer multiple of child pixel size."
                )

            if not (self.pixel_size * self.dst_width) % parent.pixel_size == 0.0:
                raise ValueError("Child domain x size not compatible with parent pixel size.")

            if not (self.pixel_size * self.dst_height) % parent.pixel_size == 0.0:
                raise ValueError("Child domain y size not compatible with parent pixel size.")

            # Child coordinates relative to parent
            if domain_config.lower_left_x is not None and domain_config.lower_left_y is not None:
                self.lower_left_x = domain_config.lower_left_x
                self.lower_left_y = domain_config.lower_left_y

                diff_lower_left_x = self.lower_left_x - parent.lower_left_x
                diff_lower_left_y = self.lower_left_y - parent.lower_left_y
                if diff_lower_left_x % parent.pixel_size != 0.0:
                    raise ValueError(
                        "Child x position not compatible with parent pixel size and position."
                    )
                if diff_lower_left_y % parent.pixel_size != 0.0:
                    raise ValueError(
                        "Child y position not compatible with parent pixel size and position."
                    )

            else:
                # Calculate lower_left_x/y from origin_*

                # Preliminay origin_x/y
                if domain_config.origin_x is not None and domain_config.origin_y is not None:
                    origin_x_prelim = domain_config.origin_x
                    origin_y_prelim = domain_config.origin_y
                elif domain_config.origin_lon is not None and domain_config.origin_lat is not None:
                    tmp_x, tmp_y = self.transform_points_from_wgs84(
                        [domain_config.origin_lon], [domain_config.origin_lat]
                    )
                    origin_x_prelim = tmp_x[0]
                    origin_y_prelim = tmp_y[0]
                else:
                    raise ValueError("Not all required input needed for geo conversion given.")

                # Preliminary lower left corner of child domain, possibly not compatible
                # with its parent
                # origin_x/y -> lower_left_x/y: rotate around root_parent.origin_x/y
                #                               with -rotation_angle
                lower_left_x_prelim, lower_left_y_prelim = self._rotate(
                    origin_x_prelim - root_parent.origin_x,
                    origin_y_prelim - root_parent.origin_y,
                    -self.rotation_angle,
                )

                # Calculate lower left corner of child domain compatible with its parent
                self.lower_left_x = (
                    np.round((lower_left_x_prelim - parent.lower_left_x) / parent.pixel_size)
                    * parent.pixel_size
                    + parent.lower_left_x
                )
                self.lower_left_y = (
                    np.round((lower_left_y_prelim - parent.lower_left_y) / parent.pixel_size)
                    * parent.pixel_size
                    + parent.lower_left_y
                )

                print(" Position relative the root parent domain:")
                print(f" lower_left_x: {self.lower_left_x}")
                print(f" lower_left_y: {self.lower_left_y}")

                if (
                    self.lower_left_x != lower_left_x_prelim
                    or self.lower_left_y != lower_left_y_prelim
                ):
                    print(" Adjusted origin_x/y or origin_lon/lat to be compatible with parent.")

            # Calculate origin_x/y from lower_left_x/y
            # lower_left_x/y -> origin_x/y: rotate with rotation_angle and
            #                               add root_parent.origin_x/y
            lower_left_x_rot, lower_left_y_rot = self._rotate(
                self.lower_left_x, self.lower_left_y, self.rotation_angle
            )
            self.origin_x = root_parent.origin_x + lower_left_x_rot
            self.origin_y = root_parent.origin_y + lower_left_y_rot

        # origin_x/y set above? not done for root domain
        if hasattr(self, "origin_x") and hasattr(self, "origin_y"):
            tmp_lon, tmp_lat = self.transform_points_to_wgs84([self.origin_x], [self.origin_y])
            self.origin_lon = tmp_lon[0]
            self.origin_lat = tmp_lat[0]
        # need to calculate origin_x/y from origin_lon/lat
        elif domain_config.origin_lon is not None and domain_config.origin_lat is not None:
            self.origin_lon = domain_config.origin_lon
            self.origin_lat = domain_config.origin_lat

            tmp_x, tmp_y = self.transform_points_from_wgs84(
                [domain_config.origin_lon], [domain_config.origin_lat]
            )
            self.origin_x = tmp_x[0]
            self.origin_y = tmp_y[0]
        else:
            raise ValueError("Not all required input needed for geo conversion given.")

        # top left corner of the domain is needed for the transform
        top_left_x = self.origin_x
        top_left_y = self.origin_y + self.dst_height * self.pixel_size

        # unrotated transform
        self.dst_transform = rasterio.transform.from_origin(
            top_left_x, top_left_y, self.pixel_size, self.pixel_size
        )

        # rotation with rotation_angle (clockwise) around
        # (x=0, y=self.dst_height) (relative to top left corner)
        rotation = affine.Affine.rotation(self.rotation_angle, (0, self.dst_height))
        self.dst_transform = self.dst_transform * rotation

    @staticmethod
    def _transform_points(
        src_crs: rasterio.CRS,
        dst_crs: rasterio.CRS,
        x: npt.ArrayLike,
        y: npt.ArrayLike,
    ) -> Tuple[List[float], List[float]]:
        """Transform points with coordinates `x` and `y` from `src_crs` to `dst_crs`."""
        x_transform, y_transform, *_ = rasterio.warp.transform(src_crs, dst_crs, x, y)
        return x_transform, y_transform

    @staticmethod
    def _rotate(x: float, y: float, angle: float) -> Tuple[float, float]:
        """Rotate point `(x, y)` by `angle` (anti-clockwise)."""
        return affine.Affine.rotation(angle) * (x, y)

    def global_palm_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
        x_global = np.arange(0.5, self.dst_width + 0.5) * self.pixel_size + self.lower_left_x
        y_global = np.arange(0.5, self.dst_height + 0.5) * self.pixel_size + self.lower_left_y
        return x_global, y_global

    def geographic_coordinates(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Calculate x and y coordinates of the domain as defined by `self.dst_transform
        and transform them from `self.dst_crs` to WGS84."""

        # create meshgrid of coordinates, flip rows to get correct order
        cols, rows = np.meshgrid(np.arange(self.dst_width), np.arange(self.dst_height))
        rows = np.flipud(rows)

        # x and y coordinates in dst_crs
        xs, ys = rasterio.transform.xy(self.dst_transform, rows, cols)
        x_coord = np.array(xs)
        y_coord = np.array(ys)

        # transform to wgs84
        lon, lat = self.transform_points_to_wgs84(x_coord.flatten(), y_coord.flatten())
        lon_coord = np.array(lon).reshape(self.dst_height, self.dst_width)
        lat_coord = np.array(lat).reshape(self.dst_height, self.dst_width)

        return x_coord, y_coord, lon_coord, lat_coord

    def transform_points_from_wgs84(
        self, x: npt.ArrayLike, y: npt.ArrayLike
    ) -> Tuple[List[float], List[float]]:
        """Transform points with coordinates `x` and `y` from WGS84 to `self.dst_crs`."""
        x_transform, y_transform, *_ = self._transform_points(
            self.crs_wgs84,
            self.dst_crs,
            x,
            y,
        )
        return x_transform, y_transform

    def transform_points_to_wgs84(
        self, x: npt.ArrayLike, y: npt.ArrayLike
    ) -> Tuple[List[float], List[float]]:
        """Transform points with coordinates `x` and `y` from `self.dst_crs` to WGS84."""
        x_transform, y_transform, *_ = self._transform_points(self.dst_crs, self.crs_wgs84, x, y)
        return x_transform, y_transform

    def dst_crs_to_cf(self) -> Dict[str, Any]:
        crs = pyproj.CRS(self.dst_crs).to_cf()
        # PALM wants the epsg_code and units, add it if not already present
        if "epsg_code" not in crs:
            crs["epsg_code"] = f"EPSG:{self.dst_crs.to_epsg()}"
        if "units" not in crs:
            unit = self.dst_crs.linear_units
            if unit in ["metre", "meter", "m"]:
                crs["units"] = "m"
            else:
                raise ValueError(f"Unknown linear unit {unit}")
        return crs
