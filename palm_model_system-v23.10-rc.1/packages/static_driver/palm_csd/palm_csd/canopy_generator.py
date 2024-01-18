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
# Canopy generator routines for creating 3D leaf and basal area densities for
# single trees and tree patches
#
# @Author Bjoern Maronga (maronga@muk.uni-hannover.de)
# ------------------------------------------------------------------------------ #
import numpy as np
import numpy.ma as ma
import math


def process_patch(dz, patch_height, patch_type_2d, vegetation_type, max_height_lad, patch_lai, alpha, beta):

    phdz = patch_height[:, :] / dz
    pch_index = ma.where(patch_height.mask, int(-1), phdz.astype(int) + 1)
    ma.masked_equal(pch_index, 0, copy=False)
    pch_index = ma.where(pch_index == -1, 0, pch_index)

    max_canopy_height = max(ma.max(patch_height), max_height_lad)

    z = np.arange(0, math.floor(max_canopy_height / dz) * dz + 2 * dz, dz)

    z[1:] = z[1:] - 0.5 * dz

    nz = len(z)
    ny = len(patch_height[:, 0])
    nx = len(patch_height[0, :])

    pre_lad = ma.zeros(nz)
    lad_loc = ma.empty((nz, ny, nx))
    lad_loc.mask = True

    for i in range(0, nx):
        for j in range(0, ny):
            int_bpdf = 0.0
            if patch_height[j, i] >= (0.5 * dz):
                for k in range(1, pch_index[j, i]):
                    int_bpdf = int_bpdf + (((z[k] / patch_height[j, i])**(alpha - 1)) * (
                                (1.0 - (z[k] / patch_height[j, i]))**(beta - 1)) * (
                                                       dz / patch_height[j, i]))

                for k in range(1, pch_index[j, i]):
                    pre_lad[k] = patch_lai[j, i] * (
                                ((dz * k / patch_height[j, i])**(alpha - 1.0)) * (
                                    (1.0 - (dz * k / patch_height[j, i]))**(
                                        beta - 1.0)) / int_bpdf) / patch_height[j, i]

                lad_loc[0, j, i] = pre_lad[0]

                for k in range(0, pch_index[j, i]):
                    lad_loc[k, j, i] = 0.5 * (pre_lad[k - 1] + pre_lad[k])

    patch_id_2d = ma.where(lad_loc.mask[0, :, :], 0, 1)
    patch_id_3d = ma.where(lad_loc.mask, 0, 1)
    patch_type_3d = ma.empty((nz, ny, nx))

    for k in range(0, nz):
        patch_id_3d[k, :, :] = ma.where((patch_id_2d != 0) & ~lad_loc.mask[k, :, :],
                                     patch_id_2d, ma.masked)
        patch_type_3d[k, :, :] = ma.where((patch_id_2d != 0) & ~lad_loc.mask[k, :, :],
                                       patch_type_2d, ma.masked)

    return lad_loc, patch_id_3d, patch_type_3d, nz, 0
