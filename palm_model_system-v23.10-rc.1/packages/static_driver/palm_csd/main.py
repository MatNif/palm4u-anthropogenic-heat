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
# Create static-driver files for the PALM model system from rastered netCDF input.
#
# @Author Tobias Gronemeier (gronemeier@muk.uni-hannover.de)
# ------------------------------------------------------------------------------ #

import sys
from palm_csd.create_driver import create_driver


def main():
    '''Start the main routine.'''

    # Get name of configuration file
    input_configuration_file = 'csd_default.config'
    for i in range(1, len(sys.argv)):
        input_configuration_file = str(sys.argv[i])

    # Start main program
    create_driver(input_configuration_file)


if __name__ == '__main__':

    main()
