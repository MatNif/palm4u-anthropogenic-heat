#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#--------------------------------------------------------------------------------#
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
# Copyright 2022-2022  pecanode GmbH
#--------------------------------------------------------------------------------#
#
# Description:
# ------------
# Processing tool for creating PIDS conform virtual measurement setup file
# from UC2 data-standard conform observational data or from prescribed input
# coordinates.
#
# @Authors Matthias Suehring (suehring@muk.uni-hannover.de)
#          Tobias Gronemeier (gronemeier@muk.uni-hannover.de)
#          Helge Knoop (helge.knoop@pecanode.com)
#
# @todo Add further feature tpyes for customized observations. At the moment only
#       timeSeries is possible.
#--------------------------------------------------------------------------------#

import sys

try:
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
except ImportError:
    sys.exit(
        'ERROR: You need argparse!\n' +
        '   install it from http://pypi.python.org/pypi/argparse\n' +
        '   or run \"pip install argparse\".'
    )

from palm_cvd.create_driver import create_driver


class PALMCVDArgumentParser(ArgumentParser):

    def __init__(self):
        super().__init__(
            description='This is the PALM-CVD virtual measurements preprocessor tool\n' +
                        'Developer Support: support@pecanode.com',
            formatter_class=RawTextHelpFormatter,
            add_help=True,
        )
        group = self.add_mutually_exclusive_group(required=True)
        group.add_argument(
            '-i', '--ini-file',
            dest='ini_file',
            action='store',
            default=None,
            required=False,
            help='Define the ini configuration input file.',
            type=str,
            metavar='FILE',
        )
        group.add_argument(
            '-g', '--geo-json',
            dest='geo_json',
            action='store',
            default=None,
            required=False,
            help='Define the geo.json input file.',
            type=str,
            metavar='FILE',
        )
        self.add_argument(
            '-x', '--origin-x',
            dest='origin_x',
            action='store',
            default=-999.9,
            required=False,
            help='Define the EUTM coordinate of the lower-left corner of the PALM domain.',
            type=float,
            metavar='FLOAT',
        )
        self.add_argument(
            '-y', '--origin-y',
            dest='origin_y',
            action='store',
            default=-999.9,
            required=False,
            help='Define the NUTM coordinate of the lower-left corner of the PALM domain.',
            type=float,
            metavar='FLOAT',
        )
        self.add_argument(
            '-z', '--origin-z',
            dest='origin_z',
            action='store',
            default=0.0,
            required=False,
            help='Define the z coordinate (in meters above sea level) of the lower-left corner of the PALM domain.',
            type=float,
            metavar='FLOAT',
        )
        self.add_argument(
            '-o', '--output',
            dest='output',
            action='store',
            default='vm_driver.nc',
            required=True,
            help='Define the path to the output filename.',
            type=str,
            metavar='FILE',
        )


if __name__ == '__main__':
    parser = PALMCVDArgumentParser()
    args = parser.parse_args()
    create_driver(
        output_filename=args.output,
        origin_x=args.origin_x,
        origin_y=args.origin_y,
        origin_z=args.origin_z,
        ini_in=args.ini_file,
        gj_in=args.geo_json,
    )
