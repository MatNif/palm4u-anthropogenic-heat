#!/usr/bin/env python3
# --------------------------------------------------------------------------------#
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
# --------------------------------------------------------------------------------#
#
# Description:
# ------------
"""Merge virtual measurement output.

Removes empty time stamps from the netCDF files and concatenates files
from several restart files into one file.

Example:
module load nco anaconda3
python3 postprocess_vm_measurements.py my_palm_simulation/OUTPUT
"""
#
# @Authors Matthias Sühring (suehring@muk.uni-hannover.de)
#          Tobias Gronemeier (gronemeier@muk.uni-hannover.de)
# --------------------------------------------------------------------------------#


import argparse
import subprocess
import os
import sys
try:
    import numpy as np
except ImportError:
    sys.exit(
        'package "numpy" is required but not installed! Run\n'
        + 'python -m pip install --user numpy\nto install it.')
try:
    from netCDF4 import Dataset
except ImportError:
    sys.exit(
        'package "netCDF4" is required but not installed! Run\n'
        + 'python -m pip install --user netCDF4\nto install it.')


def concatenate(files_per_site, sites, output_directory, overwrite_file=False):
    """Concatenate netCDF files via ncrcat.

    Concatenate a list of netCDF files using NCO command 'ncrcat'.
    Return value: output file
    """

    if not os.path.isdir(output_directory):
        mkdir = os.mkdir(output_directory)

    if output_directory[-1] != '/':
        output_directory += '/'

    for site_index, file_list in enumerate(files_per_site):

        ncrcat_command = "ncrcat"

        if overwrite_file:
            ncrcat_command += " -O"

        for file_name in file_list:
            ncrcat_command += " " + file_name

        # Check if output file already exists
        output_file = output_directory + sites[site_index]
        if not overwrite_file and os.path.isfile(output_file):
            for i in range(1000):
                output_file = output_directory + sites[site_index] + "_{:03d}".format(i)
                if not os.path.isfile(output_file):
                    break
                elif i == 999:
                    raise IOError("could not guarantee non overwriting output file: {}".format(
                            output_file))
        ncrcat_command += " " + output_file

        print(ncrcat_command)
        ncrcat_output = subprocess.run(ncrcat_command, shell=True, check=True)

    return output_file


def truncate(input_file, time_index_shift=0, overwrite_file=False):
    """Truncate netCDF files via ncrcat.

    Truncate all time dimensions of the input file and convert them to
    record dimensions. The output is saved to 'input_file.trunc' or to
    'input_file.trunc.nc' if the input_file has a '.nc' extension.
    If "overwrite_file" is true, write output directly to input_file.
    Shift the time index variables by time_index_shift.

    Return values:
        highest time index of time dimension in output file
        output-file name
    """

    # Gather information about time coordinate in file
    ncfile = Dataset(input_file, "r")
    time_dim = ncfile.dimensions["ntime"]
    time_var = ncfile.variables["time"][:, :]
    time_mask = ~np.ma.getmaskarray(time_var)
    start_index = 0

    soil = any([var == "time_soil" for var in ncfile.variables.keys()])

    end_index1 = 10E+10
    end_index2 = 10E+10
    end_index3 = 10E+10
    end_index4 = 10E+10
    if np.any(time_mask is False):
        end_ind = np.where(time_mask is False)
        end_index1 = end_ind[0][0] - 1
        cut = True
    else:
        end_index2 = time_var.shape[0]
        cut = False

    # Further check for faulty time values.
    # time should not be negative.
    if np.any(time_var <= 0):
        end_ind = np.where(time_var <= 0 )
        end_index3 = end_ind[0][0] - 1
        cut = True

    # time should not be unrealistically large.
    if np.any(time_var >= 10E6):
        end_ind = np.where(time_var >= 10E6 )
        end_index4 = end_ind[0][0] - 1
        cut = True

    end_index = min( end_index1, end_index2, end_index3, end_index4 )

    # time should increase monotonically.
    for t in range( 1, end_index ):
       if ( time_var[t,0] <= time_var[t-1,0] ):
          end_index = t-1
          break

    for att in ncfile.ncattrs():
        if (att == "site"):
            site = ncfile.getncattr(att)
        if (att == "featureType"):
            feat = ncfile.getncattr(att)

    # Compose nco commands
    ncks_command = "ncks"

    if overwrite_file:
        ncks_command += " -O"
        output_file = input_file
    else:
        # Add '.trunc' to file name before '.nc' file extension
        output_file, file_extension = os.path.splitext(input_file)
        if file_extension != '.nc':
            output_file += file_extension + '.trunc'
        else:
            output_file += '.trunc' + file_extension
        if os.path.isfile(output_file):
            raise IOError("truncated file already exists: {}".format(output_file))

    if cut:
        # set dimension limits
        ncks_command += " -d ntime,{0},{1}".format(start_index, end_index)
        if soil:
            ncks_command += " -d ntime_soil,{0},{1}".format(start_index, end_index)

    # convert time into record dimension
    time_is_limited = not time_dim.isunlimited()
    if time_is_limited:
        ncks_command += " --mk_rec_dmn"
        ncks_command += " ntime"

    if cut or time_is_limited:
        # set input and output file
        ncks_command += " {0} {1}".format(input_file, output_file)

        # execute ncks
        print(ncks_command)
        ncks_output = subprocess.run(ncks_command, shell=True, check=True, stdout=subprocess.PIPE)

        new_input_file = output_file

    else:
        new_input_file = input_file

    # If soil is present, also convert soil time to record dimension
    # (must be done separately due to NCO limitations)
    if soil:
        soil_time_is_limited = not ncfile.dimensions["ntime_soil"].isunlimited()
        if soil_time_is_limited:

            ncks_command = "ncks -O --mk_rec_dmn ntime_soil {0} {1}".format(
                    new_input_file, output_file)
            print(ncks_command)
            ncks_output = subprocess.run(
                    ncks_command, shell=True, check=True, stdout=subprocess.PIPE)

            new_input_file = output_file

    # Add time shift to ntime variables
    if time_index_shift != 0:
        ncap2_command = "ncap2 -O -s 'ntime=ntime+{}' ".format(time_index_shift)
        if soil:
            ncap2_command += " -s 'ntime_soil=ntime_soil+{}' ".format(time_index_shift)

        ncap2_command += " {0} {1}".format(new_input_file, output_file)

        print(ncap2_command)
        ncap2_output = subprocess.run(ncap2_command, shell=True, check=True)

        end_index += time_index_shift
        new_input_file = output_file

    return new_input_file, end_index


def main(base_input_directory, output_directory, overwrite_file=False):

    if base_input_directory[-1] != '/':
        base_input_directory += '/'

    if output_directory[-1] != '/':
        output_directory += '/'

    # Get directory list
    input_directory_list = [
            base_input_directory + directory + '/' for directory in
            sorted(os.listdir(base_input_directory))]

    # Obtain list of sites that need to be processed
    sites = sorted(os.listdir(input_directory_list[0]))

    files_per_site_and_directory = [[None] * len(input_directory_list) for i in range(len(sites))]

    # Truncate each file and save end index of time dimension
    for site_index, site_name in enumerate(sites):

        start_index = 0
        for dir_index, directory in enumerate(input_directory_list):
            files_per_site_and_directory[site_index][dir_index], end_index = \
                    truncate(directory + site_name, start_index, overwrite_file)
            start_index = end_index

    # Concatenate all files
    file_concatenated = concatenate(
            files_per_site_and_directory, sites, output_directory, overwrite_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Merge virtual measurement output from multiple PALM run cycles',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            'input',
            metavar='IN',
            help='PALM output directory containing virtual measurements')
    parser.add_argument(
            '--out',
            '-o',
            metavar='OUT',
            default='./merge',
            help='Output directory to store merged data')
    parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Overwrite input files with output files')

    args = parser.parse_args()

    main(args.input, output_directory=args.out, overwrite_file=args.overwrite)
