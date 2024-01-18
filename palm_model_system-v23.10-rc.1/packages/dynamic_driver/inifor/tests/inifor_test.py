#!/usr/bin/env python
#------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2017-2021 Leibniz Universitaet Hannover
# Copyright 2017-2021 Deutscher Wetterdienst Offenbach
#------------------------------------------------------------------------------#
#
# Authors:
# --------
# @author Eckhard Kadasch
#
# Description:
# ------------
# This script runs the INIFOR integration test. It runs the 'inifor' executable
# found in the test directory (the directory this script resides in) and
# compares the produced dynamic driver file with a referene file.
#------------------------------------------------------------------------------#
from netCDF4 import Dataset
from numpy import max as numpy_max
from os import path
from subprocess import run
from sys import argv, exit

EXIT_CODE_OK = 0
EXIT_CODE_FAIL = 1
ATTRIBUTES_TO_CHECK = set(['origin_lon', 'origin_lat', 'origin_z'])
DEBUGGING = 'on'
REFERENCE_FILENAME = 'reference.nc'
TEST_FILENAME = 'test.nc'
    

def main(argv):

    test_dir = directory_of_this_script()
    case_dir = test_dir + '/cases/20130720_cosmo'
    path_to_reference_file = case_dir + '/' + REFERENCE_FILENAME
    path_to_test_file = case_dir + '/' + TEST_FILENAME
    inifor_binary = test_dir + '/../build/bin/inifor'
    reference_call = f"{inifor_binary} \
--namelist {case_dir}/namelist \
--path {case_dir} \
--date 2013072012 \
--elevation 42.0 \
--averaging-angle 0.02 \
--output {path_to_test_file}"
    run(['rm', '-f', path_to_test_file], cwd=test_dir)
    print(reference_call)
    run(reference_call.split(' '), cwd=test_dir)

    try:

       with Dataset(path_to_reference_file, 'r') as reference_file, \
            Dataset(path_to_test_file, 'r') as test_file:

            print_debug('Comparing test file: %s' % path_to_test_file)
            print_debug('with reference file: %s' % path_to_reference_file)

            exit_code = compare_files(reference_file, test_file)

    except OSError as e:

        print_debug('%s: %s' % (e.strerror, e.filename.decode('utf-8')))
        exit_code = EXIT_CODE_FAIL

    print_test_result(exit_code)
    return exit_code


def compare_files(reference_file, test_file):

    test_result, missing_items = test_file_contains_reference_attributes(test_file)
    try:
        assert test_result
        print_debug('All required global attributes are present.')
    except AssertionError:
        print_debug('The following global attributes are missing:')
        print_debug(missing_items)
        return EXIT_CODE_FAIL

    try:
        assert all_attributes_match(reference_file, test_file)
        print_debug('All attributes match.')
    except AssertionError:
        print_debug('Some global attributes do not match.')
        return EXIT_CODE_FAIL

    test_result, missing_items = test_file_contains_reference_variables(reference_file, test_file)
    try:
        assert test_result
        print_debug('All variables are present.')
    except AssertionError:
        print_debug('The following variables are missing:')
        print_debug(missing_items)
        return EXIT_CODE_FAIL

    try:
        assert all_variables_match(reference_file, test_file)
        print_debug('All variables match.')
    except AssertionError:
        print_debug('Some variables do not match.')
        return EXIT_CODE_FAIL

    return EXIT_CODE_OK


def test_file_contains_reference_variables(reference_file, test_file):
    """
    Check if any netCDF variable contained in the reference file is missing
    from the test file. Additional variables in the test file are permitted.
    """
    reference_vars = set(reference_file.variables.keys())
    test_vars = set(test_file.variables.keys())

    return set_contains_reference_items(
        reference_set=reference_vars,
        test_set=test_vars
    )


def all_variables_match(file_a, file_b):
    vars_a = set(file_a.variables.keys())
    vars_b = set(file_b.variables.keys())
    shared_vars = vars_a.intersection(vars_b)
    true_if_all_match = True
    
    for var in shared_vars:
        try:
            assert (file_a.variables[var][:] == file_b.variables[var][:]).all()
            print_debug('  %s: data matches' % var)
        except AssertionError:
            max_diff = numpy_max(file_a.variables[var][:] - file_b.variables[var][:])
            print_debug('  %s: max error = %f' % (var, max_diff))
            true_if_all_match = False

    return true_if_all_match


def test_file_contains_reference_attributes(test_file):
    """
    Check if any required global netCDF attribute (listed in
    ATTRIBUTES_TO_CHECK) are missing from the test file. Additional attributes
    are permitted.
    """
    return set_contains_reference_items(
        reference_set=ATTRIBUTES_TO_CHECK,
        test_set=set(test_file.ncattrs())
    )


def all_attributes_match(reference_file, test_file):
    for attribute in ATTRIBUTES_TO_CHECK:
        reference_value = float(reference_file.getncattr(attribute))
        test_value = float(test_file.getncattr(attribute))
        try:
            assert (test_value == reference_value)
            print_debug('  %s: value matches' % attribute)
        except AssertionError:
            diff = test_value - reference_value
            print_debug('  %s: error = %f' % (attribute, diff))
            return False

    return True


def set_contains_reference_items(reference_set, test_set):
    missing_items = reference_set.difference(test_set)
    if len(missing_items) == 0:
        return True, missing_items
    else:
        return False, missing_items


def directory_of_this_script():
    return path.dirname(path.realpath(__file__))


def print_debug(message):
    if DEBUGGING == 'on':
        print('inifor_test: %s' % message)


def print_test_result(exit_code):

    if exit_code == EXIT_CODE_OK:
        print_debug('SUCCESS: INIFOR passed all tests.')
    else:
        print_debug('FAILURE: INIFOR failed the above test.')


if __name__ == '__main__':
    exit(main(argv))

