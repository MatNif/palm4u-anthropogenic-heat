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

from netCDF4 import Dataset, stringtochar
import numpy as np
import shutil
import sys
import pandas as pd
import utm
import os
import json
import configparser

global global_acronym
global global_author
global global_campaign
global global_comment
global global_contact
global global_data_content
global global_dependencies
global global_institution
global global_keywords
global global_location
global global_references
global global_site
global global_source
global global_name_longitude
global global_name_latitude
global global_featureType
global global_height
global global_palm_version
global global_data_path
global input_file_coord
global number_positions
global input_from_observations
global input_from_csv_xlsx
global coordinates
global vars_to_be_measured
global custom_coordinates
global coordinates_geojson

global_acronym          = " "
global_author           = " "
global_campaign         = " "
global_comment          = " "
global_contact          = " "
global_data_content     = " "
global_dependencies     = " "
global_institution      = " "
global_keywords         = " "
global_location         = " "
global_references       = " "
global_site             = " "
global_source           = " "
global_palm_version     = 6.0
global_name_latitude    = " "
global_name_longitude   = " "
global_height           = -999
global_featureType      = " "
global_data_path        = " "
input_file_coord        = " "
number_positions        = -999
input_from_observations = False
input_from_csv_xlsx     = False
coordinates             = []
vars_to_be_measured     = []
custom_coordinates      = False
coordinates_geojson     = []


def read_gj_file( input_file ):

   global coordinates_geojson

   # Check if the given config file exists.
   if ( os.path.isfile( input_file ) == False ):
      print ("Error. No geojson file " + input_file + " found.")
      sys.exit( " " )

   with open( input_file ) as file_in:
      data = json.load( file_in )

   # Check if this is point data.
   if data['name'] != 'points':
      print( "palm_cvd only supports input of point coordinates at the moment" )
      sys.exit( " " )

   coordinates_geojson = []
   for feature in data['features']:
      dum = feature['geometry']['coordinates'][:]
      coordinates_geojson.append( dum )


# Function to read the config file
def read_config_file( input_config ):

   global global_acronym
   global global_author
   global global_campaign
   global global_comment
   global global_contact
   global global_data_content
   global global_dependencies
   global global_institution
   global global_keywords
   global global_location
   global global_references
   global global_site
   global global_source
   global global_name_longitude
   global global_name_latitude
   global global_featureType
   global global_height
   global global_palm_version
   global global_data_path
   global input_file_coord
   global number_positions
   global input_from_observations
   global input_from_csv_xlsx
   global coordinates
   global vars_to_be_measured
   global custom_coordinates
   global coordinates_geojson
   global input_file_coord
   global number_positions

   # Allow empty settings
   config = configparser.RawConfigParser(allow_no_value=True)

   # Check if the given config file exists.
   if ( os.path.isfile( input_config ) == False ):
      print ("Error. No configuration file " + input_config + " found.")
      sys.exit( " " )

   config.read(input_config)

   for section in range( 0, len( config.sections() ) ):

      current_section = config.sections()[section]

      # Read global attributes which are written into the output file header.
      # Always check if attribute is given in the config file.
      if ( current_section == 'global' ):

         if config.has_option( current_section, 'acronym' ):
            global_acronym = config.get( current_section, 'acronym' )
         if config.has_option( current_section, 'author' ):
            global_author = config.get( current_section, 'author' )
         if config.has_option( current_section, 'campaign' ):
            global_campaign = config.get( current_section, 'campaign' )
         if config.has_option( current_section, 'comment' ):
            global_comment = config.get( current_section, 'comment' )
         if config.has_option( current_section, 'contact_person' ):
            global_contact = config.get( current_section, 'contact_person' )
         if config.has_option( current_section, 'data_content' ):
            global_data_content = config.get( current_section, 'data_content' )
         if config.has_option( current_section, 'dependencies' ):
            global_dependencies = config.get( current_section, 'dependencies' )
         if config.has_option( current_section, 'institution' ):
            global_institution = config.get( current_section, 'institution' )
         if config.has_option( current_section, 'keywords' ):
            global_keywords = config.get( current_section, 'keywords' )
         if config.has_option( current_section, 'location' ):
            global_location = config.get( current_section, 'location' )
         if config.has_option( current_section, 'references' ):
            global_references = config.get( current_section, 'references' )
         if config.has_option( current_section, 'site' ):
            global_site = config.get( current_section, 'site' )
         if config.has_option( current_section, 'source' ):
            global_source = config.get( current_section, 'source' )
         if config.has_option( current_section, 'palm_version' ):
            global_palm_version = config.get( current_section, 'palm_version' )

      # Read data input path for observational data
      elif ( current_section == 'input_from_netcdf' ):

         if config.has_option( current_section, 'data_path' ):
            global_data_path = config.get( current_section, 'data_path' )
            input_from_observations = True

      # Read customized coordinates where virtual measurements shall be taken,
      # as well as the variables that should be sampled.
      elif ( current_section == 'custom_positions' ):

         if config.has_option( current_section, 'number_positions' ):
            number_positions = config.get( current_section, 'number_positions' )

         coordinates = []
         for count in range( 0, int( number_positions ) ):
            coordinates.append( json.loads( config.get( current_section, \
                                                        "coordinates" + str( count + 1 ) ) ) )
            # If coordinates are given, set a global flag.
            custom_coordinates = True

         vars_to_be_measured=[]
         for count in range( 0, int( number_positions ) ):
            vars_to_be_measured.append( json.loads( config.get( current_section, \
                                                    "vars_to_be_measured" + str( count + 1 ) ) ) )

      # Read further information required for file processing of csv / xlsx files
      elif ( current_section == 'csv_xlsx_data' ):

         input_from_csv_xlsx = True
         if config.has_option( current_section, 'data_path' ):
            global_data_path = config.get( current_section, 'data_path' )
         if config.has_option( current_section, 'featureType' ):
            global_featureType = config.get( current_section, 'featureType'     )
         if config.has_option( current_section, 'height' ):
            global_height = float( config.get(current_section, 'height'      )  )
         if config.has_option( current_section, 'name_longitude' ):
            global_name_longitude = config.get( current_section, 'name_longitude' )
         if config.has_option( current_section, 'name_latitude' ):
            global_name_latitude = config.get( current_section, 'name_latitude'   )

         if config.has_option( current_section, 'input_file_coord' ):
            input_file_coord = config.get( current_section, 'input_file_coord' )
         # @todo following rly neccessary?
         if config.has_option( current_section, 'vars_to_be_measured' ):
            if config.get( current_section, 'vars_to_be_measured' ):
               vars_to_be_measured_tmp = config.get( current_section,'vars_to_be_measured' )
               vars_to_be_measured = vars_to_be_measured_tmp.split( ', ' )


# Function to convert an input xslx file to csv
def convert_to_csv( path_to_data ):
    # Read input .xslx file.
    try:
       xlsxfile_in = pd.read_excel(path_to_data + input_file_coord)
    except ImportError:
       print('Input of .csv file failed.\n',
             'Package openpyxl is missing. Use pip to install openpyxl')
       sys.exit( " " )

    xlsxfile_in = pd.read_excel(path_to_data + input_file_coord)

    number_positions = xlsxfile_in.shape[0]

    # Change file name to XX.csv.
    input_file_coord_tmp = input_file_coord.split('.')
    input_file_coord = input_file_coord_tmp[0] + '.csv'

    # Save .xlsx data .csv file.
    csvfile_in = xlsxfile_in.to_csv(path_to_data + input_file_coord, sep=',')
    print('Converting to .csv done. File is temporarily stored in ' + path_to_data + input_file_coord )

    return 0

# Function to create a suitable input data file out of csv input data.
def convert_csv_input_data( path_to_data ):
   # This function is used for the following case:
   # aviable measurement data is given in an csv file, geographical
   # information is aviable in WGS84 decimal which needs to be converted
   # to WGS84 UTM [m], no measurment hight given in the file.
   # Only for timeseries and trajectory!
   #
   # File structure in this case:
   # number of station (consecutively increasing), latitude, longitude

   # Get data from .csv file
   # Check if extra variable names for coordinates are given. (It could be the
   # case for a standard different from [UC]2 standard.)
   if global_name_latitude == " " and global_name_longitude == " ":
      global_name_latitude = 'N_UTM'
      global_name_longitude = 'E_UTM'
   print('Reading data: ' + path_to_data + input_file_coord)
   csvfile_in = pd.read_csv(
      path_to_data + input_file_coord, sep=',', header=0, usecols=[
      global_name_latitude, global_name_longitude], dtype=float )

   number_positions = csvfile_in.shape[0]

   # Save input data in separate arrays.
   x = csvfile_in.loc[:, global_name_longitude]
   longitude = np.array(x[:])
   shape_longitude = longitude.shape

   y = csvfile_in.loc[:, global_name_latitude]
   latitude = np.array(y[:])
   shape_latitude = latitude.shape

   # Size check of latitude and longitude.
   if shape_latitude != shape_longitude:
      print('Error. Input coordinates do not have the same shape.')
      print('Please check your input data.')
      sys.exit( " " )

   # Check if coordinates are given as WGS84 UTM coord. [m], not in [°] or
   # decimal.
   max_tmp = max(latitude)
   if max_tmp < 100.:
      coordinates_tmp = utm.from_latlon(np.array(latitude), np.array(
          longitude))
      utm_longitude = coordinates_tmp[0]  # longitude
      utm_latitude = coordinates_tmp[1]  # latitude
      print('Change coordinates to UTM [m] format finished.')
   else:
      utm_longitude = longitude
      utm_latitude = latitude

   # Save nc file separately.
#   dirtmp = './tmp/'
#   if not os.path.isdir(dirtmp):
#       os.mkdir(dirtmp)


   # Save data in .nc file
   input_file_tmp = input_file_coord.split('.')
   # no spaces in netCDF file name
   input_file_tmp = [u.replace(' ', '_') for u in input_file_tmp]
   input_file = input_file_tmp[0] + '.nc'

   # print('Attention! Do not forget to set the measurement height '
   #       'correctly in /trunk/palm_cvd as' + '\n' +
   #       'well as origin coordinates in the vm config file!\n\n')
   nc_file = Dataset(path_to_data + input_file, 'w', format='NETCDF4')
   # get global from vm config
   nc_file.site = global_site
   nc_file.featureType = global_featureType

   # number of stations/measurement points == shape_longitude
   # Write dimensions.
   nc_file.createDimension('traj', 1)
   nc_file.createDimension('ntime', number_positions)

   # Write Variables.
   E_UTM = nc_file.createVariable('E_UTM', 'f8', ('traj', 'ntime',))           # longitude axis
   N_UTM = nc_file.createVariable('N_UTM', 'f8', ('traj', 'ntime',))           # latitude axis
   z = nc_file.createVariable('z', 'f8', ('traj', 'ntime',))                   # sensor height sea lvl
   station_h = nc_file.createVariable('station_h', 'f8', ('traj', 'ntime',))   # surface height
   height = nc_file.createVariable('height', 'f8', ('traj', 'ntime',))         # sensor height above ground

   # Set heights of surface & measurement device
   z_sensor = np.zeros(shape_latitude)  # measurement height, from sea level
   z_surface = np.zeros(shape_latitude)  # surface height, from sea level
   z_surface[:] = origin_z
   h_sensor = global_height
   if global_featureType == 'trajectory':
       height[:] = h_sensor
   elif global_featureType == 'timeSeries':
       z_sensor[:] = z_surface + h_sensor
       z[:] = z_sensor[:]
       station_h[:] = z_surface[:]

   E_UTM[:] = utm_longitude[:]
   N_UTM[:] = utm_latitude[:]

   nc_file.close()

   print('Coordinates saved in netCDF file.')

   return 0

def create_driver( output_filename, origin_x, origin_y, origin_z, ini_in=None, gj_in=None ):
   # Define strings
   name_featuretype   = 'featureType'
   name_ts            = 'timeSeries'
   name_tspr          = 'timeSeriesProfile'
   name_traj          = 'trajectory'
   name_ntime         = 'ntime'
   name_time          = 'time'
   name_station       = 'station'
   name_traj_dim      = 'traj'
   name_nz            = 'nz'
   name_datacontent   = 'data_content'
   name_eutm          = 'E_UTM'
   name_nutm          = 'N_UTM'
   name_hao           = 'height'
   name_station_h     = 'station_h'
   name_z             = 'z'
   name_soil_sampling = 'soil_sample'
   name_num_stat      = 'number_of_stations'
   name_fill          = '_FillValue'
   name_site          = 'site'
   name_acro          = 'acronym'
   name_content       = 'data_content'
   name_orig_x        = 'origin_x'
   name_orig_y        = 'origin_y'
   name_orig_z        = 'origin_z'

   max_string_len     = 50

   name_measvars      = 'measured_variables'

   non_measurable_vars = ['station_name', 'time', 'time_bounds', 'crs', \
                          'vrs', 'x', 'y', 'z', 'lon', 'lat', 'ntime', 'station', 'traj', \
                          'E_UTM', 'N_UTM', 'height_above_origin', 'station_h', \
                          'traj_name', 'height', 'band_pm_size', 'bands_pm', 'bands_pm_size_bounds', \
                          'bands_pm_size', 'ancillary_detected_layer' ]

   soil_vars            = [ 't_soil', 'm_soil', 'lwc', 'lwcs', 'smp' ]

   dims_out             = [ name_eutm, name_nutm, name_hao, name_z, name_station_h ]

   # Define list of attributes which need to be of type float. In the data set this is not
   # necessarily guranteed.
   atts_float           = [ 'origin_x', 'origin_y', 'origin_z', 'origin_lon', 'origin_lat', 'rotation_angle' ]

   # Define list of default variables that shall be measured at each site
   vars_default         = [ 'u', 'v', 'w', 'theta', 'hus' ]


   #Read config file
   if ini_in is not None:
      read_config_file( ini_in )
   else:
      read_gj_file( gj_in )

   # Initialize counter variable for the number of sites
   num_sites = 0

   # Open output file
   ncfile_out = Dataset( output_filename, "w", format="NETCDF4" )

   # First, add global attributes. Some of them are only defined when
   # palm_cvd is configured from ini file.
   if ini_in is not None:

      ncfile_out.setncattr( 'acronym',        global_acronym      )
      ncfile_out.setncattr( 'author',         global_author       )
      ncfile_out.setncattr( 'campaign',       global_campaign     )
      ncfile_out.setncattr( 'comment',        global_comment      )
      ncfile_out.setncattr( 'contact_person', global_contact      )
      ncfile_out.setncattr( 'data_content',   global_data_content )
      ncfile_out.setncattr( 'dependencies',   global_dependencies )
      ncfile_out.setncattr( 'institution',    global_institution  )
      ncfile_out.setncattr( 'keywords',       global_keywords     )
      ncfile_out.setncattr( 'location',       global_location     )
      ncfile_out.setncattr( 'references',     global_references   )
      ncfile_out.setncattr( 'site',           global_site         )
      ncfile_out.setncattr( 'source',         global_source       )
      ncfile_out.setncattr( 'palm_version',   global_palm_version )

   #Origin coordinates will be always given via the command line.
   ncfile_out.setncattr( 'origin_x', origin_x )
   ncfile_out.setncattr( 'origin_y', origin_y )
   ncfile_out.setncattr( 'origin_z', origin_z )

   # Create universal dimension for the string length.
   ncfile_out.createDimension("string_len", max_string_len)

   # Branch when palm_cvd is configured with ini file.
   file_from_csv = False
   if ini_in is not None:
      data_path = global_data_path

      # Check if observational data is available. In this case,
      # obtain an alphabetically sorted list of input data. List is sorted
      # just for the sake of clarity in the resulting setup file.
      if ( input_from_observations == True  or  input_from_csv_xlsx == True ):
         # Is there a specific (csv) file given? This case, make a temporary working copy of the
         # input files in ./tmp
         if input_file_coord != ' '  and  ( '.xlsx' in input_file_coord or  '.csv' in input_file_coord ):
            if not os.path.isdir('./tmp'):
               os.mkdir('./tmp')

            os.system( 'cp ' + data_path + input_file_coord + ' ' + './tmp' )

            # Reset data path to temporary working directory
            data_path = './tmp/'

         if input_file_coord != ' '  and  '.xlsx' in input_file_coord:
             # excel_input = True
             print('Excel file in .xlsx format found, converting to .csv.')
             convert_to_csv( data_path )
             convert_csv_input_data( data_path )
             file_from_csv = True
         elif input_file_coord != ' '  and  '.csv' in input_file_coord:
             convert_csv_input_data( data_path )
             file_from_csv = True

         # Remove temporarily stored .csv or .xlsx files if necessary.
         if file_from_csv == True:
            list_files = sorted( os.listdir( data_path ) )
            for filename in list_files:
               if '.csv' in filename  or  '.xlsx' in filename:
                  os.remove( data_path + filename )

         list_input_data = sorted( os.listdir( data_path ) )

      # Set initial counter ID
      counter_id = 1

      if ( input_from_observations  or  input_from_csv_xlsx ):

         # Run loop over all listed input data. Depending on the data set, this could be
         # a list of files or a list of subdirectories.
         # This is done to reduce the number of virtual measurements in the model. Each
         # virtual measurement has an overhead and consumes memory.
         sites = []
         input_files = []
         input_files_orig = []
         for dirname in list_input_data:

             data_file = data_path + dirname

             if ( os.path.isdir(data_file) == True ):
                # Directory may contain various file versions.
                # Take the one with highest cycle number.
                highest_cycle_nr = 0
                for filename in os.listdir(data_file):
                   start_seq = len( filename ) - 6
                   end_seq   = len( filename ) - 3
                   if int( filename[start_seq:end_seq] ) > highest_cycle_nr:
                      highest_cycle_nr = int(filename[start_seq:end_seq])
                      latest_file      = filename
                input_file = data_file + "/" + latest_file
                input_file_orig = latest_file
             else:
                input_file = data_file
                input_file_orig = dirname

             # Open the NetCDF file
             ncfile_in = Dataset( input_file, "r", format="NETCDF4", encoding='ascii')
             input_files.append(input_file)
             input_files_orig.append(input_file_orig)

         # Gather all files according to their feature type and all sites for the respective feature type
         files_traj  = []
         files_ts    = []
         files_tspr  = []
         sites_traj  = []
         sites_ts    = []
         sites_tspr  = []
         for input_file in input_files:
            ncfile_in = Dataset( input_file, "r", format="NETCDF4", encoding='ascii' )

            for att in ncfile_in.ncattrs():
               if ( att == name_featuretype ):
                  feature = ncfile_in.getncattr(att)
               if ( att == name_site ):
                  site = ncfile_in.getncattr(att)

            if ( feature == name_traj ):
               files_traj.append(input_file)
            elif ( feature == name_ts ):
               files_ts.append(input_file)
            else:
               files_tspr.append(input_file)

            if ( feature == name_traj  and  site not in sites_traj ):
               sites_traj.append(site)
            if ( feature == name_ts  and  site not in sites_ts ):
               sites_ts.append(site)
            if ( feature == name_tspr  and  site not in sites_tspr ):
               sites_tspr.append(site)

            ncfile_in.close()


         for site_traj in sites_traj:
            # For the given site already define the featureTpye and site
            ncfile_out.setncattr( name_featuretype + str(counter_id), name_traj )
            ncfile_out.setncattr( name_site + str(counter_id), site_traj )

            # Define the number of coordinates for the site
            num_coords = 0

            e_utm_traj = np.array([])
            n_utm_traj = np.array([])
            h_traj     = np.array([])
            measured_variables = ['u', 'v', 'w', 'theta', 'hus']
            for input_file in files_traj:
               print( "traj", input_file, " ", site_traj )
               ncfile_in = Dataset( input_file, "r", format="NETCDF4", encoding='ascii' )
               for att in ncfile_in.ncattrs():
                  if ( att == name_site ):
                     site = ncfile_in.getncattr(att)

               if ( site == site_traj ):
                  orig_x = -999.9
                  orig_y = -999.9
                  orig_z = -999.9
                  for att in ncfile_in.ncattrs():
                     if ( att == name_orig_x ):
                        orig_x = ncfile_in.getncattr(att)
                     if ( att == name_orig_y ):
                        orig_y = ncfile_in.getncattr(att)
                     if ( att == name_orig_z ):
                        orig_z = ncfile_in.getncattr(att)

                  ntime = len( ncfile_in.dimensions[name_ntime]    )
                  ntraj = len( ncfile_in.dimensions[name_traj_dim] )

                  num_coords += ntime * ntraj
                  # Gather UTM and height coordinates and merge them into one array. Further, gather
                  # the variables that shall be sampled. Coordinates are checked to for NaN values and
                  # are tranformed to arithmetric numbers. Further, 2D input array is transformed into
                  # a 1D array.
                  for var in ncfile_in.variables.keys():
                     if ( var in dims_out  and  var == name_eutm ):
                        e_utm_traj = np.append(e_utm_traj, np.nan_to_num( ncfile_in.variables[var][:,:] ).flatten())
                        #e_utm_traj.append( np.nan_to_num( ncfile_in.variables[var][:,:] ).flatten() )
                     if ( var in dims_out  and  var == name_nutm ):
                        n_utm_traj = np.append(n_utm_traj, np.nan_to_num( ncfile_in.variables[var][:,:] ).flatten())
                     if ( var in dims_out  and  var == name_hao ):
                        h_traj = np.append(h_traj, np.nan_to_num( ncfile_in.variables[var][:,:] ).flatten())

                     if ( var not in non_measurable_vars  and  \
                          var not in vars_default         and  \
                          var not in measured_variables ):
                        measured_variables.append(var)

               ncfile_in.close()

            # After all files for the current site are processed, write the origin-coordinates for x,y,z
            ncfile_out.setncattr( name_orig_x + str(counter_id), orig_x )
            ncfile_out.setncattr( name_orig_y + str(counter_id), orig_y )
            ncfile_out.setncattr( name_orig_z + str(counter_id), orig_z )
            # Create the dimensions
            ncfile_out.createDimension( name_station + str(counter_id), num_coords )

            temp_traj = ncfile_out.createVariable( name_eutm + str(counter_id), float, name_station + str(counter_id) )
            temp_traj[:] = e_utm_traj

            temp_traj = ncfile_out.createVariable( name_nutm + str(counter_id), float, name_station + str(counter_id) )
            temp_traj[:] = n_utm_traj

            temp_traj = ncfile_out.createVariable( name_hao  + str(counter_id), float, name_station + str(counter_id) )
            temp_traj[:] = h_traj

            # Check if any of the measured variables is a soil variable. Set flag accordingly.
            soil = False
            for var in measured_variables:
               if ( var in soil_vars ):
                  soil = True
            # Write soil flag
            ncfile_out.setncattr( name_soil_sampling + str( counter_id), np.int8(soil) )

            # Create dimension for sample-variable string
            ncfile_out.createDimension( "nvar"+ str(counter_id), len( measured_variables ) )

            measured_var = ncfile_out.createVariable( 'measured_variables' + str(counter_id), 'S1', \
                                                      ("nvar" + str(counter_id), "string_len") ) # must be NC_CHAR

            # Write the variables to the file
            for counter, meas in enumerate( measured_variables ):
               measured_var[counter] = stringtochar( np.array( meas,"S%s"%(max_string_len) ) )

            # Increment the counter
            counter_id += 1


         for site_tspr in sites_tspr:
            # For the given site already define the featureTpye and site
            ncfile_out.setncattr( name_featuretype + str(counter_id), name_tspr )
            ncfile_out.setncattr( name_site + str(counter_id), site_tspr )

            # Define the number of coordinates for the site
            num_coords = 0
            e_utm_tspr     = np.array([])
            n_utm_tspr     = np.array([])
            station_h_tspr = np.array([])
            z_tspr         = np.array([])

            measured_variables = ['u', 'v', 'w', 'theta', 'hus']
            for input_file in files_tspr:
               print( "tspr", input_file, " ", site_tspr )
               ncfile_in = Dataset( input_file, "r", format="NETCDF4", encoding='ascii' )
               for att in ncfile_in.ncattrs():
                  if ( att == name_site ):
                     site = ncfile_in.getncattr(att)

               if ( site == site_tspr ):
                  for att in ncfile_in.ncattrs():
                     if ( att == name_orig_x ):
                        orig_x = ncfile_in.getncattr(att)
                     if ( att == name_orig_y ):
                        orig_y = ncfile_in.getncattr(att)
                     if ( att == name_orig_z ):
                        orig_z = ncfile_in.getncattr(att)

                  nstation = len( ncfile_in.dimensions[name_station] )
                  ntime    = len( ncfile_in.dimensions[name_ntime]   )
                  nz       = len( ncfile_in.dimensions[name_nz]   )

                  num_coords += nstation * ntime * nz
                  # Gather UTM and height coordinates and merge them into one array. Further, gather
                  # the variables that shall be sampled. Coordinates are checked to for NaN values and
                  # are tranformed to arithmetric numbers. Further, 2D input array is transformed into
                  # a 1D array.
                  for var in ncfile_in.variables.keys():
                     tspr_tmp1 = np.zeros((nstation))
                     tspr_tmp2 = np.zeros((ntime*nz))
                     if ( var in dims_out  and  var == name_eutm ):
                        tspr_tmp1 = np.nan_to_num( ncfile_in.variables[var][:] )
                        for ns in range(0,int(nstation)):
                           tspr_tmp2[:] = tspr_tmp1[ns]
                           e_utm_tspr = np.append(e_utm_tspr, tspr_tmp2)
                     if ( var in dims_out  and  var == name_nutm ):
                        tspr_tmp1 = np.nan_to_num( ncfile_in.variables[var][:] )
                        for ns in range(0,int(nstation)):
                           tspr_tmp2[:] = tspr_tmp1[ns]
                           n_utm_tspr = np.append(n_utm_tspr, tspr_tmp2)
                     if ( var in dims_out  and  var == name_z ):
                        z_tspr_tmp = np.nan_to_num( ncfile_in.variables[var][:,:,:] )
                        z_tspr = np.append(z_tspr, np.concatenate( z_tspr_tmp ))
                     if ( var in dims_out  and  var == name_station_h ):
                        tspr_tmp1 = np.nan_to_num( ncfile_in.variables[var][:] )
                        for ns in range(0,int(nstation)):
                           tspr_tmp2[:] = tspr_tmp1[ns]
                           station_h_tspr = np.append(station_h_tspr, tspr_tmp2)

                     if ( var not in non_measurable_vars  and  \
                          var not in vars_default         and  \
                          var not in measured_variables ):
                        measured_variables.append(var)

               ncfile_in.close()

            # After all files for the current site are processed, write the origin-coordinates for x,y,z
            ncfile_out.setncattr( name_orig_x + str(counter_id), orig_x )
            ncfile_out.setncattr( name_orig_y + str(counter_id), orig_y )
            ncfile_out.setncattr( name_orig_z + str(counter_id), orig_z )
            # Create the dimensions
            ncfile_out.createDimension( name_station + str(counter_id), num_coords )

            temp_tspr = ncfile_out.createVariable( name_eutm      + str(counter_id), float, name_station + str(counter_id) )
            temp_tspr[:] = e_utm_tspr

            temp_tspr = ncfile_out.createVariable( name_nutm      + str(counter_id), float, name_station + str(counter_id) )
            temp_tspr[:] = n_utm_tspr

            temp_tspr = ncfile_out.createVariable( name_z         + str(counter_id), float, name_station + str(counter_id) )
            temp_tspr[:] = z_tspr

            temp_tspr = ncfile_out.createVariable( name_station_h + str(counter_id), float, name_station + str(counter_id) )
            temp_tspr[:] = station_h_tspr

            # Check if any of the measured variables is a soil variable. Set flag accordingly.
            soil = False
            for var in measured_variables:
               if ( var in soil_vars ):
                  soil = True
            # Write soil flag
            ncfile_out.setncattr( name_soil_sampling + str( counter_id), np.int8(soil) )

            # Create dimension for sample-variable string
            ncfile_out.createDimension( "nvar"+ str(counter_id), len( measured_variables ) )

            measured_var = ncfile_out.createVariable( 'measured_variables' + str(counter_id), 'S1', \
                                                      ("nvar" + str(counter_id), "string_len") ) # must be NC_CHAR

            # Write the variables to the file
            for counter, meas in enumerate( measured_variables ):
               measured_var[counter] = stringtochar( np.array( meas,"S%s"%(max_string_len) ) )

            # Increment the counter
            counter_id += 1


         for site_ts in sites_ts:
            # For the given site already define the featureTpye and site
            ncfile_out.setncattr( name_featuretype + str(counter_id), name_ts )
            ncfile_out.setncattr( name_site + str(counter_id), site_ts )

            # Define the number of coordinates for the site
            num_coords = 0
            e_utm_ts     = np.array([])
            n_utm_ts     = np.array([])
            station_h_ts = np.array([])
            z_ts         = np.array([])

            measured_variables = ['u', 'v', 'w', 'theta', 'hus']
            for input_file in files_ts:
               print( "ts", input_file, " ", site_ts )
               ncfile_in = Dataset( input_file, "r", format="NETCDF4", encoding='ascii' )
               for att in ncfile_in.ncattrs():
                  if ( att == name_site ):
                     site = ncfile_in.getncattr(att)

               if ( site == site_ts ):

                  for att in ncfile_in.ncattrs():
                     if ( att == name_orig_x ):
                        orig_x = ncfile_in.getncattr(att)
                     if ( att == name_orig_y ):
                        orig_y = ncfile_in.getncattr(att)
                     if ( att == name_orig_z ):
                        orig_z = ncfile_in.getncattr(att)

                  nstation = len( ncfile_in.dimensions[name_station]    )
                  num_coords += nstation
                  # Gather UTM and height coordinates and merge them into one array. Further, gather
                  # the variables that shall be sampled. Coordinates are checked to for NaN values and
                  # are tranformed to arithmetric numbers.
                  for var in ncfile_in.variables.keys():
                     if ( var in dims_out  and  var == name_eutm ):
                        e_utm_ts = np.append(e_utm_ts, np.nan_to_num( ncfile_in.variables[var][:] ))
                     if ( var in dims_out  and  var == name_nutm ):
                        n_utm_ts = np.append(n_utm_ts, np.nan_to_num( ncfile_in.variables[var][:] ))
                     if ( var in dims_out  and  var == name_z ):
                        z_ts = np.append(z_ts, np.nan_to_num( ncfile_in.variables[var][:] ))
                     if ( var in dims_out  and  var == name_station_h ):
                        station_h_ts = np.append(station_h_ts, np.nan_to_num( ncfile_in.variables[var][:] ))

                     if ( var not in non_measurable_vars  and  \
                          var not in vars_default         and  \
                          var not in measured_variables ):
                        measured_variables.append(var)

               ncfile_in.close()

            # After all files for the current site are processed, write the origin-coordinates for x,y,z
            ncfile_out.setncattr( name_orig_x + str(counter_id), orig_x )
            ncfile_out.setncattr( name_orig_y + str(counter_id), orig_y )
            ncfile_out.setncattr( name_orig_z + str(counter_id), orig_z )
            # Create the dimensions
            ncfile_out.createDimension( name_station + str(counter_id), num_coords )

            temp_ts = ncfile_out.createVariable( name_eutm      + str(counter_id), float, name_station + str(counter_id) )
            temp_ts[:] = e_utm_ts

            temp_ts = ncfile_out.createVariable( name_nutm      + str(counter_id), float, name_station + str(counter_id) )
            temp_ts[:] = n_utm_ts

            temp_ts = ncfile_out.createVariable( name_z         + str(counter_id), float, name_station + str(counter_id) )
            temp_ts[:] = z_ts

            temp_ts = ncfile_out.createVariable( name_station_h + str(counter_id), float, name_station + str(counter_id) )
            temp_ts[:] = station_h_ts

            # Check if any of the measured variables is a soil variable. Set flag accordingly.
            soil = False
            for var in measured_variables:
               if ( var in soil_vars ):
                  soil = True
            # Write soil flag
            ncfile_out.setncattr( name_soil_sampling + str( counter_id), np.int8(soil) )

            # Create dimension for sample-variable string
            ncfile_out.createDimension( "nvar"+ str(counter_id), len( measured_variables ) )

            measured_var = ncfile_out.createVariable( 'measured_variables' + str(counter_id), 'S1', \
                                                      ("nvar" + str(counter_id), "string_len") ) # must be NC_CHAR

            # Write the variables to the file
            for counter, meas in enumerate( measured_variables ):
               measured_var[counter] = stringtochar( np.array( meas,"S%s"%(max_string_len) ) )

            # Increment the counter
            counter_id += 1



         # Store the number of observational sites
         num_sites = len( sites_traj ) + len( sites_ts ) + len( sites_tspr )


      # Now process the customized input data. Please note, at the moment only timeseries are
      # are possible.
      if ( custom_coordinates ):
         num_sites = counter_id - 1
         for coord in coordinates:
            # Define mandatory attributes
            ncfile_out.setncattr( name_featuretype + str(counter_id),  \
                                  name_ts )
            ncfile_out.setncattr( name_site        + str(counter_id),  \
                                  "custom"         + str(counter_id - num_sites) )
            ncfile_out.setncattr( name_orig_x      + str(counter_id),  \
                                  coord[0] )
            ncfile_out.setncattr( name_orig_y      + str(counter_id),  \
                                  coord[1] )
            ncfile_out.setncattr( name_orig_z      + str(counter_id),  \
                                  0.0 )

            # Define dimensions
            ntime = 1
            nstat = 1
            ncfile_out.createDimension( name_ntime   + str(counter_id), ntime )
            ncfile_out.createDimension( name_station + str(counter_id), nstat )

            # Define coordinate variables
            temp_ts = ncfile_out.createVariable( name_eutm      + str(counter_id), \
                                                 float,                            \
                                                 name_station   + str(counter_id) )
            temp_ts[:] = np.array( coord[0] )

            temp_ts = ncfile_out.createVariable( name_nutm      + str(counter_id), \
                                                 float,                            \
                                                 name_station   + str(counter_id) )
            temp_ts[:] = np.array( coord[1] )

            temp_ts = ncfile_out.createVariable( name_z         + str(counter_id), \
                                                 float,                            \
                                                 name_station   + str(counter_id) )
            temp_ts[:] = np.array( coord[2] )

            temp_ts = ncfile_out.createVariable( name_station_h + str(counter_id), \
                                                 float,                            \
                                                 name_station   + str(counter_id) )
            temp_ts[:] = np.array( 0.0 )


            counter_id += 1

         # Reset counter variable
         counter_id = num_sites + 1

         # check if variables are prescribed. If so, prepare final output string
         # stored in measured_variables.
         if ( vars_to_be_measured ):

            for custom_vars in vars_to_be_measured:

               measured_variables = []
               for var in vars_default:
                  measured_variables.append(var)

               # Check if given variables are already in the default variables.
               # If not, extend.
               for var in custom_vars:
                   if ( var  not in  measured_variables ):

                      measured_variables.append(var)

               ncfile_out.createDimension( "nvar"+ str(counter_id), \
                                           len( measured_variables ) )

               measured_var = ncfile_out.createVariable( 'measured_variables' + str(counter_id), 'S1', \
                                                        ("nvar" + str(counter_id), "string_len") ) # must be NC_CHAR

               # Write the variables to the file
               for counter, meas in enumerate( measured_variables ):
                  measured_var[counter] = stringtochar( np.array( meas,"S%s"%(max_string_len) ) )

               # Add soil attribute for the current measurement.
               soil = False
               if ( any( var == soil_vars for var in measured_variables) ):
                  soil = True

               # Write soil flag
               ncfile_out.setncattr( name_soil_sampling + str( counter_id), np.int8(soil) )

               # Increment counter variable
               counter_id += 1

               del ( measured_variables[:] )

         # Add the number of customized sites.
         num_sites = counter_id - 1

   # Create driver from geo-json input.
   if gj_in is not None:
      num_sites = len(coordinates_geojson)
      counter_id = 1

      for coord in coordinates_geojson:
         # Define mandatory attributes
         ncfile_out.setncattr( name_featuretype + str(counter_id), name_ts )
         ncfile_out.setncattr( name_site        + str(counter_id), "geojson" + str(counter_id) )

         ncfile_out.setncattr( name_orig_x      + str(counter_id), coord[0] )
         ncfile_out.setncattr( name_orig_y      + str(counter_id), coord[1] )
         ncfile_out.setncattr( name_orig_z      + str(counter_id), 0.0 )

         # Define dimensions
         ntime = 1
         nstat = 1
         ncfile_out.createDimension( name_ntime   + str(counter_id), ntime )
         ncfile_out.createDimension( name_station + str(counter_id), nstat )

         # Define coordinate variables
         temp_ts = ncfile_out.createVariable( name_eutm      + str(counter_id), \
                                              float,                            \
                                              name_station   + str(counter_id) )
         temp_ts[:] = np.array( coord[0] )

         temp_ts = ncfile_out.createVariable( name_nutm      + str(counter_id), \
                                              float,                            \
                                              name_station   + str(counter_id) )
         temp_ts[:] = np.array( coord[1] )

         temp_ts = ncfile_out.createVariable( name_z         + str(counter_id), \
                                              float,                            \
                                              name_station   + str(counter_id) )
         temp_ts[:] = np.array( coord[2] )

         temp_ts = ncfile_out.createVariable( name_station_h + str(counter_id), \
                                              float,                            \
                                              name_station   + str(counter_id) )
         temp_ts[:] = np.array( 0.0 )


         # Define variables that shall be sampled.
         # Note, at the moment this is simple the default list.
         ncfile_out.createDimension( "nvar"+ str(counter_id), len( vars_default ) )

         measured_var = ncfile_out.createVariable( 'measured_variables' + str(counter_id), 'S1', \
                                                   ("nvar" + str(counter_id), "string_len") ) # must be NC_CHAR

         # Write the variables to the file
         for counter, meas in enumerate( vars_default ):
            measured_var[counter] = stringtochar( np.array( meas,"S%s"%(max_string_len) ) )

         # Write soil flag
         soil = False
         ncfile_out.setncattr( name_soil_sampling + str( counter_id), np.int8(soil) )

         counter_id += 1

   # Finally, write the total number of sites to the output file
   ncfile_out.setncattr( name_num_stat, num_sites )

   # Clean-up temporary files
   if file_from_csv == True:
      shutil.rmtree('./tmp/')

   print( " " )
   print( "*** palm_cvd has been finished. You can find the output file under: " )
   print( "    " + output_filename )

