!> @virtual_measurement_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
! Copyright 2022-2022 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
!
! Authors:
! --------
! @author Matthias Suehring
! @author Tobias Gronemeier
!
! Description:
! ------------
!> The module acts as an interface between 'real-world' observations and model simulations.
!> Virtual measurements will be taken in the model at the coordinates representative for the
!> 'real-world' observation coordinates. More precisely, coordinates and measured quanties will be
!> read from a NetCDF file which contains all required information. In the model, the same
!> quantities (as long as all the required PALM-modules are switched-on) will be sampled at the
!> respective positions and output into an extra file, which allows for straight-forward comparison
!> of model results with observations.
!--------------------------------------------------------------------------------------------------!
 MODULE virtual_measurement_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               exner,                                                                              &
               heatflux_output_conversion,                                                         &
               hyp,                                                                                &
               q,                                                                                  &
               ql,                                                                                 &
               pt,                                                                                 &
               rho_air,                                                                            &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               waterflux_output_conversion,                                                        &
               zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  convert_utm_to_geographic,                                                          &
               degc_to_k,                                                                          &
               magnus,                                                                             &
               pi,                                                                                 &
               rd_d_rv

    USE biometeorology_mod,                                                                        &
        ONLY:  bio_vm_sampling

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    USE chem_modules,                                                                              &
        ONLY:  chem_species

    USE control_parameters,                                                                        &
        ONLY:  air_chemistry,                                                                      &
               biometeorology,                                                                     &
               coupling_char,                                                                      &
               debug_output,                                                                       &
               dz,                                                                                 &
               end_time,                                                                           &
               humidity,                                                                           &
               land_surface,                                                                       &
               message_string,                                                                     &
               neutral,                                                                            &
               origin_date_time,                                                                   &
               rho_surface,                                                                        &
               surface_pressure,                                                                   &
               time_since_reference_point,                                                         &
               urban_surface,                                                                      &
               virtual_measurement

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE data_output_module

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               topo_top_ind,                                                                       &
               topo_flags

    USE kinds

    USE land_surface_model_mod,                                                                    &
        ONLY:  lsm_vm_sampling

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               coord_ref_sys,                                                                      &
               crs_list,                                                                           &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               init_model,                                                                         &
               input_file_atts,                                                                    &
               input_file_vm,                                                                      &
               input_pids_static,                                                                  &
               input_pids_vm,                                                                      &
               inquire_fill_value,                                                                 &
               open_read_file,                                                                     &
               pids_id

    USE pegrid

    USE radiation_model_mod,                                                                       &
        ONLY:  radiation,                                                                          &
               radiation_vm_sampling

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    USE urban_surface_mod,                                                                         &
        ONLY:  usm_vm_sampling

    IMPLICIT NONE

    TYPE virt_general
       INTEGER(iwp) ::  nvm = 0  !< number of virtual measurements

       REAL(wp) ::  origin_x_pd = HUGE( -1.0_wp ) !< EUTM coordinate of the western PALM domain boundary
       REAL(wp) ::  origin_y_pd = HUGE( -1.0_wp ) !< NUTM coordinate of the southern PALM domain boundary
    END TYPE virt_general

    TYPE virt_var_atts
       CHARACTER(LEN=100) ::  coordinates           !< defined longname of the variable
       CHARACTER(LEN=100) ::  grid_mapping          !< defined longname of the variable
       CHARACTER(LEN=100) ::  long_name             !< defined longname of the variable
       CHARACTER(LEN=100) ::  name                  !< variable name
       CHARACTER(LEN=100) ::  standard_name         !< defined standard name of the variable
       CHARACTER(LEN=100) ::  units                 !< unit of the output variable

       LOGICAL            ::  sampled               !< flag indicating whether a variable has been sampled

       REAL(wp)           ::  fill_value = -9999.0  !< _FillValue attribute
    END TYPE virt_var_atts

    TYPE virt_mea
       CHARACTER(LEN=100) ::  feature_type                      !< type of the real-world measurement
       CHARACTER(LEN=100) ::  feature_type_out = 'timeSeries'   !< type of the virtual measurement
                                                                 !< (all will be timeSeries, even trajectories)
       CHARACTER(LEN=100) ::  nc_filename                       !< name of the NetCDF output file for the station
       CHARACTER(LEN=100) ::  site                              !< name of the measurement site

       CHARACTER(LEN=1000) ::  data_content = REPEAT(' ', 1000)  !< string of measured variables (data output only)

       INTEGER(iwp) ::  end_coord_a     = 0  !< end coordinate in NetCDF file for local atmosphere observations
       INTEGER(iwp) ::  end_coord_s     = 0  !< end coordinate in NetCDF file for local soil observations
       INTEGER(iwp) ::  file_time_index = 0  !< time index in NetCDF output file
       INTEGER(iwp) ::  ns              = 0  !< number of observation coordinates on subdomain, for atmospheric measurements
       INTEGER(iwp) ::  ns_tot          = 0  !< total number of observation coordinates, for atmospheric measurements
       INTEGER(iwp) ::  nst                  !< number of coordinate points given for a measurement
       INTEGER(iwp) ::  nmeas                !< number of measured variables (atmosphere + soil)
       INTEGER(iwp) ::  ns_soil         = 0  !< number of observation coordinates on subdomain, for soil measurements
       INTEGER(iwp) ::  ns_soil_tot     = 0  !< total number of observation coordinates, for soil measurements
       INTEGER(iwp) ::  start_coord_a   = 0  !< start coordinate in NetCDF file for local atmosphere observations
       INTEGER(iwp) ::  start_coord_s   = 0  !< start coordinate in NetCDF file for local soil observations

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  dim_t  !< number observations individual for each trajectory
                                                          !< or station that are no _FillValues

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i       !< grid index for measurement position in x-direction
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j       !< grid index for measurement position in y-direction
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k       !< grid index for measurement position in k-direction

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_soil  !< grid index for measurement position in x-direction
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_soil  !< grid index for measurement position in y-direction
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_soil  !< grid index for measurement position in k-direction

       LOGICAL ::  interpolate         = .FALSE.  !< flag indicating that interpolation onto sampling location is desired
       LOGICAL ::  soil_sampling       = .FALSE.  !< flag indicating that soil state variables are present
       LOGICAL ::  soil_var_out        = .FALSE.  !< flag indicating that soil state variables were actually sampled
       LOGICAL ::  trajectory          = .FALSE.  !< flag indicating that the observation is a mobile observation
       LOGICAL ::  timeseries          = .FALSE.  !< flag indicating that the observation is a point measurement fixed in space
       LOGICAL ::  timeseries_profile  = .FALSE.  !< flag indicating that the observation is a stationary profile measurement

       REAL(wp) ::  fill_eutm          !< fill value for UTM coordinates in case of missing values
       REAL(wp) ::  fill_nutm          !< fill value for UTM coordinates in case of missing values
       REAL(wp) ::  fill_zar           !< fill value for heigth coordinates in case of missing values
       REAL(wp) ::  fillout = -9999.0  !< fill value for output in case an observation is taken e.g. from inside a building
       REAL(wp) ::  origin_x_obs       !< origin of the observation in UTM coordiates in x-direction
       REAL(wp) ::  origin_y_obs       !< origin of the observation in UTM coordiates in y-direction

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coords  !< exact sampling location
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  depth   !< measurement depth in soil
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zar     !< measurement height above ground level

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  measured_vars       !< measured variables
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  measured_vars_soil  !< measured variables

       TYPE( virt_var_atts ), DIMENSION(:), ALLOCATABLE ::  var_atts !< variable attributes
    END TYPE virt_mea

    CHARACTER(LEN=5)  ::  char_eutm = 'E_UTM'            !< dimension name for UTM coordinate easting
    CHARACTER(LEN=11) ::  char_feature = 'featureType'   !< attribute name for feature type

    ! This need to be generalized
    CHARACTER(LEN=10)  ::  char_fill = '_FillValue'                 !< attribute name for fill value
    CHARACTER(LEN=6)   ::  char_height = 'height'                   !< attribute name indicating height above surface (for trajectories only)
    CHARACTER(LEN=9)   ::  char_long = 'long_name'                  !< attribute name for long_name
    CHARACTER(LEN=18)  ::  char_mv = 'measured_variables'           !< variable name for the array with the measured variable names
    CHARACTER(LEN=5)   ::  char_nutm = 'N_UTM'                      !< dimension name for UTM coordinate northing
    CHARACTER(LEN=18)  ::  char_numstations = 'number_of_stations'  !< attribute name for number of stations
    CHARACTER(LEN=8)   ::  char_origx = 'origin_x'                  !< attribute name for station coordinate in x
    CHARACTER(LEN=8)   ::  char_origy = 'origin_y'                  !< attribute name for station coordinate in y
    CHARACTER(LEN=4)   ::  char_site = 'site'                       !< attribute name for site name
    CHARACTER(LEN=11)  ::  char_soil = 'soil_sample'                !< attribute name for soil sampling indication
    CHARACTER(LEN=13)  ::  char_standard = 'standard_name'          !< attribute name for standard_name
    CHARACTER(LEN=9)   ::  char_station_h = 'station_h'             !< variable name indicating height of the site
    CHARACTER(LEN=5)   ::  char_unit = 'units'                      !< attribute name for standard_name
    CHARACTER(LEN=1)   ::  char_zar = 'z'                           !< attribute name indicating height above reference level
    CHARACTER(LEN=800) ::  dom_error_message                        !< error message returned by the data-output module
    CHARACTER(LEN=10)  ::  type_ts   = 'timeSeries'                 !< name of stationary point measurements
    CHARACTER(LEN=10)  ::  type_traj = 'trajectory'                 !< name of line measurements
    CHARACTER(LEN=17)  ::  type_tspr = 'timeSeriesProfile'          !< name of stationary profile measurements

    CHARACTER(LEN=6), DIMENSION(1:5) ::  soil_vars = (/ 't_soil',                                  & !< list of soil variables
                                                        'm_soil',                                  &
                                                        'lwc   ',                                  &
                                                        'lwcs  ',                                  &
                                                        'smp   ' /)

    CHARACTER(LEN=10), DIMENSION(0:1,1:9) ::  chem_vars = RESHAPE( (/ 'mcpm1     ', 'PM1       ',  &
                                                                      'mcpm2p5   ', 'PM2.5     ',  &
                                                                      'mcpm10    ', 'PM10      ',  &
                                                                      'mfno2     ', 'NO2       ',  &
                                                                      'mfno      ', 'NO        ',  &
                                                                      'mcno2     ', 'NO2       ',  &
                                                                      'mcno      ', 'NO        ',  &
                                                                      'tro3      ', 'O3        ',  &
                                                                      'ncaa      ', 'PM10      '   & !simply assume ncaa to be PM10
                                                                    /), (/ 2, 9 /) )

    INTEGER(iwp) ::  maximum_name_length = 32  !< maximum name length of station names
    INTEGER(iwp) ::  ntimesteps                !< number of timesteps defined in NetCDF output file
    INTEGER(iwp) ::  off_pr              = 0   !< number of neighboring grid points (in horizontal direction) where virtual profile
                                               !< measurements shall be taken, in addition to the given coordinates in the driver
    INTEGER(iwp) ::  off_pr_z            = 0   !< number of additional grid points (in each upwardd and downward direction) where
                                               !< virtual profile measurements shall be taken, in addition to the given z coordinates in the driver
    INTEGER(iwp) ::  off_ts              = 0   !< number of neighboring grid points (in horizontal direction) where virtual profile
                                               !< measurements shall be taken, in addition to the given coordinates in the driver
    INTEGER(iwp) ::  off_ts_z            = 0   !< number of additional grid points (in each upwardd and downward direction) where
                                               !< virtual profile measurements shall be taken, in addition to the given z coordinates in the driver
    INTEGER(iwp) ::  off_tr              = 0   !< number of neighboring grid points (in horizontal direction) where virtual profile
                                               !< measurements shall be taken, in addition to the given coordinates in the driver
    INTEGER(iwp) ::  off_tr_z            = 0   !< number of additional grid points (in each upward and downward direction) where
                                               !< virtual profile measurements shall be taken, in addition to the given z coordinates in the driver
    LOGICAL ::  interpolate_to_location      = .FALSE.  !< flag to control if either nearest neighbor or interpolation is used
    LOGICAL ::  global_attribute             = .TRUE.   !< flag indicating a global attribute
    LOGICAL ::  initial_write_coordinates    = .FALSE.  !< flag indicating a global attribute
    LOGICAL ::  warn_out_of_domain_locations = .FALSE.  !< flag to enable warnings

    REAL(wp) ::  dt_virtual_measurement      = -HUGE( 0.0_wp )  !< overall sampling interval
    REAL(wp) ::  dt_virtual_measurement_pr   = -HUGE( 0.0_wp )  !< sampling interval escpecially for profile meausrements
    REAL(wp) ::  dt_virtual_measurement_ts   = -HUGE( 0.0_wp )  !< sampling interval escpecially for point meausrements
    REAL(wp) ::  dt_virtual_measurement_tr   = -HUGE( 0.0_wp )  !< sampling interval escpecially for trajectory meausrements
    REAL(wp) ::  time_virtual_measurement_pr = 0.0_wp           !< time since last sampling of profile data
    REAL(wp) ::  time_virtual_measurement_ts = 0.0_wp           !< time since last sampling of point data
    REAL(wp) ::  time_virtual_measurement_tr = 0.0_wp           !< time since last sampling of trajectory data
    REAL(wp) ::  vm_time_start               = 0.0_wp           !< time after which sampling shall start

    TYPE( virt_general )                        ::  vmea_general  !< data structure which encompasses global variables
    TYPE( virt_mea ), DIMENSION(:), ALLOCATABLE ::  vmea          !< data structure containing station-specific variables

    INTERFACE vm_check_parameters
       MODULE PROCEDURE vm_check_parameters
    END INTERFACE vm_check_parameters

    INTERFACE vm_data_output
       MODULE PROCEDURE vm_data_output
    END INTERFACE vm_data_output

    INTERFACE vm_init
       MODULE PROCEDURE vm_init
    END INTERFACE vm_init

    INTERFACE vm_init_output
       MODULE PROCEDURE vm_init_output
    END INTERFACE vm_init_output

    INTERFACE vm_parin
       MODULE PROCEDURE vm_parin
    END INTERFACE vm_parin

    INTERFACE vm_sampling
       MODULE PROCEDURE vm_sampling
    END INTERFACE vm_sampling

    SAVE

    PRIVATE

!
!-- Public interfaces
    PUBLIC  vm_check_parameters,                                                                   &
            vm_data_output,                                                                        &
            vm_init,                                                                               &
            vm_init_output,                                                                        &
            vm_parin,                                                                              &
            vm_sampling

!
!-- Public variables
    PUBLIC  dt_virtual_measurement_pr,                                                             &
            dt_virtual_measurement_ts,                                                             &
            dt_virtual_measurement_tr,                                                             &
            time_virtual_measurement_pr,                                                           &
            time_virtual_measurement_ts,                                                           &
            time_virtual_measurement_tr,                                                           &
            vmea,                                                                                  &
            vmea_general,                                                                          &
            vm_time_start

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters for virtual measurement module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_check_parameters

    IF ( .NOT. virtual_measurement )  RETURN
!
!-- Virtual measurements require a setup file.
    IF ( .NOT. input_pids_vm )  THEN
       message_string = 'virtual measurements require a setup input a setup input file for ' //    &
                        'the site locations'
       CALL message( 'vm_check_parameters', 'VME0001', 1, 2, 0, 6, 0 )
    ENDIF

#if !defined( __netcdf4_parallel )
!
!-- In case of non-parallel NetCDF the virtual measurement output is not
!-- working. This is only designed for parallel NetCDF.
    message_string = 'virtual measurements require parallel netCDF'
    CALL message( 'vm_check_parameters', 'VME0002', 1, 2, 0, 6, 0 )
#endif
!
!-- Check if the given number of neighboring grid points do not exceed the number
!-- of ghost points.
    IF ( off_pr > nbgp - 1  .OR.  off_ts > nbgp - 1  .OR.  off_tr > nbgp - 1 )  THEN
       WRITE ( message_string,* )  'for virtual measurements the number of surrounding grid ' //   &
                                   'points must not be larger than the number of ghost points ' // &
                                   '- 1, which is: ', nbgp - 1
       CALL message( 'vm_check_parameters', 'VME0003', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set sampling intervals to general sampling interval if not set properly. Finally, check
!-- sampling intervals for consistency.
    IF ( dt_virtual_measurement_pr < 0.0_wp )  dt_virtual_measurement_pr = dt_virtual_measurement
    IF ( dt_virtual_measurement_ts < 0.0_wp )  dt_virtual_measurement_ts = dt_virtual_measurement
    IF ( dt_virtual_measurement_tr < 0.0_wp )  dt_virtual_measurement_tr = dt_virtual_measurement

    IF ( dt_virtual_measurement_pr <= 0.0  .OR.  dt_virtual_measurement_ts <= 0.0_wp  .OR.         &
         dt_virtual_measurement_tr <= 0.0 )                                                        &
    THEN
       message_string = 'dt_virtual_measurement(pr/ts/tr) must be > 0.0.'
       CALL message( 'check_parameters', 'VME0004', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE vm_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defines variable attributes according to UC2 standard. Note, later this list can be
!> moved to the data-output module where it can be re-used also for other output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_set_attributes( output_variable )

    TYPE( virt_var_atts ), INTENT(INOUT) ::  output_variable !< data structure with attributes that need to be set

    output_variable%long_name     = 'none'
    output_variable%standard_name = 'none'
    output_variable%units         = 'none'
    output_variable%coordinates   = 'lon lat E_UTM N_UTM x y z time station_name'
    output_variable%grid_mapping  = 'crs'

    SELECT CASE ( TRIM( output_variable%name ) )

       CASE ( 'u' )
          output_variable%long_name     = 'u wind component'
          output_variable%units         = 'm s-1'

       CASE ( 'ua' )
          output_variable%long_name     = 'eastward wind'
          output_variable%standard_name = 'eastward_wind'
          output_variable%units         = 'm s-1'

       CASE ( 'v' )
          output_variable%long_name     = 'v wind component'
          output_variable%units         = 'm s-1'

       CASE ( 'va' )
          output_variable%long_name     = 'northward wind'
          output_variable%standard_name = 'northward_wind'
          output_variable%units         = 'm s-1'

       CASE ( 'w' )
          output_variable%long_name     = 'w wind component'
          output_variable%standard_name = 'upward_air_velocity'
          output_variable%units         = 'm s-1'

       CASE ( 'wspeed' )
          output_variable%long_name     = 'wind speed'
          output_variable%standard_name = 'wind_speed'
          output_variable%units         = 'm s-1'

       CASE ( 'wdir' )
          output_variable%long_name     = 'wind from direction'
          output_variable%standard_name = 'wind_from_direction'
          output_variable%units         = 'degrees'

       CASE ( 'theta' )
          output_variable%long_name     = 'air potential temperature'
          output_variable%standard_name = 'air_potential_temperature'
          output_variable%units         = 'K'

       CASE ( 'utheta' )
          output_variable%long_name     = 'eastward kinematic sensible heat flux in air'
          output_variable%units         = 'K m s-1'

       CASE ( 'vtheta' )
          output_variable%long_name     = 'northward kinematic sensible heat flux in air'
          output_variable%units         = 'K m s-1'

       CASE ( 'wtheta' )
          output_variable%long_name     = 'upward kinematic sensible heat flux in air'
          output_variable%units         = 'K m s-1'

       CASE ( 'ta' )
          output_variable%long_name     = 'air temperature'
          output_variable%standard_name = 'air_temperature'
          output_variable%units         = 'degree_C'

       CASE ( 't_va' )
          output_variable%long_name     = 'virtual acoustic temperature'
          output_variable%units         = 'K'

       CASE ( 't_mrt' )
          output_variable%long_name     = 'mean radiant temperature'
          output_variable%units         = 'degree_C'

       CASE ( 't_perceived' )
          output_variable%long_name     = 'perceived temperature'
          output_variable%units         = 'degree_C'

       CASE ( 't_pet' )
          output_variable%long_name     = 'physiological equivalent temperature'
          output_variable%units         = 'degree_C'

       CASE ( 't_utci' )
          output_variable%long_name     = 'universal thermal climate index'
          output_variable%units         = 'degree_C'

       CASE ( 'haa' )
          output_variable%long_name     = 'absolute atmospheric humidity'
          output_variable%units         = 'kg m-3'

       CASE ( 'hus' )
          output_variable%long_name     = 'specific humidity'
          output_variable%standard_name = 'specific_humidity'
          output_variable%units         = 'kg kg-1'

       CASE ( 'hur' )
          output_variable%long_name     = 'relative humidity'
          output_variable%standard_name = 'relative_humidity'
          output_variable%units         = '1'

       CASE ( 'rlu' )
          output_variable%long_name     = 'upwelling longwave flux in air'
          output_variable%standard_name = 'upwelling_longwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rlus' )
          output_variable%long_name     = 'surface upwelling longwave flux in air'
          output_variable%standard_name = 'surface_upwelling_longwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rld' )
          output_variable%long_name     = 'downwelling longwave flux in air'
          output_variable%standard_name = 'downwelling_longwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rsddif' )
          output_variable%long_name     = 'diffuse downwelling shortwave flux in air'
          output_variable%standard_name = 'diffuse_downwelling_shortwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rsd' )
          output_variable%long_name     = 'downwelling shortwave flux in air'
          output_variable%standard_name = 'downwelling_shortwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rnds' )
          output_variable%long_name     = 'surface net downward radiative flux'
          output_variable%standard_name = 'surface_net_downward_radiative_flux'
          output_variable%units         = 'W m-2'

       CASE ( 'rsu' )
          output_variable%long_name     = 'upwelling shortwave flux in air'
          output_variable%standard_name = 'upwelling_shortwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rsus' )
          output_variable%long_name     = 'surface upwelling shortwave flux in air'
          output_variable%standard_name = 'surface_upwelling_shortwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'rsds' )
          output_variable%long_name     = 'surface downwelling shortwave flux in air'
          output_variable%standard_name = 'surface_downwelling_shortwave_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'hfss' )
          output_variable%long_name     = 'surface upward sensible heat flux'
          output_variable%standard_name = 'surface_upward_sensible_heat_flux'
          output_variable%units         = 'W m-2'

       CASE ( 'hfls' )
          output_variable%long_name     = 'surface upward latent heat flux'
          output_variable%standard_name = 'surface_upward_latent_heat_flux'
          output_variable%units         = 'W m-2'

       CASE ( 'ts' )
          output_variable%long_name     = 'surface temperature'
          output_variable%standard_name = 'surface_temperature'
          output_variable%units         = 'K'

       CASE ( 'thetas' )
          output_variable%long_name     = 'surface layer temperature scale'
          output_variable%units         = 'K'

       CASE ( 'us' )
          output_variable%long_name     = 'friction velocity'
          output_variable%units         = 'm s-1'

       CASE ( 'uw' )
          output_variable%long_name     = 'upward eastward kinematic momentum flux in air'
          output_variable%units         = 'm2 s-2'

       CASE ( 'vw' )
          output_variable%long_name     = 'upward northward kinematic momentum flux in air'
          output_variable%units         = 'm2 s-2'

       CASE ( 'uv' )
          output_variable%long_name     = 'eastward northward kinematic momentum flux in air'
          output_variable%units         = 'm2 s-2'

       CASE ( 'm_soil' )
          output_variable%long_name     = 'soil moisture volumetric'
          output_variable%units         = 'm3 m-3'

       CASE ( 't_soil' )
          output_variable%long_name     = 'soil temperature'
          output_variable%standard_name = 'soil_temperature'
          output_variable%units         = 'degree_C'

       CASE ( 'hfdg' )
          output_variable%long_name     = 'downward heat flux at ground level in soil'
          output_variable%standard_name = 'downward_heat_flux_at_ground_level_in_soil'
          output_variable%units         = 'W m-2'

       CASE ( 'hfds' )
          output_variable%long_name     = 'downward heat flux in soil'
          output_variable%standard_name = 'downward_heat_flux_in_soil'
          output_variable%units         = 'W m-2'

       CASE ( 'hfla' )
          output_variable%long_name     = 'upward latent heat flux in air'
          output_variable%standard_name = 'upward_latent_heat_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'hfsa' )
          output_variable%long_name     = 'upward latent heat flux in air'
          output_variable%standard_name = 'upward_sensible_heat_flux_in_air'
          output_variable%units         = 'W m-2'

       CASE ( 'jno2' )
          output_variable%long_name     = 'photolysis rate of nitrogen dioxide'
          output_variable%standard_name = 'photolysis_rate_of_nitrogen_dioxide'
          output_variable%units         = 's-1'

       CASE ( 'lwcs' )
          output_variable%long_name     = 'liquid water content of soil layer'
          output_variable%standard_name = 'liquid_water_content_of_soil_layer'
          output_variable%units         = 'kg m-2'

       CASE ( 'lwp' )
          output_variable%long_name     = 'liquid water path'
          output_variable%standard_name = 'atmosphere_mass_content_of_cloud_liquid_water'
          output_variable%units         = 'kg m-2'

       CASE ( 'ps' )
          output_variable%long_name     = 'surface air pressure'
          output_variable%standard_name = 'surface_air_pressure'
          output_variable%units         = 'hPa'

       CASE ( 'pwv' )
          output_variable%long_name     = 'water vapor partial pressure in air'
          output_variable%standard_name = 'water_vapor_partial_pressure_in_air'
          output_variable%units         = 'hPa'

       CASE ( 't_lw' )
          output_variable%long_name     = 'land water temperature'
          output_variable%units         = 'degree_C'

       CASE ( 'tb' )
          output_variable%long_name     = 'brightness temperature'
          output_variable%standard_name = 'brightness_temperature'
          output_variable%units         = 'K'

       CASE ( 'uqv' )
          output_variable%long_name     = 'eastward kinematic latent heat flux in air'
          output_variable%units         = 'g kg-1 m s-1'

       CASE ( 'vqv' )
          output_variable%long_name     = 'northward kinematic latent heat flux in air'
          output_variable%units         = 'g kg-1 m s-1'

       CASE ( 'wqv' )
          output_variable%long_name     = 'upward kinematic latent heat flux in air'
          output_variable%units         = 'g kg-1 m s-1'

       CASE ( 'mcpm1' )
          output_variable%long_name     = 'mass concentration of pm1 ambient aerosol particles in air'
          output_variable%standard_name = 'mass_concentration_of_pm1_ambient_aerosol_particles_in_air'
          output_variable%units         = 'kg m-3'

       CASE ( 'mcpm10' )
          output_variable%long_name     = 'mass concentration of pm10 ambient aerosol particles in air'
          output_variable%standard_name = 'mass_concentration_of_pm10_ambient_aerosol_particles_in_air'
          output_variable%units         = 'kg m-3'

       CASE ( 'mcpm2p5' )
          output_variable%long_name     = 'mass concentration of pm2p5 ambient aerosol particles in air'
          output_variable%standard_name = 'mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air'
          output_variable%units         = 'kg m-3'

       CASE ( 'mfno', 'mcno'  )
          output_variable%long_name     = 'mole fraction of nitrogen monoxide in air'
          output_variable%standard_name = 'mole_fraction_of_nitrogen_monoxide_in_air'
          output_variable%units         = 'ppm' !'mol mol-1'

       CASE ( 'mfno2', 'mcno2'  )
          output_variable%long_name     = 'mole fraction of nitrogen dioxide in air'
          output_variable%standard_name = 'mole_fraction_of_nitrogen_dioxide_in_air'
          output_variable%units         = 'ppm' !'mol mol-1'

       CASE ( 'ncaa' )
          output_variable%long_name     = 'number concentration of ambient aerosol particles in air'
          output_variable%standard_name = 'number_concentration_of_ambient_aerosol_particles_in_air'
          output_variable%units         = 'm-3' !'mol mol-1'

       CASE ( 'tro3'  )
          output_variable%long_name     = 'mole fraction of ozone in air'
          output_variable%standard_name = 'mole_fraction_of_ozone_in_air'
          output_variable%units         = 'ppm' !'mol mol-1'

       CASE DEFAULT

    END SELECT

 END SUBROUTINE vm_set_attributes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read namelist for the virtual measurement module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /virtual_measurement_parameters/  dt_virtual_measurement,                             &
                                               dt_virtual_measurement_pr,                          &
                                               dt_virtual_measurement_ts,                          &
                                               dt_virtual_measurement_tr,                          &
                                               interpolate_to_location,                            &
                                               off_pr,                                             &
                                               off_pr_z,                                           &
                                               off_tr,                                             &
                                               off_tr_z,                                           &
                                               off_ts,                                             &
                                               off_ts_z,                                           &
                                               switch_off_module,                                  &
                                               vm_time_start,                                      &
                                               warn_out_of_domain_locations

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, virtual_measurement_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    virtual_measurements_parameters namelist was found and read correctly. Enable the
!--    virtual measurements.
       IF ( .NOT. switch_off_module )  virtual_measurement = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    virtual_measurement_parameters namelist was found but contained errors. Print an error
!--    message including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'virtual_measurement_parameters', line )

    ENDIF

 END SUBROUTINE vm_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize virtual measurements: read coordiante arrays and measured variables, set indicies
!> indicating the measurement points, read further attributes, etc..
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_init

    CHARACTER(LEN=5)                  ::  dum                           !< dummy string indicating station id
    CHARACTER(LEN=100), DIMENSION(200) ::  measured_variables_file = ''  !< array with all measured variables read from NetCDF
    CHARACTER(LEN=100), DIMENSION(200) ::  measured_variables      = ''  !< dummy array with all measured variables that are allowed

    INTEGER(iwp) ::  i          !< grid index of virtual observation point in x-direction
    INTEGER(iwp) ::  is         !< grid index of real observation point of the respective station in x-direction
    INTEGER(iwp) ::  j          !< grid index of observation point in x-direction
    INTEGER(iwp) ::  js         !< grid index of real observation point of the respective station in y-direction
    INTEGER(iwp) ::  k          !< grid index of observation point in x-direction
    INTEGER(iwp) ::  kl         !< lower vertical index of surrounding grid points of an observation coordinate
    INTEGER(iwp) ::  ks         !< grid index of real observation point of the respective station in z-direction
    INTEGER(iwp) ::  ksurf      !< topography top index
    INTEGER(iwp) ::  ku         !< upper vertical index of surrounding grid points of an observation coordinate
    INTEGER(iwp) ::  l          !< running index over all stations
    INTEGER(iwp) ::  len_char   !< character length of single measured variables without Null character
    INTEGER(iwp) ::  ll         !< running index over all measured variables in file
    INTEGER(iwp) ::  m          !< running index for surface elements
#if defined( __netcdf4_parallel )
    INTEGER(iwp) ::  n          !< running index over all cores
#endif
    INTEGER(iwp) ::  nofill     !< dummy for nofill return value (not used)
    INTEGER(iwp) ::  ns         !< counter variable for number of observation points on subdomain
    INTEGER(iwp) ::  off        !< number of horizontally surrounding grid points to be sampled
    INTEGER(iwp) ::  off_z      !< number of vertically surrounding grid points to be sampled
    INTEGER(iwp) ::  t          !< running index over number of trajectories

    INTEGER(ibp)                                ::  soil_dum  !< dummy variable to input a soil flag

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE     ::  ns_all  !< dummy array used to sum-up the number of observation coordinates

#if defined( __netcdf4_parallel )
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  ns_atmos  !< number of observation points for each station on each mpi rank
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  ns_soil   !< number of observation points for each station on each mpi rank
#endif

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  meas_flag  !< mask array indicating measurement positions

    LOGICAL  ::  on_pe  !< flag indicating that the respective measurement coordinate is on subdomain

    REAL(wp) ::  fill_height !< _FillValue for height coordinate (for trajectories)
    REAL(wp) ::  fill_eutm   !< _FillValue for coordinate array E_UTM
    REAL(wp) ::  fill_nutm   !< _FillValue for coordinate array N_UTM
    REAL(wp) ::  fill_zar    !< _FillValue for zar coordinate

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e_utm      !< easting UTM coordinate, temporary variable
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e_utm_tmp  !< EUTM coordinate before rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  height     !< observation height above ground (for trajectories)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  n_utm      !< northing UTM coordinate, temporary variable
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  n_utm_tmp  !< NUTM coordinate before rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  station_h  !< station height above reference
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zar        !< observation height above reference

    IF ( debug_output )  CALL debug_message( 'vm_init', 'start' )
#if defined( __netcdf )
!
!-- Open the input file.
    CALL open_read_file( TRIM( input_file_vm ) // TRIM( coupling_char ), pids_id )
!
!-- Obtain number of sites.
    CALL get_attribute( pids_id, char_numstations, vmea_general%nvm, global_attribute )
!
!-- Allocate data structure which encompasses all required information, such as  grid points indicies,
!-- absolute UTM coordinates, the measured quantities, etc. .
    ALLOCATE( vmea(1:vmea_general%nvm) )
!
!-- Allocate flag array. This dummy array is used to identify grid points where virtual measurements
!-- should be taken. In order to include also the surrounding grid points of the original coordinate,
!-- ghost points are required.
    ALLOCATE( meas_flag(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    meas_flag = 0
!
!-- Determine reference coordinate of the lower-left corner of the model domain. This will be used
!-- to determine the discrete sampling locations. Either this coordinate is taken from the
!-- static input file, which is the preference, or it is read from the measurement input file. In
!-- case it is not provided anyway, an error message will be given.
!-- Take it from the static input file.
    IF ( input_pids_static )  THEN
      vmea_general%origin_x_pd = init_model%origin_x
      vmea_general%origin_y_pd = init_model%origin_y
!
!-- Read from measurement file.
    ELSE
       CALL get_attribute( pids_id, char_origx, vmea_general%origin_x_pd, global_attribute,        &
                           ignore_error = .FALSE. )
       CALL get_attribute( pids_id, char_origy, vmea_general%origin_y_pd, global_attribute,        &
                           ignore_error = .FALSE. )
    ENDIF

!
!-- Check if the reference coordinate is given.
    IF ( vmea_general%origin_x_pd < 0.0_wp  .OR.  vmea_general%origin_y_pd < 0.0_wp )  THEN
       message_string = "no reference coordinate for the lower-left model domain corner are " //   &
                        "provided, neither in the static driver nor in the measurement input file"
       CALL message( 'vm_init', 'VME0005', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Loop over all sites in the setup file.
    DO  l = 1, vmea_general%nvm
!
!--    Determine suffix which contains the ID, ordered according to the number of measurements.
       IF( l < 10 )  THEN
          WRITE( dum, '(I1)')  l
       ELSEIF( l < 100 )  THEN
          WRITE( dum, '(I2)')  l
       ELSEIF( l < 1000 )  THEN
          WRITE( dum, '(I3)')  l
       ELSEIF( l < 10000 )  THEN
          WRITE( dum, '(I4)')  l
       ELSEIF( l < 100000 )  THEN
          WRITE( dum, '(I5)')  l
       ENDIF
!
!--    Read the origin site coordinates (UTM).
       CALL get_attribute( pids_id, char_origx // TRIM( dum ), vmea(l)%origin_x_obs, global_attribute )
       CALL get_attribute( pids_id, char_origy // TRIM( dum ), vmea(l)%origin_y_obs, global_attribute )
!
!--    Read site name.
       CALL get_attribute( pids_id, char_site // TRIM( dum ), vmea(l)%site, global_attribute )
!
!--    Read a flag which indicates that also soil quantities are take at the respective site
!--    (is part of the virtual measurement driver).
       CALL get_attribute( pids_id, char_soil // TRIM( dum ), soil_dum, global_attribute )
!
!--    Set flag indicating soil-sampling.
       IF ( soil_dum == 1 )  vmea(l)%soil_sampling = .TRUE.
!
!--    Read type of the measurement (trajectory, profile, timeseries).
       CALL get_attribute( pids_id, char_feature // TRIM( dum ), vmea(l)%feature_type,             &
                           global_attribute )
!
!---   Set logicals depending on the type of the measurement
       IF ( INDEX( vmea(l)%feature_type, type_tspr     ) /= 0 )  THEN
          vmea(l)%timeseries_profile = .TRUE.
       ELSEIF ( INDEX( vmea(l)%feature_type, type_ts   ) /= 0 )  THEN
          vmea(l)%timeseries = .TRUE.
       ELSEIF ( INDEX( vmea(l)%feature_type, type_traj ) /= 0 )  THEN
          vmea(l)%trajectory = .TRUE.
!
!--    Give error message in case the type matches none of the pre-defined types.
       ELSE
          message_string = 'attribute feature_type = "' // TRIM( vmea(l)%feature_type ) //         &
                           '" is not allowed'
          CALL message( 'vm_init', 'VME0006', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read string with all measured variables at this site.
       measured_variables_file = ''
       CALL get_variable( pids_id, char_mv // TRIM( dum ), measured_variables_file )
!
!--    Count the number of measured variables.
!--    For some NetCDF internal reasons, characters end with a NULL, i.e. also empty
!--    characters contain a NULL. Therefore, check the strings for a NULL to get the correct
!--    character length in order to compare them with the list of allowed variables.
       vmea(l)%nmeas  = 1
       DO  ll = 1, SIZE( measured_variables_file )
          IF ( measured_variables_file(ll)(1:1) /= CHAR(0)  .AND.                                  &
               measured_variables_file(ll)(1:1) /= ' ')  THEN
!
!--          Obtain character length of the character
             len_char = 1
             DO  WHILE ( measured_variables_file(ll)(len_char:len_char) /= CHAR(0)  .AND.          &
                 measured_variables_file(ll)(len_char:len_char) /= ' ' )
                len_char = len_char + 1
             ENDDO
             len_char = len_char - 1

             measured_variables(vmea(l)%nmeas) = measured_variables_file(ll)(1:len_char)
             vmea(l)%nmeas = vmea(l)%nmeas + 1

          ENDIF
       ENDDO
       vmea(l)%nmeas = vmea(l)%nmeas - 1
!
!--    Allocate data-type array for the measured variables names and attributes at the respective
!--    site.
       ALLOCATE( vmea(l)%var_atts(1:vmea(l)%nmeas) )
!
!--    Store the variable names in a data structure, which assigns further attributes to this name.
!--    Further, for data output reasons, create a string of output variables, which will be written
!--    into the attribute data_content.
       DO  ll = 1, vmea(l)%nmeas
          vmea(l)%var_atts(ll)%name = TRIM( measured_variables(ll) )

!           vmea(l)%data_content = TRIM( vmea(l)%data_content ) // " " //                            &
!                                  TRIM( vmea(l)%var_atts(ll)%name )
       ENDDO
!
!--    Read all the UTM coordinates for the site. Based on the coordinates, define the grid-index
!--    space on each subdomain where virtual measurements should be taken. The entire
!--    coordinate array (on the entire model domain) won't be stored as this would exceed memory
!--    requirements, particularly for trajectories.
       IF ( vmea(l)%nmeas > 0 )  THEN
!
!--       For stationary and mobile measurements UTM coordinates are just one value and its dimension
!--       is "station".
          CALL get_dimension_length( pids_id, vmea(l)%nst, "station" // TRIM( dum ) )
!
!--       Allocate temporary arrays for UTM and height coordinates. On file, UTM coordinates
!--       can be 1D or 2D variables
          ALLOCATE( e_utm(1:vmea(l)%nst)       )
          ALLOCATE( n_utm(1:vmea(l)%nst)       )
          ALLOCATE( station_h(1:vmea(l)%nst)   )
          ALLOCATE( zar(1:vmea(l)%nst)         )
          IF ( vmea(l)%trajectory )  ALLOCATE( height(1:vmea(l)%nst) )
          e_utm     = 0.0_wp
          n_utm     = 0.0_wp
          station_h = 0.0_wp
          zar       = 0.0_wp
          IF ( vmea(l)%trajectory )  height = 0.0_wp

          ALLOCATE( e_utm_tmp(1:vmea(l)%nst) )
          ALLOCATE( n_utm_tmp(1:vmea(l)%nst) )
!
!--       Read UTM and height coordinates for all trajectories and times. In case
!--       they contain any missing values, replace them with default _FillValues.
          CALL inquire_fill_value( pids_id, char_eutm   // TRIM( dum ), nofill, fill_eutm    )
          CALL inquire_fill_value( pids_id, char_nutm   // TRIM( dum ), nofill, fill_nutm    )
          CALL inquire_fill_value( pids_id, char_zar    // TRIM( dum ), nofill, fill_zar     )
          IF ( vmea(l)%trajectory )                                                                &
             CALL inquire_fill_value( pids_id, char_height // TRIM( dum ), nofill, fill_height  )
!
!--       Further line is just to avoid compiler warnings. nofill might be used in future.
          IF ( nofill == 0  .OR.  nofill /= 0 )  CONTINUE
!
!--       Read observation coordinates. For trajectories the observation height is stored directly
!--       in z, while for timeSeries it is stored in z - station_h, according to UC2-standard.
          IF ( vmea(l)%trajectory )  THEN
             CALL get_variable( pids_id, char_eutm      // TRIM( dum ), e_utm(:)     )
             CALL get_variable( pids_id, char_nutm      // TRIM( dum ), n_utm(:)     )
             CALL get_variable( pids_id, char_height    // TRIM( dum ), height(:)    )
          ELSE
             CALL get_variable( pids_id, char_eutm      // TRIM( dum ), e_utm(:)     )
             CALL get_variable( pids_id, char_nutm      // TRIM( dum ), n_utm(:)     )
             CALL get_variable( pids_id, char_station_h // TRIM( dum ), station_h(:) )
             CALL get_variable( pids_id, char_zar       // TRIM( dum ), zar(:)       )
          ENDIF

          e_utm = MERGE( e_utm, vmea(l)%fillout, e_utm /= fill_eutm )
          n_utm = MERGE( n_utm, vmea(l)%fillout, n_utm /= fill_nutm )
          zar   = MERGE( zar,   vmea(l)%fillout, zar   /= fill_zar  )
          IF ( vmea(l)%trajectory )                                                                &
             height = MERGE( height, vmea(l)%fillout, height /= fill_height )
!
!--       Compute observation height above ground. Note, for trajectory measurements the height
!--       above the surface is actually stored in 'height'.
          IF ( vmea(l)%trajectory )  THEN
             zar      = height
             fill_zar = fill_height
          ELSE
             zar = zar - station_h
          ENDIF
!
!--       Based on UTM coordinates, check if the measurement station or parts of the trajectory are
!--       on subdomain. This case, setup grid index space sample these quantities.
          meas_flag = 0
          DO  t = 1, vmea(l)%nst
!
!--          First, compute relative x- and y-coordinates with respect to the lower-left origin of
!--          the model domain, which is the difference between UTM coordinates. If the origin
!--          is not correct, the virtual sites will be misplaced. Further, in case of a rotated
!--          model domain, the UTM coordinates must also be rotated.
             e_utm_tmp(t) = e_utm(t) - vmea_general%origin_x_pd
             n_utm_tmp(t) = n_utm(t) - vmea_general%origin_y_pd
             e_utm(t) = COS( init_model%rotation_angle * pi / 180.0_wp ) * e_utm_tmp(t)            &
                      - SIN( init_model%rotation_angle * pi / 180.0_wp ) * n_utm_tmp(t)
             n_utm(t) = SIN( init_model%rotation_angle * pi / 180.0_wp ) * e_utm_tmp(t)            &
                      + COS( init_model%rotation_angle * pi / 180.0_wp ) * n_utm_tmp(t)
!
!--          Compute grid indices relative to origin and check if these are on the subdomain. Note,
!--          virtual measurements will be taken also at grid points surrounding the station, hence,
!--          check also for these grid points. The number of surrounding grid points is set
!--          according to the featureType.
             IF ( vmea(l)%timeseries_profile )  THEN
                off   = off_pr
                off_z = off_pr_z
             ELSEIF ( vmea(l)%timeseries )  THEN
                off   = off_ts
                off_z = off_ts_z
             ELSEIF ( vmea(l)%trajectory )  THEN
                off   = off_tr
                off_z = off_tr_z
             ENDIF

             is = INT( e_utm(t) * ddx, KIND = iwp )
             js = INT( n_utm(t) * ddy, KIND = iwp )

!
!--          For timeseries output, check if the sampling location lies within the model domain.
             IF ( warn_out_of_domain_locations  .AND.  vmea(l)%timeseries )  THEN
                IF ( is < 0  .OR.  is > nx+1  .OR.  js < 0  .OR.  js > ny+1  .OR.                  &
                     zar(t) < 0.0_wp  .OR.  zar(t) >= zw(nzt) )                                    &
                THEN
                   WRITE( message_string, * ) 'sampling coordinates of site "',                    &
                                              TRIM( vmea(l)%site ), '" do not lie within the ',    &
                                              'PALM domain'
                   CALL message( 'vm_init', 'VME0007', 0, 1, 0, 6, 0 )
                ENDIF
             ENDIF
!
!--          Is the observation point on subdomain?
             on_pe = ( is >= nxl  .AND.  is <= nxr  .AND.  js >= nys  .AND.  js <= nyn )
!
!--          Check if observation coordinate is on subdomain.
             IF ( on_pe )  THEN
!
!--             Determine vertical index which corresponds to the observation height.
                ksurf = topo_top_ind(js,is,0)
                ks = MINLOC( ABS( zu - zw(ksurf) - zar(t) ), DIM = 1 ) - 1
!
!--             Set mask array at the observation coordinates. Also, flag the surrounding
!--             coordinate points, but first check whether the surrounding coordinate points are
!--             on the subdomain.
                kl = MERGE( ks-off_z, ksurf, ks-off_z >= nzb  .AND. ks-off_z >= ksurf )
                ku = MERGE( ks+off_z, nzt,   ks+off_z < nzt+1 )

                DO  i = is-off, is+off
                   DO  j = js-off, js+off
                      DO  k = kl, ku
                         meas_flag(k,j,i) = MERGE( IBSET( meas_flag(k,j,i), 0 ), 0,                &
                                                   BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--          In case of timeseries output with one single point, exact interpolation onto the
!--          sampling location is possible. This case, store the respective relative coordinate.
!--          For trajectories and profiles interpolation onto the exact sampling location is
!--          not possible at the moment. If this is desired at a later point in time, the following
!--          block needs to be moved out of the loop over t.
             IF ( vmea(l)%timeseries  .AND.  interpolate_to_location  .AND.  off_ts == 0  .AND.    &
                  off_ts_z == 0 )                                                                  &
             THEN
                ALLOCATE( vmea(l)%coords(3) )

                IF ( on_pe )  THEN
                   vmea(l)%coords(1) = e_utm(t)
                   vmea(l)%coords(2) = n_utm(t)
                   vmea(l)%coords(3) = zar(t) + zw(topo_top_ind(js,is,0))

                   vmea(l)%interpolate = .TRUE.
                ENDIF
             ENDIF

          ENDDO
!
!--       Based on the flag array, count the number of sampling coordinates. Sampling coordinates in
!--       atmosphere and soil may be different, as within the soil all levels will be measured.
!--       Hence, count individually. Start with atmoshere.
          ns = 0
          DO  i = nxl-off, nxr+off
             DO  j = nys-off, nyn+off
                DO  k = nzb, nzt+1
                   ns = ns + MERGE( 1, 0, BTEST( meas_flag(k,j,i), 0 ) )
                ENDDO
             ENDDO
          ENDDO

!
!--       Store number of observation points on subdomain and allocate index arrays as well as array
!--       containing height information.
          vmea(l)%ns = ns

          ALLOCATE( vmea(l)%i(1:vmea(l)%ns) )
          ALLOCATE( vmea(l)%j(1:vmea(l)%ns) )
          ALLOCATE( vmea(l)%k(1:vmea(l)%ns) )
          ALLOCATE( vmea(l)%zar(1:vmea(l)%ns) )
!
!--       Based on the flag array store the grid indices which correspond to the observation
!--       coordinates.
          ns = 0
          DO  i = nxl-off, nxr+off
             DO  j = nys-off, nyn+off
                DO  k = nzb, nzt+1
                   IF ( BTEST( meas_flag(k,j,i), 0 ) )  THEN
                      ns = ns + 1
                      vmea(l)%i(ns) = i
                      vmea(l)%j(ns) = j
                      vmea(l)%k(ns) = k
                      vmea(l)%zar(ns) = zu(k) - zw(topo_top_ind(j,i,0))
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
!
!--       Same for the soil. Based on the flag array, count the number of sampling coordinates in
!--       soil. Sample at all soil levels in this case. Please note, soil variables can only be
!--       sampled on subdomains, not on ghost layers.
          IF ( vmea(l)%soil_sampling )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( ANY( BTEST( meas_flag(:,j,i), 0 ) ) )  THEN
                      IF ( surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )  THEN
                         vmea(l)%ns_soil = vmea(l)%ns_soil +                                       &
                                           surf_lsm%nzt_soil - surf_lsm%nzb_soil + 1
                      ENDIF
                      IF ( surf_usm%start_index(j,i) <= surf_usm%end_index(j,i) )  THEN
                         vmea(l)%ns_soil = vmea(l)%ns_soil +                                       &
                                           surf_usm%nzt_wall - surf_usm%nzb_wall + 1
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
!
!--       Allocate index arrays as well as array containing height information for soil.
          IF ( vmea(l)%soil_sampling )  THEN
             ALLOCATE( vmea(l)%i_soil(1:vmea(l)%ns_soil) )
             ALLOCATE( vmea(l)%j_soil(1:vmea(l)%ns_soil) )
             ALLOCATE( vmea(l)%k_soil(1:vmea(l)%ns_soil) )
             ALLOCATE( vmea(l)%depth(1:vmea(l)%ns_soil)  )
          ENDIF
!
!--       For soil, store the grid indices.
          ns = 0
          IF ( vmea(l)%soil_sampling )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( ANY( BTEST( meas_flag(:,j,i), 0 ) ) )  THEN
                      IF ( surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )  THEN
                         m = surf_lsm%start_index(j,i)
                         DO  k = surf_lsm%nzb_soil, surf_lsm%nzt_soil
                            ns = ns + 1
                            vmea(l)%i_soil(ns) = i
                            vmea(l)%j_soil(ns) = j
                            vmea(l)%k_soil(ns) = k
                            vmea(l)%depth(ns)  = -surf_lsm%zs(k)
                         ENDDO
                      ENDIF

                      IF ( surf_usm%start_index(j,i) <= surf_usm%end_index(j,i) )  THEN
                         m = surf_usm%start_index(j,i)
                         DO  k = surf_usm%nzb_wall, surf_usm%nzt_wall
                            ns = ns + 1
                            vmea(l)%i_soil(ns) = i
                            vmea(l)%j_soil(ns) = j
                            vmea(l)%k_soil(ns) = k
                            vmea(l)%depth(ns)  = -surf_usm%zw(k,m)
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
!
!--       Allocate array to save the sampled values.
          ALLOCATE( vmea(l)%measured_vars(1:vmea(l)%ns,1:vmea(l)%nmeas) )

          IF ( vmea(l)%soil_sampling )  THEN
             ALLOCATE( vmea(l)%measured_vars_soil(1:vmea(l)%ns_soil,1:vmea(l)%nmeas) )
          ENDIF
!
!--       Initialize with _FillValues.
          vmea(l)%measured_vars(1:vmea(l)%ns,1:vmea(l)%nmeas) = vmea(l)%fillout
          IF ( vmea(l)%soil_sampling )  THEN
             vmea(l)%measured_vars_soil(1:vmea(l)%ns_soil,1:vmea(l)%nmeas) = vmea(l)%fillout
          ENDIF
!
!--       Deallocate temporary coordinate arrays.
          IF ( ALLOCATED( e_utm )     )  DEALLOCATE( e_utm )
          IF ( ALLOCATED( n_utm )     )  DEALLOCATE( n_utm )
          IF ( ALLOCATED( e_utm_tmp ) )  DEALLOCATE( e_utm_tmp )
          IF ( ALLOCATED( n_utm_tmp ) )  DEALLOCATE( n_utm_tmp )
          IF ( ALLOCATED( n_utm )     )  DEALLOCATE( n_utm )
          IF ( ALLOCATED( zar  )      )  DEALLOCATE( zar  )
          IF ( ALLOCATED( height )    )  DEALLOCATE( height )
          IF ( ALLOCATED( station_h ) )  DEALLOCATE( station_h )

       ENDIF
    ENDDO
!
!-- Dellocate flag array
    DEALLOCATE( meas_flag )
!
!-- Close input file for virtual measurements.
    CALL close_input_file( pids_id )
!
!-- Sum-up the number of observation coordiates, for atmosphere first.
!-- This is actually only required for data output.
    ALLOCATE( ns_all(1:vmea_general%nvm) )
    ns_all = 0
#if defined( __parallel )
    CALL MPI_ALLREDUCE( vmea(:)%ns, ns_all(:), vmea_general%nvm, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    ns_all(:) = vmea(:)%ns
#endif
    vmea(:)%ns_tot = ns_all(:)
!
!-- Now for soil
    ns_all = 0
#if defined( __parallel )
    CALL MPI_ALLREDUCE( vmea(:)%ns_soil, ns_all(:), vmea_general%nvm, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    ns_all(:) = vmea(:)%ns_soil
#endif
    vmea(:)%ns_soil_tot = ns_all(:)

    DEALLOCATE( ns_all )
!
!-- In case of parallel NetCDF the start coordinate for each mpi rank needs to be defined, so that
!-- each processor knows where to write the data.
#if defined( __netcdf4_parallel )
    ALLOCATE( ns_atmos(0:numprocs-1,1:vmea_general%nvm) )
    ALLOCATE( ns_soil(0:numprocs-1,1:vmea_general%nvm)  )
    ns_atmos = 0
    ns_soil  = 0

    DO  l = 1, vmea_general%nvm
       ns_atmos(myid,l) = vmea(l)%ns
       ns_soil(myid,l)  = vmea(l)%ns_soil
    ENDDO

#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ns_atmos, numprocs * vmea_general%nvm,                       &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ns_soil, numprocs * vmea_general%nvm,                        &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    ns_atmos(0,:) = vmea(:)%ns
    ns_soil(0,:)  = vmea(:)%ns_soil
#endif

!
!-- Determine the start coordinate in NetCDF file for the local arrays. Start coordinates are
!-- initialized with zero for sake of simplicity in summation. However, in NetCDF the start
!-- coordinates must be >= 1, so that a one needs to be added at the end.
    DO  l = 1, vmea_general%nvm
       DO  n  = 0, myid - 1
          vmea(l)%start_coord_a = vmea(l)%start_coord_a + ns_atmos(n,l)
          vmea(l)%start_coord_s = vmea(l)%start_coord_s + ns_soil(n,l)
       ENDDO
!
!--    Start coordinate in NetCDF starts always at one not at 0.
       vmea(l)%start_coord_a = vmea(l)%start_coord_a + 1
       vmea(l)%start_coord_s = vmea(l)%start_coord_s + 1
!
!--    Determine the local end coordinate
       vmea(l)%end_coord_a = vmea(l)%start_coord_a + vmea(l)%ns - 1
       vmea(l)%end_coord_s = vmea(l)%start_coord_s + vmea(l)%ns_soil - 1
    ENDDO

    DEALLOCATE( ns_atmos )
    DEALLOCATE( ns_soil  )

#endif

#endif
    IF ( debug_output )  CALL debug_message( 'vm_init', 'end' )

 END SUBROUTINE vm_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize output using data-output module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_init_output

    CHARACTER(LEN=100) ::  variable_name  !< name of output variable

    INTEGER(iwp) ::  l             !< loop index
    INTEGER(iwp) ::  n             !< loop index
    INTEGER      ::  return_value  !< returned status value of called function

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ndim !< dummy to write dimension

    REAL(wp) ::  dt_vm    !< sampling interval of specific measurement
    REAL(wp) ::  dum_lat  !< transformed geographical coordinate (latitude)
    REAL(wp) ::  dum_lon  !< transformed geographical coordinate (longitude)

!
!-- Create directory where output files will be stored.
    CALL local_system( 'mkdir -p VM_OUTPUT' // TRIM( coupling_char ) )

!
!-- Loop over all sites.
    DO  l = 1, vmea_general%nvm
!
!--    Skip if no observations will be taken for this site.
       IF ( vmea(l)%ns_tot == 0  .AND.  vmea(l)%ns_soil_tot == 0 )  CYCLE
!
!--    Determine the number of output timesteps, which depend on the type of the measurement.
       IF ( vmea(l)%timeseries )  THEN
          dt_vm = dt_virtual_measurement_ts
       ELSEIF ( vmea(l)%timeseries_profile )  THEN
          dt_vm = dt_virtual_measurement_pr
       ELSE
          dt_vm = dt_virtual_measurement_tr
       ENDIF

       ntimesteps = CEILING( ( end_time - MAX( vm_time_start, time_since_reference_point )         &
                             ) / dt_vm )

!
!--    Define output file.
       WRITE( vmea(l)%nc_filename, '(A,I4.4)' ) 'VM_OUTPUT' // TRIM( coupling_char ) // '/' //     &
              'site', l

       return_value = dom_def_file( vmea(l)%nc_filename, 'netcdf4-parallel' )
!
!--    Define global attributes.
!--    Before, transform UTM into geographical coordinates.
       CALL convert_utm_to_geographic( crs_list, vmea(l)%origin_x_obs, vmea(l)%origin_y_obs,       &
                                       dum_lon, dum_lat )

       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'site',                   &
                                   value = TRIM( vmea(l)%site ) )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'title',                  &
                                   value = 'Virtual measurement output')
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'source',                 &
                                   value = 'PALM-4U')
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'institution',            &
                                   value = input_file_atts%institution )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'acronym',                &
                                   value = input_file_atts%acronym )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'author',                 &
                                   value = input_file_atts%author )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'contact_person',         &
                                   value = input_file_atts%author )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'iop',                    &
                                   value = input_file_atts%campaign )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'campaign',               &
                                   value = 'PALM-4U' )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'origin_time ',           &
                                   value = origin_date_time)
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'location',               &
                                   value = input_file_atts%location )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'origin_x',               &
                                   value = vmea(l)%origin_x_obs )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'origin_y',               &
                                   value = vmea(l)%origin_y_obs )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'origin_lon',             &
                                   value = dum_lon )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'origin_lat',             &
                                   value = dum_lat )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'origin_z',               &
                                   value = init_model%origin_z )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'rotation_angle',         &
                                   value = input_file_atts%rotation_angle )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'featureType',            &
                                   value = TRIM( vmea(l)%feature_type_out ) )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'data_content',           &
                                   value = TRIM( vmea(l)%data_content ) )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'creation_time',          &
                                   value = input_file_atts%creation_time )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'version', value = 1 ) !input_file_atts%version
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'Conventions',            &
                                   value = input_file_atts%conventions )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'dependencies',           &
                                   value = input_file_atts%dependencies )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'history',                &
                                   value = input_file_atts%history )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'references',             &
                                   value = input_file_atts%references )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'comment',                &
                                   value = input_file_atts%comment )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'keywords',               &
                                   value = input_file_atts%keywords )
       return_value = dom_def_att( vmea(l)%nc_filename, attribute_name = 'licence',                &
                                   value = '[UC]2 Open Licence; see [UC]2 ' //                     &
                                           'data policy available at ' //                          &
                                           'www.uc2-program.org/uc2_data_policy.pdf' )
!
!--    Define dimensions.
!--    station
       ALLOCATE( ndim(1:vmea(l)%ns_tot) )
       DO  n = 1, vmea(l)%ns_tot
          ndim(n) = n
       ENDDO
       return_value = dom_def_dim( vmea(l)%nc_filename, dimension_name = 'station',                &
                                   output_type = 'int32', bounds = (/1_iwp, vmea(l)%ns_tot/),      &
                                   values_int32 = ndim )
       DEALLOCATE( ndim )
!
!--    ntime
       ALLOCATE( ndim(1:ntimesteps) )
       DO  n = 1, ntimesteps
          ndim(n) = n
       ENDDO

       return_value = dom_def_dim( vmea(l)%nc_filename, dimension_name = 'ntime',                  &
                                   output_type = 'int32', bounds = (/1_iwp, ntimesteps/),          &
                                   values_int32 = ndim )
       DEALLOCATE( ndim )
!
!--    nv
       ALLOCATE( ndim(1:2) )
       DO  n = 1, 2
          ndim(n) = n
       ENDDO

       return_value = dom_def_dim( vmea(l)%nc_filename, dimension_name = 'nv',                     &
                                   output_type = 'int32', bounds = (/1_iwp, 2_iwp/),               &
                                   values_int32 = ndim )
       DEALLOCATE( ndim )
!
!--    maximum name length
       ALLOCATE( ndim(1:maximum_name_length) )
       DO  n = 1, maximum_name_length
          ndim(n) = n
       ENDDO

       return_value = dom_def_dim( vmea(l)%nc_filename, dimension_name = 'max_name_len',           &
                                   output_type = 'int32',                                          &
                                   bounds = (/1_iwp, maximum_name_length /), values_int32 = ndim )
       DEALLOCATE( ndim )
!
!--    Define coordinate variables.
!--    time
       variable_name = 'time'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/ 'station  ', 'ntime    '/),                &
                                   write_mode = 'independent', output_type = 'real32' )
!
!--    station_name
       variable_name = 'station_name'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/ 'max_name_len', 'station     ' /),         &
                                   write_mode = 'independent', output_type = 'char' )
!
!--    vrs (vertical reference system)
       variable_name = 'vrs'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/ 'station' /), write_mode = 'independent',  &
                                   output_type = 'int8' )
!
!--    crs (coordinate reference system)
       variable_name = 'crs'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/ 'station' /), write_mode = 'independent',  &
                                   output_type = 'int8' )
!
!--    z
       variable_name = 'z'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    station_h
       variable_name = 'station_h'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    x
       variable_name = 'x'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    y
       variable_name = 'y'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    E-UTM
       variable_name = 'E_UTM'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    N-UTM
       variable_name = 'N_UTM'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    latitude
       variable_name = 'lat'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    longitude
       variable_name = 'lon'
       return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,             &
                                   dimension_names = (/'station'/), write_mode = 'independent',    &
                                   output_type = 'real32' )
!
!--    Set attributes for the coordinate variables. Note, not all coordinates have the same number
!--    of attributes.
!--    Units
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time',                    &
                                   attribute_name = char_unit, value = 'seconds since ' //         &
                                   origin_date_time )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z',                       &
                                   attribute_name = char_unit, value = 'm' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_h',               &
                                   attribute_name = char_unit, value = 'm' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'x',                       &
                                   attribute_name = char_unit, value = 'm' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'y',                       &
                                   attribute_name = char_unit, value = 'm' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'E_UTM',                   &
                                   attribute_name = char_unit, value = 'm' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'N_UTM',                   &
                                   attribute_name = char_unit, value = 'm' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lat',                     &
                                   attribute_name = char_unit, value = 'degrees_north' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lon',                     &
                                   attribute_name = char_unit, value = 'degrees_east' )
!
!--    Long name
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_name',            &
                                   attribute_name = char_long, value = 'station name')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time',                    &
                                   attribute_name = char_long, value = 'time')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z',                       &
                                   attribute_name = char_long, value = 'height above origin' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_h',               &
                                   attribute_name = char_long, value = 'surface altitude' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'x',                       &
                                   attribute_name = char_long,                                     &
                                   value = 'distance to origin in x-direction')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'y',                       &
                                   attribute_name = char_long,                                     &
                                   value = 'distance to origin in y-direction')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'E_UTM',                   &
                                   attribute_name = char_long, value = 'easting' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'N_UTM',                   &
                                   attribute_name = char_long, value = 'northing' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lat',                     &
                                   attribute_name = char_long, value = 'latitude' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lon',                     &
                                   attribute_name = char_long, value = 'longitude' )
!
!--    Standard name
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_name',            &
                                   attribute_name = char_standard, value = 'platform_name')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time',                    &
                                   attribute_name = char_standard, value = 'time')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z',                       &
                                   attribute_name = char_standard,                                 &
                                   value = 'height_above_mean_sea_level' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_h',               &
                                   attribute_name = char_standard, value = 'surface_altitude' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'E_UTM',                   &
                                   attribute_name = char_standard,                                 &
                                   value = 'projection_x_coordinate' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'N_UTM',                   &
                                   attribute_name = char_standard,                                 &
                                   value = 'projection_y_coordinate' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lat',                     &
                                   attribute_name = char_standard, value = 'latitude' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lon',                     &
                                   attribute_name = char_standard, value = 'longitude' )
!
!--    Axis
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time',                    &
                                   attribute_name = 'axis', value = 'T')
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z',                       &
                                   attribute_name = 'axis', value = 'Z' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'x',                       &
                                   attribute_name = 'axis', value = 'X' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'y',                       &
                                   attribute_name = 'axis', value = 'Y' )
!
!--    Set further individual attributes for the coordinate variables.
!--    For station name
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_name',            &
                                   attribute_name = 'cf_role', value = 'timeseries_id' )
!
!--    For time
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time',                    &
                                   attribute_name = 'calendar', value = 'proleptic_gregorian' )
!        return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time',                    &
!                                    attribute_name = 'bounds', value = 'time_bounds' )
!
!--    For vertical reference system
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'vrs',                     &
                                   attribute_name = char_long, value = 'vertical reference system' )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'vrs',                     &
                                   attribute_name = 'system_name', value = 'DHHN2016' )
!
!--    For z
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z',                       &
                                   attribute_name = 'positive', value = 'up' )
!
!--    For coordinate reference system
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'epsg_code', value = coord_ref_sys%epsg_code )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'false_easting',                               &
                                   value = coord_ref_sys%false_easting )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'false_northing',                              &
                                   value = coord_ref_sys%false_northing )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'grid_mapping_name',                           &
                                   value = coord_ref_sys%grid_mapping_name )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'inverse_flattening',                          &
                                   value = coord_ref_sys%inverse_flattening )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'latitude_of_projection_origin',               &
                                   value = coord_ref_sys%latitude_of_projection_origin )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = char_long, value = coord_ref_sys%long_name )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'longitude_of_central_meridian',               &
                                   value = coord_ref_sys%longitude_of_central_meridian )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'longitude_of_prime_meridian',                 &
                                   value = coord_ref_sys%longitude_of_prime_meridian )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'scale_factor_at_central_meridian',            &
                                   value = coord_ref_sys%scale_factor_at_central_meridian )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = 'semi_major_axis',                             &
                                   value = coord_ref_sys%semi_major_axis )
       return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'crs',                     &
                                   attribute_name = char_unit, value = coord_ref_sys%units )
!
!--    In case of sampled soil quantities, define further dimensions and coordinates.
       IF ( vmea(l)%soil_sampling )  THEN
!
!--       station for soil
          ALLOCATE( ndim(1:vmea(l)%ns_soil_tot) )
          DO  n = 1, vmea(l)%ns_soil_tot
             ndim(n) = n
          ENDDO

          return_value = dom_def_dim( vmea(l)%nc_filename, dimension_name = 'station_soil',        &
                                      output_type = 'int32',                                       &
                                      bounds = (/1_iwp,vmea(l)%ns_soil_tot/), values_int32 = ndim )
          DEALLOCATE( ndim )
!
!--       ntime for soil
          ALLOCATE( ndim(1:ntimesteps) )
          DO  n = 1, ntimesteps
             ndim(n) = n
          ENDDO

          return_value = dom_def_dim( vmea(l)%nc_filename, dimension_name = 'ntime_soil',          &
                                      output_type = 'int32', bounds = (/1_iwp,ntimesteps/),        &
                                      values_int32 = ndim )
          DEALLOCATE( ndim )
!
!--       time for soil
          variable_name = 'time_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil', 'ntime_soil  '/),        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       station_name for soil
          variable_name = 'station_name_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/ 'max_name_len', 'station_soil' /),      &
                                      write_mode = 'independent', output_type = 'char' )
!
!--       z
          variable_name = 'z_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       station_h for soil
          variable_name = 'station_h_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent',  output_type = 'real32' )
!
!--       x soil
          variable_name = 'x_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!-        y soil
          variable_name = 'y_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       E-UTM soil
          variable_name = 'E_UTM_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       N-UTM soil
          variable_name = 'N_UTM_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       latitude soil
          variable_name = 'lat_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       longitude soil
          variable_name = 'lon_soil'
          return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      dimension_names = (/'station_soil'/),                        &
                                      write_mode = 'independent', output_type = 'real32' )
!
!--       Set attributes for the coordinate variables. Note, not all coordinates have the same
!--       number of attributes.
!--       Units
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time_soil',            &
                                      attribute_name = char_unit, value = 'seconds since ' //      &
                                      origin_date_time )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z_soil',               &
                                      attribute_name = char_unit, value = 'm' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_h_soil',       &
                                      attribute_name = char_unit, value = 'm' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'x_soil',               &
                                      attribute_name = char_unit, value = 'm' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'y_soil',               &
                                      attribute_name = char_unit, value = 'm' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'E_UTM_soil',           &
                                      attribute_name = char_unit, value = 'm' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'N_UTM_soil',           &
                                      attribute_name = char_unit, value = 'm' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lat_soil',             &
                                      attribute_name = char_unit, value = 'degrees_north' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lon_soil',             &
                                      attribute_name = char_unit, value = 'degrees_east' )
!
!--       Long name
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_name_soil',    &
                                      attribute_name = char_long, value = 'station name')
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time_soil',            &
                                      attribute_name = char_long, value = 'time')
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z_soil',               &
                                      attribute_name = char_long, value = 'height above origin' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_h_soil',       &
                                      attribute_name = char_long, value = 'surface altitude' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'x_soil',               &
                                      attribute_name = char_long,                                  &
                                      value = 'distance to origin in x-direction' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'y_soil',               &
                                      attribute_name = char_long,                                  &
                                      value = 'distance to origin in y-direction' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'E_UTM_soil',           &
                                      attribute_name = char_long, value = 'easting' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'N_UTM_soil',           &
                                      attribute_name = char_long, value = 'northing' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lat_soil',             &
                                      attribute_name = char_long, value = 'latitude' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lon_soil',             &
                                      attribute_name = char_long, value = 'longitude' )
!
!--       Standard name
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_name_soil',    &
                                      attribute_name = char_standard, value = 'platform_name')
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time_soil',            &
                                      attribute_name = char_standard, value = 'time')
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z_soil',               &
                                      attribute_name = char_standard,                              &
                                      value = 'height_above_mean_sea_level' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_h_soil',       &
                                      attribute_name = char_standard, value = 'surface_altitude' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'E_UTM_soil',           &
                                      attribute_name = char_standard,                              &
                                      value = 'projection_x_coordinate' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'N_UTM_soil',           &
                                      attribute_name = char_standard,                              &
                                      value = 'projection_y_coordinate' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lat_soil',             &
                                      attribute_name = char_standard, value = 'latitude' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'lon_soil',             &
                                      attribute_name = char_standard, value = 'longitude' )
!
!--       Axis
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time_soil',            &
                                      attribute_name = 'axis', value = 'T')
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z_soil',               &
                                      attribute_name = 'axis', value = 'Z' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'x_soil',               &
                                      attribute_name = 'axis', value = 'X' )
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'y_soil',               &
                                      attribute_name = 'axis', value = 'Y' )
!
!--       Set further individual attributes for the coordinate variables.
!--       For station name soil
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'station_name_soil',    &
                                      attribute_name = 'cf_role', value = 'timeseries_id' )
!
!--       For time soil
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time_soil',            &
                                      attribute_name = 'calendar', value = 'proleptic_gregorian' )
!           return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'time_soil',            &
!                                       attribute_name = 'bounds', value = 'time_bounds' )
!
!--       For z soil
          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = 'z_soil',               &
                                      attribute_name = 'positive', value = 'up' )
       ENDIF
!
!--    Define variables that shall be sampled.
       DO  n = 1, vmea(l)%nmeas
          variable_name = TRIM( vmea(l)%var_atts(n)%name )
!
!--       In order to link the correct dimension names, atmosphere and soil variables need to be
!--       distinguished.
          IF ( vmea(l)%soil_sampling  .AND.                                                        &
               ANY( TRIM( vmea(l)%var_atts(n)%name) == soil_vars ) )  THEN

             return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,       &
                                         dimension_names = (/'station_soil', 'ntime_soil  '/),     &
                                         write_mode = 'independent', output_type = 'real32' )
          ELSE

             return_value = dom_def_var( vmea(l)%nc_filename, variable_name = variable_name,       &
                                         dimension_names = (/'station', 'ntime  '/),               &
                                         write_mode = 'independent', output_type = 'real32' )
          ENDIF
!
!--       Set variable attributes. For some variables not all attributes are defined,
!--       e.g. standard_name for the horizontal wind components.
          CALL vm_set_attributes( vmea(l)%var_atts(n) )

          IF ( vmea(l)%var_atts(n)%long_name /= 'none' )  THEN
             return_value = dom_def_att( vmea(l)%nc_filename,  variable_name = variable_name,      &
                                         attribute_name = char_long,                               &
                                         value = TRIM( vmea(l)%var_atts(n)%long_name ) )
          ENDIF
          IF ( vmea(l)%var_atts(n)%standard_name /= 'none' )  THEN
             return_value = dom_def_att( vmea(l)%nc_filename, variable_name = variable_name,       &
                                         attribute_name = char_standard,                           &
                                         value = TRIM( vmea(l)%var_atts(n)%standard_name ) )
          ENDIF
          IF ( vmea(l)%var_atts(n)%units /= 'none' )  THEN
             return_value = dom_def_att( vmea(l)%nc_filename, variable_name = variable_name,       &
                                         attribute_name = char_unit,                               &
                                         value = TRIM( vmea(l)%var_atts(n)%units ) )
          ENDIF

          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      attribute_name = 'grid_mapping',                             &
                                      value = TRIM( vmea(l)%var_atts(n)%grid_mapping ) )

          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      attribute_name = 'coordinates',                              &
                                      value = TRIM( vmea(l)%var_atts(n)%coordinates ) )

          return_value = dom_def_att( vmea(l)%nc_filename, variable_name = variable_name,          &
                                      attribute_name = char_fill,                                  &
                                      value = REAL( vmea(l)%var_atts(n)%fill_value, KIND=4 ) )

       ENDDO  ! loop over variables per site

    ENDDO  ! loop over sites

!
!-- Check if DOM reported any error
    dom_error_message = dom_get_error_message()
    IF ( TRIM( dom_error_message ) /= '' )  THEN
       message_string = 'error while defining output: "' // TRIM( dom_error_message ) // '"'
       CALL message( 'vm_init_output', 'VME0008', 0, 1, 0, 6, 0 )
    ENDIF

 END SUBROUTINE vm_init_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parallel NetCDF output via data-output module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_data_output( output_pr, output_ts, output_tr )

    CHARACTER(LEN=100) ::  variable_name  !< name of output variable

    CHARACTER(LEN=maximum_name_length), DIMENSION(:), ALLOCATABLE :: station_name  !< string for station name, consecutively ordered

    CHARACTER(LEN=1), DIMENSION(:,:), ALLOCATABLE, TARGET ::  output_values_2d_char_target  !< target for output name arrays
    CHARACTER(LEN=1), DIMENSION(:,:), POINTER             ::  output_values_2d_char_pointer !< pointer for output name arrays

    INTEGER(iwp)       ::  l             !< loop index for the number of sites
    INTEGER(iwp)       ::  n             !< loop index for observation points
    INTEGER(iwp)       ::  nn            !< loop index for number of characters in a name
    INTEGER            ::  return_value  !< returned status value of called function
    INTEGER(iwp)       ::  t_ind         !< time index

    LOGICAL ::  output_pr  !< flag indicating profile data to be output
    LOGICAL ::  output_ts  !< flag indicating timeseries data to be output
    LOGICAL ::  output_tr  !< flag indicating trajectory data to be output

    REAL(wp) ::  eutm  !< easting (UTM)
    REAL(wp) ::  nutm  !< northing (UTM)

    REAL(wp), DIMENSION(:), ALLOCATABLE           ::  dum_lat                   !< transformed geographical coordinate (latitude)
    REAL(wp), DIMENSION(:), ALLOCATABLE           ::  dum_lon                   !< transformed geographical coordinate (longitude)
    REAL(wp), DIMENSION(:), ALLOCATABLE           ::  oro_rel                   !< relative altitude of model surface
    REAL(sp), DIMENSION(:), POINTER               ::  output_values_1d_pointer  !< pointer for 1d output array
    REAL(sp), DIMENSION(:), ALLOCATABLE, TARGET   ::  output_values_1d_target   !< target for 1d output array
    REAL(sp), DIMENSION(:,:), POINTER             ::  output_values_2d_pointer  !< pointer for 2d output array
    REAL(sp), DIMENSION(:,:), ALLOCATABLE, TARGET ::  output_values_2d_target   !< target for 2d output array


    CALL cpu_log( log_point_s(26), 'VM output', 'start' )

!
!-- At the first call of this routine write the spatial coordinates.
    IF ( .NOT. initial_write_coordinates )  THEN
!
!--    Write spatial coordinates.
       DO  l = 1, vmea_general%nvm
!
!--       Skip if no observations were taken.
          IF ( vmea(l)%ns_tot == 0  .AND.  vmea(l)%ns_soil_tot == 0 )  CYCLE

          ALLOCATE( output_values_1d_target(vmea(l)%start_coord_a:vmea(l)%end_coord_a) )
!
!--       Output of Easting coordinate. Before output, recalculate EUTM. Note, in case of single-
!--       point timeseries output and exact sampling location output, the sampling coordinates
!--       are already available in the %coords array.
          IF ( vmea(l)%interpolate  .AND.  SIZE( output_values_1d_target ) == 1 )  THEN
             output_values_1d_target = vmea_general%origin_x_pd                                    &
                                    + vmea(l)%coords(1) *                                          &
                                      COS( init_model%rotation_angle * pi / 180.0_wp )             &
                                    + vmea(l)%coords(2) *                                          &
                                      SIN( init_model%rotation_angle * pi / 180.0_wp )
          ELSE
             output_values_1d_target = vmea_general%origin_x_pd                                    &
                                    + REAL( vmea(l)%i(1:vmea(l)%ns) + 0.5_wp, KIND = wp ) * dx *   &
                                      COS( init_model%rotation_angle * pi / 180.0_wp )             &
                                    + REAL( vmea(l)%j(1:vmea(l)%ns) + 0.5_wp, KIND = wp ) * dy *   &
                                      SIN( init_model%rotation_angle * pi / 180.0_wp )
          ENDIF

          output_values_1d_pointer => output_values_1d_target
          return_value = dom_write_var( vmea(l)%nc_filename, 'E_UTM',                              &
                                        values_real32_1d = output_values_1d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_a/),                  &
                                        bounds_end   = (/vmea(l)%end_coord_a  /) )
!
!--       Output of Northing coordinate. Before output, recalculate NUTM. Note, in case of single-
!--       point timeseries output and exact sampling location output, the sampling coordinates
!--       are already available in the %coords array.
          IF ( vmea(l)%interpolate  .AND.  SIZE( output_values_1d_target ) == 1 )  THEN
             output_values_1d_target = vmea_general%origin_y_pd                                    &
                                    - vmea(l)%coords(1) *                                          &
                                      SIN( init_model%rotation_angle * pi / 180.0_wp )             &
                                    + vmea(l)%coords(2) *                                          &
                                      COS( init_model%rotation_angle * pi / 180.0_wp )
          ELSE
             output_values_1d_target = vmea_general%origin_y_pd                                    &
                                    - REAL( vmea(l)%i(1:vmea(l)%ns) + 0.5_wp, KIND = wp ) * dx *   &
                                      SIN( init_model%rotation_angle * pi / 180.0_wp )             &
                                    + REAL( vmea(l)%j(1:vmea(l)%ns) + 0.5_wp, KIND = wp ) * dy *   &
                                      COS( init_model%rotation_angle * pi / 180.0_wp )
          ENDIF

          output_values_1d_pointer => output_values_1d_target
          return_value = dom_write_var( vmea(l)%nc_filename, 'N_UTM',                              &
                                        values_real32_1d = output_values_1d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_a/),                  &
                                        bounds_end   = (/vmea(l)%end_coord_a  /) )
!
!--       Output of longitude and latitude coordinate. Before output, convert it.
          ALLOCATE( dum_lat(1:vmea(l)%ns) )
          ALLOCATE( dum_lon(1:vmea(l)%ns) )
!
!--       In case of single-point timeseries output and exact sampling location output, the
!--       sampling coordinates are already available in the %coords array.
          IF ( vmea(l)%interpolate  .AND.  SIZE( output_values_1d_target ) == 1 )  THEN
             DO  n = 1, vmea(l)%ns
                eutm = vmea_general%origin_x_pd + vmea(l)%coords(1) *                              &
                                                  COS( init_model%rotation_angle * pi / 180.0_wp ) &
                                                + vmea(l)%coords(2) *                              &
                                                  SIN( init_model%rotation_angle * pi / 180.0_wp )
                nutm = vmea_general%origin_y_pd - vmea(l)%coords(1) *                              &
                                                  SIN( init_model%rotation_angle * pi / 180.0_wp ) &
                                                + vmea(l)%coords(2) *                              &
                                                  COS( init_model%rotation_angle * pi / 180.0_wp )
                CALL convert_utm_to_geographic( crs_list, eutm, nutm, dum_lon(n), dum_lat(n) )
             ENDDO
          ELSE
             DO  n = 1, vmea(l)%ns
                eutm = vmea_general%origin_x_pd + REAL( vmea(l)%i(n) + 0.5_wp, KIND = wp ) * dx *  &
                                                  COS( init_model%rotation_angle * pi / 180.0_wp ) &
                                                + REAL( vmea(l)%j(n) + 0.5_wp, KIND = wp ) * dy *  &
                                                  SIN( init_model%rotation_angle * pi / 180.0_wp )
                nutm = vmea_general%origin_y_pd - REAL( vmea(l)%i(n) + 0.5_wp, KIND = wp ) * dx *  &
                                                  SIN( init_model%rotation_angle * pi / 180.0_wp ) &
                                                + REAL( vmea(l)%j(n) + 0.5_wp, KIND = wp ) * dy *  &
                                                  COS( init_model%rotation_angle * pi / 180.0_wp )
                CALL convert_utm_to_geographic( crs_list, eutm, nutm, dum_lon(n), dum_lat(n) )
             ENDDO
          ENDIF

          output_values_1d_target = dum_lat
          output_values_1d_pointer => output_values_1d_target
          return_value = dom_write_var( vmea(l)%nc_filename, 'lat',                                &
                                        values_real32_1d = output_values_1d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_a/),                  &
                                        bounds_end   = (/vmea(l)%end_coord_a  /) )

          output_values_1d_target = dum_lon
          output_values_1d_pointer => output_values_1d_target
          return_value = dom_write_var( vmea(l)%nc_filename, 'lon',                                &
                                        values_real32_1d = output_values_1d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_a/),                  &
                                        bounds_end   = (/vmea(l)%end_coord_a  /) )
          DEALLOCATE( dum_lat )
          DEALLOCATE( dum_lon )
!
!--       Output of relative height coordinate.
!--       Before this is output, first define the relative orography height and add this to z.
          ALLOCATE( oro_rel(1:vmea(l)%ns) )
          DO  n = 1, vmea(l)%ns
             oro_rel(n) = zw(topo_top_ind(vmea(l)%j(n),vmea(l)%i(n),3))
          ENDDO
!
!--       In case of single-point timeseries output and exact sampling location output, the
!--       sampling coordinates are already available in the %coords array.
          IF ( vmea(l)%interpolate  .AND.  SIZE( output_values_1d_target ) == 1 )  THEN
             output_values_1d_target = vmea(l)%coords(3)
          ELSE
             output_values_1d_target = vmea(l)%zar(1:vmea(l)%ns) + oro_rel(:)
          ENDIF

          output_values_1d_pointer => output_values_1d_target
          return_value = dom_write_var( vmea(l)%nc_filename, 'z',                                  &
                                        values_real32_1d = output_values_1d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_a/),                  &
                                        bounds_end   = (/vmea(l)%end_coord_a  /) )
!
!--       Write surface altitude for the station. Note, since z is already a relative observation
!--       height, station_h must be zero, in order to obtain the observation level.
          output_values_1d_target = oro_rel(:)
          output_values_1d_pointer => output_values_1d_target
          return_value = dom_write_var( vmea(l)%nc_filename, 'station_h',                          &
                                        values_real32_1d = output_values_1d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_a/),                  &
                                        bounds_end   = (/vmea(l)%end_coord_a  /) )

          DEALLOCATE( oro_rel )
          DEALLOCATE( output_values_1d_target )
!
!--       Write station name
          ALLOCATE ( station_name(vmea(l)%start_coord_a:vmea(l)%end_coord_a) )
          ALLOCATE ( output_values_2d_char_target(vmea(l)%start_coord_a:vmea(l)%end_coord_a,       &
                                                  1:maximum_name_length) )

          DO  n = vmea(l)%start_coord_a, vmea(l)%end_coord_a
             station_name(n) = REPEAT( ' ', maximum_name_length )
             WRITE( station_name(n), '(A,I10.10)') "station", n
             DO  nn = 1, maximum_name_length
                output_values_2d_char_target(n,nn) = station_name(n)(nn:nn)
             ENDDO
          ENDDO

          output_values_2d_char_pointer => output_values_2d_char_target

          return_value = dom_write_var( vmea(l)%nc_filename, 'station_name',                       &
                                        values_char_2d = output_values_2d_char_pointer,            &
                                        bounds_start = (/ 1, vmea(l)%start_coord_a /),             &
                                        bounds_end   = (/ maximum_name_length,                     &
                                        vmea(l)%end_coord_a /) )

          DEALLOCATE( station_name )
          DEALLOCATE( output_values_2d_char_target )
!
!--       In case of sampled soil quantities, output also the respective coordinate arrays.
          IF ( vmea(l)%soil_sampling )  THEN
             ALLOCATE( output_values_1d_target(vmea(l)%start_coord_s:vmea(l)%end_coord_s) )
!
!--          Output of Easting coordinate. Before output, recalculate EUTM.
             output_values_1d_target = vmea_general%origin_x_pd                                    &
                            + REAL( vmea(l)%i_soil(1:vmea(l)%ns_soil) + 0.5_wp, KIND = wp ) * dx * &
                              COS( init_model%rotation_angle * pi / 180.0_wp )                     &
                            + REAL( vmea(l)%j_soil(1:vmea(l)%ns_soil) + 0.5_wp, KIND = wp ) * dy * &
                              SIN( init_model%rotation_angle * pi / 180.0_wp )
             output_values_1d_pointer => output_values_1d_target
             return_value = dom_write_var( vmea(l)%nc_filename, 'E_UTM_soil',                      &
                                           values_real32_1d = output_values_1d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s/),               &
                                           bounds_end   = (/vmea(l)%end_coord_s  /) )
!
!--          Output of Northing coordinate. Before output, recalculate NUTM.
             output_values_1d_target = vmea_general%origin_y_pd                                    &
                            - REAL( vmea(l)%i_soil(1:vmea(l)%ns_soil) + 0.5_wp, KIND = wp ) * dx * &
                              SIN( init_model%rotation_angle * pi / 180.0_wp )                     &
                            + REAL( vmea(l)%j_soil(1:vmea(l)%ns_soil) + 0.5_wp, KIND = wp ) * dy * &
                              COS( init_model%rotation_angle * pi / 180.0_wp )

             output_values_1d_pointer => output_values_1d_target
             return_value = dom_write_var( vmea(l)%nc_filename, 'N_UTM_soil',                      &
                                           values_real32_1d = output_values_1d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s/),               &
                                           bounds_end   = (/vmea(l)%end_coord_s  /) )
!
!--          Output of longitude and latitude coordinate. Before output, convert it.
             ALLOCATE( dum_lat(1:vmea(l)%ns_soil) )
             ALLOCATE( dum_lon(1:vmea(l)%ns_soil) )

             DO  n = 1, vmea(l)%ns_soil
                eutm = vmea_general%origin_x_pd                           &
                                            + REAL( vmea(l)%i_soil(n) + 0.5_wp, KIND = wp ) * dx * &
                                              COS( init_model%rotation_angle * pi / 180.0_wp )     &
                                            + REAL( vmea(l)%j_soil(n) + 0.5_wp, KIND = wp ) * dy * &
                                              SIN( init_model%rotation_angle * pi / 180.0_wp )
                nutm = vmea_general%origin_y_pd                           &
                                            - REAL( vmea(l)%i_soil(n) + 0.5_wp, KIND = wp ) * dx * &
                                              SIN( init_model%rotation_angle * pi / 180.0_wp )     &
                                            + REAL( vmea(l)%j_soil(n) + 0.5_wp, KIND = wp ) * dy * &
                                              COS( init_model%rotation_angle * pi / 180.0_wp )
                CALL convert_utm_to_geographic( crs_list, eutm, nutm, dum_lon(n), dum_lat(n) )
             ENDDO

             output_values_1d_target = dum_lat
             output_values_1d_pointer => output_values_1d_target
             return_value = dom_write_var( vmea(l)%nc_filename, 'lat_soil',                        &
                                           values_real32_1d = output_values_1d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s/),               &
                                           bounds_end   = (/vmea(l)%end_coord_s  /) )

             output_values_1d_target = dum_lon
             output_values_1d_pointer => output_values_1d_target
             return_value = dom_write_var( vmea(l)%nc_filename, 'lon_soil',                        &
                                           values_real32_1d = output_values_1d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s/),               &
                                           bounds_end   = (/vmea(l)%end_coord_s  /) )
             DEALLOCATE( dum_lat )
             DEALLOCATE( dum_lon )
!
!--          Output of relative height coordinate.
!--          Before this is output, first define the relative orographie height and add this to z.
             ALLOCATE( oro_rel(1:vmea(l)%ns_soil) )
             DO  n = 1, vmea(l)%ns_soil
                oro_rel(n) = zw(topo_top_ind(vmea(l)%j_soil(n),vmea(l)%i_soil(n),3))
             ENDDO

             output_values_1d_target = vmea(l)%depth(1:vmea(l)%ns_soil) + oro_rel(:)
             output_values_1d_pointer => output_values_1d_target
             return_value = dom_write_var( vmea(l)%nc_filename, 'z_soil',                          &
                                           values_real32_1d = output_values_1d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s/),               &
                                           bounds_end   = (/vmea(l)%end_coord_s  /) )
!
!--          Write surface altitude for the station. Note, since z is already a relative observation
!--          height, station_h must be zero, in order to obtain the observation level.
             output_values_1d_target = oro_rel(:)
             output_values_1d_pointer => output_values_1d_target
             return_value = dom_write_var( vmea(l)%nc_filename, 'station_h_soil',                  &
                                           values_real32_1d = output_values_1d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s/),               &
                                           bounds_end   = (/vmea(l)%end_coord_s  /) )

             DEALLOCATE( oro_rel )
             DEALLOCATE( output_values_1d_target )
!
!--          Write station name
             ALLOCATE ( station_name(vmea(l)%start_coord_s:vmea(l)%end_coord_s) )
             ALLOCATE ( output_values_2d_char_target(vmea(l)%start_coord_s:vmea(l)%end_coord_s,    &
                                                     1:maximum_name_length) )

             DO  n = vmea(l)%start_coord_s, vmea(l)%end_coord_s
                station_name(n) = REPEAT( ' ', maximum_name_length )
                WRITE( station_name(n), '(A,I10.10)') "station", n
                DO  nn = 1, maximum_name_length
                   output_values_2d_char_target(n,nn) = station_name(n)(nn:nn)
                ENDDO
             ENDDO
             output_values_2d_char_pointer => output_values_2d_char_target

             return_value = dom_write_var( vmea(l)%nc_filename, 'station_name_soil',               &
                                           values_char_2d = output_values_2d_char_pointer,         &
                                           bounds_start = (/ 1, vmea(l)%start_coord_s /),          &
                                           bounds_end   = (/ maximum_name_length,                  &
                                           vmea(l)%end_coord_s   /) )

             DEALLOCATE( station_name )
             DEALLOCATE( output_values_2d_char_target )
!
!--          At the beginning, check if also time coordinate for soil variables need to be written.
!--          This is required when the input flag soil_sampling is .T. and also a soil-related
!--          variable has been sampled. The latter is not necessarily true if the driver input
!--          does not indicate any soil variable. Store the result in a dedicated variable.
             DO  n = 1, vmea(l)%nmeas
                IF ( ANY( TRIM( vmea(l)%var_atts(n)%name ) == soil_vars )  .AND.                   &
                     vmea(l)%var_atts(n)%sampled )  vmea(l)%soil_var_out = .TRUE.
             ENDDO
          ENDIF

       ENDDO  ! loop over sites

       initial_write_coordinates = .TRUE.

    ENDIF
!
!-- Loop over all sites.
    DO  l = 1, vmea_general%nvm
!
!--    Skip the measurement if its feature type does not match with the type requested to output.
       IF ( ( vmea(l)%timeseries          .AND.  output_ts )  .OR.                                 &
            ( vmea(l)%timeseries_profile  .AND.  output_pr )  .OR.                                 &
            ( vmea(l)%trajectory          .AND.  output_tr ) )                                     &
       THEN
          CONTINUE
       ELSE
          CYCLE
       ENDIF
!
!--    Skip if no observations were taken.
       IF ( vmea(l)%ns_tot == 0  .AND.  vmea(l)%ns_soil_tot == 0 )  CYCLE
!
!--    Determine time index in file.
       t_ind = vmea(l)%file_time_index + 1
!
!--    Write time coordinate to file.
       variable_name = 'time'

       ALLOCATE( output_values_2d_target(t_ind:t_ind,vmea(l)%start_coord_a:vmea(l)%end_coord_a) )
       output_values_2d_target(t_ind,:) = time_since_reference_point
       output_values_2d_pointer => output_values_2d_target

       return_value = dom_write_var( vmea(l)%nc_filename, variable_name,                           &
                                     values_real32_2d = output_values_2d_pointer,                  &
                                     bounds_start = (/vmea(l)%start_coord_a, t_ind/),              &
                                     bounds_end   = (/vmea(l)%end_coord_a, t_ind/) )

       DEALLOCATE( output_values_2d_target )
!
!--    Check if also time coordinate for soil variables need to be written.
       IF( vmea(l)%soil_var_out )  THEN
          variable_name = 'time_soil'

          ALLOCATE( output_values_2d_target(t_ind:t_ind,vmea(l)%start_coord_s:vmea(l)%end_coord_s) )
          output_values_2d_target(t_ind,:) = time_since_reference_point
          output_values_2d_pointer => output_values_2d_target

          return_value = dom_write_var( vmea(l)%nc_filename, variable_name,                        &
                                        values_real32_2d = output_values_2d_pointer,               &
                                        bounds_start = (/vmea(l)%start_coord_s, t_ind/),           &
                                        bounds_end   = (/vmea(l)%end_coord_s, t_ind /) )

          DEALLOCATE( output_values_2d_target )
       ENDIF
!
!--    Write output variables. Distinguish between atmosphere and soil variables.
       DO  n = 1, vmea(l)%nmeas
!
!--       Skip writing if the variable was not sampled.
          IF ( .NOT. vmea(l)%var_atts(n)%sampled )  CYCLE
!
!--       Output soil quantities.
          IF ( vmea(l)%soil_sampling  .AND.                                                        &
               ANY( TRIM( vmea(l)%var_atts(n)%name ) == soil_vars ) )  THEN

             variable_name = TRIM( vmea(l)%var_atts(n)%name )

             ALLOCATE( output_values_2d_target(t_ind:t_ind,vmea(l)%start_coord_s:vmea(l)%end_coord_s) )
             output_values_2d_target(t_ind,:) = vmea(l)%measured_vars_soil(:,n)
             output_values_2d_pointer => output_values_2d_target
             return_value = dom_write_var( vmea(l)%nc_filename, variable_name,                     &
                                           values_real32_2d = output_values_2d_pointer,            &
                                           bounds_start = (/vmea(l)%start_coord_s, t_ind/),        &
                                           bounds_end   = (/vmea(l)%end_coord_s, t_ind  /) )

             DEALLOCATE( output_values_2d_target )
!
!--       Output atmosphere quantities.
          ELSE

             variable_name = TRIM( vmea(l)%var_atts(n)%name )
             ALLOCATE( output_values_2d_target(t_ind:t_ind,vmea(l)%start_coord_a:vmea(l)%end_coord_a) )
             output_values_2d_target(t_ind,:) = vmea(l)%measured_vars(:,n)
             output_values_2d_pointer => output_values_2d_target
             return_value = dom_write_var( vmea(l)%nc_filename, variable_name,                     &
                                           values_real32_2d = output_values_2d_pointer,            &
                                           bounds_start = (/ vmea(l)%start_coord_a, t_ind /),      &
                                           bounds_end   = (/ vmea(l)%end_coord_a, t_ind /) )

             DEALLOCATE( output_values_2d_target )

          ENDIF
       ENDDO
!
!--    Update number of written time indices
       vmea(l)%file_time_index = t_ind

    ENDDO  ! loop over sites

!
!-- Check if DOM reported any error
    dom_error_message = dom_get_error_message()
    IF ( TRIM( dom_error_message ) /= '' )  THEN
       message_string = 'error while writing output: "' // TRIM( dom_error_message ) // '"'
       CALL message( 'vm_data_output', 'VME0009', 0, 1, 0, 6, 0 )
    ENDIF

    CALL cpu_log( log_point_s(26), 'VM output', 'stop' )

 END SUBROUTINE vm_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sampling of the actual quantities along the observation coordinates
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE vm_sampling( sample_pr, sample_ts, sample_tr )


     INTEGER(iwp) ::  i         !< grid index in x-direction
     INTEGER(iwp) ::  j         !< grid index in y-direction
     INTEGER(iwp) ::  k         !< grid index in z-direction
     INTEGER(iwp) ::  ind_chem  !< dummy index to identify chemistry variable and translate it from (UC)2 standard to interal naming
     INTEGER(iwp) ::  l         !< running index over the number of stations
     INTEGER(iwp) ::  m         !< running index over all virtual observation coordinates
     INTEGER(iwp) ::  mm        !< index of surface element which corresponds to the virtual observation coordinate
     INTEGER(iwp) ::  n         !< running index over all measured variables at a station
     INTEGER(iwp) ::  nn        !< running index over the number of chemcal species

     LOGICAL ::  sample_pr  !< flag indicating profile data to be sampled
     LOGICAL ::  sample_ts  !< flag indicating timeseries data to be sampled
     LOGICAL ::  sample_tr  !< flag indicating trajectory data to be sampled

     REAL(wp) ::  e_s   !< saturation water vapor pressure
     REAL(wp) ::  q_s   !< saturation mixing ratio
     REAL(wp) ::  q_wv  !< mixing ratio


     CALL cpu_log( log_point_s(27), 'VM sampling', 'start' )
!
!--  Loop over all sites.
     DO  l = 1, vmea_general%nvm
!
!--     Skip the measurement if its feature type does not match with the type requested to sample.
        IF ( ( vmea(l)%timeseries          .AND.  sample_ts )  .OR.                                &
             ( vmea(l)%timeseries_profile  .AND.  sample_pr )  .OR.                                &
             ( vmea(l)%trajectory          .AND.  sample_tr ) )                                    &
        THEN
           CONTINUE
        ELSE
           CYCLE
        ENDIF
!
!--     At the beginning, set _FillValues
        IF ( ALLOCATED( vmea(l)%measured_vars      ) ) vmea(l)%measured_vars      = vmea(l)%fillout
        IF ( ALLOCATED( vmea(l)%measured_vars_soil ) ) vmea(l)%measured_vars_soil = vmea(l)%fillout
!
!--     Loop over all variables measured at this site.
        DO  n = 1, vmea(l)%nmeas
!
!--        Set sampling flag. Only variables that have been actually sampled are output.
!--        Variables that are not treated by the model are skipped.
           vmea(l)%var_atts(n)%sampled = .FALSE.

           SELECT CASE ( TRIM( vmea(l)%var_atts(n)%name ) )

              CASE ( 'theta' ) ! potential temperature
                 IF ( .NOT. neutral )  THEN
                    IF ( vmea(l)%interpolate )  THEN
                       DO  m = 1, vmea(l)%ns
                          k = vmea(l)%k(m)
                          j = vmea(l)%j(m)
                          i = vmea(l)%i(m)
                          vmea(l)%measured_vars(m,n) =                                             &
                                    interpolate_to_sampling_location( pt, k, j, i, 0.5_wp, 0.5_wp, &
                                                                      zu, vmea(l)%coords(:) )
                       ENDDO
                    ELSE
                       DO  m = 1, vmea(l)%ns
                          k = vmea(l)%k(m)
                          j = vmea(l)%j(m)
                          i = vmea(l)%i(m)
                          vmea(l)%measured_vars(m,n) = pt(k,j,i)
                       ENDDO
                    ENDIF
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'ta' ) ! absolute temperature
                 IF ( .NOT. neutral )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = pt(k,j,i) * exner( k ) - degc_to_k
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'hus' ) ! mixing ratio
                 IF ( humidity )  THEN
                    IF ( vmea(l)%interpolate )  THEN
                       DO  m = 1, vmea(l)%ns
                          k = vmea(l)%k(m)
                          j = vmea(l)%j(m)
                          i = vmea(l)%i(m)
                          vmea(l)%measured_vars(m,n) =                                             &
                                    interpolate_to_sampling_location( q, k, j, i, 0.5_wp, 0.5_wp,  &
                                                                      zu, vmea(l)%coords(:) )
                       ENDDO
                    ELSE
                       DO  m = 1, vmea(l)%ns
                          k = vmea(l)%k(m)
                          j = vmea(l)%j(m)
                          i = vmea(l)%i(m)
                          vmea(l)%measured_vars(m,n) = q(k,j,i)
                       ENDDO
                    ENDIF
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'haa' ) ! absolute humidity
                 IF ( humidity )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = ( q(k,j,i) / ( 1.0_wp - q(k,j,i) ) ) *         &
                                                    rho_air(k)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'pwv' ) ! water vapor partial pressure
                 IF ( humidity )  THEN
!                     DO  m = 1, vmea(l)%ns
!                        k = vmea(l)%k(m)
!                        j = vmea(l)%j(m)
!                        i = vmea(l)%i(m)
!                        vmea(l)%measured_vars(m,n) = ( q(k,j,i) / ( 1.0_wp - q(k,j,i) ) )          &
!                                                     * rho_air(k)
!                     ENDDO
                 ENDIF

              CASE ( 'hur' ) ! relative humidity
                 IF ( humidity )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
!
!--                    Calculate actual temperature, water vapor saturation pressure and, based on
!--                    this, the saturation mixing ratio.
                       e_s  = magnus( exner(k) * pt(k,j,i) )
                       q_s  = rd_d_rv * e_s / ( hyp(k) - e_s )
                       q_wv = ( q(k,j,i) / ( 1.0_wp - q(k,j,i) ) ) * rho_air(k)

                       vmea(l)%measured_vars(m,n) = q_wv / ( q_s + 1E-10_wp )
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'u', 'ua' ) ! u-component
                 IF ( vmea(l)%interpolate )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) =                                                &
                                    interpolate_to_sampling_location( u, k, j, i, 0.5_wp, 0.0_wp,  &
                                                                      zu, vmea(l)%coords(:) )
                    ENDDO
                 ELSE
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
                    ENDDO
                 ENDIF
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'v', 'va' ) ! v-component
                 IF ( vmea(l)%interpolate )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) =                                                &
                                    interpolate_to_sampling_location( v, k, j, i, 0.0_wp, 0.5_wp,  &
                                                                      zu, vmea(l)%coords(:) )
                    ENDDO
                 ELSE
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )
                    ENDDO
                 ENDIF
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'w' ) ! w-component
                 IF ( vmea(l)%interpolate )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) =                                                &
                                    interpolate_to_sampling_location( w, k, j, i, 0.5_wp, 0.5_wp,  &
                                                                      zw, vmea(l)%coords(:) )
                    ENDDO
                 ELSE
                    DO  m = 1, vmea(l)%ns
                       k = MAX ( 1, vmea(l)%k(m) )
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                    ENDDO
                 ENDIF
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'wspeed' ) ! horizontal wind speed
                 DO  m = 1, vmea(l)%ns
                    k = vmea(l)%k(m)
                    j = vmea(l)%j(m)
                    i = vmea(l)%i(m)
                    vmea(l)%measured_vars(m,n) = SQRT(   ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) )**2 &
                                                       + ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) )**2 &
                                                     )
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'wdir' ) ! wind direction
                 DO  m = 1, vmea(l)%ns
                    k = vmea(l)%k(m)
                    j = vmea(l)%j(m)
                    i = vmea(l)%i(m)
                    vmea(l)%measured_vars(m,n) = 180.0_wp + 180.0_wp / pi * ATAN2(                 &
                                                               0.5_wp * ( v(k,j,i) + v(k,j+1,i) ), &
                                                               0.5_wp * ( u(k,j,i) + u(k,j,i+1) )  &
                                                                                 )
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'utheta' )
                 IF ( .NOT. neutral )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) * pt(k,j,i)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'vtheta' )
                 IF ( .NOT. neutral )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) * pt(k,j,i)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'wtheta' )
                 IF ( .NOT. neutral )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = MAX ( 1, vmea(l)%k(m) )
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( w(k-1,j,i) + w(k,j,i) ) * pt(k,j,i)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'uqv' )
                 IF ( humidity )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) * q(k,j,i)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'vqv' )
                 IF ( humidity )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = vmea(l)%k(m)
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) * q(k,j,i)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'wqv' )
                 IF ( humidity )  THEN
                    DO  m = 1, vmea(l)%ns
                       k = MAX ( 1, vmea(l)%k(m) )
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)
                       vmea(l)%measured_vars(m,n) = 0.5_wp * ( w(k-1,j,i) + w(k,j,i) ) * q(k,j,i)
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'uw' )
                 DO  m = 1, vmea(l)%ns
                    k = MAX ( 1, vmea(l)%k(m) )
                    j = vmea(l)%j(m)
                    i = vmea(l)%i(m)
                    vmea(l)%measured_vars(m,n) = 0.25_wp * ( w(k-1,j,i) + w(k,j,i) ) *             &
                                                           ( u(k,j,i)   + u(k,j,i+1) )
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'vw' )
                 DO  m = 1, vmea(l)%ns
                    k = MAX ( 1, vmea(l)%k(m) )
                    j = vmea(l)%j(m)
                    i = vmea(l)%i(m)
                    vmea(l)%measured_vars(m,n) = 0.25_wp * ( w(k-1,j,i) + w(k,j,i) ) *             &
                                                           ( v(k,j,i)   + v(k,j+1,i) )
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'uv' )
                 DO  m = 1, vmea(l)%ns
                    k = vmea(l)%k(m)
                    j = vmea(l)%j(m)
                    i = vmea(l)%i(m)
                    vmea(l)%measured_vars(m,n) = 0.25_wp * ( u(k,j,i)   + u(k,j,i+1) ) *           &
                                                           ( v(k,j,i)   + v(k,j+1,i) )
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'us' ) ! friction velocity
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points. Hence,
!--                 limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_def%start_index(j,i), surf_def%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_def%us(mm),                        &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_def%upward(mm) )
                    ENDDO
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%us(mm),                        &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%us(mm),                        &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'thetas' ) ! scaling parameter temperature
                 IF ( .NOT. neutral )  THEN
                    DO  m = 1, vmea(l)%ns
!
!--                    Surface data is only available on inner subdomains, not on ghost points.
!-                     Hence, limit the indices.
                       j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                       j = MERGE( j           , nyn, j            < nyn )
                       i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                       i = MERGE( i           , nxr, i            < nxr )
!
!--                    Take only values from horizontally-upward facing surfaces.
                       DO  mm = surf_def%start_index(j,i), surf_def%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_def%ts(mm),                     &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_def%upward(mm) )
                       ENDDO
                       DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%ts(mm),                     &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_lsm%upward(mm) )
                       ENDDO
                       DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_usm%ts(mm),                     &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_usm%upward(mm) )
                       ENDDO
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'hfls' ) ! surface latent heat flux
                 IF ( humidity )  THEN
                    DO  m = 1, vmea(l)%ns
!
!--                    Surface data is only available on inner subdomains, not on ghost points.
!--                    Hence, limit the indices.
                       j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                       j = MERGE( j           , nyn, j            < nyn )
                       i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                       i = MERGE( i           , nxr, i            < nxr )
                       k = vmea(l)%k(m)
!
!--                    Take only values from horizontally-upward facing surfaces.
                       DO  mm = surf_def%start_index(j,i), surf_def%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_def%qsws(mm) *                  &
                                                              waterflux_output_conversion(k),      &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_def%upward(mm) )
                       ENDDO
                       DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%qsws(mm) *                  &
                                                              waterflux_output_conversion(k),      &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_lsm%upward(mm) )
                       ENDDO
                       DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_usm%qsws(mm) *                  &
                                                              waterflux_output_conversion(k),      &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_usm%upward(mm) )
                       ENDDO
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'hfss' ) ! surface sensible heat flux
                 IF ( .NOT. neutral )  THEN
                    DO  m = 1, vmea(l)%ns
!
!--                    Surface data is only available on inner subdomains, not on ghost points.
!--                    Hence, limit the indices.
                       j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                       j = MERGE( j           , nyn, j            < nyn )
                       i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                       i = MERGE( i           , nxr, i            < nxr )
                       k = vmea(l)%k(m)
!
!--                    Take only values from horizontally-upward facing surfaces.
                       DO  mm = surf_def%start_index(j,i), surf_def%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_def%shf(mm) *                   &
                                                              heatflux_output_conversion(k),       &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_def%upward(mm) )
                       ENDDO
                       DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%shf(mm) *                   &
                                                              heatflux_output_conversion(k),       &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_lsm%upward(mm) )
                       ENDDO
                       DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                          vmea(l)%measured_vars(m,n) = MERGE( surf_usm%shf(mm) *                   &
                                                              heatflux_output_conversion(k),       &
                                                              vmea(l)%measured_vars(m,n),          &
                                                              surf_usm%upward(mm) )
                       ENDDO
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'rnds' ) ! surface net radiation
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points.
!--                 Hence, limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%rad_net(mm),                   &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%rad_net(mm),                   &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.
!
!--           Surface shortwave out. Note surface flux and flux in the air are considered the same
!--           at the moment. The air flux is not available from the RTM.
              CASE ( 'rsus', 'rsu' )
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points.
!--                 Hence, limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%rad_sw_out(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%rad_sw_out(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.
!
!--           Surface shortwave in. Note surface flux and flux in the air are considered the same
!--           at the moment. The air flux is not available from the RTM.
              CASE ( 'rsds', 'rsd' )
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points.
!--                 Hence, limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%rad_sw_in(mm),                 &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%rad_sw_in(mm),                 &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.
!
!--           Surface longwave out. Note surface flux and flux in the air are considered the same
!--           at the moment. The air flux is not available from the RTM.
              CASE ( 'rlus', 'rlu' )
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points.
!--                 Hence, limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%rad_lw_out(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%rad_lw_out(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.
!
!--           Surface longwave in. Note surface flux and flux in the air are considered the same
!--           at the moment. The air flux is not available from the RTM.
              CASE ( 'rlds', 'rld' )
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points.
!--                 Hence, limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%rad_lw_in(mm),                 &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%rad_lw_in(mm),                 &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'ts', 'tb' ) ! surface temperature and brighness temperature
                 DO  m = 1, vmea(l)%ns
!
!--                 Surface data is only available on inner subdomains, not on ghost points. Hence,
!--                 limit the indices.
                    j = MERGE( vmea(l)%j(m), nys, vmea(l)%j(m) > nys )
                    j = MERGE( j           , nyn, j            < nyn )
                    i = MERGE( vmea(l)%i(m), nxl, vmea(l)%i(m) > nxl )
                    i = MERGE( i           , nxr, i            < nxr )
!
!--                 Take only values from horizontally-upward facing surfaces.
                    DO  mm = surf_def%start_index(j,i), surf_def%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_def%pt_surface(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_def%upward(mm) )
                    ENDDO
                    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_lsm%pt_surface(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_lsm%upward(mm) )
                    ENDDO
                    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                       vmea(l)%measured_vars(m,n) = MERGE( surf_usm%pt_surface(mm),                &
                                                           vmea(l)%measured_vars(m,n),             &
                                                           surf_usm%upward(mm) )
                    ENDDO
                 ENDDO
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE ( 'lwp' ) ! liquid water path
                 IF ( ASSOCIATED( ql ) )  THEN
                    DO  m = 1, vmea(l)%ns
                       j = vmea(l)%j(m)
                       i = vmea(l)%i(m)

                       vmea(l)%measured_vars(m,n) = SUM( ql(nzb:nzt,j,i) * dzw(1:nzt+1) )          &
                                                    * rho_surface
                    ENDDO
                    vmea(l)%var_atts(n)%sampled = .TRUE.
                 ENDIF

              CASE ( 'ps' ) ! surface pressure
                 vmea(l)%measured_vars(:,n) = surface_pressure
                 vmea(l)%var_atts(n)%sampled = .TRUE.

              CASE DEFAULT
!
!--              Check for chemistry variables. Since the list of possible chemistry variables
!--              depends on the chosen mechanism, no fixed pre-defined list of samplable variables
!--              can be defined, so it cannot be well considered within the CASE statements.
                 IF ( air_chemistry )  THEN
!
!--                 First, search for the measured variable in the chem_vars list, in order to
!--                 obtain the internal name of the variable.
                    ind_chem = -999
                    DO  nn = 1, UBOUND( chem_vars, 2 )
                       IF ( TRIM( vmea(l)%var_atts(n)%name ) == TRIM( chem_vars(0,nn) ) )          &
                          ind_chem = nn
                    ENDDO
!
!--                 Run loop over all chemical species, if the measured variable matches the interal
!--                 name, sample the variable. Note, nvar as a chemistry-module variable. Also
!--                 check if ind_chem lies in an allowed range.
                    IF ( ind_chem >= LBOUND( chem_vars(0,:), DIM = 1 )  .AND.                      &
                         ind_chem <= UBOUND( chem_vars(0,:), DIM = 1 ) )  THEN
                       DO  nn = 1, nvar
                          IF ( TRIM( chem_vars(1,ind_chem) ) == TRIM( chem_species(nn)%name ) )    &
                          THEN
                             DO  m = 1, vmea(l)%ns
                                k = vmea(l)%k(m)
                                j = vmea(l)%j(m)
                                i = vmea(l)%i(m)
                                vmea(l)%measured_vars(m,n) = chem_species(nn)%conc(k,j,i)
                             ENDDO
                             vmea(l)%var_atts(n)%sampled = .TRUE.
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDIF
!
!--              Check modules for matching variables.
!--              This is needs to be done outside of the module-interface framework.
!--              Since virtual_measurements_mod is already embedded into the module-interface and
!--              thus module_interface_mod already depends on virtual_measurements_mod, a call to
!--              a wrapper routine located in module_interface_mod is not possible here because
!--              this would create a circular dependency.
                 IF ( biometeorology )  CALL bio_vm_sampling( TRIM( vmea(l)%var_atts(n)%name ),    &
                                                              vmea(l)%measured_vars(:,n),          &
                                                              vmea(l)%i, vmea(l)%j, vmea(l)%k,     &
                                                              vmea(l)%ns,                          &
                                                              vmea(l)%measured_vars_soil(:,n),     &
                                                              vmea(l)%i_soil, vmea(l)%j_soil,      &
                                                              vmea(l)%k_soil, vmea(l)%ns_soil,     &
                                                              vmea(l)%var_atts(n)%sampled )

                 IF ( land_surface  )  CALL lsm_vm_sampling( TRIM( vmea(l)%var_atts(n)%name ),     &
                                                             vmea(l)%measured_vars(:,n),           &
                                                             vmea(l)%i, vmea(l)%j, vmea(l)%k,      &
                                                             vmea(l)%ns,                           &
                                                             vmea(l)%measured_vars_soil(:,n),      &
                                                             vmea(l)%i_soil, vmea(l)%j_soil,       &
                                                             vmea(l)%k_soil, vmea(l)%ns_soil,      &
                                                             vmea(l)%var_atts(n)%sampled )

                 IF ( radiation )  CALL radiation_vm_sampling( TRIM( vmea(l)%var_atts(n)%name ),   &
                                                               vmea(l)%measured_vars(:,n),         &
                                                               vmea(l)%i, vmea(l)%j, vmea(l)%k,    &
                                                               vmea(l)%ns,                         &
                                                               vmea(l)%measured_vars_soil(:,n),    &
                                                               vmea(l)%i_soil, vmea(l)%j_soil,     &
                                                               vmea(l)%k_soil, vmea(l)%ns_soil,    &
                                                               vmea(l)%var_atts(n)%sampled )

                 IF ( urban_surface )  CALL usm_vm_sampling( TRIM( vmea(l)%var_atts(n)%name ),     &
                                                             vmea(l)%measured_vars(:,n),           &
                                                             vmea(l)%i, vmea(l)%j, vmea(l)%k,      &
                                                             vmea(l)%ns,                           &
                                                             vmea(l)%measured_vars_soil(:,n),      &
                                                             vmea(l)%i_soil, vmea(l)%j_soil,       &
                                                             vmea(l)%k_soil, vmea(l)%ns_soil,      &
                                                             vmea(l)%var_atts(n)%sampled )

           END SELECT

        ENDDO

     ENDDO

     CALL cpu_log( log_point_s(27), 'VM sampling', 'stop' )

 END SUBROUTINE vm_sampling

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Bi-linear interpolation onto sampling location.
!--------------------------------------------------------------------------------------------------!
 FUNCTION interpolate_to_sampling_location( var, k, j, i, off_y, off_x, z, coords )                &
          RESULT( var_int )

    INTEGER(iwp) ::  i   !< nearest grid index of sampling location in x
    INTEGER(iwp) ::  ii  !< interpolation grid index along x-dimension
    INTEGER(iwp) ::  j   !< nearest grid index of sampling location in y
    INTEGER(iwp) ::  jj  !< interpolation grid index along y-dimension
    INTEGER(iwp) ::  k   !< nearest grid index of sampling location in z
    INTEGER(iwp) ::  kk  !< lower interpolation grid index along z-dimension
    INTEGER(iwp) ::  kkp !< upper interpolation grid index along z-dimension

    REAL(wp) ::  off_x    !< offset value in x to indicate variable position on staggered grid
    REAL(wp) ::  off_y    !< offset value in y to indicate variable position on staggered grid
    REAL(wp) ::  var_int  !< interpolated value
    REAL(wp) ::  vl       !< interpolated value at kk-level
    REAL(wp) ::  vly1     !< interpolated value along x at kk and jj
    REAL(wp) ::  vly2     !< interpolated value along x at kk and jj+1
    REAL(wp) ::  vu       !< interpolated value at kkp-level
    REAL(wp) ::  vuy1     !< interpolated value along x at kkp and jj
    REAL(wp) ::  vuy2     !< interpolated value along x at kkp and jj+1

    REAL(wp), DIMENSION(3) ::  coords     !< sampling coordinates

    REAL(wp), DIMENSION(nzb:nzt+1) ::  z  !< z-dimension of passed variable

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var !< passed variable

!
!-- Depending on the sampling point, determine the interpolation grid-point indices.
!-- In 2D, for instance, if the sampling point is located in the south-left corner of grid
!-- point (j,i), the interpolated variable is obtained from grid point (j-1,j,i-1,i).
!-- If it is in the north-right corner, it will be obtained from (j,j+1,i,i+1).
    IF ( ( i + off_x ) * dx - coords(1) >= 0.0_wp  )  THEN
       ii = i - 1
    ELSE
       ii = i
    ENDIF

    IF ( ( j + off_y ) * dy - coords(2) >= 0.0_wp )  THEN
       jj = j - 1
    ELSE
       jj = j
    ENDIF
!
!-- Determine the vertical interpolation indices and limit them to values between nzb and nzt+1.
    IF ( z(k) - coords(3) >= 0.0_wp )  THEN
       kk = MAX( nzb, k - 1 )
    ELSE
       kk = k
    ENDIF

    IF ( z(k) - coords(3) >= 0.0_wp )  THEN
       kkp = k
    ELSE
       kkp = MIN( k + 1, nzt + 1 )
    ENDIF

!
!-- Bi-linearly interpolate the required variable onto x-y sampling location at discrete
!-- height levels kk and kkp (which is kk+1). Therefore, first interpolate the variable along
!-- x at its discrete y-locations jj and (jj+1). In a second step interpolate onto the x-y
!-- sampling location along y.
    vly1 = var(kk,jj,ii)   + ( var(kk,jj,ii+1)   - var(kk,jj,ii)   ) * ddx                         &
                           * ( coords(1) - ( ii + off_x ) * dx )
    vly2 = var(kk,jj+1,ii) + ( var(kk,jj+1,ii+1) - var(kk,jj+1,ii) ) * ddx                         &
                           * ( coords(1) - ( ii + off_x ) * dx )
    vl   = vly1 + ( vly2 - vly1 ) * ddy * ( coords(2) - ( jj + off_y ) * dy )

    vuy1 = var(kkp,jj,ii)   + ( var(kkp,jj,ii+1)   - var(kkp,jj,ii)   ) * ddx                      &
                           * ( coords(1) - ( ii + off_x ) * dx )
    vuy2 = var(kkp,jj+1,ii) + ( var(kkp,jj+1,ii+1) - var(kkp,jj+1,ii) ) * ddx                      &
                           * ( coords(1) - ( ii + off_x ) * dx )
    vu   = vuy1 + ( vuy2 - vuy1 ) * ddy * ( coords(2) - ( jj + off_y ) * dy )
!
!-- Now interpolate linearly onto the z-position.
    var_int = vl + ( vu - vl ) / ( z(kkp+1) - z(kk) ) * ( coords(3) - z(kk) )

 END FUNCTION interpolate_to_sampling_location


 END MODULE virtual_measurement_mod
