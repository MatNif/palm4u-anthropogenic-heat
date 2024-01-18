!> @file src/inifor_grid.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2017-2021 Leibniz Universitaet Hannover
! Copyright 2017-2021 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Authors:
! --------
!> @author Eckhard Kadasch (Deutscher Wetterdienst, Offenbach)
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> The grid module contains all variables and routines to specify and work with
!> the numerical grids in INIFOR. By convention, all angles are stored in
!> radians.
!------------------------------------------------------------------------------!

#if defined ( __netcdf )
 MODULE inifor_grid

    USE inifor_control
    USE inifor_defs,                                                           &
        ONLY:  CFG_INIT_PROFILE, CFG_INIT_SOIL_VOLUME, CFG_INIT_SOIL_PROFILE,  &
               CFG_FORCING_HETERO, CFG_FORCING_HOMO, CFG_FORCING_NUDGING,      &
               DATE, EARTH_RADIUS, TO_RADIANS, TO_DEGREES, PI,                 &
               SNAME, LNAME, PATH, FORCING_STEP, FILL_ITERATIONS,              &
               BETA, P_SL, T_SL, BETA, RD, RV, G, P_REF, RD_PALM, CP_PALM,     &
               RHO_L, OMEGA, HECTO, wp, iwp,                                   &
               PIDS_ORIGIN_LON,                                                &
               PIDS_ORIGIN_LAT,                                                &
               PIDS_ORIGIN_Z
    USE inifor_io,                                                             &
        ONLY:  file_is_present, get_cosmo_grid, get_input_file_list,           &
               get_netcdf_attribute, get_netcdf_dim_vector,                    &
               get_netcdf_variable, set_palm_origin,                           &
               netcdf_variable_present_in_file, parse_command_line_arguments,  &
               validate_config, validate_dataset, has_surface_value
    USE inifor_transform,                                                      &
        ONLY:  average_2d, centre_velocities,                                  &
               compute_horizontal_interp_weights, compute_sigma,               &
               fill_water_cells, find_horizontal_neighbours,                   &
               find_vertical_neighbours_and_weights_interp,                    &
               find_vertical_neighbours_and_weights_average,                   &
               gamma_from_hemisphere, interpolate_2d, lamc_to_lamn,            &
               map_terrain_driver,                                             &
               phic_to_phin, phirot2phi, phi2phirot, project,                  &
               rla2rlarot, rlarot2rla, rotate_to_cosmo, uv2uvrot
    USE inifor_types
    USE inifor_util
    USE netcdf,                                                                &
        ONLY:  NF90_MAX_NAME, NF90_MAX_VAR_DIMS
    
    IMPLICIT NONE
    
    SAVE
    
    REAL(wp) ::  averaging_angle   = 0.0_wp       !< latitudal and longitudal width of averaging regions [rad]
    REAL(wp) ::  averaging_width_ns = 0.0_wp      !< longitudal width of averaging regions [m]
    REAL(wp) ::  averaging_width_ew = 0.0_wp      !< latitudal width of averaging regions [m]
    REAL(wp) ::  blending_dz       = 250.0_wp     !< height above mesoscale terrain where S blending function for the terrain mapping starts [m]
    REAL(wp) ::  blending_z_ubound = 1000.0_wp    !< height above sea level where terrain mapping has no effect (end of S blending function) [m]
    REAL(wp) ::  phi_equat         = 0.0_wp       !< latitude of rotated equator of COSMO-DE grid [rad]
    REAL(wp) ::  phi_n             = 0.0_wp       !< latitude of rotated pole of COSMO-DE grid [rad]
    REAL(wp) ::  lambda_n          = 0.0_wp       !< longitude of rotaded pole of COSMO-DE grid [rad]
    REAL(wp) ::  phi_c             = 0.0_wp       !< rotated-grid latitude of the center of the PALM domain [rad]
    REAL(wp) ::  lambda_c          = 0.0_wp       !< rotated-grid longitude of the centre of the PALM domain [rad]
    REAL(wp) ::  phi_cn            = 0.0_wp       !< latitude of the rotated pole relative to the COSMO-DE grid [rad]
    REAL(wp) ::  lambda_cn         = 0.0_wp       !< longitude of the rotated pole relative to the COSMO-DE grid [rad]
    REAL(wp) ::  lam_centre        = 0.0_wp       !< longitude of the PLAM domain centre in the source (COSMO rotated-pole) system [rad]
    REAL(wp) ::  phi_centre        = 0.0_wp       !< latitude of the PLAM domain centre in the source (COSMO rotated-pole) system [rad]
    REAL(wp) ::  lam_east          = 0.0_wp       !< longitude of the east central-averaging-domain boundary in the source (COSMO rotated-pole) system [rad]
    REAL(wp) ::  lam_west          = 0.0_wp       !< longitude of the west central-averaging-domain boundary in the source (COSMO rotated-pole) system [rad]
    REAL(wp) ::  phi_north         = 0.0_wp       !< latitude of the north central-averaging-domain boundary in the source (COSMO rotated-pole) system [rad]
    REAL(wp) ::  phi_south         = 0.0_wp       !< latitude of the south central-averaging-domain boundary in the source (COSMO rotated-pole) system [rad]
    REAL(wp) ::  gam               = 0.0_wp       !< angle for working around phirot2phi/rlarot2rla bug
    REAL(wp) ::  dx                = 0.0_wp       !< PALM-4U grid spacing in x direction [m]
    REAL(wp) ::  dy                = 0.0_wp       !< PALM-4U grid spacing in y direction [m]
    REAL(wp) ::  dz(10)            = -1.0_wp      !< PALM-4U grid spacing in z direction [m]
    REAL(wp) ::  dz_max            = 1000.0_wp    !< maximum vertical grid spacing [m]
    REAL(wp) ::  dz_stretch_factor = 1.08_wp      !< factor for vertical grid stretching [m]
    REAL(wp) ::  dz_stretch_level  = -9999999.9_wp!< height above which the vertical grid will be stretched [m]
    REAL(wp) ::  dz_stretch_level_start(9) = -9999999.9_wp !< namelist parameter
    REAL(wp) ::  dz_stretch_level_end(9) = 9999999.9_wp !< namelist parameter
    REAL(wp) ::  dz_stretch_factor_array(9) = 1.08_wp !< namelist parameter
    REAL(wp) ::  dxi               = 0.0_wp       !< inverse PALM-4U grid spacing in x direction [m^-1]
    REAL(wp) ::  dyi               = 0.0_wp       !< inverse PALM-4U grid spacing in y direction [m^-1]
    REAL(wp) ::  dzi               = 0.0_wp       !< inverse PALM-4U grid spacing in z direction [m^-1]
    REAL(wp) ::  f3                = 0.0_wp       !< Coriolis parameter
    REAL(wp) ::  lx                = 0.0_wp       !< PALM-4U domain size in x direction [m]
    REAL(wp) ::  ly                = 0.0_wp       !< PALM-4U domain size in y direction [m]
    REAL(wp) ::  p0                = 0.0_wp       !< PALM-4U surface pressure, at z0 [Pa]
    REAL(wp) ::  x0                = 0.0_wp       !< x coordinate of PALM-4U Earth tangent [m] 
    REAL(wp) ::  y0                = 0.0_wp       !< y coordinate of PALM-4U Earth tangent [m] 
    REAL(wp) ::  z0                = 0.0_wp       !< Elevation of the PALM-4U domain above sea level [m]
    REAL(wp) ::  z_top(1)          = 0.0_wp       !< height of the scalar top boundary [m]
    REAL(wp) ::  zw_top(1)         = 0.0_wp       !< height of the vertical velocity top boundary [m]
    REAL(wp) ::  lonmin_cosmo      = 0.0_wp       !< Minimunm longitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp) ::  lonmax_cosmo      = 0.0_wp       !< Maximum longitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp) ::  latmin_cosmo      = 0.0_wp       !< Minimunm latitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp) ::  latmax_cosmo      = 0.0_wp       !< Maximum latitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp) ::  lonmin_palm       = 0.0_wp       !< Minimunm longitude of PALM grid [COSMO rotated-pole rad]
    REAL(wp) ::  lonmax_palm       = 0.0_wp       !< Maximum longitude of PALM grid [COSMO rotated-pole rad]
    REAL(wp) ::  latmin_palm       = 0.0_wp       !< Minimunm latitude of PALM grid [COSMO rotated-pole rad]
    REAL(wp) ::  latmax_palm       = 0.0_wp       !< Maximum latitude of PALM grid [COSMO rotated-pole rad]
    REAL(wp) ::  lonmin_tot        = 0.0_wp       !< Minimunm longitude of required COSMO data [COSMO rotated-pole rad]
    REAL(wp) ::  lonmax_tot        = 0.0_wp       !< Maximum longitude of required COSMO data [COSMO rotated-pole rad]
    REAL(wp) ::  latmin_tot        = 0.0_wp       !< Minimunm latitude of required COSMO data [COSMO rotated-pole rad]
    REAL(wp) ::  latmax_tot        = 0.0_wp       !< Maximum latitude of required COSMO data [COSMO rotated-pole rad]
    REAL(wp) ::  latitude          = 0.0_wp       !< geographical latitude of the PALM-4U origin, from inipar namelist [deg]
    REAL(wp) ::  longitude         = 0.0_wp       !< geographical longitude of the PALM-4U origin, from inipar namelist [deg]
    REAL(wp) ::  origin_lat        = 0.0_wp       !< geographical latitude of the PALM-4U origin, from static driver netCDF file [deg]
    REAL(wp) ::  origin_lon        = 0.0_wp       !< geographical longitude of the PALM-4U origin, from static driver netCDF file [deg]
    REAL(wp) ::  rotation_angle    = 0.0_wp       !< clockwise angle the PALM-4U north is rotated away from geographical north [deg]
    REAL(wp) ::  end_time          = 0.0_wp       !< PALM-4U simulation time [s]

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  hhl             !< heights of half layers (cell faces) above sea level in COSMO-DE, read in from external file
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  hfl             !< heights of full layers (cell centres) above sea level in COSMO-DE, computed as arithmetic average of hhl
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  zs              !< PALM terrain height ('z surface') above sea level at scalar points  [m]
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  zs_u            !< PALM terrain height ('z surface') above sea level at u points  [m]
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  zs_v            !< PALM terrain height ('z surface') above sea level at v points  [m]
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  depths          !< COSMO-DE's TERRA-ML soil layer depths
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  d_depth         !< COSMO-DE's TERRA-ML soil layer thicknesses
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  d_depth_rho_inv !< inverted soil water mass
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  rlon            !< longitudes of COSMO-DE's rotated-pole grid
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  rlat            !< latitudes of COSMO-DE's rotated-pole grid
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  time            !< output times 
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  x               !< base palm grid x coordinate vector pointed to by grid_definitions
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  xu              !< base palm grid xu coordinate vector pointed to by grid_definitions
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  y               !< base palm grid y coordinate vector pointed to by grid_definitions
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  yv              !< base palm grid yv coordinate vector pointed to by grid_definitions
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  z_column        !< base palm grid z coordinate vector including the top boundary coordinate (entire column)
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  zw_column       !< base palm grid zw coordinate vector including the top boundary coordinate (entire column)
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  z               !< base palm grid z coordinate vector pointed to by grid_definitions
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET     ::  zw              !< base palm grid zw coordinate vector pointed to by grid_definitions

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  soiltyp     !< COSMO-DE soil type map
    INTEGER(iwp) ::  dz_stretch_level_end_index(9)               !< vertical grid level index until which the vertical grid spacing is stretched
    INTEGER(iwp) ::  dz_stretch_level_start_index(9)             !< vertical grid level index above which the vertical grid spacing is stretched
    INTEGER(iwp) ::  iostat !< return status of READ statement
    INTEGER(iwp) ::  nt    !< number of output time steps
    INTEGER(iwp) ::  nx    !< number of PALM-4U grid points in x direction
    INTEGER(iwp) ::  ny    !< number of PALM-4U grid points in y direction
    INTEGER(iwp) ::  nz    !< number of PALM-4U grid points in z direction
    INTEGER(iwp) ::  nlon  !< number of longitudal points in target grid (COSMO-DE)
    INTEGER(iwp) ::  nlat  !< number of latitudal points in target grid (COSMO-DE)
    INTEGER(iwp) ::  nlev  !< number of levels in target grid (COSMO-DE)
    INTEGER(iwp) ::  ndepths !< number of COSMO-DE soil layers
    INTEGER(iwp) ::  start_hour_flow          !< start of flow forcing in number of hours relative to start_date
    INTEGER(iwp) ::  start_hour_radiation     !< start of radiation forcing in number of hours relative to start_date, 0 to 2 hours before start_hour_flow to reconstruct hourly averages from one- to three hourly averages of the input data
    INTEGER(iwp) ::  start_hour_precipitation !< start of forcing for precipitaiton forcing in number of hours relative to start_date
    INTEGER(iwp) ::  end_hour  !< simulation time in hours
    INTEGER(iwp) ::  step_hour !< number of hours between forcing time steps

    LOGICAL ::  init_variables_required       !< flag controlling whether init variables are to be processed
    LOGICAL ::  boundary_variables_required   !< flag controlling whether boundary grids are to be allocated and boundary variables are to be computed
    LOGICAL ::  ls_forcing_variables_required !< flag controlling whether large-scale forcing variables are to be computed
    LOGICAL ::  surface_forcing_required      !< flag controlling whether surface forcing variables are to be computed
    LOGICAL ::  palm_domain_outside_cosmo     !< indicates whether COSMO grid covers the PALM domain and the geostrophic averaging domains

    TYPE(nc_var), ALLOCATABLE, TARGET ::  input_var_table(:)  !< table of input variables
    TYPE(nc_var), ALLOCATABLE, TARGET ::  output_var_table(:) !< table of input variables
    TYPE(nc_var) ::  dummy_var                                !< dummy variable, used for reading HHL, rlon, rlat, and PALM terrain

    TYPE(grid_definition), TARGET ::  palm_grid                       !< PALM-4U grid in the target system (COSMO-DE rotated-pole)
    TYPE(grid_definition), TARGET ::  palm_intermediate               !< PALM-4U grid with coarse vertical grid wiht levels interpolated from COSMO-DE grid
    TYPE(grid_definition), TARGET ::  cosmo_grid                      !< target system (COSMO-DE rotated-pole)
    TYPE(grid_definition), TARGET ::  scalars_east_grid               !< grid for eastern scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_west_grid               !< grid for western scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_north_grid              !< grid for northern scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_south_grid              !< grid for southern scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_top_grid                !< grid for top scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_east_intermediate       !< intermediate grid for eastern scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_west_intermediate       !< intermediate grid for western scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_north_intermediate      !< intermediate grid for northern scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_south_intermediate      !< intermediate grid for southern scalar boundary condition
    TYPE(grid_definition), TARGET ::  scalars_top_intermediate        !< intermediate grid for top scalar boundary condition
    TYPE(grid_definition), TARGET ::  u_initial_grid                  !< grid for u initial condition
    TYPE(grid_definition), TARGET ::  u_east_grid                     !< grid for eastern u boundary condition
    TYPE(grid_definition), TARGET ::  u_west_grid                     !< grid for western u boundary condition
    TYPE(grid_definition), TARGET ::  u_north_grid                    !< grid for northern u boundary condition
    TYPE(grid_definition), TARGET ::  u_south_grid                    !< grid for southern u boundary condition
    TYPE(grid_definition), TARGET ::  u_top_grid                      !< grid for top u boundary condition
    TYPE(grid_definition), TARGET ::  u_initial_intermediate          !< intermediate grid for u initial condition
    TYPE(grid_definition), TARGET ::  u_east_intermediate             !< intermediate grid for eastern u boundary condition
    TYPE(grid_definition), TARGET ::  u_west_intermediate             !< intermediate grid for western u boundary condition
    TYPE(grid_definition), TARGET ::  u_north_intermediate            !< intermediate grid for northern u boundary condition
    TYPE(grid_definition), TARGET ::  u_south_intermediate            !< intermediate grid for southern u boundary condition
    TYPE(grid_definition), TARGET ::  u_top_intermediate              !< intermediate grid for top u boundary condition
    TYPE(grid_definition), TARGET ::  v_initial_grid                  !< grid for v initial condition
    TYPE(grid_definition), TARGET ::  v_east_grid                     !< grid for eastern v boundary condition
    TYPE(grid_definition), TARGET ::  v_west_grid                     !< grid for western v boundary condition
    TYPE(grid_definition), TARGET ::  v_north_grid                    !< grid for northern v boundary condition
    TYPE(grid_definition), TARGET ::  v_south_grid                    !< grid for southern v boundary condition
    TYPE(grid_definition), TARGET ::  v_top_grid                      !< grid for top v boundary condition
    TYPE(grid_definition), TARGET ::  v_initial_intermediate          !< intermediate grid for v initial condition
    TYPE(grid_definition), TARGET ::  v_east_intermediate             !< intermediate grid for eastern v boundary condition
    TYPE(grid_definition), TARGET ::  v_west_intermediate             !< intermediate grid for western v boundary condition
    TYPE(grid_definition), TARGET ::  v_north_intermediate            !< intermediate grid for northern v boundary condition
    TYPE(grid_definition), TARGET ::  v_south_intermediate            !< intermediate grid for southern v boundary condition
    TYPE(grid_definition), TARGET ::  v_top_intermediate              !< intermediate grid for top v boundary condition
    TYPE(grid_definition), TARGET ::  w_initial_grid                  !< grid for w initial condition
    TYPE(grid_definition), TARGET ::  w_east_grid                     !< grid for eastern w boundary condition
    TYPE(grid_definition), TARGET ::  w_west_grid                     !< grid for western w boundary condition
    TYPE(grid_definition), TARGET ::  w_north_grid                    !< grid for northern w boundary condition
    TYPE(grid_definition), TARGET ::  w_south_grid                    !< grid for southern w boundary condition
    TYPE(grid_definition), TARGET ::  w_top_grid                      !< grid for top w boundary condition
    TYPE(grid_definition), TARGET ::  w_initial_intermediate          !< intermediate grid for w initial condition
    TYPE(grid_definition), TARGET ::  w_east_intermediate             !< intermediate grid for eastern w boundary condition
    TYPE(grid_definition), TARGET ::  w_west_intermediate             !< intermediate grid for western w boundary condition
    TYPE(grid_definition), TARGET ::  w_north_intermediate            !< intermediate grid for northern w boundary condition
    TYPE(grid_definition), TARGET ::  w_south_intermediate            !< intermediate grid for southern w boundary condition
    TYPE(grid_definition), TARGET ::  w_top_intermediate              !< intermediate grid for top w boundary condition
    TYPE(grid_definition), TARGET ::  north_geostrophic_scalar_profile!< grid of the northern geostrophic scalar averaging region
    TYPE(grid_definition), TARGET ::  south_geostrophic_scalar_profile!< grid of the southern geostrophic scalar averaging region
    TYPE(grid_definition), TARGET ::  west_geostrophic_scalar_profile !< grid of the western geostrophic scalar averaging region
    TYPE(grid_definition), TARGET ::  east_geostrophic_scalar_profile !< grid of the eastern geostrophic scalar averaging region
    TYPE(grid_definition), TARGET ::  geostrophic_scalar_profile      !< grid of the central geostrophic scalar averaging region
    TYPE(grid_definition), TARGET ::  geostrophic_w_profile           !< grid of the central geostrophic w-velocity averaging region
    TYPE(grid_definition), TARGET ::  averaged_soil_profile           !< averaging grid for initial soil profiles
    TYPE(grid_definition), TARGET ::  averaged_scalar_profile         !< averaging grid for initial and boundary condition scalar profiles
    TYPE(grid_definition), TARGET ::  averaged_w_profile              !< averaging grid for initial and boundary condition scalar profiles
    TYPE(grid_definition), TARGET ::  averaged_scalar_top_point       !< averaging grid for top scalar boundary conditions for homogeneous forcing mode
    TYPE(grid_definition), TARGET ::  averaged_w_top_point            !< averaging grid for top w boundary condition for homogeneous forcing mode

    TYPE(io_group), ALLOCATABLE, TARGET ::  io_group_list(:)  !< List of I/O groups, which group together output variables that share the same input variable
 
    NAMELIST /inipar/ nx, ny, nz, dx, dy, dz, longitude, latitude,             &
                      dz_max, dz_stretch_factor, dz_stretch_level,             &
                      dz_stretch_level_start, dz_stretch_level_end,            &
                      blending_dz, blending_z_ubound
    NAMELIST /d3par/  end_time
    
    CHARACTER(LEN=LNAME) ::  nc_source_text = ''  !< Text describing the source of the output data, e.g. 'COSMO-DE analysis from ...'

    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  flow_files      !< list of atmospheric input files (<prefix>YYYYMMDDHH-flow.nc)
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  precip_files    !< list of precipitation input files (<prefix>YYYYMMDDHH-precip.nc) 
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  soil_files      !< list of soil input files (temperature, moisture, <prefix>YYYYMMDDHH-soil.nc)
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  radiation_files !< list of radiation input files (<prefix>YYYYMMDDHH-rad.nc)

    CHARACTER(LEN=SNAME) ::  input_prefix          !< prefix of input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  flow_prefix           !< prefix of flow input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  soil_prefix           !< prefix of soil input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  radiation_prefix      !< prefix of radiation input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  precipitation_prefix  !< prefix of input files for precipitation forcing, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  flow_suffix           !< suffix of flow input files, e.g. 'flow'
    CHARACTER(LEN=SNAME) ::  soil_suffix           !< suffix of soil input files, e.g. 'soil'
    CHARACTER(LEN=SNAME) ::  radiation_suffix      !< suffix of radiation input files, e.g. 'radiation'
    CHARACTER(LEN=SNAME) ::  precipitation_suffix  !< suffix of input files for precipition forcing, e.g. 'precip'
                          
    TYPE(nc_file) ::  output_file !< metadata of the dynamic driver 

    TYPE(inifor_config) ::  cfg !< container of the INIFOR command-line configuration

 CONTAINS
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine initializes INIFOR. This includes parsing command-line options,
!> setting the names of the input and output files, reading INIFOR's namelist
!> file as well as reading in and setting grid parameters for defining
!> interplation grids later in setup_grids().
!------------------------------------------------------------------------------!
 SUBROUTINE setup_parameters()

!
!------------------------------------------------------------------------------
! Section 1: Define default parameters
!------------------------------------------------------------------------------
    cfg%start_date = '2013072100'

!
!-- Default coordinates of the PALM origin in the geographical reference system
    origin_lat = 52.325079_wp * TO_RADIANS ! south-west of Berlin, origin used for the Dec 2017 showcase simulation
    origin_lon = 13.082744_wp * TO_RADIANS
    cfg%z0 = 35.0_wp

!
!-- Default atmospheric parameters
    cfg%ug = 0.0_wp
    cfg%vg = 0.0_wp
    cfg%p0 = P_SL

!
!-- Parameters for file names
    start_hour_flow = 0
    start_hour_radiation = 0
    start_hour_precipitation = start_hour_flow 
    step_hour = FORCING_STEP

    input_prefix = 'laf'
    cfg%flow_prefix = input_prefix
    cfg%input_prefix = input_prefix
    cfg%soil_prefix = input_prefix
    cfg%radiation_prefix = input_prefix
    cfg%precipitation_prefix  = input_prefix

    flow_suffix = '-flow'
    soil_suffix = '-soil'
    radiation_suffix = '-rad'
    precipitation_suffix = '-precip'

    cfg%debug = .FALSE.
    cfg%averaging_angle = 2.0_wp
!
!------------------------------------------------------------------------------
! Section 2: Read command-line arguments, namelist, and grid parameters
!------------------------------------------------------------------------------

!
!-- Set default paths and modes
    cfg%input_path         = './'
    cfg%hhl_file           = ''
    cfg%soiltyp_file       = ''
    cfg%namelist_file      = './namelist'
    cfg%static_driver_file = ''
    cfg%output_file = './palm-4u-input.nc'
    cfg%ic_mode = CFG_INIT_PROFILE
    cfg%isc_mode = CFG_INIT_SOIL_VOLUME
    cfg%bc_mode = CFG_FORCING_HETERO
    cfg%averaging_mode = 'level'

!
!-- Overwrite defaults with user configuration
    CALL parse_command_line_arguments( cfg )
    CALL report('main_loop', 'Running INIFOR version ' // VERSION)

    flow_prefix = TRIM(cfg%input_prefix)
    radiation_prefix = TRIM(cfg%input_prefix)
    soil_prefix = TRIM(cfg%input_prefix)
    precipitation_prefix = TRIM(cfg%input_prefix)
    IF (cfg%flow_prefix_is_set)  flow_prefix = TRIM(cfg%flow_prefix)
    IF (cfg%radiation_prefix_is_set)  radiation_prefix = TRIM(cfg%radiation_prefix)
    IF (cfg%soil_prefix_is_set)  soil_prefix = TRIM(cfg%soil_prefix)
    IF (cfg%precipitation_prefix_is_set)  precipitation_prefix = TRIM(cfg%precipitation_prefix)

    output_file%name = cfg%output_file

    init_variables_required = .TRUE.
    boundary_variables_required = (                                            &
       ( TRIM( cfg%bc_mode ) == CFG_FORCING_HETERO )  .OR.                     &
       ( TRIM( cfg%bc_mode ) == CFG_FORCING_HOMO )                             &
    )
    ls_forcing_variables_required = TRIM( cfg%bc_mode ) == CFG_FORCING_NUDGING
    surface_forcing_required = .TRUE.

    IF ( ls_forcing_variables_required )  THEN
       message = "Averaging of large-scale forcing profiles " //            &
                 "has not been implemented, yet."
       CALL inifor_abort('setup_parameters', message)
    ENDIF

!
!-- Set default file paths, if not specified by user.
    CALL normalize_path(cfg%input_path)
    IF (TRIM(cfg%hhl_file) == '')  cfg%hhl_file = TRIM(cfg%input_path) // 'hhl.nc'
    IF (TRIM(cfg%soiltyp_file) == '')  cfg%soiltyp_file = TRIM(cfg%input_path) // 'soil.nc'

    CALL validate_config( cfg ) 

    CALL report('setup_parameters', "atmosphere initialization mode: " // TRIM( cfg%ic_mode ) )
    CALL report('setup_parameters', "      soil initialization mode: " // TRIM( cfg%isc_mode ) )
    CALL report('setup_parameters', "                  forcing mode: " // TRIM( cfg%bc_mode ) )
    CALL report('setup_parameters', "                averaging mode: " // TRIM( cfg%averaging_mode ) )
    CALL report('setup_parameters', "               averaging angle: " // TRIM( real_to_str( cfg%averaging_angle ) )// " deg" )
    CALL report_toggle('setup_parameters', "               terrain mapping:", cfg%map_terrain )
    CALL report('setup_parameters', "                     data path: " // TRIM( cfg%input_path ) )
    CALL report('setup_parameters', "                      hhl file: " // TRIM( cfg%hhl_file ) )
    CALL report('setup_parameters', "                  soiltyp file: " // TRIM( cfg%soiltyp_file ) )
    CALL report('setup_parameters', "                 namelist file: " // TRIM( cfg%namelist_file ) )
    CALL report('setup_parameters', "              output data file: " // TRIM( output_file%name ) )
    CALL report_toggle('setup_parameters', "                 precipitation:", cfg%process_precipitation)
    CALL report_toggle('setup_parameters', "                debugging mode:", cfg%debug )

    CALL log_runtime('time', 'init')
!
!-- Read in namelist parameters
    OPEN(10, FILE=cfg%namelist_file, STATUS='old', IOSTAT=iostat)
    IF ( iostat /= 0 )  THEN
       message = "Failed to open file '" //             &
                 TRIM( cfg%namelist_file ) // "'. "
       CALL inifor_abort( 'setup_parameters', message )
    ENDIF

    READ(10, NML=inipar, IOSTAT=iostat) ! nx, ny, nz, dx, dy, dz
    IF ( iostat > 0 )  THEN      
       message = "Failed to read namelist 'inipar' from file '" //             &
                 TRIM( cfg%namelist_file ) // "'. "
       CALL inifor_abort( 'setup_parameters', message )
       CLOSE(10)
    ENDIF

    READ(10, NML=d3par, IOSTAT=iostat)  ! end_time
    IF ( iostat > 0 )  THEN
       message = "Failed to read namelist 'd3par' from file '" //              &
                 TRIM( cfg%namelist_file ) // "'. "
       CALL inifor_abort( 'setup_parameters', message )
       CLOSE(10)
    ENDIF
    CLOSE(10)
    
    CALL log_runtime('time', 'read')

    end_hour = CEILING( end_time / 3600.0 * step_hour )

!
!-- Generate input file lists
    CALL get_input_file_list(                                                  &
       cfg%start_date, start_hour_flow, end_hour, step_hour,                   &
       cfg%input_path, flow_prefix, flow_suffix, flow_files)
    CALL get_input_file_list(                                                  &
       cfg%start_date, start_hour_flow, start_hour_flow, step_hour,            &
       cfg%input_path, soil_prefix, soil_suffix, soil_files)
    CALL get_input_file_list(                                                  &
       cfg%start_date, start_hour_radiation, end_hour, step_hour,              &
       cfg%input_path, radiation_prefix, radiation_suffix, radiation_files)
    CALL get_input_file_list(                                                  &
       cfg%start_date, start_hour_flow, end_hour, step_hour,                   &
       cfg%input_path, precipitation_prefix, precipitation_suffix, precip_files )

!
!------------------------------------------------------------------------------
! Section 3: Check for consistency
!------------------------------------------------------------------------------
    !CALL validate_dataset( flow_files, cfg%hhl_file )

!
!------------------------------------------------------------------------------
! Section 4: Compute additional parameters
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Section 4.1: COSMO-DE parameters
!------------------------------------------------------------------------------


    CALL log_runtime('time', 'init')
!
!-- Read COSMO soil type map
    dummy_var%name = 'SOILTYP'
    CALL get_netcdf_variable(cfg%soiltyp_file, dummy_var, soiltyp)

    CALL set_palm_origin(cfg, longitude, latitude, origin_lon, origin_lat, z0)
    p0 = cfg%p0

    IF ( cfg%map_terrain )  THEN
       dummy_var%name = 'zt'
       CALL get_netcdf_variable( cfg%static_driver_file, dummy_var, zs)
       zs(:,:,:) = zs(:,:,:) + z0

       ALLOCATE( zs_u(1:nx,  1:ny+1,1) )
       ALLOCATE( zs_v(1:nx+1,1:ny,  1) )
       zs_u(:,:,:) = 0.5_wp * ( zs(1:nx, :,    :) + zs(2:nx+1, :,      :) )
       zs_v(:,:,:) = 0.5_wp * ( zs(:,    1:ny, :) + zs(:,      2:ny+1, :) )
    ENDIF

    CALL log_runtime('time', 'read')

    CALL get_cosmo_grid( cfg%hhl_file, soil_files(1), rlon, rlat, hhl, hfl,    &
                         depths, d_depth, d_depth_rho_inv, phi_n, lambda_n,    &
                         phi_equat,                                            &
                         lonmin_cosmo, lonmax_cosmo,                           &
                         latmin_cosmo, latmax_cosmo,                           &
                         nlon, nlat, nlev, ndepths )


!------------------------------------------------------------------------------
! Section 4.2: PALM-4U parameters
!------------------------------------------------------------------------------
!
!-- PALM-4U domain extents
    lx = (nx+1) * dx
    ly = (ny+1) * dy
    
!
!-- PALM-4U point of Earth tangency
    x0 = 0.0_wp
    y0 = 0.0_wp

!
!-- time vector
    nt = CEILING(end_time / (step_hour * 3600.0_wp)) + 1
    ALLOCATE( time(nt) )
    CALL linspace(0.0_wp, 3600.0_wp * (nt-1), time)
    output_file%time => time
    CALL log_runtime('time', 'init')

!
!-- Convert the PALM-4U origin coordinates to COSMO's rotated-pole grid
    phi_c    = TO_RADIANS *                                                 &
               phi2phirot( origin_lat * TO_DEGREES, origin_lon * TO_DEGREES,&
                           phi_n * TO_DEGREES, lambda_n * TO_DEGREES )
    lambda_c = TO_RADIANS *                                                 &
               rla2rlarot( origin_lat * TO_DEGREES, origin_lon * TO_DEGREES,&
                           phi_n * TO_DEGREES, lambda_n * TO_DEGREES,     &
                           0.0_wp )

!
!-- Set gamma according to whether PALM domain is in the northern or southern
!-- hemisphere of the COSMO rotated-pole system. Gamma assumes either the
!-- value 0 or PI and is needed to work around around a bug in the
!-- rotated-pole coordinate transformations.
    gam = gamma_from_hemisphere(origin_lat, phi_equat)

!
!-- Compute the north pole of the rotated-pole grid centred at the PALM-4U
!-- domain centre. The resulting (phi_cn, lambda_cn) are coordinates in
!-- COSMO-DE's rotated-pole grid.
    phi_cn    = phic_to_phin(phi_c) 
    lambda_cn = lamc_to_lamn(phi_c, lambda_c) 

    message =   "PALM-4U origin:" // NEW_LINE('') // &
       "           lon (lambda) = " // &
       TRIM(real_to_str_f(origin_lon * TO_DEGREES)) // " deg"// NEW_LINE(' ') //&
       "           lat (phi   ) = " // &
       TRIM(real_to_str_f(origin_lat * TO_DEGREES)) // " deg (geographical)" // NEW_LINE(' ') //&
       "           lon (lambda) = " // &
       TRIM(real_to_str_f(lambda_c * TO_DEGREES)) // " deg" // NEW_LINE(' ') // &
       "           lat (phi   ) = " // &
       TRIM(real_to_str_f(phi_c * TO_DEGREES)) // " deg (COSMO-DE rotated-pole)"
    CALL report ('setup_parameters', message)

    message = "COSMO rotated north pole:" // NEW_LINE(' ') // &
       "           lon (lambda) = " // &
       TRIM(real_to_str_f(lambda_n * TO_DEGREES)) // " deg" // NEW_LINE(' ') //&
       "           lat (phi   ) = " // &
       TRIM(real_to_str_f(phi_n * TO_DEGREES)) // " deg (geographical)"
    CALL report ('setup_parameters', message)
       
    message = "North pole of the rotated palm system:" // NEW_LINE(' ') // &
       "           lon (lambda) = " // &
       TRIM(real_to_str_f(lambda_cn * TO_DEGREES)) // " deg" // NEW_LINE(' ') // &
       "           lat (phi   ) = " // &
       TRIM(real_to_str_f(phi_cn * TO_DEGREES)) // " deg (COSMO-DE rotated-pole)"
    CALL report ('setup_parameters', message)

    CALL log_runtime('time', 'comp')

!------------------------------------------------------------------------------
! Section 4.3: INIFOR averaging domains
!------------------------------------------------------------------------------

!
!-- Compute coordiantes of the PALM centre in the source (COSMO) system
    phi_centre = phirot2phi(                                                   &
       phirot = project(0.5_wp*ly, y0, EARTH_RADIUS) * TO_DEGREES,             &
       rlarot = project(0.5_wp*lx, x0, EARTH_RADIUS) * TO_DEGREES,             &
       polphi = phi_cn * TO_DEGREES,                                           &
       polgam = gam * TO_DEGREES                                               &
    ) * TO_RADIANS

    lam_centre = rlarot2rla(                                                   &
       phirot = project(0.5_wp*ly, y0, EARTH_RADIUS) * TO_DEGREES,             &
       rlarot = project(0.5_wp*lx, x0, EARTH_RADIUS) * TO_DEGREES,             &
       polphi = phi_cn * TO_DEGREES, pollam = lambda_cn * TO_DEGREES,          &
       polgam = gam * TO_DEGREES                                               &
    ) * TO_RADIANS

    message = "PALM-4U centre:" // NEW_LINE('') // &
       "           lon (lambda) = " // &
       TRIM(real_to_str_f(lam_centre * TO_DEGREES)) // " deg" // NEW_LINE(' ') // &
       "           lat (phi   ) = " // &
       TRIM(real_to_str_f(phi_centre * TO_DEGREES)) // " deg (COSMO-DE rotated-pole)"
    CALL report( 'setup_parameters', message )

!
!-- Compute boundaries of the central averaging box
    averaging_angle = cfg%averaging_angle * TO_RADIANS
    lam_east = lam_centre + 0.5_wp * averaging_angle
    lam_west = lam_centre - 0.5_wp * averaging_angle
    phi_north = phi_centre + 0.5_wp * averaging_angle
    phi_south = phi_centre - 0.5_wp * averaging_angle
    averaging_width_ew = averaging_angle * COS(phi_centre) * EARTH_RADIUS
    averaging_width_ns = averaging_angle * EARTH_RADIUS

!
!-- Coriolis parameter
    f3 = 2.0_wp * OMEGA * SIN(                                                 &
       TO_RADIANS*phirot2phi( phi_centre * TO_DEGREES, lam_centre * TO_DEGREES,&
                              phi_n * TO_DEGREES,                              &
                              gam * TO_DEGREES )                               &
    )

 END SUBROUTINE setup_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Defines the COSMO, PALM-4U, PALM-4U boundary grids, in particular their
!> coordinates and interpolation weights
!------------------------------------------------------------------------------!
 SUBROUTINE setup_grids()
    CHARACTER ::  interp_mode

!------------------------------------------------------------------------------
! Section 0: Define base PALM-4U grid coordinate vectors
!------------------------------------------------------------------------------
!
!-- palm x y z, we allocate the column to nz+1 in order to include the top
!-- scalar boundary. The interpolation grids will be associated with
!-- a shorter column that omits the top element.
    ALLOCATE( x(0:nx), y(0:ny), z(1:nz), z_column(1:nz+1) )
    CALL linspace(0.5_wp * dx, lx - 0.5_wp * dx, x)
    CALL linspace(0.5_wp * dy, ly - 0.5_wp * dy, y)
    CALL stretched_z(z_column, dz, dz_max=dz_max,                           &
                     dz_stretch_factor=dz_stretch_factor,                   &
                     dz_stretch_level=dz_stretch_level,                     &
                     dz_stretch_level_start=dz_stretch_level_start,         &
                     dz_stretch_level_end=dz_stretch_level_end,             &
                     dz_stretch_factor_array=dz_stretch_factor_array)
    z(1:nz)  = z_column(1:nz)
    z_top(:) = z_column(nz+1)

!
!-- palm xu yv zw, compared to the scalar grid, velocity coordinates
!-- contain one element less.
    ALLOCATE( xu(1:nx),  yv(1:ny), zw(1:nz-1), zw_column(1:nz))
    CALL linspace(dx, lx - dx, xu)
    CALL linspace(dy, ly - dy, yv)
    CALL midpoints(z_column, zw_column)
    zw(1:nz-1) = zw_column(1:nz-1)
    zw_top(:)  = zw_column(nz)


!------------------------------------------------------------------------------
! Section 1: Define initialization and boundary grids
!------------------------------------------------------------------------------
    CALL init_grid_definition('palm', grid=palm_grid,                       &
            xmin=0.0_wp, xmax=lx,                                           &
            ymin=0.0_wp, ymax=ly,                                           &
            x0=x0, y0=y0, z0=z0,                                            &
            nx=nx, ny=ny, nz=nz, z=z, zw=zw, ic_mode=cfg%ic_mode)

!
!-- Subtracting 1 because arrays will be allocated with nlon + 1 elements.
    CALL init_grid_definition('cosmo-de', grid=cosmo_grid,                  &
            xmin=lonmin_cosmo, xmax=lonmax_cosmo,                           &
            ymin=latmin_cosmo, ymax=latmax_cosmo,                           &
            x0=x0, y0=y0, z0=0.0_wp,                                        &
            nx=nlon-1, ny=nlat-1, nz=nlev-1)

!
!-- Define intermediate grid. This is the same as palm_grid except with a
!-- much coarser vertical grid. The vertical levels are interpolated in each
!-- PALM column from COSMO's secondary levels. The main levels are then
!-- computed as the averages of the bounding secondary levels.
    CALL init_grid_definition('palm intermediate', grid=palm_intermediate,  &
            xmin=0.0_wp, xmax=lx,                                           &
            ymin=0.0_wp, ymax=ly,                                           &
            x0=x0, y0=y0, z0=z0,                                            &
            nx=nx, ny=ny, nz=nlev-2)

    CALL init_grid_definition('boundary', grid=u_initial_grid,              &
            xmin = dx, xmax = lx - dx,                                      &
            ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
            x0=x0, y0=y0, z0 = z0,                                          &
            nx = nx-1, ny = ny, nz = nz,                                    &
            z=z, ic_mode=cfg%ic_mode)

    CALL init_grid_definition('boundary', grid=v_initial_grid,              &
            xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
            ymin = dy, ymax = ly - dy,                                      &
            x0=x0, y0=y0, z0 = z0,                                          &
            nx = nx, ny = ny-1, nz = nz,                                    &
            z=z, ic_mode=cfg%ic_mode)

    CALL init_grid_definition('boundary', grid=w_initial_grid,              &
            xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
            ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
            x0=x0, y0=y0, z0 = z0,                                          &
            nx = nx, ny = ny, nz = nz-1,                                    &
            z=zw, ic_mode=cfg%ic_mode)

    CALL init_grid_definition('boundary intermediate', grid=u_initial_intermediate,      &
            xmin = dx, xmax = lx - dx,                                      &
            ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
            x0=x0, y0=y0, z0 = z0,                                          &
            nx = nx-1, ny = ny, nz = nlev - 2)

    CALL init_grid_definition('boundary intermediate', grid=v_initial_intermediate,      &
            xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
            ymin = dy, ymax = ly - dy,                                      &
            x0=x0, y0=y0, z0 = z0,                                          &
            nx = nx, ny = ny-1, nz = nlev - 2)

    CALL init_grid_definition('boundary intermediate', grid=w_initial_intermediate,      &
            xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
            ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
            x0=x0, y0=y0, z0 = z0,                                          &
            nx = nx, ny = ny, nz = nlev - 1)

    IF (boundary_variables_required)  THEN
!
!------------------------------------------------------------------------------
! Section 2: Define PALM-4U boundary grids
!------------------------------------------------------------------------------
       CALL init_grid_definition('boundary', grid=scalars_east_grid,           &
               xmin = lx + 0.5_wp * dx, xmax = lx + 0.5_wp * dx,               &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=scalars_west_grid,           &
               xmin = -0.5_wp * dx, xmax = -0.5_wp * dx,                       &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=scalars_north_grid,          &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = ly + 0.5_wp * dy, ymax = ly + 0.5_wp * dy,               &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=scalars_south_grid,          &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin = -0.5_wp * dy, ymax = -0.5_wp * dy,                       &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=scalars_top_grid,            &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny, nz = 1_iwp, z=z_top)

       CALL init_grid_definition('boundary', grid=u_east_grid,                 &
               xmin = lx, xmax = lx,                                           &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=u_west_grid,                 &
               xmin = 0.0_wp, xmax = 0.0_wp,                                   &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=u_north_grid,                &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = ly + 0.5_wp * dy, ymax = ly + 0.5_wp * dy,               &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = 0_iwp, nz = nz, z=z)
    
       CALL init_grid_definition('boundary', grid=u_south_grid,                &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = -0.5_wp * dy, ymax = -0.5_wp * dy,                       &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = 0_iwp, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=u_top_grid,                  &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = ny, nz = 1_iwp, z=z_top)

       CALL init_grid_definition('boundary', grid=v_east_grid,                 &
               xmin = lx + 0.5_wp * dx, xmax = lx + 0.5_wp * dx,               &
               ymin = dy, ymax = ly - dy,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny-1, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=v_west_grid,                 &
               xmin = -0.5_wp * dx, xmax = -0.5_wp * dx,                       &
               ymin = dy, ymax = ly - dy,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny-1, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=v_north_grid,                &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = ly, ymax = ly,                                           &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=v_south_grid,                &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = 0.0_wp, ymax = 0.0_wp,                                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nz, z=z)

       CALL init_grid_definition('boundary', grid=v_top_grid,                  &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = dy, ymax = ly - dy,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny-1, nz = 1_iwp, z=z_top)

       CALL init_grid_definition('boundary', grid=w_east_grid,                 &
               xmin = lx + 0.5_wp * dx, xmax = lx + 0.5_wp * dx,               &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nz - 1, z=zw)

       CALL init_grid_definition('boundary', grid=w_west_grid,                 &
               xmin = -0.5_wp * dx, xmax = -0.5_wp * dx,                       &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nz - 1, z=zw)

       CALL init_grid_definition('boundary', grid=w_north_grid,                &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = ly + 0.5_wp * dy, ymax = ly + 0.5_wp * dy,               &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nz - 1, z=zw)

       CALL init_grid_definition('boundary', grid=w_south_grid,                &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin = -0.5_wp * dy, ymax = -0.5_wp * dy,                       &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nz - 1, z=zw)

       CALL init_grid_definition('boundary', grid=w_top_grid,                  &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny, nz = 1_iwp, z=zw_top)

       CALL init_grid_definition('boundary intermediate', grid=scalars_east_intermediate,   &
               xmin = lx + 0.5_wp * dx, xmax = lx + 0.5_wp * dx,               &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=scalars_west_intermediate,   &
               xmin = -0.5_wp * dx, xmax = -0.5_wp * dx,                       &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=scalars_north_intermediate,  &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = ly + 0.5_wp * dy, ymax = ly + 0.5_wp * dy,               &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=scalars_south_intermediate,  &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin = -0.5_wp * dy, ymax = -0.5_wp * dy,                       &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=scalars_top_intermediate,    &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=u_east_intermediate,         &
               xmin = lx, xmax = lx,                                           &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=u_west_intermediate,         &
               xmin = 0.0_wp, xmax = 0.0_wp,                                   &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=u_north_intermediate,        &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = ly + 0.5_wp * dy, ymax = ly + 0.5_wp * dy,               &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = 0_iwp, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=u_south_intermediate,        &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = -0.5_wp * dy, ymax = -0.5_wp * dy,                       &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = 0_iwp, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=u_top_intermediate,          &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = 0.5_wp * dy, ymax = ly - 0.5_wp * dy,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = ny, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=v_east_intermediate,         &
               xmin = lx + 0.5_wp * dx, xmax = lx + 0.5_wp * dx,               &
               ymin = dy, ymax = ly - dy,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny-1, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=v_west_intermediate,         &
               xmin = -0.5_wp * dx, xmax = -0.5_wp * dx,                       &
               ymin = dy, ymax = ly - dy,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny-1, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=v_north_intermediate,        &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = ly, ymax = ly,                                           &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=v_south_intermediate,        &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = 0.0_wp, ymax = 0.0_wp,                                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=v_top_intermediate,          &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = dy, ymax = ly - dy,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny-1, nz = nlev - 2)

       CALL init_grid_definition('boundary intermediate', grid=w_east_intermediate,         &
               xmin = lx + 0.5_wp * dx, xmax = lx + 0.5_wp * dx,               &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nlev - 1)

       CALL init_grid_definition('boundary intermediate', grid=w_west_intermediate,         &
               xmin = -0.5_wp * dx, xmax = -0.5_wp * dx,                       &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = 0_iwp, ny = ny, nz = nlev - 1)

       CALL init_grid_definition('boundary intermediate', grid=w_north_intermediate,        &
               xmin = 0.5_wp * dx, xmax = lx - 0.5_wp * dx,                    &
               ymin = ly + 0.5_wp * dy, ymax = ly + 0.5_wp * dy,               &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nlev - 1)

       CALL init_grid_definition('boundary intermediate', grid=w_south_intermediate,        &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin = -0.5_wp * dy, ymax = -0.5_wp * dy,                       &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = 0_iwp, nz = nlev - 1)

       CALL init_grid_definition('boundary intermediate', grid=w_top_intermediate,          &
               xmin =  0.5_wp * dx, xmax = lx - 0.5_wp * dx,                   &
               ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                   &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny, nz = nlev - 1)

    ENDIF

!                                                                              
!------------------------------------------------------------------------------
! Section 3: Define profile grids
!------------------------------------------------------------------------------

    lonmin_palm = MINVAL(palm_intermediate%clon)
    lonmax_palm = MAXVAL(palm_intermediate%clon)
    latmin_palm = MINVAL(palm_intermediate%clat)
    latmax_palm = MAXVAL(palm_intermediate%clat)

    lonmin_tot = MIN(lam_centre - averaging_angle, lonmin_palm)
    lonmax_tot = MAX(lam_centre + averaging_angle, lonmax_palm)
    latmin_tot = MIN(phi_centre - averaging_angle, latmin_palm)
    latmax_tot = MAX(phi_centre + averaging_angle, latmax_palm)

    palm_domain_outside_cosmo = ANY(                                           &
       (/ lonmin_tot,   -lonmax_tot,   latmin_tot,   -latmax_tot/) .LT.        &
       (/ lonmin_cosmo, -lonmax_cosmo, latmin_cosmo, -latmax_cosmo/)           &
    )

    IF ( palm_domain_outside_cosmo )  THEN
       message = 'PALM domain or geostrophic averaging domains extend ' //     &
                 'outside COSMO domain.'
       CALL inifor_abort( 'setup_grids', message )
    ENDIF

    CALL init_averaging_grid(averaged_soil_profile, cosmo_grid,                &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = depths, z0 = z0,             &
            lonmin = lonmin_palm, lonmax = lonmax_palm,                        &
            latmin = latmin_palm, latmax = latmax_palm,                        &
            kind='scalar', name='averaged soil profile')

    CALL init_averaging_grid(averaged_scalar_profile, cosmo_grid,              &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z, z0 = z0,                  &
            lonmin = lonmin_palm, lonmax = lonmax_palm,                        &
            latmin = latmin_palm, latmax = latmax_palm,                        &
            kind='scalar', name='averaged scalar profile')

    CALL init_averaging_grid(averaged_w_profile, cosmo_grid,                   &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = zw, z0 = z0,                 &
            lonmin = lonmin_palm, lonmax = lonmax_palm,                        &
            latmin = latmin_palm, latmax = latmax_palm,                        &
            kind='w', name='averaged w profile')

    CALL init_averaging_grid(averaged_scalar_top_point, cosmo_grid,            &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z_top, z0 = z0,              &
            lonmin = lonmin_palm, lonmax = lonmax_palm,                        &
            latmin = latmin_palm, latmax = latmax_palm,                        &
            kind='scalar', name='averaged scalar top point')

    CALL init_averaging_grid(averaged_w_top_point, cosmo_grid,                 &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = zw_top, z0 = z0,             &
            lonmin = lonmin_palm, lonmax = lonmax_palm,                        &
            latmin = latmin_palm, latmax = latmax_palm,                        &
            kind='w', name='averaged w top point')

    ! Also used for computing pressure profiles
    CALL init_averaging_grid(geostrophic_scalar_profile, cosmo_grid,           &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z, z0 = z0,                  &
            lonmin = lam_west, lonmax = lam_east,                              &
            latmin = phi_south, latmax = phi_north,                            &
            kind='scalar', name='centre geostrophic scalar profile')
!
!-- Initialize averaging grid for internal variables.
!-- Note, at the moment the if conditions is commented to avoid divisions by
!-- zero.
!     IF ( ls_forcing_variables_required )  THEN

    CALL init_averaging_grid(geostrophic_w_profile, cosmo_grid,                &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = zw, z0 = z0,                 &
            lonmin = lam_west, lonmax = lam_east,                              &
            latmin = phi_south, latmax = phi_north,                            &
            kind='w', name='centre geostrophic w profile')

    CALL init_averaging_grid(south_geostrophic_scalar_profile, cosmo_grid,     &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z, z0 = z0,                  &
            lonmin = lam_west, lonmax = lam_east,                              &
            latmin = phi_centre - averaging_angle,                             &
            latmax = phi_centre,                                               &
            kind='scalar', name='south geostrophic scalar profile')

    CALL init_averaging_grid(north_geostrophic_scalar_profile, cosmo_grid,     &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z, z0 = z0,                  &
            lonmin = lam_west, lonmax = lam_east,                              &
            latmin = phi_centre,                                               &
            latmax = phi_centre + averaging_angle,                             &
            kind='scalar', name='north geostrophic scalar profile')

    CALL init_averaging_grid(west_geostrophic_scalar_profile, cosmo_grid,      &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z, z0 = z0,                  &
            lonmin = lam_centre - averaging_angle,                             &
            lonmax = lam_centre,                                               &
            latmin = phi_south, latmax = phi_north,                            &
            kind='scalar', name='west geostrophic scalar profile')

    CALL init_averaging_grid(east_geostrophic_scalar_profile, cosmo_grid,      &
            x = 0.5_wp * lx, y = 0.5_wp * ly, z = z, z0 = z0,                  &
            lonmin = lam_centre,                                               &
            lonmax = lam_centre + averaging_angle,                             &
            latmin = phi_south, latmax = phi_north,                            &
            kind='scalar', name='east geostrophic scalar profile')

!     END IF


!                                                                              
!------------------------------------------------------------------------------
! Section 4: Precompute neighbours and weights for interpolation              
!------------------------------------------------------------------------------
    IF (cfg%map_terrain)  THEN
       CALL setup_terrain_mapping( palm_grid, zs = zs(:,:,:) )
       CALL setup_terrain_mapping( u_initial_grid, zs = zs_u(:,:,:) )
       CALL setup_terrain_mapping( v_initial_grid, zs = zs_v(:,:,:) )
       CALL setup_terrain_mapping( w_initial_grid, zs = zs(:,:,:) )
    ENDIF

    IF (cfg%map_terrain .AND. boundary_variables_required)  THEN
       CALL setup_terrain_mapping( scalars_east_grid, zs = zs(nx+1:nx+1,1:ny+1,:) )
       CALL setup_terrain_mapping( scalars_west_grid, zs = zs(1:1,1:ny+1,:) )
       CALL setup_terrain_mapping( scalars_north_grid, zs = zs(1:nx+1,ny+1:ny+1,:) )
       CALL setup_terrain_mapping( scalars_south_grid, zs = zs(1:nx+1,1:1,:) )
       CALL setup_terrain_mapping( scalars_top_grid, zs = zs(:,:,:) )
       CALL setup_terrain_mapping( u_east_grid, zs = zs(nx+1:nx+1,1:ny+1,:) )
       CALL setup_terrain_mapping( u_west_grid, zs = zs(1:1,1:ny+1,:) )
       CALL setup_terrain_mapping( u_north_grid, zs = zs(1:nx+1,ny+1:ny+1,:) )
       CALL setup_terrain_mapping( u_south_grid, zs = zs(1:nx+1,1:1,:) )
       CALL setup_terrain_mapping( u_top_grid, zs = zs(:,:,:) )
       CALL setup_terrain_mapping( v_east_grid, zs = zs(nx+1:nx+1,1:ny+1,:) )
       CALL setup_terrain_mapping( v_west_grid, zs = zs(1:1,1:ny+1,:) )
       CALL setup_terrain_mapping( v_north_grid,zs = zs(1:nx+1,ny+1:ny+1,:) )
       CALL setup_terrain_mapping( v_south_grid, zs = zs(1:nx+1,1:1,:) )
       CALL setup_terrain_mapping( v_top_grid, zs = zs(:,:,:) )
       CALL setup_terrain_mapping( w_east_grid, zs = zs(nx+1:nx+1,1:ny+1,:) )
       CALL setup_terrain_mapping( w_west_grid, zs = zs(1:1,1:ny+1,:) )
       CALL setup_terrain_mapping( w_north_grid, zs = zs(1:nx+1,ny+1:ny+1,:) )
       CALL setup_terrain_mapping( w_south_grid, zs = zs(1:nx+1,1:1,:) )
       CALL setup_terrain_mapping( w_top_grid, zs = zs(:,:,:) )
    ENDIF

    interp_mode = 's'
    CALL setup_interpolation(cosmo_grid, palm_grid, palm_intermediate, interp_mode, ic_mode=cfg%ic_mode)
    IF (boundary_variables_required)  THEN
       CALL setup_interpolation(cosmo_grid, scalars_east_grid, scalars_east_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, scalars_west_grid, scalars_west_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, scalars_north_grid, scalars_north_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, scalars_south_grid, scalars_south_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, scalars_top_grid, scalars_top_intermediate, interp_mode)
    ENDIF

    interp_mode = 'u'
    CALL setup_interpolation(cosmo_grid, u_initial_grid, u_initial_intermediate, interp_mode, ic_mode=cfg%ic_mode)
    IF (boundary_variables_required)  THEN
       CALL setup_interpolation(cosmo_grid, u_east_grid, u_east_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, u_west_grid, u_west_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, u_north_grid, u_north_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, u_south_grid, u_south_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, u_top_grid, u_top_intermediate, interp_mode)
    ENDIF

    interp_mode = 'v'
    CALL setup_interpolation(cosmo_grid, v_initial_grid, v_initial_intermediate, interp_mode, ic_mode=cfg%ic_mode)
    IF (boundary_variables_required)  THEN
       CALL setup_interpolation(cosmo_grid, v_east_grid, v_east_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, v_west_grid, v_west_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, v_north_grid, v_north_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, v_south_grid, v_south_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, v_top_grid, v_top_intermediate, interp_mode)
    ENDIF

    interp_mode = 'w'
    CALL setup_interpolation(cosmo_grid, w_initial_grid, w_initial_intermediate, interp_mode, ic_mode=cfg%ic_mode)
    IF (boundary_variables_required)  THEN
       CALL setup_interpolation(cosmo_grid, w_east_grid, w_east_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, w_west_grid, w_west_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, w_north_grid, w_north_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, w_south_grid, w_south_intermediate, interp_mode)
       CALL setup_interpolation(cosmo_grid, w_top_grid, w_top_intermediate, interp_mode)
    ENDIF

    IF (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)  THEN
        !TODO: remove this conditional if not needed. 
    ENDIF

 END SUBROUTINE setup_grids


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Driver for computing neighbour indices and weights for horizontal and
!> vertical interpolation.
!------------------------------------------------------------------------------!
 SUBROUTINE setup_interpolation(cosmo_grid, grid, intermediate_grid, kind, ic_mode)

    TYPE(grid_definition), INTENT(IN), TARGET    ::  cosmo_grid
    TYPE(grid_definition), INTENT(INOUT), TARGET ::  grid, intermediate_grid
    CHARACTER, INTENT(IN)                        ::  kind
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL       ::  ic_mode

    REAL(wp), DIMENSION(:), POINTER     ::  cosmo_lat, cosmo_lon
    REAL(wp), DIMENSION(:,:,:), POINTER ::  cosmo_h

    LOGICAL :: setup_volumetric

!------------------------------------------------------------------------------
! Section 1: Horizontal interpolation                                        
!------------------------------------------------------------------------------
!
!-- Select horizontal coordinates according to kind of points (s/w, u, v)
    SELECT CASE(kind)

!
!-- scalars
    CASE('s')

       cosmo_lat => cosmo_grid%lat
       cosmo_lon => cosmo_grid%lon
       cosmo_h   => cosmo_grid%hfl
!
!-- vertical velocity
    CASE('w')

       cosmo_lat => cosmo_grid%lat
       cosmo_lon => cosmo_grid%lon
       cosmo_h   => cosmo_grid%hhl
!
!-- x velocity
    CASE('u')

       cosmo_lat => cosmo_grid%lat
       cosmo_lon => cosmo_grid%lonu
       cosmo_h   => cosmo_grid%hfl

!
!-- y velocity
    CASE('v')

       cosmo_lat => cosmo_grid%latv
       cosmo_lon => cosmo_grid%lon
       cosmo_h   => cosmo_grid%hfl

    CASE DEFAULT

       message = "Interpolation quantity '" // kind // "' is not supported."
       CALL inifor_abort('setup_interpolation', message)

    END SELECT

    CALL find_horizontal_neighbours(cosmo_lat, cosmo_lon,                      &
       intermediate_grid%clat, intermediate_grid%clon,                         &
       intermediate_grid%ii, intermediate_grid%jj)

    CALL compute_horizontal_interp_weights(cosmo_lat, cosmo_lon,               &
       intermediate_grid%clat, intermediate_grid%clon,                         &
       intermediate_grid%ii, intermediate_grid%jj,                             &
       intermediate_grid%w_horiz)

!------------------------------------------------------------------------------
! Section 2: Vertical interpolation
!------------------------------------------------------------------------------

!
!-- If profile initialization is chosen, we--somewhat counterintuitively--
!-- don't need to compute vertical interpolation weights. At least, we
!-- don't need them on the intermediate grid, which fills the entire PALM
!-- domain volume. Instead we need vertical weights for the intermediate
!-- profile grids, which get computed in setup_averaging().
    setup_volumetric = .TRUE.
    IF (PRESENT(ic_mode))  THEN
       IF (TRIM(ic_mode) == CFG_INIT_PROFILE)  setup_volumetric = .FALSE.
    ENDIF

    IF (setup_volumetric)  THEN
       ALLOCATE( intermediate_grid%intermediate_h(0:intermediate_grid%nx,      &
                                                  0:intermediate_grid%ny,      &
                                                  0:intermediate_grid%nz) ) 
       intermediate_grid%intermediate_h(:,:,:) = - EARTH_RADIUS

!
!--    For w points, use hhl, for scalars use hfl
!--    compute the full heights for the intermediate grids
       CALL interpolate_2d(cosmo_h, intermediate_grid%intermediate_h, intermediate_grid)

       IF (cfg%map_terrain)  THEN

!
!--       Interpolate mesoscale terrain height (zs) on  every intermediate grid
!--       column.
          CALL interpolate_2d(                                                 &
             invar = cosmo_grid%hhl(:,:,1:1),                                  &
             outvar = intermediate_grid%zs,                                    &
             outgrid = intermediate_grid                                       &
          )

!
!--       Save mesoscale model top height (zt).
          intermediate_grid%zt = cosmo_grid%hhl(1,1,UBOUND( cosmo_grid%hhl, 3))

          CALL compute_sigma(                                                  &
             intermediate_grid%intermediate_h(:,:,:),                          &
             intermediate_grid%zs(:,:,:),                                      &
             intermediate_grid%zt,                                             &
             intermediate_grid%sigma(:,:,:)                                    &
          )

          CALL map_terrain_driver( intermediate_grid, grid%zs_palm,            &
                                   blending_dz, blending_z_ubound )

       ENDIF

       CALL find_vertical_neighbours_and_weights_interp(grid, intermediate_grid)
    ENDIF
       
 END SUBROUTINE setup_interpolation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes grid_definition-type variables.
!> 
!> Input parameters:
!> -----------------
!> kind : Grid kind, distinguishes between PALM-4U and COSMO-DE grids
!>    as well as grids covering the boundary surfaces. Valid kinds are:
!>       - 'palm'
!>       - 'cosmo-de'
!>       - 'eastwest-scalar'
!> 
!> <xyx>min, <xyz>max : Domain minima and maxima in x, y, and z direction. Note
!>    that these values do not necessarily translate to the outmost coordinates
!>    of the generated grid but rather refer to the extent of the underlying
!>    PALM-4U computational domain (i.e. the outer cell faces). The coordinates
!>    of the generated grid will be inferred from this information taking into
!>    account the initialization mode ic_mode. For example, the coordinates of a
!>    boundary grid initialized using mode 'eastwest-scalar' will be located in
!>    planes one half grid point outwards of xmin and xmax.
!>
!> z0 : Elevation of the PALM-4U domain above sea level [m]
!>
!> n<xyz> : Number of grod points in x, y, and z direction
!>
!> Output parameters:
!> ------------------
!> grid : Grid variable to be initialized.
!------------------------------------------------------------------------------!
 SUBROUTINE init_grid_definition( kind, xmin, xmax, ymin, ymax,                &
                                  x0, y0, z0, nx, ny, nz, z, zw, grid, ic_mode )
    CHARACTER(LEN=*), INTENT(IN)           ::  kind
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL ::  ic_mode
    INTEGER(iwp), INTENT(IN)               ::  nx, ny, nz
    REAL(wp), INTENT(IN)                   ::  xmin, xmax, ymin, ymax
    REAL(wp), INTENT(IN)                   ::  x0, y0, z0
    REAL(wp), INTENT(IN), TARGET, OPTIONAL ::  z(:)
    REAL(wp), INTENT(IN), TARGET, OPTIONAL ::  zw(:)
    TYPE(grid_definition), INTENT(INOUT)   ::  grid

    grid%nx = nx
    grid%ny = ny
    grid%nz = nz

    grid%lx = xmax - xmin
    grid%ly = ymax - ymin

    grid%x0 = x0
    grid%y0 = y0
    grid%z0 = z0

    SELECT CASE( TRIM(kind) )

       CASE('boundary')
       
          IF (.NOT.PRESENT(z))  THEN
             message = "z has not been passed but is required for 'boundary' grids"
             CALL inifor_abort('init_grid_definition', message)
          ENDIF
       
          ALLOCATE( grid%x(0:nx) )
          CALL linspace(xmin, xmax, grid%x)
       
          ALLOCATE( grid%y(0:ny) )
          CALL linspace(ymin, ymax, grid%y)
       
          grid%z => z
!
!--       Allocate neighbour indices and weights
          IF (TRIM(ic_mode) .NE. CFG_INIT_PROFILE)  THEN
             ALLOCATE( grid%kk(0:nx, 0:ny, 1:nz, 2) )
             grid%kk(:,:,:,:) = -1
       
             ALLOCATE( grid%w_verti(0:nx, 0:ny, 1:nz, 2) )
             grid%w_verti(:,:,:,:) = 0.0_wp
          ENDIF
       
       CASE('boundary intermediate')
       
          ALLOCATE( grid%x(0:nx) )
          CALL linspace(xmin, xmax, grid%x)
       
          ALLOCATE( grid%y(0:ny) )
          CALL linspace(ymin, ymax, grid%y)
       
          ALLOCATE( grid%clon(0:nx, 0:ny), grid%clat(0:nx, 0:ny)  )
       
          CALL rotate_to_cosmo(                                                &
             phir = project( grid%y, y0, EARTH_RADIUS ) ,                      & ! = plate-carree latitude
             lamr = project( grid%x, x0, EARTH_RADIUS ) ,                      & ! = plate-carree longitude
             phip = phi_cn, lamp = lambda_cn,                                  &
             phi  = grid%clat,                                                 &
             lam  = grid%clon,                                                 &
             gam  = gam                                                        &
          )

!
!--       For intermediate grids, set mesoscale terrain height. Initializing
!--       with zero; will be filled with horizontally interpolated mesoscale
!--       terrain in setup_interpolation().
          IF (cfg%map_terrain)  THEN
             ALLOCATE( grid%zs(0:nx,0:ny,1) )
             ALLOCATE( grid%sigma(0:nx,0:ny,0:nz) )
             grid%zs(:,:,1) = 0.0_wp
             grid%sigma(:,:,:) = 0.0_wp
          ENDIF
       
!      
!--       Allocate neighbour indices and weights
          ALLOCATE( grid%ii(0:nx, 0:ny, 4),                                    &
                    grid%jj(0:nx, 0:ny, 4) )
          grid%ii(:,:,:)   = -1
          grid%jj(:,:,:)   = -1
       
          ALLOCATE( grid%w_horiz(0:nx, 0:ny, 4) )
          grid%w_horiz(:,:,:)   = 0.0_wp
       
!      
!--    This mode initializes a Cartesian PALM-4U grid and adds the
!--    corresponding latitudes and longitudes of the rotated pole grid.
       CASE('palm')
       
          IF (.NOT.PRESENT(z))  THEN
             message = "z has not been passed but is required for 'palm' grids"
             CALL inifor_abort('init_grid_definition', message)
          ENDIF
       
          IF (.NOT.PRESENT(zw))  THEN
             message = "zw has not been passed but is required for 'palm' grids"
             CALL inifor_abort('init_grid_definition', message)
          ENDIF
       
          grid%name(1) = 'x and lon'
          grid%name(2) = 'y and lat'
          grid%name(3) = 'z'
       
!      
!--       TODO: Remove use of global dx, dy, dz variables. Consider
!--       TODO: associating global x,y, and z arrays.
          ALLOCATE( grid%x(0:nx),   grid%y(0:ny) )
          ALLOCATE( grid%xu(1:nx),  grid%yv(1:ny) )
          CALL linspace(xmin + 0.5_wp* dx, xmax - 0.5_wp* dx, grid%x)
          CALL linspace(ymin + 0.5_wp* dy, ymax - 0.5_wp* dy, grid%y)
          grid%z => z
          CALL linspace(xmin +  dx, xmax -  dx, grid%xu)
          CALL linspace(ymin +  dy, ymax -  dy, grid%yv)
          grid%zw => zw
       
          grid%depths => depths
       
!      
!--       Allocate neighbour indices and weights
          IF (TRIM(ic_mode) .NE. CFG_INIT_PROFILE)  THEN
             ALLOCATE( grid%kk(0:nx, 0:ny, 1:nz, 2) )
             grid%kk(:,:,:,:) = -1
       
             ALLOCATE( grid%w_verti(0:nx, 0:ny, 1:nz, 2) )
             grid%w_verti(:,:,:,:) = 0.0_wp
          ENDIF
       
       CASE('palm intermediate')
       
          grid%name(1) = 'x and lon'
          grid%name(2) = 'y and lat'
          grid%name(3) = 'interpolated hhl or hfl'
       
!      
!--       TODO: Remove use of global dx, dy, dz variables. Consider
!--       TODO: associating global x,y, and z arrays.
          ALLOCATE( grid%x(0:nx),   grid%y(0:ny) )
          ALLOCATE( grid%xu(1:nx),  grid%yv(1:ny) )
          CALL linspace(xmin + 0.5_wp*dx, xmax - 0.5_wp*dx, grid%x)
          CALL linspace(ymin + 0.5_wp*dy, ymax - 0.5_wp*dy, grid%y)
          CALL linspace(xmin + dx, xmax - dx, grid%xu)
          CALL linspace(ymin + dy, ymax - dy, grid%yv)
       
          grid%depths => depths
       
!      
!--       Allocate rotated-pole coordinates, clon is for (c)osmo-de (lon)gitude
          ALLOCATE( grid%clon(0:nx, 0:ny),   grid%clat(0:nx, 0:ny)  )
          ALLOCATE( grid%clonu(1:nx, 0:ny),  grid%clatu(1:nx, 0:ny) )
          ALLOCATE( grid%clonv(0:nx, 1:ny),  grid%clatv(0:nx, 1:ny) )

!      
!--       Compute rotated-pole coordinates of...
!--       ... PALM-4U centres
          CALL rotate_to_cosmo(                                                &
             phir = project( grid%y, y0, EARTH_RADIUS ) , & ! = plate-carree latitude
             lamr = project( grid%x, x0, EARTH_RADIUS ) , & ! = plate-carree longitude
             phip = phi_cn, lamp = lambda_cn,                                  &
             phi  = grid%clat,                                                 &
             lam  = grid%clon,                                                 &
             gam  = gam                                                        &
          )
       
!      
!--       ... PALM-4U u winds
          CALL rotate_to_cosmo(                                                &
             phir = project( grid%y,  y0, EARTH_RADIUS ), & ! = plate-carree latitude
             lamr = project( grid%xu, x0, EARTH_RADIUS ), & ! = plate-carree longitude
             phip = phi_cn, lamp = lambda_cn,                                  &
             phi  = grid%clatu,                                                &
             lam  = grid%clonu,                                                &
             gam  = gam                                                        &
          )
       
!      
!--       ... PALM-4U v winds
          CALL rotate_to_cosmo(                                                &
             phir = project( grid%yv, y0, EARTH_RADIUS ), & ! = plate-carree latitude
             lamr = project( grid%x,  x0, EARTH_RADIUS ), & ! = plate-carree longitude
             phip = phi_cn, lamp = lambda_cn,                                  &
             phi  = grid%clatv,                                                &
             lam  = grid%clonv,                                                &
             gam  = gam                                                        &
          )

!
!--       For intermediate grids, set mesoscale terrain height. Initializing
!--       with zero; will be filled with horizontally interpolated mesoscale
!--       terrain in setup_interpolation().
          IF (cfg%map_terrain)  THEN
             ALLOCATE( grid%zs(0:nx,0:ny,1) )
             ALLOCATE( grid%sigma(0:nx,0:ny,0:nz) )
             grid%zs(:,:,1) = 0.0_wp
             grid%sigma(:,:,:) = 0.0_wp
          ENDIF

!      
!--       Allocate neighbour indices and weights
          ALLOCATE( grid%ii(0:nx, 0:ny, 4),                                    &
                    grid%jj(0:nx, 0:ny, 4) )
          grid%ii(:,:,:)   = -1
          grid%jj(:,:,:)   = -1
       
          ALLOCATE( grid%w_horiz(0:nx, 0:ny, 4) )
          grid%w_horiz(:,:,:)   = 0.0_wp
       
       CASE('cosmo-de')
          grid%name(1) = 'rlon'         ! of COMSO-DE cell centres (scalars)
          grid%name(2) = 'rlat'         ! of COMSO-DE cell centres (scalars)
          grid%name(3) = 'height'
       
          ALLOCATE( grid%lon(0:nx),   grid%lat(0:ny)  )
          ALLOCATE( grid%lonu(0:nx),  grid%latv(0:ny) )
       
          CALL linspace(xmin, xmax, grid%lon)
          CALL linspace(ymin, ymax, grid%lat)
          grid%lonu(:) = grid%lon + 0.5_wp * (grid%lx / grid%nx)
          grid%latv(:) = grid%lat + 0.5_wp * (grid%ly / grid%ny)
       
!      
!--       Point to heights of half levels (hhl) and compute heights of full
!--       levels (hfl) as arithmetic averages
          grid%hhl => hhl
          grid%hfl => hfl
          grid%depths => depths
       
       CASE DEFAULT
           message = "Grid kind '" // TRIM(kind) // "' is not recognized."
           CALL inifor_abort('init_grid_definition', message)

    END SELECT

 END SUBROUTINE init_grid_definition


!------------------------------------------------------------------------------!
! Description:
! ------------
!> ....
!------------------------------------------------------------------------------!
 SUBROUTINE setup_terrain_mapping( boundary_grid, zs)

    TYPE(grid_definition), INTENT(INOUT) :: boundary_grid
    REAL(wp), TARGET, INTENT(IN)         :: zs(0:,0:,1:)

    boundary_grid%zs_palm(0:,0:,1:) => zs(:,:,:)

 END SUBROUTINE setup_terrain_mapping

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes averaging grids
!> 
!> Input parameters:
!> -----------------
!>
!> cosmo_grid : grid_definition-type varieble initialized as 'cosmo-de' grid
!>    providing longitudes and latitudes of the parent grid
!>
!> x, y : location of the profile to be averaged in the PALM system [m]
!>
!> z : vertical grid coordinates of the profile in the PALM system [m]
!>
!> <lat/lon><min/max>: boundaries of the averaging region in the parent system,
!>    i.e. the COSMO-DE rotated-pole system. [rad]
!>
!> kind : kind of quantity to be averaged using this averaging grid.
!>    Destinguishes COSMO-DE scalar and w-velocity levels. Note that finding the
!>    parent/COSMO columns for the region in get_latlon_averaging_region() is
!>    independent of 'kind' b/ecause INIFOR uses column-centred u and v velocity
!>    components, which are computed in the preprocessing step.
!>
!> Output parameters:
!> ------------------
!> avg_grid : averagin grid to be initialized
!------------------------------------------------------------------------------!
 SUBROUTINE init_averaging_grid(avg_grid, cosmo_grid, x, y, z, z0,             &
    lonmin, lonmax, latmin, latmax, kind, name)

    TYPE(grid_definition), INTENT(INOUT) ::  avg_grid
    TYPE(grid_definition), INTENT(IN)    ::  cosmo_grid
    REAL(wp), INTENT(IN)                 ::  x, y, z0
    REAL(wp), INTENT(IN), TARGET         ::  z(:)
    REAL(wp), INTENT(IN)                 ::  lonmin !< lower longitude bound of the averaging grid region [COSMO rotated-pole rad]
    REAL(wp), INTENT(IN)                 ::  lonmax !< upper longitude bound of the averaging grid region [COSMO rotated-pole rad]
    REAL(wp), INTENT(IN)                 ::  latmin !< lower latitude bound of the averaging grid region [COSMO rotated-pole rad]
    REAL(wp), INTENT(IN)                 ::  latmax !< lower latitude bound of the averaging grid region [COSMO rotated-pole rad]

    CHARACTER(LEN=*), INTENT(IN)         ::  kind
    CHARACTER(LEN=*), INTENT(IN)         ::  name

    LOGICAL                              ::  level_based_averaging

    ALLOCATE( avg_grid%x(1) )
    ALLOCATE( avg_grid%y(1) )
    avg_grid%x(1) = x
    avg_grid%y(1) = y
    avg_grid%z => z
    avg_grid%z0 = z0

    avg_grid%nz = SIZE(z, 1)

    ALLOCATE( avg_grid%lon(2) )
    ALLOCATE( avg_grid%lat(2) )
    avg_grid%lon(1:2) = (/lonmin, lonmax/)
    avg_grid%lat(1:2) = (/latmin, latmax/)

    avg_grid%kind = TRIM(kind)
    avg_grid%name(1) = TRIM(name)

!
!-- Find and store COSMO columns that fall into the coordinate range
!-- given by avg_grid%clon, %clat
    CALL get_latlon_averaging_region(avg_grid, cosmo_grid)

    ALLOCATE (avg_grid%kkk(avg_grid%n_columns, avg_grid%nz, 2) )
    ALLOCATE (avg_grid%w(avg_grid%n_columns, avg_grid%nz, 2) )
!
!-- Compute average COSMO levels in the averaging region
    SELECT CASE(avg_grid%kind)

       CASE('scalar', 'u', 'v')
          avg_grid%cosmo_h => cosmo_grid%hfl

       CASE('w')
          avg_grid%cosmo_h => cosmo_grid%hhl

       CASE DEFAULT
          message = "Averaging grid kind '" // TRIM(avg_grid%kind) // &
                    "' is not supported. Use 'scalar', 'u', or 'v'."
          CALL inifor_abort('init_averaging_grid', message)

    END SELECT

!
!-- For level-based averaging, compute average heights
    level_based_averaging = ( TRIM(cfg%averaging_mode) == 'level' )
    IF (level_based_averaging)  THEN
       ALLOCATE(avg_grid%intermediate_h(1,1,SIZE(avg_grid%cosmo_h, 3)) )
 
       CALL average_2d(avg_grid%cosmo_h, avg_grid%intermediate_h(1,1,:),       &
                       avg_grid%iii, avg_grid%jjj)

    ENDIF

!
!-- Compute vertical weights and neighbours
    CALL find_vertical_neighbours_and_weights_average(                         &
       avg_grid, level_based_averaging                                         &
    )

 END SUBROUTINE init_averaging_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_latlon_averaging_region() finds all mesocsale columns within the
!> latitude-longitude reactagle given by the four values in avg_grid%lon(1:2)
!> and %lat(1:2). The total number of all found columns is stored in
!> avg_grid%n_columns, and their indices are stored in the sequential lists
!> avg_grid%iii(:) and %jjj(:).
!------------------------------------------------------------------------------!
 SUBROUTINE get_latlon_averaging_region(avg_grid, cosmo_grid)
    TYPE(grid_definition), INTENT(INOUT)         ::  avg_grid
    TYPE(grid_definition), TARGET, INTENT(IN)    ::  cosmo_grid

    REAL(wp), DIMENSION(:), POINTER              ::  cosmo_lon, cosmo_lat
    REAL(wp)                                     ::  dlon, dlat

    INTEGER(iwp) ::  i, j, imin, imax, jmin, jmax, l, nx, ny


    SELECT CASE( TRIM(avg_grid%kind) )

       CASE('scalar', 'w')
          cosmo_lon => cosmo_grid%lon
          cosmo_lat => cosmo_grid%lat

       CASE('u')
          cosmo_lon => cosmo_grid%lonu
          cosmo_lat => cosmo_grid%lat

       CASE('v')
          cosmo_lon => cosmo_grid%lon
          cosmo_lat => cosmo_grid%latv

       CASE DEFAULT
          message = "Averaging grid kind '" // TRIM(avg_grid%kind) // &
                    "' is not supported. Use 'scalar', 'u', or 'v'."
          CALL inifor_abort('get_latlon_averaging_region', message)

    END SELECT

    dlon = cosmo_lon(1) - cosmo_lon(0)
    dlat = cosmo_lat(1) - cosmo_lat(0)

    imin = FLOOR  ( (avg_grid%lon(1) - cosmo_lon(0)) / dlon )
    imax = CEILING( (avg_grid%lon(2) - cosmo_lon(0)) / dlon )

    jmin = FLOOR  ( (avg_grid%lat(1) - cosmo_lat(0)) / dlat )
    jmax = CEILING( (avg_grid%lat(2) - cosmo_lat(0)) / dlat )
    
    message = "Grid " // TRIM(avg_grid%name(1)) // " averages over " //        &
              TRIM(str(imin)) // " <= i <= " // TRIM(str(imax)) //             &
              " and " //                                                       &
              TRIM(str(jmin)) // " <= j <= " // TRIM(str(jmax))
    CALL report( 'get_latlon_averaging_region', message )

    nx = imax - imin + 1
    ny = jmax - jmin + 1
    avg_grid%n_columns = nx * ny

    ALLOCATE( avg_grid%iii(avg_grid%n_columns),                                &
              avg_grid%jjj(avg_grid%n_columns) )

    l = 0
    DO  j = jmin, jmax
    DO  i = imin, imax
       l = l + 1
       avg_grid%iii(l) = i
       avg_grid%jjj(l) = j
    ENDDO
    ENDDO

 END SUBROUTINE get_latlon_averaging_region


!------------------------------------------------------------------------------!
! Description:
! ------------
!> PALM's stretched vertical grid generator. Forked from PALM revision 3139, see
!> https://palm.muk.uni-hannover.de/trac/browser/palm/trunk/SOURCE/init_grid.f90?rev=3139
!> 
!> This routine computes the levels of scalar points. The levels of the velocity
!> points are then obtained as the midpoints inbetween using the INIFOR routine
!> 'modpoints'.
!------------------------------------------------------------------------------!
 SUBROUTINE stretched_z(z, dz, dz_max, dz_stretch_factor, dz_stretch_level,    &
                        dz_stretch_level_start, dz_stretch_level_end,          &
                        dz_stretch_factor_array)

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  z, dz, dz_stretch_factor_array
    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  dz_stretch_level_start, dz_stretch_level_end
    REAL(wp), INTENT(IN) ::  dz_max, dz_stretch_factor, dz_stretch_level

    INTEGER(iwp) ::  number_stretch_level_start        !< number of user-specified start levels for stretching
    INTEGER(iwp) ::  number_stretch_level_end          !< number of user-specified end levels for stretching

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  min_dz_stretch_level_end
    REAL(wp) ::  dz_level_end, dz_stretched

    INTEGER(iwp) ::  dz_stretch_level_end_index(9)      !< vertical grid level index until which the vertical grid spacing is stretched
    INTEGER(iwp) ::  dz_stretch_level_start_index(9)    !< vertical grid level index above which the vertical grid spacing is stretched
    INTEGER(iwp) ::  dz_stretch_level_index = 0
    INTEGER(iwp) ::  k, n, number_dz

!
!-- Compute height of u-levels from constant grid length and dz stretch factors
    IF ( dz(1) == -1.0_wp )  THEN
       message = 'missing dz'
       CALL inifor_abort( 'stretched_z', message)
    ELSEIF ( dz(1) <= 0.0_wp )  THEN
       WRITE( message, * ) 'dz=', dz(1),' <= 0.0'
       CALL inifor_abort( 'stretched_z', message)
    ENDIF

!
!-- Initialize dz_stretch_level_start with the value of dz_stretch_level
!-- if it was set by the user
    IF ( dz_stretch_level /= -9999999.9_wp )  THEN
       dz_stretch_level_start(1) = dz_stretch_level
    ENDIF
       
!
!-- Determine number of dz values and stretching levels specified by the
!-- user to allow right controlling of the stretching mechanism and to
!-- perform error checks. The additional requirement that dz /= dz_max
!-- for counting number of user-specified dz values is necessary. Otherwise
!-- restarts would abort if the old stretching mechanism with dz_stretch_level
!-- is used (Attention: The user is not allowed to specify a dz value equal
!-- to the default of dz_max = 999.0). 
    number_dz = COUNT( dz /= -1.0_wp .AND. dz /= dz_max )
    number_stretch_level_start = COUNT( dz_stretch_level_start /=              &
                                        -9999999.9_wp )
    number_stretch_level_end = COUNT( dz_stretch_level_end /=                  &
                                      9999999.9_wp )

!
!-- The number of specified end levels +1 has to be the same than the number 
!-- of specified dz values
    IF ( number_dz /= number_stretch_level_end + 1 )  THEN
       WRITE( message, * ) 'The number of values for dz = ',                   &
                           number_dz, 'has to be the same than ',              &
                           'the number of values for ',                        &
                           'dz_stretch_level_end + 1 = ',                      &
                           number_stretch_level_end+1
       CALL inifor_abort( 'stretched_z', message)
    ENDIF
    
!
!-- The number of specified start levels has to be the same or one less than 
!-- the number of specified dz values
    IF ( number_dz /= number_stretch_level_start + 1 .AND.                     &
         number_dz /= number_stretch_level_start )  THEN
       WRITE( message, * ) 'The number of values for dz = ',                   &
                           number_dz, 'has to be the same or one ',            &
                           'more than& the number of values for ',             &
                           'dz_stretch_level_start = ',                        &
                           number_stretch_level_start
       CALL inifor_abort( 'stretched_z', message)
    ENDIF
    
!-- The number of specified start levels has to be the same or one more than 
!-- the number of specified end levels
    IF ( number_stretch_level_start /= number_stretch_level_end + 1 .AND.      &
         number_stretch_level_start /= number_stretch_level_end )  THEN
       WRITE( message, * ) 'The number of values for ',                        &
                           'dz_stretch_level_start = ',                        &
                           dz_stretch_level_start, 'has to be the ',           &
                           'same or one more than& the number of ',            &
                           'values for dz_stretch_level_end = ',               &
                           number_stretch_level_end
       CALL inifor_abort( 'stretched_z', message)
    ENDIF

!
!-- Initialize dz for the free atmosphere with the value of dz_max
    IF ( dz(number_stretch_level_start+1) == -1.0_wp .AND.                     &
         number_stretch_level_start /= 0 )  THEN 
       dz(number_stretch_level_start+1) = dz_max
    ENDIF
       
!
!-- Initialize the stretching factor if (infinitely) stretching in the free 
!-- atmosphere is desired (dz_stretch_level_end was not specified for the 
!-- free atmosphere)
    IF ( number_stretch_level_start == number_stretch_level_end + 1 )  THEN 
       dz_stretch_factor_array(number_stretch_level_start) =                   &
       dz_stretch_factor
    ENDIF

!-- Allocation of arrays for stretching
    ALLOCATE( min_dz_stretch_level_end(number_stretch_level_start) )

!
!-- The stretching region has to be large enough to allow for a smooth
!-- transition between two different grid spacings
    DO  n = 1, number_stretch_level_start
       min_dz_stretch_level_end(n) = dz_stretch_level_start(n) +               &
                                     4 * MAX( dz(n),dz(n+1) )
    ENDDO

    IF ( ANY( min_dz_stretch_level_end(1:number_stretch_level_start) >         &
              dz_stretch_level_end(1:number_stretch_level_start) ) )  THEN
          message = 'Each dz_stretch_level_end has to be larger '  //          &
                    'than its corresponding value for ' //                     &
                    'dz_stretch_level_start + 4*MAX(dz(n),dz(n+1)) '//         &
                    'to allow for smooth grid stretching'
          CALL inifor_abort('stretched_z', message)
    ENDIF
    
!
!-- Stretching must not be applied within the prandtl_layer 
!-- (first two grid points). For the default case dz_stretch_level_start 
!-- is negative. Therefore the absolut value is checked here.
    IF ( ANY( ABS( dz_stretch_level_start ) < dz(1) * 1.5_wp ) )  THEN
       WRITE( message, * ) 'Eeach dz_stretch_level_start has to be ',          &
                           'larger than ', dz(1) * 1.5
          CALL inifor_abort( 'stretched_z', message)
    ENDIF

!
!-- The stretching has to start and end on a grid level. Therefore 
!-- user-specified values have to ''interpolate'' to the next lowest level
    IF ( number_stretch_level_start /= 0 )  THEN
       dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) -           &
                                         dz(1)/2.0) / dz(1) )                  &
                                   * dz(1) + dz(1)/2.0
    ENDIF
    
    IF ( number_stretch_level_start > 1 )  THEN
       DO  n = 2, number_stretch_level_start
          dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) /         &
                                           dz(n) ) * dz(n)
       ENDDO
    ENDIF
    
    IF ( number_stretch_level_end /= 0 )  THEN
       DO  n = 1, number_stretch_level_end
          dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) /             &
                                         dz(n+1) ) * dz(n+1)
       ENDDO
    ENDIF
 
!
!-- Determine stretching factor if necessary
    IF ( number_stretch_level_end >= 1 )  THEN 
       CALL calculate_stretching_factor( number_stretch_level_end, dz,         & 
                                         dz_stretch_factor_array,              &   
                                         dz_stretch_level_end,                 &
                                         dz_stretch_level_start )
    ENDIF

    z(1) = dz(1) * 0.5_wp
!
    dz_stretch_level_index = n
    dz_stretched = dz(1)
    DO  k = 2, n

       IF ( dz_stretch_level <= z(k-1)  .AND.  dz_stretched < dz_max )  THEN

          dz_stretched = dz_stretched * dz_stretch_factor
          dz_stretched = MIN( dz_stretched, dz_max )

          IF ( dz_stretch_level_index == n )  dz_stretch_level_index = k-1

       ENDIF

       z(k) = z(k-1) + dz_stretched

    ENDDO
!-- Determine u and v height levels considering the possibility of grid
!-- stretching in several heights.
    n = 1
    dz_stretch_level_start_index(:) = UBOUND(z, 1)
    dz_stretch_level_end_index(:) = UBOUND(z, 1)
    dz_stretched = dz(1)

!-- The default value of dz_stretch_level_start is negative, thus the first
!-- condition is always true. Hence, the second condition is necessary.
    DO  k = 2, UBOUND(z, 1)
       IF ( dz_stretch_level_start(n) <= z(k-1) .AND.                          &
            dz_stretch_level_start(n) /= -9999999.9_wp )  THEN
          dz_stretched = dz_stretched * dz_stretch_factor_array(n)
          
          IF ( dz(n) > dz(n+1) )  THEN
             dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
          ELSE
             dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
          ENDIF
          
          IF ( dz_stretch_level_start_index(n) == UBOUND(z, 1) )               &
          dz_stretch_level_start_index(n) = k-1
          
       ENDIF
       
       z(k) = z(k-1) + dz_stretched
       
!
!--    Make sure that the stretching ends exactly at dz_stretch_level_end 
       dz_level_end = ABS( z(k) - dz_stretch_level_end(n) ) 
       
       IF ( dz_level_end < dz(n+1)/3.0 )  THEN
          z(k) = dz_stretch_level_end(n)
          dz_stretched = dz(n+1)
          dz_stretch_level_end_index(n) = k
          n = n + 1             
       ENDIF
    ENDDO

    DEALLOCATE( min_dz_stretch_level_end )

 END SUBROUTINE stretched_z


!------------------------------------------------------------------------------!
! Description: [PALM subroutine]
! -----------------------------------------------------------------------------!
!> Calculation of the stretching factor through an iterative method. Ideas were 
!> taken from the paper "Regional stretched grid generation and its application
!> to the NCAR RegCM (1999)". Normally, no analytic solution exists because the
!> system of equations has two variables (r,l) but four requirements 
!> (l=integer, r=[0,88;1,2], Eq(6), Eq(5) starting from index j=1) which
!> results into an overdetermined system. 
!------------------------------------------------------------------------------!
 SUBROUTINE calculate_stretching_factor( number_end, dz,                       &
                                         dz_stretch_factor_array,              &   
                                         dz_stretch_level_end,                 &
                                         dz_stretch_level_start )
 
    REAL(wp), DIMENSION(:), INTENT(IN)    ::  dz
    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  dz_stretch_factor_array
    REAL(wp), DIMENSION(:), INTENT(IN)    ::  dz_stretch_level_end, dz_stretch_level_start
  
    INTEGER(iwp) ::  iterations  !< number of iterations until stretch_factor_lower/upper_limit is reached  
    INTEGER(iwp) ::  l_rounded   !< after l_rounded grid levels dz(n) is strechted to dz(n+1) with stretch_factor_2 
    INTEGER(iwp) ::  n           !< loop variable for stretching
    
    INTEGER(iwp), INTENT(IN) ::  number_end !< number of user-specified end levels for stretching
        
    REAL(wp) ::  delta_l               !< absolute difference between l and l_rounded
    REAL(wp) ::  delta_stretch_factor  !< absolute difference between stretch_factor_1 and stretch_factor_2
    REAL(wp) ::  delta_total_new       !< sum of delta_l and delta_stretch_factor for the next iteration (should be as small as possible) 
    REAL(wp) ::  delta_total_old       !< sum of delta_l and delta_stretch_factor for the last iteration 
    REAL(wp) ::  distance              !< distance between dz_stretch_level_start and dz_stretch_level_end (stretching region)
    REAL(wp) ::  l                     !< value that fulfil Eq. (5) in the paper mentioned above together with stretch_factor_1 exactly
    REAL(wp) ::  numerator             !< numerator of the quotient
    REAL(wp) ::  stretch_factor_1      !< stretching factor that fulfil Eq. (5) togehter with l exactly
    REAL(wp) ::  stretch_factor_2      !< stretching factor that fulfil Eq. (6) togehter with l_rounded exactly
    
    REAL(wp) ::  dz_stretch_factor_array_2(9) = 1.08_wp  !< Array that contains all stretch_factor_2 that belongs to stretch_factor_1 
    
    REAL(wp), PARAMETER ::  stretch_factor_interval = 1.0E-06  !< interval for sampling possible stretching factors
    REAL(wp), PARAMETER ::  stretch_factor_lower_limit = 0.88  !< lowest possible stretching factor
    REAL(wp), PARAMETER ::  stretch_factor_upper_limit = 1.12  !< highest possible stretching factor
 
 
    l = 0
    DO  n = 1, number_end
    
       iterations = 1
       stretch_factor_1 = 1.0 
       stretch_factor_2 = 1.0
       delta_total_old = 1.0
       
       IF ( dz(n) > dz(n+1) )  THEN
          DO  WHILE ( stretch_factor_1 >= stretch_factor_lower_limit ) 
             
             stretch_factor_1 = 1.0 - iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                         &
                        dz_stretch_level_start(n) ) 
             numerator = distance*stretch_factor_1/dz(n) +                     &
                         stretch_factor_1 - distance/dz(n)
             
             IF ( numerator > 0.0 )  THEN
                l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0
                l_rounded = NINT( l )
                delta_l = ABS( l_rounded - l ) / l
             ENDIF
             
             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )
             
             delta_stretch_factor = ABS( stretch_factor_1 -                    &
                                         stretch_factor_2 ) /                  &
                                    stretch_factor_2
             
             delta_total_new = delta_l + delta_stretch_factor

!
!--                stretch_factor_1 is taken to guarantee that the stretching
!--                procedure ends as close as possible to dz_stretch_level_end.
!--                stretch_factor_2 would guarantee that the stretched dz(n) is
!--                equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old)  THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF
             
             iterations = iterations + 1
            
          ENDDO
             
       ELSEIF ( dz(n) < dz(n+1) )  THEN 
          DO  WHILE ( stretch_factor_1 <= stretch_factor_upper_limit )
                     
             stretch_factor_1 = 1.0 + iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                         &
                        dz_stretch_level_start(n) ) 
             numerator = distance*stretch_factor_1/dz(n) +                     &
                         stretch_factor_1 - distance/dz(n)
             
             l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0
             l_rounded = NINT( l )
             delta_l = ABS( l_rounded - l ) / l
             
             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 -                    &
                                        stretch_factor_2 ) /                   &
                                        stretch_factor_2
             
             delta_total_new = delta_l + delta_stretch_factor
             
!
!--          stretch_factor_1 is taken to guarantee that the stretching
!--          procedure ends as close as possible to dz_stretch_level_end.
!--          stretch_factor_2 would guarantee that the stretched dz(n) is
!--          equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old)  THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF
             
             iterations = iterations + 1
          ENDDO
          
       ELSE
          message = 'Two adjacent values of dz must be different'
          CALL inifor_abort( 'calculate_stretching_factor', message)
       ENDIF

!
!--    Check if also the second stretching factor fits into the allowed
!--    interval. If not, print a warning for the user.
       IF ( dz_stretch_factor_array_2(n) < stretch_factor_lower_limit .OR.     & 
            dz_stretch_factor_array_2(n) > stretch_factor_upper_limit )  THEN
          WRITE( message, * ) 'stretch_factor_2 = ',                           &
                                     dz_stretch_factor_array_2(n), ' which is',&
                                     ' responsible for exactly reaching& dz =',&
                                      dz(n+1), 'after a specific amount of',   & 
                                     ' grid levels& exceeds the upper',        &
                                     ' limit =', stretch_factor_upper_limit,   &
                                     ' &or lower limit = ',                    &
                                     stretch_factor_lower_limit
          CALL inifor_abort( 'calculate_stretching_factor', message )
            
       ENDIF
    ENDDO
        
 END SUBROUTINE calculate_stretching_factor

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the midpoint between two neighbouring coordinates in the given
!> coordinate vector 'z' and stores it in 'zw'.
!------------------------------------------------------------------------------!
 SUBROUTINE midpoints(z, zw)

     REAL(wp), INTENT(IN)  ::  z(0:)
     REAL(wp), INTENT(OUT) ::  zw(1:)

     INTEGER(iwp) ::  k

     DO  k = 1, UBOUND(zw, 1)
        zw(k) = 0.5_wp * (z(k-1) + z(k))
     ENDDO

 END SUBROUTINE midpoints

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Defines INFOR's IO groups.
!------------------------------------------------------------------------------!
 SUBROUTINE setup_io_groups()

    INTEGER(iwp) ::  ngroups

    ngroups = 17
    ALLOCATE( io_group_list(ngroups) )

!
!-- soil temp
    io_group_list(1) = init_io_group(                                          &
       in_files = soil_files,                                                  &
       out_vars = output_var_table(1:1),                                       &
       in_var_list = input_var_table(1:1),                                     &
       kind = 'soil-temperature'                                               &
    )

!
!-- soil water
    io_group_list(2) = init_io_group(                                          &
       in_files = soil_files,                                                  &
       out_vars = output_var_table(2:2),                                       &
       in_var_list = input_var_table(2:2),                                     &
       kind = 'soil-water'                                                     &
    )

!
!-- potential temperature, surface pressure, specific humidity including
!-- nudging and subsidence, and geostrophic winds ug, vg
    io_group_list(3) = init_io_group(                                          &
       in_files = flow_files,                                                  &
       out_vars = [output_var_table(56:64),                                    & ! internal averaged density and pressure profiles
                   output_var_table(3:8), output_var_table(9:14),              &
                   output_var_table(42:42), output_var_table(43:44),           &
                   output_var_table(49:51), output_var_table(52:54)],          &
       in_var_list = (/input_var_table(3), input_var_table(17),                & ! T, P, QV
                       input_var_table(4) /),                                  &
       kind = 'thermodynamics',                                                &
       n_output_quantities = 4_iwp                                             & ! P, Theta, Rho, qv
    )

!
!-- Moved to therodynamic io_group
    !io_group_list(4) = init_io_group(                                       &
    !   in_files = flow_files,                                               &
    !   out_vars = [output_var_table(9:14), output_var_table(52:54)],        &
    !   in_var_list = input_var_table(4:4),                                  &
    !   kind = 'scalar'                                                      &
    !)

!
!-- u and v velocity
    io_group_list(5) = init_io_group(                                          &
       in_files = flow_files,                                                  &
       out_vars = [output_var_table(15:26), output_var_table(45:46)],          &
       in_var_list = input_var_table(5:6),                                     &
       kind = 'velocities'                                                     &
    )
   
!
!-- v velocity, deprecated!
    !io_group_list(6) = init_io_group(                                       &
    !   in_files = flow_files,                                               &
    !   out_vars = output_var_table(21:26),                                  &
    !   in_var_list = input_var_table(6:6),                                  &
    !   kind = 'horizontal velocity'                                         &
    !)
    !io_group_list(6)%to_be_processed = .FALSE.
   
!
!-- w velocity and subsidence and w nudging
    io_group_list(7) = init_io_group(                                          &
       in_files = flow_files,                                                  &
       out_vars = [output_var_table(27:32), output_var_table(47:48)],          &
       in_var_list = input_var_table(7:7),                                     &
       kind = 'scalar'                                                         &
    )
!
!-- rain
    io_group_list(8) = init_io_group(                                          &
       in_files = precip_files,                                                &
       out_vars = output_var_table(33:33),                                     &
       in_var_list = input_var_table(8:8),                                     &
       kind = 'surface'                                                        &
    )
    io_group_list(8)%to_be_processed = cfg%process_precipitation
!
!-- snow
    io_group_list(9) = init_io_group(                                          &
       in_files = precip_files,                                                &
       out_vars = output_var_table(34:34),                                     &
       in_var_list = input_var_table(9:9),                                     &
       kind = 'surface'                                                        &
    )
    io_group_list(9)%to_be_processed = cfg%process_precipitation
!
!-- graupel
    io_group_list(10) = init_io_group(                                         &
       in_files = precip_files,                                                &
       out_vars = output_var_table(35:35),                                     &
       in_var_list = input_var_table(10:10),                                   &
       kind = 'surface'                                                        &
    )
    io_group_list(10)%to_be_processed = cfg%process_precipitation
!
!-- evapotranspiration
    io_group_list(11) = init_io_group(                                         &
       in_files = precip_files,                                                &
       out_vars = output_var_table(37:37),                                     &
       in_var_list = input_var_table(11:11),                                   &
       kind = 'accumulated'                                                    &
    )
    io_group_list(11)%to_be_processed = .FALSE.
!
!-- 2m air temperature
    io_group_list(12) = init_io_group(                                         &
       in_files = precip_files,                                                &
       out_vars = output_var_table(36:36),                                     &
       in_var_list = input_var_table(12:12),                                   &
       kind = 'surface'                                                        &
    )
    io_group_list(12)%to_be_processed = .FALSE.
!
!-- incoming diffusive sw flux
    io_group_list(13) = init_io_group(                                         &
       in_files = radiation_files,                                             &
       out_vars = output_var_table(38:38),                                     &
       in_var_list = input_var_table(13:13),                                   &
       kind = 'running average'                                                &
    )
    io_group_list(13)%to_be_processed = .FALSE.
!
!-- incoming direct sw flux
    io_group_list(14) = init_io_group(                                         &
       in_files = radiation_files,                                             &
       out_vars = output_var_table(39:39),                                     &
       in_var_list = input_var_table(14:14),                                   &
       kind = 'running average'                                                &
    )
    io_group_list(14)%to_be_processed = .FALSE.
!
!-- sw radiation balance
    io_group_list(15) = init_io_group(                                         &
       in_files = radiation_files,                                             &
       out_vars = output_var_table(40:40),                                     &
       in_var_list = input_var_table(15:15),                                   &
       kind = 'running average'                                                &
    )
    io_group_list(15)%to_be_processed = .FALSE.
!
!-- lw radiation balance
    io_group_list(16) = init_io_group(                                         &
       in_files = radiation_files,                                             &
       out_vars = output_var_table(41:41),                                     &
       in_var_list = input_var_table(16:16),                                   &
       kind = 'running average'                                                &
    )
    io_group_list(16)%to_be_processed = .FALSE.

!-- fog and cloud water
    io_group_list(17) = init_io_group(                                         &
       in_files = flow_files,                                                  &
       out_vars = output_var_table(65:94),                                     & !qc, qi, qr, qs, qg
       in_var_list = input_var_table(11:15),                                   & !QC, QI, QR, QS, QG
       kind = 'scalar'                                                         &
    )
    io_group_list(17)%to_be_processed = .TRUE.

    CALL validate_io_groups( io_group_list )

 END SUBROUTINE setup_io_groups


 SUBROUTINE validate_io_groups( io_group_list )
    TYPE(io_group), TARGET, INTENT(INOUT) ::  io_group_list(:)

    INTEGER(iwp)                 ::  group_idx, file_idx, outvar_idx
    TYPE(io_group), POINTER      ::  group
    CHARACTER(LEN=PATH), POINTER ::  filename
    TYPE(nc_var), POINTER        ::  outvar, invar

    DO group_idx = 1, SIZE( io_group_list )
       group => io_group_list(group_idx)

       IF (group%to_be_processed)  THEN
          DO file_idx = 1, SIZE( group%in_files )
             filename => group%in_files(file_idx)

             IF (.NOT. file_is_present( filename, 'input', message ))  THEN
                IF (group%is_optional)  THEN

!--                deactivate group, warn and continue
                   group%to_be_processed = .FALSE.
                   message = "Optional IO group '" // TRIM( group%kind ) // "' was deactivated" // &
                      " because the required file '" // TRIM( filename ) // "' was not found."
                   CALL warn( 'validate_io_groups', message )
                ELSE
!--                abort due to missing input file
                   message = "The input file '" // TRIM( filename ) // "' of the" // &
                      " mandatory IO group '" // TRIM( group%kind ) // "' was not found."
                   CALL inifor_abort( 'validate_io_groups', message  )
                ENDIF
             ELSE
                message = "Set up input file name '" // TRIM(filename) // "'"
                CALL report('validate_io_groups', message)
             ENDIF

          ENDDO
       
!--       deactivate all optional output variables that have missing inputs
          DO outvar_idx = 1, SIZE(group%out_vars)
             outvar => group%out_vars(outvar_idx)

             IF (outvar%is_optional)  THEN

!--             only check first input file, assuming files at later times have the
!--             same contents
                filename => group%in_files(LBOUND( group%in_files, 1))
                invar => group%in_var_list(outvar%input_id)

!--             if corresponding invar not present, deactivate output_var
                IF ( .NOT.netcdf_variable_present_in_file( invar%name, filename ) )  THEN
                      ! deactivate output variable variable, warn and continue
                      invar%to_be_processed = .FALSE.
                      outvar%to_be_processed = .FALSE.
                      message = "Skipping optional output variable '" // TRIM( outvar%name ) // &
                         ", corresponding input variable not available in file '" // TRIM( filename ) // "'."
                      CALL warn( 'validate_io_groups', message )
                ENDIF
             ENDIF

!--       outvar loop
          ENDDO

!--   group%to_be_processed conditional
      ENDIF

!-- groups loop
    ENDDO
    
 END SUBROUTINE validate_io_groups



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes an IO group with the list of input files, the lists of
!> input and output variables, and sets the number of output quantities.
!>
!> In this context, output quantities refers to the physical quantities required
!> to interpolate netCDF output variables, not the output variables itsleft. For
!> instance, the output variables ls_forcing_left/_right/_north/_south all rely
!> on the output quantity Theta.
!------------------------------------------------------------------------------!
 FUNCTION init_io_group(in_files, out_vars, in_var_list, kind,                 &
                        n_output_quantities) RESULT(group)
    CHARACTER(LEN=PATH), INTENT(IN) ::  in_files(:)
    CHARACTER(LEN=*), INTENT(IN)    ::  kind
    TYPE(nc_var), INTENT(IN)        ::  out_vars(:)
    TYPE(nc_var), INTENT(IN)        ::  in_var_list(:)
    INTEGER(iwp), OPTIONAL          ::  n_output_quantities

    TYPE(io_group)                  ::  group

    group%nt = SIZE(in_files)
    group%nv = SIZE(out_vars)
    group%n_inputs = SIZE(in_var_list)
    group%kind = TRIM(kind)

!
!-- Set the number of output quantities, which is used to allocate input buffers
!-- in the preprocess() routine. For passive scalars, there is a one-to-one
!-- correspondance of input and output quantities. For instance, for every water
!-- phase input variable (e.g. cloud water QC), there is a corresponding output
!-- quantity (qc, which lead to a number of output netCDF variables:
!-- init_atmosphere_qc, ls_forcing_left_qc, etc.)
!--     If more output quantities than input quantities are to be produced,
!-- n_output_quantities can be set through the optional subroutine parameter.
!-- For the 'thermodynamics' IO group, one quantity more than input variables is
!-- needed to compute all output variables of the IO group.  Concretely, in
!-- preprocess() the density is computed from T,P or PP,QV in adddition to the
!-- variables Theta, p, qv. In read_input_variables(), n_output_quantities is
!-- used to allocate the correct number of input buffers.

    IF  ( PRESENT(n_output_quantities) )  THEN
       group%n_output_quantities = n_output_quantities
    ELSE
       group%n_output_quantities = group%n_inputs
    ENDIF

    ALLOCATE(group%in_var_list(group%n_inputs))
    ALLOCATE(group%in_files(group%nt))
    ALLOCATE(group%out_vars(group%nv))

    group%in_var_list = in_var_list
    group%in_files = in_files
    group%out_vars = out_vars
    group%to_be_processed = .TRUE.

 END FUNCTION init_io_group


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocates all allocated variables.
!------------------------------------------------------------------------------!
 SUBROUTINE fini_grids()

    CALL report('fini_grids', 'Deallocating grids', cfg%debug)
    
    DEALLOCATE(x, y, z, xu, yv, zw, z_column, zw_column)

    DEALLOCATE(palm_grid%x,  palm_grid%y,  palm_grid%z,                        &
               palm_grid%xu, palm_grid%yv, palm_grid%zw,                       &
               palm_grid%clon,  palm_grid%clat,                                &
               palm_grid%clonu, palm_grid%clatu)

    DEALLOCATE(palm_intermediate%x,  palm_intermediate%y,  palm_intermediate%z, &
               palm_intermediate%xu, palm_intermediate%yv, palm_intermediate%zw,&
               palm_intermediate%clon,  palm_intermediate%clat,                &  
               palm_intermediate%clonu, palm_intermediate%clatu)

    DEALLOCATE(cosmo_grid%lon,  cosmo_grid%lat,                                &
               cosmo_grid%lonu, cosmo_grid%latv,                               &
               cosmo_grid%hfl)

 END SUBROUTINE fini_grids


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes the variable list.
!------------------------------------------------------------------------------!
 SUBROUTINE setup_variable_tables
    INTEGER(iwp)                    ::  n_invar = 0  !< number of variables in the input variable table
    INTEGER(iwp)                    ::  n_outvar = 0 !< number of variables in the output variable table
    TYPE(nc_var), POINTER           ::  var

    IF (TRIM(cfg%start_date) == '')  THEN
       message = 'Simulation start date has not been set.'
       CALL inifor_abort('setup_variable_tables', message)
    ENDIF

    nc_source_text = 'COSMO analysis from ' // TRIM(cfg%start_date)

    n_invar = 17
    n_outvar = 94
    ALLOCATE( input_var_table(n_invar) )
    ALLOCATE( output_var_table(n_outvar) )

!
!------------------------------------------------------------------------------
!- Section 1: NetCDF input variables
!------------------------------------------------------------------------------

!
!-- COSMO's soil temperature T_SO may contain the surface temperature at
!-- depth=0, which is a redundant copy the first soil layer. INIFOR uses this
!-- flag to decide if this level is present and has to be ignored in
!-- inifor_io:get_netcdf_start_and_count().
    input_var_table(:)%has_redundant_first_level = .FALSE.

    var => input_var_table(1)
    var%name = 'T_SO'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .FALSE.
    input_var_table(1)%has_redundant_first_level = has_surface_value( var, soil_files(1) )

    var => input_var_table(2)
    var%name = 'W_SO'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .FALSE.

    var => input_var_table(3)
    var%name = 'T'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(4)
    var%name = 'QV'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(5)
    var%name = 'U'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(6)
    var%name = 'V'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(7)
    var%name = 'W'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(8)
    var%name = 'RAIN_GSP'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .FALSE.

    var => input_var_table(9)
    var%name = 'SNOW_GSP'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .FALSE.

    var => input_var_table(10)
    var%name = 'GRAU_GSP'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .FALSE.

    var => input_var_table(11)
    var%name = 'QC'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(12)
    var%name = 'QI'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(13)
    var%name = 'QR'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(14)
    var%name = 'QS'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(15)
    var%name = 'QG'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

    var => input_var_table(16)
    var%name = '--'
    var%to_be_processed = .FALSE.
    var%is_upside_down = .FALSE.

    var => input_var_table(17)
    var%name = 'P'
    var%to_be_processed = .TRUE.
    var%is_upside_down = .TRUE.

!
!------------------------------------------------------------------------------
!- Section 2: NetCDF output variables
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
! Section 2.1: Realistic forcings, i.e. 3D initial and boundary conditions
!------------------------------------------------------------------------------
    output_var_table(1) = init_nc_var(                                      &
       name              = 'init_soil_t',                                   &
       std_name          = "",                                              &
       long_name         = "initial soil temperature",                      &
       units             = "K",                                             &
       kind              = "init soil",                                     &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%isc_mode) == CFG_INIT_SOIL_PROFILE)           &
    )

    output_var_table(2) = init_nc_var(                                      &
       name              = 'init_soil_m',                                   &
       std_name          = "",                                              &
       long_name         = "initial soil moisture",                         &
       units             = "m^3/m^3",                                       &
       kind              = "init soil",                                     &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%isc_mode) == CFG_INIT_SOIL_PROFILE)           &
    )

    output_var_table(3) = init_nc_var(                                      &
       name              = 'init_atmosphere_pt',                            &
       std_name          = "",                                              &
       long_name         = "initial potential temperature",                 &
       units             = "K",                                             &
       kind              = "init scalar",                                   &
       input_id          = 1_iwp,                                           & ! first in (T, p) IO group
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )

    output_var_table(4) = init_nc_var(                                      &
       name              = 'ls_forcing_left_pt',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the potential temperature", &
       units             = "K",                                             &
       kind              = "left scalar",                                   &
       input_id          = 1_iwp,                                           &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO),                &
       output_file = output_file                                            &
    )

    output_var_table(5) = init_nc_var(                                      &
       name              = 'ls_forcing_right_pt',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the potential temperature", &
       units             = "K",                                             &
       kind              = "right scalar",                                  &
       input_id          = 1_iwp,                                           &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO),                &
       output_file = output_file                                            &
    )

    output_var_table(6) = init_nc_var(                                      &
       name              = 'ls_forcing_north_pt',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the potential temperature", &
       units             = "K",                                             &
       kind              = "north scalar",                                  &
       input_id          = 1_iwp,                                           &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO),                &
       output_file = output_file                                            &
    )

    output_var_table(7) = init_nc_var(                                      &
       name              = 'ls_forcing_south_pt',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the potential temperature", &
       units             = "K",                                             &
       kind              = "south scalar",                                  &
       input_id          = 1_iwp,                                           &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO),                &
       output_file = output_file                                            &
    )

    output_var_table(8) = init_nc_var(                                      &
       name              = 'ls_forcing_top_pt',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the potential temperature", &
       units             = "K",                                             &
       kind              = "top scalar",                                    &
       input_id          = 1_iwp,                                           &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO),                &
       output_file = output_file                                            &
    )

    output_var_table(9) = init_nc_var(                                      &
       name              = 'init_atmosphere_qv',                            &
       std_name          = "",                                              &
       long_name         = "initial specific humidity",                     &
       units             = "kg/kg",                                         &
       kind              = "init scalar",                                   &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )

    output_var_table(10) = init_nc_var(                                     &
       name              = 'ls_forcing_left_qv',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the specific humidity", &
       units             = "kg/kg",                                         &
       kind              = "left scalar",                                   &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(11) = init_nc_var(                                     &
       name              = 'ls_forcing_right_qv',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the specific humidity", &
       units             = "kg/kg",                                         &
       kind              = "right scalar",                                  &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(12) = init_nc_var(                                     &
       name              = 'ls_forcing_north_qv',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the specific humidity", &
       units             = "kg/kg",                                         &
       kind              = "north scalar",                                  &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(13) = init_nc_var(                                     &
       name              = 'ls_forcing_south_qv',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the specific humidity", &
       units             = "kg/kg",                                         &
       kind              = "south scalar",                                  &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(14) = init_nc_var(                                     &
       name              = 'ls_forcing_top_qv',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the specific humidity", &
       units             = "kg/kg",                                         &
       kind              = "top scalar",                                    &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(15) = init_nc_var(                                     &
       name              = 'init_atmosphere_u',                             &
       std_name          = "",                                              &
       long_name         = "initial wind component in x direction",         &
       units             = "m/s",                                           &
       kind              = "init u",                                        &
       input_id          = 1_iwp,                                           & ! first in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = u_initial_grid,                                  &
       intermediate_grid = u_initial_intermediate,                          &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )

    output_var_table(16) = init_nc_var(                                     &
       name              = 'ls_forcing_left_u',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the wind component in x direction", &
       units             = "m/s",                                           &
       kind              = "left u",                                        &
       input_id          = 1_iwp,                                           & ! first in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = u_west_grid,                                     &
       intermediate_grid = u_west_intermediate,                             &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(17) = init_nc_var(                                     &
       name              = 'ls_forcing_right_u',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the wind component in x direction", &
       units             = "m/s",                                           &
       kind              = "right u",                                       &
       input_id          = 1_iwp,                                           & ! first in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = u_east_grid,                                     &
       intermediate_grid = u_east_intermediate,                             &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(18) = init_nc_var(                                     &
       name              = 'ls_forcing_north_u',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the wind component in x direction", &
       units             = "m/s",                                           &
       kind              = "north u",                                       &
       input_id          = 1_iwp,                                           & ! first in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = u_north_grid,                                    &
       intermediate_grid = u_north_intermediate,                            &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(19) = init_nc_var(                                     &
       name              = 'ls_forcing_south_u',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the wind component in x direction", &
       units             = "m/s",                                           &
       kind              = "south u",                                       &
       input_id          = 1_iwp,                                           & ! first in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = u_south_grid,                                    &
       intermediate_grid = u_south_intermediate,                            &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(20) = init_nc_var(                                     &
       name              = 'ls_forcing_top_u',                              &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the wind component in x direction", &
       units             = "m/s",                                           &
       kind              = "top u",                                         &
       input_id          = 1_iwp,                                           & ! first in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = u_top_grid,                                      &
       intermediate_grid = u_top_intermediate,                              &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(21) = init_nc_var(                                     &
       name              = 'init_atmosphere_v',                             &
       std_name          = "",                                              &
       long_name         = "initial wind component in y direction",         &
       units             = "m/s",                                           &
       kind              = "init v",                                        &
       input_id          = 2_iwp,                                           & ! second in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = v_initial_grid,                                  &
       intermediate_grid = v_initial_intermediate,                          &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )

    output_var_table(22) = init_nc_var(                                     &
       name              = 'ls_forcing_left_v',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the wind component in y direction", &
       units             = "m/s",                                           &
       kind              = "right v",                                       &
       input_id          = 2_iwp,                                           & ! second in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = v_west_grid,                                     &
       intermediate_grid = v_west_intermediate,                             &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(23) = init_nc_var(                                     &
       name              = 'ls_forcing_right_v',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the wind component in y direction", &
       units             = "m/s",                                           &
       kind              = "right v",                                       &
       input_id          = 2_iwp,                                           & ! second in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = v_east_grid,                                     &
       intermediate_grid = v_east_intermediate,                             &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(24) = init_nc_var(                                     &
       name              = 'ls_forcing_north_v',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the wind component in y direction", &
       units             = "m/s",                                           &
       kind              = "north v",                                       &
       input_id          = 2_iwp,                                           & ! second in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = v_north_grid,                                    &
       intermediate_grid = v_north_intermediate,                            &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(25) = init_nc_var(                                     &
       name              = 'ls_forcing_south_v',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the wind component in y direction", &
       units             = "m/s",                                           &
       kind              = "south v",                                       &
       input_id          = 2_iwp,                                           & ! second in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = v_south_grid,                                    &
       intermediate_grid = v_south_intermediate,                            &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(26) = init_nc_var(                                     &
       name              = 'ls_forcing_top_v',                              &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the wind component in y direction", &
       units             = "m/s",                                           &
       kind              = "top v",                                         &
       input_id          = 2_iwp,                                           & ! second in (U, V) I/O group
       output_file       = output_file,                                     &
       grid              = v_top_grid,                                      &
       intermediate_grid = v_top_intermediate,                              &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(27) = init_nc_var(                                     &
       name              = 'init_atmosphere_w',                             &
       std_name          = "",                                              &
       long_name         = "initial wind component in z direction",         &
       units             = "m/s",                                           &
       kind              = "init w",                                        &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = w_initial_grid,                                  &
       intermediate_grid = w_initial_intermediate,                          &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )

    output_var_table(28) = init_nc_var(                                     &
       name              = 'ls_forcing_left_w',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the wind component in z direction", &
       units             = "m/s",                                           &
       kind              = "left w",                                        &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = w_west_grid,                                     &
       intermediate_grid = w_west_intermediate,                             &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(29) = init_nc_var(                                     &
       name              = 'ls_forcing_right_w',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the wind component in z direction", &
       units             = "m/s",                                           &
       kind              = "right w",                                       &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = w_east_grid,                                     &
       intermediate_grid = w_east_intermediate,                             &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(30) = init_nc_var(                                     &
       name              = 'ls_forcing_north_w',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the wind component in z direction", &
       units             = "m/s",                                           &
       kind              = "north w",                                       &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = w_north_grid,                                    &
       intermediate_grid = w_north_intermediate,                            &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(31) = init_nc_var(                                     &
       name              = 'ls_forcing_south_w',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the wind component in z direction", &
       units             = "m/s",                                           &
       kind              = "south w",                                       &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = w_south_grid,                                    &
       intermediate_grid = w_south_intermediate,                            &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(32) = init_nc_var(                                     &
       name              = 'ls_forcing_top_w',                              &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the wind component in z direction", &
       units             = "m/s",                                           &
       kind              = "top w",                                         &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = w_top_grid,                                      &
       intermediate_grid = w_top_intermediate,                              &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )

    output_var_table(33) = init_nc_var(                                     &
       name              = 'ls_forcing_soil_rain',                          &
       std_name          = "",                                              &
       long_name         = "large-scale forcing rain",                      &
       units             = "kg/m2",                                         &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(34) = init_nc_var(                                     &
       name              = 'ls_forcing_soil_snow',                          &
       std_name          = "",                                              &
       long_name         = "large-scale forcing snow",                      &
       units             = "kg/m2",                                         &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(35) = init_nc_var(                                     &
       name              = 'ls_forcing_soil_graupel',                       &
       std_name          = "",                                              &
       long_name         = "large-scale forcing graupel",                   &
       units             = "kg/m2",                                         &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(36) = init_nc_var(                                     &
       name              = 'ls_forcing_soil_t_2m',                          &
       std_name          = "",                                              &
       long_name         = "large-scale forcing 2m air temperature",        &
       units             = "kg/m2",                                         &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(37) = init_nc_var(                                     &
       name              = 'ls_forcing_soil_evap',                          &
       std_name          = "",                                              &
       long_name         = "large-scale forcing evapo-transpiration",       &
       units             = "kg/m2",                                         &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(38) = init_nc_var(                                     &
       name              = 'rad_swd_dif_0',                                 &
       std_name          = "",                                              &
       long_name         = "incoming diffuse shortwave radiative flux at the surface", &
       units             = "W/m2",                                          &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(39) = init_nc_var(                                     &
       name              = 'rad_swd_dir_0',                                 &
       std_name          = "",                                              &
       long_name         = "incoming direct shortwave radiative flux at the surface", &
       units             = "W/m2",                                          &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(40) = init_nc_var(                                     &
       name              = 'rad_sw_bal_0',                                  &
       std_name          = "",                                              &
       long_name         = "shortwave radiation balance at the surface",    &
       units             = "W/m2",                                          &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )

    output_var_table(41) = init_nc_var(                                     &
       name              = 'rad_lw_bal_0',                                  &
       std_name          = "",                                              &
       long_name         = "longwave radiation balance at the surface",     &
       units             = "W/m2",                                          &
       kind              = "surface forcing",                               &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )
!
!------------------------------------------------------------------------------
! Section 2.2: Idealized large-scale forcings
!------------------------------------------------------------------------------
    output_var_table(42) = init_nc_var(                                     &
       name              = 'surface_forcing_surface_pressure',              &
       std_name          = "",                                              &
       long_name         = "surface pressure",                              &
       units             = "Pa",                                            &
       kind              = "time series",                                   &
       input_id          = 2_iwp,                                           & ! second in (T, p) I/O group
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate                                &
    )
    output_var_table(42)%averaging_grid => geostrophic_scalar_profile

    output_var_table(43) = init_nc_var(                                     &
       name              = 'ls_forcing_ug',                                 &
       std_name          = "",                                              &
       long_name         = "geostrophic wind (u component)",                &
       units             = "m/s",                                           &
       kind              = "geostrophic",                                   &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )

    output_var_table(44) = init_nc_var(                                     &
       name              = 'ls_forcing_vg',                                 &
       std_name          = "",                                              &
       long_name         = "geostrophic wind (v component)",                &
       units             = "m/s",                                           &
       kind              = "geostrophic",                                   &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )

    output_var_table(45) = init_nc_var(                                     &
       name              = 'nudging_u',                                     &
       std_name          = "",                                              &
       long_name         = "wind component in x direction",                 &
       units             = "m/s",                                           &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )

    output_var_table(46) = init_nc_var(                                     &
       name              = 'nudging_v',                                     &
       std_name          = "",                                              &
       long_name         = "wind component in y direction",                 &
       units             = "m/s",                                           &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )

    output_var_table(47) = init_nc_var(                                     &
       name              = 'ls_forcing_sub_w',                              &
       std_name          = "",                                              &
       long_name         = "subsidence velocity of w",                      &
       units             = "m/s",                                           &
       kind              = "large-scale w forcing",                         &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )

    output_var_table(48) = init_nc_var(                                     &
       name              = 'nudging_w',                                     &
       std_name          = "",                                              &
       long_name         = "wind component in w direction",                 &
       units             = "m/s",                                           &
       kind              = "large-scale w forcing",                         &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_w_profile,                              &
       intermediate_grid = geostrophic_w_profile                               &
    )
    output_var_table(48)%to_be_processed = ls_forcing_variables_required


    output_var_table(49) = init_nc_var(                                     &
       name              = 'ls_forcing_adv_pt',                             &
       std_name          = "",                                              &
       long_name         = "advection of potential temperature",            &
       units             = "K/s",                                           &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(49)%to_be_processed = ls_forcing_variables_required

    output_var_table(50) = init_nc_var(                                     &
       name              = 'ls_forcing_sub_pt',                             &
       std_name          = "",                                              &
       long_name         = "subsidence velocity of potential temperature",  &
       units             = "K/s",                                           &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(50)%to_be_processed = ls_forcing_variables_required

    output_var_table(51) = init_nc_var(                                     &
       name              = 'nudging_pt',                                    &
       std_name          = "",                                              &
       long_name         = "potential temperature",                         &
       units             = "K",                                             &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(51)%to_be_processed = ls_forcing_variables_required

    output_var_table(52) = init_nc_var(                                     &
       name              = 'ls_forcing_adv_qv',                             &
       std_name          = "",                                              &
       long_name         = "advection of specific humidity",                &
       units             = "kg/kg/s",                                       &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(52)%to_be_processed = ls_forcing_variables_required


    output_var_table(53) = init_nc_var(                                     &
       name              = 'ls_forcing_sub_qv',                             &
       std_name          = "",                                              &
       long_name         = "subsidence velocity of specific humidity",      &
       units             = "kg/kg/s",                                       &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(53)%to_be_processed = ls_forcing_variables_required

    output_var_table(54) = init_nc_var(                                     &
       name              = 'nudging_qv',                                    &
       std_name          = "",                                              &
       long_name         = "specific humidity",                             &
       units             = "kg/kg",                                         &
       kind              = "large-scale scalar forcing",                    &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(54)%to_be_processed = ls_forcing_variables_required

    output_var_table(55) = init_nc_var(                                     &
       name              = 'nudging_tau',                                   &
       std_name          = "",                                              &
       long_name         = "nudging relaxation time scale",                 &
       units             = "s",                                             &
       kind              = "constant scalar profile",                       &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(55)%to_be_processed = ls_forcing_variables_required


    output_var_table(56) = init_nc_var(                                     &
       name              = 'internal_density_centre',                       &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = geostrophic_scalar_profile,                         &
       intermediate_grid = geostrophic_scalar_profile                          &
    )
    output_var_table(56)%averaging_grid => geostrophic_scalar_profile
!-- Except for this one, all 'internal profile'-kind output variables are only
!-- needed for the computation of geostrophic wind components and are all
!-- enabled in init_nc_var() if they are to be computed. The averaged density
!-- profile, however, is also needed for the computation of the extrapolated
!-- surface pressure, which is why we enable it here for all INIFOR settings.
    output_var_table(56)%to_be_processed = .TRUE.


    output_var_table(57) = init_nc_var(                                     &
       name              = 'internal_density_north',                        &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = north_geostrophic_scalar_profile,                   &
       intermediate_grid = north_geostrophic_scalar_profile                    &
    )
    output_var_table(57)%averaging_grid => north_geostrophic_scalar_profile
    output_var_table(57)%to_be_processed = .NOT. cfg%ug_defined_by_user


    output_var_table(58) = init_nc_var(                                     &
       name              = 'internal_density_south',                        &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = south_geostrophic_scalar_profile,                   &
       intermediate_grid = south_geostrophic_scalar_profile                    &
    )
    output_var_table(58)%averaging_grid => south_geostrophic_scalar_profile
    output_var_table(58)%to_be_processed = .NOT. cfg%ug_defined_by_user


    output_var_table(59) = init_nc_var(                                     &
       name              = 'internal_density_east',                         &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = east_geostrophic_scalar_profile,                    &
       intermediate_grid = east_geostrophic_scalar_profile                     &
    )
    output_var_table(59)%averaging_grid => east_geostrophic_scalar_profile
    output_var_table(59)%to_be_processed = .NOT. cfg%ug_defined_by_user


    output_var_table(60) = init_nc_var(                                     &
       name              = 'internal_density_west',                         &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = west_geostrophic_scalar_profile,                    &
       intermediate_grid = west_geostrophic_scalar_profile                     &
    )
    output_var_table(60)%averaging_grid => west_geostrophic_scalar_profile
    output_var_table(60)%to_be_processed = .NOT. cfg%ug_defined_by_user

    output_var_table(61) = init_nc_var(                                     &
       name              = 'internal_pressure_north',                       &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = north_geostrophic_scalar_profile,                   &
       intermediate_grid = north_geostrophic_scalar_profile                    &
    )
    output_var_table(61)%averaging_grid => north_geostrophic_scalar_profile
    output_var_table(61)%to_be_processed = .NOT. cfg%ug_defined_by_user


    output_var_table(62) = init_nc_var(                                     &
       name              = 'internal_pressure_south',                       &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = south_geostrophic_scalar_profile,                   &
       intermediate_grid = south_geostrophic_scalar_profile                    &
    )
    output_var_table(62)%averaging_grid => south_geostrophic_scalar_profile
    output_var_table(62)%to_be_processed = .NOT. cfg%ug_defined_by_user


    output_var_table(63) = init_nc_var(                                     &
       name              = 'internal_pressure_east',                        &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = east_geostrophic_scalar_profile,                    &
       intermediate_grid = east_geostrophic_scalar_profile                     &
    )
    output_var_table(63)%averaging_grid => east_geostrophic_scalar_profile
    output_var_table(63)%to_be_processed = .NOT. cfg%ug_defined_by_user


    output_var_table(64) = init_nc_var(                                     &
       name              = 'internal_pressure_west',                        &
       std_name          = "",                                              &
       long_name         = "",                                              &
       units             = "",                                              &
       kind              = "internal profile",                              &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = west_geostrophic_scalar_profile,                    &
       intermediate_grid = west_geostrophic_scalar_profile                     &
    )
    output_var_table(64)%averaging_grid => west_geostrophic_scalar_profile
    output_var_table(64)%to_be_processed = .NOT. cfg%ug_defined_by_user

    output_var_table(65) = init_nc_var(                                     &
       name              = 'init_atmosphere_qc',                            &
       std_name          = "",                                              &
       long_name         = "initial cloud water mixture fraction",          &
       units             = "kg/kg",                                         &
       kind              = "init scalar",                                   &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )
    output_var_table(65)%is_optional = .TRUE.

    output_var_table(66) = init_nc_var(                                     &
       name              = 'ls_forcing_left_qc',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the cloud water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "left scalar",                                   &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(66)%is_optional = .TRUE.

    output_var_table(67) = init_nc_var(                                     &
       name              = 'ls_forcing_right_qc',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the cloud water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "right scalar",                                  &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(67)%is_optional = .TRUE.

    output_var_table(68) = init_nc_var(                                     &
       name              = 'ls_forcing_north_qc',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the cloud water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "north scalar",                                  &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(68)%is_optional = .TRUE.

    output_var_table(69) = init_nc_var(                                     &
       name              = 'ls_forcing_south_qc',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the cloud water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "south scalar",                                  &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(69)%is_optional = .TRUE.

    output_var_table(70) = init_nc_var(                                     &
       name              = 'ls_forcing_top_qc',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the cloud water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "top scalar",                                    &
       input_id          = 1_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(70)%is_optional = .TRUE.

    output_var_table(71) = init_nc_var(                                     &
       name              = 'init_atmosphere_qi',                            &
       std_name          = "",                                              &
       long_name         = "initial cloud ice mixture fraction",            &
       units             = "kg/kg",                                         &
       kind              = "init scalar",                                   &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )
    output_var_table(71)%is_optional = .TRUE.

    output_var_table(72) = init_nc_var(                                     &
       name              = 'ls_forcing_left_qi',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the cloud ice mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "left scalar",                                   &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(72)%is_optional = .TRUE.

    output_var_table(73) = init_nc_var(                                     &
       name              = 'ls_forcing_right_qi',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the cloud ice mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "right scalar",                                  &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(73)%is_optional = .TRUE.

    output_var_table(74) = init_nc_var(                                     &
       name              = 'ls_forcing_north_qi',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the cloud ice mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "north scalar",                                  &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(74)%is_optional = .TRUE.

    output_var_table(75) = init_nc_var(                                     &
       name              = 'ls_forcing_south_qi',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the cloud ice mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "south scalar",                                  &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(75)%is_optional = .TRUE.

    output_var_table(76) = init_nc_var(                                     &
       name              = 'ls_forcing_top_qi',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the cloud ice mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "top scalar",                                    &
       input_id          = 2_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(76)%is_optional = .TRUE.

    output_var_table(77) = init_nc_var(                                     &
       name              = 'init_atmosphere_qr',                            &
       std_name          = "",                                              &
       long_name         = "initial rain water mixture fraction",           &
       units             = "kg/kg",                                         &
       kind              = "init scalar",                                   &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )
    output_var_table(77)%is_optional = .TRUE.

    output_var_table(78) = init_nc_var(                                     &
       name              = 'ls_forcing_left_qr',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the rain water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "left scalar",                                   &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(78)%is_optional = .TRUE.

    output_var_table(79) = init_nc_var(                                     &
       name              = 'ls_forcing_right_qr',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the rain water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "right scalar",                                  &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(79)%is_optional = .TRUE.

    output_var_table(80) = init_nc_var(                                     &
       name              = 'ls_forcing_north_qr',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the rain water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "north scalar",                                  &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(80)%is_optional = .TRUE.

    output_var_table(81) = init_nc_var(                                     &
       name              = 'ls_forcing_south_qr',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the rain water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "south scalar",                                  &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(81)%is_optional = .TRUE.

    output_var_table(82) = init_nc_var(                                     &
       name              = 'ls_forcing_top_qr',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the rain water mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "top scalar",                                    &
       input_id          = 3_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(82)%is_optional = .TRUE.

    output_var_table(83) = init_nc_var(                                     &
       name              = 'init_atmosphere_qs',                            &
       std_name          = "",                                              &
       long_name         = "initial snow mixture fraction",                 &
       units             = "kg/kg",                                         &
       kind              = "init scalar",                                   &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )
    output_var_table(83)%is_optional = .TRUE.

    output_var_table(84) = init_nc_var(                                     &
       name              = 'ls_forcing_left_qs',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the snow mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "left scalar",                                   &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(84)%is_optional = .TRUE.

    output_var_table(85) = init_nc_var(                                     &
       name              = 'ls_forcing_right_qs',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the snow mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "right scalar",                                  &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(85)%is_optional = .TRUE.

    output_var_table(86) = init_nc_var(                                     &
       name              = 'ls_forcing_north_qs',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the snow mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "north scalar",                                  &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(86)%is_optional = .TRUE.

    output_var_table(87) = init_nc_var(                                     &
       name              = 'ls_forcing_south_qs',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the snow mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "south scalar",                                  &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(87)%is_optional = .TRUE.

    output_var_table(88) = init_nc_var(                                     &
       name              = 'ls_forcing_top_qs',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the snow mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "top scalar",                                    &
       input_id          = 4_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(88)%is_optional = .TRUE.

    output_var_table(89) = init_nc_var(                                     &
       name              = 'init_atmosphere_qg',                            &
       std_name          = "",                                              &
       long_name         = "initial graupel mixture fraction",              &
       units             = "kg/kg",                                         &
       kind              = "init scalar",                                   &
       input_id          = 5_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = palm_grid,                                       &
       intermediate_grid = palm_intermediate,                               &
       is_profile = (TRIM(cfg%ic_mode) == CFG_INIT_PROFILE)                 &
    )
    output_var_table(89)%is_optional = .TRUE.

    output_var_table(90) = init_nc_var(                                     &
       name              = 'ls_forcing_left_qg',                            &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for left model boundary for the graupel mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "left scalar",                                   &
       input_id          = 5_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_west_grid,                               &
       intermediate_grid = scalars_west_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(90)%is_optional = .TRUE.

    output_var_table(91) = init_nc_var(                                     &
       name              = 'ls_forcing_right_qg',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for right model boundary for the graupel mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "right scalar",                                  &
       input_id          = 5_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_east_grid,                               &
       intermediate_grid = scalars_east_intermediate,                       &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(91)%is_optional = .TRUE.

    output_var_table(92) = init_nc_var(                                     &
       name              = 'ls_forcing_north_qg',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for north model boundary for the graupel mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "north scalar",                                  &
       input_id          = 5_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_north_grid,                              &
       intermediate_grid = scalars_north_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(92)%is_optional = .TRUE.

    output_var_table(93) = init_nc_var(                                     &
       name              = 'ls_forcing_south_qg',                           &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for south model boundary for the graupel mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "south scalar",                                  &
       input_id          = 5_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_south_grid,                              &
       intermediate_grid = scalars_south_intermediate,                      &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(93)%is_optional = .TRUE.

    output_var_table(94) = init_nc_var(                                     &
       name              = 'ls_forcing_top_qg',                             &
       std_name          = "",                                              &
       long_name         = "large-scale forcing for top model boundary for the graupel mixture fraction", &
       units             = "kg/kg",                                         &
       kind              = "top scalar",                                    &
       input_id          = 5_iwp,                                           &
       output_file       = output_file,                                     &
       grid              = scalars_top_grid,                                &
       intermediate_grid = scalars_top_intermediate,                        &
       is_profile = (TRIM(cfg%bc_mode) == CFG_FORCING_HOMO)                 &
    )
    output_var_table(94)%is_optional = .TRUE.

!
!-- Attributes shared among all variables
    output_var_table(:)%source = nc_source_text


 END SUBROUTINE setup_variable_tables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes an nc_var varible with the given parameters. The 'kind'
!> parameter is used to infer the correct netCDF IDs and the level of detail,
!> 'lod', as defined by the PALM-4U input data standard.
!------------------------------------------------------------------------------!
 FUNCTION init_nc_var(name, std_name, long_name, units, kind, input_id,     &
                      grid, intermediate_grid, output_file, is_profile)     &
    RESULT(var)

    CHARACTER(LEN=*), INTENT(IN)      ::  name, std_name, long_name, units, kind
    INTEGER(iwp), INTENT(IN)          ::  input_id
    TYPE(grid_definition), INTENT(IN), TARGET ::  grid, intermediate_grid
    TYPE(nc_file), INTENT(IN)         ::  output_file
    LOGICAL, INTENT(IN), OPTIONAL     ::  is_profile

    CHARACTER(LEN=LNAME)              ::  out_var_kind 
    TYPE(nc_var)                      ::  var

    out_var_kind = TRIM(kind)

    IF (PRESENT(is_profile))  THEN
       IF (is_profile)  out_var_kind = TRIM(kind) // ' profile'
    ENDIF

    var%name              = name
    var%standard_name     = std_name
    var%long_name         = long_name
    var%units             = units
    var%kind              = TRIM(out_var_kind)
    var%input_id          = input_id
    var%nt                = SIZE (output_file%time)
    var%grid              => grid
    var%intermediate_grid => intermediate_grid

    SELECT CASE( TRIM(out_var_kind) )

!
!--    TODO: Using global module variables 'init_variables_required' and
!--    TODO: 'boundary_variables_required'. Encapsulate in settings type
!--    TODO: and pass into init_nc_var.
       CASE( 'init soil' )
          var%nt              = 1
          var%lod             = 2
          var%ndim            = 3
          var%dimids(1:3)     = output_file%dimids_soil
          var%dimvarids(1:3)  = output_file%dimvarids_soil
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_2d"

       CASE( 'init scalar' )
          var%nt              = 1
          var%lod             = 2
          var%ndim            = 3
          var%dimids(1:3)     = output_file%dimids_scl
          var%dimvarids(1:3)  = output_file%dimvarids_scl
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'init u' )
          var%nt              = 1
          var%lod             = 2
          var%ndim            = 3
          var%dimids(1)       = output_file%dimids_vel(1)
          var%dimids(2)       = output_file%dimids_scl(2)
          var%dimids(3)       = output_file%dimids_scl(3)
          var%dimvarids(1)    = output_file%dimvarids_vel(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(2)
          var%dimvarids(3)    = output_file%dimvarids_scl(3)
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'init v' )
          var%nt              = 1
          var%lod             = 2
          var%ndim            = 3
          var%dimids(1)       = output_file%dimids_scl(1)
          var%dimids(2)       = output_file%dimids_vel(2)
          var%dimids(3)       = output_file%dimids_scl(3)
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_vel(2)
          var%dimvarids(3)    = output_file%dimvarids_scl(3)
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'init w' )
          var%nt              = 1
          var%lod             = 2
          var%ndim            = 3
          var%dimids(1)       = output_file%dimids_scl(1)
          var%dimids(2)       = output_file%dimids_scl(2)
          var%dimids(3)       = output_file%dimids_vel(3)
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(2)
          var%dimvarids(3)    = output_file%dimvarids_vel(3)
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'init soil profile' )
          var%nt              = 1
          var%lod             = 1
          var%ndim            = 1
          var%dimids(1)       = output_file%dimids_soil(3)
          var%dimvarids(1)    = output_file%dimvarids_soil(3)
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average levels"
          var%averaging_grid  => averaged_soil_profile

       CASE( 'init scalar profile', 'init u profile', 'init v profile' )
          var%nt              = 1
          var%lod             = 1
          var%ndim            = 1
          var%dimids(1)       = output_file%dimids_scl(3)    !z
          var%dimvarids(1)    = output_file%dimvarids_scl(3) !z
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average profile"
          var%averaging_grid  => averaged_scalar_profile

       CASE( 'init w profile')
          var%nt              = 1
          var%lod             = 1
          var%ndim            = 1
          var%dimids(1)       = output_file%dimids_vel(3)    !z
          var%dimvarids(1)    = output_file%dimvarids_vel(3) !z
          var%to_be_processed = init_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average profile"
          var%averaging_grid  => averaged_w_profile

       CASE( 'surface forcing' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1:2)     = output_file%dimids_soil(1:2)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1:2)  = output_file%dimvarids_soil(1:2)
          var%to_be_processed = surface_forcing_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_2d"

       CASE( 'left scalar', 'right scalar')
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1)       = output_file%dimids_scl(2)
          var%dimids(2)       = output_file%dimids_scl(3)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(2)
          var%dimvarids(2)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'north scalar', 'south scalar')
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1)       = output_file%dimids_scl(1)
          var%dimids(2)       = output_file%dimids_scl(3)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'top scalar', 'top w' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1)       = output_file%dimids_scl(1)
          var%dimids(2)       = output_file%dimids_scl(2)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(2)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'left u', 'right u' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1)       = output_file%dimids_scl(2)
          var%dimids(2)       = output_file%dimids_scl(3)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(2)
          var%dimvarids(2)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'north u', 'south u' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_vel(1) !x
          var%dimids(2)       = output_file%dimids_scl(3) !z
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_vel(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'top u' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_vel(1) !x
          var%dimids(2)       = output_file%dimids_scl(2) !z
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_vel(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(2)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'left v', 'right v' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1)       = output_file%dimids_vel(2)
          var%dimids(2)       = output_file%dimids_scl(3)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_vel(2)
          var%dimvarids(2)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'north v', 'south v' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(1) !x
          var%dimids(2)       = output_file%dimids_scl(3) !z
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'top v' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(1) !x
          var%dimids(2)       = output_file%dimids_vel(2) !z
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_vel(2)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'left w', 'right w')
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time
          var%dimids(1)       = output_file%dimids_scl(2)
          var%dimids(2)       = output_file%dimids_vel(3)
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(2)
          var%dimvarids(2)    = output_file%dimvarids_vel(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'north w', 'south w' )
          var%lod             = 2
          var%ndim            = 3
          var%dimids(3)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(1) !x
          var%dimids(2)       = output_file%dimids_vel(3) !z
          var%dimvarids(3)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(1)
          var%dimvarids(2)    = output_file%dimvarids_vel(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "interpolate_3d"

       CASE( 'left scalar profile', 'right scalar profile',                    &
             'north scalar profile', 'south scalar profile',                   &
             'left u profile', 'right u profile',                              &
             'north u profile', 'south u profile',                             &
             'left v profile', 'right v profile',                              &
             'north v profile', 'south v profile' )
          var%lod             = 1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average profile"
          var%averaging_grid  => averaged_scalar_profile

       CASE( 'top scalar profile', 'top u profile', 'top v profile' )
          var%lod             = 1
          var%ndim            = 1
          var%dimids(1)       = output_file%dimid_time    !t
          var%dimvarids(1)    = output_file%dimvarid_time
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average profile"
          var%averaging_grid  => averaged_scalar_top_point

       CASE( 'left w profile', 'right w profile',                              &
             'north w profile', 'south w profile' )
          var%lod             = 1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_vel(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_vel(3)
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average profile"
          var%averaging_grid  => averaged_w_profile

       CASE( 'top w profile' )
          var%lod             = 1
          var%ndim            = 1
          var%dimids(1)       = output_file%dimid_time    !t
          var%dimvarids(1)    = output_file%dimvarid_time
          var%to_be_processed = boundary_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average profile"
          var%averaging_grid  => averaged_w_top_point

       CASE( 'time series' )
          var%lod             = 0
          var%ndim            = 1
          var%dimids(1)       = output_file%dimid_time    !t
          var%dimvarids(1)    = output_file%dimvarid_time
          var%to_be_processed = .TRUE.
          var%is_internal     = .FALSE.
          var%task            = "average profile"

       CASE( 'constant scalar profile' )
          var%lod             = 1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(3)
          var%to_be_processed = .TRUE.
          var%is_internal     = .FALSE.
          var%task            = "set profile"

       CASE( 'large-scale scalar forcing' )
          var%lod             = 1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(3)
          var%to_be_processed = ls_forcing_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average large-scale profile"

       CASE( 'geostrophic' )
          var%lod             = 1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(3)
          var%to_be_processed = ls_forcing_variables_required .OR. cfg%ug_defined_by_user
          var%is_internal     = .FALSE.
          var%task            = "geostrophic winds"

       CASE( 'large-scale w forcing' )
          var%lod             = 1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_vel(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_vel(3)
          var%to_be_processed = ls_forcing_variables_required
          var%is_internal     = .FALSE.
          var%task            = "average large-scale profile"

       CASE( 'internal profile' )
          var%lod             = -1
          var%ndim            = 2
          var%dimids(2)       = output_file%dimid_time    !t
          var%dimids(1)       = output_file%dimids_scl(3) !z
          var%dimvarids(2)    = output_file%dimvarid_time
          var%dimvarids(1)    = output_file%dimvarids_scl(3)
          var%to_be_processed = ls_forcing_variables_required
          var%is_internal     = .TRUE.
          var%task            = "internal profile"

       CASE DEFAULT
           message = "Variable kind '" // TRIM(out_var_kind) // "' not recognized."
           CALL inifor_abort ('init_nc_var', message)

    END SELECT

 END FUNCTION init_nc_var


 SUBROUTINE fini_variables()

    CALL report('fini_variables', 'Deallocating variable table', cfg%debug)
    DEALLOCATE( input_var_table )

 END SUBROUTINE fini_variables


 SUBROUTINE fini_io_groups()

    CALL report('fini_io_groups', 'Deallocating IO groups', cfg%debug)
    DEALLOCATE( io_group_list )

 END SUBROUTINE fini_io_groups


 SUBROUTINE fini_file_lists()
    
    CALL report('fini_file_lists', 'Deallocating file lists', cfg%debug)
    DEALLOCATE( flow_files, soil_files, radiation_files, precip_files )

 END SUBROUTINE fini_file_lists


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Carries out any physical conversion of the quantities in the given input
!> buffer needed to obtain the quantity required by PALM-4U. For instance, 
!> velocities are rotated to the PALM-4U coordinate system and the potential
!> temperature is computed from the absolute temperature and pressure.
!>
!> Note, that the preprocessing does not include any grid change. The result
!> array will match a COSMO-DE scalar array.
!------------------------------------------------------------------------------!
 SUBROUTINE preprocess( group, input_buffer, cosmo_grid )

    TYPE(io_group), INTENT(INOUT), TARGET       ::  group
    TYPE(container), INTENT(INOUT), ALLOCATABLE ::  input_buffer(:)
    TYPE(grid_definition), INTENT(IN)           ::  cosmo_grid
    
    REAL(wp), ALLOCATABLE                       ::  basic_state_pressure(:)
    TYPE(container), ALLOCATABLE                ::  preprocess_buffer(:)
    INTEGER(iwp)                                ::  i, j, k
    INTEGER(iwp)                                ::  nx, ny, nz
    
    input_buffer(:)%is_preprocessed = .FALSE.
     
    SELECT CASE( group%kind )
       
       CASE( 'velocities' )
!
!--       Allocate a compute buffer with the same number of arrays as the input
          ALLOCATE( preprocess_buffer( SIZE(input_buffer) ) )

!
!--       Allocate u and v arrays with scalar dimensions
          nx = SIZE(input_buffer(1)%array, 1)
          ny = SIZE(input_buffer(1)%array, 2)
          nz = SIZE(input_buffer(1)%array, 3)
          ALLOCATE( preprocess_buffer(1)%array(nx, ny, nz) ) ! u buffer
          ALLOCATE( preprocess_buffer(2)%array(nx, ny, nz) ) ! v buffer

          CALL log_runtime('time', 'alloc')

!
!--       interpolate U and V to centres
          CALL centre_velocities( u_face = input_buffer(1)%array,            &
                                  v_face = input_buffer(2)%array,            &
                                  u_centre = preprocess_buffer(1)%array,     &
                                  v_centre = preprocess_buffer(2)%array )
          
          cfg%rotation_method = 'rotated-pole'
          SELECT CASE(cfg%rotation_method)

             CASE('rotated-pole')
!            
!--             rotate U and V to PALM-4U orientation and overwrite U and V with
!--             rotated velocities
                DO  k = 1, nz
                DO  j = 1, ny
                DO  i = 1, nx
                   CALL uv2uvrot( urot = preprocess_buffer(1)%array(i,j,k),     &
                                  vrot = preprocess_buffer(2)%array(i,j,k),     &
                                  rlat = cosmo_grid%lat(j-1),                   &
                                  rlon = cosmo_grid%lon(i-1),                   &
                                  pollat = phi_cn,                                &
                                  pollon = lambda_cn,                             &
                                  u = input_buffer(1)%array(i,j,k),             &
                                  v = input_buffer(2)%array(i,j,k) )
                ENDDO
                ENDDO
                ENDDO
             
             CASE DEFAULT
                message = "Rotation method '" // TRIM(cfg%rotation_method) //   &
                   "' not recognized."
                CALL inifor_abort('preprocess', message)

          END SELECT

          input_buffer(1)%array(1,:,:) = 0.0_wp
          input_buffer(2)%array(1,:,:) = 0.0_wp
          input_buffer(1)%array(:,1,:) = 0.0_wp
          input_buffer(2)%array(:,1,:) = 0.0_wp

          input_buffer(1:2)%is_preprocessed = .TRUE.
          CALL log_runtime('time', 'comp')

          DEALLOCATE( preprocess_buffer )
          CALL log_runtime('time', 'alloc')

          message = "Input buffers for group '" // TRIM(group%kind) // "'"//&
             " preprocessed sucessfully."
          CALL report('preprocess', message)
       
       CASE( 'thermodynamics' ) ! T, P, QV
          nx = SIZE(input_buffer(1)%array, 1)
          ny = SIZE(input_buffer(1)%array, 2)
          nz = SIZE(input_buffer(1)%array, 3)

!
!--       Compute absolute pressure if presure perturbation has been read in.
          IF ( TRIM(group%in_var_list(2)%name) == 'PP' )  THEN
             message = "Absolute pressure, P, not available, " //              &
                       "computing from pressure preturbation PP."
             CALL report('preprocess', message)

             ALLOCATE( basic_state_pressure(1:nz) )
             CALL log_runtime('time', 'alloc')

             DO  j = 1, ny
             DO  i = 1, nx

                CALL get_basic_state( cosmo_grid%hfl(i,j,:), BETA, P_SL, T_SL, &
                                      RD, G, basic_state_pressure )

!
!--             Overwrite pressure perturbation with absolute pressure. HECTO
!--             converts pressure perturbation from hPa to Pa.
                input_buffer (2)%array(i,j,:) =                              &
                   HECTO * input_buffer (2)%array(i,j,:) +                   &
                   basic_state_pressure(:)

             ENDDO
             ENDDO
             CALL log_runtime('time', 'comp')

             DEALLOCATE( basic_state_pressure )
             CALL log_runtime('time', 'alloc')

             group%in_var_list(2)%name = 'P'

          ENDIF
!
!--       mark pressure as preprocessed
          input_buffer(2)%is_preprocessed = .TRUE.

!
!--       Copy temperature to the last input buffer array
          ALLOCATE(                                                            &
              input_buffer( group%n_output_quantities )%array (nx, ny, nz) &
          )
          CALL log_runtime('time', 'alloc')
          input_buffer(group%n_output_quantities)%array(:,:,:) =           &
              input_buffer(1)%array(:,:,:)

!
!--       Convert absolute in place to potential temperature
          CALL potential_temperature(                                          &
             t = input_buffer(1)%array(:,:,:),                               &
             p = input_buffer(2)%array(:,:,:),                               &
             p_ref = P_REF,                                                    &
             r = RD_PALM,                                                      &
             cp = CP_PALM                                                      &
          )

!
!--       mark potential temperature as preprocessed
          input_buffer(1)%is_preprocessed = .TRUE.

!
!--       Convert temperature copy to density
          CALL moist_density(                                                  &
             t_rho = input_buffer(group%n_output_quantities)%array(:,:,:), &
             p = input_buffer(2)%array(:,:,:),                               &
             qv = input_buffer(3)%array(:,:,:),                              &
             rd = RD,                                                          &
             rv = RV                                                           &
          )

!
!--       mark qv as preprocessed
          input_buffer(3)%is_preprocessed = .TRUE.

!
!--       mark density as preprocessed
          input_buffer(group%n_output_quantities)%is_preprocessed = .TRUE.


          message = "Input buffers for group '" // TRIM(group%kind) // "'"//&
             " preprocessed sucessfully."
          CALL report('preprocess', message)
       
       CASE( 'scalar' ) ! S or W
          input_buffer(:)%is_preprocessed = .TRUE.

       CASE( 'soil-temperature' ) ! 
          
          CALL fill_water_cells(soiltyp, input_buffer(1)%array, &
                                SIZE(input_buffer(1)%array, 3, kind=iwp), &
                                FILL_ITERATIONS)
          input_buffer(:)%is_preprocessed = .TRUE.

       CASE( 'soil-water' ) ! 

          CALL fill_water_cells(soiltyp, input_buffer(1)%array, &
                                SIZE(input_buffer(1)%array, 3, kind=iwp), &
                                FILL_ITERATIONS)

          nx = SIZE(input_buffer(1)%array, 1)
          ny = SIZE(input_buffer(1)%array, 2)
          nz = SIZE(input_buffer(1)%array, 3)

          DO  k = 1, nz
          DO  j = 1, ny
          DO  i = 1, nx
             input_buffer(1)%array(i,j,k) =                                  &
                 input_buffer(1)%array(i,j,k) * d_depth_rho_inv(k)
          ENDDO
          ENDDO
          ENDDO

          message = "Converted soil water from [kg/m^2] to [m^3/m^3]"
          CALL report('preprocess', message)

          input_buffer(:)%is_preprocessed = .TRUE.

       CASE( 'surface' ) ! 
          input_buffer(:)%is_preprocessed = .TRUE.

       CASE DEFAULT
          message = "IO group kind '" // TRIM(group%kind) // "' is not supported."
          CALL inifor_abort('prerpocess', message)

    END SELECT
    CALL log_runtime('time', 'comp')

 END SUBROUTINE preprocess


 END MODULE inifor_grid
#endif
