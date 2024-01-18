!> @file src/inifor_types.f90
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
! Description:
! ------------
!> The types module provides derived data types used in INIFOR.
!------------------------------------------------------------------------------!
 MODULE inifor_types
 
 USE inifor_defs,                                                              &
    ONLY:  DATE, PATH, SNAME, LNAME, iwp, wp

#if defined ( __netcdf )
 USE netcdf,                                                                   &
    ONLY:  NF90_MAX_VAR_DIMS, NF90_MAX_NAME
#endif

 IMPLICIT NONE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Contaner for the INIFOR command-line configuration
!------------------------------------------------------------------------------!
 TYPE inifor_config
    CHARACTER(LEN=DATE)  ::  start_date           !< String of the FORMAT YYYYMMDDHH indicating the start of the intended PALM-4U simulation

    CHARACTER(LEN=PATH)  ::  input_path           !< Path to the input data file directory
    CHARACTER(LEN=PATH)  ::  hhl_file             !< Path to the file containing the COSMO-DE HHL variable (height of half layers, i.e. vertical cell faces)
    CHARACTER(LEN=PATH)  ::  namelist_file        !< Path to the PALM-4U namelist file
    CHARACTER(LEN=PATH)  ::  output_file          !< Path to the INIFOR output file (i.e. PALM-4U dynamic driver')
    CHARACTER(LEN=PATH)  ::  soiltyp_file         !< Path to the file containing the COSMO-DE SOILTYP variable (map of COSMO-DE soil types)
    CHARACTER(LEN=PATH)  ::  static_driver_file   !< Path to the file containing the COSMO-DE SOILTYP variable (map of COSMO-DE soil types)

    CHARACTER(LEN=SNAME) ::  flow_prefix          !< Prefix of flow input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  input_prefix         !< Prefix of all input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  radiation_prefix     !< Prefix of radiation input files, e.g 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  soil_prefix          !< Prefix of soil input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  precipitation_prefix !< Prefix of input files for precipitation forcing, e.g 'laf' for COSMO-DE analyses

    CHARACTER(LEN=SNAME) ::  averaging_mode       !< destinguishes between level-based and heigh-based averaging
    CHARACTER(LEN=SNAME) ::  bc_mode              !< destinguishes realistic and idealistic forcing
    CHARACTER(LEN=SNAME) ::  ic_mode              !< destinguishes volume and profile initialization
    CHARACTER(LEN=SNAME) ::  isc_mode             !< destinguishes volume and profile soil initialization
    CHARACTER(LEN=SNAME) ::  rotation_method      !< selects method for velocity rotation

    REAL(wp)             ::  p0                   !< manually specified surface pressure [Pa]
    REAL(wp)             ::  ug                   !< manually spefied geostrophic wind component in x direction [m/s]
    REAL(wp)             ::  vg                   !< manually spefied geostrophic wind component in y direction [m/s]
    REAL(wp)             ::  z0                   !< elevation of the PALM-4U domain above sea level [m]
    REAL(wp)             ::  averaging_angle      !< latitudal and longitudal width of averaging regions [deg]
    
    LOGICAL              ::  debug                       !< indicates whether --debug option was given
    LOGICAL              ::  flow_prefix_is_set          !< indicates whether the flow prefix was set manually
    LOGICAL              ::  input_prefix_is_set         !< indicates whether the input prefix was set manually
    LOGICAL              ::  map_terrain                 !< indicates whether mesoscale grid should be mapped to microscale terrain
    LOGICAL              ::  p0_is_set                   !< indicates whether p0 was set manually
    LOGICAL              ::  radiation_prefix_is_set     !< indicates whether the radiation prefix was set manually
    LOGICAL              ::  soil_prefix_is_set          !< indicates whether the soil prefix was set manually
    LOGICAL              ::  precipitation_prefix_is_set !< indicates whether the precipitation prefix was set manually
    LOGICAL              ::  process_precipitation       !< indicates whether precipitation should be processed
    LOGICAL              ::  static_driver_is_set        !< indicates whether a static driver was given
    LOGICAL              ::  ug_defined_by_user          !< indicates whether ug was set manually
    LOGICAL              ::  vg_defined_by_user          !< indicates whether vg was set manually
    LOGICAL              ::  z0_is_set                   !< indicates whether z0 was set manually
 END TYPE inifor_config


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Container for grid data, in partucular coordinates, interpolation neighbours
!> and weights
!------------------------------------------------------------------------------!
 TYPE grid_definition
    CHARACTER(LEN=SNAME)  ::  name(3)       !< names of the grid dimensions, e.g. (/'x', 'y', 'z'/) or (/'latitude', 'longitude', 'height'/)
    CHARACTER(LEN=SNAME)  ::  kind          !< names of the grid dimensions, e.g. (/'x', 'y', 'z'/) or (/'latitude', 'longitude', 'height'/)
    INTEGER(iwp)               ::  k_min         !< Index of lowest PALM grid level that is not cut by local COSMO orography; vertically separates interpolation and extrapolation region.
    INTEGER(iwp)               ::  nx            !< number of gridpoints in the first dimension
    INTEGER(iwp)               ::  ny            !< number of gridpoints in the second dimension
    INTEGER(iwp)               ::  nz            !< number of gridpoints in the third dimension, used for PALM points
    INTEGER(iwp)               ::  nlev          !< number of COSMO grid levels
    INTEGER(iwp)               ::  n_columns     !< number of averaging columns of the source grid 
    INTEGER(iwp), ALLOCATABLE  ::  ii(:,:,:)     !< Given a point (i,j,k) in the PALM-4U grid, ii(i,j,l) gives the x index of the l'th horizontl neighbour on the COSMO-DE grid.
    INTEGER(iwp), ALLOCATABLE  ::  jj(:,:,:)     !< Given a point (i,j,k) in the PALM-4U grid, jj(i,j,l) gives the y index of the l'th horizontl neighbour on the COSMO-DE grid.
    INTEGER(iwp), ALLOCATABLE  ::  kk(:,:,:,:)   !< Given a point (i,j,k) in the PALM-4U grid, kk(i,j,k,l) gives the z index of the l'th vertical neighbour in the intermediate grid.
    INTEGER(iwp), ALLOCATABLE  ::  iii(:)        !< profile averaging neighbour indices 
    INTEGER(iwp), ALLOCATABLE  ::  jjj(:)        !< profile averaging neighbour indices 
    INTEGER(iwp), ALLOCATABLE  ::  kkk(:,:,:)    !< indices of vertical interpolation neightbours, kkk(<source column>, <PALM k level>, <neighbour index>)
    REAL(wp)              ::  lx            !< domain length in the first dimension [m]
    REAL(wp)              ::  ly            !< domain length in the second dimension [m]
    REAL(wp)              ::  x0            !< x coordinate of PALM-4U domain projection centre, i.e. location of zero distortion
    REAL(wp)              ::  y0            !< y coordinate of PALM-4U domain projection centre, i.e. location of zwro distortion
    REAL(wp)              ::  z0            !< displacement of the coordinate origin above sea level [m]
    REAL(wp), ALLOCATABLE ::  x(:)          !< coordinates of cell centers in x direction [m]
    REAL(wp), ALLOCATABLE ::  y(:)          !< coordinates of cell centers in y direction [m]
    REAL(wp), POINTER     ::  z(:)          !< coordinates of cell centers in z direction [m]
    REAL(wp), ALLOCATABLE ::  intermediate_h(:,:,:) !< heights grid point for intermediate grids [m]
    REAL(wp), ALLOCATABLE ::  sigma(:,:,:)  !< sigma-z coordinate of the original (unmapped) mesoscale grid, = (z - zt)/(zs - zt) [-]
    REAL(wp), POINTER     ::  cosmo_h(:,:,:)!< pointer to appropriate COSMO level heights (scalar/w) [m]
    REAL(wp), POINTER     ::  hhl(:,:,:)    !< heights of half layers (cell faces) above sea level in COSMO-DE, read in from 
    REAL(wp), POINTER     ::  hfl(:,:,:)    !< heights of full layers (cell centres) above sea level in COSMO-DE, computed as arithmetic average of hhl
    REAL(wp), POINTER     ::  depths(:)     !< depths of output soil layers, equal the depths of the source model (e.g. COSMO-DE)
    REAL(wp), ALLOCATABLE ::  xu(:)         !< coordinates of cell faces in x direction [m]
    REAL(wp), ALLOCATABLE ::  yv(:)         !< coordinates of cell faces in y direction [m]
    REAL(wp), POINTER     ::  zw(:)         !< coordinates of cell faces in z direction [m]
    REAL(wp)              ::  zt            !< mesoscale model top height [m]
    REAL(wp), POINTER     ::  zs_palm(:,:,:)!< PALM terrain height ('z surface') [m]
    REAL(wp), ALLOCATABLE ::  zs(:,:,:)     !< interpolated mesoscale terrain height [m]
    REAL(wp), ALLOCATABLE ::  lat(:)        !< rotated-pole latitudes of scalars (cell centers) of the COSMO-DE grid [rad]
    REAL(wp), ALLOCATABLE ::  lon(:)        !< rotated-pole longitudes of scalars (cell centres) of the COSMO-DE grid [rad]
    REAL(wp), ALLOCATABLE ::  latv(:)       !< rotated-pole latitudes of v winds (face centres in latitudal/y direction) [rad]
    REAL(wp), ALLOCATABLE ::  lonu(:)       !< rotated-pole latitudes of u winds (face centres in longitudal/x direction) [rad]
    REAL(wp), ALLOCATABLE ::  clat(:,:)     !< latitudes of PALM-4U cell centres in COSMO-DE's rotated-pole grid [rad]
    REAL(wp), ALLOCATABLE ::  clon(:,:)     !< longitudes of PALM-4U scalars (cell centres) in COSMO-DE's rotated-pole grid [rad]
    REAL(wp), ALLOCATABLE ::  clatu(:,:)    !< latitudes of PALM-4U u winds (cell faces in u direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(wp), ALLOCATABLE ::  clonu(:,:)    !< longitudes of PALM-4U u winds (cell faces in u direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(wp), ALLOCATABLE ::  clatv(:,:)    !< latitudes of PALM-4U v winds (cell faces in v direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(wp), ALLOCATABLE ::  clonv(:,:)    !< longitudes of PALM-4U v winds (cell faces in v direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(wp), ALLOCATABLE ::  w_horiz(:,:,:)   !< weights for bilinear horizontal interpolation
    REAL(wp), ALLOCATABLE ::  w_verti(:,:,:,:) !< weights for linear vertical interpolation
    REAL(wp), ALLOCATABLE ::  w(:,:,:)      !< vertical interpolation weights, w(<source_column>, <PALM k level>, <neighbour index>) [-]
 END TYPE grid_definition


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Container for name and dimensions of the netCDF output file
!------------------------------------------------------------------------------!
 TYPE nc_file
    CHARACTER(LEN=PATH)   ::  name              !< file name
    INTEGER               ::  dimid_time        !< NetCDF IDs of the time dimension
    INTEGER               ::  dimids_scl(3)     !< NetCDF IDs of the grid dimensions for scalar points x, y, z 
    INTEGER               ::  dimids_vel(3)     !< NetCDF IDs of the grid dimensions for velocity points xu, yu, zu
    INTEGER               ::  dimids_soil(3)    !< NetCDF IDs of the grid dimensions for soil points x, y, depth
    INTEGER               ::  dimvarid_time     !< NetCDF IDs of the time variable
    INTEGER               ::  dimvarids_scl(3)  !< NetCDF IDs of the grid coordinates of scalars x, y, z 
    INTEGER               ::  dimvarids_vel(3)  !< NetCDF IDs of the grid coordinates of velocities xu, yu, zu. Note that velocities are located at mix of both coordinates, e.g. u(xu, y, z).
    INTEGER               ::  dimvarids_soil(3) !< NetCDF IDs of the grid coordinates for soil points x, y, depth 
    REAL(wp), POINTER     ::  time(:)           !< vector of output time steps
 END TYPE nc_file


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Metadata container for netCDF variables
!------------------------------------------------------------------------------!
#if defined ( __netcdf )
 TYPE nc_var
    INTEGER                               ::  varid     !< NetCDF ID of the variable
    INTEGER                               ::  input_id  !< ID of the correpsonding input variables, only valid for output variables
    INTEGER                               ::  ndim      !< number of NetCDF dimensions
    INTEGER(iwp)                          ::  nt        !< number of output time steps
    INTEGER                               ::  lod       !< NetCDF attribute indicating the PALM-4U level of detail
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) ::  dimids    !< NetCDF IDs of the dimensions
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) ::  dimvarids !< IDs of NetCDF dimension variables
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) ::  dimlen    !< length of NetCDF dimensions
    CHARACTER(LEN=NF90_MAX_NAME), DIMENSION(NF90_MAX_VAR_DIMS) ::  dimname !< names of NetCDF dimensions
    CHARACTER(LEN=SNAME)                  ::  name                      !< NetCDF short name of the variable
    CHARACTER(LEN=LNAME)                  ::  standard_name             !< NetCDF standard name of the variable
    CHARACTER(LEN=LNAME)                  ::  long_name                 !< NetCDF long name of the variable
    CHARACTER(LEN=LNAME)                  ::  source                    !< NetCDF attribute indicating the data source for the output
    CHARACTER(LEN=SNAME)                  ::  units                     !< NetCDF units of the variable
    CHARACTER(LEN=SNAME)                  ::  kind                      !< Kind of grid
    CHARACTER(LEN=SNAME)                  ::  task                      !< Processing task that generates this variable, e.g. 'interpolate_2d' or 'average profile'
    LOGICAL                               ::  to_be_processed = .FALSE. !< INIFOR flag indicating whether variable shall be processed
    LOGICAL                               ::  is_internal = .FALSE.     !< INIFOR flag indicating whether variable shall be written to netCDF file (.FALSE.) or kept for later (.TRUE.)
    LOGICAL                               ::  is_optional = .FALSE.     !< Flag indicating whether INIFOR may continue if the the netCDF variable cannot be processed, e.g. if files are missing
    LOGICAL                               ::  is_read = .FALSE.         !< INIFOR flag indicating whether variable has been read
    LOGICAL                               ::  is_upside_down  = .FALSE. !< INIFOR flag indicating whether vertical dimension is reversed (typically the case with COSMO-DE atmospheric fields)
    LOGICAL                               ::  has_redundant_first_level !< INIFOR flag inidicating whether a soil variable has a redundant first level (e.g. COSMO's T_SO may contain the surface temperature at depth=0, which is a redundant copy the first model layer)
    TYPE(grid_definition), POINTER        ::  grid                      !< Pointer to the corresponding output grid
    TYPE(grid_definition), POINTER        ::  intermediate_grid         !< Pointer to the corresponding intermediate grid
    TYPE(grid_definition), POINTER        ::  averaging_grid            !< Pointer to the corresponding intermediate grid
 END TYPE nc_var


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Input/Output group, groups together nc_var-type output variabels that share 
!> input variables as well as lists of the netCDF files they are stored in.
!> For instance, all boundary surfaces and initialization fields of the
!> potential temperature are base on the input netCDF variables T and P.
!------------------------------------------------------------------------------!
 TYPE io_group
    INTEGER(iwp)                     ::  nt                  !< maximum number of output time steps across all output variables
    INTEGER(iwp)                     ::  nv                  !< number of netCDF output variables
    INTEGER(iwp)                     ::  n_inputs            !< number of input variables
    INTEGER(iwp)                     ::  n_output_quantities !< number of physical quantities required for computing netCDF output variables
    CHARACTER(LEN=SNAME)             ::  name                !< name of I/O group
    CHARACTER(LEN=SNAME)             ::  kind                !< kind of I/O group
    CHARACTER(LEN=PATH), ALLOCATABLE ::  in_files(:)         !< list of nt input files
    TYPE(nc_var), ALLOCATABLE        ::  out_vars(:)         !< list of output variables 
    TYPE(nc_var), ALLOCATABLE        ::  in_var_list(:)      !< list of input variables
    LOGICAL                          ::  to_be_processed = .FALSE. !< Inifor flag indicating whether I/O group shall be processed
    LOGICAL                          ::  is_accumulated = .FALSE.  !< Flag indicating whether this I/O group contains accumulated variables
    LOGICAL                          ::  is_optional = .FALSE.     !< Flag indicating whether INIFOR may continue if group cannot be processed, e.g. if files are missing
    LOGICAL                          ::  is_preprocessed = .FALSE. !< Inifor flag indicating whether the I/O group has been preprocessed
 END TYPE io_group 
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Container for input data arrays. read_input_variables() allocates a
!> one-dimensional array of containers, to accomodate all inputs of the given
!> IO group in one variable.
!------------------------------------------------------------------------------!
 TYPE container
   REAL(wp), ALLOCATABLE ::  array(:,:,:)               !< generic data array
   LOGICAL               ::  is_preprocessed = .FALSE.  !< flag indicating whether input array has been preprocessed
 END TYPE container

 END MODULE inifor_types
