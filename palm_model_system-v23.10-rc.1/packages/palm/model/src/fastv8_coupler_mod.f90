!> @file fastv8_coupler_mod.f90
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
! Copyright 2017-2023 Carl von Ossietzky Universität Oldenburg
! Copyright 2022-2023 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Coupler to multiple FASTv8 wind turbine Simulations using TCP/IP sockets
!--------------------------------------------------------------------------------------------------!
 MODULE fastv8_coupler_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               zu

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  pi

    USE control_parameters,                                                                        &
        ONLY:  current_timestep_number,                                                            &
               debug_output,                                                                       &
               dt_3d,                                                                              &
               dz,                                                                                 &
               end_time,                                                                           &
               initializing_actions,                                                               &
               interpolate_to_grid_center,                                                         &
               length,                                                                             &
               message_string,                                                                     &
               restart_data_format_output,                                                         &
               restart_string,                                                                     &
               simulated_time,                                                                     &
               terminate_run

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               dx,                                                                                 &
               ddy,                                                                                &
               dy


    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  comm2d,                                                                             &
               ierr,                                                                               &
               myid


    IMPLICIT NONE


    LOGICAL ::  fastv8_coupler_enabled = .FALSE.  !<
    LOGICAL ::  first_ts_with_wt = .TRUE.         !<

!
!-- Number of turbines
    INTEGER(iwp), PARAMETER ::  nturbines_max = 300_iwp  !< maximum number of turbines allowed (can be increased if required)
    INTEGER(iwp) ::  nturbines = 0_iwp                   !< namelist parameter
!
!-- Number of time steps in FAST
    INTEGER(iwp) ::  n_dt_fast = 0.0_wp  !<
!
!-- Time step size in FAST and PALM derived from the FAST dt
    REAL(wp) ::  dt_fast       !<
    REAL(wp) ::  dt_palm       !<
    REAL(wp) ::  current_time_fast_min  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  current_time_fast  !<
!
!-- Number of blades (1 value per wind turbine, value is received from FAST)
    INTEGER(iwp), DIMENSION(1:nturbines_max) ::  fast_n_blades = 0_iwp  !<

!
!-- Maximum allowed number of blades per turbine
    INTEGER(iwp) ::  fast_n_blades_max = 0_iwp  !< namelist parameter

!
!-- Number of blade elements per blade in FAST (1 value per wind turbine, value is received from FAST)
    INTEGER(iwp), DIMENSION(1:nturbines_max) ::  fast_n_blade_elem = 0_iwp  !<

!
!-- Maximum allowed number of blade elements per blade
    INTEGER(iwp) ::  fast_n_blade_elem_max = 0_iwp  !< namelist parameter

!
!-- Rotational velocity
    INTEGER(iwp) ::  n_sector_max  !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  n_sector  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rotspeed          !< rotational velocity from FAST
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  sector_angle      !< angle (in radiance) of sector, the rotor coveres in a PALM timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  sector_angle_deg  !< angle (in degree) of sector, the rotor coveres in a PALM timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  alpha_rot         !< angle (in radiance) of sector, the rotor coveres in a FAST timestep

    REAL(wp), DIMENSION(:),     ALLOCATABLE ::  shaft_height_fast  !< shaft height from FAST
    REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::  shaft_coordinates  !< shaft coordinates from FAST
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  r_n                !< rotation matrx

!
!-- Radius of turbines [ToDo]
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_n_radius = 0.0_wp  !<

    CHARACTER (LEN=15)                             ::  fast_host_addr_default = '127.0.0.1'  !<
    CHARACTER (LEN=15), DIMENSION(1:nturbines_max) ::  fast_host_addr = REPEAT(' ', 15)      !<
    CHARACTER (LEN=5),  DIMENSION(1:nturbines_max) ::  fast_host_port = REPEAT(' ', 5)       !<

!
!-- Tower base reference position in FAST (three components)
!-- Dimension has to be larger than / equal to nturbines
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_tower_ref_pos_x = 0.0_wp  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_tower_ref_pos_y = 0.0_wp  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_tower_ref_pos_z = 0.0_wp  !<

!
!-- Tower base reference position in PALM (three components, namelist parameter)
!-- Dimension has to be larger than / equal to nturbines
    REAL(wp), DIMENSION(1:nturbines_max) ::  palm_tower_ref_pos_x = - HUGE(1.0_wp)  !< namelist parameter
    REAL(wp), DIMENSION(1:nturbines_max) ::  palm_tower_ref_pos_y = - HUGE(1.0_wp)  !< namelist parameter
    REAL(wp), DIMENSION(1:nturbines_max) ::  palm_tower_ref_pos_z = 0.0_wp          !< namelist parameter

!
!-- Sum of tower base reference position in PALM and FAST
!-- Dimension has to be larger than / equal to nturbines
    REAL(wp), DIMENSION(1:nturbines_max) ::  tower_ref_pos_x  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  tower_ref_pos_y  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  tower_ref_pos_z  !<

!
!-- Hub center position (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fast_hub_center_pos  !<

!
!-- New 081012:
!-- In the calculations in PALM information on the positions
!-- is required at two different points in time
!-- Therefore, a second array for the positions of the hub center
!-- has to be added
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fast_hub_center_pos_old  !<

!
!-- Spinner aerodynamic force (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3
!    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fast_spinner_force

!
!-- Hub center velocity (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  palm_hub_center_vel  !<

    !-- Position of rotor blade elements (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions
!-- to be allocated nturbines, fast_n_blades, fast_n_blade_elem, 3
!-- Note: fast_n_blades has to be set to the maximum number of blades
!-- for all simulated wind turbines; the missing value has to be
!-- entered for turbines with less blades than the turbine with the
!-- maximum number of blades
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_position_blade  !<

!
!-- In the calculations in PALM information on the positions
!-- is required at two different points in time
!-- Therefore, a second array for the positions of the blades
!-- has to be added
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_position_blade_old  !<

!
!-- Aerodynamic forces on blade segments (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions
!-- to be allocated nturbines, fast_n_blades, fast_n_blade_elem, 3
!-- Note: fast_n_blades has to be set to the maximum number of blades
!-- for all simulated wind turbines; the missing value has to be
!-- entered for turbines with less blades than the turbine with the
!-- maximum number of blades
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_force_blade      !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_force_blade_old  !<

!
!-- Blade velocities (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions
!-- to be allocated nturbines, fast_n_blades, fast_n_blade_elem, 3
!-- Note: fast_n_blades has to be set to the maximum number of blades
!-- for all simulated wind turbines; the missing value has to be
!-- entered for turbines with less blades than the turbine with the
!-- maximum number of blades
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  palm_vel_blade  !<

!
!-- Time at which the wind turbine is switched on in the
!-- simulation (by default it would be switched on at the beginning
!-- of the run, namelist parameter)
    REAL(wp) ::  time_turbine_on = 0.0_wp  !<

!
!-- Simulated time in FAST (about simulated_time-time_turbine_on)
    !REAL(wp)    :: fast_simulated_time

!
!-- 3D arrays of forces caused by wind turbines in the model domain
!-- (forces smeared with a regularization kernel)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  force_x  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  force_y  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  force_z  !<

!
!-- Regularization factor for force distribution (epsilon = reg_fac * dx)
    REAL(wp) ::  reg_fac = 2.0_wp  !< namelist parameter

!
!-- Force distribution cut-off factor
!-- Edge length of box around turbine in which forces are distributed (edge length = edgelen_fac * turbine radius)
!-- Requirements for the size of the box:
!-- Around every point of attack there should be the at least the distance of 3 * grid spacing * reg_fac
!-- for force distribution available
    REAL(wp) ::  fbox_fac = 1.5_wp  !< namelist parameter

!-- Force distribution boxes around turbines
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3, 2
!-- number of turbines; x_min, x_max; y_min, y_max; z_min, z_max
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  fboxcorners  !<

!
!-- Field with values of the exponential function
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  expvalue  !<

!
!-- x-, y- and z-components of the turbine induced forces:
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  thrx  !< thrx_av
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tory  !< tory_av
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  torz  !< torz_av


!
! -- Parameters and Variables for the calculation of tower area and thrust applied by the tower
    REAL(wp) ::  thrust_tower_x  !<
    REAL(wp) ::  thrust_tower_y  !<

    REAL(wp), DIMENSION(1:nturbines_max) ::  dtow         = 0.0_wp  !< tower diameter [m]
    REAL(wp), DIMENSION(1:nturbines_max) ::  htow         = 0.0_wp  !< tower height [m]
    REAL(wp), DIMENSION(1:nturbines_max) ::  turb_C_d_tow = 1.2_wp  !< drag coefficient for tower

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tower_area_x  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tower_area_y  !<

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tower_area_x_4d  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tower_area_y_4d  !<

    SAVE

    PRIVATE
!
!-- Public functions
    PUBLIC                                                                                         &
       f8c_actions,                                                                                &
       f8c_check_data_output,                                                                      &
       f8c_check_parameters,                                                                       &
       f8c_data_output_2d,                                                                         &
       f8c_data_output_3d,                                                                         &
       f8c_define_netcdf_grid,                                                                     &
       f8c_header,                                                                                 &
       f8c_init,                                                                                   &
       f8c_init_arrays,                                                                            &
       f8c_parin
!
!-- Public parameters, constants and initial values
   PUBLIC                                                                                          &
      current_time_fast,                                                                           &
      dt_fast,                                                                                     &
      fast_force_blade,                                                                            &
      fast_position_blade,                                                                         &
      fast_n_blades,                                                                               &
      fast_n_blade_elem,                                                                           &
      fast_n_radius,                                                                               &
      fast_hub_center_pos,                                                                         &
      nturbines,                                                                                   &
      palm_hub_center_vel,                                                                         &
      palm_vel_blade,                                                                              &
      rotspeed,                                                                                    &
      shaft_height_fast,                                                                           &
      fastv8_coupler_enabled,                                                                      &
      shaft_coordinates


    INTERFACE f8c_parin
       MODULE PROCEDURE f8c_parin
    END INTERFACE f8c_parin

    INTERFACE f8c_check_parameters
       MODULE PROCEDURE f8c_check_parameters
    END INTERFACE f8c_check_parameters

    INTERFACE f8c_check_data_output
       MODULE PROCEDURE f8c_check_data_output
    END INTERFACE f8c_check_data_output

    INTERFACE f8c_define_netcdf_grid
       MODULE PROCEDURE f8c_define_netcdf_grid
    END INTERFACE f8c_define_netcdf_grid

    INTERFACE f8c_init
       MODULE PROCEDURE f8c_init
    END INTERFACE f8c_init

    INTERFACE f8c_init_arrays
       MODULE PROCEDURE f8c_init_arrays
    END INTERFACE f8c_init_arrays

    INTERFACE f8c_header
       MODULE PROCEDURE f8c_header
    END INTERFACE f8c_header

    INTERFACE f8c_actions
       MODULE PROCEDURE f8c_actions
       MODULE PROCEDURE f8c_actions_ij
    END INTERFACE f8c_actions

    INTERFACE f8c_data_output_2d
       MODULE PROCEDURE f8c_data_output_2d
    END INTERFACE f8c_data_output_2d

    INTERFACE f8c_data_output_3d
       MODULE PROCEDURE f8c_data_output_3d
    END INTERFACE f8c_data_output_3d


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &fastv8_coupler_parameters for user module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_parin

    CHARACTER (LEN=80) ::  line  !< string containing the last line read from namelist file

    INTEGER(iwp) ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /fastv8_coupler_parameters/                                                           &
       switch_off_module,                                                                          &
       nturbines,                                                                                  &
       time_turbine_on,                                                                            &
       fast_n_blades_max,                                                                          &
       fast_n_blade_elem_max,                                                                      &
       fast_host_addr_default,                                                                     &
       fast_host_addr,                                                                             &
       fast_host_port,                                                                             &
       palm_tower_ref_pos_x,                                                                       &
       palm_tower_ref_pos_y,                                                                       &
       palm_tower_ref_pos_z,                                                                       &
       reg_fac,                                                                                    &
       fbox_fac,                                                                                   &
       dtow,                                                                                       &
       htow,                                                                                       &
       turb_C_d_tow


!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND ( 11 )
    READ( 11, fastv8_coupler_parameters, IOSTAT=io_status )

!
!-- Actions depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    fastv8_coupler_parameters namelist was found and read correctly. Set flag that
!--    fastv8_coupler_mod is switched on.
       IF ( .NOT. switch_off_module )  fastv8_coupler_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    User namelist was found, but contained errors. Print an error message containing the line
!--    that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'fastv8_coupler_parameters', line )

    ENDIF

 END SUBROUTINE f8c_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check &fastv8_coupler_parameters control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_check_parameters

    INTEGER(iwp) ::  i  !< loop variable

#ifndef __fastv8
    IF ( fastv8_coupler_enabled )  THEN
       WRITE( message_string, * )  'Using FASTv8 Coupler requires cpp macro "__fastv8" during build'
       CALL message( 'f8c_check_parameters', 'F8C0001', 1, 2, 0, 6, 0 )
    ENDIF
#endif

    IF ( nturbines <= 0_iwp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"nturbines" must set to a value greater than 0 ',                                       &
          '(currently it is set to ', nturbines, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0002', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( nturbines > nturbines_max )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter "nturbines" ',    &
          'must be set to a value less or equal to ', nturbines_max,                               &
          '(currently it is set to ', nturbines, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0003', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( fast_n_blades_max <= 0_iwp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"fast_n_blades_max" must set to a value greater than 0 ',                               &
          '(currently it is set to ', fast_n_blades_max, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0004', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( fast_n_blade_elem_max <= 0_iwp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"fast_n_blade_elem_max" must set to a value greater than 0 ',                           &
          '(currently it is set to ', fast_n_blade_elem_max, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0005', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM(fast_host_addr_default) == '' )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"fast_host_addr_default" must be set to a valid IP adress e.g 127.0.0.1',               &
          '(currently it is not set)'
       CALL message( 'f8c_check_parameters', 'F8C0006', 1, 2, 0, 6, 0 )
    ENDIF

    DO i = 1, nturbines
       IF ( TRIM(fast_host_port(i)) == '' )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"fast_host_port(', i, ')" must be set to a valid port number ',                      &
             '(currently it is not set)'
          CALL message( 'f8c_check_parameters', 'F8C0007', 1, 2, 0, 6, 0 )
       ENDIF
    END DO

    DO i = 1, nturbines
       IF ( palm_tower_ref_pos_x(i) < 0.0_wp  .OR.  palm_tower_ref_pos_x(i) > (nx+1) * dx )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"palm_tower_ref_pos_x(', i, ')" must be within the model domain. ',                  &
             '(currently it is set to ', palm_tower_ref_pos_x(i), ' )'
          CALL message( 'f8c_check_parameters', 'F8C0008', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( palm_tower_ref_pos_y(i) < 0.0_wp  .OR.  palm_tower_ref_pos_y(i) > (ny+1) * dy )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"palm_tower_ref_pos_y(', i, ')" must be within the model domain.',                   &
             '(currently it is set to ', palm_tower_ref_pos_y(i), ' )'
          CALL message( 'f8c_check_parameters', 'F8C0009', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( palm_tower_ref_pos_z(i) < 0.0_wp  .OR.  palm_tower_ref_pos_z(i) > zu(nzt) )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"palm_tower_ref_pos_z(', i, ')" must be within the model domain.',                   &
             '(currently it is set to ', palm_tower_ref_pos_z(i), ' )'
          CALL message( 'f8c_check_parameters', 'F8C0010', 1, 2, 0, 6, 0 )
       ENDIF
    END DO

    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
      WRITE( message_string, * )  'Currently initializing_actions = "read_restart_data" is ',      &
          'not implemented while using the FASTv8 Coupler'
      CALL message( 'f8c_check_parameters', 'F8C0011', 1, 2, 0, 6, 0 )
   ENDIF


 END SUBROUTINE f8c_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined output quantities. For those variables not recognized by the user,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!> is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_check_data_output( variable, unit )

    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<


    SELECT CASE ( TRIM( variable ) )

       CASE ( 'thrx' )
          unit = 'm/s2'

       CASE ( 'tory' )
          unit = 'm/s2'

       CASE ( 'torz' )
          unit = 'm/s2'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE f8c_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize user-defined arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_init_arrays

!
!-- Allocation of the arrays required for the output of additional 2D and 3D data.
    ALLOCATE( thrx(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tory(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( torz(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Allocation of position, force and velocity arrays that are
!-- required for the exchange of information between PALM and FAST.
!-- Arrays with dimensions dependent on the number of turbines:
!-- current_time_fast: first dimension number of turbines
    ALLOCATE( current_time_fast(1:nturbines) )
!
!-- n_sector: first dimension number of turbines
    ALLOCATE( n_sector(1:nturbines) )
!
!-- rotspeed: first dimension number of turbines
    ALLOCATE( rotspeed(1:nturbines) )
!
!-- sector_angle: first dimension number of turbines
    ALLOCATE( sector_angle(1:nturbines) )
!
!-- sector_angle_deg: first dimension number of turbines
    ALLOCATE( sector_angle_deg(1:nturbines) )
!
!-- alpha_rot: first dimension number of turbines
    ALLOCATE( alpha_rot(1:nturbines) )
!
!-- shaft_height_fast: first dimension number of turbines
    ALLOCATE( shaft_height_fast(1:nturbines) )
!
!-- shaft_coordinates: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( shaft_coordinates(1:nturbines,1:3) )
!
!-- r_n: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( r_n(1:nturbines,1:3,1:3) )
!
!-- fast_hub_center_pos: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( fast_hub_center_pos(1:nturbines,1:3) )
!
!-- fast_hub_center_pos_old: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( fast_hub_center_pos_old(1:nturbines,1:3) )
!
!-- palm_hub_center_vel: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( palm_hub_center_vel(1:nturbines,1:3) )
!
!-- Arrays with dimension dependent on the number of blades
!-- and the number of blade elements per blade
!
!-- fast_position_blade: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( fast_position_blade(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- New: 081012
!-- fast_position_blade_old: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( fast_position_blade_old(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- fast_force_blade: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( fast_force_blade(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
    ALLOCATE( fast_force_blade_old(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- palm_vel_blade: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( palm_vel_blade(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- Allocation of the arrays of forces (one in x-, y- and z-direction)
    ALLOCATE( force_x(nzb:nzt,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( force_y(nzb:nzt,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( force_z(nzb:nzt,nysg:nyng,nxlg:nxrg) )

!-- fboxcorners: first dimension number of turbines,
!-- second dimension three direction in space, third dimension min and max value
    ALLOCATE( fboxcorners(1:nturbines,1:3,1:2) )

!
!-- Allocate the field with values of the exponential function
    ALLOCATE( expvalue(0:10000000) )

!
!-- First of all the arrays that contain the information on the overlapping
!-- areas for each LES grid box have to be allocated; afterwards they are set to
!-- a basic value of 0 (no overlap between the rotor area and the LES grid box):
    ALLOCATE( tower_area_x(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tower_area_y(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!-- Additionally 4D-arrays with the number of the wind turbines as fourth
!-- dimension are required:
    ALLOCATE( tower_area_x_4d(1:nturbines,nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tower_area_y_4d(1:nturbines,nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

 END SUBROUTINE f8c_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined initializing actions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_init

    USE ISO_C_BINDING,                                                                             &
       ONLY: C_CHAR,                                                                               &
             C_INT,                                                                                &
             C_LOC,                                                                                &
             C_NULL_CHAR,                                                                          &
             C_NULL_PTR,                                                                           &
             C_PTR


    INTEGER(iwp) ::  i  !< loop variable
    INTEGER(iwp) ::  j  !< loop variable
    INTEGER(iwp) ::  k  !< loop variable
    INTEGER(iwp) ::  o  !< loop variable

    INTEGER(iwp) ::  iflag    !<
    INTEGER(iwp) ::  endflag  !<

    INTEGER(iwp) ::  i_hub  !<
    INTEGER(iwp) ::  j_hub  !<
    INTEGER(iwp) ::  k_hub  !<

    INTEGER(iwp) ::  tower_l  !<
    INTEGER(iwp) ::  tower_r  !<
    INTEGER(iwp) ::  tower_n  !<
    INTEGER(iwp) ::  tower_s  !<

    REAL(wp) ::  rcx  !<
    REAL(wp) ::  rcy  !<
    REAL(wp) ::  rcz  !<

    TYPE(C_PTR) ::  fast_host_addr_c(SIZE(fast_host_addr) + 1)  !<
    TYPE(C_PTR) ::  fast_host_port_c(SIZE(fast_host_port) + 1)  !<

    TYPE string
      CHARACTER(LEN=:,KIND=C_CHAR), ALLOCATABLE ::  item  !<
    END TYPE string

    TYPE(string), TARGET :: tmp_addr(SIZE(fast_host_addr))  !<
    TYPE(string), TARGET :: tmp_port(SIZE(fast_host_port))  !<

#if defined( __fastv8 )
    INTERFACE

       FUNCTION init_comm( ifnumturb, host_addr, host_port ) BIND(C,NAME='init_comm' )
          USE ISO_C_BINDING,                                                                       &
              ONLY:  C_INT,                                                                        &
                     C_PTR
          INTEGER(kind=C_INT) ::  init_comm
          INTEGER(kind=C_INT), VALUE ::  ifnumturb
          TYPE(C_PTR), INTENT(IN) ::  host_addr(*)
          TYPE(C_PTR), INTENT(IN) ::  host_port(*)
       END FUNCTION init_comm

    END INTERFACE
#endif

!
!-- Compute values of the exponential function.
    expvalue(:) = 0.0_wp
    DO  i = 0, 10000000
       expvalue(i) = EXP(-0.001_wp*i)
       IF ( expvalue(i) < 0.00000001 )  THEN
          EXIT
       ENDIF
    ENDDO

    endflag = 0

!
!-- Only do if turbines are supposed to be simulated.
    IF ( nturbines > 0 )  THEN

       IF ( myid == 0 )  THEN

          DO  i = 1, nturbines
             IF ( TRIM( fast_host_addr(i) ) == '' )  THEN
               fast_host_addr(i) = fast_host_addr_default
             ENDIF
          ENDDO

          DO  i = 1, SIZE( fast_host_addr )
!
!--          This may involve kind conversion.
             tmp_addr(i)%item = TRIM( fast_host_addr(i) ) // C_NULL_CHAR
             fast_host_addr_c(i) = C_LOC( tmp_addr(i)%item )
          ENDDO
          fast_host_addr_c(SIZE( fast_host_addr )) = C_NULL_PTR

          DO  i = 1, SIZE( fast_host_port )
!
!--          This may involve kind conversion.
             tmp_port(i)%item = TRIM( fast_host_port(i) ) // C_NULL_CHAR
             fast_host_port_c(i) = C_LOC( tmp_port(i)%item )
          END DO
          fast_host_port_c(SIZE( fast_host_port )) = C_NULL_PTR

!
!--       Initialisieren der Verbindung mit FAST
#if defined( __fastv8 )
          iflag = init_comm( nturbines, fast_host_addr_c, fast_host_port_c )
#else
          iflag = -1
#endif
          IF ( iflag /= 0 )  THEN
             WRITE( message_string, * )  'Unable to initalize communication with FAST server(s)'
             CALL message( 'f8c_actions', 'F8C0012', 2, 2, 0, 6, 0 )
          ENDIF

       ENDIF

    ENDIF !(nturbines > 0)

    current_time_fast = 0.0_wp

!
!-- Determine the area within each grid cell that overlaps with the area of the nacelle and the
!-- tower (needed for calculation of the forces).
!-- Note: so far this is only a 2D version, in that the mean flow is
!-- perpendicular to the rotor area.
    tower_area_x_4d(:,:,:,:) = 0.0_wp
    tower_area_y_4d(:,:,:,:) = 0.0_wp

    IF ( myid == 0  .AND.  debug_output )  THEN
       WRITE(9,*) "f8c_init - before determination of tower_area"
       FLUSH( 9 )
    ENDIF

!
!-- Loop over all turbines.
    DO  o = 1, nturbines

       rcx = palm_tower_ref_pos_x(o)
       rcy = palm_tower_ref_pos_y(o)
       rcz = palm_tower_ref_pos_z(o) + htow(o)

       tower_area_x(:,:,:) = 0.0_wp
       tower_area_y(:,:,:) = 0.0_wp

!
!--    Determine the indices of the hub height.
       i_hub = INT( rcx / dx )
       j_hub = INT( ( rcy + 0.5_wp * dy ) / dy )
       k_hub = INT( ( rcz + 0.5_wp * dz(1) ) / dz(1) )

!
!--    Determine the indices of the grid boxes containing the left and
!--    the right boundaries of the tower.
       tower_n = INT( ( rcy + 0.5_wp * dtow(o) - 0.5_wp * dy ) / dy )
       tower_s = INT( ( rcy - 0.5_wp * dtow(o) - 0.5_wp * dy ) / dy )
       tower_r = INT( ( rcx + 0.5_wp * dtow(o) - 0.5_wp * dx ) / dx )
       tower_l = INT( ( rcx - 0.5_wp * dtow(o) - 0.5_wp * dx ) / dx )

       IF ( myid == 0  .AND.  debug_output )  THEN
          WRITE(9,*) 'id: ', i, ' tower_n ', tower_n
          WRITE(9,*) 'id: ', i, ' tower_s ', tower_s
          WRITE(9,*) 'id: ', i, ' tower_r ', tower_r
          WRITE(9,*) 'id: ', i, ' tower_l ', tower_l
          FLUSH( 9 )
       ENDIF
!
!--    Determine the fraction of the grid box area overlapping with the tower area.
       IF ( ( nxlg <= i_hub ) .AND. ( nxrg >= i_hub ) .AND.                                        &
            ( nysg <= j_hub ) .AND. ( nyng >= j_hub ) )                                            &
       THEN
!
!--       Loop from the southernmost grid index of the tower to the northernmost.
          DO  k = nzb, k_hub
!
!--          Loop from the southernmost grid index of the tower to the northernmost.
             DO  j = tower_s, tower_n
!
!--             If tower not completely inside one grid box.
                IF ( tower_n - tower_s >= 1 )  THEN

                   IF ( j == tower_s )  THEN
                      tower_area_x(k,j,i_hub) =                                                    &
                                                ! extension in z-direction
                                                MIN(  rcz - ( k * dz(1) - 0.5_wp * dz(1) ),        &
                                                      dz(1) ) *                                    &
                                                ! extension in y-direction
                                                ( ( tower_s + 1 + 0.5_wp ) * dy                    &
                                                - ( rcy - 0.5_wp * dtow(o) ) )
                   ELSEIF ( j == tower_n )  THEN
                      tower_area_x(k,j,i_hub) =                                                    &
                                                ! extension in z-direction
                                                MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ) ,        &
                                                     dz(1) ) *                                     &
                                                ! extension in y-direction
                                                ( rcy + 0.5_wp * dtow(o)                           &
                                                - ( tower_n + 0.5_wp ) * dy )
                   ELSE
!
!--                   Grid boxes inbetween (where tower_area = grid box area):
                      tower_area_x(k,j,i_hub) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) * dy
                   ENDIF
                ELSE
!
!--                Tower lies completely within one grid box.
                   tower_area_x(k,j,i_hub) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),            &
                                                  dz(1) ) * dtow(o)
                ENDIF

             ENDDO

!
!--          Loop from the left grid index of the tower to the right index
             DO  i = tower_l, tower_r
!
!--             If tower not completely inside one grid box
                IF ( tower_l - tower_r >= 1 )  THEN
!--                leftmost and rightmost grid box:
                   IF ( i == tower_l )  THEN

                      tower_area_y(k,j_hub,i) =                                                    &
                                                ! extension in z-direction
                                                MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) *                                     &
                                                ! extension in y-direction
                                                ( ( tower_l + 1 ) * dx                             &
                                                - ( rcx - 0.5_wp * dtow(o) ) )
                   ELSEIF ( i == tower_r )  THEN
                      tower_area_y(k,j_hub,i) =                                                    &
                                                ! extension in z-direction
                                                MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) *                                     &
                                                ! extension in y-direction
                                                ( rcx + 0.5 * dtow(o) - tower_r * dx )
                   ELSE
!
!--                   Grid boxes inbetween (where tower_area = grid box area).
                      tower_area_y(k,j_hub,i) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) * dy
                   ENDIF
                ELSE
!
!--                Tower lies completely within one grid box.
                   tower_area_y(k,j_hub,i) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),            &
                                                  dz(1) ) * dtow(o)
                ENDIF
             ENDDO
          ENDDO
       ENDIF !( ( nxlg <= i_hub ) .AND. ( nxrg >= i_hub ) .AND. ( nysg <= j_hub ) .AND. ( nyng >= j_hub ) )

       CALL exchange_horiz( tower_area_x, nbgp )
       CALL exchange_horiz( tower_area_y, nbgp )
!
!--    Restore tower area field in 4D-field containing each turbine.
       tower_area_x_4d(o,:,:,:) = tower_area_x(:,:,:)
       tower_area_y_4d(o,:,:,:) = tower_area_y(:,:,:)
!
!--    Tabulate the points on the circle that are required in the following for
!--    the calculation of the Riemann integral (node points; they are called
!--    circle_points in the following).

    ENDDO  ! end of loop over turbines

!
!-- Initializing 3D force arrays
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
         DO  k = nzb, nzt
             force_x(k,j,i) = 0.0_wp
             force_y(k,j,i) = 0.0_wp
             force_z(k,j,i) = 0.0_wp
          ENDDO
       ENDDO
   ENDDO

 END SUBROUTINE f8c_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which output quantities are defined.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    CHARACTER(LEN=*), INTENT(OUT) ::  grid_x  !< x grid of output variable
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_y  !< y grid of output variable
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_z  !< z grid of output variable
    CHARACTER(LEN=*), INTENT(IN)  ::  var     !< name of output variable

    LOGICAL, INTENT(OUT) ::  found   !< flag if output variable is found

    found  = .TRUE.
!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'thrx', 'thrx_xy', 'thrx_xz', 'thrx_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       u grid
          ELSE
             grid_x = 'xu'
             grid_y = 'y'
             grid_z = 'zu'
          ENDIF

       CASE ( 'tory', 'tory_xy', 'tory_xz', 'tory_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       v grid
          ELSE
             grid_x = 'x'
             grid_y = 'yv'
             grid_z = 'zu'
          ENDIF

       CASE ( 'torz', 'torz_xy', 'torz_xz', 'torz_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       w grid
          ELSE
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zw'
          ENDIF

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE f8c_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print a header with user-defined information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_header( io )

    INTEGER(iwp) ::  io  !<


    WRITE ( io, 1 )
    WRITE ( io, 20 )
    WRITE ( io, 21 ) nturbines

!
!-- Format-descriptors
1   FORMAT ( //' FASTv8 coupler information:'/ ' ------------------------------------------'/ )

20  FORMAT ( '--> Essential parameters:' )
21  FORMAT ( '       number of active turnines    :   nturbines   = ', I3)

 END SUBROUTINE f8c_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_actions( location )

    CHARACTER(LEN=*) ::  location  !<
!
!-- loop variables
    INTEGER(iwp) ::  i    !<
    INTEGER(iwp) ::  j    !<
    INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  l    !<
    INTEGER(iwp) ::  o    !<
    INTEGER(iwp) ::  m    !<
    INTEGER(iwp) ::  n    !<
    INTEGER(iwp) ::  n_s  !<
!
!-- grid point numbers
    INTEGER(iwp) ::  ii  !<
    INTEGER(iwp) ::  jj  !<
    INTEGER(iwp) ::  kk  !<
!
!-- run control
    INTEGER(iwp) ::  iflag  !<

    INTEGER(iwp) ::  expargument  !<

#if defined( __parallel )
    INTEGER(iwp) ::  size_array_type_1  !<
#endif
!
!-- Velocity interpolation and force distribution
    REAL(wp) ::  x  !<
    REAL(wp) ::  y  !<
    REAL(wp) ::  z  !<

    REAL(wp) ::  aa  !<
    REAL(wp) ::  bb  !<
    REAL(wp) ::  cc  !<
    REAL(wp) ::  dd  !<
    REAL(wp) ::  gg  !<

    REAL(wp) ::  u_int_l  !<
    REAL(wp) ::  u_int_u  !<
    REAL(wp) ::  v_int_l  !<
    REAL(wp) ::  v_int_u  !<
    REAL(wp) ::  w_int_l  !<
    REAL(wp) ::  w_int_u  !<

    REAL(wp) ::  eps_kernel    !<
    REAL(wp) ::  eps_kernel_2  !<
    REAL(wp) ::  nenner        !<
    REAL(wp) ::  kernel_value  !<
    REAL(wp) ::  distance_hub  !<
    REAL(wp) ::  fboxelen      !<

    REAL(wp), DIMENSION(3) ::  tmp_res  !<

    REAL(wp), DIMENSION(:,:),     ALLOCATABLE ::  u_int_h  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  u_int_b  !<

    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  x_prime_sec  !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  y_prime_sec  !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  z_prime_sec  !<

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  x_prime  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  y_prime  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  z_prime  !<

    LOGICAL ::  palm_dt_3d = .FALSE.

#if defined( __fastv8 )
    INTERFACE

       FUNCTION commclient(target, message) BIND(C,NAME='commclient')
          use ISO_C_BINDING
          INTEGER(kind=C_INT) ::  commclient
          INTEGER(kind=C_INT), value ::  target
          INTEGER(kind=C_INT), value ::  message
       END FUNCTION commclient

       FUNCTION time_data(simtime, simdt, simstep, simturbon, simendtime) BIND(C,NAME='time_data')
          use ISO_C_BINDING
          INTEGER(kind=C_INT) ::  time_data
          REAL(kind=c_double), value ::  simtime
          REAL(kind=c_double), value ::  simdt
          INTEGER(kind=C_INT), value ::  simstep
          REAL(kind=c_double), value ::  simturbon
          REAL(kind=c_double), value ::  simendtime
       END FUNCTION time_data

    END INTERFACE
#endif

    CALL cpu_log( log_point(24), 'f8c_actions', 'start' )

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( myid == 0 )  THEN

             IF ( debug_output )  THEN
                WRITE(9,*) 'starting new time step', simulated_time, current_timestep_number
                FLUSH(9)
             ENDIF
!
!--          update information between fortran and c code
#if defined( __fastv8 )
             IF ( simulated_time < time_turbine_on )  THEN
                iflag = time_data(simulated_time, dt_3d, current_timestep_number,                  &
                                  time_turbine_on, end_time)
             ELSE
!
!--             Hier wird die Zeit von PALM an die .c Datei geschickt
                iflag = time_data(simulated_time, dt_fast, current_timestep_number,                &
                                  time_turbine_on, end_time)
             ENDIF !( simulated_time < time_turbine_on )
#else
             iflag = -1
#endif

             IF ( iflag /= 0 )  THEN
                WRITE( message_string, * )  'Unable to update information between ',       &
                   'FORTRAN and C:', simulated_time, ':', current_timestep_number
                CALL message( 'f8c_actions', 'F8C0013', 2, 2, 0, 6, 0 )
             ENDIF !( iflag /= 0 )

          ENDIF !( myid == 0 )

!--       Check if connection to all FAST instances can be established and exchange some data for initialisation.
          IF ( myid == 0 )  THEN
             IF ( current_timestep_number == 0 )  THEN

               CALL location_message( 'Checking connection to FAST server(s)', 'start' )
#if defined( __fastv8 )
                iflag = commclient(0, 1)
#else
                iflag = -1
#endif

                IF ( iflag < 0 )  THEN
                   WRITE( message_string, * )                                                      &
                          'Unable to establish connection with all FAST server(s)'
                   CALL message( 'f8c_actions', 'F8C0014', 2, 2, 0, 6, 0 )
                ELSEIF ( iflag == 1 )  THEN
                   WRITE( message_string, * )  'At least one FAST server finished the simulation'
                   CALL message( 'f8c_actions', 'F8C0015', 2, 2, 0, 6, 0 )
                ELSEIF ( iflag == 0 )  THEN
                   CALL location_message(                                                          &
                           'Connection to FAST server(s) successfully established', '' )
                ENDIF !( iflag < 0 )
                CALL location_message( 'Checking connection to FAST server(s)',    'finished' )
             ENDIF !( current_timestep_number == 0 )

          ENDIF !( myid == 0 )

!--       Errechnen wie viele FAST sub-Zeitschritte es geben muss
          IF ( myid == 0 )  THEN
             dt_palm = dt_3d
             palm_dt_3d = .TRUE.
             n_dt_fast = FLOOR( dt_palm / dt_fast )
          ENDIF !( myid == 0 )

#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( dt_fast, 1, MPI_REAL, 0, comm2d, ierr )
          CALL MPI_BCAST( dt_palm, 1, MPI_REAL, 0, comm2d, ierr )
          CALL MPI_BCAST( fast_n_radius, nturbines, MPI_REAL, 0, comm2d, ierr )  ! [todo]
          CALL MPI_BARRIER( comm2d, ierr )
#endif

          IF ( simulated_time >= time_turbine_on )  THEN

             dt_3d = dt_palm
!
!--          Besondere Aktionen in Zusammenhang mit der Kopplung von
!--          PALM und FAST sind nur dann durchzufuehren, wenn die
!--          Windenergieanlagen eingeschaltet sind; der Zeitpunkt des
!--          Einschaltens der Anlagen ist ueber den Parameter
!--          time_turbine_on festgelegt, der in der
!--          Namelistparameterdatei unter userpar vorzugeben ist
!
!--          Beim ersten Zeitschritt ist der Ablauf etwas anders
!--          als an den anderen Zeitschritten. Zum ersten Zeitschritt,
!--          anders als zu den anderen Zeitschritten, erhaelt PALM noch
!--          keine Informationen ueber Kraefte aus FAST, sondern
!--          lediglich Informationen ueber die Positionen der
!--          Blattelemente
!
!--          first time step with wind turbine: Sending of signal to FAST to show that PALM is ready
!--          Receive blade element positions from FAST
!--          all other time steps with wind turbine: asking for positions
             IF ( myid == 0 )  THEN

                IF ( first_ts_with_wt )  THEN

#if defined( __fastv8 )
                   iflag = commclient(0, 2)
#else
                   iflag = -1
#endif

                   IF ( iflag < 0 )  THEN
                      WRITE( message_string, * )  'Unable to send start signale to FAST server(s)'
                      CALL message( 'f8c_actions', 'F8C0016', 2, 2, 0, 6, 0 )
                   ELSEIF ( iflag == 1 )  THEN
                      WRITE( message_string, * )  'At least one FAST server finished the simulation'
                      CALL message( 'f8c_actions', 'F8C0017', 2, 2, 0, 6, 0 )
                   ELSEIF ( iflag == 0 )  THEN
                      CALL location_message(                                                       &
                              'Start signal to FAST servers successfully broadcasted', '' )
                   ENDIF !( iflag < 0 )
!
!--                Determination of the corners of the boxes around every turbine for force distribution
                   IF ( iflag == 0 )  THEN

                      DO  i = 1, nturbines

                         fboxelen = fast_n_radius(i) * fbox_fac

                         fboxcorners(i, 1, 1) = FLOOR(                                             &
                            ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - fboxelen ) * ddx                          &
                                                     )

                         fboxcorners(i, 1, 2) = CEILING(                                           &
                            ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) + fboxelen ) * ddx                          &
                                                       )

                         fboxcorners(i, 2, 1) = FLOOR(                                             &
                            ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - fboxelen ) * ddy                          &
                                                     )

                         fboxcorners(i, 2, 2) = CEILING(                                           &
                            ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) + fboxelen ) * ddy                          &
                                                       )

                         fboxcorners(i, 3, 1) = FLOOR(                                             &
                            ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - fboxelen ) / dz(1)                        &
                                                     )

                         fboxcorners(i, 3, 2) = CEILING(                                           &
                            ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) + fboxelen ) / dz(1)                        &
                                                       )
                      ENDDO

                   ENDIF !( iflag == 0)

                   first_ts_with_wt = .FALSE.

                ENDIF !( first_ts_with_wt )

             ENDIF !( myid == 0 )

#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
!
!--          sending of fboxcorners
             CALL MPI_BCAST( fboxcorners(1,1,1), nturbines*3*2, MPI_REAL, 0, comm2d, ierr )
!
!--          Senden von fast_n_blades (nturbines) -> senden Anzahl der Blätter für Turbine 1
             CALL MPI_BCAST( fast_n_blades(1), nturbines, MPI_REAL, 0, comm2d, ierr )
!
!--          Senden von fast_n_blade_elem (nturbines) -> senden Anzahl der Blattelemente von Turbine 1
             CALL MPI_BCAST( fast_n_blade_elem(1), nturbines, MPI_REAL, 0, comm2d, ierr )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
!
!--          Counter for  number of FAST timesteps for a single PALM timestep
             l = 0
             IF ( myid == 0 )  THEN
                current_time_fast_min = MINVAL(current_time_fast(:))

                DO i = 1, nturbines
                   IF ( current_time_fast_min < current_time_fast(i) )  THEN
                      WRITE( message_string, * )  'Different times on different FAST instances: ', &
                         'i, current_time_fast_min, current_time_fast(i) ',                        &
                          i, current_time_fast_min, current_time_fast(i)
                      CALL message( 'f8c_actions', 'F8C0018', 2, 2, 0, 6, 0 )
                   ENDIF
                END DO

             ENDIF
!
!--          Loop for the section method: while FAST is catching up to PALM,
!--          PALM is frozen while the communication with FAST continues as before.
             DO WHILE ( current_time_fast(1) <= simulated_time + dt_palm )

#if defined( __parallel )
!
!--             Kommunizieren der Informationen über Positionen, die
!--             von PE0 erhalten wurden, sind auch an alle anderen PEs
!--             Dies laeuft ab wie in der Kopplung mit FAST
!
!--             Größe der Felder, deren Werte von PE0 allen anderen
!--             Prozessoren mitgeteilt werden
                size_array_type_1 = nturbines * fast_n_blades_max * fast_n_blade_elem_max * 3
                CALL MPI_BARRIER( comm2d, ierr )
                CALL MPI_BCAST( fast_tower_ref_pos_x(1), nturbines,                                &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_tower_ref_pos_y(1), nturbines,                                &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_tower_ref_pos_z(1), nturbines,                                &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_hub_center_pos(1,1), nturbines*3,                             &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_position_blade(1,1,1,1), size_array_type_1,                   &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_force_blade(1,1,1,1), size_array_type_1,                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( shaft_height_fast, nturbines,                                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BARRIER( comm2d, ierr )
#endif
!
!--             Wie im Code für die Kopplung von PALM und FAST:
!--             Interpolation des in PALM berechneten Geschwindigkeitsfeldes auf
!--             die Blattpositionen, die von FAST übergeben worden sind
!
!--             Allokierung temporaerer Felder für die interpolierten Geschwindigkeiten
!--             1. u_int_h: Geschwindigkeit an der Nabenposition, für jede Anlage drei Komponenten
                ALLOCATE(u_int_h(1:nturbines,1:3))
!
!--             2. u_int_b: Geschwindigkeit am Rotorblatt, für jedes Blattsegment eines der n Blätter einer Anlage drei Komponenten
                ALLOCATE(u_int_b(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3))
!
!--             Initialisierung der interpolierten Felder mit Nullwerten; dies wird eigentlich
!--             nicht benoetigt, da diese Felder auf jeden Fall auf jedem Prozessor im Laufe
!--             des Interpolationsverfahrens mit einem Wert belegt werden (wenn Punkt, auf den
!--             interpoliert werden soll, außerhalb des Teilmodellgebiets eines Prozessors
!--             liegt, wird dem entsprechenden Feldelement der Wert 0 zugeordnet)
                u_int_h(:,:) = 0.0_wp
                u_int_b(:,:,:,:) = 0.0_wp

                palm_hub_center_vel(:,:) = 0.0_wp
!
!--             count intermediate FAST time step
                l = l+1

                DO  i = 1, nturbines

                   IF ( myid == 0 .AND. debug_output )  THEN
                      WRITE(9,*) "rotspeed: ", i, rotspeed(i)
                      FLUSH(9)
                   ENDIF
!
!--                Bilineares Interpolationsverfahren:
!--                1. Interpolation des Geschwindigkeitsfeldes auf die Position des Hubs --> palm_hub_center_vel
!--                Nabenwindgeschwindigkeit an Turbine i, Interpolation von u:
                   ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                &
                          fast_hub_center_pos(i,1) ) * ddx )
                   jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                &
                          fast_hub_center_pos(i,2) - 0.5_wp * dy) * ddy )
                   kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                &
                          fast_hub_center_pos(i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                Auf einem Prozessorelement findet nur dann eine Interpolation statt,
!--                wenn alle Stützpunkte der Interpolation auf dem Prozessorelement zur
!--                Verfügung stehen
                   IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                         x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - (ii * dx)
                         y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - 0.5_wp * dy - jj * dy
                         z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - kk * dz(1) + 0.5_wp * dz(1)

                         aa = abs((dx-x)*(dy-y))
                         bb = abs((x)*(dy-y))
                         cc = abs((dx-x)*(y))
                         dd = x*y
                         gg = dx*dy

                         u_int_l = ( ( aa ) * u(kk,jj,ii)     +    &
                                     ( bb ) * u(kk,jj,ii+1)   +    &
                                     ( cc ) * u(kk,jj+1,ii)   +    &
                                     ( dd ) * u(kk,jj+1,ii+1) ) /  &
                                     ( gg )

                         u_int_u = ( ( aa ) * u(kk+1,jj,ii)     +    &
                                     ( bb ) * u(kk+1,jj,ii+1)   +    &
                                     ( cc ) * u(kk+1,jj+1,ii)   +    &
                                     ( dd ) * u(kk+1,jj+1,ii+1) ) /  &
                                     ( gg )

                         u_int_h(i,1) = (1/dz(1)) * ((dz(1)-z)*u_int_l + z*u_int_u)

                      ELSE !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                         u_int_h(i,1) = 0.0_wp
                      ENDIF !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                   ELSE !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
                      u_int_h(i,1) = 0.0_wp
                   ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                Nabenwindgeschwindigkeit an Turbine i, Interpolation von v:
                   ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                &
                          fast_hub_center_pos(i,1) - 0.5_wp * dx)  * ddx )
                   jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                &
                          fast_hub_center_pos(i,2) ) * ddy )
                   kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                &
                          fast_hub_center_pos(i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                Auf einem Prozessorelement findet nur dann eine Interpolation statt,
!--                wenn alle Stützpunkte der Interpolation auf dem Prozessorelement zur
!--                Verfügung stehen
                   IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                         x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - ii * dx - 0.5_wp * dx
                         y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - jj * dy
                         z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - kk * dz(1) + 0.5_wp * dz(1)

                         aa = abs((dx-x)*(dy-y))
                         bb = abs((x)*(dy-y))
                         cc = abs((dx-x)*(y))
                         dd = x*y
                         gg = dx*dy

                         v_int_l = ( ( aa ) * v(kk,jj,ii)     +    &
                                     ( bb ) * v(kk,jj,ii+1)   +    &
                                     ( cc ) * v(kk,jj+1,ii)   +    &
                                     ( dd ) * v(kk,jj+1,ii+1) ) /  &
                                     ( gg )

                         v_int_u = ( ( aa ) * v(kk+1,jj,ii)     +    &
                                     ( bb ) * v(kk+1,jj,ii+1)   +    &
                                     ( cc ) * v(kk+1,jj+1,ii)   +    &
                                     ( dd ) * v(kk+1,jj+1,ii+1) ) /  &
                                     ( gg )
!
!--                      ToDo: check, why the computed values are not used and why u_int_h(i,2) is set here
                         u_int_h(i,2) = 0.0_wp!(1/dz) * ((dz-z)*v_int_l + z*v_int_u)

                      ELSE !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                        u_int_h(i,2) = 0.0_wp
                      ENDIF !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                   ELSE !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
                      u_int_h(i,2) = 0.0_wp
                   ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                Nabenwindgeschwindigkeit an Turbine i, Interpolation von w:
                   ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                &
                        fast_hub_center_pos(i,1) - 0.5_wp * dx) * ddx )
                   jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                &
                        fast_hub_center_pos(i,2) - 0.5_wp * dy) * ddy )
                   kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                &
                        fast_hub_center_pos(i,3) ) / dz(1) )
!
!--                Auf einem Prozessorelement findet nur dann eine Interpolation statt,
!--                wenn alle Stützpunkte der Interpolation auf dem Prozessorelement zur
!--                Verfügung stehen
                   IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                         x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - ii * dx - 0.5_wp * dx
                         y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - jj * dy - 0.5_wp * dy
                         z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - kk * dz(1)

                         aa = abs((dx-x)*(dy-y))
                         bb = abs((x)*(dy-y))
                         cc = abs((dx-x)*(y))
                         dd = x*y
                         gg = dx*dy

                         w_int_l = ( ( aa ) * w(kk,jj,ii)     +    &
                                     ( bb ) * w(kk,jj,ii+1)   +    &
                                     ( cc ) * w(kk,jj+1,ii)   +    &
                                     ( dd ) * w(kk,jj+1,ii+1) ) /  &
                                     ( gg )

                         w_int_u = ( ( aa ) * w(kk+1,jj,ii)     +    &
                                     ( bb ) * w(kk+1,jj,ii+1)   +    &
                                     ( cc ) * w(kk+1,jj+1,ii)   +    &
                                     ( dd ) * w(kk+1,jj+1,ii+1) ) /  &
                                     ( gg )
!
!--                      ToDo: check, why the computed values are not used and why u_int_h(i,3) is set here
                         u_int_h(i,3) = 0.0_wp!(1/dz) * ((dz-z)*w_int_l + z*w_int_u)

                      ELSE !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                         u_int_h(i,3) = 0.0_wp
                      ENDIF !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                   ELSE !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
                      u_int_h(i,3) = 0.0_wp
                   ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )

                ENDDO !i = 1, nturbines
!
!--             2. Interpolation des Geschwindigkeitsfeldes auf die Position der Blattelemente --> palm_vel_blade

                palm_vel_blade(:,:,:,:) = 0.0_wp

                DO  i = 1, nturbines
                   DO  j = 1, fast_n_blades(i)
                      DO  k = 1, fast_n_blade_elem(i)
!
!--                      Geschwindigkeit an Element k des Blattes j der Turbine i, Interpolation von u:
                         ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +          &
                                fast_position_blade(k,j,i,1)  - 4 * fast_n_radius(i)) * ddx )! wind speed 2D in front of turb
                         jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +          &
                                fast_position_blade(k,j,i,2) - 0.5_wp * dy) * ddy )
                         kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +          &
                                fast_position_blade(k,j,i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                      Auf einem Prozessorelement findet nur dann eine Interpolation
!--                      statt, wenn alle Stützpunkte der Interpolation auf dem
!--                      Prozessorelement zur Verfügung stehen
                         IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                            IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                               x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +            &
                                    fast_position_blade(k,j,i,1) - 4 * fast_n_radius(i) - (ii * dx) ! wind speed 2D in front of turb
                               y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +            &
                                    fast_position_blade(k,j,i,2) - jj * dy - 0.5_wp * dy
                               z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +            &
                                    fast_position_blade(k,j,i,3) - kk * dz(1) + 0.5_wp * dz(1)

                               aa = abs((dx-x)*(dy-y))
                               bb = abs((x)*(dy-y))
                               cc = abs((dx-x)*(y))
                               dd = x*y
                               gg = dx*dy

                               u_int_l = ( ( aa ) * u(kk,jj,ii)     +    &
                                           ( bb ) * u(kk,jj,ii+1)   +    &
                                           ( cc ) * u(kk,jj+1,ii)   +    &
                                           ( dd ) * u(kk,jj+1,ii+1) ) /  &
                                         ( gg )

                               u_int_u = ( ( aa ) * u(kk+1,jj,ii)     +    &
                                           ( bb ) * u(kk+1,jj,ii+1)   +    &
                                           ( cc ) * u(kk+1,jj+1,ii)   +    &
                                           ( dd ) * u(kk+1,jj+1,ii+1) ) /  &
                                         ( gg )

                               u_int_b(k,j,i,1) = (1/dz(1)) * ((dz(1)-z)*u_int_l + z*u_int_u)

                            ELSE
                               u_int_b(k,j,i,1) = 0.0_wp
                            ENDIF
                         ELSE
                            u_int_b(k,j,i,1) = 0.0_wp
                         ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                      Geschwindigkeit an Element k des Blattes j der Turbine i, Interpolation von v:
                         ii = FLOOR( ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +         &
                                       fast_position_blade(k,j,i,1) -                              &
                                       4 * fast_n_radius(i) - 0.5_wp * dx)  * ddx )  ! wind speed 2D in front of turb

                         jj = FLOOR( ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +         &
                                       fast_position_blade(k,j,i,2) ) * ddy )

                         kk = FLOOR( ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +         &
                                       fast_position_blade(k,j,i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                      Auf einem Prozessorelement findet nur dann eine Interpolation
!--                      statt, wenn alle Stützpunkte der Interpolation auf dem
!--                      Prozessorelement zur Verfügung stehen
                         IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                            IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                               x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +            &
                                    fast_position_blade(k,j,i,1) -                                 &
                                    4 * fast_n_radius(i) - ii * dx - 0.5_wp * dx   ! wind speed 2D in front of turb

                               y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +            &
                                    fast_position_blade(k,j,i,2) - jj * dy

                               z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +            &
                                    fast_position_blade(k,j,i,3) - kk * dz(1) + 0.5_wp * dz(1)

                               aa = abs((dx-x)*(dy-y))
                               bb = abs((x)*(dy-y))
                               cc = abs((dx-x)*(y))
                               dd = x*y
                               gg = dx*dy

                               v_int_l = ( ( aa ) * v(kk,jj,ii)     +                              &
                                           ( bb ) * v(kk,jj,ii+1)   +                              &
                                           ( cc ) * v(kk,jj+1,ii)   +                              &
                                           ( dd ) * v(kk,jj+1,ii+1) ) /                            &
                                           ( gg )

                               v_int_u = ( ( aa ) * v(kk+1,jj,ii)     +                            &
                                           ( bb ) * v(kk+1,jj,ii+1)   +                            &
                                           ( cc ) * v(kk+1,jj+1,ii)   +                            &
                                           ( dd ) * v(kk+1,jj+1,ii+1) ) /                          &
                                         ( gg )

                               u_int_b(k,j,i,2) = (1/dz(1)) * ((dz(1)-z)*v_int_l + z*v_int_u)


                            ELSE
                               u_int_b(k,j,i,2) = 0.0_wp
                            ENDIF
                         ELSE
                            u_int_b(k,j,i,2) = 0.0_wp
                         ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                      Geschwindigkeit im Element k des Blattes k der Turbine i,
!--                      Interpolation von w:
                         ii = FLOOR( ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +         &
                                       fast_position_blade(k,j,i,1) -                              &
                                       4 * fast_n_radius(i) - 0.5_wp * dx) * ddx )  ! wind speed 2D in front of turb

                         jj = FLOOR( ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +         &
                                       fast_position_blade(k,j,i,2) - 0.5_wp * dy) * ddy )

                         kk = FLOOR( ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +         &
                                       fast_position_blade(k,j,i,3) ) / dz(1) )
!
!--                      Auf einem Prozessorelement findet nur dann eine Interpolation
!--                      statt, wenn alle Stützpunkte der Interpolation auf dem
!--                      Prozessorelement zur Verfügung stehen
                         IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                            IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                               x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +            &
                                    fast_position_blade(k,j,i,1) -                                 &
                                    4 * fast_n_radius(i) - ii * dx - 0.5_wp * dx   ! wind speed 2D in front of turb

                               y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +            &
                                    fast_position_blade(k,j,i,2) - jj * dy - 0.5_wp * dy

                               z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +            &
                                    fast_position_blade(k,j,i,3) - kk * dz(1)

                               aa = abs((dx-x)*(dy-y))
                               bb = abs((x)*(dy-y))
                               cc = abs((dx-x)*(y))
                               dd = x*y
                               gg = dx*dy

                               w_int_l = ( ( aa ) * w(kk,jj,ii)     +                              &
                                           ( bb ) * w(kk,jj,ii+1)   +                              &
                                           ( cc ) * w(kk,jj+1,ii)   +                              &
                                           ( dd ) * w(kk,jj+1,ii+1) ) /                            &
                                         ( gg )

                               w_int_u = ( ( aa ) * w(kk+1,jj,ii)     +                            &
                                           ( bb ) * w(kk+1,jj,ii+1)   +                            &
                                           ( cc ) * w(kk+1,jj+1,ii)   +                            &
                                           ( dd ) * w(kk+1,jj+1,ii+1) ) /                          &
                                         ( gg )

                               u_int_b(k,j,i,3) = (1/dz(1)) * ((dz(1)-z)*w_int_l + z*w_int_u)


                            ELSE
                               u_int_b(k,j,i,3) = 0.0_wp
                            ENDIF
                         ELSE
                            u_int_b(k,j,i,3) = 0.0_wp
                         ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )

                      ENDDO !k = 1, fast_n_blade_elem(i)
                   ENDDO !j = 1, fast_n_blades(i)
                ENDDO !i = 1, nturbines


#if defined( __parallel )
!
!--             In Vorbereitung auf den Austausch der Informationen ueber die
!--             Geschwindigkeiten an den aus FAST erhaltenen Informationen, werden
!--             die Groessen der per MPI-Befehl gesendeten Felder bestimmt
                size_array_type_1 = nturbines*fast_n_blades_max*fast_n_blade_elem_max*3
!
!--             Austausch der Informationen über die Geschwindigkeiten an den
!--             von FAST erhaltenen Positionen mittels MPI-Befehlen
!--             Wichtig ist, dass PE0 alle Informationen erhält, da ueber dieses
!--             Prozessorelement dann spaeter die Kommunikation mit dem FAST-Programm ablaeuft
!--             Sicherheitshalber werden vor dem Austausch von Informationen zwischen
!--             aller Prozessorelementen noch einmal alle Prozessorelemente synchronisiert
                CALL MPI_BARRIER( comm2d, ierr )
!
!--             Austausch der Komponenten der Geschwindigkeit an der Rotornabe
!--             Informationen darueber stehen danach im Feld palm_hub_center_vel zur Verfuegung
                CALL MPI_ALLREDUCE( u_int_h(1,1), palm_hub_center_vel(1,1), nturbines*3,           &
                                    MPI_REAL, MPI_SUM, comm2d, ierr )
!
!--             Austausch der Komponenten der Geschwindigkeiten an den Blattelementen
!--             Informationen darueber stehen danach im Feld palm_vel_blade zur Verfuegung
                CALL MPI_ALLREDUCE( u_int_b(1,1,1,1), palm_vel_blade(1,1,1,1), size_array_type_1,  &
                                    MPI_REAL, MPI_SUM, comm2d, ierr )
!
!--             Erneute Synchronisation der Prozessorelemente
                CALL MPI_BARRIER( comm2d, ierr )
#else
                palm_hub_center_vel = u_int_h
                palm_vel_blade = u_int_b
#endif
!
!--             Da die temporaeren Felder mit den interpolierten Geschwindigkeitswerten nicht
!--             mehr benoetigt werden, koennen diese nunmehr deallokiert werden
                DEALLOCATE( u_int_h, u_int_b )
!
!--             PALM-FAST-Interaktion:
!--             Die aus den Geschwindigkeitsfeldern von PALM interpolierten Geschwindigkeiten an
!--             den Positionen, die aus FAST erhalten wurden, werden nun mittels Socket
!--             von PALM an FAST übergeben
!
!--             Im Wesentlichen zwei Schritte:
!--             1. Schritt: Daten zusammenpacken in das Format, das FAST erwartet
!--             2. Schritt: PALM-Client starten, der sich mit FAST-Server verbindet und so die
!--                         Daten bereitstellt
!--             Da der Austausch mit FAST wieder nur ueber PE0 erfolgt, wird im Folgenden
!--             wieder eine Fallunterscheidung nach Prozessorelementen getroffen
                IF ( myid == 0 )  THEN
                   IF ( .NOT. terminate_run )  THEN

                      IF ( debug_output )  THEN
                         WRITE(9,*) 'Sending velocities to FAST server(s)...'
                         FLUSH(9)
                      ENDIF

#if defined( __fastv8 )
                      iflag = commclient(0, 3)
#else
                      iflag = -1
#endif

                      IF ( iflag < 0 )  THEN
                         WRITE( message_string, * )  'Unable to send velocities to FAST server(s)',&
                                '. Please check if any of the FAST servers finished the simulation.'
                         CALL message( 'f8c_actions', 'F8C0019', 0, 2, 0, 6, 0 )
                         terminate_run = .TRUE.
                      ELSEIF ( iflag == 1 )  THEN
                         WRITE( message_string, * )                                                &
                                'At least one FAST server finished the simulation'
                         CALL message( 'f8c_actions', 'F8C0020', 0, 2, 0, 6, 0 )
                         terminate_run = .TRUE.
                      ELSEIF ( iflag == 0 )  THEN
                         IF ( debug_output )  THEN
                            WRITE(9,*) 'Velocities successfully broadcasted'
                            FLUSH(9)
                         ENDIF
                      ENDIF !( iflag < 0 )

                   ENDIF !(.NOT. terminate_run )
                ENDIF !( myid == 0 )

#if defined( __parallel )
                   CALL MPI_BARRIER( comm2d, ierr )
                   CALL MPI_BCAST( terminate_run, 1, MPI_LOGICAL, 0, comm2d, ierr )
                   CALL MPI_BARRIER( comm2d, ierr )
#endif
                   IF ( terminate_run )  THEN
                      RETURN
                   ENDIF
!
!--             Restoring the positions and forces to have the starting position within
!--             the sector and the forces of the line in the middle of the sector.
                IF ( myid == 0 )  THEN
                   IF ( l == 1 )  THEN

                      DO i = 1, nturbines

                         fast_hub_center_pos_old(i,1) = fast_hub_center_pos(i,1)
                         fast_hub_center_pos_old(i,2) = fast_hub_center_pos(i,2)
                         fast_hub_center_pos_old(i,3) = fast_hub_center_pos(i,3)

                         DO m = 1, fast_n_blades(i)
                            DO n = 1, fast_n_blade_elem(i)

                               fast_position_blade_old(n,m,i,1) = fast_position_blade(n,m,i,1)
                               fast_position_blade_old(n,m,i,2) = fast_position_blade(n,m,i,2)
                               fast_position_blade_old(n,m,i,3) = fast_position_blade(n,m,i,3)

                            ENDDO
                         ENDDO

                      ENDDO

                   ENDIF

                   IF ( l == ceiling( (dt_palm/dt_fast ) / 2 ) )  THEN

                      DO i = 1, nturbines
                         DO m = 1, fast_n_blades(i)
                            DO n = 1, fast_n_blade_elem(i)
                               fast_force_blade_old(n,m,i,1) = fast_force_blade(n,m,i,1)
                               fast_force_blade_old(n,m,i,2) = fast_force_blade(n,m,i,2)
                               fast_force_blade_old(n,m,i,3) = fast_force_blade(n,m,i,3)
                            ENDDO
                         ENDDO
                      ENDDO

                   ENDIF

                ENDIF !(myid == 0)

#if defined( __parallel )
!
!--             Wegen unterschiedlicher Aktivitaeten der einzelnen Prozessorelemente
!--             im Vorfeld werden die Prozessorelemente an dieser Stelle
!--             sicherheitshalber synchronisiert
                CALL MPI_BARRIER( comm2d, ierr )
!
!--             Nach der Kommunikation mit FAST sind Informationen ueber
!--             Kraefte und Positionen zunaechst nur auf PE0 bekannt. Mittels
!--             MPI_BCAST sollen die Informationen auch auf den anderen
!--             Prozessorelementen bekanntgemacht werden
!--             Vorher ist noch die Groesse der mittels MPI_BCAST
!--             auszutauschenden Felder zu bestimmen.
                size_array_type_1 = nturbines * fast_n_blades_max * fast_n_blade_elem_max * 3
!
!--             Broadcasting von fast_hub_center_pos (nturbines,3)
!--             (Position der Rotornabe (x-, y- und z-Komponente))
                CALL MPI_BCAST( fast_hub_center_pos(1,1), nturbines*3,                             &
                                MPI_REAL, 0, comm2d, ierr )
!
!--             Broadcasting von fast_position_blade
!--             (fast_n_blade_elem_max,fast_n_blades_max,nturbines,3)
!--             (Positionen der Blattelemente (x-, y- und z-Komponenten))
                CALL MPI_BCAST( fast_position_blade(1,1,1,1), size_array_type_1,                   &
                                MPI_REAL, 0, comm2d, ierr )
!
!--             Broadcasting von fast_force_blade
!--             (fast_n_blade_elem_max,fast_n_blades_max,nturbines,3)
!--             (Kraefte an den Blattelementen (x-, y- und z-Komponenten))
                CALL MPI_BCAST( fast_force_blade(1,1,1,1), size_array_type_1,                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_position_blade_old(1,1,1,1), size_array_type_1,               &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_force_blade_old(1,1,1,1), size_array_type_1,                  &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( current_time_fast, nturbines,                                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDDO ! DO WHILE ( current_time_fast(1) <= simulated_time + dt_palm )
!
!--          Calculating the angle of the sector and the rotation matrix,
!--          to rotate the positions starting from the first line in the sector
             IF ( myid == 0 )  THEN

                DO  i = 1, nturbines

                   IF ( palm_dt_3d )  THEN
                      sector_angle(i) = rotspeed(i) * dt_palm
                      n_sector(i) =  FLOOR( dt_palm/dt_fast )
                      alpha_rot(i) =  sector_angle(i) / n_sector(i)
                   ELSE
                      sector_angle_deg(i) = 120
                      n_sector(i) = MAX(n_dt_fast, 1)
                      alpha_rot(i) = ( sector_angle_deg(i) / n_dt_fast ) * pi / 180.0_wp
                   ENDIF
!
!--                Drehmatrix:
                   r_n(i,1,:) = (/                                                                 &
                      (shaft_coordinates(i,1)**2)*(1-COS(alpha_rot(i)))+COS(alpha_rot(i)),         &
                      shaft_coordinates(i,1)*shaft_coordinates(i,2)*(1-COS(alpha_rot(i))) -        &
                         shaft_coordinates(i,3)*sin(alpha_rot(i)),                                 &
                      shaft_coordinates(i,1)*shaft_coordinates(i,3)*(1-COS(alpha_rot(i))) +        &
                         shaft_coordinates(i,2)*sin(alpha_rot(i)) &
                                /)
                   r_n(i,2,:) = (/                                                                 &
                      shaft_coordinates(i,2)*shaft_coordinates(i,1)*(1-COS(alpha_rot(i))) +        &
                         shaft_coordinates(i,3)*sin(alpha_rot(i)),                                 &
                      (shaft_coordinates(i,2)**2)*(1-COS(alpha_rot(i)))+COS(alpha_rot(i)),         &
                      shaft_coordinates(i,2)*shaft_coordinates(i,3)*(1-COS(alpha_rot(i))) -        &
                         shaft_coordinates(i,1)*sin(alpha_rot(i))                                  &
                                /)
                   r_n(i,3,:) = (/                                                                 &
                      shaft_coordinates(i,3)*shaft_coordinates(i,1)*(1-COS(alpha_rot(i))) -        &
                         shaft_coordinates(i,2)*sin(alpha_rot(i)),                                 &
                      shaft_coordinates(i,3)*shaft_coordinates(i,2)*(1-COS(alpha_rot(i))) +        &
                         shaft_coordinates(i,1)*sin(alpha_rot(i)),                                 &
                      (shaft_coordinates(i,3)**2)*(1-COS(alpha_rot(i)))+COS(alpha_rot(i))          &
                                /)

                ENDDO !DO  i = 1, nturbines

             ENDIF !( myid == 0 )

#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
             CALL MPI_BCAST( n_sector, nturbines, MPI_INTEGER, 0, comm2d, ierr )
             CALL MPI_BCAST( r_n, nturbines*9, MPI_REAL, 0, comm2d, ierr )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             n_sector_max = MAXVAL( n_sector(:) )

             ALLOCATE( x_prime(1:fast_n_blade_elem_max,                                            &
                               1:fast_n_blades_max,                                                &
                               1:nturbines,                                                        &
                               1:n_sector_max) )
             ALLOCATE( y_prime(1:fast_n_blade_elem_max,                                            &
                               1:fast_n_blades_max,                                                &
                               1:nturbines,                                                        &
                               1:n_sector_max) )
             ALLOCATE( z_prime(1:fast_n_blade_elem_max,                                            &
                               1:fast_n_blades_max,                                                &
                               1:nturbines,                                                        &
                               1:n_sector_max) )

             x_prime(:,:,:,:) = 0.0_wp
             y_prime(:,:,:,:) = 0.0_wp
             z_prime(:,:,:,:) = 0.0_wp

             ALLOCATE( x_prime_sec(1:fast_n_blade_elem_max,                                        &
                                   1:n_sector_max*fast_n_blades_max,                               &
                                   1:nturbines) )
             ALLOCATE( y_prime_sec(1:fast_n_blade_elem_max,                                        &
                                   1:n_sector_max*fast_n_blades_max,                               &
                                   1:nturbines) )
             ALLOCATE( z_prime_sec(1:fast_n_blade_elem_max,                                        &
                                   1:n_sector_max*fast_n_blades_max,                               &
                                   1:nturbines) )

             x_prime_sec(:,:,:) = 0.0_wp
             y_prime_sec(:,:,:) = 0.0_wp
             z_prime_sec(:,:,:) = 0.0_wp
!
!--          Calculating the positions of the points on the line samples distributed over the three blade sectors
             DO  i = 1, nturbines
                DO  j = 1, fast_n_blades(i)
                   DO  k = 1, fast_n_blade_elem(i)

                      x = fast_tower_ref_pos_x(i) +                                                &
                          fast_position_blade_old(k,j,i,1)
                      y = fast_tower_ref_pos_y(i) +                                                &
                          fast_position_blade_old(k,j,i,2)
                      z = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                      &
                          fast_position_blade_old(k,j,i,3) - shaft_height_fast(i)
!
!--                   Multiplication of the points with the rotation matrix -> positions of the lines within the sector
                      x_prime(k,j,i,1) = x
                      y_prime(k,j,i,1) = y
                      z_prime(k,j,i,1) = z
                      DO m = 2, n_sector(i)

                         tmp_res = MATMUL( r_n(i,:,:),                                             &
                            (/ x_prime(k,j,i,m-1), y_prime(k,j,i,m-1), z_prime(k,j,i,m-1) /)       &
                                         )

                         x_prime(k,j,i,m) = tmp_res(1)
                         y_prime(k,j,i,m) = tmp_res(2)
                         z_prime(k,j,i,m) = tmp_res(3)

                      ENDDO !DO m = 2, n_sector(i)

                   ENDDO !DO  k = 1, fast_n_blade_elem(i)
                ENDDO !DO  j = 1, fast_n_blades(i)
!
!--             Positions in the PALM grid
                x_prime(:,:,i,:) = x_prime(:,:,i,:) + palm_tower_ref_pos_x(i)
                y_prime(:,:,i,:) = y_prime(:,:,i,:) + palm_tower_ref_pos_y(i)
                z_prime(:,:,i,:) = z_prime(:,:,i,:) + shaft_height_fast(i)

             ENDDO !DO i = 1, nturbines

             DO  i = 1, nturbines
                DO  j = 1, fast_n_blades(i)
                   DO  k = 1, fast_n_blade_elem(i)
                      DO m = 1, n_sector(i)
                         x_prime_sec(k,(m-1)*fast_n_blades(i) + j,i) = x_prime(k,j,i,m)
                         y_prime_sec(k,(m-1)*fast_n_blades(i) + j,i) = y_prime(k,j,i,m)
                         z_prime_sec(k,(m-1)*fast_n_blades(i) + j,i) = z_prime(k,j,i,m)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO !DO i = 1, nturbines

             !-Verschmieren und Aufsummieren der Kraefte

             !-Einfuehren eines Regulierungskernels
             eps_kernel = reg_fac * dx
             eps_kernel_2 = -reg_fac * reg_fac * dx * dx

             nenner = (eps_kernel**3.0) * (pi**(3.0/2.0))
!
!--          Informationen ueber die zusaetzlichen durch die Windenergieanlage verursachten
!--          Kraefte werden an jedem Gitterpunkt des Modellgebiets benoetigt
!--          In einer Schleife ueber alle Gitterpunkte des Teilmodellgebiets wird im Folgenden
!--          die Verschmierung der Kraefte durchgefuehrt, so dass jedem Gitterpunkt eine
!--          zusaetzliche Kraft zugeordnet werden kann
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      force_x(k,j,i) = 0.0_wp
                      force_y(k,j,i) = 0.0_wp
                      force_z(k,j,i) = 0.0_wp
                   ENDDO
                ENDDO
             ENDDO

             DO  o = 1, nturbines
!
!--             check if coordinates are even in the limit of the processor
                IF ( .NOT. ( nxlg >= fboxcorners(o, 1, 2) )  .OR.                                  &
                           ( nxrg <= fboxcorners(o, 1, 1) ) )  THEN

                   IF ( .NOT. ( nysg >= fboxcorners(o, 2, 2) )  .OR.                               &
                              ( nyng <= fboxcorners(o, 2, 1) ) )  THEN

                      IF ( .NOT. ( nzb >= fboxcorners(o, 3, 2) )  .OR.                             &
                                 ( nzt <= fboxcorners(o, 3, 1) ) )  THEN

                         tower_ref_pos_x(o) = palm_tower_ref_pos_x(o) + fast_tower_ref_pos_x(o)
                         tower_ref_pos_y(o) = palm_tower_ref_pos_y(o) + fast_tower_ref_pos_y(o)
                         tower_ref_pos_z(o) = palm_tower_ref_pos_z(o) + fast_tower_ref_pos_z(o)

                         DO  i = nxlg, nxrg
                            IF ( ( i >= fboxcorners(o, 1, 1) )  .AND.                              &
                                 ( i <= fboxcorners(o, 1, 2) ) )  THEN

                               DO  j = nysg, nyng
                                  IF ( ( j >= fboxcorners(o, 2, 1) )  .AND.                        &
                                       ( j <= fboxcorners(o, 2, 2) ) )  THEN

                                     DO  k = nzb, nzt
                                        IF ( ( k >= fboxcorners(o, 3, 1) )  .AND.                  &
                                             ( k <= fboxcorners(o, 3, 2) ) )  THEN

                                           DO  m = 1, fast_n_blades(o)
                                              DO  n = 1, fast_n_blade_elem(o)
                                                 DO n_s = 1, n_sector(o)
!
!--                                                 Bestimmung der Distanz zwischen Blattelement und dem Gitterpunkt des u-grids
                                                    x  = x_prime(n,m,o,n_s) - i * dx
                                                    y  = y_prime(n,m,o,n_s) - j * dy - 0.5_wp * dy
                                                    z  = z_prime(n,m,o,n_s) - k * dz(1)            &
                                                                            + 0.5_wp * dz(1)

                                                    distance_hub = x**2_iwp + y**2_iwp + z**2_iwp

                                                    expargument = INT(                             &
                                                       - distance_hub / eps_kernel_2 * 1000.0_wp   &
                                                                     )

                                                    IF ( expargument < 10000001_iwp )  THEN
                                                       kernel_value = expvalue(expargument)/nenner
                                                    ELSE
                                                       kernel_value = 0.0_wp
                                                    ENDIF
!
!--                                                 Berechnung der Kraft in x-Richtung auf dem u-grid
                                                    IF ( mod(m,3) == 0 )  THEN
                                                       l = 3
                                                    ELSE
                                                       l = mod(m,3)
                                                    ENDIF

                                                    force_x(k,j,i) = force_x(k,j,i) -              &
                                                                     ( 1/float(n_sector(o)) ) *    &
                                                                     fast_force_blade_old(n,m,o,1) &
                                                                     * kernel_value
!
!--                                                 Bestimmung der Distanz zwischen Blattelement und dem Gitterpunkt des v-grids
                                                    x  = x_prime(n,m,o,n_s) - i * dx - 0.5_wp * dx
                                                    y  = y_prime(n,m,o,n_s) - j * dy
                                                    z  = z_prime(n,m,o,n_s) - k * dz(1)            &
                                                                            + 0.5_wp * dz(1)

                                                    distance_hub = x**2_iwp + y**2_iwp + z**2_iwp

                                                    expargument = INT(                             &
                                                       - distance_hub / eps_kernel_2 * 1000.0_wp   &
                                                                     )

                                                    IF ( expargument < 10000001_iwp )  THEN
                                                       kernel_value = expvalue(expargument)/nenner
                                                    ELSE
                                                       kernel_value = 0.0_wp
                                                    ENDIF
!
!--                                                 Berechnung der Kraft in y-Richtung auf dem v-grid
                                                    force_y(k,j,i) = force_y(k,j,i) -              &
                                                                     ( 1/float(n_sector(o)) ) *    &
                                                                     fast_force_blade_old(n,m,o,2) &
                                                                     * kernel_value
!
!--                                                 Bestimmung der Distanz zwischen Blattelement und dem Gitterpunkt des w-grids
                                                    x  = x_prime(n,m,o,n_s) - i * dx - 0.5_wp * dx
                                                    y  = y_prime(n,m,o,n_s) - j * dy - 0.5_wp * dy
                                                    z  = z_prime(n,m,o,n_s) - k * dz(1)

                                                    distance_hub = x**2_iwp + y**2_iwp + z**2_iwp

                                                    expargument = INT(                             &
                                                       - distance_hub / eps_kernel_2 * 1000.0_wp   &
                                                                     )

                                                    IF ( expargument < 10000001_iwp )  THEN
                                                       kernel_value = expvalue(expargument)/nenner
                                                    ELSE
                                                       kernel_value = 0.0_wp
                                                    ENDIF
!
!--                                                 Berechnung der Kraft in z-Richtung auf dem w-grid
                                                    force_z(k,j,i) = force_z(k,j,i) -              &
                                                                     ( 1/float(n_sector(o)) ) *    &
                                                                     fast_force_blade_old(n,m,o,3) &
                                                                     * kernel_value
                                                 ENDDO
                                              ENDDO
                                           ENDDO

                                           IF ( u(k,j,i) >= 0.0_wp )  THEN
                                              thrust_tower_x   = 0.5_wp * turb_C_d_tow(o) *        &
                                                              tower_area_x_4d(o,k,j,i) /           &
                                                              ( dx * dy * dz(1) ) * u(k,j,i)**2
                                              thrust_tower_y   = 0.5_wp * turb_C_d_tow(o) *        &
                                                              tower_area_y_4d(o,k,j,i) /           &
                                                              ( dx * dy * dz(1) ) * v(k,j,i)**2
                                           ELSE
                                              thrust_tower_x   = 0.0_wp
                                              thrust_tower_y   = 0.0_wp
                                           ENDIF

                                           force_x(k,j,i) = force_x(k,j,i) - thrust_tower_x
                                           force_y(k,j,i) = force_y(k,j,i) - thrust_tower_y

                                        ENDIF !((k >= fboxcorners(o, 3, 1)) .AND. (k <= fboxcorners(o, 3, 2)))
                                     ENDDO !DO k = nzb, nzt

                                  ENDIF !((j >= fboxcorners(o, 2, 1)) .AND. (j <= fboxcorners(o, 2, 2)))
                               ENDDO !DO j = nysg, nyng

                            ENDIF !((i >= fboxcorners(o, 1, 1)) .AND. (i <= fboxcorners(o, 1, 2)))
                         ENDDO !DO  i = nxlg, nxrg

                      ENDIF !(.NOT.((nzb >= fboxcorners(o, 3, 2)) .OR. (nzt <= fboxcorners(o, 3, 1))))
                   ENDIF !(.NOT.((nysg >= fboxcorners(o, 2, 2)) .OR. (nyng <= fboxcorners(o, 2, 1))))
                ENDIF !(.NOT.((nxlg >= fboxcorners(o, 1, 2)) .OR. (nxrg <= fboxcorners(o, 1, 1))))

             ENDDO !DO o = 1, nturbines

#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif

             DEALLOCATE( x_prime, y_prime, z_prime, stat=ierr )
             DEALLOCATE( x_prime_sec, y_prime_sec, z_prime_sec, stat=ierr )

          ENDIF !( simulated_time >= time_turbine_on )


       CASE ( 'after_integration' )
!
!--       Total thrust and torque components calculated above:
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb, nzt !nzb_u_inner(j,i)+1, nzt
                   thrx(k,j,i) = force_x(k,j,i)
                   tory(k,j,i) = force_y(k,j,i)
                   torz(k,j,i) = force_z(k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'u-tendency' )

          IF ( simulated_time >= time_turbine_on )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt !nzb_u_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + force_x(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'v-tendency' )

          IF ( simulated_time >= time_turbine_on )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt !nzb_v_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + force_y(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'w-tendency' )

          IF ( simulated_time >= time_turbine_on )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt !nzb_w_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + force_z(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

    CALL cpu_log( log_point(24), 'f8c_actions', 'stop' )

 END SUBROUTINE f8c_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_actions_ij( i, j, location )


    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

!
!-- Here the user-defined actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

          IF ( simulated_time >= time_turbine_on )  THEN
             DO  k = nzb, nzt !nzb_u_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + force_x(k,j,i)
             ENDDO
          ENDIF

       CASE ( 'v-tendency' )

          IF ( simulated_time >= time_turbine_on )  THEN
             DO  k = nzb, nzt !nzb_v_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + force_y(k,j,i)
             ENDDO
          ENDIF

       CASE ( 'w-tendency' )

          IF ( simulated_time >= time_turbine_on )  THEN
             DO  k = nzb, nzt !nzb_w_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + force_z(k,j,i)
             ENDDO
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE f8c_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k) and sets the grid on which it is defined. Allowed values for grid are "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )


    CHARACTER(LEN=*), INTENT(INOUT) ::  grid      !< name of vertical grid
    CHARACTER(LEN=*), INTENT(IN)    ::  mode      !< either 'xy', 'xz' or 'yz'
    CHARACTER(LEN=*), INTENT(IN)    ::  variable  !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av      !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp), INTENT(IN) ::  nzb_do  !< lower limit of the domain (usually nzb)
    INTEGER(iwp), INTENT(IN) ::  nzt_do  !< upper limit of the domain (usually nzt+1)

    INTEGER(iwp) ::  i  !< grid index along x-direction
    INTEGER(iwp) ::  j  !< grid index along y-direction
    INTEGER(iwp) ::  k  !< grid index along z-direction

    LOGICAL, INTENT(INOUT) ::  found  !<
    LOGICAL, INTENT(INOUT) ::  two_d  !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.  local_pf(nxl,nys,nzb_do) == 0.0_wp  .OR.  two_d )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'thrx_xy', 'thrx_xz', 'thrx_yz' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   local_pf(i,j,k) = thrx(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'tory_xy', 'tory_xz', 'tory_yz' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   local_pf(i,j,k) = tory(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'torz_xy', 'torz_xz', 'torz_yz' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   local_pf(i,j,k) = torz(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          IF ( mode == 'xy' ) grid = 'zu'


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE f8c_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av     !<
    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  k      !<
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found  !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<


!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.  local_pf(nxl,nys,nzb_do) == 0.0_wp )  CONTINUE

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'thrx' )
          IF ( av == 0 )  THEN

             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = 0.5_wp * ( thrx(k,j,i) + thrx(k,j,i+1) )
                     ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = thrx(k,j,i)
                     ENDDO
                   ENDDO
                ENDDO
             ENDIF

          ENDIF
          found = .TRUE.

       CASE ( 'tory' )
          IF ( av == 0 )  THEN

             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = 0.5_wp * ( tory(k,j,i) + tory(k,j+1,i) )
                     ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = tory(k,j,i)
                     ENDDO
                   ENDDO
                ENDDO
             ENDIF

          ENDIF
          found = .TRUE.

       CASE ( 'torz' )
          IF ( av == 0 )  THEN

             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = 0.5_wp * ( torz(k,j,i) + torz(k-1,j,i) )
                     ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = torz(k,j,i)
                     ENDDO
                   ENDDO
                ENDDO
             ENDIF

          ENDIF
          found = .TRUE.

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE f8c_data_output_3d


 END MODULE fastv8_coupler_mod
