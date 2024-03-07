!> @file wind_turbine_model_mod.f90
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
! Copyright 2009-2021 Carl von Ossietzky Universitaet Oldenburg
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> This module introduces anthropogenic heat emissions from external sources into PALM. It considers
!> heat from buildings, traffic, and other point sources. 
!> AH-profiles for BUILDINGS are typically obtained from a building energy model (BEM) 
!> that replicates PALM's buildings in terms of geometry and materials, but applies a different 
!> heating or cooling system to them.
!> TRAFFIC AH-profiles are meant to be obtained from a traffic model that replicates PALM's streets 
!> and incorporates traffic flow data. 
!> POINT SORUCES are a more generic data fromat that can be used to represent any other source of 
!> anthropogenic heat: e.g. incineration plants, power plants, or industrial complexes.
!
!--------------------------------------------------------------------------------------------------!
 MODULE anthropogenic_heat_mod

#if defined( __parallel )
    USE MPI
#endif
    
    USE arrays_3d,                                                                                     &
       ONLY:  tend,                                                                                    &
              u,                                                                                       &
              v,                                                                                       &
              w,                                                                                       &
              zu,                                                                                      &
              zw
    
    USE basic_constants_and_equations_mod,                                                             &
       ONLY:  pi
    
    USE control_parameters,                                                                            &
       ONLY:  coupling_char,                                                                           &
              debug_output,                                                                            &
              dt_3d,                                                                                   &
              dz,                                                                                      &
              end_time,                                                                                &
              external_anthropogenic_heat,                                                             &
              initializing_actions,                                                                    &
              message_string,                                                                          &
              origin_date_time,                                                                        &
              restart_data_format_output,                                                              &
              time_since_reference_point,                                                              &
              wind_turbine
    
    USE cpulog,                                                                                        &
       ONLY:  cpu_log,                                                                                 &
              log_point_s
    
    USE data_output_module
   
    USE grid_variables,                                                                                &
       ONLY:  ddx,                                                                                     &
              dx,                                                                                      &
              ddy,                                                                                     &
              dy
    
    USE indices,                                                                                       &
       ONLY:  nbgp,                                                                                    &
              nx,                                                                                      &
              nxl,                                                                                     &
              nxlg,                                                                                    &
              nxr,                                                                                     &
              nxrg,                                                                                    &
              ny,                                                                                      &
              nyn,                                                                                     &
              nyng,                                                                                    &
              nys,                                                                                     &
              nysg,                                                                                    &
              nz,                                                                                      &
              nzb,                                                                                     &
              nzt,                                                                                     &
              topo_flags
    
    USE kinds
    
    USE netcdf_data_input_mod,                                                                         &
       ONLY:  char_fill,                                                                               &
              check_existence,                                                                         &
              close_input_file,                                                                        &
              get_variable,                                                                            &
              get_attribute,                                                                           &
              get_dimension_length,                                                                    &
              init_model,                                                                              &
              input_pids_wtm,                                                                          &
              input_pids_ah,                                                                           &
              inquire_num_variables,                                                                   &
              inquire_variable_names,                                                                  &
              input_file_wtm,                                                                          &
              input_file_ah,                                                                           &
              num_var_pids,                                                                            &
              open_read_file,                                                                          &
              pids_id,                                                                                 &
              vars_pids
    
    USE pegrid
    
    USE restart_data_mpi_io_mod,                                                                       &
       ONLY:  rrd_mpi_io_global_array,                                                                 &
              wrd_mpi_io_global_array

    USE surface_mod,                                                                                   &
       ONLY:  surf_usm,                                                                                &
              surf_lsm
    

    !
    !-- Define data types
    TYPE int_dimension
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: val      !< 1d integer-value array
       LOGICAL ::  from_file = .FALSE.                     !< flag indicating whether an input variable is available and read from file 
                                                           !< or default value is used
    END TYPE int_dimension
    
    TYPE real_2d_matrix
       REAL(wp) ::  fill                                  !< fill value
       REAL(wp), DIMENSION(:,:), ALLOCATABLE :: val       !< 2d real-value matrix
       LOGICAL ::  from_file = .FALSE.                    !< flag indicating whether an input variable is available and read from file 
                                                          !< or default value is used
    END TYPE real_2d_matrix

    TYPE point_coordinates
       REAL(wp) ::  x_fill                                !< fill value for x-coordinates
       REAL(wp) ::  y_fill                                !< fill value for y-coordinates
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  x          !< x-coordinates of point sources
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  y          !< y-coordinates of point sources
       LOGICAL ::  from_file = .FALSE.                    !< flag indicating whether an input variable is available and read from file 
                                                          !< or default value is used
    END TYPE point_coordinates

       
    IMPLICIT NONE
    
    PRIVATE
    
    ! ------- NEW ANTHROPOGENIC HEAT MODEL VARIABLES ------- !

    TYPE(int_dimension) :: building_ids          !< ids of buildings with anthropogenic heat profiles
    TYPE(int_dimension) :: street_ids            !< ids of streets with anthropogenic heat profiles
    TYPE(int_dimension) :: point_ids             !< ids of point sources with anthropogenic heat profiles
    REAL(wp), DIMENSION(:), ALLOCATABLE :: t     !< time steps


    TYPE(real_2d_matrix) :: building_ah       !< anthropogenic heat profiles for buildings
    TYPE(real_2d_matrix) :: street_ah         !< anthropogenic heat profiles for streets
    TYPE(real_2d_matrix) :: point_ah          !< anthropogenic heat profiles for point sources
    

    TYPE(point_coordinates) :: point_coords      !< exact coordinates of point sources


    ! ------- OLD WIND TURBINE MODEL VARIABLES ------- !
    
        CHARACTER(LEN=800) ::  dom_error_message  !< error message returned by the data-output module
        CHARACTER(LEN=100) ::  variable_name      !< name of output variable
        CHARACTER(LEN=30)  ::  nc_filename        !<
    
    
        INTEGER(iwp), PARAMETER ::  n_turbines_max = 1E4  !< maximum number of turbines (for array allocation)
    
    
    !
    !-- Variables specified in the namelist wind_turbine_par
        INTEGER(iwp) ::  n_airfoils = 8  !< number of airfoils of the used turbine model (for ADM-R and ALM)
        INTEGER(iwp) ::  n_turbines = 1  !< number of turbines
    
        LOGICAL ::  pitch_control       = .FALSE.  !< switch for use of pitch controller
        LOGICAL ::  speed_control       = .FALSE.  !< switch for use of speed controller
        LOGICAL ::  tip_loss_correction = .FALSE.  !< switch for use of tip loss correct.
        LOGICAL ::  yaw_control         = .FALSE.  !< switch for use of yaw controller
    
    
        LOGICAL ::  initial_write_coordinates = .FALSE.  !<
    
    
        REAL(wp), DIMENSION(:), POINTER   ::  output_values_1d_pointer  !< pointer for 2d output array
        REAL(wp), POINTER                 ::  output_values_0d_pointer  !< pointer for 2d output array
        REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET ::  output_values_1d_target  !< pointer for 2d output array
        REAL(wp), TARGET                  ::  output_values_0d_target   !< pointer for 2d output array
    
    
        REAL(wp) ::  dt_data_output_wtm = 0.0_wp  !< data output interval
        REAL(wp) ::  time_wtm           = 0.0_wp  !< time since last data output
    
    
        REAL(wp) ::  segment_length_tangential  = 1.0_wp  !< length of the segments, the rotor area is divided into
                                                          !< (in tangential direction, as factor of MIN(dx,dy,dz))
        REAL(wp) ::  segment_width_radial       = 0.5_wp  !< width of the segments, the rotor area is divided into
                                                          !< (in radial direction, as factor of MIN(dx,dy,dz))
    
        REAL(wp) ::  smearing_kernel_size       = 2.0_wp  !< size of the smearing kernel as multiples of dx
    
        REAL(wp) ::  time_turbine_on            = 0.0_wp  !< time at which turbines are started
        REAL(wp) ::  tilt_angle                 = 0.0_wp  !< vertical tilt_angle of the rotor [degree] ( positive = backwards )
    
                                                                                     !< ( clockwise, 0 = facing west )
        REAL(wp), DIMENSION(1:n_turbines_max) ::  hub_x               = 9999999.9_wp !< position of hub in x-direction
        REAL(wp), DIMENSION(1:n_turbines_max) ::  hub_y               = 9999999.9_wp !< position of hub in y-direction
        REAL(wp), DIMENSION(1:n_turbines_max) ::  hub_z               = 9999999.9_wp !< position of hub in z-direction
        REAL(wp), DIMENSION(1:n_turbines_max) ::  nacelle_cd          = 0.85_wp      !< drag coefficient for nacelle
        REAL(wp), DIMENSION(1:n_turbines_max) ::  nacelle_radius      = 0.0_wp       !< nacelle radius [m]
        REAL(wp), DIMENSION(1:n_turbines_max) ::  tower_cd            = 1.2_wp       !< drag coefficient for tower
        REAL(wp), DIMENSION(1:n_turbines_max) ::  pitch_angle         = 0.0_wp       !< constant pitch angle
        REAL(wp), DIMENSION(1:n_turbines_max) ::  rotor_radius        = 63.0_wp      !< rotor radius [m]
        REAL(wp), DIMENSION(1:n_turbines_max), TARGET ::  rotor_speed = 0.9_wp       !< inital or constant rotor speed [rad/s]
        REAL(wp), DIMENSION(1:n_turbines_max) ::  tower_diameter      = 0.0_wp       !< tower diameter [m]
        REAL(wp), DIMENSION(1:n_turbines_max) ::  yaw_angle           = 0.0_wp       !< yaw angle [degree]
    
    
    !
    !-- Variables specified in the namelist for speed controller
    !-- Default values are from the NREL 5MW research turbine (Jonkman, 2008)
        REAL(wp) ::  air_density               = 1.225_wp       !< Air density to convert to W [kg/m3]
        REAL(wp) ::  gear_efficiency           = 1.0_wp         !< Loss between rotor and generator
        REAL(wp) ::  gear_ratio                = 97.0_wp        !< Gear ratio from rotor to generator
        REAL(wp) ::  generator_power_rated     = 5296610.0_wp   !< rated turbine power [W]
        REAL(wp) ::  generator_inertia         = 534.116_wp     !< Inertia of the generator [kg*m2]
        REAL(wp) ::  generator_efficiency      = 0.944_wp       !< Electric efficiency of the generator
        REAL(wp) ::  generator_speed_rated     = 121.6805_wp    !< Rated generator speed [rad/s]
        REAL(wp) ::  generator_torque_max      = 47402.91_wp    !< Maximum of the generator torque [Nm]
        REAL(wp) ::  generator_torque_rate_max = 15000.0_wp     !< Max generator torque increase [Nm/s]
        REAL(wp) ::  pitch_rate                = 8.0_wp         !< Max pitch rate [degree/s]
        REAL(wp) ::  region_2_slope            = 2.332287_wp    !< Slope constant for region 2
        REAL(wp) ::  region_2_min              = 91.21091_wp    !< Lower generator speed boundary of region 2 [rad/s]
        REAL(wp) ::  region_15_min             = 70.16224_wp    !< Lower generator speed boundary of region 1.5 [rad/s]
        REAL(wp) ::  rotor_inertia             = 34784179.0_wp  !< Inertia of the rotor [kg*m2]
    
    
    
    !
    !-- Variables specified in the namelist for yaw control:
        REAL(wp) ::  yaw_misalignment_max = 0.08726_wp   !< maximum tolerated yaw missalignment [rad]
        REAL(wp) ::  yaw_misalignment_min = 0.008726_wp  !< minimum yaw missalignment for which the actuator stops [rad]
        REAL(wp) ::  yaw_speed            = 0.005236_wp  !< speed of the yaw actuator [rad/s]
    
    !
    !-- Variables for initialization of the turbine model:
        INTEGER(iwp) ::  inot        !< turbine loop index (turbine id)
        INTEGER(iwp) ::  nsegs_max   !< maximum number of segments (all turbines, required for allocation of arrays)
        INTEGER(iwp) ::  nrings_max  !< maximum number of rings (all turbines, required for allocation of arrays)
        INTEGER(iwp) ::  ring        !< ring loop index (ring number)
        INTEGER(iwp) ::  upper_end   !<
    
        INTEGER(iwp), DIMENSION(1) ::  lct  !<
    
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_hub    !< index belonging to x-position of the turbine
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_smear  !< index defining the area for the smearing of the forces (x-direction)
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_hub    !< index belonging to y-position of the turbine
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_smear  !< index defining the area for the smearing of the forces (y-direction)
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_hub    !< index belonging to hub height
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_smear  !< index defining the area for the smearing of the forces (z-direction)
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nrings   !< number of rings per turbine
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nsegs_total  !< total number of segments per turbine
    
        INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nsegs  !< number of segments per ring and turbine
    
    !
    !-  Parameters for the smearing from the rotor to the cartesian grid:
        REAL(wp) ::  delta_t_factor  !<
        REAL(wp) ::  eps_factor      !<
        REAL(wp) ::  eps_min         !<
        REAL(wp) ::  eps_min2        !<
        REAL(wp) ::  pol_a           !< parameter for the polynomial smearing fct
        REAL(wp) ::  pol_b           !< parameter for the polynomial smearing fct
    
    !
    !-- Variables for the calculation of lift and drag coefficients:
        REAL(wp), DIMENSION(:), ALLOCATABLE  ::  ard      !<
        REAL(wp), DIMENSION(:), ALLOCATABLE  ::  crd      !<
        REAL(wp), DIMENSION(:), ALLOCATABLE  ::  delta_r  !< radial segment length
        REAL(wp), DIMENSION(:), ALLOCATABLE  ::  lrd      !<
    
        REAL(wp) ::  accu_cl_cd_tab = 0.1_wp  !< Accuracy of the interpolation of the lift and drag coeff [deg]
    
        REAL(wp), DIMENSION(:,:), ALLOCATABLE :: turb_cd_tab  !< table of the blade drag coefficient
        REAL(wp), DIMENSION(:,:), ALLOCATABLE :: turb_cl_tab  !< table of the blade lift coefficient
    
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  nac_cd_surf  !< 3d field of the nacelle drag coefficient
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tow_cd_surf  !< 3d field of the tower drag coefficient
    
    !
    !-- Variables for the calculation of the forces:
        REAL(wp) ::  cur_r       !<
        REAL(wp) ::  phi_rotor   !<
        REAL(wp) ::  pre_factor  !<
        REAL(wp) ::  torque_seg  !<
        REAL(wp) ::  u_int_l     !<
        REAL(wp) ::  u_int_u     !<
        REAL(wp) ::  u_rot       !<
        REAL(wp) ::  v_int_l     !<
        REAL(wp) ::  v_int_u     !<
        REAL(wp) ::  w_int_l     !<
        REAL(wp) ::  w_int_u     !<
        !$OMP THREADPRIVATE (cur_r, phi_rotor, pre_factor, torque_seg, u_int_l, u_int_u, u_rot, &
        !$OMP&               v_int_l, v_int_u, w_int_l, w_int_u)
    !
    !-  Tendencies from the nacelle and tower thrust:
        REAL(wp) ::  tend_nac_x = 0.0_wp  !<
        REAL(wp) ::  tend_nac_y = 0.0_wp  !<
        REAL(wp) ::  tend_tow_x = 0.0_wp  !<
        REAL(wp) ::  tend_tow_y = 0.0_wp  !<
        !$OMP THREADPRIVATE (tend_nac_x, tend_tow_x, tend_nac_y, tend_tow_y)
    
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  alpha_attack  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  chord         !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  phi_rel       !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  torque_total  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  thrust_rotor  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  turb_cd       !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  turb_cl       !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  vrel          !<
        REAL(wp), DIMENSION(:), ALLOCATABLE ::  vtheta        !<
    
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rbx, rby, rbz     !< coordinates of the blade elements
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rotx, roty, rotz  !< normal vectors to the rotor coordinates
    
    !
    !-  Fields for the interpolation of velocities on the rotor grid:
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_int      !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_int_1_l  !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_int      !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_int_1_l  !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_int      !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_int_1_l  !<
    
    !
    !-  Rotor tendencies on the segments:
        REAL(wp), DIMENSION(:), ALLOCATABLE :: thrust_seg    !<
        REAL(wp), DIMENSION(:), ALLOCATABLE :: torque_seg_y  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE :: torque_seg_z  !<
    
    !
    !-  Rotor tendencies on the rings:
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  thrust_ring    !<
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  torque_ring_y  !<
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  torque_ring_z  !<
    
    !
    !-  Rotor tendencies on rotor grids for all turbines:
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  thrust     !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  torque_y   !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  torque_z   !<
    
    !
    !-  Rotor tendencies on coordinate grid:
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_x  !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_y  !<
        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_z  !<
    
    !
    !-  Variables for the rotation of the rotor coordinates:
        REAL(wp), DIMENSION(1:n_turbines_max,1:3,1:3) ::  rot_coord_trans  !< matrix for rotation of rotor coordinates
    
        REAL(wp), DIMENSION(1:3) ::  rot_eigen_rad  !<
        REAL(wp), DIMENSION(1:3) ::  rot_eigen_azi  !<
        REAL(wp), DIMENSION(1:3) ::  rot_eigen_nor  !<
        REAL(wp), DIMENSION(1:3) ::  re             !<
        REAL(wp), DIMENSION(1:3) ::  rea            !<
        REAL(wp), DIMENSION(1:3) ::  ren            !<
        REAL(wp), DIMENSION(1:3) ::  rote           !<
        REAL(wp), DIMENSION(1:3) ::  rota           !<
        REAL(wp), DIMENSION(1:3) ::  rotn           !<
    
    !
    !-- Fixed variables for the speed controller:
        LOGICAL  ::  start_up = .TRUE.  !<
    
        REAL(wp) ::  fcorner           !< corner freq for the controller low pass filter
        REAL(wp) ::  om_rate            !< rotor speed change
        REAL(wp) ::  region_25_min      !< min region 2.5
        REAL(wp) ::  region_25_slope    !< slope in region 2.5
        REAL(wp) ::  slope15            !< slope in region 1.5
        REAL(wp) ::  trq_rate           !< torque change
        REAL(wp) ::  vs_sysp            !<
        REAL(wp) ::  lp_coeff           !< coeff for the controller low pass filter
    
        REAL(wp), DIMENSION(n_turbines_max) :: rotor_speed_l = 0.0_wp  !< local rot speed [rad/s]
    
    !
    !-- Fixed variables for the yaw controller:
        INTEGER(iwp)                          ::  wdlon           !<
        INTEGER(iwp)                          ::  wdsho            !<
    
        LOGICAL,  DIMENSION(1:n_turbines_max) ::  doyaw = .FALSE.  !<
    
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  u_inflow         !< wind speed at hub
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  u_inflow_l       !<
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wdir             !< wind direction at hub
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wdir_l           !<
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wd2_l            !< local (cpu) short running avg of the wd
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wd30_l           !< local (cpu) long running avg of the wd
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  yawdir           !< direction to yaw
        REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  yaw_angle_l      !< local (cpu) yaw angle
    
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wd2              !<
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wd30             !<
    
    
    !
    !-- Variables that have to be saved in the binary file for restarts:
        REAL(wp), DIMENSION(1:n_turbines_max) ::  pitch_angle_old         = 0.0_wp  !< old constant pitch angle
        REAL(wp), DIMENSION(1:n_turbines_max) ::  generator_speed         = 0.0_wp  !< curr. generator speed
        REAL(wp), DIMENSION(1:n_turbines_max) ::  generator_speed_f       = 0.0_wp  !< filtered generator speed
        REAL(wp), DIMENSION(1:n_turbines_max) ::  generator_speed_old     = 0.0_wp  !< last generator speed
        REAL(wp), DIMENSION(1:n_turbines_max) ::  generator_speed_f_old   = 0.0_wp  !< last filtered generator speed
        REAL(wp), DIMENSION(1:n_turbines_max) ::  torque_gen              = 0.0_wp  !< generator torque
        REAL(wp), DIMENSION(1:n_turbines_max) ::  torque_gen_old          = 0.0_wp  !< last generator torque
    
    
        SAVE
    
    
        INTERFACE wtm_parin
           MODULE PROCEDURE wtm_parin
        END INTERFACE wtm_parin
    
        INTERFACE wtm_check_parameters
           MODULE PROCEDURE wtm_check_parameters
        END INTERFACE wtm_check_parameters
    
        INTERFACE wtm_data_output
           MODULE PROCEDURE wtm_data_output
        END INTERFACE wtm_data_output
    
        INTERFACE wtm_init_arrays
           MODULE PROCEDURE wtm_init_arrays
        END INTERFACE wtm_init_arrays
    
        INTERFACE wtm_init
           MODULE PROCEDURE wtm_init
        END INTERFACE wtm_init
    
        INTERFACE wtm_init_output
           MODULE PROCEDURE wtm_init_output
        END INTERFACE wtm_init_output
    
        INTERFACE wtm_actions
           MODULE PROCEDURE wtm_actions
           MODULE PROCEDURE wtm_actions_ij
        END INTERFACE wtm_actions
    
        INTERFACE wtm_rrd_global
           MODULE PROCEDURE wtm_rrd_global_ftn
           MODULE PROCEDURE wtm_rrd_global_mpi
        END INTERFACE wtm_rrd_global
    
        INTERFACE wtm_wrd_global
           MODULE PROCEDURE wtm_wrd_global
        END INTERFACE wtm_wrd_global
    
    
        PUBLIC                                                                                         &
               dt_data_output_wtm,                                                                     &
               time_wtm,                                                                               &
               wind_turbine
    
        PUBLIC                                                                                         &
               wtm_parin,                                                                              &
               wtm_check_parameters,                                                                   &
               wtm_data_output,                                                                        &
               wtm_init_arrays,                                                                        &
               wtm_init_output,                                                                        &
               wtm_init,                                                                               &
               wtm_actions,                                                                            &
               wtm_rrd_global,                                                                         &
               wtm_wrd_global
    
    
     CONTAINS


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Parin for &anthropogenic_heat_par for external anthropogenic heat sources model
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE ah_parin
    
      IMPLICIT NONE
  
      CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file
  
      INTEGER(iwp) ::  io_status   !< status after reading the namelist file
  
      LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                               !< although the respective module namelist appears in
                                               !< the namelist file
  
      NAMELIST /anthropogenic_heat_parameters/  switch_off_module                                  
  
  !
  !-- Move to the beginning of the namelist file and try to find and read the namelist.
      REWIND( 11 )
      READ( 11, anthropogenic_heat_parameters, IOSTAT=io_status )
  
  !
  !-- Action depending on the READ status
      IF ( io_status == 0 )  THEN
  !
  !--    anthropogenic_heat_parameters namelist was found and read correctly. Enable the
  !--    external anthropogenic heat module.
         IF ( .NOT. switch_off_module )  external_anthropogenic_heat = .TRUE.
  
      ELSEIF ( io_status > 0 )  THEN
  !
  !--    anthropogenic_heat_parameters namelist was found but contained errors. Print an error message
  !--    including the line that caused the problem.
         BACKSPACE( 11 )
         READ( 11 , '(A)' ) line
         CALL parin_fail_message( 'anthropogenic_heat_turbines', line )
  
      ENDIF
  
   END SUBROUTINE ah_parin

    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Parin for &wind_turbine_par for wind turbine model
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_parin
    
        IMPLICIT NONE
    
        CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file
    
        INTEGER(iwp) ::  io_status   !< status after reading the namelist file
    
        LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                                 !< although the respective module namelist appears in
                                                 !< the namelist file
    
        NAMELIST /wind_turbine_parameters/  air_density,                                               &
                                            dt_data_output_wtm,                                        &
                                            gear_efficiency,                                           &
                                            gear_ratio,                                                &
                                            generator_efficiency,                                      &
                                            generator_inertia,                                         &
                                            generator_power_rated,                                     &
                                            generator_speed_rated,                                     &
                                            generator_torque_max,                                      &
                                            generator_torque_rate_max,                                 &
                                            hub_x,                                                     &
                                            hub_y,                                                     &
                                            hub_z,                                                     &
                                            nacelle_cd,                                                &
                                            nacelle_radius,                                            &
                                            n_airfoils,                                                &
                                            n_turbines,                                                &
                                            pitch_angle,                                               &
                                            pitch_control,                                             &
                                            pitch_rate,                                                &
                                            region_15_min,                                             &
                                            region_2_min,                                              &
                                            region_2_slope,                                            &
                                            rotor_inertia,                                             &
                                            rotor_radius,                                              &
                                            rotor_speed,                                               &
                                            segment_length_tangential,                                 &
                                            segment_width_radial,                                      &
                                            smearing_kernel_size,                                      &
                                            speed_control,                                             &
                                            switch_off_module,                                         &
                                            tilt_angle,                                                &
                                            time_turbine_on,                                           &
                                            tip_loss_correction,                                       &
                                            tower_cd,                                                  &
                                            tower_diameter,                                            &
                                            yaw_angle,                                                 &
                                            yaw_control,                                               &
                                            yaw_misalignment_max,                                      &
                                            yaw_misalignment_min,                                      &
                                            yaw_speed
    
    !
    !-- Move to the beginning of the namelist file and try to find and read the namelist.
        REWIND( 11 )
        READ( 11, wind_turbine_parameters, IOSTAT=io_status )
    
    !
    !-- Action depending on the READ status
        IF ( io_status == 0 )  THEN
    !
    !--    wind_turbine_parameters namelist was found and read correctly. Enable the
    !--    wind turbine module.
           IF ( .NOT. switch_off_module )  wind_turbine = .TRUE.
    
        ELSEIF ( io_status > 0 )  THEN
    !
    !--    wind_turbine_parameters namelist was found but contained errors. Print an error message
    !--    including the line that caused the problem.
           BACKSPACE( 11 )
           READ( 11 , '(A)' ) line
           CALL parin_fail_message( 'wind_turbine_parameters', line )
    
        ENDIF
    
     END SUBROUTINE wtm_parin
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine writes the respective restart data.
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_wrd_global
    
    
        IMPLICIT NONE
    
    
       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
    
           CALL wrd_write_string( 'generator_speed' )
           WRITE ( 14 )  generator_speed
    
           CALL wrd_write_string( 'generator_speed_f' )
           WRITE ( 14 )  generator_speed_f
    
           CALL wrd_write_string( 'generator_speed_f_old' )
           WRITE ( 14 )  generator_speed_f_old
    
           CALL wrd_write_string( 'generator_speed_old' )
           WRITE ( 14 )  generator_speed_old
    
           CALL wrd_write_string( 'rotor_speed' )
           WRITE ( 14 )  rotor_speed
    
           CALL wrd_write_string( 'yaw_angle' )
           WRITE ( 14 )  yaw_angle
    
           CALL wrd_write_string( 'pitch_angle' )
           WRITE ( 14 )  pitch_angle
    
           CALL wrd_write_string( 'pitch_angle_old' )
           WRITE ( 14 )  pitch_angle_old
    
           CALL wrd_write_string( 'torque_gen' )
           WRITE ( 14 )  torque_gen
    
           CALL wrd_write_string( 'torque_gen_old' )
           WRITE ( 14 )  torque_gen_old
    
        ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
    
           CALL wrd_mpi_io_global_array( 'generator_speed', generator_speed )
           CALL wrd_mpi_io_global_array( 'generator_speed_f', generator_speed_f )
           CALL wrd_mpi_io_global_array( 'generator_speed_f_old', generator_speed_f_old )
           CALL wrd_mpi_io_global_array( 'generator_speed_old', generator_speed_old )
           CALL wrd_mpi_io_global_array( 'rotor_speed', rotor_speed )
           CALL wrd_mpi_io_global_array( 'yaw_angle', yaw_angle )
           CALL wrd_mpi_io_global_array( 'pitch_angle', pitch_angle )
           CALL wrd_mpi_io_global_array( 'pitch_angle_old', pitch_angle_old )
           CALL wrd_mpi_io_global_array( 'torque_gen', torque_gen )
           CALL wrd_mpi_io_global_array( 'torque_gen_old', torque_gen_old )
    
        ENDIF
    
     END SUBROUTINE wtm_wrd_global
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Read module-specific global restart data (Fortran binary format).
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_rrd_global_ftn( found )
    
    
        USE control_parameters,                                                                        &
            ONLY:  length,                                                                             &
                   restart_string
    
        IMPLICIT NONE
    
        LOGICAL, INTENT(OUT) ::  found
    
    
        found = .TRUE.
    
        SELECT CASE ( restart_string(1:length) )
    
           CASE ( 'generator_speed' )
              READ ( 13 )  generator_speed
           CASE ( 'generator_speed_f' )
              READ ( 13 )  generator_speed_f
           CASE ( 'generator_speed_f_old' )
              READ ( 13 )  generator_speed_f_old
           CASE ( 'generator_speed_old' )
              READ ( 13 )  generator_speed_old
           CASE ( 'rotor_speed' )
              READ ( 13 )  rotor_speed
           CASE ( 'yaw_angle' )
              READ ( 13 )  yaw_angle
           CASE ( 'pitch_angle' )
              READ ( 13 )  pitch_angle
           CASE ( 'pitch_angle_old' )
              READ ( 13 )  pitch_angle_old
           CASE ( 'torque_gen' )
              READ ( 13 )  torque_gen
           CASE ( 'torque_gen_old' )
              READ ( 13 )  torque_gen_old
    
           CASE DEFAULT
    
              found = .FALSE.
    
        END SELECT
    
    
     END SUBROUTINE wtm_rrd_global_ftn
    
    
    !------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Read module-specific global restart data (MPI-IO).
    !------------------------------------------------------------------------------!
     SUBROUTINE wtm_rrd_global_mpi
    
        CALL rrd_mpi_io_global_array( 'generator_speed', generator_speed )
        CALL rrd_mpi_io_global_array( 'generator_speed_f', generator_speed_f )
        CALL rrd_mpi_io_global_array( 'generator_speed_f_old', generator_speed_f_old )
        CALL rrd_mpi_io_global_array( 'generator_speed_old', generator_speed_old )
        CALL rrd_mpi_io_global_array( 'rotor_speed', rotor_speed )
        CALL rrd_mpi_io_global_array( 'yaw_angle', yaw_angle )
        CALL rrd_mpi_io_global_array( 'pitch_angle', pitch_angle )
        CALL rrd_mpi_io_global_array( 'pitch_angle_old', pitch_angle_old )
        CALL rrd_mpi_io_global_array( 'torque_gen', torque_gen )
        CALL rrd_mpi_io_global_array( 'torque_gen_old', torque_gen_old )
    
     END SUBROUTINE wtm_rrd_global_mpi
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Check namelist parameter
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_check_parameters
    
        IMPLICIT NONE
    
        IF ( .NOT. input_pids_wtm )  THEN
           IF ( ( .NOT.speed_control ) .AND. pitch_control )  THEN
              message_string = 'pitch_control = .TRUE. requires speed_control = .TRUE.'
              CALL message( 'wtm_check_parameters', 'WTM0001', 1, 2, 0, 6, 0 )
           ENDIF
    
           IF ( ANY( rotor_speed(1:n_turbines) < 0.0 ) )  THEN
              message_string = 'rotor_speed < 0.0'
              CALL message( 'wtm_check_parameters', 'WTM0002', 1, 2, 0, 6, 0 )
           ENDIF
    
           IF ( ANY( hub_x(1:n_turbines) == 9999999.9_wp ) .OR.                                        &
                ANY( hub_y(1:n_turbines) == 9999999.9_wp ) .OR.                                        &
                ANY( hub_z(1:n_turbines) == 9999999.9_wp ) )  THEN
    
              message_string = 'hub_x, hub_y, hub_z have to be given for each turbine'
              CALL message( 'wtm_check_parameters', 'WTM0003', 1, 2, 0, 6, 0 )
           ENDIF
        ENDIF
    
     END SUBROUTINE wtm_check_parameters

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads anthropogenic heat profiles emitted from building and ground surfaces from a NetCDF file.
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_anthro_heat_profiles
        
       IMPLICIT NONE

       INTEGER(iwp) ::  id_netcdf                                    !< NetCDF id of input file
       INTEGER(iwp) ::  num_vars                                     !< number of variables in input file
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names   !< variable names in static input file   
       
       
       INTEGER(iwp) ::  n_buildings = 0           !< number of buildings (for array allocation)
       INTEGER(iwp) ::  n_streets = 0             !< number of streets (for array allocation)
       INTEGER(iwp) ::  n_points = 0              !< number of point sources (for array allocation)
       INTEGER(iwp) ::  n_timesteps = 0           !< number of time steps (for array allocation)
        
    !-- If no anthropogenic heat input file is available, skip this routine
       IF ( .NOT. input_pids_ah )  RETURN
    !
    !-- Measure CPU time
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'start' )
    !
    !-- Skip the following if no expternal anthropogenic heat profiles are to be considered in the calculations.
       IF ( .NOT. external_anthropogenic_heat )  RETURN
        
    #if defined ( __netcdf )
    !
    !-- Open file in read-only mode
       CALL open_read_file( TRIM( input_file_ah ) // TRIM( coupling_char ) , id_netcdf )
    !
    !-- Inquire all variable names.
    !-- This will be used to check whether an optional input variable exists or not.
       CALL inquire_num_variables( id_netcdf, num_vars )
        
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_netcdf, var_names )
    !
    !-- Read dimensions from file
       CALL get_dimension_length( id_netcdf, n_buildings, 'building_id' )
       CALL get_dimension_length( id_netcdf, n_streets, 'street_id' )
       CALL get_dimension_length( id_netcdf, n_points, 'point_id' )
       CALL get_dimension_length( id_netcdf, n_timesteps, 'time' )
    !
    !-- Read building ids from file 
       IF ( n_buildings > 0 .AND. check_existence( var_names, 'building_id' ) )  THEN
          building_ids%from_file = .TRUE.
        
          ALLOCATE( building_ids%val(1:n_buildings) ) 
        
          CALL get_variable( id_netcdf, 'building_id', building_ids%val )
       ELSE
          building_ids%from_file = .FALSE.
       ENDIF
    !
    !-- Read streed ids from file 
       IF ( n_streets > 0 .AND. check_existence( var_names, 'street_id' ) )  THEN
          street_ids%from_file = .TRUE.
      
          ALLOCATE( street_ids%val(1:n_streets) ) 
      
          CALL get_variable( id_netcdf, 'street_id', street_ids%val )
       ELSE
          street_ids%from_file = .FALSE.
       ENDIF
    !
    !-- Read point ids from file 
       IF ( n_points > 0 .AND. check_existence( var_names, 'point_id' ) )  THEN
         point_ids%from_file = .TRUE.
     
         ALLOCATE( point_ids%val(1:n_points) ) 
     
         CALL get_variable( id_netcdf, 'point_id', point_ids%val )
      ELSE
         point_ids%from_file = .FALSE.
      ENDIF
    !
    !-- Read timesteps from file
       IF ( n_timesteps > 0 .AND. check_existence( var_names, 'time' ) .AND. 
           ( building_ids%from_file .OR. street_ids%from_file .OR. point_ids%from_file ) )  THEN   
          ALLOCATE( t(1:n_timesteps) ) 
         
          CALL get_variable( id_netcdf, 'time', t)
       ENDIF

    !
    !-- Read anthrpogenic heat profiles from buildings from file
       IF ( building_ids%from_file .AND. check_existence( var_names, 'building_ah') ) THEN
          building_ah%from_file = .TRUE.
          CALL get_attribute( id_netcdf, char_fill, building_ah%fill, .FALSE., 'building_ah', .FALSE. )
        
          ALLOCATE( building_ah%val(1:n_buildings,1:n_timesteps) )
        
          CALL get_variable( id_netcdf, 'building_ah', building_ah%val )
          !-- TODO:(DONE) Output warning if the anthropogenic heat profile is not fully defined (use message.f90 module, use message_string from control_parameters under modules.f90)
          CALL ah_check_input_profiles('buildings', building_ah%val, building_ah%fill)
       ELSE
          building_ah%from_file = .FALSE.
       ENDIF
    !
    !-- Read anthrpogenic heat profiles from streets from file
       IF ( street_ids%from_file .AND. check_existence( var_names, 'street_ah') ) THEN
          street_ah%from_file = .TRUE.
          CALL get_attribute( id_netcdf, char_fill, street_ah%fill, .FALSE., 'street_ah', .FALSE. )
       
          ALLOCATE( street_ah%val(1:n_buildings,1:n_timesteps) )
       
          CALL get_variable( id_netcdf, 'street_ah', street_ah%val )
          !-- TODO:(DONE) Output warning if the anthropogenic heat profile is not fully defined (use message.f90 module, use message_string from control_parameters under modules.f90)
          CALL ah_check_input_profiles('streets', street_ah%val, street_ah%fill)
       ELSE
          street_ah%from_file = .FALSE.
       ENDIF
    !
    !-- Read anthropogenic heat profiles from points from file
       IF ( point_ids%from_file .AND. check_existence( var_names, 'point_ah') ) THEN
          point_ah%from_file = .TRUE.
          CALL get_attribute( id_netcdf, char_fill, point_ah%fill, .FALSE., 'point_ah', .FALSE. )
       
          ALLOCATE( point_ah%val(1:n_buildings,1:n_timesteps) )
       
          CALL get_variable( id_netcdf, 'point_ah', point_ah%val )
          !-- TODO:(DONE) Output warning if the anthropogenic heat profile is not fully defined (use message.f90 module)
          CALL ah_check_input_profiles('points', point_ah%val, point_ah%fill)
       ELSE
          point_ah%from_file = .FALSE.
       ENDIF
   
    !
    !-- Read coordinates of point sources from file
       IF ( point_ah%from_file .AND. check_existence( var_names, 'point_x') .AND. check_existence( var_names, 'point_y') ) THEN
         point_coords%from_file = .TRUE.
         CALL get_attribute( id_netcdf, char_fill, point_coords%x_fill, .FALSE., 'point_x', .FALSE. )
         CALL get_attribute( id_netcdf, char_fill, point_coords%y_fill, .FALSE., 'point_y', .FALSE. )
      
         ALLOCATE( point_coords%x(1:n_points) )
         ALLOCATE( point_coords%y(1:n_points) )
      
         CALL get_variable( id_netcdf, 'point_x', point_coords%x )
         CALL get_variable( id_netcdf, 'point_y', point_coords%y )

         !-- TODO: Implement error handling if point coorinates are not fully defined (use message.f90 module)
         
      ELSE
         point_coords%from_file = .FALSE.
      ENDIF
    
    !
    !-- Finally, close input file
       CALL close_input_file( id_netcdf )
    #endif
    !
    !-- End of CPU measurement
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'stop' )
    
    END SUBROUTINE netcdf_data_input_anthro_heat_profiles
    

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Check anthropogenic heat profiles from input
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_check_input_profiles(source, profile, fill_value)

       IMPLICIT NONE

       CHARACTER(len=20), INTENT(IN) :: source           !< source of the anthropogenic heat profile
       REAL(wp), DIMENSION(:,:), INTENT(IN) :: profile   !< anthropogenic heat profile
       REAL(wp), INTENT(IN) :: fill_value                !< fill value of the anthropogenic heat profile

       !-- Check if the anthropogenic heat profile contains only non-negative values
       IF ( ANY( ( profile < 0.0_wp ) .AND. ( profile > fill_value )) ) THEN

         message_string = 'Anthropogenic heat profiles from ' // source // ' must be non-negative.'
         CALL message( 'netcdf_data_input_anthro_heat_profiles', 'WTM0001', 1, 2, 0, 6, 0 )
       
       !-- Check if all values of the input anthropogenic heat profile have been defined
       ELSE IF ( ANY( profile == fill_value ) ) THEN

         message_string = 'Some timesteps in the anthropogenic heat profile from ' // source // ' are not fully defined.'
         CALL message( 'netcdf_data_input_anthro_heat_profiles', 'WTM0001', 0, 1, 0, 6, 0 ) 
         
       ENDIF 

    END SUBROUTINE ah_check_input_profiles


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Check anthropogenic heat profiles from input
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_check_point_source_locations(point_coordinates)

       IMPLICIT NONE

       TYPE(point_coords_t), INTENT(IN) :: point_coordinates  !< point source coordinates

       !-- Check if the point source coordinates contain only non-negative values


    END SUBROUTINE ah_check_point_source_locations
    
    !
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine applies the imported anthropogenic heat profiles to the the 
    !> corresponding urban surfaces.
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE apply_anthro_heat_profiles_to_surfaces

       IMPLICIT NONE

    END SUBROUTINE apply_anthro_heat_profiles_to_surfaces



    !
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine fetches the roof surface tiles based on building ids
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_building_id_to_surfaces(building_id, surf_indexes)

       USE indoor_model_mod, ONLY: buildings

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: building_id                       !< building id
       INTEGER(iwp), DIMENSION(:), INTENT(OUT) :: surf_indexes       !< surface index

       DO b = lbound(buildings), ubound(buildings)
          IF ( buildings(b)%id == building_id ) THEN
             ALLOCATE( surf_indexes(lbound(buildings(b)%m):ubound(buildings(b)%m)) )
             surf_indexes = buildings(b)%m
             EXIT
          ENDIF
       END DO

    END SUBROUTINE ah_building_id_to_surfaces


    !
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine fetches the street surface tiles based on street ids
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_street_id_to_surfaces(street_id)

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: street_id  !< street id

       !TODO: Implement body of this routine (look in chem_modules.f90
       !      line 1297 ff. chem_emissions_mod.f90)

    END SUBROUTINE ah_street_id_to_surfaces


    !
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine fetches the ground surface tiles based on point-source ids
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_point_id_to_surfaces(point_id, surf_index)

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: point_id     !< point id
       INTEGER(iwp), INTENT(OUT) :: surf_index  !< surface index

       IF ( point_coords%from_file ) THEN
          !-- Retrieve UTM coordinates of point sources from file
          DO p = 1, n_points
             IF ( point_coords%val(p) == point_id ) THEN
                x_coord_abs = point_coords%x(p)
                y_coord_abs = point_coords%y(p)
                EXIT
             ENDIF
          END DO
         
          !
          ! -- offset and rotate UTM coordinates to match the model grid
          x_coord_rel(t) = x_coord_abs(t) - init_model%origin_x
          y_coord_rel(t) = y_coord_abs(t) - init_model%origin_y
          x_coord(t) = COS( init_model%rotation_angle * pi / 180.0_wp ) * x_coord_rel(t)            &
                     - SIN( init_model%rotation_angle * pi / 180.0_wp ) * y_coord_rel(t)
          y_coord(t) = SIN( init_model%rotation_angle * pi / 180.0_wp ) * x_coord_rel(t)            &
                     + COS( init_model%rotation_angle * pi / 180.0_wp ) * y_coord_rel(t)

          !
          ! -- Then, compute the indices of the grid cell in which the point source is located.
          is = INT( x_coord(t) * ddx, KIND = iwp )
          js = INT( y_coord(t) * ddy, KIND = iwp )

          ! TODO:(DONE) Finish writing comments
          ! -- Find the surface index of the grid cell in which the point source is located.
          ! -- First, search the land surfaces (lsm) for a match with the coordinates of the point source.
          DO m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             IF ( i == is .AND. j == js ) THEN
                surf_index = m
                found = .TRUE.
                EXIT
             ENDIF
          END DO

          ! -- If no match was found in the land surfaces, search the urban surfaces (usm).
          IF (.NOT. found) THEN
             DO m = 1, surf_usm%ns   ! TODO: Limit this looping to horizontal surfaces
                i = surf_usm%i(m)
                j = surf_usm%j(m)
                IF ( i == is .AND. j == js ) THEN
                   surf_index = m
                   found = .TRUE.
                   EXIT
                ENDIF
             END DO
          ENDIF
         
       ENDIF

    END SUBROUTINE ah_point_id_to_surfaces


    !
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Allocate wind turbine model arrays
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_init_arrays
    
        IMPLICIT NONE
    
        REAL(wp) ::  delta_r_factor  !<
        REAL(wp) ::  delta_r_init    !<
    
    #if defined( __netcdf )
    !
    ! Read wtm input file (netcdf) if it exists:
        IF ( input_pids_wtm )  THEN
    
    !
    !--    Open the wtm  input file:
           CALL open_read_file( TRIM( input_file_wtm ) //                                              &
                                TRIM( coupling_char ), pids_id )
    
           CALL inquire_num_variables( pids_id, num_var_pids )
    
    !
    !--    Allocate memory to store variable names and read them:
           ALLOCATE( vars_pids(1:num_var_pids) )
           CALL inquire_variable_names( pids_id, vars_pids )
    
    !
    !--    Input of all wtm parameters:
           IF ( check_existence( vars_pids, 'tower_diameter' ) )  THEN
              CALL get_variable( pids_id, 'tower_diameter', tower_diameter(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'rotor_speed' ) )  THEN
              CALL get_variable( pids_id, 'rotor_speed', rotor_speed(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'pitch_angle' ) )  THEN
              CALL get_variable( pids_id, 'pitch_angle', pitch_angle(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'yaw_angle' ) )  THEN
              CALL get_variable( pids_id, 'yaw_angle', yaw_angle(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'hub_x' ) )  THEN
              CALL get_variable( pids_id, 'hub_x', hub_x(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'hub_y' ) )  THEN
              CALL get_variable( pids_id, 'hub_y', hub_y(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'hub_z' ) )  THEN
              CALL get_variable( pids_id, 'hub_z', hub_z(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'nacelle_radius' ) )  THEN
              CALL get_variable( pids_id, 'nacelle_radius', nacelle_radius(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'rotor_radius' ) )  THEN
              CALL get_variable( pids_id, 'rotor_radius', rotor_radius(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'nacelle_cd' ) )  THEN
              CALL get_variable( pids_id, 'nacelle_cd', nacelle_cd(1:n_turbines) )
           ENDIF
    
           IF ( check_existence( vars_pids, 'tower_cd' ) )  THEN
              CALL get_variable( pids_id, 'tower_cd', tower_cd(1:n_turbines) )
           ENDIF
    !
    !--    Close wtm input file:
           CALL close_input_file( pids_id )
    
        ENDIF
    #endif
    
    !
    !-- To be able to allocate arrays with dimension of rotor rings and segments,
    !-- the maximum possible numbers of rings and segments have to be calculated:
        ALLOCATE( nrings(1:n_turbines) )
        ALLOCATE( delta_r(1:n_turbines) )
    
        nrings(:)  = 0
        delta_r(:) = 0.0_wp
    
    !
    !-- Thickness (radial) of each ring and length (tangential) of each segment:
        delta_r_factor = segment_width_radial
        delta_t_factor = segment_length_tangential
        delta_r_init   = delta_r_factor * MIN( dx, dy, dz(1) )
    
        DO  inot = 1, n_turbines
    !
    !--    Determine number of rings:
           nrings(inot) = NINT( rotor_radius(inot) / delta_r_init )
    
           delta_r(inot) = rotor_radius(inot) / nrings(inot)
    
        ENDDO
    
        nrings_max = MAXVAL( nrings )
    
        ALLOCATE( nsegs(1:nrings_max,1:n_turbines) )
        ALLOCATE( nsegs_total(1:n_turbines) )
    
        nsegs(:,:)     = 0
        nsegs_total(:) = 0
    
    
        DO  inot = 1, n_turbines
           DO  ring = 1, nrings(inot)
    !
    !--       Determine number of segments for each ring:
              nsegs(ring,inot) = MAX( 8, CEILING( delta_r_factor * pi * ( 2.0_wp * ring - 1.0_wp ) /   &
                                                  delta_t_factor ) )
           ENDDO
    !
    !--    Total sum of all rotor segments:
           nsegs_total(inot) = SUM( nsegs(:,inot) )
        ENDDO
    
    !
    !-- Maximum number of segments per ring:
        nsegs_max = MAXVAL( nsegs )    
    
    !
    !-- Allocate 1D arrays (dimension = number of turbines):
        ALLOCATE( i_hub(1:n_turbines) )
        ALLOCATE( i_smear(1:n_turbines) )
        ALLOCATE( j_hub(1:n_turbines) )
        ALLOCATE( j_smear(1:n_turbines) )
        ALLOCATE( k_hub(1:n_turbines) )
        ALLOCATE( k_smear(1:n_turbines) )
        ALLOCATE( torque_total(1:n_turbines) )
        ALLOCATE( thrust_rotor(1:n_turbines) )
    
    !
    !-- Allocation of the 1D arrays for yaw control:
        ALLOCATE( yawdir(1:n_turbines) )
        ALLOCATE( u_inflow(1:n_turbines) )
        ALLOCATE( wdir(1:n_turbines) )
        ALLOCATE( u_inflow_l(1:n_turbines) )
        ALLOCATE( wdir_l(1:n_turbines) )
        ALLOCATE( yaw_angle_l(1:n_turbines) )
    
    !
    !-- Allocate 1D arrays (dimension = number of rotor segments):
        ALLOCATE( alpha_attack(1:nsegs_max) )
        ALLOCATE( chord(1:nsegs_max) )
        ALLOCATE( phi_rel(1:nsegs_max) )
        ALLOCATE( thrust_seg(1:nsegs_max) )
        ALLOCATE( torque_seg_y(1:nsegs_max) )
        ALLOCATE( torque_seg_z(1:nsegs_max) )
        ALLOCATE( turb_cd(1:nsegs_max) )
        ALLOCATE( turb_cl(1:nsegs_max) )
        ALLOCATE( vrel(1:nsegs_max) )
        ALLOCATE( vtheta(1:nsegs_max) )
    
    !
    !-- Allocate 2D arrays (dimension = number of rotor rings and segments):
        ALLOCATE( rbx(1:nrings_max,1:nsegs_max) )
        ALLOCATE( rby(1:nrings_max,1:nsegs_max) )
        ALLOCATE( rbz(1:nrings_max,1:nsegs_max) )
        ALLOCATE( thrust_ring(1:nrings_max,1:nsegs_max) )
        ALLOCATE( torque_ring_y(1:nrings_max,1:nsegs_max) )
        ALLOCATE( torque_ring_z(1:nrings_max,1:nsegs_max) )
    
    !
    !-- Allocate additional 2D arrays:
        ALLOCATE( rotx(1:n_turbines,1:3) )
        ALLOCATE( roty(1:n_turbines,1:3) )
        ALLOCATE( rotz(1:n_turbines,1:3) )
    
    !
    !-- Allocate 3D arrays (dimension = number of grid points):
        ALLOCATE( nac_cd_surf(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( rot_tend_x(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( rot_tend_y(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( rot_tend_z(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( thrust(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( torque_y(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( torque_z(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
        ALLOCATE( tow_cd_surf(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    
    !
    !-- Allocate additional 3D arrays:
        ALLOCATE( u_int(1:n_turbines,1:nrings_max,1:nsegs_max) )
        ALLOCATE( u_int_1_l(1:n_turbines,1:nrings_max,1:nsegs_max) )
        ALLOCATE( v_int(1:n_turbines,1:nrings_max,1:nsegs_max) )
        ALLOCATE( v_int_1_l(1:n_turbines,1:nrings_max,1:nsegs_max) )
        ALLOCATE( w_int(1:n_turbines,1:nrings_max,1:nsegs_max) )
        ALLOCATE( w_int_1_l(1:n_turbines,1:nrings_max,1:nsegs_max) )
    
    !
    !-- All of the arrays are initialized with a value of zero:
        i_hub(:)                 = 0
        i_smear(:)               = 0
        j_hub(:)                 = 0
        j_smear(:)               = 0
        k_hub(:)                 = 0
        k_smear(:)               = 0
    
        torque_total(:)          = 0.0_wp
        thrust_rotor(:)          = 0.0_wp
    
        IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
           generator_speed(:)             = 0.0_wp
           generator_speed_old(:)         = 0.0_wp
           generator_speed_f(:)           = 0.0_wp
           generator_speed_f_old(:)       = 0.0_wp
           pitch_angle_old(:)             = 0.0_wp
           torque_gen(:)                  = 0.0_wp
           torque_gen_old(:)              = 0.0_wp
        ENDIF
    
        yawdir(:)                = 0.0_wp
        wdir_l(:)                = 0.0_wp
        wdir(:)                  = 0.0_wp
        u_inflow(:)              = 0.0_wp
        u_inflow_l(:)            = 0.0_wp
        yaw_angle_l(:)           = 0.0_wp
    
    !
    !-- Allocate 1D arrays (dimension = number of rotor segments):
        alpha_attack(:)          = 0.0_wp
        chord(:)                 = 0.0_wp
        phi_rel(:)               = 0.0_wp
        thrust_seg(:)            = 0.0_wp
        torque_seg_y(:)          = 0.0_wp
        torque_seg_z(:)          = 0.0_wp
        turb_cd(:)               = 0.0_wp
        turb_cl(:)               = 0.0_wp
        vrel(:)                  = 0.0_wp
        vtheta(:)                = 0.0_wp
    
        rbx(:,:)                 = 0.0_wp
        rby(:,:)                 = 0.0_wp
        rbz(:,:)                 = 0.0_wp
        thrust_ring(:,:)         = 0.0_wp
        torque_ring_y(:,:)       = 0.0_wp
        torque_ring_z(:,:)       = 0.0_wp
    
        rotx(:,:)                = 0.0_wp
        roty(:,:)                = 0.0_wp
        rotz(:,:)                = 0.0_wp
    
        nac_cd_surf(:,:,:)       = 0.0_wp
        rot_tend_x(:,:,:)        = 0.0_wp
        rot_tend_y(:,:,:)        = 0.0_wp
        rot_tend_z(:,:,:)        = 0.0_wp
        thrust(:,:,:)            = 0.0_wp
        torque_y(:,:,:)          = 0.0_wp
        torque_z(:,:,:)          = 0.0_wp
        tow_cd_surf(:,:,:)       = 0.0_wp
    
        u_int(:,:,:)             = 0.0_wp
        u_int_1_l(:,:,:)         = 0.0_wp
        v_int(:,:,:)             = 0.0_wp
        v_int_1_l(:,:,:)         = 0.0_wp
        w_int(:,:,:)             = 0.0_wp
        w_int_1_l(:,:,:)         = 0.0_wp
    
    
     END SUBROUTINE wtm_init_arrays
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Initialization of the wind turbine model
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_init
    
        USE control_parameters,                                                                        &
            ONLY:  dz_stretch_level_start
    
        USE exchange_horiz_mod,                                                                        &
            ONLY:  exchange_horiz
    
        IMPLICIT NONE
    
        INTEGER(iwp) ::  i  !< running index
        INTEGER(iwp) ::  j  !< running index
        INTEGER(iwp) ::  k  !< running index
    
    !
    !-- Help variables for the smearing function:
        REAL(wp) ::  eps_kernel  !<
    !
    !-- Help variables for calculation of the tower drag:
        INTEGER(iwp) ::  tower_n  !<
        INTEGER(iwp) ::  tower_s  !<
    !
    !-- Help variables for the calculation of the nacelle drag:
        INTEGER(iwp) ::  i_ip         !<
        INTEGER(iwp) ::  i_ipg        !<
        
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_nacb  !<
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_nacl  !<
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_nacr  !<
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_nact  !<
    
        REAL(wp) ::  dy_int    !<
        REAL(wp) ::  dz_int    !<
        REAL(wp) ::  sqrt_arg  !<
        REAL(wp) ::  yvalue    !<
    
        REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  circle_points  !<
    
    
        IF ( debug_output )  CALL debug_message( 'wtm_init', 'start' )
    
        ALLOCATE( index_nacb(1:n_turbines) )
        ALLOCATE( index_nacl(1:n_turbines) )
        ALLOCATE( index_nacr(1:n_turbines) )
        ALLOCATE( index_nact(1:n_turbines) )
    
    !
    !-- Calculation of parameters for the regularization kernel (smearing of the forces)
    !
    !-- In the following, some of the required parameters for the smearing will be calculated:
    !-- The kernel is set equal to twice the grid spacing which has turned out to be a reasonable
    !-- value (see e.g. Troldborg et al. (2013), Wind Energy, DOI: 10.1002/we.1608).
        eps_kernel = smearing_kernel_size * dx
    !
    !-- The zero point (eps_min) of the polynomial function must be the following if the integral of
    !-- the polynomial function (for values < eps_min) shall be equal to the integral of the Gaussian
    !-- function used before.
        eps_min = ( 105.0_wp / 32.0_wp )**( 1.0_wp / 3.0_wp ) * pi**( 1.0_wp / 6.0_wp ) * eps_kernel
    !
    !-- Stretching (non-uniform grid spacing) is not considered in the wind turbine model.
    !-- Therefore, vertical stretching has to be applied above the area where the wtm is active.
    !-- ABS (...) is required because the default value of dz_stretch_level_start is -9999999.9_wp
    !-- (negative).
        IF ( ABS( dz_stretch_level_start(1) ) <=                                                       &
             MAXVAL( hub_z(1:n_turbines) ) + MAXVAL( rotor_radius(1:n_turbines) ) + eps_min )          &
        THEN
           WRITE( message_string, * ) 'vertical grid stretching only allowed above ',                  &
                                      MAXVAL( hub_z(1:n_turbines) ) +                                  &
                                      MAXVAL( rotor_radius(1:n_turbines) ) + eps_min, ' m'
           CALL message( 'wtm_init', 'WTM0004', 1, 2, 0, 6, 0 )
        ENDIF
    
        eps_min2 = eps_min**2
    !
    !-- Parameters in the polynomial function.
        pol_a = 1.0_wp / eps_min**4
        pol_b = 2.0_wp / eps_min**2
    !
    !-- Normalization factor which is the inverse of the integral of the smearing function.
        eps_factor = 105.0_wp / ( 32.0_wp * pi * eps_min**3 )
    
    !-- Change tilt angle to rad.
        tilt_angle = tilt_angle * pi / 180.0_wp
    
    !
    !-- Change yaw angle to rad.
        IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
           yaw_angle(:) = yaw_angle(:) * pi / 180.0_wp
        ENDIF
    
    
        DO  inot = 1, n_turbines
    !
    !--    Rotate the rotor coordinates in case yaw and tilt are defined.
           CALL wtm_rotate_rotor( inot )
    
    !
    !--    Determine the indices of the hub height.
           i_hub(inot) = INT(   hub_x(inot)                    / dx )
           j_hub(inot) = INT( ( hub_y(inot) + 0.5_wp * dy )    / dy )
           k_hub(inot) = INT( ( hub_z(inot) + 0.5_wp * dz(1) ) / dz(1) )
    
    !
    !--    Determining the area to which the smearing of the forces is applied.
    !--    As smearing now is effectively applied only for distances smaller than eps_min, the
    !--    smearing area can be further limited and regarded as a function of eps_min.
           i_smear(inot) = CEILING( ( rotor_radius(inot) + eps_min ) / dx )
           j_smear(inot) = CEILING( ( rotor_radius(inot) + eps_min ) / dy )
           k_smear(inot) = CEILING( ( rotor_radius(inot) + eps_min ) / dz(1) )
    
        ENDDO
    
    !
    !-- Call the wtm_init_speed_control subroutine and calculate the local rotor_speed for the
    !-- respective processor.
        IF ( speed_control)  THEN
    
           CALL wtm_init_speed_control
    
           IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
    
              DO  inot = 1, n_turbines
    
                 IF ( nxl > i_hub(inot) )  THEN
                    torque_gen(inot)        = 0.0_wp
                    generator_speed_f(inot) = 0.0_wp
                    rotor_speed_l(inot)     = 0.0_wp
                 ENDIF
    
                 IF ( nxr < i_hub(inot) )  THEN
                    torque_gen(inot)        = 0.0_wp
                    generator_speed_f(inot) = 0.0_wp
                    rotor_speed_l(inot)     = 0.0_wp
                 ENDIF
    
                 IF ( nys > j_hub(inot) )  THEN
                    torque_gen(inot)        = 0.0_wp
                    generator_speed_f(inot) = 0.0_wp
                    rotor_speed_l(inot)     = 0.0_wp
                 ENDIF
    
                 IF ( nyn < j_hub(inot) )  THEN
                    torque_gen(inot)        = 0.0_wp
                    generator_speed_f(inot) = 0.0_wp
                    rotor_speed_l(inot)     = 0.0_wp
                 ENDIF
    
                 IF ( ( nxl <= i_hub(inot) ) .AND. ( nxr >= i_hub(inot) ) )  THEN
                    IF ( ( nys <= j_hub(inot) ) .AND. ( nyn >= j_hub(inot) ) )  THEN
    
                       rotor_speed_l(inot) = generator_speed(inot) / gear_ratio
    
                    ENDIF
                 ENDIF
    
              END DO
    
           ENDIF
    
        ENDIF
    
    !
    !-- Determine the area within each grid cell that overlaps with the area of the nacelle and the
    !-- tower (needed for calculation of the forces).
    !-- Note: so far this is only a 2D version, in that the mean flow is perpendicular to the rotor
    !--       area.
    !
    !-- Allocation of the array containing information on the intersection points between rotor disk
    !-- and the numerical grid.
        upper_end = ( ny + 1 ) * 10000
        ALLOCATE( circle_points(1:2,1:upper_end) )
        circle_points(:,:) = 0.0_wp
    
        DO  inot = 1, n_turbines
    !
    !--    Determine the grid index (u-grid) that corresponds to the location of the rotor center
    !--    (reduces the amount of calculations in the case that the mean flow is perpendicular to the
    !--    rotor area).
           i = i_hub(inot)
    !
    !--    Determine the left and the right edge of the nacelle (corresponding grid point indices).
           index_nacl(inot) = INT( ( hub_y(inot) - nacelle_radius(inot) + 0.5_wp * dy ) / dy )
           index_nacr(inot) = INT( ( hub_y(inot) + nacelle_radius(inot) + 0.5_wp * dy ) / dy )
    !
    !--    Determine the bottom and the top edge of the nacelle (corresponding grid point indices).
    !--    The grid point index has to be increased by 1, as the first level for the u-component
    !--    (index 0) is situated below the surface. All points between z=0 and z=dz/s would already
    !--    be contained in grid box 1.
           index_nacb(inot) = INT( ( hub_z(inot) - nacelle_radius(inot) ) / dz(1) ) + 1
           index_nact(inot) = INT( ( hub_z(inot) + nacelle_radius(inot) ) / dz(1) ) + 1
    !
    !--    Determine the indices of the grid boxes containing the left and the right boundaries of
    !--    the tower.
           tower_n = ( hub_y(inot) + 0.5_wp * tower_diameter(inot) - 0.5_wp * dy ) / dy
           tower_s = ( hub_y(inot) - 0.5_wp * tower_diameter(inot) - 0.5_wp * dy ) / dy
    !
    !--    Determine the fraction of the grid box area overlapping with the tower area and multiply
    !--    it with the drag of the tower.
           IF ( ( nxlg <= i )  .AND.  ( nxrg >= i ) )  THEN
    
              DO  j = nys, nyn
    !
    !--          Loop from south to north boundary of tower.
                 IF ( ( j >= tower_s )  .AND.  ( j <= tower_n ) )  THEN
    
                    DO  k = nzb, nzt
    
                       IF ( k == k_hub(inot) )  THEN
                          IF ( tower_n - tower_s >= 1 )  THEN
    !
    !--                      Leftmost and rightmost grid box.
                             IF ( j == tower_s )  THEN
                                tow_cd_surf(k,j,i) =                                                   &
                                          ( hub_z(inot) - ( k_hub(inot) * dz(1) - 0.5_wp * dz(1) ) ) * &
                                          ( ( tower_s + 1.0_wp + 0.5_wp ) * dy    -                    &
                                            ( hub_y(inot) - 0.5_wp * tower_diameter(inot) )            &
                                          ) * tower_cd(inot)
    
                             ELSEIF ( j == tower_n )  THEN
                                tow_cd_surf(k,j,i) =                                                   &
                                          ( hub_z(inot) - ( k_hub(inot) * dz(1) - 0.5_wp * dz(1) ) ) * &
                                          ( ( hub_y(inot) + 0.5_wp * tower_diameter(inot) ) -          &
                                            ( tower_n + 0.5_wp ) * dy                                  &
                                          ) * tower_cd(inot)
    !
    !--                      Grid boxes inbetween (where tow_cd_surf = grid box area).
                             ELSE
                                tow_cd_surf(k,j,i) = ( hub_z(inot) -                                   &
                                                       ( k_hub(inot) * dz(1) - 0.5_wp * dz(1) )        &
                                                     ) * dy * tower_cd(inot)
                             ENDIF
    !
    !--                   Tower lies completely within one grid box.
                          ELSE
                             tow_cd_surf(k,j,i) = ( hub_z(inot) -                                      &
                                                    ( k_hub(inot) * dz(1) - 0.5_wp * dz(1) )           &
                                                  ) * tower_diameter(inot) * tower_cd(inot)
                          ENDIF
    !
    !--                In case that k is smaller than k_hub the following actions are carried out.
                       ELSEIF ( k < k_hub(inot) )  THEN
    
                          IF ( ( tower_n - tower_s ) >= 1 )  THEN
    !
    !--                      Leftmost and rightmost grid box.
                             IF ( j == tower_s )  THEN
                                tow_cd_surf(k,j,i) = dz(1) * ( ( tower_s + 1 + 0.5_wp ) * dy -         &
                                                       ( hub_y(inot) - 0.5_wp * tower_diameter(inot) ) &
                                                             ) * tower_cd(inot)
                             ELSEIF ( j == tower_n )  THEN
                                tow_cd_surf(k,j,i) = dz(1) * (                                         &
                                                       ( hub_y(inot) + 0.5_wp * tower_diameter(inot) ) &
                                                               - ( tower_n + 0.5_wp ) * dy             &
                                                             ) * tower_cd(inot)
    !
    !--                      Grid boxes inbetween (where tow_cd_surf = grid box area).
                             ELSE
                                tow_cd_surf(k,j,i) = dz(1) * dy * tower_cd(inot)
                             ENDIF
    !
    !--                   Tower lies completely within one grid box.
                          ELSE
    
                             tow_cd_surf(k,j,i) = dz(1) * tower_diameter(inot) * tower_cd(inot)
    
                          ENDIF ! larger than grid box
    
                       ENDIF ! k == k_hub
    
                    ENDDO ! over k
    
                 ENDIF ! inside north and south boundary of tower
    
              ENDDO ! over j
    
           ENDIF ! hub inside domain + ghostpoints
    
    
           CALL exchange_horiz( tow_cd_surf, nbgp )
    
    !
    !--    Calculation of the nacelle area.
    !--    Tabulate the points on the circle that are required in the following for the calculation
    !--    of the Riemann integral (node points; they are called circle_points in the following).
           dy_int = dy / 10000.0_wp
    
           IF ( ( nacelle_cd(inot) /= 0.0_wp ) .AND. ( nacelle_radius(inot) /= 0.0_wp ) ) THEN
    
              DO  i_ip = 1, upper_end
    
                 yvalue   = dy_int * ( i_ip - 0.5_wp ) + 0.5_wp * dy
                 sqrt_arg = nacelle_radius(inot)**2 - ( yvalue - hub_y(inot) )**2
    
                 IF ( sqrt_arg >= 0.0_wp )  THEN
    !
    !--             Bottom intersection point.
                    circle_points(1,i_ip) = hub_z(inot) - SQRT( sqrt_arg )
    !
    !--             Top intersection point.
                    circle_points(2,i_ip) = hub_z(inot) + SQRT( sqrt_arg )
                 ELSE
                    circle_points(:,i_ip) = -111111
                 ENDIF
    
             ENDDO
    
             DO  j = nys, nyn
    !        
    !--         In case that the grid box is located completely outside the nacelle (y) it can
    !--         automatically be stated that there is no overlap between the grid box and the nacelle
    !--         and consequently we can set nac_cd_surf(:,j,i) = 0.0.
                IF ( ( j >= index_nacl(inot) )  .AND.  ( j <= index_nacr(inot) ) )  THEN
    
                   DO  k = nzb+1, nzt
    !        
    !--               In case that the grid box is located completely outside the nacelle (z) it can
    !--               automatically be stated that there is no overlap between the grid box and the
    !--               nacelle and consequently we can set nac_cd_surf(k,j,i) = 0.0.
                      IF ( ( k >= index_nacb(inot) )  .OR.  ( k <= index_nact(inot) ) )  THEN
    !        
    !--                  For all other cases Riemann integrals are calculated. Here, the points on the
    !--                  circle that have been determined above are used in order to calculate the
    !--                  overlap between the gridbox and the nacelle area (area approached by 10000
    !--                  rectangulars dz_int * dy_int).
                         DO  i_ipg = 1, 10000
    
                            dz_int = dz(k)
                            i_ip = j * 10000 + i_ipg
    !        
    !--                     Determine the vertical extension dz_int of the circle within the current
    !--                     grid box.
                            IF ( ( circle_points(2,i_ip) < zw(k) )  .AND.                              &
                                 ( circle_points(2,i_ip) >= zw(k-1) ) )                                &
                            THEN
                               dz_int = dz_int - ( zw(k) - circle_points(2,i_ip) )
                            ENDIF
    
                            IF ( ( circle_points(1,i_ip) <= zw(k) )  .AND.                             &
                                 ( circle_points(1,i_ip) > zw(k-1) ) )                                 &
                            THEN
                               dz_int = dz_int - ( circle_points(1,i_ip) - zw(k-1) )
                            ENDIF
    
                            IF ( zw(k-1) > circle_points(2,i_ip) )  dz_int = 0.0_wp
                            IF ( zw(k)   < circle_points(1,i_ip) )  dz_int = 0.0_wp
    
                            IF ( ( nxlg <= i )  .AND.  ( nxrg >= i ) )  THEN
                               nac_cd_surf(k,j,i) = nac_cd_surf(k,j,i) +                               &
                                                    dy_int * dz_int * nacelle_cd(inot)
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
    
             ENDDO
    
             CALL exchange_horiz( nac_cd_surf, nbgp )
    
           ENDIF  
    
        ENDDO ! over turbines
    !
    !-- Normalize tower and nacelle drag.
        tow_cd_surf = tow_cd_surf / ( dx * dy * dz(1) )
        nac_cd_surf = nac_cd_surf / ( dx * dy * dz(1) )
    
        CALL wtm_read_blade_tables
    
        IF ( debug_output )  CALL debug_message( 'wtm_init', 'end' )
    
     END SUBROUTINE wtm_init
    
    
     SUBROUTINE wtm_init_output
    
    
    !    INTEGER(iwp) ::  ntimesteps              !< number of timesteps defined in NetCDF output file
    !    INTEGER(iwp) ::  ntimesteps_max = 80000  !< number of maximum timesteps defined in NetCDF output file
        INTEGER(iwp) ::  return_value            !< returned status value of called function
    
        INTEGER(iwp) ::  n  !< running index
    
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ndim !< dummy to write dimension
    
    
    !
    !-- Create NetCDF output file:
    #if defined( __netcdf4 )
        nc_filename = 'DATA_1D_TS_WTM_NETCDF'
        return_value = dom_def_file( nc_filename, 'netcdf4-serial' )
    #else
        message_string = 'Wind turbine model output requires netCDF version 4. &' //                   &
                         'No output file will be created.'
        CALL message( 'wtm_init_output', 'WTM0005', 0, 1, 0, 6, 0 )
    #endif
        IF ( myid == 0 )  THEN
    !
    !--    Define dimensions in output file:
           ALLOCATE( ndim(1:n_turbines) )
           DO  n = 1, n_turbines
              ndim(n) = n
           ENDDO
           return_value = dom_def_dim( nc_filename,                                                    &
                                       dimension_name = 'turbine',                                     &
                                       output_type = 'int32',                                          &
                                       bounds = (/1_iwp, n_turbines/),                                 &
                                       values_int32 = ndim )
           DEALLOCATE( ndim )
    
           return_value = dom_def_dim( nc_filename,                                                    &
                                       dimension_name = 'time',                                        &
                                       output_type = 'real32',                                         &
                                       bounds = (/1_iwp/),                                             &
                                       values_realwp = (/0.0_wp/) )
    
           variable_name = 'x'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/'turbine'/),                                &
                                       output_type = 'real32' )
    
           variable_name = 'y'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/'turbine'/),                                &
                                       output_type = 'real32' )
    
           variable_name = 'z'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/'turbine'/),                                &
                                       output_type = 'real32' )
    
    
           variable_name = 'turbine_name'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/'turbine'/),                                &
                                       output_type = 'real32' )
    
    
           variable_name = 'rotor_diameter'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/'turbine'/),                                &
                                       output_type = 'real32' )
    
           variable_name = 'tower_diameter'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/'turbine'/),                                &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'time',                                         &
                                       attribute_name = 'units',                                       &
                                       value = 'seconds since ' // origin_date_time )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'x',                                            &
                                       attribute_name = 'units',                                       &
                                       value = 'm' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'y',                                            &
                                       attribute_name = 'units',                                       &
                                       value = 'm' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'z',                                            &
                                       attribute_name = 'units',                                       &
                                       value = 'm' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'rotor_diameter',                               &
                                       attribute_name = 'units',                                       &
                                       value = 'm' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'tower_diameter',                               &
                                       attribute_name = 'units',                                       &
                                       value = 'm' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'x',                                            &
                                       attribute_name = 'long_name',                                   &
                                       value = 'x location of rotor center' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'y',                                            &
                                       attribute_name = 'long_name',                                   &
                                       value = 'y location of rotor center' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'z',                                            &
                                       attribute_name = 'long_name',                                   &
                                       value = 'z location of rotor center' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'turbine_name',                                 &
                                       attribute_name = 'long_name',                                   &
                                       value = 'turbine name')
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'rotor_diameter',                               &
                                       attribute_name = 'long_name',                                   &
                                       value = 'rotor diameter')
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'tower_diameter',                               &
                                       attribute_name = 'long_name',                                   &
                                       value = 'tower diameter')
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'time',                                         &
                                       attribute_name = 'standard_name',                               &
                                       value = 'time')
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'time',                                         &
                                       attribute_name = 'axis',                                        &
                                       value = 'T')
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'x',                                            &
                                       attribute_name = 'axis',                                        &
                                       value = 'X' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = 'y',                                            &
                                       attribute_name = 'axis',                                        &
                                       value = 'Y' )
    
    
           variable_name = 'generator_power'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'W' )
    
           variable_name = 'generator_speed'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'rad/s' )
    
           variable_name = 'generator_torque'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'Nm' )
    
           variable_name = 'pitch_angle'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'degrees' )
    
           variable_name = 'rotor_power'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'W' )
    
           variable_name = 'rotor_speed'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'rad/s' )
    
           variable_name = 'rotor_thrust'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'N' )
    
           variable_name = 'rotor_torque'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'Nm' )
    
           variable_name = 'wind_direction'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'degrees' )
    
           variable_name = 'yaw_angle'
           return_value = dom_def_var( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       dimension_names = (/ 'turbine', 'time   ' /),                   &
                                       output_type = 'real32' )
    
           return_value = dom_def_att( nc_filename,                                                    &
                                       variable_name = variable_name,                                  &
                                       attribute_name = 'units',                                       &
                                       value = 'degrees' )
    
    !
    !--    Check if DOM reported any error
           dom_error_message = dom_get_error_message()
           IF ( TRIM( dom_error_message ) /= '' )  THEN
              message_string = 'error while defining output: "' // TRIM( dom_error_message ) // '"'
              CALL message( 'wtm_init_output', 'WTM0006', 0, 1, 0, 6, 0 )
           ENDIF
    
        ENDIF
    
     END SUBROUTINE wtm_init_output
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Read in layout of the rotor blade , the lift and drag tables and the distribution of lift and
    !> drag tables along the blade
    !--------------------------------------------------------------------------------------------------!
    !
     SUBROUTINE wtm_read_blade_tables
    
    
        IMPLICIT NONE
    
        CHARACTER(200) :: chmess  !< Read in string
    
        INTEGER(iwp) ::  ii  !< running index
        INTEGER(iwp) ::  jj  !< running index
    
        INTEGER(iwp) ::  ierrn  !<
    
        INTEGER(iwp) ::  dlen        !< no. rows of local table
        INTEGER(iwp) ::  dlenbl      !< no. rows of cd, cl table
        INTEGER(iwp) ::  dlenbl_int  !< no. rows of interpolated cd, cl tables
        INTEGER(iwp) ::  ialpha      !< table position of current alpha value
        INTEGER(iwp) ::  iialpha     !<
        INTEGER(iwp) ::  iir         !<
        INTEGER(iwp) ::  radres      !< radial resolution
        INTEGER(iwp) ::  t1          !< no. of airfoil
        INTEGER(iwp) ::  t2          !< no. of airfoil
        INTEGER(iwp) ::  trow        !<
    
    
        REAL(wp) :: alpha_attack_i  !<
        REAL(wp) :: weight_a        !<
        REAL(wp) :: weight_b        !<
    
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: ttoint1  !<
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: ttoint2  !<
    
        REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cd_sel1  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cd_sel2  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cl_sel1  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cl_sel2  !<
        REAL(wp), DIMENSION(:), ALLOCATABLE :: read_cl_cd    !< read in var array
    
        REAL(wp), DIMENSION(:),   ALLOCATABLE :: alpha_attack_tab  !<
        REAL(wp), DIMENSION(:),   ALLOCATABLE :: trad1             !<
        REAL(wp), DIMENSION(:),   ALLOCATABLE :: trad2             !<
        REAL(wp), DIMENSION(:,:), ALLOCATABLE :: turb_cd_table     !<
        REAL(wp), DIMENSION(:,:), ALLOCATABLE :: turb_cl_table     !<
    
        ALLOCATE ( read_cl_cd(1:2 * n_airfoils + 1) )
    
    !
    !-- Read in the distribution of lift and drag tables along the blade, the layout of the rotor
    !-- blade and the lift and drag tables:
        OPEN ( 90, FILE='WTM_DATA', STATUS='OLD', FORM='FORMATTED', IOSTAT=ierrn )
    
        IF ( ierrn /= 0 )  THEN
           message_string = 'file WTM_DATA does not exist'
           CALL message( 'wtm_init', 'WTM0007', 1, 2, 0, 6, 0 )
        ENDIF
    !
    !-- Read distribution table:
        dlen = 0
    
        READ ( 90, '(3/)' )
    
        rloop3: DO
           READ ( 90, *, IOSTAT=ierrn ) chmess
           IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop3
           dlen = dlen + 1
        ENDDO rloop3
    
        ALLOCATE( trad1(1:dlen), trad2(1:dlen), ttoint1(1:dlen), ttoint2(1:dlen) )
    
        DO  jj = 1, dlen+1
           BACKSPACE ( 90, IOSTAT=ierrn )
        ENDDO
    
        DO  jj = 1, dlen
           READ ( 90, * ) trad1(jj), trad2(jj), ttoint1(jj), ttoint2(jj)
        ENDDO
    
    !
    !-- Read layout table:
        dlen = 0
    
        READ ( 90, '(3/)')
    
        rloop1: DO
           READ ( 90, *, IOSTAT=ierrn ) chmess
           IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop1
           dlen = dlen + 1
        ENDDO rloop1
    
        ALLOCATE( lrd(1:dlen), ard(1:dlen), crd(1:dlen) )
        DO  jj = 1, dlen + 1
           BACKSPACE ( 90, IOSTAT=ierrn )
        ENDDO
        DO  jj = 1, dlen
           READ ( 90, * ) lrd(jj), ard(jj), crd(jj)
        ENDDO
    
    !
    !-- Read tables (turb_cl(alpha),turb_cd(alpha) for the different profiles:
        dlen = 0
    
        READ ( 90, '(3/)' )
    
        rloop2: DO
           READ ( 90, *, IOSTAT=ierrn ) chmess
           IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop2
           dlen = dlen + 1
        ENDDO rloop2
    
        ALLOCATE( alpha_attack_tab(1:dlen), turb_cl_table(1:dlen,1:n_airfoils), &
                  turb_cd_table(1:dlen,1:n_airfoils) )
    
        DO  jj = 1, dlen + 1
           BACKSPACE ( 90, IOSTAT=ierrn )
        ENDDO
    
        DO  jj = 1, dlen
           READ ( 90, * ) read_cl_cd
           alpha_attack_tab(jj) = read_cl_cd(1)
           DO  ii = 1, n_airfoils
              turb_cl_table(jj,ii) = read_cl_cd(ii * 2)
              turb_cd_table(jj,ii) = read_cl_cd(ii * 2 + 1)
           ENDDO
    
        ENDDO
    
        dlenbl = dlen
    
        CLOSE ( 90 )
    
    !
    !-- For each possible radial position (resolution: 0.1 m --> 631 values if rotor_radius(1)=63m)
    !-- and each possible angle of attack (resolution: 0.1 degrees --> 3601 values!) determine the
    !-- lift and drag coefficient by interpolating between the tabulated values of each table
    !-- (interpolate to current angle of attack) and between the tables (interpolate to current
    !-- radial position):
        ALLOCATE( turb_cl_sel1(1:dlenbl) )
        ALLOCATE( turb_cl_sel2(1:dlenbl) )
        ALLOCATE( turb_cd_sel1(1:dlenbl) )
        ALLOCATE( turb_cd_sel2(1:dlenbl) )
    
        radres     = INT( rotor_radius(1) * 10.0_wp ) + 1_iwp
        dlenbl_int = INT( 360.0_wp / accu_cl_cd_tab ) + 1_iwp
    
        ALLOCATE( turb_cl_tab(1:dlenbl_int,1:radres) )
        ALLOCATE( turb_cd_tab(1:dlenbl_int,1:radres) )
    
        DO  iir = 1, radres ! loop over radius
    
           cur_r = ( iir - 1_iwp ) * 0.1_wp
    !
    !--    Find position in table 1:
           lct = MINLOC( ABS( trad1 - cur_r ) )
    !          lct(1) = lct(1)
    
           IF ( ( trad1(lct(1)) - cur_r ) > 0.0 )  THEN
              lct(1) = lct(1) - 1
           ENDIF
    
           trow = lct(1)
    !
    !--    Calculate weights for radius interpolation:
           weight_a = ( trad2(trow) - cur_r ) / ( trad2(trow) - trad1(trow) )
           weight_b = ( cur_r - trad1(trow) ) / ( trad2(trow) - trad1(trow) )
           t1 = ttoint1(trow)
           t2 = ttoint2(trow)
    
           IF ( t1 == t2 )  THEN  ! if both are the same, the weights are NaN
              weight_a = 0.5_wp   ! then do interpolate in between same twice
              weight_b = 0.5_wp   ! using 0.5 as weight
           ENDIF
    
           IF ( t1 == 0 .AND. t2 == 0 )  THEN
              turb_cd_sel1 = 0.0_wp
              turb_cd_sel2 = 0.0_wp
              turb_cl_sel1 = 0.0_wp
              turb_cl_sel2 = 0.0_wp
    
              turb_cd_tab(1,iir) = 0.0_wp  ! For -180 degrees (iialpha=1) the values
              turb_cl_tab(1,iir) = 0.0_wp  ! for each radius has to be set
                                           ! explicitly
           ELSE
              turb_cd_sel1 = turb_cd_table(:,t1)
              turb_cd_sel2 = turb_cd_table(:,t2)
              turb_cl_sel1 = turb_cl_table(:,t1)
              turb_cl_sel2 = turb_cl_table(:,t2)
    !
    !--       For -180 degrees (iialpha=1) the values for each radius has to be set explicitly:
              turb_cd_tab(1,iir) = ( weight_a * turb_cd_table(1,t1) + weight_b * turb_cd_table(1,t2) )
              turb_cl_tab(1,iir) = ( weight_a * turb_cl_table(1,t1) + weight_b * turb_cl_table(1,t2) )
           ENDIF
    
           DO  iialpha = 2, dlenbl_int  ! loop over angles
    
              alpha_attack_i = -180.0_wp + REAL( iialpha-1 ) * accu_cl_cd_tab
              ialpha = 1
    
              DO WHILE  ( ( alpha_attack_i > alpha_attack_tab(ialpha) ) .AND. ( ialpha < dlen ) )
                 ialpha = ialpha + 1
              ENDDO
    
    !
    !--       Interpolation of lift and drag coefficiencts on fine grid of radius segments and angles
    !--       of attack:
              turb_cl_tab(iialpha,iir) = ( alpha_attack_tab(ialpha) - alpha_attack_i ) /               &
                                         ( alpha_attack_tab(ialpha) - alpha_attack_tab(ialpha-1) ) *   &
                                         ( weight_a * turb_cl_sel1(ialpha-1) +                         &
                                           weight_b * turb_cl_sel2(ialpha-1) ) +                       &
                                         ( alpha_attack_i - alpha_attack_tab(ialpha-1) ) /             &
                                         ( alpha_attack_tab(ialpha) - alpha_attack_tab(ialpha-1) ) *   &
                                         ( weight_a * turb_cl_sel1(ialpha) +                           &
                                           weight_b * turb_cl_sel2(ialpha) )
              turb_cd_tab(iialpha,iir) = ( alpha_attack_tab(ialpha) - alpha_attack_i ) /               &
                                         ( alpha_attack_tab(ialpha) - alpha_attack_tab(ialpha-1) ) *   &
                                         ( weight_a * turb_cd_sel1(ialpha-1) +                         &
                                           weight_b * turb_cd_sel2(ialpha-1) ) +                       &
                                         ( alpha_attack_i - alpha_attack_tab(ialpha-1) ) /             &
                                         ( alpha_attack_tab(ialpha) - alpha_attack_tab(ialpha-1) ) *   &
                                         ( weight_a * turb_cd_sel1(ialpha) +                           &
                                           weight_b * turb_cd_sel2(ialpha) )
    
           ENDDO   ! end loop over angles of attack
    
        ENDDO   ! end loop over radius
    
    
     END SUBROUTINE wtm_read_blade_tables
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> The projection matrix for the coordinate system of therotor disc in respect to the yaw and tilt
    !> angle of the rotor is calculated
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_rotate_rotor( inot )
    
    
        IMPLICIT NONE
    
        INTEGER(iwp) :: inot
    !
    !-- Calculation of the rotation matrix for the application of the tilt to the rotors:
        rot_eigen_rad(1) = SIN( yaw_angle(inot) )   ! x-component of the radial eigenvector
        rot_eigen_rad(2) = COS( yaw_angle(inot) )   ! y-component of the radial eigenvector
        rot_eigen_rad(3) = 0.0_wp                   ! z-component of the radial eigenvector
    
        rot_eigen_azi(1) = 0.0_wp                   ! x-component of the azimuth eigenvector
        rot_eigen_azi(2) = 0.0_wp                   ! y-component of the azimuth eigenvector
        rot_eigen_azi(3) = 1.0_wp                   ! z-component of the azimuth eigenvector
    
        rot_eigen_nor(1) =  COS( yaw_angle(inot) )  ! x-component of the normal eigenvector
        rot_eigen_nor(2) = -SIN( yaw_angle(inot) )  ! y-component of the normal eigenvector
        rot_eigen_nor(3) = 0.0_wp                   ! z-component of the normal eigenvector
    
    !
    !-- Calculation of the coordinate transformation matrix to apply a tilt to the rotor.
    !-- If tilt = 0, rot_coord_trans is a unit matrix.
        rot_coord_trans(inot,1,1) = rot_eigen_rad(1)**2                   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) ) + COS( tilt_angle )
        rot_coord_trans(inot,1,2) = rot_eigen_rad(1) * rot_eigen_rad(2)   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) )              -                      &
                                    rot_eigen_rad(3) * SIN( tilt_angle )
        rot_coord_trans(inot,1,3) = rot_eigen_rad(1) * rot_eigen_rad(3)   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) )              +                      &
                                    rot_eigen_rad(2) * SIN( tilt_angle )
        rot_coord_trans(inot,2,1) = rot_eigen_rad(2) * rot_eigen_rad(1)   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) )              +                      &
                                    rot_eigen_rad(3) * SIN( tilt_angle )
        rot_coord_trans(inot,2,2) = rot_eigen_rad(2)**2                   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) ) + COS( tilt_angle )
        rot_coord_trans(inot,2,3) = rot_eigen_rad(2) * rot_eigen_rad(3)   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) )              -                      &
                                    rot_eigen_rad(1) * SIN( tilt_angle )
        rot_coord_trans(inot,3,1) = rot_eigen_rad(3) * rot_eigen_rad(1)   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) )              -                      &
                                    rot_eigen_rad(2) * SIN( tilt_angle )
        rot_coord_trans(inot,3,2) = rot_eigen_rad(3) * rot_eigen_rad(2)   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) )              +                      &
                                    rot_eigen_rad(1) * SIN( tilt_angle )
        rot_coord_trans(inot,3,3) = rot_eigen_rad(3)**2                   *                            &
                                    ( 1.0_wp - COS( tilt_angle ) ) + COS( tilt_angle )
    
    !
    !-- Vectors for the Transformation of forces from the rotor's spheric coordinate system to the
    !-- cartesian coordinate system:
        rotx(inot,:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_nor )
        roty(inot,:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_rad )
        rotz(inot,:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_azi )
    
     END SUBROUTINE wtm_rotate_rotor
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Calculation of the forces generated by the wind turbine
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_forces
    
    
        IMPLICIT NONE
    
    
        INTEGER(iwp) ::  i, j, k       !< loop indices
        INTEGER(iwp) ::  ii, jj, kk    !<
        INTEGER(iwp) ::  inot          !< turbine loop index (turbine id)
        INTEGER(iwp) ::  iialpha, iir  !<
        INTEGER(iwp) ::  rseg          !<
        INTEGER(iwp) ::  ring          !<
    
        REAL(wp)     ::  flag              !< flag to mask topography grid points
        REAL(wp)     ::  sin_rot, cos_rot  !<
        REAL(wp)     ::  sin_yaw, cos_yaw  !<
    
        REAL(wp) ::  aa, bb, cc, dd  !< interpolation distances
        REAL(wp) ::  gg              !< interpolation volume var
    
        REAL(wp) ::  dist_u_3d, dist_v_3d, dist_w_3d  !< smearing distances
    
    
    !
    !-- Variables for pitch control:
        INTEGER(iwp), DIMENSION(1) ::  lct = 0
    
        LOGICAL ::  pitch_sw = .FALSE.
    
        REAL(wp), DIMENSION(1)     ::  rad_d = 0.0_wp
    
        REAL(wp) :: tl_factor  !< factor for tip loss correction
    
    
        CALL cpu_log( log_point_s(61), 'wtm_forces', 'start' )
    
    
        IF ( time_since_reference_point >= time_turbine_on )   THEN
    
    !
    !--    Set forces to zero for each new time step:
           thrust(:,:,:)         = 0.0_wp
           torque_y(:,:,:)       = 0.0_wp
           torque_z(:,:,:)       = 0.0_wp
           torque_total(:)       = 0.0_wp
           rot_tend_x(:,:,:)     = 0.0_wp
           rot_tend_y(:,:,:)     = 0.0_wp
           rot_tend_z(:,:,:)     = 0.0_wp
           thrust_rotor(:)       = 0.0_wp
    !
    !--    Loop over number of turbines:
           DO  inot = 1, n_turbines
    
              cos_yaw = COS(yaw_angle(inot))
              sin_yaw = SIN(yaw_angle(inot))
    !
    !--       Loop over rings of each turbine:
              !$OMP PARALLEL PRIVATE (ring, rseg, thrust_seg, torque_seg_y, torque_seg_z, sin_rot,  &
              !$OMP&                  cos_rot, re, rbx, rby, rbz, ii, jj, kk, aa, bb, cc, dd, gg)
              !$OMP DO
              DO  ring = 1, nrings(inot)
    
                 thrust_seg(:)   = 0.0_wp
                 torque_seg_y(:) = 0.0_wp
                 torque_seg_z(:) = 0.0_wp
    !
    !--          Determine distance between each ring (center) and the hub:
                 cur_r = (ring - 0.5_wp) * delta_r(inot)
    
    !
    !--          Loop over segments of each ring of each turbine:
                 DO  rseg = 1, nsegs(ring,inot)
    !
    !--             !----------------------------------------------------------------------------------!
    !--             !-- Determine coordinates of the ring segments                                   --!
    !--             !----------------------------------------------------------------------------------!
    !
    !--             Determine angle of ring segment towards zero degree angle of rotor system
    !--             (at zero degree rotor direction vectors aligned with y-axis):
                    phi_rotor = rseg * 2.0_wp * pi / nsegs(ring,inot)
                    cos_rot   = COS( phi_rotor )
                    sin_rot   = SIN( phi_rotor )
    
    !--             Now the direction vectors can be determined with respect to the yaw and tilt angle:
                    re(1) = cos_rot * sin_yaw
                    re(2) = cos_rot * cos_yaw
                    re(3) = sin_rot
    
                    rote = MATMUL( rot_coord_trans(inot,:,:), re )
    !
    !--             Coordinates of the single segments (center points):
                    rbx(ring,rseg) = hub_x(inot) + cur_r * rote(1)
                    rby(ring,rseg) = hub_y(inot) + cur_r * rote(2)
                    rbz(ring,rseg) = hub_z(inot) + cur_r * rote(3)
    
    !--             !----------------------------------------------------------------------------------!
    !--             !-- Interpolation of the velocity components from the cartesian grid point to    --!
    !--             !-- the coordinates of each ring segment (follows a method used in the           --!
    !--             !-- particle model)                                                              --!
    !--             !----------------------------------------------------------------------------------!
    
                    u_int(inot,ring,rseg)     = 0.0_wp
                    u_int_1_l(inot,ring,rseg) = 0.0_wp
    
                    v_int(inot,ring,rseg)     = 0.0_wp
                    v_int_1_l(inot,ring,rseg) = 0.0_wp
    
                    w_int(inot,ring,rseg)     = 0.0_wp
                    w_int_1_l(inot,ring,rseg) = 0.0_wp
    
    !
    !--             Interpolation of the u-component:
                    ii =   rbx(ring,rseg) * ddx
                    jj = ( rby(ring,rseg) - 0.5_wp * dy ) * ddy
                    kk = ( rbz(ring,rseg) + 0.5_wp * dz(1) ) / dz(1)
    !
    !--             Interpolate only if all required information is available on the current PE:
                    IF ( ( ii >= nxl )  .AND.  ( ii <= nxr ) )  THEN
                       IF ( ( jj >= nys )  .AND.  ( jj <= nyn ) )  THEN
    
                          aa = ( ( ii + 1          ) * dx - rbx(ring,rseg) ) *                         &
                               ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                          bb = ( rbx(ring,rseg) - ii * dx )                  *                         &
                               ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                          cc = ( ( ii+1            ) * dx - rbx(ring,rseg) ) *                         &
                               ( rby(ring,rseg) - ( jj + 0.5_wp ) * dy )
                          dd = ( rbx(ring,rseg) -              ii * dx )     *                         &
                               ( rby(ring,rseg) - ( jj + 0.5_wp ) * dy )
                          gg = dx * dy
    
                          u_int_l = ( aa * u(kk,jj,ii) + bb * u(kk,jj,ii+1)   + cc * u(kk,jj+1,ii)     &
                                    +  dd * u(kk,jj+1,ii+1) ) / gg
    
                          u_int_u = ( aa * u(kk+1,jj,ii) + bb * u(kk+1,jj,ii+1) + cc * u(kk+1,jj+1,ii) &
                                    + dd * u(kk+1,jj+1,ii+1) ) / gg
    
                          u_int_1_l(inot,ring,rseg) = u_int_l + ( rbz(ring,rseg) - zu(kk) ) / dz(1) *  &
                                                    ( u_int_u - u_int_l )
    
                       ELSE
                          u_int_1_l(inot,ring,rseg) = 0.0_wp
                       ENDIF
    
                    ELSE
                       u_int_1_l(inot,ring,rseg) = 0.0_wp
                    ENDIF
    
    
    !
    !--             Interpolation of the v-component:
                    ii = ( rbx(ring,rseg) - 0.5_wp * dx ) * ddx
                    jj =   rby(ring,rseg)                 * ddy
                    kk = ( rbz(ring,rseg) + 0.5_wp * dz(1) ) / dz(1)
    !
    !--             Interpolate only if all required information is available on the current PE:
                    IF ( ( ii >= nxl )  .AND.  ( ii <= nxr ) )  THEN
                       IF ( ( jj >= nys )  .AND.  ( jj <= nyn ) )  THEN
    
                          aa = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *                         &
                               ( ( jj + 1 )          * dy - rby(ring,rseg) )
                          bb = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *                         &
                               ( ( jj + 1 ) * dy          - rby(ring,rseg) )
                          cc = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *                         &
                               ( rby(ring,rseg)           -        jj * dy )
                          dd = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *                         &
                               ( rby(ring,rseg)           -        jj * dy )
                          gg = dx * dy
    
                          v_int_l = ( aa * v(kk,jj,ii) + bb * v(kk,jj,ii+1) + cc * v(kk,jj+1,ii)       &
                                    + dd * v(kk,jj+1,ii+1) ) / gg
    
                          v_int_u = ( aa * v(kk+1,jj,ii) + bb * v(kk+1,jj,ii+1) + cc * v(kk+1,jj+1,ii) &
                                    + dd * v(kk+1,jj+1,ii+1) ) / gg
    
                          v_int_1_l(inot,ring,rseg) = v_int_l + ( rbz(ring,rseg) - zu(kk) ) / dz(1) *  &
                                                    ( v_int_u - v_int_l )
    
                       ELSE
                          v_int_1_l(inot,ring,rseg) = 0.0_wp
                       ENDIF
    
                    ELSE
                       v_int_1_l(inot,ring,rseg) = 0.0_wp
                    ENDIF
    
    
    !
    !--             Interpolation of the w-component:
                    ii = ( rbx(ring,rseg) - 0.5_wp * dx ) * ddx
                    jj = ( rby(ring,rseg) - 0.5_wp * dy ) * ddy
                    kk =   rbz(ring,rseg)                 / dz(1)
    !
    !--             Interpolate only if all required information is available on the current PE:
                    IF ( ( ii >= nxl )  .AND.  ( ii <= nxr ) )  THEN
                       IF ( ( jj >= nys )  .AND.  ( jj <= nyn ) )  THEN
    
                          aa = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *                         &
                               ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                          bb = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *                         &
                               ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                          cc = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *                         &
                               ( rby(ring,rseg)     - ( jj + 0.5_wp ) * dy )
                          dd = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *                         &
                               ( rby(ring,rseg)     - ( jj + 0.5_wp ) * dy )
                          gg = dx * dy
    
                          w_int_l = ( aa * w(kk,jj,ii) + bb * w(kk,jj,ii+1) + cc * w(kk,jj+1,ii)       &
                                    + dd * w(kk,jj+1,ii+1) ) / gg
    
                          w_int_u = ( aa * w(kk+1,jj,ii) + bb * w(kk+1,jj,ii+1) + cc * w(kk+1,jj+1,ii) &
                                    + dd * w(kk+1,jj+1,ii+1) ) / gg
    
                          w_int_1_l(inot,ring,rseg) = w_int_l + ( rbz(ring,rseg) - zw(kk) ) / dz(1) *  &
                                                    ( w_int_u - w_int_l )
                       ELSE
                          w_int_1_l(inot,ring,rseg) = 0.0_wp
                       ENDIF
    
                    ELSE
                       w_int_1_l(inot,ring,rseg) = 0.0_wp
                    ENDIF
    
                 ENDDO
    
              ENDDO
              !$OMP END PARALLEL
    
           ENDDO
    
    !
    !--    Exchange between PEs (information required on each PE):
    #if defined( __parallel )
           CALL MPI_ALLREDUCE( u_int_1_l, u_int, n_turbines * MAXVAL(nrings) *                         &
                               MAXVAL(nsegs), MPI_REAL, MPI_SUM, comm2d, ierr )
           CALL MPI_ALLREDUCE( v_int_1_l, v_int, n_turbines * MAXVAL(nrings) *                         &
                               MAXVAL(nsegs), MPI_REAL, MPI_SUM, comm2d, ierr )
           CALL MPI_ALLREDUCE( w_int_1_l, w_int, n_turbines * MAXVAL(nrings) *                         &
                               MAXVAL(nsegs), MPI_REAL, MPI_SUM, comm2d, ierr )
    #else
           u_int = u_int_1_l
           v_int = v_int_1_l
           w_int = w_int_1_l
    #endif
    
    
    !
    !--    Loop over number of turbines:
           DO  inot = 1, n_turbines
    pit_loop: DO
    
              IF ( pitch_sw )  THEN
                 torque_total(inot) = 0.0_wp
                 thrust_rotor(inot) = 0.0_wp
                 pitch_angle(inot)    = pitch_angle(inot) + 0.25_wp
    !              IF ( myid == 0 ) PRINT*, 'Pitch', inot, pitch_angle(inot)
              ELSE
                 cos_yaw = COS( yaw_angle(inot) )
                 sin_yaw = SIN( yaw_angle(inot) )
                 IF ( pitch_control )  THEN
                    pitch_angle(inot) = MAX( pitch_angle_old(inot) - pitch_rate * dt_3d , 0.0_wp )
                 ENDIF
              ENDIF
    
    !
    !--       Loop over rings of each turbine:
              !$OMP PARALLEL PRIVATE (ring, rseg, sin_rot, cos_rot, re, rea, ren, rote, rota, rotn, &
              !$OMP&                  vtheta, phi_rel, lct, rad_d, alpha_attack, vrel,              &
              !$OMP&                  chord, iialpha, iir, turb_cl, tl_factor, thrust_seg,          &
              !$OMP&                  torque_seg_y, turb_cd, torque_seg_z, thrust_ring,             &
              !$OMP&                  torque_ring_y, torque_ring_z)
              !$OMP DO
              DO  ring = 1, nrings(inot)
    !
    !--          Determine distance between each ring (center) and the hub:
                 cur_r = ( ring - 0.5_wp ) * delta_r(inot)
    !
    !--          Loop over segments of each ring of each turbine:
                 DO  rseg = 1, nsegs(ring,inot)
    !
    !--             Determine angle of ring segment towards zero degree angle of rotor system
    !--             (at zero degree rotor direction vectors aligned with y-axis):
                    phi_rotor = rseg * 2.0_wp * pi / nsegs(ring,inot)
                    cos_rot   = COS( phi_rotor )
                    sin_rot   = SIN( phi_rotor )
    !
    !--             Now the direction vectors can be determined with respect to the yaw and tilt angle:
                    re(1) = cos_rot * sin_yaw
                    re(2) = cos_rot * cos_yaw
                    re(3) = sin_rot
    
    !               The current unit vector in azimuthal direction:
                    rea(1) = - sin_rot * sin_yaw
                    rea(2) = - sin_rot * cos_yaw
                    rea(3) =   cos_rot
    
    !
    !--             To respect the yawing angle for the calculations of velocities and forces the
    !--             unit vectors perpendicular to the rotor area in direction of the positive yaw
    !--             angle are defined:
                    ren(1) =   cos_yaw
                    ren(2) = - sin_yaw
                    ren(3) = 0.0_wp
    !
    !--             Multiplication with the coordinate transformation matrix gives the final unit
    !--             vector with consideration of the rotor tilt:
                    rote = MATMUL( rot_coord_trans(inot,:,:), re )
                    rota = MATMUL( rot_coord_trans(inot,:,:), rea )
                    rotn = MATMUL( rot_coord_trans(inot,:,:), ren )
    !
    !--             Coordinates of the single segments (center points):
                    rbx(ring,rseg) = hub_x(inot) + cur_r * rote(1)
    
                    rby(ring,rseg) = hub_y(inot) + cur_r * rote(2)
    
                    rbz(ring,rseg) = hub_z(inot) + cur_r * rote(3)
    
    !
    !--             !----------------------------------------------------------------------------------!
    !--             !-- Calculation of various angles and relative velocities                        --!
    !--             !----------------------------------------------------------------------------------!
    !
    !--             In the following the 3D-velocity field is projected its components perpendicular
    !--             and parallel to the rotor area.
    !--             The calculation of forces will be done in the rotor-coordinates y' and z.
    !--             The yaw angle will be reintroduced when the force is applied on the hydrodynamic
    !--             equations.
    !
    !--             Projection of the xy-velocities relative to the rotor area:
    !
    !--             Velocity perpendicular to the rotor area:
                    u_rot = u_int(inot,ring,rseg) * rotn(1) + v_int(inot,ring,rseg)*rotn(2) +          &
                            w_int(inot,ring,rseg)*rotn(3)
    !
    
    !--             Projection of the 3D-velocity vector in the azimuthal direction:
                    vtheta(rseg) = rota(1) * u_int(inot,ring,rseg) + rota(2) * v_int(inot,ring,rseg) + &
                                   rota(3) * w_int(inot,ring,rseg)
    !
    
    !--             Determination of the angle phi_rel between the rotor plane and the direction of the
    !--             flow relative to the rotor:
                    phi_rel(rseg) = ATAN2( u_rot , ( rotor_speed(inot) * cur_r - vtheta(rseg) ) )
    
    !
    !--             Interpolation of the local pitch angle from tabulated values to the current
    !--             radial position:
                    lct = minloc( ABS( cur_r-lrd ) )
                    rad_d = cur_r-lrd(lct)
    
                    IF ( cur_r == 0.0_wp )  THEN
                       alpha_attack(rseg) = 0.0_wp
    
                    ELSE IF ( cur_r >= lrd(size(ard)) )  THEN
                       alpha_attack(rseg) = ( ard(size(ard) ) + ard(size(ard) -1 ) ) / 2.0_wp
    
                    ELSE
                       alpha_attack(rseg) = ( ard( lct(1) ) * ( ( lrd( lct(1) + 1 ) - cur_r ) /        &
                                            ( lrd( lct(1) + 1 ) - lrd( lct(1) ) ) ) ) +                &
                                            ( ard( lct(1) + 1 ) * ( ( cur_r - lrd( lct(1) ) ) /        &
                                            ( lrd( lct(1) + 1 ) - lrd( lct(1) ) ) ) )
                    ENDIF
    
    !
    !--             In Fortran radian instead of degree is used as unit for all angles.
    !--             Therefore, a transformation from angles given in degree to angles given in radian
    !--             is necessary here:
                    alpha_attack(rseg) = alpha_attack(rseg) * ( ( 2.0_wp * pi ) / 360.0_wp )
    !
    !--             Substraction of the local pitch angle to obtain the local angle of attack:
                    alpha_attack(rseg) = phi_rel(rseg) - alpha_attack(rseg)
    !
    !--             Preliminary transformation back from angles given in radian to angles given in
    !--             degree:
                    alpha_attack(rseg) = alpha_attack(rseg) * ( 360.0_wp / ( 2.0_wp * pi ) )
    !
    !--             Correct with collective pitch angle:
                    alpha_attack(rseg) = alpha_attack(rseg) - pitch_angle(inot)
    
    !
    !--             Determination of the magnitude of the flow velocity relative to the rotor:
                    vrel(rseg) = SQRT( u_rot**2 + ( rotor_speed(inot) * cur_r - vtheta(rseg) )**2 )
    
    !
    !--             !----------------------------------------------------------------------------------!
    !--             !-- Interpolation of chord as well as lift and drag                              --!
    !--             !-- coefficients from tabulated values                                           --!
    !--             !----------------------------------------------------------------------------------!
    
    !
    !--             Interpolation of the chord_length from tabulated values to the current radial
    !--             position:
                    IF ( cur_r == 0.0_wp )  THEN
                       chord(rseg) = 0.0_wp
    
                    ELSE IF ( cur_r >= lrd( size(crd) ) )  THEN
                       chord(rseg) = ( crd( size(crd) ) + ard( size(crd) - 1 ) ) / 2.0_wp
    
                    ELSE
                       chord(rseg) = ( crd( lct(1) ) * ( ( lrd( lct(1) + 1 ) - cur_r ) /               &
                                     ( lrd( lct(1) + 1 ) - lrd( lct(1) ) ) ) ) + ( crd( lct(1) + 1 )   &
                                     * ( ( cur_r - lrd( lct(1) ) ) / ( lrd( lct(1) + 1 ) -             &
                                     lrd( lct(1) ) ) ) )
                    ENDIF
    
    !
    !--             Determine index of current angle of attack, needed for finding the appropriate
    !--             interpolated values of the lift and drag coefficients
    !--             (-180.0 degrees = 1, +180.0 degrees = 3601, so one index every 0.1 degrees):
                    iialpha = CEILING( ( alpha_attack(rseg) + 180.0_wp )                               &
                            * ( 1.0_wp / accu_cl_cd_tab ) ) + 1.0_wp
    !
    !--             Determine index of current radial position, needed for finding the appropriate
    !--             interpolated values of the lift and drag coefficients (one index every 0.1 m):
                    iir = CEILING( cur_r * 10.0_wp )
    !
    !--             Read in interpolated values of the lift and drag coefficients for the current
    !--             radial position and angle of attack:
                    turb_cl(rseg) = turb_cl_tab(iialpha,iir)
                    turb_cd(rseg) = turb_cd_tab(iialpha,iir)
    
    !
    !--             Final transformation back from angles given in degree to angles given in radian:
                    alpha_attack(rseg) = alpha_attack(rseg) * ( ( 2.0_wp * pi ) / 360.0_wp )
    
                    IF ( tip_loss_correction )  THEN
    !
    !--               Tip loss correction following Schito.
    !--               Schito applies the tip loss correction only to the lift force.
    !--               Therefore, the tip loss correction is only applied to the lift coefficient and
    !--               not to the drag coefficient in our case.
                      IF ( phi_rel(rseg) == 0.0_wp )  THEN
                         tl_factor = 1.0_wp
    
                      ELSE
                         tl_factor = ( 2.0 / pi ) *                                                    &
                              ACOS( EXP( -1.0 * ( 3.0 * ( rotor_radius(inot) - cur_r ) /               &
                              ( 2.0 * cur_r * ABS( SIN( phi_rel(rseg) ) ) ) ) ) )
                      ENDIF
    
                      turb_cl(rseg)  = tl_factor * turb_cl(rseg)
    
                    ENDIF
    !
    !--             !----------------------------------------------------------------------------------!
    !--             !-- Calculation of the forces                                                    --!
    !--             !----------------------------------------------------------------------------------!
    
    !
    !--             Calculate the pre_factor for the thrust and torque forces:
                    pre_factor = 0.5_wp * ( vrel(rseg)**2 ) * 3.0_wp * chord(rseg) * delta_r(inot)     &
                               / nsegs(ring,inot)
    
    !
    !--             Calculate the thrust force (x-component of the total force) for each ring segment:
                    thrust_seg(rseg) = pre_factor * ( turb_cl(rseg) * COS (phi_rel(rseg) ) +           &
                                       turb_cd(rseg) * SIN( phi_rel(rseg) ) )
    
    !
    !--             Determination of the second of the additional forces acting on the flow in the
    !--             azimuthal direction: force vector as basis for torque (torque itself would be the
    !--             vector product of the radius vector and the force vector):
                    torque_seg = pre_factor * ( turb_cl(rseg) * SIN (phi_rel(rseg) ) -                 &
                                 turb_cd(rseg) * COS( phi_rel(rseg) ) )
    !
    !--             Decomposition of the force vector into two parts: One acting along the
    !--             y-direction and one acting along the z-direction of the rotor coordinate system:
                    torque_seg_y(rseg) = -torque_seg * sin_rot
                    torque_seg_z(rseg) =  torque_seg * cos_rot
    
    !
    !--             Add the segment thrust to the thrust of the whole rotor:
                    !$OMP CRITICAL
                    thrust_rotor(inot) = thrust_rotor(inot) + thrust_seg(rseg)
    
    
                    torque_total(inot) = torque_total(inot) + (torque_seg * cur_r)
                    !$OMP END CRITICAL
    
                 ENDDO   !-- end of loop over ring segments
    
    !
    !--          Restore the forces into arrays containing all the segments of each ring:
                 thrust_ring(ring,:)   = thrust_seg(:)
                 torque_ring_y(ring,:) = torque_seg_y(:)
                 torque_ring_z(ring,:) = torque_seg_z(:)
    
    
              ENDDO   !-- end of loop over rings
              !$OMP END PARALLEL
    
    
              CALL cpu_log( log_point_s(62), 'wtm_controller', 'start' )
    
    
              IF ( speed_control )  THEN
    !
    !--          Calculation of the current generator speed for rotor speed control:
    
    !
    !--          The acceleration of the rotor speed is calculated from the force balance of the
    !--          accelerating torque and the torque of the rotating rotor and generator:
                 om_rate = ( torque_total(inot) * air_density * gear_efficiency -                      &
                             gear_ratio * torque_gen_old(inot) ) / ( rotor_inertia +                   &
                             gear_ratio * gear_ratio * generator_inertia ) * dt_3d
    
    !
    !--          The generator speed is given by the product of gear gear_ratio and rotor speed:
                 generator_speed(inot) = gear_ratio * ( rotor_speed(inot) + om_rate )
    
              ENDIF
    
              IF ( pitch_control )  THEN
    
    !
    !--          If the current generator speed is above rated, the pitch is not saturated and the
    !--          change from the last time step is within the maximum pitch rate, then the pitch loop
    !--          is repeated with a pitch gain:
                 IF ( ( generator_speed(inot)  > generator_speed_rated )  .AND.                        &
                      ( pitch_angle(inot) < 25.0_wp )  .AND.                                           &
                      ( pitch_angle(inot) < pitch_angle_old(inot) + pitch_rate * dt_3d ) )  THEN
                    pitch_sw = .TRUE.
    !
    !--             Go back to beginning of pit_loop:
                    CYCLE pit_loop
                 ENDIF
    
    !
    !--          The current pitch is saved for the next time step:
                 pitch_angle_old(inot) = pitch_angle(inot)
                 pitch_sw = .FALSE.
              ENDIF
              EXIT pit_loop
           ENDDO pit_loop ! Recursive pitch control loop
    
    
    !
    !--       Call the rotor speed controller:
              IF ( speed_control )  THEN
    !
    !--          Find processor at i_hub, j_hub:
                 IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )  THEN
                    IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )  THEN
                       CALL wtm_speed_control( inot )
                    ENDIF
                 ENDIF
    
              ENDIF
    
    
              CALL cpu_log( log_point_s(62), 'wtm_controller', 'stop' )
    
              CALL cpu_log( log_point_s(63), 'wtm_smearing', 'start' )
    
    
    !--       !----------------------------------------------------------------------------------------!
    !--       !--                  Regularization kernel                                             --!
    !--       !-- Smearing of the forces and interpolation to cartesian grid                         --!
    !--       !----------------------------------------------------------------------------------------!
    !
    !--       The aerodynamic blade forces need to be distributed smoothly on several mesh points in
    !--       order to avoid singular behaviour.
    !
    !--       Summation over sum of weighted forces. The weighting factor (calculated in user_init)
    !--       includes information on the distance between the center of the grid cell and the rotor
    !--       segment under consideration.
    !
    !--       To save computing time, apply smearing only for the relevant part of the model domain:
    !
    !--
    !--       Calculation of the boundaries:
              i_smear(inot) = CEILING( ( rotor_radius(inot) * ABS( roty(inot,1) ) + eps_min ) / dx )
              j_smear(inot) = CEILING( ( rotor_radius(inot) * ABS( roty(inot,2) ) + eps_min ) / dy )
    
              !$OMP PARALLEL PRIVATE (i, j, k, ring, rseg, flag, dist_u_3d, dist_v_3d, dist_w_3d)
              !$OMP DO
              DO  i = MAX( nxl, i_hub(inot) - i_smear(inot) ), MIN( nxr, i_hub(inot) + i_smear(inot) )
                 DO  j = MAX( nys, j_hub(inot) - j_smear(inot) ), MIN( nyn, j_hub(inot) + j_smear(inot) )
    !                 DO  k = MAX( nzb_u_inner(j,i)+1, k_hub(inot) - k_smear(inot) ),                   &
    !                              k_hub(inot) + k_smear(inot)
                    DO  k = nzb + 1, k_hub(inot) + k_smear(inot)
                       DO  ring = 1, nrings(inot)
                          DO  rseg = 1, nsegs(ring,inot)
    !
    !--                      Predetermine flag to mask topography:
                             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
    
    !
    !--                      Determine the square of the distance between the current grid point and
    !--                      each rotor area segment:
                             dist_u_3d = ( i * dx               - rbx(ring,rseg) )**2 +                &
                                         ( j * dy + 0.5_wp * dy - rby(ring,rseg) )**2 +                &
                                         ( k * dz(1) - 0.5_wp * dz(1) - rbz(ring,rseg) )**2
    
                             dist_v_3d = ( i * dx + 0.5_wp * dx - rbx(ring,rseg) )**2 +                &
                                         ( j * dy               - rby(ring,rseg) )**2 +                &
                                         ( k * dz(1) - 0.5_wp * dz(1) - rbz(ring,rseg) )**2
    
                             dist_w_3d = ( i * dx + 0.5_wp * dx - rbx(ring,rseg) )**2 +                &
                                         ( j * dy + 0.5_wp * dy - rby(ring,rseg) )**2 +                &
                                         ( k * dz(1)               - rbz(ring,rseg) )**2
    
    !
    !--                      3D-smearing of the forces with a polynomial function (much faster than
    !--                      the old Gaussian function), using some parameters that have been
    !--                      calculated in user_init. The function is only similar to Gaussian
    !--                      function for squared distances <= eps_min2:
                             IF ( dist_u_3d <= eps_min2 )  THEN
                                thrust(k,j,i) = thrust(k,j,i) + thrust_ring(ring,rseg) *               &
                                                ( ( pol_a * dist_u_3d - pol_b ) *                      &
                                                 dist_u_3d + 1.0_wp ) * eps_factor * flag
                             ENDIF
    
                             IF ( dist_v_3d <= eps_min2 )  THEN
                                torque_y(k,j,i) = torque_y(k,j,i) + torque_ring_y(ring,rseg) *         &
                                                  ( ( pol_a * dist_v_3d - pol_b ) *                    &
                                                   dist_v_3d + 1.0_wp ) * eps_factor * flag
                             ENDIF
    
                             IF ( dist_w_3d <= eps_min2 )  THEN
                                torque_z(k,j,i) = torque_z(k,j,i) + torque_ring_z(ring,rseg) *         &
                                                  ( ( pol_a * dist_w_3d - pol_b ) *                    &
                                                   dist_w_3d + 1.0_wp ) * eps_factor * flag
                             ENDIF
    
                          ENDDO  ! End of loop over rseg
                       ENDDO     ! End of loop over ring
    
    !
    !--                Rotation of force components:
                       rot_tend_x(k,j,i) = rot_tend_x(k,j,i) + ( thrust(k,j,i) * rotx(inot,1) +        &
                                           torque_y(k,j,i) * roty(inot,1) + torque_z(k,j,i) *          &
                                           rotz(inot,1) ) * flag
    
                       rot_tend_y(k,j,i) = rot_tend_y(k,j,i) + ( thrust(k,j,i) * rotx(inot,2) +        &
                                           torque_y(k,j,i) * roty(inot,2) + torque_z(k,j,i) *          &
                                           rotz(inot,2) ) * flag
    
                       rot_tend_z(k,j,i) = rot_tend_z(k,j,i) + ( thrust(k,j,i) * rotx(inot,3) +        &
                                           torque_y(k,j,i) * roty(inot,3) +  torque_z(k,j,i) *         &
                                           rotz(inot,3) ) * flag
    
                    ENDDO  ! End of loop over k
                 ENDDO     ! End of loop over j
              ENDDO        ! End of loop over i
              !$OMP END PARALLEL
    
              CALL cpu_log( log_point_s(63), 'wtm_smearing', 'stop' )
    
           ENDDO                  !-- end of loop over turbines
    
    
           IF ( yaw_control )  THEN
    !
    !--       Allocate arrays for yaw control at first call. Can't be allocated before dt_3d is set.
              IF ( start_up )  THEN
                 wdlon= MAX( 1 , NINT( 30.0_wp / dt_3d ) )  ! 30s running mean array
                 ALLOCATE( wd30(1:n_turbines,1:WDLON) )
                 wd30 = 999.0_wp                             ! Set to dummy value
                 ALLOCATE( wd30_l(1:WDLON) )
    
                 wdsho = MAX( 1 , NINT( 2.0_wp / dt_3d ) )   ! 2s running mean array
                 ALLOCATE( wd2(1:n_turbines,1:wdsho) )
                 wd2 = 999.0_wp                              ! Set to dummy value
                 ALLOCATE( wd2_l(1:wdsho) )
                 start_up = .FALSE.
              ENDIF
    
    !
    !--       Calculate the inflow wind speed:
    !--
    !--       Loop over number of turbines:
              DO  inot = 1, n_turbines
    !
    !--          Find processor at i_hub, j_hub:
                 IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )  THEN
                    IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )  THEN
                       u_inflow_l(inot) = u(k_hub(inot),j_hub(inot),i_hub(inot))
                       wdir_l(inot)     = -1.0_wp * ATAN2( 0.5_wp *                                    &
                                        ( v(k_hub(inot), j_hub(inot), i_hub(inot) + 1) +               &
                                          v(k_hub(inot), j_hub(inot), i_hub(inot)) ) , 0.5_wp *        &
                                        ( u(k_hub(inot), j_hub(inot) + 1, i_hub(inot)) +               &
                                          u(k_hub(inot), j_hub(inot), i_hub(inot)) ) )
    
                       CALL wtm_yawcontrol( inot )
    
                       yaw_angle_l(inot) = yaw_angle(inot)
    
                    ENDIF
                 ENDIF
    
              ENDDO  ! end of loop over turbines
    
    !
    !--       Transfer of information to the other cpus:
    #if defined( __parallel )
              CALL MPI_ALLREDUCE( u_inflow_l, u_inflow, n_turbines, MPI_REAL,                          &
                                  MPI_SUM, comm2d, ierr )
              CALL MPI_ALLREDUCE( wdir_l, wdir, n_turbines, MPI_REAL, MPI_SUM,                         &
                                  comm2d, ierr )
              CALL MPI_ALLREDUCE( yaw_angle_l, yaw_angle, n_turbines, MPI_REAL,                        &
                                  MPI_SUM, comm2d, ierr )
    #else
              u_inflow   = u_inflow_l
              wdir       = wdir_l
              yaw_angle  = yaw_angle_l
    
    #endif
              DO  inot = 1, n_turbines
    !
    !--          Update rotor orientation:
                 CALL wtm_rotate_rotor( inot )
    
              ENDDO ! End of loop over turbines
    
           ENDIF  ! end of yaw control
    
           IF ( speed_control )  THEN
    !
    !--       Transfer of information to the other cpus:
    !           CALL MPI_ALLREDUCE( generator_speed, generator_speed_old, n_turbines,                   &
    !                               MPI_REAL,MPI_SUM, comm2d, ierr )
    #if defined( __parallel )
              CALL MPI_ALLREDUCE( torque_gen, torque_gen_old, n_turbines,                              &
                                  MPI_REAL, MPI_SUM, comm2d, ierr )
              CALL MPI_ALLREDUCE( rotor_speed_l, rotor_speed, n_turbines,                              &
                                  MPI_REAL, MPI_SUM, comm2d, ierr )
              CALL MPI_ALLREDUCE( generator_speed_f, generator_speed_f_old, n_turbines,                &
                                  MPI_REAL, MPI_SUM, comm2d, ierr )
    #else
              torque_gen_old        = torque_gen
              rotor_speed           = rotor_speed_l
              generator_speed_f_old = generator_speed_f
    #endif
    
           ENDIF
    
        ENDIF
    
        CALL cpu_log( log_point_s(61), 'wtm_forces', 'stop' )
    
    
     END SUBROUTINE wtm_forces
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Yaw controller for the wind turbine model
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_yawcontrol( inot )
    
        USE kinds
    
        IMPLICIT NONE
    
        INTEGER(iwp)             :: inot
        INTEGER(iwp)             :: i_wd_30
    
        REAL(wp)                 :: missal
    
        i_wd_30 = 0_iwp
    
    !
    !-- The yaw controller computes a 30s running mean of the wind direction. If the difference
    !-- between turbine alignment and wind direction exceeds 5 degrees, the turbine is yawed. The
    !-- mechanism stops as soon as the 2s-running mean of the missalignment is smaller than 0.5
    !-- degrees. Attention: If the timestep during the simulation changes significantly the lengths
    !-- of the running means change and it does not correspond to 30s/2s anymore.
    !-- ! Needs to be modified for these situations !
    !-- For wind from the east, the averaging of the wind direction could cause problems and the yaw
    !-- controller is probably flawed. -> Routine for averaging needs to be improved!
    !
    !-- Check if turbine is not yawing:
        IF ( .NOT. doyaw(inot) )  THEN
    !
    !--    Write current wind direction into array:
           wd30_l    = wd30(inot,:)
           wd30_l    = CSHIFT( wd30_l, SHIFT = -1 )
           wd30_l(1) = wdir(inot)
    !
    !--    Check if array is full ( no more dummies ):
           IF ( .NOT. ANY( wd30_l == 999.) )  THEN
    
              missal = SUM( wd30_l ) / SIZE( wd30_l ) - yaw_angle(inot)
    !
    !--       Check if turbine is missaligned by more than yaw_misalignment_max:
              IF ( ABS( missal ) > yaw_misalignment_max )  THEN
    !
    !--          Check in which direction to yaw:
                 yawdir(inot) = SIGN( 1.0_wp, missal )
    !
    !--          Start yawing of turbine:
                 yaw_angle(inot) = yaw_angle(inot) + yawdir(inot) * yaw_speed * dt_3d
                 doyaw(inot)     = .TRUE.
                 wd30_l          = 999.  ! fill with dummies again
              ENDIF
           ENDIF
    
           wd30(inot,:) = wd30_l
    
    !
    !-- If turbine is already yawing:
    !-- Initialize 2 s running mean and yaw until the missalignment is smaller than
    !-- yaw_misalignment_min
    
        ELSE
    !
    !--    Initialize 2 s running mean:
    
           wd2_l = wd2(inot,:)
           wd2_l = CSHIFT( wd2_l, SHIFT = -1 )
           wd2_l(1) = wdir(inot)
    !
    !--    Check if array is full ( no more dummies ):
           IF ( .NOT. ANY( wd2_l == 999.0_wp ) ) THEN
    !
    !--       Calculate missalignment of turbine:
              missal = SUM( wd2_l - yaw_angle(inot) ) / SIZE( wd2_l )
    !
    !--       Check if missalignment is still larger than 0.5 degree and if the yaw direction is
    !--       still right:
              IF ( ( ABS( missal ) > yaw_misalignment_min )  .AND.                                     &
                   ( yawdir(inot) == SIGN( 1.0_wp, missal ) ) )  THEN
    !
    !--          Continue yawing:
                 yaw_angle(inot) = yaw_angle(inot) + yawdir(inot) * yaw_speed * dt_3d
    
              ELSE
    !
    !--          Stop yawing:
                 doyaw(inot) = .FALSE.
                 wd2_l       = 999.0_wp  ! fill with dummies again
              ENDIF
           ELSE
    !
    !--       Continue yawing:
              yaw_angle(inot) = yaw_angle(inot) + yawdir(inot) * yaw_speed * dt_3d
           ENDIF
    
           wd2(inot,:) = wd2_l
    
        ENDIF
    
     END SUBROUTINE wtm_yawcontrol
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Initialization of the speed control
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_init_speed_control
    
    
        IMPLICIT NONE
    
    !
    !-- If speed control is set, remaining variables and control_parameters for the control algorithm
    !-- are calculated:
    !
    !-- Calculate slope constant for region 15:
        slope15         = ( region_2_slope * region_2_min * region_2_min ) /                           &
                          ( region_2_min - region_15_min )
    !
    !-- Calculate upper limit of slipage region:
        vs_sysp         = generator_speed_rated / 1.1_wp
    !
    !-- Calculate slope of slipage region:
        region_25_slope = ( generator_power_rated / generator_speed_rated ) /                          &
                        ( generator_speed_rated - vs_sysp )
    !
    !-- Calculate lower limit of slipage region:
        region_25_min   = ( region_25_slope - SQRT( region_25_slope *                                  &
                          ( region_25_slope - 4.0_wp * region_2_slope * vs_sysp ) ) ) /                &
                          ( 2.0_wp * region_2_slope )
    !
    !-- Frequency for the simple low pass filter:
        fcorner      = 0.25_wp
    !
    !-- At the first timestep the torque is set to its maximum to prevent an overspeeding of the rotor:
        IF ( TRIM( initializing_actions ) /= 'read_restart_data' ) THEN
           torque_gen_old(:) = generator_torque_max
        ENDIF
    
     END SUBROUTINE wtm_init_speed_control
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Simple controller for the regulation of the rotor speed
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_speed_control( inot )
    
    
        IMPLICIT NONE
    
        INTEGER(iwp)::  inot
    
    
    !
    !-- The controller is based on the fortran script from Jonkman et al. 2009 "Definition of a 5 MW
    !-- Reference Wind Turbine for offshore system developement"
    
    !
    !-- The generator speed is filtered by a low pass filter  for the control of the generator torque:
        lp_coeff = EXP( -2.0_wp * 3.14_wp * dt_3d * fcorner)
        generator_speed_f(inot) = ( 1.0_wp - lp_coeff ) * generator_speed(inot) + lp_coeff *           &
                                   generator_speed_f_old(inot)
    
        IF ( generator_speed_f(inot) <= region_15_min )  THEN
    !
    !--    Region 1: Generator torque is set to zero to accelerate the rotor:
           torque_gen(inot) = 0
    
        ELSEIF ( generator_speed_f(inot) <= region_2_min )  THEN
    !
    !--    Region 1.5: Generator torque is increasing linearly with rotor speed:
           torque_gen(inot) = slope15 * ( generator_speed_f(inot) - region_15_min )
    
        ELSEIF ( generator_speed_f(inot) <= region_25_min )  THEN
    !
    !--    Region 2: Generator torque is increased by the square of the generator speed to keep the
    !--              TSR optimal:
           torque_gen(inot) = region_2_slope * generator_speed_f(inot) * generator_speed_f(inot)
    
        ELSEIF ( generator_speed_f(inot) < generator_speed_rated )  THEN
    !
    !--    Region 2.5: Slipage region between 2 and 3:
           torque_gen(inot) = region_25_slope * ( generator_speed_f(inot) - vs_sysp )
    
        ELSE
    !
    !--    Region 3: Generator torque is antiproportional to the rotor speed to keep the power
    !--              constant:
           torque_gen(inot) = generator_power_rated / generator_speed_f(inot)
    
        ENDIF
    !
    !-- Calculate torque rate and confine with a max:
        trq_rate = ( torque_gen(inot) - torque_gen_old(inot) ) / dt_3d
        trq_rate = MIN( MAX( trq_rate, -1.0_wp * generator_torque_rate_max ), generator_torque_rate_max )
    !
    !-- Calculate new gen torque and confine with max torque:
        torque_gen(inot) = torque_gen_old(inot) + trq_rate * dt_3d
        torque_gen(inot) = MIN( torque_gen(inot), generator_torque_max )
    !
    !-- Overwrite values for next timestep:
        rotor_speed_l(inot) = generator_speed(inot) / gear_ratio
    
    
     END SUBROUTINE wtm_speed_control
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Application of the additional forces generated by the wind turbine on the flow components
    !> (tendency terms)
    !> Call for all grid points
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_actions( location )
    
    
        CHARACTER(LEN=*) ::  location  !<
    
        INTEGER(iwp) ::  i  !< running index
        INTEGER(iwp) ::  j  !< running index
        INTEGER(iwp) ::  k  !< running index
    
    
        SELECT CASE ( location )
    
        CASE ( 'before_timestep' )
    
           CALL wtm_forces
    
        CASE ( 'u-tendency' )
    !
    
    !--    Apply the x-component of the force to the u-component of the flow:
           IF ( time_since_reference_point >= time_turbine_on )  THEN
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    DO  k = nzb + 1, MAXVAL( k_hub ) + MAXVAL( k_smear )
    !
    !--                Calculate the thrust generated by the nacelle and the tower:
                       tend_nac_x  = 0.5_wp * nac_cd_surf(k,j,i) * SIGN( u(k,j,i)**2 , u(k,j,i) )
                       tend_tow_x  = 0.5_wp * tow_cd_surf(k,j,i) * SIGN( u(k,j,i)**2 , u(k,j,i) )
                       tend(k,j,i) = tend(k,j,i) + ( - rot_tend_x(k,j,i) - tend_nac_x - tend_tow_x ) * &
                                     MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
    
        CASE ( 'v-tendency' )
    !
    !--    Apply the y-component of the force to the v-component of the flow:
           IF ( time_since_reference_point >= time_turbine_on )  THEN
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    DO  k = nzb+1, MAXVAL(k_hub) + MAXVAL(k_smear)
                       tend_nac_y  = 0.5_wp * nac_cd_surf(k,j,i) * SIGN( v(k,j,i)**2 , v(k,j,i) )
                       tend_tow_y  = 0.5_wp * tow_cd_surf(k,j,i) * SIGN( v(k,j,i)**2 , v(k,j,i) )
                       tend(k,j,i) = tend(k,j,i) + ( - rot_tend_y(k,j,i) - tend_nac_y - tend_tow_y ) * &
                                     MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
    
        CASE ( 'w-tendency' )
    !
    !--    Apply the z-component of the force to the w-component of the flow:
           IF ( time_since_reference_point >= time_turbine_on )  THEN
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    DO  k = nzb+1,  MAXVAL(k_hub) + MAXVAL(k_smear)
                       tend(k,j,i) = tend(k,j,i) - rot_tend_z(k,j,i) *                                 &
                                     MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
    
    
        CASE DEFAULT
           CONTINUE
    
        END SELECT
    
    
     END SUBROUTINE wtm_actions
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Application of the additional forces generated by the wind turbine on the flow components
    !> (tendency terms)
    !> Call for grid point i,j
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE wtm_actions_ij( i, j, location )
    
    
        CHARACTER (LEN=*) ::  location  !<
    
        INTEGER(iwp) ::  i  !< running index
        INTEGER(iwp) ::  j  !< running index
        INTEGER(iwp) ::  k  !< running index
    
        SELECT CASE ( location )
    
        CASE ( 'before_timestep' )
    
           CALL wtm_forces
    
        CASE ( 'u-tendency' )
    !
    !--    Apply the x-component of the force to the u-component of the flow:
           IF ( time_since_reference_point >= time_turbine_on )  THEN
              DO  k = nzb+1,  MAXVAL(k_hub) + MAXVAL(k_smear)
    !
    !--          Calculate the thrust generated by the nacelle and the tower:
                 tend_nac_x  = 0.5_wp * nac_cd_surf(k,j,i) * SIGN( u(k,j,i)**2 , u(k,j,i) )
                 tend_tow_x  = 0.5_wp * tow_cd_surf(k,j,i) * SIGN( u(k,j,i)**2 , u(k,j,i) )
                 tend(k,j,i) = tend(k,j,i) + ( - rot_tend_x(k,j,i) - tend_nac_x - tend_tow_x ) *       &
                               MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
              ENDDO
           ENDIF
    
        CASE ( 'v-tendency' )
    !
    !--    Apply the y-component of the force to the v-component of the flow:
           IF ( time_since_reference_point >= time_turbine_on )  THEN
              DO  k = nzb+1,  MAXVAL(k_hub) + MAXVAL(k_smear)
                 tend_nac_y  = 0.5_wp * nac_cd_surf(k,j,i) * SIGN( v(k,j,i)**2 , v(k,j,i) )
                 tend_tow_y  = 0.5_wp * tow_cd_surf(k,j,i) * SIGN( v(k,j,i)**2 , v(k,j,i) )
                 tend(k,j,i) = tend(k,j,i) + ( - rot_tend_y(k,j,i) - tend_nac_y - tend_tow_y )  *      &
                               MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                 ENDDO
           ENDIF
    
        CASE ( 'w-tendency' )
    !
    !--    Apply the z-component of the force to the w-component of the flow:
           IF ( time_since_reference_point >= time_turbine_on )  THEN
              DO  k = nzb+1,  MAXVAL(k_hub) + MAXVAL(k_smear)
                 tend(k,j,i) = tend(k,j,i) - rot_tend_z(k,j,i) *                                       &
                               MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
              ENDDO
           ENDIF
    
    
        CASE DEFAULT
           CONTINUE
    
        END SELECT
    
    
     END SUBROUTINE wtm_actions_ij
    
    
     SUBROUTINE wtm_data_output
    
    
        INTEGER(iwp) ::  return_value  !< returned status value of called function
        INTEGER(iwp) ::  t_ind = 0     !< time index
    
        IF ( myid == 0 )  THEN
    
    !
    !--    At the first call of this routine write the spatial coordinates:
           IF ( .NOT. initial_write_coordinates )  THEN
              ALLOCATE ( output_values_1d_target(1:n_turbines) )
              output_values_1d_target = hub_x(1:n_turbines)
              output_values_1d_pointer => output_values_1d_target
              return_value = dom_write_var( nc_filename,                                               &
                                         'x',                                                          &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1/),                                     &
                                         bounds_end       = (/n_turbines/) )
    
              output_values_1d_target = hub_y(1:n_turbines)
              output_values_1d_pointer => output_values_1d_target
              return_value = dom_write_var( nc_filename,                                               &
                                         'y',                                                          &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1/),                                     &
                                         bounds_end       = (/n_turbines/) )
    
              output_values_1d_target = hub_z(1:n_turbines)
              output_values_1d_pointer => output_values_1d_target
              return_value = dom_write_var( nc_filename,                                               &
                                         'z',                                                          &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1/),                                     &
                                         bounds_end       = (/n_turbines/) )
    
              output_values_1d_target = rotor_radius(1:n_turbines) * 2.0_wp
              output_values_1d_pointer => output_values_1d_target
              return_value = dom_write_var( nc_filename,                                               &
                                         'rotor_diameter',                                             &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1/),                                     &
                                         bounds_end       = (/n_turbines/) )
    
              output_values_1d_target = tower_diameter(1:n_turbines)
              output_values_1d_pointer => output_values_1d_target
              return_value = dom_write_var( nc_filename,                                               &
                                         'tower_diameter',                                             &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1/),                                     &
                                         bounds_end       = (/n_turbines/) )
    
              initial_write_coordinates = .TRUE.
              DEALLOCATE ( output_values_1d_target )
           ENDIF
    
           t_ind = t_ind + 1
    
           ALLOCATE ( output_values_1d_target(1:n_turbines) )
           output_values_1d_target = rotor_speed(:)
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'rotor_speed',                                                &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = generator_speed(:)
           output_values_1d_pointer => output_values_1d_target
           return_value = dom_write_var( nc_filename,                                                  &
                                         'generator_speed',                                            &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = torque_gen_old(:)
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'generator_torque',                                           &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = torque_total(:)
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'rotor_torque',                                               &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = pitch_angle(:)
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'pitch_angle',                                                &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = torque_gen_old(:) * generator_speed(:) * generator_efficiency
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'generator_power',                                            &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           DO  inot = 1, n_turbines
              output_values_1d_target(inot) = torque_total(inot) * rotor_speed(inot) * air_density
           ENDDO
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'rotor_power',                                                &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = thrust_rotor(:)
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'rotor_thrust',                                               &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = wdir(:)*180.0_wp/pi
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'wind_direction',                                             &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_1d_target = (yaw_angle(:)) * 180.0_wp / pi
           output_values_1d_pointer => output_values_1d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'yaw_angle',                                                  &
                                         values_realwp_1d = output_values_1d_pointer,                  &
                                         bounds_start     = (/1, t_ind/),                              &
                                         bounds_end       = (/n_turbines, t_ind /) )
    
           output_values_0d_target = time_since_reference_point
           output_values_0d_pointer => output_values_0d_target
    
           return_value = dom_write_var( nc_filename,                                                  &
                                         'time',                                                       &
                                         values_realwp_0d = output_values_0d_pointer,                  &
                                         bounds_start     = (/t_ind/),                                 &
                                         bounds_end       = (/t_ind/) )
    
           DEALLOCATE ( output_values_1d_target )
    
    !
    !--    Check if DOM reported any error
           dom_error_message = dom_get_error_message()
           IF ( TRIM( dom_error_message ) /= '' )  THEN
              message_string = 'error while writing output: "' // TRIM( dom_error_message ) // '"'
              CALL message( 'wtm_data_output', 'WTM0008', 0, 1, 0, 6, 0 )
           ENDIF
    
        ENDIF
    
     END SUBROUTINE wtm_data_output
    
     END MODULE wind_turbine_model_mod
    