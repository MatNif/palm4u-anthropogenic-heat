!> @file read_restart_data_mod.f90
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
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Reads restart data from restart-file(s) (binary format).
!>
!> @todo: Revise max_pr_cs (profiles for chemistry)
!> @todo: Modularize reading of restart data for diagnostic quantities, which
!>        is not possible with the current module-interface structure
!--------------------------------------------------------------------------------------------------!
 MODULE read_restart_data_mod

    USE arrays_3d,                                                                                 &
        ONLY:  mean_inflow_profiles,                                                               &
               pt_init,                                                                            &
               q_init,                                                                             &
               ref_state,                                                                          &
               sa_init,                                                                            &
               s_init,                                                                             &
               u_init,                                                                             &
               ug,                                                                                 &
               v_init,                                                                             &
               vg,                                                                                 &
               e,                                                                                  &
               kh,                                                                                 &
               km,                                                                                 &
               p,                                                                                  &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               s,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               vpt,                                                                                &
               w

    USE averaging

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE chem_modules,                                                                              &
        ONLY:  max_pr_cs

    USE control_parameters

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE gust_mod,                                                                                  &
        ONLY:  gust_module_enabled

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nx_on_file,                                                                         &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               ny_on_file,                                                                         &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE indoor_model_mod,                                                                          &
        ONLY:  time_indoor

    USE kinds

    USE model_1d_mod,                                                                              &
        ONLY:  damp_level_1d,                                                                      &
               dt_pr_1d,                                                                           &
               dt_run_control_1d,                                                                  &
               end_time_1d

    USE module_interface,                                                                          &
        ONLY:  module_interface_rrd_global,                                                        &
               module_interface_rrd_local,                                                         &
               module_interface_rrd_local_spinup

    USE netcdf_interface,                                                                          &
        ONLY:  netcdf_precision,                                                                   &
               output_for_t0

    USE particle_attributes,                                                                       &
        ONLY:  first_call_lpm,                                                                     &
               particle_advection,                                                                 &
               time_dopts

    USE pegrid

    USE radiation_model_mod,                                                                       &
        ONLY:  time_radiation

    USE random_function_mod,                                                                       &
        ONLY:  random_iv,                                                                          &
               random_iy

    USE random_generator_parallel,                                                                 &
        ONLY:  id_random_array,                                                                    &
               seq_random_array

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rd_mpi_io_close,                                                                    &
               rd_mpi_io_open,                                                                     &
               rrd_mpi_io,                                                                         &
               rrd_mpi_io_global_array

    USE spectra_mod,                                                                               &
        ONLY:  average_count_sp,                                                                   &
               spectrum_x,                                                                         &
               spectrum_y

    USE surface_data_output_mod,                                                                   &
        ONLY:  surface_data_output_rrd_global,                                                     &
               surface_data_output_rrd_local

    USE surface_mod,                                                                               &
        ONLY:  surface_rrd_local

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               hom_sum,                                                                            &
               pr_max,                                                                             &
               statistic_regions,                                                                  &
               u_max,                                                                              &
               u_max_ijk,                                                                          &
               v_max,                                                                              &
               v_max_ijk,                                                                          &
               w_max,                                                                              &
               w_max_ijk,                                                                          &
               z_i

    USE user,                                                                                      &
        ONLY:  user_module_enabled

    USE virtual_measurement_mod,                                                                   &
        ONLY:  time_virtual_measurement_pr,                                                        &
               time_virtual_measurement_ts,                                                        &
               time_virtual_measurement_tr


    IMPLICIT NONE


    INTERFACE rrd_global
       MODULE PROCEDURE rrd_global
    END INTERFACE rrd_global

    INTERFACE rrd_global_spinup
       MODULE PROCEDURE rrd_global_spinup
    END INTERFACE rrd_global_spinup

    INTERFACE rrd_read_parts_of_global
       MODULE PROCEDURE rrd_read_parts_of_global
    END INTERFACE rrd_read_parts_of_global

    INTERFACE rrd_local
       MODULE PROCEDURE rrd_local
    END INTERFACE rrd_local

    INTERFACE rrd_local_spinup
       MODULE PROCEDURE rrd_local_spinup
    END INTERFACE rrd_local_spinup

    INTERFACE rrd_skip_global
       MODULE PROCEDURE rrd_skip_global
    END INTERFACE rrd_skip_global


    PUBLIC rrd_global,                                                                             &
           rrd_global_spinup,                                                                      &
           rrd_read_parts_of_global,                                                               &
           rrd_local,                                                                              &
           rrd_local_spinup,                                                                       &
           rrd_skip_global


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads values of global control variables from restart-file (binary format)
!> created by PE0 of the previous run
!------------------------------------------------------------------------------!
 SUBROUTINE rrd_global


    CHARACTER(LEN=10) ::  binary_version_global  !<
    CHARACTER(LEN=10) ::  version_on_file        !<
    CHARACTER(LEN=20) ::  tmp_name               !< temporary variable

    INTEGER ::  i                                !< loop index

    LOGICAL ::  array_found                      !<
    LOGICAL ::  found                            !<


    CALL location_message( 'read global restart data', 'start' )

!
!-- Caution: When any of the read instructions have been changed, the
!-- -------  version number stored in the variable binary_version_global has
!--          to be increased. The same changes must also be done in wrd_write_global.
    binary_version_global = '22.04'

    IF ( TRIM( restart_data_format_input ) == 'fortran_binary' )  THEN
!
!--    Input in Fortran binary format
       CALL check_open( 13 )
!
!--    Make version number check first
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file

       IF ( TRIM( version_on_file ) /= TRIM( binary_version_global ) )  THEN
          WRITE( message_string, * ) 'version mismatch concerning ',           &
                                     'binary_version_global:',                 &
                                     '&version on file    = "',                &
                                     TRIM( version_on_file ), '"',             &
                                     '&version in program = "',                &
                                     TRIM( binary_version_global ), '"'
          CALL message( 'rrd_global', 'PAC0268', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Read number of PEs and horizontal index bounds of all PEs used in the
!--    previous run
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( TRIM( restart_string(1:length) ) /= 'numprocs' )  THEN
          WRITE( message_string, * ) 'numprocs not found in data from prior run on PE ', myid
          CALL message( 'rrd_global', 'PAC0269', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  numprocs_previous_run

       IF ( .NOT. ALLOCATED( hor_index_bounds_previous_run ) )  THEN
          ALLOCATE( hor_index_bounds_previous_run(4,0:numprocs_previous_run-1) )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'hor_index_bounds' )  THEN
          WRITE( message_string, * ) 'hor_index_bounds not found in data ',    &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_global', 'PAC0270', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  hor_index_bounds_previous_run

!
!--    Read vertical number of gridpoints and number of different areas used
!--    for computing statistics. Allocate arrays depending on these values,
!--    which are needed for the following read instructions.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'nz' )  THEN
          WRITE( message_string, * ) 'nz not found in data from prior run on PE ', myid
          CALL message( 'rrd_global', 'PAC0271', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  nz

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'max_pr_user' )  THEN
          WRITE( message_string, * ) 'max_pr_user not found in data from ',    &
                                     'prior run on PE ', myid
          CALL message( 'rrd_global', 'PAC0272', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  max_pr_user    ! This value is checked against the number of
                                   ! user profiles given for the current run
                                   ! in routine user_parin (it has to match)

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'statistic_regions' )  THEN
          WRITE( message_string, * ) 'statistic_regions not found in data ',   &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_global', 'PAC0273', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  statistic_regions

!
!--    The following global arrays (better to say, they have the same size and values on each
!--    subdomain) are by default allocated in routine parin, but not in case of restarts!
       IF ( .NOT. ALLOCATED( ug ) )  THEN
          ALLOCATE( ug(0:nz+1), u_init(0:nz+1), vg(0:nz+1),                                        &
                    v_init(0:nz+1), pt_init(0:nz+1), q_init(0:nz+1),                               &
                    ref_state(0:nz+1), s_init(0:nz+1), sa_init(0:nz+1),                            &
                    hom(0:nz+1,2,pr_max,0:statistic_regions),                                      &
                    hom_sum(0:nz+1,pr_max,0:statistic_regions) )
          hom = 0.0_wp
          hom_sum = 0.0_wp
       ENDIF

!
!--    Now read all control parameters:
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO WHILE ( restart_string(1:length) /= 'binary_version_local' )

          found = .FALSE.

          SELECT CASE ( restart_string(1:length) )

             CASE ( 'advected_distance_x' )
                READ ( 13 )  advected_distance_x
             CASE ( 'advected_distance_y' )
                READ ( 13 )  advected_distance_y
             CASE ( 'allow_negative_scalar_values' )
                READ ( 13 )  allow_negative_scalar_values
             CASE ( 'alpha_surface' )
                READ ( 13 )  alpha_surface
             CASE ( 'average_count_pr' )
                READ ( 13 )  average_count_pr
             CASE ( 'average_count_sp' )
                READ ( 13 )  average_count_sp
             CASE ( 'average_count_3d' )
                READ ( 13 )  average_count_3d
             CASE ( 'bc_e_b' )
                READ ( 13 )  bc_e_b
             CASE ( 'bc_lr' )
                READ ( 13 )  bc_lr
             CASE ( 'bc_ns' )
                READ ( 13 )  bc_ns
             CASE ( 'bc_p_b' )
                READ ( 13 )  bc_p_b
             CASE ( 'bc_p_t' )
                READ ( 13 )  bc_p_t
             CASE ( 'bc_pt_b' )
                READ ( 13 )  bc_pt_b
             CASE ( 'bc_pt_t' )
                READ ( 13 )  bc_pt_t
             CASE ( 'bc_pt_t_val' )
                READ ( 13 )  bc_pt_t_val
             CASE ( 'bc_q_b' )
                READ ( 13 )  bc_q_b
             CASE ( 'bc_q_t' )
                READ ( 13 )  bc_q_t
             CASE ( 'bc_q_t_val' )
                READ ( 13 )  bc_q_t_val
             CASE ( 'bc_s_b' )
                READ ( 13 )  bc_s_b
             CASE ( 'bc_s_t' )
                READ ( 13 )  bc_s_t
             CASE ( 'bc_uv_b' )
                READ ( 13 )  bc_uv_b
             CASE ( 'bc_uv_t' )
                READ ( 13 )  bc_uv_t
             CASE ( 'building_height' )
                READ ( 13 )  building_height
             CASE ( 'building_length_x' )
                READ ( 13 )  building_length_x
             CASE ( 'building_length_y' )
                READ ( 13 )  building_length_y
             CASE ( 'building_wall_left' )
                READ ( 13 )  building_wall_left
             CASE ( 'building_wall_south' )
                READ ( 13 )  building_wall_south
             CASE ( 'call_psolver_at_all_substeps' )
                READ ( 13 )  call_psolver_at_all_substeps
             CASE ( 'canyon_height' )
                READ ( 13 )  canyon_height
             CASE ( 'canyon_wall_left' )
                READ ( 13 )  canyon_wall_left
             CASE ( 'canyon_wall_south' )
                READ ( 13 )  canyon_wall_south
             CASE ( 'canyon_width_x' )
                READ ( 13 )  canyon_width_x
             CASE ( 'canyon_width_y' )
                READ ( 13 )  canyon_width_y
             CASE ( 'cfl_factor' )
                READ ( 13 )  cfl_factor
             CASE ( 'cloud_droplets' )
                READ ( 13 )  cloud_droplets
             CASE ( 'collective_wait' )
                READ ( 13 )  collective_wait
             CASE ( 'conserve_volume_flow' )
                READ ( 13 )  conserve_volume_flow
             CASE ( 'conserve_volume_flow_mode' )
                READ ( 13 )  conserve_volume_flow_mode
             CASE ( 'constant_flux_layer' )
                READ ( 13 )  constant_flux_layer
             CASE ( 'coupling_start_time' )
                READ ( 13 )  coupling_start_time
             CASE ( 'current_timestep_number' )
                READ ( 13 )  current_timestep_number
             CASE ( 'cycle_mg' )
                READ ( 13 )  cycle_mg
             CASE ( 'damp_level_1d' )
                READ ( 13 )  damp_level_1d
             CASE ( 'origin_date_time' )
                READ ( 13 )  origin_date_time
             CASE ( 'dissipation_1d' )
                READ ( 13 )  dissipation_1d
             CASE ( 'dp_external' )
                READ ( 13 )  dp_external
             CASE ( 'dp_level_b' )
                READ ( 13 )  dp_level_b
             CASE ( 'dp_smooth' )
                READ ( 13 )  dp_smooth
             CASE ( 'dpdxy' )
                READ ( 13 )  dpdxy
             CASE ( 'dt_3d' )
                READ ( 13 )  dt_3d
             CASE ( 'dt_pr_1d' )
                READ ( 13 )  dt_pr_1d
             CASE ( 'dt_run_control_1d' )
                READ ( 13 )  dt_run_control_1d
             CASE ( 'dx' )
                READ ( 13 )  dx
             CASE ( 'dy' )
                READ ( 13 )  dy
             CASE ( 'dz' )
                READ ( 13 )  dz
             CASE ( 'dz_max' )
                READ ( 13 )  dz_max
             CASE ( 'dz_stretch_factor' )
                READ ( 13 )  dz_stretch_factor
             CASE ( 'dz_stretch_factor_array' )
                READ ( 13 )  dz_stretch_factor_array
             CASE ( 'dz_stretch_level' )
                READ ( 13 )  dz_stretch_level
             CASE ( 'dz_stretch_level_end' )
                READ ( 13 )  dz_stretch_level_end
             CASE ( 'dz_stretch_level_start' )
                READ ( 13 )  dz_stretch_level_start
             CASE ( 'e_min' )
                READ ( 13 )  e_min
             CASE ( 'end_time_1d' )
                READ ( 13 )  end_time_1d
             CASE ( 'fft_method' )
                READ ( 13 )  fft_method
             CASE ( 'first_call_lpm' )
                READ ( 13 )  first_call_lpm
             CASE ( 'galilei_transformation' )
                READ ( 13 )  galilei_transformation
             CASE ( 'hom' )
                READ ( 13 )  hom
             CASE ( 'hom_sum' )
                READ ( 13 )  hom_sum
             CASE ( 'homogenize_surface_temperature' )
                READ ( 13 )  homogenize_surface_temperature
             CASE ( 'humidity' )
                READ ( 13 )  humidity
             CASE ( 'inflow_disturbance_begin' )
                READ ( 13 )  inflow_disturbance_begin
             CASE ( 'inflow_disturbance_end' )
                READ ( 13 )  inflow_disturbance_end
             CASE ( 'km_constant' )
                READ ( 13 )  km_constant
             CASE ( 'large_scale_forcing' )
                READ ( 13 )  large_scale_forcing
             CASE ( 'large_scale_subsidence' )
                READ ( 13 )  large_scale_subsidence
             CASE ( 'latitude' )
                READ ( 13 )  latitude
             CASE ( 'longitude' )
                READ ( 13 )  longitude
             CASE ( 'loop_optimization' )
                READ ( 13 )  loop_optimization
             CASE ( 'masking_method' )
                READ ( 13 )  masking_method
             CASE ( 'mean_inflow_profiles' )
                IF ( .NOT. ALLOCATED( mean_inflow_profiles ) )  THEN
                   ALLOCATE( mean_inflow_profiles(0:nz+1,1:num_mean_inflow_profiles) )
                ENDIF
                READ ( 13 )  mean_inflow_profiles
             CASE ( 'mg_cycles' )
                READ ( 13 )  mg_cycles
             CASE ( 'mg_switch_to_pe0_level' )
                READ ( 13 )  mg_switch_to_pe0_level
             CASE ( 'mixing_length_1d' )
                READ ( 13 )  mixing_length_1d
             CASE ( 'momentum_advec' )
                READ ( 13 )  momentum_advec
             CASE ( 'netcdf_precision' )
                READ ( 13 )  netcdf_precision
             CASE ( 'neutral' )
                READ ( 13 )  neutral
             CASE ( 'ngsrb' )
                READ ( 13 )  ngsrb
             CASE ( 'nsor' )
                READ ( 13 )  nsor
             CASE ( 'nsor_ini' )
                READ ( 13 )  nsor_ini
             CASE ( 'nudging' )
                READ ( 13 )  nudging
             CASE ( 'num_leg' )
                READ ( 13 )  num_leg
             CASE ( 'nx' )
                READ ( 13 )  nx
                nx_on_file = nx
             CASE ( 'ny' )
                READ ( 13 )  ny
                ny_on_file = ny
             CASE ( 'ocean_mode' )
                READ ( 13 )  ocean_mode
             CASE ( 'omega' )
                READ ( 13 )  omega
             CASE ( 'omega_sor' )
                READ ( 13 )  omega_sor
             CASE ( 'outflow_damping_factor' )
                READ ( 13 )  outflow_damping_factor
             CASE ( 'outflow_damping_width' )
                READ ( 13 )  outflow_damping_width
             CASE ( 'output_for_t0' )
                READ (13)    output_for_t0
             CASE ( 'passive_scalar' )
                READ ( 13 )  passive_scalar
             CASE ( 'prandtl_number' )
                READ ( 13 )  prandtl_number
             CASE ( 'psolver' )
                READ ( 13 )  psolver
             CASE ( 'pt_damping_factor' )
                READ ( 13 )  pt_damping_factor
             CASE ( 'pt_damping_width' )
                READ ( 13 )  pt_damping_width
             CASE ( 'pt_init' )
                READ ( 13 )  pt_init
             CASE ( 'pt_reference' )
                READ ( 13 )  pt_reference
             CASE ( 'pt_surface' )
                READ ( 13 )  pt_surface
             CASE ( 'pt_surface_heating_rate' )
                READ ( 13 )  pt_surface_heating_rate
             CASE ( 'pt_vertical_gradient' )
                READ ( 13 )  pt_vertical_gradient
             CASE ( 'pt_vertical_gradient_level' )
                READ ( 13 )  pt_vertical_gradient_level
             CASE ( 'pt_vertical_gradient_level_ind' )
                READ ( 13 )  pt_vertical_gradient_level_ind
             CASE ( 'q_init' )
                READ ( 13 )  q_init
             CASE ( 'q_surface' )
                READ ( 13 )  q_surface
             CASE ( 'q_vertical_gradient' )
                READ ( 13 )  q_vertical_gradient
             CASE ( 'q_vertical_gradient_level' )
                READ ( 13 )  q_vertical_gradient_level
             CASE ( 'q_vertical_gradient_level_ind' )
                READ ( 13 )  q_vertical_gradient_level_ind
             CASE ( 'random_generator' )
                READ ( 13 )  random_generator
             CASE ( 'random_heatflux' )
                READ ( 13 )  random_heatflux
             CASE ( 'rans_mode' )
                READ ( 13 )  rans_mode
             CASE ( 'rayleigh_damping_factor' )
                READ ( 13 )  rayleigh_damping_factor
             CASE ( 'rayleigh_damping_height' )
                READ ( 13 )  rayleigh_damping_height
             CASE ( 'ref_state' )
                READ ( 13 )  ref_state
             CASE ( 'reference_state' )
                READ ( 13 )  reference_state
             CASE ( 'residual_limit' )
                READ ( 13 )  residual_limit
             CASE ( 'roughness_length' )
                READ ( 13 )  roughness_length
             CASE ( 'runnr' )
                READ ( 13 )  runnr
             CASE ( 's_init' )
                READ ( 13 )  s_init
             CASE ( 's_surface' )
                READ ( 13 )  s_surface
             CASE ( 's_vertical_gradient' )
                READ ( 13 )  s_vertical_gradient
             CASE ( 's_vertical_gradient_level' )
                READ ( 13 )  s_vertical_gradient_level
             CASE ( 's_vertical_gradient_level_ind' )
                READ ( 13 )  s_vertical_gradient_level_ind
             CASE ( 'scalar_advec' )
                READ ( 13 )  scalar_advec
             CASE ( 'simulated_time' )
                READ ( 13 )  simulated_time
             CASE ( 'spectrum_x' )
                IF ( .NOT. ALLOCATED( spectrum_x ) )  THEN
                   ALLOCATE( spectrum_x( 1:nx/2, 1:100, 1:10 ) )
                ENDIF
                READ ( 13 )  spectrum_x
             CASE ( 'spectrum_y' )
                IF ( .NOT. ALLOCATED( spectrum_y ) )  THEN
                   ALLOCATE( spectrum_y( 1:ny/2, 1:100, 1:10 ) )
                ENDIF
                READ ( 13 )  spectrum_y
             CASE ( 'spinup_time' )
                READ ( 13 )  spinup_time
             CASE ( 'subs_vertical_gradient' )
                READ ( 13 )  subs_vertical_gradient
             CASE ( 'subs_vertical_gradient_level' )
                READ ( 13 )  subs_vertical_gradient_level
             CASE ( 'subs_vertical_gradient_level_i' )
                READ ( 13 )  subs_vertical_gradient_level_i
             CASE ( 'surface_heatflux' )
                READ ( 13 )  surface_heatflux
             CASE ( 'surface_pressure' )
                READ ( 13 )  surface_pressure
             CASE ( 'surface_scalarflux' )
                READ ( 13 )  surface_scalarflux
             CASE ( 'surface_waterflux' )
                READ ( 13 )  surface_waterflux
             CASE ( 'time_coupling' )
                READ ( 13 )  time_coupling
             CASE ( 'time_disturb' )
                READ ( 13 )  time_disturb
             CASE ( 'time_do2d_xy' )
                READ ( 13 )  time_do2d_xy
             CASE ( 'time_do2d_xz' )
                READ ( 13 )  time_do2d_xz
             CASE ( 'time_do2d_yz' )
                READ ( 13 )  time_do2d_yz
             CASE ( 'time_do3d' )
                READ ( 13 )  time_do3d
             CASE ( 'time_do_av' )
                READ ( 13 )  time_do_av
             CASE ( 'time_do_sla' )
                READ ( 13 )  time_do_sla
             CASE ( 'time_domask' )
                READ ( 13 )  time_domask
             CASE ( 'time_dopr' )
                READ ( 13 )  time_dopr
             CASE ( 'time_dopr_av' )
                READ ( 13 )  time_dopr_av
             CASE ( 'time_dopr_listing' )
                READ ( 13 )  time_dopr_listing
             CASE ( 'time_dopts' )
                READ ( 13 )  time_dopts
             CASE ( 'time_dosp' )
                READ ( 13 )  time_dosp
             CASE ( 'time_dots' )
                READ ( 13 )  time_dots
             CASE ( 'time_indoor' )
                READ ( 13 )  time_indoor
             CASE ( 'time_radiation' )
                READ ( 13 )  time_radiation
             CASE ( 'time_restart' )
                READ ( 13 )  time_restart
             CASE ( 'time_run_control' )
                READ ( 13 )  time_run_control
             CASE ( 'time_since_reference_point' )
                READ ( 13 )  time_since_reference_point
             CASE ( 'time_virtual_measurement_pr' )
                READ ( 13 )  time_virtual_measurement_pr
             CASE ( 'time_virtual_measurement_ts' )
                READ ( 13 )  time_virtual_measurement_ts
             CASE ( 'time_virtual_measurement_tr' )
                READ ( 13 )  time_virtual_measurement_tr
             CASE ( 'timestep_scheme' )
                READ ( 13 )  timestep_scheme
             CASE ( 'top_heatflux' )
                READ ( 13 )  top_heatflux
             CASE ( 'top_momentumflux_u' )
                READ ( 13 )  top_momentumflux_u
             CASE ( 'top_momentumflux_v' )
                READ ( 13 )  top_momentumflux_v
             CASE ( 'top_scalarflux' )
                READ ( 13 )  top_scalarflux
             CASE ( 'topography' )
                READ ( 13 )  topography
             CASE ( 'topography_grid_convention' )
                READ ( 13 )  topography_grid_convention
             CASE ( 'tsc' )
                READ ( 13 )  tsc
             CASE ( 'tunnel_height' )
                READ ( 13 )  tunnel_height
             CASE ( 'tunnel_length' )
                READ ( 13 )  tunnel_length
             CASE ( 'tunnel_wall_depth' )
                READ ( 13 )  tunnel_wall_depth
             CASE ( 'tunnel_width_x' )
                READ ( 13 )  tunnel_width_x
             CASE ( 'tunnel_width_y' )
                READ ( 13 )  tunnel_width_y
             CASE ( 'turbulence_closure' )
                READ ( 13 )  turbulence_closure
             CASE ( 'u_bulk' )
                READ ( 13 )  u_bulk
             CASE ( 'u_init' )
                READ ( 13 )  u_init
             CASE ( 'u_max' )
                READ ( 13 )  u_max
             CASE ( 'u_max_ijk' )
                READ ( 13 )  u_max_ijk
             CASE ( 'ug' )
                READ ( 13 )  ug
             CASE ( 'ug_surface' )
                READ ( 13 )  ug_surface
             CASE ( 'ug_vertical_gradient' )
                READ ( 13 )  ug_vertical_gradient
             CASE ( 'ug_vertical_gradient_level' )
                READ ( 13 )  ug_vertical_gradient_level
             CASE ( 'ug_vertical_gradient_level_ind' )
                READ ( 13 )  ug_vertical_gradient_level_ind
             CASE ( 'use_surface_fluxes' )
                READ ( 13 )  use_surface_fluxes
             CASE ( 'use_top_fluxes' )
                READ ( 13 )  use_top_fluxes
             CASE ( 'use_ug_for_galilei_tr' )
                READ ( 13 )  use_ug_for_galilei_tr
             CASE ( 'use_upstream_for_tke' )
                READ ( 13 )  use_upstream_for_tke
             CASE ( 'v_bulk' )
                READ ( 13 )  v_bulk
             CASE ( 'v_init' )
                READ ( 13 )  v_init
             CASE ( 'v_max' )
                READ ( 13 )  v_max
             CASE ( 'v_max_ijk' )
                READ ( 13 )  v_max_ijk
             CASE ( 'vg' )
                READ ( 13 )  vg
             CASE ( 'vg_surface' )
                READ ( 13 )  vg_surface
             CASE ( 'vg_vertical_gradient' )
                READ ( 13 )  vg_vertical_gradient
             CASE ( 'vg_vertical_gradient_level' )
                READ ( 13 )  vg_vertical_gradient_level
             CASE ( 'vg_vertical_gradient_level_ind' )
                READ ( 13 )  vg_vertical_gradient_level_ind
             CASE ( 'virtual_flight' )
                READ ( 13 )  virtual_flight
             CASE ( 'volume_flow_area' )
                READ ( 13 )  volume_flow_area
             CASE ( 'volume_flow_initial' )
                READ ( 13 )  volume_flow_initial
             CASE ( 'w_max' )
                READ ( 13 )  w_max
             CASE ( 'w_max_ijk' )
                READ ( 13 )  w_max_ijk
             CASE ( 'wall_adjustment' )
                READ ( 13 )  wall_adjustment
             CASE ( 'wall_heatflux' )
                READ ( 13 )  wall_heatflux
             CASE ( 'wall_humidityflux' )
                READ ( 13 )  wall_humidityflux
             CASE ( 'wall_scalarflux' )
                READ ( 13 )  wall_scalarflux
             CASE ( 'y_shift' )
                READ ( 13 )  y_shift
             CASE ( 'z0h_factor' )
                READ ( 13 )  z0h_factor
             CASE ( 'zeta_max' )
                READ ( 13 )  zeta_max
             CASE ( 'zeta_min' )
                READ ( 13 )  zeta_min
             CASE ( 'z_i' )
                READ ( 13 )  z_i

             CASE DEFAULT
!
!--             Read global variables from of other modules
                CALL module_interface_rrd_global( found )
!
!--             Read global variables from surface_data_output_mod
                IF ( .NOT. found )  CALL surface_data_output_rrd_global( found )

                IF ( .NOT. found )  THEN
                   WRITE( message_string, * ) 'unknown variable named "',      &
                                              restart_string(1:length),           &
                                              '" found in global data from prior run on PE ', myid
                CALL message( 'rrd_global', 'PAC0274', 1, 2, 0, 6, 0 )

                ENDIF

          END SELECT
!
!--       Read next string
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO  ! End of loop for reading the restart string

       CALL close_file( 13 )

    ELSEIF ( restart_data_format_input(1:3) == 'mpi' )  THEN
!
!--    Read global restart data using MPI-IO
!--    ATTENTION: Arrays need to be read with routine rrd_mpi_io_global_array!

!
!--    Open the MPI-IO restart file.
       CALL rd_mpi_io_open( 'read', 'BININ' // TRIM( coupling_char ),                              &
                            open_for_global_io_only = .TRUE. )

!
!--    Make version number check first
       CALL rrd_mpi_io( 'binary_version_global',  version_on_file )

       IF ( TRIM( version_on_file ) /= TRIM( binary_version_global ) )  THEN
          WRITE( message_string, * ) 'version mismatch concerning binary_version_global:',         &
                                     '&version on file    = "', TRIM( version_on_file ), '"',      &
                                     '&version in program = "', TRIM( binary_version_global ), '"'
          CALL message( 'rrd_global', 'PAC0268', 1, 2, 0, 6, 0 )
       ENDIF

       CALL rrd_mpi_io( 'numprocs',  numprocs_previous_run )
       CALL rrd_mpi_io( 'nz' , nz )
       CALL rrd_mpi_io( 'max_pr_user',  max_pr_user )
       CALL rrd_mpi_io( 'statistic_regions', statistic_regions )

!
!--    The following global arrays (better to say, they have the same size and values on each
!--    subdomain) are by default allocated in routine parin, but not in case of restarts!
       IF ( .NOT. ALLOCATED( ug ) )  THEN
          ALLOCATE( ug(0:nz+1), u_init(0:nz+1), vg(0:nz+1),                                        &
                    v_init(0:nz+1), pt_init(0:nz+1), q_init(0:nz+1),                               &
                    ref_state(0:nz+1), s_init(0:nz+1), sa_init(0:nz+1),                            &
                    hom(0:nz+1,2,pr_max,0:statistic_regions),                                      &
                    hom_sum(0:nz+1,pr_max,0:statistic_regions) )
          hom = 0.0_wp
          hom_sum = 0.0_wp
       ENDIF

       CALL rrd_mpi_io( 'advected_distance_x',  advected_distance_x )
       CALL rrd_mpi_io( 'advected_distance_y', advected_distance_y )
       CALL rrd_mpi_io( 'allow_negative_scalar_values', allow_negative_scalar_values )
       CALL rrd_mpi_io( 'alpha_surface', alpha_surface )
       CALL rrd_mpi_io( 'average_count_pr', average_count_pr )
       CALL rrd_mpi_io( 'average_count_sp', average_count_sp )
       CALL rrd_mpi_io( 'average_count_3d', average_count_3d )
       CALL rrd_mpi_io( 'bc_e_b', bc_e_b )
       CALL rrd_mpi_io( 'bc_lr', bc_lr )
       CALL rrd_mpi_io( 'bc_ns', bc_ns )
       CALL rrd_mpi_io( 'bc_p_b', bc_p_b )
       CALL rrd_mpi_io( 'bc_p_t', bc_p_t )
       CALL rrd_mpi_io( 'bc_pt_b', bc_pt_b )
       CALL rrd_mpi_io( 'bc_pt_t', bc_pt_t )
       CALL rrd_mpi_io( 'bc_pt_t_val', bc_pt_t_val )
       CALL rrd_mpi_io( 'bc_q_b', bc_q_b )
       CALL rrd_mpi_io( 'bc_q_t', bc_q_t )
       CALL rrd_mpi_io( 'bc_q_t_val', bc_q_t_val )
       CALL rrd_mpi_io( 'bc_s_b', bc_s_b )
       CALL rrd_mpi_io( 'bc_s_t', bc_s_t )
       CALL rrd_mpi_io( 'bc_uv_b', bc_uv_b )
       CALL rrd_mpi_io( 'bc_uv_t', bc_uv_t )
       CALL rrd_mpi_io( 'biometeorology', biometeorology )
       CALL rrd_mpi_io( 'building_height', building_height )
       CALL rrd_mpi_io( 'building_length_x', building_length_x )
       CALL rrd_mpi_io( 'building_length_y', building_length_y )
       CALL rrd_mpi_io( 'building_wall_left', building_wall_left )
       CALL rrd_mpi_io( 'building_wall_south', building_wall_south )
       CALL rrd_mpi_io( 'bulk_cloud_model', bulk_cloud_model )
       CALL rrd_mpi_io( 'call_psolver_at_all_substeps', call_psolver_at_all_substeps )
       CALL rrd_mpi_io( 'canyon_height', canyon_height )
       CALL rrd_mpi_io( 'canyon_wall_left', canyon_wall_left )
       CALL rrd_mpi_io( 'canyon_wall_south', canyon_wall_south )
       CALL rrd_mpi_io( 'canyon_width_x',  canyon_width_x )
       CALL rrd_mpi_io( 'canyon_width_y', canyon_width_y )
       CALL rrd_mpi_io( 'cfl_factor', cfl_factor )
       CALL rrd_mpi_io( 'cloud_droplets',  cloud_droplets )
       CALL rrd_mpi_io( 'collective_wait', collective_wait )
       CALL rrd_mpi_io( 'conserve_volume_flow', conserve_volume_flow )
       CALL rrd_mpi_io( 'conserve_volume_flow_mode', conserve_volume_flow_mode )
       CALL rrd_mpi_io( 'constant_flux_layer', constant_flux_layer )
       CALL rrd_mpi_io( 'coupling_start_time', coupling_start_time )
       CALL rrd_mpi_io( 'current_timestep_number', current_timestep_number )
       CALL rrd_mpi_io( 'cycle_mg', cycle_mg )
       CALL rrd_mpi_io( 'damp_level_1d', damp_level_1d )
       CALL rrd_mpi_io( 'dissipation_1d', dissipation_1d )
       CALL rrd_mpi_io( 'dp_external', dp_external )
       CALL rrd_mpi_io( 'dp_level_b', dp_level_b )
       CALL rrd_mpi_io( 'dp_smooth', dp_smooth )
       CALL rrd_mpi_io_global_array( 'dpdxy', dpdxy )
       CALL rrd_mpi_io( 'dt_3d', dt_3d )
       CALL rrd_mpi_io( 'dt_pr_1d', dt_pr_1d )
       CALL rrd_mpi_io( 'dt_run_control_1d', dt_run_control_1d )
       CALL rrd_mpi_io( 'dx', dx )
       CALL rrd_mpi_io( 'dy', dy )
       CALL rrd_mpi_io_global_array( 'dz', dz )
       CALL rrd_mpi_io( 'dz_max', dz_max )
       CALL rrd_mpi_io( 'dz_stretch_factor', dz_stretch_factor )
       CALL rrd_mpi_io_global_array( 'dz_stretch_factor_array', dz_stretch_factor_array )
       CALL rrd_mpi_io( 'dz_stretch_level', dz_stretch_level )
       CALL rrd_mpi_io_global_array( 'dz_stretch_level_end', dz_stretch_level_end )
       CALL rrd_mpi_io_global_array( 'dz_stretch_level_start', dz_stretch_level_start )
       CALL rrd_mpi_io( 'e_min', e_min )
       CALL rrd_mpi_io( 'end_time_1d', end_time_1d )
       CALL rrd_mpi_io( 'fft_method', fft_method )
       CALL rrd_mpi_io( 'first_call_lpm', first_call_lpm )
       CALL rrd_mpi_io( 'galilei_transformation', galilei_transformation )
       CALL rrd_mpi_io( 'gust_module_enabled', gust_module_enabled )
       CALL rrd_mpi_io_global_array( 'hom', hom )
       CALL rrd_mpi_io_global_array( 'hom_sum', hom_sum )
       CALL rrd_mpi_io( 'homogenize_surface_temperature', homogenize_surface_temperature )
       CALL rrd_mpi_io( 'humidity', humidity )
       CALL rrd_mpi_io( 'inflow_disturbance_begin', inflow_disturbance_begin )
       CALL rrd_mpi_io( 'inflow_disturbance_end', inflow_disturbance_end )
       CALL rrd_mpi_io( 'km_constant', km_constant )
       CALL rrd_mpi_io( 'large_scale_forcing', large_scale_forcing )
       CALL rrd_mpi_io( 'large_scale_subsidence', large_scale_subsidence )
       CALL rrd_mpi_io( 'latitude', latitude )
       CALL rrd_mpi_io( 'longitude', longitude )
       CALL rrd_mpi_io( 'loop_optimization', loop_optimization )
       CALL rrd_mpi_io( 'masking_method', masking_method )
       CALL rd_mpi_io_check_array( 'mean_inflow_profiles', found = array_found )
       IF ( array_found)  THEN
          IF ( .NOT. ALLOCATED( mean_inflow_profiles ) )  THEN
             ALLOCATE( mean_inflow_profiles(0:nz+1,7) )
          ENDIF
          CALL rrd_mpi_io_global_array( 'mean_inflow_profiles', mean_inflow_profiles )
       ENDIF
       CALL rrd_mpi_io( 'mg_cycles', mg_cycles )
       CALL rrd_mpi_io( 'mg_switch_to_pe0_level', mg_switch_to_pe0_level )
       CALL rrd_mpi_io( 'mixing_length_1d', mixing_length_1d )
       CALL rrd_mpi_io( 'momentum_advec', momentum_advec )
!
!--    There is no module procedure for CHARACTER arrays
       DO  i = 1, SIZE( netcdf_precision , 1 )
          WRITE( tmp_name, '(A,I2.2)' )  'netcdf_precision_', i
          CALL rrd_mpi_io( TRIM( tmp_name ), netcdf_precision(i) )
       ENDDO
       CALL rrd_mpi_io( 'neutral', neutral )
       CALL rrd_mpi_io( 'ngsrb', ngsrb )
       CALL rrd_mpi_io( 'nsor', nsor )
       CALL rrd_mpi_io( 'nsor_ini', nsor_ini )
       CALL rrd_mpi_io( 'nudging', nudging )
       CALL rrd_mpi_io( 'num_leg', num_leg )
       CALL rrd_mpi_io( 'nx', nx )
       nx_on_file = nx
       CALL rrd_mpi_io( 'ny', ny )
       ny_on_file = ny
       CALL rrd_mpi_io( 'ocean_mode', ocean_mode )
       CALL rrd_mpi_io( 'omega', omega )
       CALL rrd_mpi_io( 'omega_sor', omega_sor )
       CALL rrd_mpi_io( 'origin_date_time', origin_date_time )       
       CALL rrd_mpi_io( 'outflow_damping_factor', outflow_damping_factor )
       CALL rrd_mpi_io( 'outflow_damping_width', outflow_damping_width )
       CALL rrd_mpi_io( 'output_for_t0', output_for_t0 )
       CALL rrd_mpi_io( 'particle_advection', particle_advection )
       CALL rrd_mpi_io( 'passive_scalar', passive_scalar )
       CALL rrd_mpi_io( 'prandtl_number', prandtl_number )
       CALL rrd_mpi_io( 'psolver', psolver )
       CALL rrd_mpi_io( 'pt_damping_factor', pt_damping_factor )
       CALL rrd_mpi_io( 'pt_damping_width', pt_damping_width )
       CALL rrd_mpi_io_global_array( 'pt_init', pt_init )
       CALL rrd_mpi_io( 'pt_reference', pt_reference )
       CALL rrd_mpi_io( 'pt_surface', pt_surface )
       CALL rrd_mpi_io( 'pt_surface_heating_rate', pt_surface_heating_rate )
       CALL rrd_mpi_io_global_array( 'pt_vertical_gradient', pt_vertical_gradient )
       CALL rrd_mpi_io_global_array( 'pt_vertical_gradient_level', pt_vertical_gradient_level )
       CALL rrd_mpi_io_global_array( 'pt_vertical_gradient_level_ind', pt_vertical_gradient_level_ind )
       CALL rrd_mpi_io_global_array( 'q_init', q_init )
       CALL rrd_mpi_io( 'q_surface', q_surface )
       CALL rrd_mpi_io_global_array( 'q_vertical_gradient', q_vertical_gradient )
       CALL rrd_mpi_io_global_array( 'q_vertical_gradient_level', q_vertical_gradient_level )
       CALL rrd_mpi_io_global_array( 'q_vertical_gradient_level_ind', q_vertical_gradient_level_ind )
       CALL rrd_mpi_io( 'random_generator', random_generator )
       CALL rrd_mpi_io( 'random_heatflux', random_heatflux )
       CALL rrd_mpi_io( 'rans_mode', rans_mode )
       CALL rrd_mpi_io( 'rayleigh_damping_factor', rayleigh_damping_factor )
       CALL rrd_mpi_io( 'rayleigh_damping_height', rayleigh_damping_height )
       CALL rrd_mpi_io_global_array( 'ref_state', ref_state )
       CALL rrd_mpi_io( 'reference_state', reference_state )
       CALL rrd_mpi_io( 'residual_limit', residual_limit )
       CALL rrd_mpi_io( 'roughness_length', roughness_length )
       CALL rrd_mpi_io( 'runnr', runnr )
       CALL rrd_mpi_io_global_array( 's_init', s_init )
       CALL rrd_mpi_io( 's_surface', s_surface )
       CALL rrd_mpi_io_global_array( 's_vertical_gradient', s_vertical_gradient )
       CALL rrd_mpi_io_global_array( 's_vertical_gradient_level', s_vertical_gradient_level )
       CALL rrd_mpi_io_global_array( 's_vertical_gradient_level_ind', s_vertical_gradient_level_ind )
       CALL rrd_mpi_io( 'scalar_advec', scalar_advec )
       CALL rrd_mpi_io( 'simulated_time', simulated_time )
       CALL rd_mpi_io_check_array( 'spectrum_x', found = array_found )
       IF (array_found )  THEN
           IF ( .NOT. ALLOCATED( spectrum_x ) )  THEN
              ALLOCATE( spectrum_x( 1:nx/2, 1:100, 1:10 ) )
           ENDIF
           CALL rrd_mpi_io_global_array( 'spectrum_x', spectrum_x )
       ENDIF
       CALL rd_mpi_io_check_array( 'spectrum_y', found = array_found )
       IF ( array_found )  THEN
           IF ( .NOT. ALLOCATED( spectrum_y ) )  THEN
              ALLOCATE( spectrum_y( 1:ny/2, 1:100, 1:10 ) )
           ENDIF
           CALL rrd_mpi_io_global_array( 'spectrum_y', spectrum_y )
       ENDIF
       CALL rrd_mpi_io( 'spinup_time ', spinup_time )
       CALL rrd_mpi_io_global_array( 'subs_vertical_gradient', subs_vertical_gradient )
       CALL rrd_mpi_io_global_array( 'subs_vertical_gradient_level', subs_vertical_gradient_level )
       CALL rrd_mpi_io_global_array( 'subs_vertical_gradient_level_i', subs_vertical_gradient_level_i )
       CALL rrd_mpi_io( 'surface_heatflux', surface_heatflux )
       CALL rrd_mpi_io( 'surface_pressure', surface_pressure )
       CALL rrd_mpi_io( 'surface_output', surface_output )
       CALL rrd_mpi_io( 'surface_scalarflux', surface_scalarflux )
       CALL rrd_mpi_io( 'surface_waterflux', surface_waterflux )
       CALL rrd_mpi_io( 'time_coupling', time_coupling )
       CALL rrd_mpi_io( 'time_disturb', time_disturb )
       CALL rrd_mpi_io( 'time_do2d_xy', time_do2d_xy )
       CALL rrd_mpi_io( 'time_do2d_xz', time_do2d_xz )
       CALL rrd_mpi_io( 'time_do2d_yz', time_do2d_yz )
       CALL rrd_mpi_io( 'time_do3d', time_do3d )
       CALL rrd_mpi_io( 'time_do_av', time_do_av )
       CALL rrd_mpi_io( 'time_do_sla', time_do_sla )
       CALL rrd_mpi_io_global_array( 'time_domask', time_domask )
       CALL rrd_mpi_io( 'time_dopr', time_dopr )
       CALL rrd_mpi_io( 'time_dopr_av', time_dopr_av )
       CALL rrd_mpi_io( 'time_dopr_listing', time_dopr_listing )
       CALL rrd_mpi_io( 'time_dopts', time_dopts )
       CALL rrd_mpi_io( 'time_dosp', time_dosp )
       CALL rrd_mpi_io( 'time_dots', time_dots )
       CALL rrd_mpi_io( 'time_indoor', time_indoor )
       CALL rrd_mpi_io( 'time_radiation', time_radiation )
       CALL rrd_mpi_io( 'time_restart', time_restart )
       CALL rrd_mpi_io( 'time_run_control', time_run_control )
       CALL rrd_mpi_io( 'time_since_reference_point', time_since_reference_point )
       CALL rrd_mpi_io( 'time_virtual_measurement_pr', time_virtual_measurement_pr )
       CALL rrd_mpi_io( 'time_virtual_measurement_ts', time_virtual_measurement_ts )
       CALL rrd_mpi_io( 'time_virtual_measurement_tr', time_virtual_measurement_tr )
       CALL rrd_mpi_io( 'timestep_scheme', timestep_scheme )
       CALL rrd_mpi_io( 'top_heatflux', top_heatflux )
       CALL rrd_mpi_io( 'top_momentumflux_u', top_momentumflux_u )
       CALL rrd_mpi_io( 'top_momentumflux_v', top_momentumflux_v )
       CALL rrd_mpi_io( 'top_scalarflux', top_scalarflux )
       CALL rrd_mpi_io( 'topography', topography )
       CALL rrd_mpi_io( 'topography_grid_convention', topography_grid_convention )
       CALL rrd_mpi_io_global_array( 'tsc', tsc )
       CALL rrd_mpi_io( 'tunnel_height', tunnel_height )
       CALL rrd_mpi_io( 'tunnel_length', tunnel_length )
       CALL rrd_mpi_io( 'tunnel_wall_depth', tunnel_wall_depth )
       CALL rrd_mpi_io( 'tunnel_width_x', tunnel_width_x )
       CALL rrd_mpi_io( 'tunnel_width_y', tunnel_width_y )
       CALL rrd_mpi_io( 'turbulence_closure', turbulence_closure )
       CALL rrd_mpi_io( 'u_bulk', u_bulk )
       CALL rrd_mpi_io_global_array( 'u_init', u_init )
       CALL rrd_mpi_io( 'u_max', u_max )
       CALL rrd_mpi_io_global_array( 'u_max_ijk', u_max_ijk )
       CALL rrd_mpi_io_global_array( 'ug', ug )
       CALL rrd_mpi_io( 'ug_surface', ug_surface )
       CALL rrd_mpi_io_global_array( 'ug_vertical_gradient', ug_vertical_gradient )
       CALL rrd_mpi_io_global_array( 'ug_vertical_gradient_level', ug_vertical_gradient_level )
       CALL rrd_mpi_io_global_array( 'ug_vertical_gradient_level_ind', ug_vertical_gradient_level_ind )
       CALL rrd_mpi_io( 'use_surface_fluxes', use_surface_fluxes )
       CALL rrd_mpi_io( 'use_top_fluxes', use_top_fluxes )
       CALL rrd_mpi_io( 'use_ug_for_galilei_tr', use_ug_for_galilei_tr )
       CALL rrd_mpi_io( 'use_upstream_for_tke', use_upstream_for_tke )
       CALL rrd_mpi_io( 'user_module_enabled', user_module_enabled )
       CALL rrd_mpi_io( 'v_bulk', v_bulk )
       CALL rrd_mpi_io_global_array( 'v_init', v_init )
       CALL rrd_mpi_io( 'v_max', v_max )
       CALL rrd_mpi_io_global_array( 'v_max_ijk', v_max_ijk )
       CALL rrd_mpi_io_global_array( 'vg', vg )
       CALL rrd_mpi_io( 'vg_surface', vg_surface )
       CALL rrd_mpi_io_global_array( 'vg_vertical_gradient', vg_vertical_gradient )
       CALL rrd_mpi_io_global_array( 'vg_vertical_gradient_level', vg_vertical_gradient_level )
       CALL rrd_mpi_io_global_array( 'vg_vertical_gradient_level_ind', vg_vertical_gradient_level_ind )
       CALL rrd_mpi_io( 'virtual_flight', virtual_flight )
       CALL rrd_mpi_io_global_array( 'volume_flow_area', volume_flow_area )
       CALL rrd_mpi_io_global_array( 'volume_flow_initial', volume_flow_initial )
       CALL rrd_mpi_io( 'w_max', w_max )
       CALL rrd_mpi_io_global_array( 'w_max_ijk', w_max_ijk )
       CALL rrd_mpi_io( 'wall_adjustment', wall_adjustment )
       CALL rrd_mpi_io_global_array( 'wall_heatflux', wall_heatflux )
       CALL rrd_mpi_io_global_array( 'wall_humidityflux', wall_humidityflux )
       CALL rrd_mpi_io_global_array( 'wall_scalarflux', wall_scalarflux )
       CALL rrd_mpi_io( 'y_shift', y_shift )
       CALL rrd_mpi_io( 'z0h_factor', z0h_factor )
       CALL rrd_mpi_io( 'zeta_max', zeta_max )
       CALL rrd_mpi_io( 'zeta_min', zeta_min )
       CALL rrd_mpi_io_global_array( 'z_i', z_i )
!
!--    Read global variables from surface_data_output_mod
       IF ( surface_output )  CALL surface_data_output_rrd_global
!
!--    Read global variables from of other modules
       CALL module_interface_rrd_global

!
!--    Close restart file
       CALL rd_mpi_io_close

    ENDIF

    CALL location_message( 'read global restart data', 'finished' )

 END SUBROUTINE rrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads values of global control variables from spinup restart-file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_global_spinup

!
!-- Open the MPI-IO restart file.
    CALL rd_mpi_io_open( 'read', 'SPINUPIN' // TRIM( coupling_char ),                              &
                         open_for_global_io_only = .TRUE. )

    CALL rrd_mpi_io( 'nx', nx )
    CALL rrd_mpi_io( 'ny', ny )
!
!-- Close restart file
    CALL rd_mpi_io_close
!
!-- Set x-/y-on-file dimensions to enable cyclic fill capability.
    nx_on_file = nx
    ny_on_file = ny

 END SUBROUTINE rrd_global_spinup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Skipping the global control variables from restart-file (binary format)
!> except some information which is required when reading restart data from a previous
!> run which used a smaller total domain or/and a different domain decomposition
!> (initializing_actions  == 'cyclic_fill').
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_read_parts_of_global


    CHARACTER (LEN=10) ::  version_on_file
    CHARACTER (LEN=20) ::  bc_lr_on_file
    CHARACTER (LEN=20) ::  bc_ns_on_file
    CHARACTER (LEN=20) ::  momentum_advec_check
    CHARACTER (LEN=20) ::  scalar_advec_check
    CHARACTER (LEN=1)  ::  cdum

    INTEGER(iwp) ::  max_pr_user_on_file
    INTEGER(iwp) ::  nz_on_file
    INTEGER(iwp) ::  statistic_regions_on_file
    INTEGER(iwp) ::  tmp_sr

    LOGICAL ::  neutral_check

    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  hom_sum_on_file
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  hom_on_file


    IF ( TRIM( restart_data_format_input ) == 'fortran_binary' )  THEN
!
!--    Input in Fortran binary format
       CALL check_open( 13 )

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file

!
!--    Read number of PEs and horizontal index bounds of all PEs used in previous run
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'numprocs' )  THEN
          WRITE( message_string, * ) 'numprocs not found in data from prior run on PE ', myid
          CALL message( 'rrd_read_parts_of_global', 'PAC0269', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  numprocs_previous_run

       IF ( .NOT. ALLOCATED( hor_index_bounds_previous_run ) )  THEN
          ALLOCATE( hor_index_bounds_previous_run(4,0:numprocs_previous_run-1) )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'hor_index_bounds' )  THEN
          WRITE( message_string, * ) 'hor_index_bounds not found in data ',    &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_read_parts_of_global', 'PAC0270', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  hor_index_bounds_previous_run

!
!--    Read vertical number of gridpoints and number of different areas used for computing
!--    statistics. Allocate arrays depending on these values, which are needed for the following
!--    read instructions.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'nz' )  THEN
          message_string = 'nz not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PAC0347', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  nz_on_file
       IF ( nz_on_file /= nz )  THEN
          WRITE( message_string, * ) 'mismatch concerning number of ',         &
                                     'gridpoints along z:',                    &
                                     '&nz on file    = "', nz_on_file, '"',    &
                                     '&nz from run   = "', nz, '"'
          CALL message( 'rrd_read_parts_of_global', 'PAC0275', 1, 2, 0, 6, 0 )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'max_pr_user' )  THEN
          message_string = 'max_pr_user not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PAC0276', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  max_pr_user_on_file
       IF ( max_pr_user_on_file /= max_pr_user )  THEN
          WRITE( message_string, * ) 'number of user profiles on res',         &
                                     'tart data file differs from the ',       &
                                     'current run:&max_pr_user on file    = "',&
                                     max_pr_user_on_file, '"',                 &
                                     '&max_pr_user from run   = "',            &
                                     max_pr_user, '"'
          CALL message( 'rrd_read_parts_of_global', 'PAC0277', 0, 0, 0, 6, 0 )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'statistic_regions' )  THEN
          message_string = 'statistic_regions not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PAC0278', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  statistic_regions_on_file
       IF ( statistic_regions_on_file /= statistic_regions )  THEN
          WRITE( message_string, * ) 'statistic regions on restart data file ',&
                                     'differ from the current run:',           &
                                     '&statistic regions on file    = "',      &
                                     statistic_regions_on_file, '"',           &
                                     '&statistic regions from run   = "',      &
                                      statistic_regions, '"',                  &
                                     '&statistic data may be lost!'
          CALL message( 'rrd_read_parts_of_global', 'PAC0279', 0, 1, 0, 6, 0 )
          tmp_sr = MIN( statistic_regions_on_file, statistic_regions )
       ELSE
          tmp_sr = statistic_regions
       ENDIF

!
!--    Now read and check some control parameters and skip the rest
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO  WHILE ( restart_string(1:length) /= 'binary_version_local' )

          SELECT CASE ( restart_string(1:length) )

             CASE ( 'average_count_pr' )
                READ ( 13 )  average_count_pr
                IF ( average_count_pr /= 0 )  THEN
                   WRITE( message_string, * ) 'Inflow profiles not temporally averaged.',          &
                                              '&Averaging will be done now using',                 &
                                              average_count_pr, ' samples.'
                   CALL message( 'rrd_read_parts_of_global', 'PAC0280',         &
                                 0, 1, 0, 6, 0 )
                ENDIF

             CASE ( 'bc_lr' )
                READ ( 13 )  bc_lr_on_file
                IF ( TRIM( bc_lr_on_file ) /= 'cyclic' )  THEN
                   message_string = 'illegal setting of bc_lr = "cyclic" in the prerun'
                   CALL message( 'rrd_read_parts_of_global', 'PAC0281', 1, 2, 0, 6, 0 )
                ENDIF

             CASE ( 'bc_ns' )
                READ ( 13 )  bc_ns_on_file
                IF ( TRIM( bc_ns_on_file ) /= 'cyclic' )  THEN
                   message_string = 'illegal setting of bc_lr = "cyclic" in the prerun'
                   CALL message( 'rrd_read_parts_of_global', 'PAC0281', 1, 2, 0, 6, 0 )
                ENDIF

             CASE ( 'hom' )
                ALLOCATE( hom_on_file(0:nz+1,2,pr_max,0:statistic_regions_on_file) )
                READ ( 13 )  hom_on_file
                hom(:,:,:,0:tmp_sr) = hom_on_file(:,:,:,0:tmp_sr)
                DEALLOCATE( hom_on_file )

             CASE ( 'hom_sum' )
                ALLOCATE( hom_sum_on_file(0:nz+1,pr_max,0:statistic_regions_on_file) )
                READ ( 13 )  hom_sum_on_file
                hom_sum(:,:,0:tmp_sr) = hom_sum_on_file(:,:,0:tmp_sr)
                DEALLOCATE( hom_sum_on_file )

             CASE ( 'momentum_advec' )
                momentum_advec_check = momentum_advec
                READ ( 13 )  momentum_advec
                IF ( TRIM( momentum_advec_check ) /= TRIM( momentum_advec ) )  &
                THEN
                   message_string = 'momentum_advec = "' // TRIM( momentum_advec_check ) //        &
                                    '" of the restart run differs from momentum_advec = "' //      &
                                    TRIM( momentum_advec) // '" of the initial run'
                   CALL message( 'rrd_read_parts_of_global', 'PAC0282',         &
                                 1, 2, 0, 6, 0 )
                ENDIF

             CASE ( 'neutral' )
                neutral_check = neutral
                READ (13 )  neutral
                IF ( neutral_check  .NEQV.  neutral )  THEN
                   message_string = 'setting of neutral in pre-run differs from setting in main run'
                   CALL message( 'rrd_read_parts_of_global', 'PAC0283', 1, 2, 0, 6, 0 )
                ENDIF

             CASE ( 'nx' )
                READ ( 13 )  nx_on_file

             CASE ( 'ny' )
                READ ( 13 )  ny_on_file

             CASE ( 'ref_state' )
                READ ( 13 )  ref_state

             CASE ( 'scalar_advec' )
                scalar_advec_check = scalar_advec
                READ ( 13 )  scalar_advec
                IF ( TRIM( scalar_advec_check ) /= TRIM( scalar_advec ) )      &
                THEN
                   message_string = 'scalar_advec = "' // TRIM( scalar_advec_check ) //            &
                                    '" of the restart run differs from scalar_advec = "' //        &
                                    TRIM( scalar_advec) // '" of the initial run'
                   CALL message( 'rrd_read_parts_of_global', 'PAC0284', 1, 2, 0, 6, 0 )
                ENDIF

             CASE DEFAULT

                READ ( 13 )  cdum

          END SELECT

          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO

       CALL close_file( 13 )

    ELSEIF ( restart_data_format_input(1:3) == 'mpi' )  THEN
!
!--    Open the MPI-IO restart file.
       CALL rd_mpi_io_open( 'read', 'BININ' // TRIM( coupling_char ),                              &
                            open_for_global_io_only = .TRUE. )

!
!--    Read vertical number of gridpoints and number of different areas used for computing
!--    statistics. Allocate arrays depending on these values, which are required for the following
!--    read instructions.
       CALL rrd_mpi_io( 'nz', nz_on_file )
       IF ( nz_on_file /= nz )  THEN
          WRITE( message_string, * ) 'mismatch concerning number of gridpoints along z:',          &
                                     '&nz on file    = "', nz_on_file, '"',                        &
                                     '&nz from run   = "', nz, '"'
          CALL message( 'rrd_read_parts_of_global', 'PAC0275', 1, 2, 0, 6, 0 )
       ENDIF

       CALL rrd_mpi_io( 'max_pr_user', max_pr_user_on_file )
       IF ( max_pr_user_on_file /= max_pr_user )  THEN
          WRITE( message_string, * ) 'number of user profiles on restart data file differs from ', &
                                     'the current run:&max_pr_user on file    = "',                &
                                     max_pr_user_on_file, '" &max_pr_user from run   = "',         &
                                     max_pr_user, '"'
          CALL message( 'rrd_read_parts_of_global', 'PAC0277', 0, 0, 0, 6, 0 )
       ENDIF

       CALL rrd_mpi_io( 'statistic_regions', statistic_regions_on_file )
       IF ( statistic_regions_on_file /= statistic_regions )  THEN
          WRITE( message_string, * ) 'statistic regions on restart data file ',&
                                     'differ from the current run:',           &
                                     '&statistic regions on file    = "',      &
                                     statistic_regions_on_file, '"',           &
                                     '&statistic regions from run   = "',      &
                                      statistic_regions, '"',                  &
                                     '&statistic data may be lost!'
          CALL message( 'rrd_read_parts_of_global', 'PAC0279', 0, 1, 0, 6, 0 )
          tmp_sr = MIN( statistic_regions_on_file, statistic_regions )
       ELSE
          tmp_sr = statistic_regions
       ENDIF

!
!--    Now read and check some control parameters and skip the rest.
       CALL rrd_mpi_io( 'average_count_pr', average_count_pr )
       IF ( average_count_pr /= 0 )  THEN
          WRITE( message_string, * ) 'Inflow profiles not temporally averaged.&Averaging ',        &
                                     'will be done now using', average_count_pr, ' samples.'
          CALL message( 'rrd_read_parts_of_global', 'PAC0280', 0, 1, 0, 6, 0 )
       ENDIF

       ALLOCATE( hom_on_file(0:nz+1,2,pr_max,0:statistic_regions_on_file) )
       CALL rrd_mpi_io_global_array( 'hom', hom_on_file )
       hom(:,:,:,0:tmp_sr) = hom_on_file(:,:,:,0:tmp_sr)
       DEALLOCATE( hom_on_file )

       ALLOCATE( hom_sum_on_file(0:nz+1,pr_max,0:statistic_regions_on_file) )
       CALL rrd_mpi_io_global_array( 'hom_sum', hom_sum_on_file )
       hom_sum(:,:,0:tmp_sr) = hom_sum_on_file(:,:,0:tmp_sr)
       DEALLOCATE( hom_sum_on_file )

       momentum_advec_check = momentum_advec
       CALL rrd_mpi_io( 'momentum_advec', momentum_advec )
       IF ( TRIM( momentum_advec_check ) /= TRIM( momentum_advec ) )  THEN
          message_string = 'momentum_advec = "' // TRIM( momentum_advec_check ) //                 &
                           '" of the restart run differs from momentum_advec = "' //               &
                           TRIM( momentum_advec) // '" of the initial run'
          CALL message( 'rrd_read_parts_of_global', 'PAC0282', 1, 2, 0, 6, 0 )
       ENDIF

       CALL rrd_mpi_io( 'bc_lr', bc_lr_on_file )
       IF ( TRIM( bc_lr_on_file ) /= 'cyclic' )  THEN
          message_string = 'illegal setting of bc_lr = "cyclic" in the prerun'
          CALL message( 'rrd_read_parts_of_global', 'PAC0281', 1, 2, 0, 6, 0 )
       ENDIF
       CALL rrd_mpi_io( 'bc_ns', bc_ns_on_file )
       IF ( TRIM( bc_ns_on_file ) /= 'cyclic'  )  THEN
          message_string = 'illegal setting of bc_ns = "cyclic" in the prerun'
          CALL message( 'rrd_read_parts_of_global', 'PAC0281', 1, 2, 0, 6, 0 )
       ENDIF

       scalar_advec_check = scalar_advec
       CALL rrd_mpi_io( 'scalar_advec', scalar_advec )
       IF ( TRIM( scalar_advec_check ) /= TRIM( scalar_advec ) )  THEN
          message_string = 'scalar_advec = "' // TRIM( scalar_advec_check ) //                     &
                           '" of the restart run differs from scalar_advec = "' //                 &
                           TRIM( momentum_advec) // '" of the initial run'
          CALL message( 'rrd_read_parts_of_global', 'PAC0284', 1, 2, 0, 6, 0 )
       ENDIF

       neutral_check = neutral
       CALL rrd_mpi_io( 'neutral', neutral )
       IF ( neutral_check  .NEQV.  neutral )  THEN
          message_string = 'setting of neutral in pre-run differs from setting in main run'
          CALL message( 'rrd_read_parts_of_global', 'PAC0283', 1, 2, 0, 6, 0 )
       ENDIF

       CALL rrd_mpi_io( 'nx', nx_on_file )
       CALL rrd_mpi_io( 'ny', ny_on_file )
       CALL rrd_mpi_io_global_array( 'ref_state', ref_state )

!
!--    Close restart file
       CALL rd_mpi_io_close

    ENDIF

!
!-- Calculate the temporal average of vertical profiles, if neccessary
    IF ( average_count_pr /= 0 )  THEN
       hom_sum = hom_sum / REAL( average_count_pr, KIND=wp )
    ENDIF

 END SUBROUTINE rrd_read_parts_of_global


! Description:
! ------------
!> Reads processor (subdomain) specific data of variables and arrays from restart file
!> (binary format).
!------------------------------------------------------------------------------!
 SUBROUTINE rrd_local


    CHARACTER (LEN=7)  ::  myid_char_save
    CHARACTER (LEN=10) ::  binary_version_local
    CHARACTER (LEN=10) ::  version_on_file
    CHARACTER (LEN=20) ::  tmp_name               !< temporary variable

    INTEGER(iwp) ::  files_to_be_opened  !<
    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  myid_on_file        !<
    INTEGER(iwp) ::  numprocs_on_file    !<
    INTEGER(iwp) ::  nxlc                !<
    INTEGER(iwp) ::  nxlf                !<
    INTEGER(iwp) ::  nxlpr               !<
    INTEGER(iwp) ::  nxl_on_file         !<
    INTEGER(iwp) ::  nxrc                !<
    INTEGER(iwp) ::  nxrf                !<
    INTEGER(iwp) ::  nxrpr               !<
    INTEGER(iwp) ::  nxr_on_file         !<
    INTEGER(iwp) ::  nync                !<
    INTEGER(iwp) ::  nynf                !<
    INTEGER(iwp) ::  nynpr               !<
    INTEGER(iwp) ::  nyn_on_file         !<
    INTEGER(iwp) ::  nysc                !<
    INTEGER(iwp) ::  nysf                !<
    INTEGER(iwp) ::  nyspr               !<
    INTEGER(iwp) ::  nys_on_file         !<
    INTEGER(iwp) ::  nzb_on_file         !<
    INTEGER(iwp) ::  nzt_on_file         !<
    INTEGER(iwp) ::  offset_x            !<
    INTEGER(iwp) ::  offset_y            !<
    INTEGER(iwp) ::  shift_x             !<
    INTEGER(iwp) ::  shift_y             !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  file_list       !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  overlap_count   !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nxlfa      !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nxrfa      !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nynfa      !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nysfa      !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  offset_xa  !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  offset_ya  !<

    INTEGER(isp), DIMENSION(:,:),   ALLOCATABLE ::  tmp_2d_id_random   !< temporary array for storing random generator data
    INTEGER(isp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp_2d_seq_random  !< temporary array for storing random generator data

    LOGICAL ::  array_found                      !<
    LOGICAL ::  found                            !<

    REAL(wp), DIMENSION(:,:),   ALLOCATABLE   ::  tmp_2d         !< temporary array for storing 2D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3d         !< temporary array for storing 3D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3d_non_standard !< temporary array for storing 3D data
                                                                 !< with non standard dimensions

!
!-- Read data from previous model run.
    CALL cpu_log( log_point_s(14), 'read-restart-data-local', 'start' )

    CALL location_message( 'reading local restart data', 'start' )

    IF ( TRIM( restart_data_format_input ) == 'fortran_binary' )  THEN
!
!--    Input in Fortran binary format

!
!--    Allocate temporary buffer arrays. In previous versions, they were
!--    declared as automated arrays, causing memory problems when these
!--    were allocate on stack.
       ALLOCATE( nxlfa(numprocs_previous_run,1000) )
       ALLOCATE( nxrfa(numprocs_previous_run,1000) )
       ALLOCATE( nynfa(numprocs_previous_run,1000) )
       ALLOCATE( nysfa(numprocs_previous_run,1000) )
       ALLOCATE( offset_xa(numprocs_previous_run,1000) )
       ALLOCATE( offset_ya(numprocs_previous_run,1000) )

!
!--    Check which of the restart files contain data needed for the subdomain
!--    of this PE
       files_to_be_opened = 0

       DO  i = 1, numprocs_previous_run
!
!--       Store array bounds of the previous run ("pr") in temporary scalars
          nxlpr = hor_index_bounds_previous_run(1,i-1)
          nxrpr = hor_index_bounds_previous_run(2,i-1)
          nyspr = hor_index_bounds_previous_run(3,i-1)
          nynpr = hor_index_bounds_previous_run(4,i-1)

!
!--       Determine the offsets. They may be non-zero in case that the total domain
!--       on file is smaller than the current total domain.
          offset_x = ( nxl / ( nx_on_file + 1 ) ) * ( nx_on_file + 1 )
          offset_y = ( nys / ( ny_on_file + 1 ) ) * ( ny_on_file + 1 )

!
!--       Start with this offset and then check, if the subdomain on file
!--       matches another time(s) in the current subdomain by shifting it
!--       for nx_on_file+1, ny_on_file+1 respectively

          shift_y = 0
          j       = 0
          DO WHILE (  nyspr+shift_y <= nyn-offset_y )

             IF ( nynpr+shift_y >= nys-offset_y ) THEN

                shift_x = 0
                DO WHILE ( nxlpr+shift_x <= nxr-offset_x )

                   IF ( nxrpr+shift_x >= nxl-offset_x ) THEN
                      j = j +1
                      IF ( j > 1000 )  THEN
!
!--                      Array bound exceeded
                         message_string = 'data from subdomain of previous run mapped more ' //    &
                                          'than 1000 times'
                         CALL message( 'rrd_local', 'PAC0285', 2, 2, -1, 6, 1 )
                      ENDIF

                      IF ( j == 1 )  THEN
                         files_to_be_opened = files_to_be_opened + 1
                         file_list(files_to_be_opened) = i-1
                      ENDIF

                      offset_xa(files_to_be_opened,j) = offset_x + shift_x
                      offset_ya(files_to_be_opened,j) = offset_y + shift_y
!
!--                   Index bounds of overlapping data
                      nxlfa(files_to_be_opened,j) = MAX( nxl-offset_x-shift_x, nxlpr )
                      nxrfa(files_to_be_opened,j) = MIN( nxr-offset_x-shift_x, nxrpr )
                      nysfa(files_to_be_opened,j) = MAX( nys-offset_y-shift_y, nyspr )
                      nynfa(files_to_be_opened,j) = MIN( nyn-offset_y-shift_y, nynpr )

                   ENDIF

                   shift_x = shift_x + ( nx_on_file + 1 )
                ENDDO

             ENDIF

             shift_y = shift_y + ( ny_on_file + 1 )
          ENDDO

          IF ( j > 0 )  overlap_count(files_to_be_opened) = j

       ENDDO

!
!--    Save the id-string of the current process, since myid_char may now be used
!--    to open files created by PEs with other id.
       myid_char_save = myid_char

       IF ( files_to_be_opened /= 1  .OR.  numprocs /= numprocs_previous_run )  THEN
          message_string = 'number of PEs or virtual PE-grid changed in restart run'
          CALL message( 'rrd_local', 'PAC0286', 0, 0, 0, 6, 0 )
          IF ( debug_output )  THEN
             IF ( files_to_be_opened <= 120 )  THEN
                WRITE( debug_string, '(2A,1X,120(I6.6,1X))' )                                      &
                  'number of PEs or virtual PE-grid changed in restart run.  PE will read from ',  &
                  'files ', file_list(1:files_to_be_opened)
             ELSE
                WRITE( debug_string, '(3A,1X,120(I6.6,1X),A)' )                                    &
                  'number of PEs or virtual PE-grid changed in restart run.  PE will read from ',  &
                  'files ', file_list(1:120), '... and more'
             ENDIF
             CALL debug_message( 'rrd_local', 'info' )
          ENDIF
       ENDIF

!
!--    Read data from all restart files determined above
       DO  i = 1, files_to_be_opened

          j = file_list(i)
!
!--       Set the filename (underscore followed by four digit processor id)
          WRITE (myid_char,'(''_'',I6.6)')  j

!
!--       Open the restart file. If this file has been created by PE0 (_000000),
!--       the global variables at the beginning of the file have to be skipped
!--       first.
          CALL check_open( 13 )
          IF ( j == 0 )  CALL rrd_skip_global

!
!--       First compare the version numbers
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)
          READ ( 13 )  version_on_file

          binary_version_local = '22.04'
          IF ( TRIM( version_on_file ) /= TRIM( binary_version_local ) )  THEN
             WRITE( message_string, * ) 'version mismatch concerning ',                            &
                                        'binary_version_local:',                                   &
                                        '&version on file    = "', TRIM( version_on_file ), '"',   &
                                        '&version in program = "', TRIM( binary_version_local ), '"'
             CALL message( 'rrd_local', 'PAC0287', 1, 2, 0, 6, 0 )
          ENDIF

!
!--       Read number of processors, processor-id, and array ranges.
!--       Compare the array ranges with those stored in the index bound array.
          READ ( 13 )  numprocs_on_file, myid_on_file, nxl_on_file, nxr_on_file, nys_on_file,      &
                       nyn_on_file, nzb_on_file, nzt_on_file

          IF ( nxl_on_file /= hor_index_bounds_previous_run(1,j) )  THEN
             WRITE( message_string, * ) 'problem with index bound nxl on ',                        &
                                        'restart file "', myid_char, '"',                          &
                                        '&nxl = ', nxl_on_file, ' but it should be',               &
                                        '&= ', hor_index_bounds_previous_run(1,j),                 &
                                        '&from the index bound information array'
             CALL message( 'rrd_local', 'PAC0288', 2, 2, -1, 6, 1 )
          ENDIF

          IF ( nxr_on_file /= hor_index_bounds_previous_run(2,j) )  THEN
              WRITE( message_string, * ) 'problem with index bound nxr on ',                       &
                                         'restart file "', myid_char, '"'  ,                       &
                                         ' nxr = ', nxr_on_file, ' but it should be',              &
                                         ' = ', hor_index_bounds_previous_run(2,j),                &
                                         ' from the index bound information array'
             CALL message( 'rrd_local', 'PAC0289', 2, 2, -1, 6, 1 )

          ENDIF

          IF ( nys_on_file /= hor_index_bounds_previous_run(3,j) )  THEN
             WRITE( message_string, * ) 'problem with index bound nys on ',                        &
                                        'restart file "', myid_char, '"',                          &
                                        '&nys = ', nys_on_file, ' but it should be',               &
                                        '&= ', hor_index_bounds_previous_run(3,j),                 &
                                        '&from the index bound information array'
             CALL message( 'rrd_local', 'PAC0290', 2, 2, -1, 6, 1 )
          ENDIF

          IF ( nyn_on_file /= hor_index_bounds_previous_run(4,j) )  THEN
             WRITE( message_string, * ) 'problem with index bound nyn on ',                        &
                                        'restart file "', myid_char, '"',                          &
                                        '&nyn = ', nyn_on_file, ' but it should be',               &
                                        '&= ', hor_index_bounds_previous_run(4,j),                 &
                                        '&from the index bound information array'
             CALL message( 'rrd_local', 'PAC0291', 2, 2, -1, 6, 1 )
          ENDIF

          IF ( nzb_on_file /= nzb )  THEN
             WRITE( message_string, * ) 'mismatch between actual data and data ',                  &
                                        'from prior run on PE ', myid,                             &
                                        '&nzb on file = ', nzb_on_file,                            &
                                        '&nzb         = ', nzb
             CALL message( 'rrd_local', 'PAC0292', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( nzt_on_file /= nzt )  THEN
             WRITE( message_string, * ) 'mismatch between actual data and data ',                  &
                                        'from prior run on PE ', myid,                             &
                                        '&nzt on file = ', nzt_on_file,                            &
                                        '&nzt         = ', nzt
             CALL message( 'rrd_local', 'PAC0293', 1, 2, 0, 6, 0 )
          ENDIF

!
!--       Allocate temporary arrays sized as the arrays on the restart file
          ALLOCATE( tmp_2d(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp),      &
                    tmp_3d(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,                               &
                           nxl_on_file-nbgp:nxr_on_file+nbgp) )

!
!--       Read arrays
!--       ATTENTION: If the following read commands have been altered, the
!--       ---------- version number of the variable binary_version_local must
!--                  be altered, too. Furthermore, the output list of arrays in
!--                  wrd_write_local must also be altered
!--                  accordingly.
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)


!
!--       Loop over processor specific field data
          DO  WHILE ( restart_string(1:length) /= '*** end ***' )

!
!--          Map data on file as often as needed (data are read only for k=1)
             DO  k = 1, overlap_count(i)

                found = .FALSE.

!
!--             Get the index range of the subdomain on file which overlap with
!--             the current subdomain
                nxlf = nxlfa(i,k)
                nxlc = nxlfa(i,k) + offset_xa(i,k)
                nxrf = nxrfa(i,k)
                nxrc = nxrfa(i,k) + offset_xa(i,k)
                nysf = nysfa(i,k)
                nysc = nysfa(i,k) + offset_ya(i,k)
                nynf = nynfa(i,k)
                nync = nynfa(i,k) + offset_ya(i,k)


                SELECT CASE ( restart_string(1:length) )

                   CASE ( 'ghf_av' )
                      IF ( .NOT. ALLOCATED( ghf_av ) )  THEN
                         ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      ghf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'e' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      e(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'e_av' )
                      IF ( .NOT. ALLOCATED( e_av ) )  THEN
                         ALLOCATE( e_av(nzb:nzt+1,nys-nbgp:nyn+nbgp,                               &
                                        nxl-nbgp:nxr+nbgp) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      e_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'kh' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      kh(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'kh_av' )
                      IF ( .NOT. ALLOCATED( kh_av ) )  THEN
                         ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      kh_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'km' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      km(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'km_av' )
                      IF ( .NOT. ALLOCATED( km_av ) )  THEN
                         ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      km_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'lpt_av' )
                      IF ( .NOT. ALLOCATED( lpt_av ) )  THEN
                         ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      lpt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'lwp_av' )
                      IF ( .NOT. ALLOCATED( lwp_av ) )  THEN
                         ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      lwp_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'ol_av' )
                      IF ( .NOT. ALLOCATED( ol_av ) )  THEN
                         ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      ol_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'p' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      p(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'p_av' )
                      IF ( .NOT. ALLOCATED( p_av ) )  THEN
                         ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      p_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'pres_drag_x_av' )
                      IF ( .NOT. ALLOCATED( pres_drag_x_av ) )  THEN
                         ALLOCATE( pres_drag_x_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      pres_drag_x_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                    &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'pres_drag_y_av' )
                      IF ( .NOT. ALLOCATED( pres_drag_y_av ) )  THEN
                         ALLOCATE( pres_drag_y_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      pres_drag_y_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                    &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'pt' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      pt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'pt_av' )
                      IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                         ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      pt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'q' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      q(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'q_av' )
                      IF ( .NOT. ALLOCATED( q_av ) )  THEN
                         ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      q_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'ql' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      ql(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'ql_av' )
                      IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                         ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      ql_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'qsurf_av' )
                      IF ( .NOT. ALLOCATED( qsurf_av ) )  THEN
                         ALLOCATE( qsurf_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      qsurf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'qsws_av' )
                      IF ( .NOT. ALLOCATED( qsws_av ) )  THEN
                         ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      qsws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                          &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'qv_av' )
                      IF ( .NOT. ALLOCATED( qv_av ) )  THEN
                         ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      qv_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'r_a_av' )
                      IF ( .NOT. ALLOCATED( r_a_av ) )  THEN
                         ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      r_a_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'random_iv' )  ! still unresolved issue
                      IF ( k == 1 )  READ ( 13 )  random_iv
                      IF ( k == 1 )  READ ( 13 )  random_iy

                   CASE ( 'seq_random_array' )
                      ALLOCATE( tmp_2d_id_random(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )
                      ALLOCATE( tmp_2d_seq_random(5,nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )
                      IF ( .NOT. ALLOCATED( id_random_array ) )  THEN
                         ALLOCATE( id_random_array(nys:nyn,nxl:nxr) )
                      ENDIF
                      IF ( .NOT. ALLOCATED( seq_random_array ) )  THEN
                         ALLOCATE( seq_random_array(5,nys:nyn,nxl:nxr) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d_id_random
                      IF ( k == 1 )  READ ( 13 )  tmp_2d_seq_random
                      id_random_array(nysc:nync,nxlc:nxrc) = tmp_2d_id_random(nysf:nynf,nxlf:nxrf)
                      seq_random_array(:,nysc:nync,nxlc:nxrc) = tmp_2d_seq_random(:,nysf:nynf,nxlf:nxrf)
                      DEALLOCATE( tmp_2d_id_random, tmp_2d_seq_random )

                   CASE ( 's' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      s(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 's_av' )
                      IF ( .NOT. ALLOCATED( s_av ) )  THEN
                         ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg))
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      s_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'shf_av' )
                      IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                         ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      shf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                           &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'ssurf_av' )
                      IF ( .NOT. ALLOCATED( ssurf_av ) )  THEN
                         ALLOCATE( ssurf_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      ssurf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'ssws_av' )
                      IF ( .NOT. ALLOCATED( ssws_av ) )  THEN
                         ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      ssws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                          &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'ts_av' )
                      IF ( .NOT. ALLOCATED( ts_av ) )  THEN
                         ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      ts_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                            &
                           tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'tsurf_av' )
                      IF ( .NOT. ALLOCATED( tsurf_av ) )  THEN
                         ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      tsurf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'u' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      u(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'u_av' )
                      IF ( .NOT. ALLOCATED( u_av ) )  THEN
                         ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      u_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'us_av' )
                      IF ( .NOT. ALLOCATED( us_av ) )  THEN
                         ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      us_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                            &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'v' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      v(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'v_av' )
                      IF ( .NOT. ALLOCATED( v_av ) )  THEN
                         ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'vpt' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      vpt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'vpt_av' )
                      IF ( .NOT. ALLOCATED( vpt_av ) )  THEN
                         ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      vpt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'w' )
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      w(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'w_av' )
                      IF ( .NOT. ALLOCATED( w_av ) )  THEN
                         ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_3d
                      w_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                         tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'z0_av' )
                      IF ( .NOT. ALLOCATED( z0_av ) )  THEN
                         ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      z0_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                            &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'z0h_av' )
                      IF ( .NOT. ALLOCATED( z0h_av ) )  THEN
                         ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      z0h_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                           &
                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE ( 'z0q_av' )
                      IF ( .NOT. ALLOCATED( z0q_av ) )  THEN
                         ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                      ENDIF
                      IF ( k == 1 )  READ ( 13 )  tmp_2d
                      z0q_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                           &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                   CASE DEFAULT

!
!--                   Read restart data of surfaces
                      IF ( .NOT. found )  CALL surface_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf,&
                                                                  nxr_on_file, nynf, nyn_on_file,  &
                                                                  nysf, nysc, nys_on_file, found )
!
!--                   Read restart data of surface_data_output_mod.  Surface data do not need
!--                   overlap data, so do not pass these information.
                      IF ( .NOT. found )  CALL surface_data_output_rrd_local( found )
!
!--                   Read restart data of other modules
                      IF ( .NOT. found ) CALL module_interface_rrd_local(                          &
                                                                  k, nxlf, nxlc, nxl_on_file, nxrf,&
                                                                  nxrc, nxr_on_file, nynf, nync,   &
                                                                  nyn_on_file, nysf, nysc,         &
                                                                  nys_on_file, tmp_2d, tmp_3d,     &
                                                                  found )


                      IF ( .NOT. found )  THEN
                         WRITE( message_string, * ) 'unknown variable named "',                    &
                                                    restart_string(1:length),                      &
                                                   '" found in subdomain data ',                   &
                                                   'from prior run on PE ', myid
                         CALL message( 'rrd_local', 'PAC0274', 1, 2, 0, 6, 0 )

                      ENDIF

                END SELECT

             ENDDO ! overlaploop

!
!--          Deallocate non standard array needed for specific variables only
             IF ( ALLOCATED( tmp_3d_non_standard ) )  DEALLOCATE( tmp_3d_non_standard )

!
!--          Read next character string
             READ ( 13 )  length
             READ ( 13 )  restart_string(1:length)

          ENDDO ! dataloop
!
!--       Close the restart file
          CALL close_file( 13 )

          DEALLOCATE( tmp_2d, tmp_3d )

       ENDDO  ! loop over restart files
!
!--    Deallocate temporary buffer arrays
       DEALLOCATE( nxlfa )
       DEALLOCATE( nxrfa )
       DEALLOCATE( nynfa )
       DEALLOCATE( nysfa )
       DEALLOCATE( offset_xa )
       DEALLOCATE( offset_ya )
!
!--    Restore the original filename for the restart file to be written
       myid_char = myid_char_save


    ELSEIF ( restart_data_format_input(1:3) == 'mpi' )  THEN

!
!--    Read local restart data using MPI-IO
!
!--    Open the MPI-IO restart file.
       CALL rd_mpi_io_open( 'read', 'BININ' // TRIM( coupling_char ) )

!
!--    Note, restart input of time-averaged quantities is skipped in case of cyclic-fill
!--    initialization. This case, input of time-averaged data is useless and can lead to faulty
!--    averaging.
       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'ghf_av' , found = array_found )
          IF ( array_found )  THEN
             IF (.NOT. ALLOCATED( ghf_av ) )  ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ghf_av', ghf_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'e', e )

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'e_av' , found = array_found )
          IF ( array_found  )  THEN
             IF ( .NOT. ALLOCATED( e_av ) )  ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'e_av', e_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'kh', kh )

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'kh_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( kh_av ) )  ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'kh_av', kh_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'km' , km)

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'km_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( km_av ) )  ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'km_av', km_av )
          ENDIF
!
          CALL rd_mpi_io_check_array( 'lpt_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( lpt_av ) )  ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'lpt_av', lpt_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'lwp_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( lwp_av ) )  ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'lwp_av', lwp_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'ol_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( ol_av ) )  ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ol_av', ol_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'p', p)

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'p_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( p_av ) )  ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'p_av', p_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'pres_drag_x_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( pres_drag_x_av ) )  ALLOCATE( pres_drag_x_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'pres_drag_x_av', pres_drag_x_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'pres_drag_y_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( pres_drag_y_av ) )  ALLOCATE( pres_drag_y_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'pres_drag_y_av', pres_drag_y_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'pt', pt)

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'pt_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( pt_av ) )  ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'pt_av', pt_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'q' , found = array_found )
       IF ( array_found )  THEN
          CALL rrd_mpi_io( 'q', q )
       ENDIF

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'q_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( q_av ) )  ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'q_av', q_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'ql' , found = array_found )
       IF ( array_found )  THEN
          CALL rrd_mpi_io( 'ql', ql )
       ENDIF

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'ql_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( ql_av ) )  ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ql_av', ql_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'qsurf_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( qsurf_av ) )  ALLOCATE( qsurf_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'qsurf_av', qsurf_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'qsws_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( qsws_av ) )  ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'qsws_av', qsws_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'qv_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( qv_av ) )  ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'qv_av', qv_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'r_a_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( r_a_av ) )  ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'r_a_av', r_a_av )
          ENDIF
       ENDIF

!
!--    ATTENTION: The random seeds are global data! If independent values for every PE are required,
!--    the general approach of PE indendent restart will be lost. That means that in general the
!--    parallel random number generator in random_generator_parallel_mod should be used!
       CALL rrd_mpi_io_global_array( 'random_iv', random_iv )
       CALL rrd_mpi_io( 'random_iy', random_iy )

       CALL rd_mpi_io_check_array( 'id_random_array' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( id_random_array ) )  ALLOCATE( id_random_array(nys:nyn,nxl:nxr) )
          IF ( .NOT. ALLOCATED( seq_random_array ) )  ALLOCATE( seq_random_array(5,nys:nyn,nxl:nxr) )
          CALL rrd_mpi_io( 'id_random_array', id_random_array)
          DO  i = 1, SIZE( seq_random_array, 1 )
             WRITE( tmp_name, '(A,I2.2)' )  'seq_random_array', i
             CALL rrd_mpi_io( TRIM(tmp_name), seq_random_array(i,:,:) )
          ENDDO
       ENDIF

       CALL rd_mpi_io_check_array( 's' , found = array_found )
       IF ( array_found )  THEN
          CALL rrd_mpi_io( 's', s )
       ENDIF

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 's_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( s_av ) )  ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 's_av', s_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'shf_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( shf_av ) )  ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'shf_av', shf_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'ssurf_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( ssurf_av ) )  ALLOCATE( ssurf_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ssurf_av', ssurf_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'ssws_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( ssws_av ) )  ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ssws_av', ssws_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'ts_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( ts_av ) )  ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ts_av', ts_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'tsurf_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( tsurf_av ) )  ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'tsurf_av', tsurf_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'u', u)

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'u_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( u_av ) )  ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'u_av', u_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'us_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( us_av ) )  ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'us_av', us_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'v', v )

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'v_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( v_av ) )  ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'v_av', v_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'vpt' , found = array_found )
       IF ( array_found )  THEN
          CALL rrd_mpi_io( 'vpt',  vpt)
       ENDIF

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'vpt_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( vpt_av ) )  ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'vpt_av', vpt_av )
          ENDIF
       ENDIF

       CALL rrd_mpi_io( 'w', w)

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'w_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( w_av ) )  ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'w_av', w_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'z0_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( z0_av ) )  ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'z0_av', z0_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'z0h_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( z0h_av ) )  ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'z0h_av', z0h_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'z0q_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( z0q_av ) )  ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'z0q_av', z0q_av )
          ENDIF
       ENDIF

!
!--    Read restart data of surfaces
       CALL surface_rrd_local
!
!--    Read restart data of surface_data_output_mod
       IF ( surface_output )  CALL surface_data_output_rrd_local
!
!--    Read restart data of other modules
       CALL module_interface_rrd_local

!
!--    Close restart file
       CALL rd_mpi_io_close

    ENDIF

    CALL location_message( 'reading local restart data', 'finished' )
!
!-- End of time measuring for reading binary data
    CALL cpu_log( log_point_s(14), 'read-restart-data-local', 'stop' )

 END SUBROUTINE rrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads local spinup data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_local_spinup
!
!-- Open the MPI-IO restart file.
    CALL rd_mpi_io_open( 'read', 'SPINUPIN' // TRIM( coupling_char ) )
!
!-- Read restart data of surfaces
    CALL surface_rrd_local
!
!-- Read restart data of other modules
    CALL module_interface_rrd_local_spinup
!
!-- Close restart file
    CALL rd_mpi_io_close

 END SUBROUTINE rrd_local_spinup

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Skipping the global control variables from restart-file (binary format)
!------------------------------------------------------------------------------!
 SUBROUTINE rrd_skip_global


    CHARACTER (LEN=1) ::  cdum


    READ ( 13 )  length
    READ ( 13 )  restart_string(1:length)

    DO  WHILE ( restart_string(1:length) /= 'binary_version_local' )

       READ ( 13 )  cdum
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

    ENDDO

    BACKSPACE ( 13 )
    BACKSPACE ( 13 )


 END SUBROUTINE rrd_skip_global


 END MODULE read_restart_data_mod
