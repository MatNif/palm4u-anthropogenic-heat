!> @file write_restart_data_mod.f90
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
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Writes restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 MODULE write_restart_data_mod

    USE arrays_3d,                                                                                 &
        ONLY:  e,                                                                                  &
               kh,                                                                                 &
               km,                                                                                 &
               mean_inflow_profiles,                                                               &
               p,                                                                                  &
               pt,                                                                                 &
               pt_init,                                                                            &
               q,                                                                                  &
               q_init,                                                                             &
               ql,                                                                                 &
               ref_state,                                                                          &
               s,                                                                                  &
               s_init,                                                                             &
               u,                                                                                  &
               u_init,                                                                             &
               ug,                                                                                 &
               v,                                                                                  &
               v_init,                                                                             &
               vg,                                                                                 &
               vpt,                                                                                &
               w

    USE averaging

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE gust_mod,                                                                                  &
        ONLY:  gust_module_enabled

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               ny,                                                                                 &
               nys,                                                                                &
               nyn,                                                                                &
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
        ONLY:  module_interface_wrd_global,                                                        &
               module_interface_wrd_local,                                                         &
               module_interface_wrd_local_spinup

    USE netcdf_interface,                                                                          &
        ONLY:  netcdf_precision,                                                                   &
               output_for_t0

    USE particle_attributes,                                                                       &
        ONLY:  first_call_lpm,                                                                     &
               particle_advection,                                                                 &
               time_dopts

    USE pegrid,                                                                                    &
        ONLY:  collective_wait,                                                                    &
               hor_index_bounds,                                                                   &
               myid,                                                                               &
               numprocs

    USE radiation_model_mod,                                                                       &
        ONLY:  time_radiation

    USE random_function_mod,                                                                       &
        ONLY:  random_iv,                                                                          &
               random_iy

    USE random_generator_parallel,                                                                 &
        ONLY:  id_random_array,                                                                    &
               seq_random_array

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  wrd_mpi_io,                                                                         &
               wrd_mpi_io_global_array

    USE spectra_mod,                                                                               &
        ONLY:  average_count_sp,                                                                   &
               spectrum_x,                                                                         &
               spectrum_y

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               hom_sum,                                                                            &
               statistic_regions,                                                                  &
               u_max,                                                                              &
               u_max_ijk,                                                                          &
               v_max,                                                                              &
               v_max_ijk,                                                                          &
               w_max,                                                                              &
               w_max_ijk,                                                                          &
               z_i

    USE surface_data_output_mod,                                                                   &
        ONLY:  surface_data_output_wrd_global,                                                     &
               surface_data_output_wrd_local

    USE surface_mod,                                                                               &
        ONLY:  surface_wrd_local

    USE user,                                                                                      &
        ONLY:  user_module_enabled

    USE virtual_measurement_mod,                                                                   &
        ONLY:  time_virtual_measurement_pr,                                                        &
               time_virtual_measurement_ts,                                                        &
               time_virtual_measurement_tr


    IMPLICIT NONE


    INTERFACE wrd_global
       MODULE PROCEDURE wrd_global
    END INTERFACE wrd_global

    INTERFACE wrd_global_spinup
       MODULE PROCEDURE wrd_global_spinup
    END INTERFACE wrd_global_spinup

    INTERFACE wrd_local
       MODULE PROCEDURE wrd_local
    END INTERFACE wrd_local

    INTERFACE wrd_local_spinup
       MODULE PROCEDURE wrd_local_spinup
    END INTERFACE wrd_local_spinup


    PUBLIC wrd_global,                                                                             &
           wrd_global_spinup,                                                                      &
           wrd_local,                                                                              &
           wrd_local_spinup


 CONTAINS


! Description:
! ------------
!> Global data of control variables and arrays is written out for restarts.
!> In case of output in Fortran binary format, this information is only written to the file
!> opened by PE0.
!--------------------------------------------------------------------------------------------------!


 SUBROUTINE wrd_global


    CHARACTER(LEN=10) ::  binary_version_global  !<
    CHARACTER(LEN=20) ::  tmp_name               !< temporary variable

    INTEGER ::  i                                !< loop index

    binary_version_global = '22.04'

    IF ( restart_data_format_output == 'fortran_binary' )  THEN
!
!--    Output in Fortran binary format
       CALL wrd_write_string( 'binary_version_global' )
       WRITE ( 14 )  binary_version_global

       CALL wrd_write_string( 'numprocs' )
       WRITE ( 14 )  numprocs

       CALL wrd_write_string( 'hor_index_bounds' )
       WRITE ( 14 )  hor_index_bounds

       CALL wrd_write_string( 'nz' )
       WRITE ( 14 )  nz

       CALL wrd_write_string( 'max_pr_user' )
       WRITE ( 14 )  max_pr_user

       CALL wrd_write_string( 'statistic_regions' )
       WRITE ( 14 )  statistic_regions

!
!--    Caution: After changes in the following parameter-list, the version number stored in the
!--    -------- variable binary_version_global has to be increased. The same changes must also be
!--             done in the parameter-list in rrd_global.
       CALL wrd_write_string( 'advected_distance_x' )
       WRITE ( 14 )  advected_distance_x

       CALL wrd_write_string( 'advected_distance_y' )
       WRITE ( 14 )  advected_distance_y

       CALL wrd_write_string( 'allow_negative_scalar_values' )
       WRITE ( 14 )  allow_negative_scalar_values

       CALL wrd_write_string( 'alpha_surface' )
       WRITE ( 14 )  alpha_surface

       CALL wrd_write_string( 'average_count_pr' )
       WRITE ( 14 )  average_count_pr

       CALL wrd_write_string( 'average_count_sp' )
       WRITE ( 14 )  average_count_sp

       CALL wrd_write_string( 'average_count_3d' )
       WRITE ( 14 )  average_count_3d

       CALL wrd_write_string( 'bc_e_b' )
       WRITE ( 14 )  bc_e_b

       CALL wrd_write_string( 'bc_lr' )
       WRITE ( 14 )  bc_lr

       CALL wrd_write_string( 'bc_ns' )
       WRITE ( 14 )  bc_ns

       CALL wrd_write_string( 'bc_p_b' )
       WRITE ( 14 )  bc_p_b

       CALL wrd_write_string( 'bc_p_t' )
       WRITE ( 14 )  bc_p_t

       CALL wrd_write_string( 'bc_pt_b' )
       WRITE ( 14 )  bc_pt_b

       CALL wrd_write_string( 'bc_pt_t' )
       WRITE ( 14 )  bc_pt_t

       CALL wrd_write_string( 'bc_pt_t_val' )
       WRITE ( 14 )  bc_pt_t_val

       CALL wrd_write_string( 'bc_q_b' )
       WRITE ( 14 )  bc_q_b

       CALL wrd_write_string( 'bc_q_t' )
       WRITE ( 14 )  bc_q_t

       CALL wrd_write_string( 'bc_q_t_val' )
       WRITE ( 14 )  bc_q_t_val

       CALL wrd_write_string( 'bc_s_b' )
       WRITE ( 14 )  bc_s_b

       CALL wrd_write_string( 'bc_s_t' )
       WRITE ( 14 )  bc_s_t

       CALL wrd_write_string( 'bc_uv_b' )
       WRITE ( 14 )  bc_uv_b

       CALL wrd_write_string( 'bc_uv_t' )
       WRITE ( 14 )  bc_uv_t

       CALL wrd_write_string( 'building_height' )
       WRITE ( 14 )  building_height

       CALL wrd_write_string( 'building_length_x' )
       WRITE ( 14 )  building_length_x

       CALL wrd_write_string( 'building_length_y' )
       WRITE ( 14 )  building_length_y

       CALL wrd_write_string( 'building_wall_left' )
       WRITE ( 14 )  building_wall_left

       CALL wrd_write_string( 'building_wall_south' )
       WRITE ( 14 )  building_wall_south

       CALL wrd_write_string( 'call_psolver_at_all_substeps' )
       WRITE ( 14 )  call_psolver_at_all_substeps

       CALL wrd_write_string( 'canyon_height' )
       WRITE ( 14 )  canyon_height

       CALL wrd_write_string( 'canyon_wall_left' )
       WRITE ( 14 )  canyon_wall_left

       CALL wrd_write_string( 'canyon_wall_south' )
       WRITE ( 14 )  canyon_wall_south

       CALL wrd_write_string( 'canyon_width_x' )
       WRITE ( 14 )  canyon_width_x

       CALL wrd_write_string( 'canyon_width_y' )
       WRITE ( 14 )  canyon_width_y

       CALL wrd_write_string( 'cfl_factor' )
       WRITE ( 14 )  cfl_factor

       CALL wrd_write_string( 'cloud_droplets' )
       WRITE ( 14 )  cloud_droplets

       CALL wrd_write_string( 'collective_wait' )
       WRITE ( 14 )  collective_wait

       CALL wrd_write_string( 'conserve_volume_flow' )
       WRITE ( 14 )  conserve_volume_flow

       CALL wrd_write_string( 'conserve_volume_flow_mode' )
       WRITE ( 14 )  conserve_volume_flow_mode

       CALL wrd_write_string( 'constant_flux_layer' )
       WRITE ( 14 )  constant_flux_layer

       CALL wrd_write_string( 'coupling_start_time' )
       WRITE ( 14 )  coupling_start_time

       CALL wrd_write_string( 'current_timestep_number' )
       WRITE ( 14 )  current_timestep_number

       CALL wrd_write_string( 'cycle_mg' )
       WRITE ( 14 )  cycle_mg

       CALL wrd_write_string( 'damp_level_1d' )
       WRITE ( 14 )  damp_level_1d

       CALL wrd_write_string( 'origin_date_time' )
       WRITE ( 14 )  origin_date_time

       CALL wrd_write_string( 'dissipation_1d' )
       WRITE ( 14 )  dissipation_1d

       CALL wrd_write_string( 'dp_external' )
       WRITE ( 14 )  dp_external

       CALL wrd_write_string( 'dp_level_b' )
       WRITE ( 14 )  dp_level_b

       CALL wrd_write_string( 'dp_smooth' )
       WRITE ( 14 )  dp_smooth

       CALL wrd_write_string( 'dpdxy' )
       WRITE ( 14 )  dpdxy

       CALL wrd_write_string( 'dt_3d' )
       WRITE ( 14 )  dt_3d

    CALL wrd_write_string( 'dt_pr_1d' )
    WRITE ( 14 )  dt_pr_1d

       CALL wrd_write_string( 'dt_run_control_1d' )
       WRITE ( 14 )  dt_run_control_1d

       CALL wrd_write_string( 'dx' )
       WRITE ( 14 )  dx

       CALL wrd_write_string( 'dy' )
       WRITE ( 14 )  dy

       CALL wrd_write_string( 'dz' )
       WRITE ( 14 )  dz

       CALL wrd_write_string( 'dz_max' )
       WRITE ( 14 )  dz_max

       CALL wrd_write_string( 'dz_stretch_factor' )
       WRITE ( 14 )  dz_stretch_factor

       CALL wrd_write_string( 'dz_stretch_factor_array' )
       WRITE ( 14 )  dz_stretch_factor_array

       CALL wrd_write_string( 'dz_stretch_level' )
       WRITE ( 14 )  dz_stretch_level

       CALL wrd_write_string( 'dz_stretch_level_end' )
       WRITE ( 14 )  dz_stretch_level_end

       CALL wrd_write_string( 'dz_stretch_level_start' )
       WRITE ( 14 )  dz_stretch_level_start

       CALL wrd_write_string( 'e_min' )
       WRITE ( 14 )  e_min

       CALL wrd_write_string( 'end_time_1d' )
       WRITE ( 14 )  end_time_1d

       CALL wrd_write_string( 'fft_method' )
       WRITE ( 14 )  fft_method

       CALL wrd_write_string( 'first_call_lpm' )
       WRITE ( 14 )  first_call_lpm

       CALL wrd_write_string( 'galilei_transformation' )
       WRITE ( 14 )  galilei_transformation

       CALL wrd_write_string( 'hom' )
       WRITE ( 14 )  hom

       CALL wrd_write_string( 'hom_sum' )
       WRITE ( 14 )  hom_sum

       CALL wrd_write_string( 'homogenize_surface_temperature' )
       WRITE ( 14 )  homogenize_surface_temperature

       CALL wrd_write_string( 'humidity' )
       WRITE ( 14 )  humidity

       CALL wrd_write_string( 'inflow_disturbance_begin' )
       WRITE ( 14 )  inflow_disturbance_begin

       CALL wrd_write_string( 'inflow_disturbance_end' )
       WRITE ( 14 )  inflow_disturbance_end

       CALL wrd_write_string( 'km_constant' )
       WRITE ( 14 )  km_constant

       CALL wrd_write_string( 'large_scale_forcing' )
       WRITE ( 14 )  large_scale_forcing

       CALL wrd_write_string( 'large_scale_subsidence' )
       WRITE ( 14 )  large_scale_subsidence

       CALL wrd_write_string( 'latitude' )
       WRITE ( 14 )  latitude

       CALL wrd_write_string( 'longitude' )
       WRITE ( 14 )  longitude

       CALL wrd_write_string( 'loop_optimization' )
       WRITE ( 14 )  loop_optimization

       CALL wrd_write_string( 'masking_method' )
       WRITE ( 14 )  masking_method

       IF ( ALLOCATED( mean_inflow_profiles ) )  THEN
          CALL wrd_write_string( 'mean_inflow_profiles' )
          WRITE ( 14 )  mean_inflow_profiles
       ENDIF

       CALL wrd_write_string( 'mg_cycles' )
       WRITE ( 14 )  mg_cycles

       CALL wrd_write_string( 'mg_switch_to_pe0_level' )
       WRITE ( 14 )  mg_switch_to_pe0_level

       CALL wrd_write_string( 'mixing_length_1d' )
       WRITE ( 14 )  mixing_length_1d

       CALL wrd_write_string( 'momentum_advec' )
       WRITE ( 14 )  momentum_advec

       CALL wrd_write_string( 'netcdf_precision' )
       WRITE ( 14 )  netcdf_precision

       CALL wrd_write_string( 'neutral' )
       WRITE ( 14 )  neutral

       CALL wrd_write_string( 'ngsrb' )
       WRITE ( 14 )  ngsrb

       CALL wrd_write_string( 'nsor' )
       WRITE ( 14 )  nsor

       CALL wrd_write_string( 'nsor_ini' )
       WRITE ( 14 )  nsor_ini

       CALL wrd_write_string( 'nudging' )
       WRITE ( 14 )  nudging

       CALL wrd_write_string( 'num_leg' )
       WRITE ( 14 )  num_leg

       CALL wrd_write_string( 'nx' )
       WRITE ( 14 )  nx

       CALL wrd_write_string( 'ny' )
       WRITE ( 14 )  ny

       CALL wrd_write_string( 'ocean_mode' )
       WRITE ( 14 )  ocean_mode

       CALL wrd_write_string( 'omega' )
       WRITE ( 14 )  omega

       CALL wrd_write_string( 'omega_sor' )
       WRITE ( 14 )  omega_sor

       CALL wrd_write_string( 'outflow_damping_factor' )
       WRITE ( 14 )  outflow_damping_factor

       CALL wrd_write_string( 'outflow_damping_width' )
       WRITE ( 14 )  outflow_damping_width

       CALL wrd_write_string( 'output_for_t0' )
       WRITE ( 14 )  output_for_t0

       CALL wrd_write_string( 'passive_scalar' )
       WRITE ( 14 )  passive_scalar

       CALL wrd_write_string( 'prandtl_number' )
       WRITE ( 14 )  prandtl_number

       CALL wrd_write_string( 'psolver' )
       WRITE ( 14 )  psolver

       CALL wrd_write_string( 'pt_damping_factor' )
       WRITE ( 14 )  pt_damping_factor

       CALL wrd_write_string( 'pt_damping_width' )
       WRITE ( 14 )  pt_damping_width

       CALL wrd_write_string( 'pt_init' )
       WRITE ( 14 )  pt_init

       CALL wrd_write_string( 'pt_reference' )
       WRITE ( 14 )  pt_reference

       CALL wrd_write_string( 'pt_surface' )
       WRITE ( 14 )  pt_surface

       CALL wrd_write_string( 'pt_surface_heating_rate' )
       WRITE ( 14 )  pt_surface_heating_rate

       CALL wrd_write_string( 'pt_vertical_gradient' )
       WRITE ( 14 )  pt_vertical_gradient

       CALL wrd_write_string( 'pt_vertical_gradient_level' )
       WRITE ( 14 )  pt_vertical_gradient_level

       CALL wrd_write_string( 'pt_vertical_gradient_level_ind' )
       WRITE ( 14 )  pt_vertical_gradient_level_ind

       CALL wrd_write_string( 'q_init' )
       WRITE ( 14 )  q_init

       CALL wrd_write_string( 'q_surface' )
       WRITE ( 14 )  q_surface

       CALL wrd_write_string( 'q_vertical_gradient' )
       WRITE ( 14 )  q_vertical_gradient

       CALL wrd_write_string( 'q_vertical_gradient_level' )
       WRITE ( 14 )  q_vertical_gradient_level

       CALL wrd_write_string( 'q_vertical_gradient_level_ind' )
       WRITE ( 14 )  q_vertical_gradient_level_ind

       CALL wrd_write_string( 'random_generator' )
       WRITE ( 14 )  random_generator

       CALL wrd_write_string( 'random_heatflux' )
       WRITE ( 14 )  random_heatflux

       CALL wrd_write_string( 'rans_mode' )
       WRITE ( 14 )  rans_mode

       CALL wrd_write_string( 'rayleigh_damping_factor' )
       WRITE ( 14 )  rayleigh_damping_factor

       CALL wrd_write_string( 'rayleigh_damping_height' )
       WRITE ( 14 )  rayleigh_damping_height

       CALL wrd_write_string( 'ref_state' )
       WRITE ( 14 )  ref_state

       CALL wrd_write_string( 'reference_state' )
       WRITE ( 14 )  reference_state

       CALL wrd_write_string( 'residual_limit' )
       WRITE ( 14 )  residual_limit

       CALL wrd_write_string( 'roughness_length' )
       WRITE ( 14 )  roughness_length

       CALL wrd_write_string( 'runnr' )
       WRITE ( 14 )  runnr

       CALL wrd_write_string( 's_init' )
       WRITE ( 14 )  s_init

       CALL wrd_write_string( 's_surface' )
       WRITE ( 14 )  s_surface

       CALL wrd_write_string( 's_vertical_gradient' )
       WRITE ( 14 )  s_vertical_gradient

       CALL wrd_write_string( 's_vertical_gradient_level' )
       WRITE ( 14 )  s_vertical_gradient_level

       CALL wrd_write_string( 's_vertical_gradient_level_ind' )
       WRITE ( 14 )  s_vertical_gradient_level_ind

       CALL wrd_write_string( 'scalar_advec' )
       WRITE ( 14 )  scalar_advec

       CALL wrd_write_string( 'simulated_time' )
       WRITE ( 14 )  simulated_time

       IF ( ALLOCATED( spectrum_x ) )  THEN
          CALL wrd_write_string( 'spectrum_x' )
          WRITE ( 14 )  spectrum_x
          CALL wrd_write_string( 'spectrum_y' )
          WRITE ( 14 )  spectrum_y
       ENDIF

       CALL wrd_write_string( 'spinup_time ' )
       WRITE ( 14 )  spinup_time

       CALL wrd_write_string( 'subs_vertical_gradient' )
       WRITE ( 14 )  subs_vertical_gradient

       CALL wrd_write_string( 'subs_vertical_gradient_level' )
       WRITE ( 14 )  subs_vertical_gradient_level

       CALL wrd_write_string( 'subs_vertical_gradient_level_i' )
       WRITE ( 14 )  subs_vertical_gradient_level_i

       CALL wrd_write_string( 'surface_heatflux' )
       WRITE ( 14 )  surface_heatflux

       CALL wrd_write_string( 'surface_pressure' )
       WRITE ( 14 )  surface_pressure

       CALL wrd_write_string( 'surface_scalarflux' )
       WRITE ( 14 )  surface_scalarflux

       CALL wrd_write_string( 'surface_waterflux' )
       WRITE ( 14 )  surface_waterflux

       CALL wrd_write_string( 'time_coupling' )
       WRITE ( 14 )  time_coupling

       CALL wrd_write_string( 'time_disturb' )
       WRITE ( 14 )  time_disturb

       CALL wrd_write_string( 'time_do2d_xy' )
       WRITE ( 14 )  time_do2d_xy

       CALL wrd_write_string( 'time_do2d_xz' )
       WRITE ( 14 )  time_do2d_xz

       CALL wrd_write_string( 'time_do2d_yz' )
       WRITE ( 14 )  time_do2d_yz

       CALL wrd_write_string( 'time_do3d' )
       WRITE ( 14 )  time_do3d

       CALL wrd_write_string( 'time_do_av' )
       WRITE ( 14 )  time_do_av

       CALL wrd_write_string( 'time_do_sla' )
       WRITE ( 14 )  time_do_sla

       CALL wrd_write_string( 'time_domask' )
       WRITE ( 14 )  time_domask

       CALL wrd_write_string( 'time_dopr' )
       WRITE ( 14 )  time_dopr

       CALL wrd_write_string( 'time_dopr_av' )
       WRITE ( 14 )  time_dopr_av

       CALL wrd_write_string( 'time_dopr_listing' )
       WRITE ( 14 )  time_dopr_listing

       CALL wrd_write_string( 'time_dopts' )
       WRITE ( 14 )  time_dopts

       CALL wrd_write_string( 'time_dosp' )
       WRITE ( 14 )  time_dosp

       CALL wrd_write_string( 'time_dots' )
       WRITE ( 14 )  time_dots

       CALL wrd_write_string( 'time_indoor' )
       WRITE ( 14 )  time_indoor

       CALL wrd_write_string( 'time_radiation' )
       WRITE ( 14 )  time_radiation

       CALL wrd_write_string( 'time_restart' )
       WRITE ( 14 )  time_restart

       CALL wrd_write_string( 'time_run_control' )
       WRITE ( 14 )  time_run_control

       CALL wrd_write_string( 'time_since_reference_point' )
       WRITE ( 14 )  time_since_reference_point

       CALL wrd_write_string( 'time_virtual_measurement_pr' )
       WRITE ( 14 )  time_virtual_measurement_pr

       CALL wrd_write_string( 'time_virtual_measurement_ts' )
       WRITE ( 14 )  time_virtual_measurement_ts

       CALL wrd_write_string( 'time_virtual_measurement_tr' )
       WRITE ( 14 )  time_virtual_measurement_tr

       CALL wrd_write_string( 'timestep_scheme' )
       WRITE ( 14 )  timestep_scheme

       CALL wrd_write_string( 'top_heatflux' )
       WRITE ( 14 )  top_heatflux

       CALL wrd_write_string( 'top_momentumflux_u' )
       WRITE ( 14 )  top_momentumflux_u

       CALL wrd_write_string( 'top_momentumflux_v' )
       WRITE ( 14 )  top_momentumflux_v

       CALL wrd_write_string( 'top_scalarflux' )
       WRITE ( 14 )  top_scalarflux

       CALL wrd_write_string( 'topography' )
       WRITE ( 14 )  topography

       CALL wrd_write_string( 'topography_grid_convention' )
       WRITE ( 14 )  topography_grid_convention

       CALL wrd_write_string( 'tsc' )
       WRITE ( 14 )  tsc

       CALL wrd_write_string( 'tunnel_height' )
       WRITE ( 14 )  tunnel_height

       CALL wrd_write_string( 'tunnel_length' )
       WRITE ( 14 )  tunnel_length

       CALL wrd_write_string( 'tunnel_wall_depth' )
       WRITE ( 14 )  tunnel_wall_depth

       CALL wrd_write_string( 'tunnel_width_x' )
       WRITE ( 14 )  tunnel_width_x

       CALL wrd_write_string( 'tunnel_width_y' )
       WRITE ( 14 )  tunnel_width_y

       CALL wrd_write_string( 'turbulence_closure' )
       WRITE ( 14 )  turbulence_closure

       CALL wrd_write_string( 'u_bulk' )
       WRITE ( 14 )  u_bulk

       CALL wrd_write_string( 'u_init' )
       WRITE ( 14 )  u_init

       CALL wrd_write_string( 'u_max' )
       WRITE ( 14 )  u_max

       CALL wrd_write_string( 'u_max_ijk' )
       WRITE ( 14 )  u_max_ijk

       CALL wrd_write_string( 'ug' )
       WRITE ( 14 )  ug

       CALL wrd_write_string( 'ug_surface' )
       WRITE ( 14 )  ug_surface

       CALL wrd_write_string( 'ug_vertical_gradient' )
       WRITE ( 14 )  ug_vertical_gradient

       CALL wrd_write_string( 'ug_vertical_gradient_level' )
       WRITE ( 14 )  ug_vertical_gradient_level

       CALL wrd_write_string( 'ug_vertical_gradient_level_ind' )
       WRITE ( 14 )  ug_vertical_gradient_level_ind

       CALL wrd_write_string( 'use_surface_fluxes' )
       WRITE ( 14 )  use_surface_fluxes

       CALL wrd_write_string( 'use_top_fluxes' )
       WRITE ( 14 )  use_top_fluxes

       CALL wrd_write_string( 'use_ug_for_galilei_tr' )
       WRITE ( 14 )  use_ug_for_galilei_tr

       CALL wrd_write_string( 'use_upstream_for_tke' )
       WRITE ( 14 )  use_upstream_for_tke

       CALL wrd_write_string( 'v_bulk' )
       WRITE ( 14 )  v_bulk

       CALL wrd_write_string( 'v_init' )
       WRITE ( 14 )  v_init

       CALL wrd_write_string( 'v_max' )
       WRITE ( 14 )  v_max

       CALL wrd_write_string( 'v_max_ijk' )
       WRITE ( 14 )  v_max_ijk

       CALL wrd_write_string( 'vg' )
       WRITE ( 14 )  vg

       CALL wrd_write_string( 'vg_surface' )
       WRITE ( 14 )  vg_surface

       CALL wrd_write_string( 'vg_vertical_gradient' )
       WRITE ( 14 ) vg_vertical_gradient

       CALL wrd_write_string( 'vg_vertical_gradient_level' )
       WRITE ( 14 )  vg_vertical_gradient_level

       CALL wrd_write_string( 'vg_vertical_gradient_level_ind' )
       WRITE ( 14 )  vg_vertical_gradient_level_ind

       CALL wrd_write_string( 'virtual_flight' )
       WRITE ( 14 )  virtual_flight

       CALL wrd_write_string( 'volume_flow_area' )
       WRITE ( 14 )  volume_flow_area

       CALL wrd_write_string( 'volume_flow_initial' )
       WRITE ( 14 )  volume_flow_initial

       CALL wrd_write_string( 'w_max' )
       WRITE ( 14 )  w_max

       CALL wrd_write_string( 'w_max_ijk' )
       WRITE ( 14 )  w_max_ijk

       CALL wrd_write_string( 'wall_adjustment' )
       WRITE ( 14 )  wall_adjustment

       CALL wrd_write_string( 'wall_heatflux' )
       WRITE ( 14 )  wall_heatflux

       CALL wrd_write_string( 'wall_humidityflux' )
       WRITE ( 14 )  wall_humidityflux

       CALL wrd_write_string( 'wall_scalarflux' )
       WRITE ( 14 )  wall_scalarflux

       CALL wrd_write_string( 'y_shift' )
       WRITE ( 14 )  y_shift

       CALL wrd_write_string( 'z0h_factor' )
       WRITE ( 14 )  z0h_factor

       CALL wrd_write_string( 'zeta_max' )
       WRITE ( 14 )  zeta_max

       CALL wrd_write_string( 'zeta_min' )
       WRITE ( 14 )  zeta_min

       CALL wrd_write_string( 'z_i' )
       WRITE ( 14 )  z_i

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--    Write global restart data using MPI-IO
!--    ATTENTION: Arrays need to be written with routine wrd_mpi_io_global_array!
       CALL wrd_mpi_io( 'binary_version_global',  binary_version_global )
       CALL wrd_mpi_io( 'numprocs',  numprocs )
       CALL wrd_mpi_io( 'nz' , nz )
       CALL wrd_mpi_io( 'max_pr_user',  max_pr_user )
       CALL wrd_mpi_io( 'statistic_regions', statistic_regions )

!
!--    Caution: After changes in the following parameter-list, the version number stored in the
!--    -------- variable binary_version_global has to be increased. The same changes must also be
!--             done in the parameter-list in rrd_global.
       CALL wrd_mpi_io( 'advected_distance_x',  advected_distance_x )
       CALL wrd_mpi_io( 'advected_distance_y', advected_distance_y )
       CALL wrd_mpi_io( 'allow_negative_scalar_values', allow_negative_scalar_values )
       CALL wrd_mpi_io( 'alpha_surface', alpha_surface )
       CALL wrd_mpi_io( 'average_count_pr', average_count_pr )
       CALL wrd_mpi_io( 'average_count_sp', average_count_sp )
       CALL wrd_mpi_io( 'average_count_3d', average_count_3d )
       CALL wrd_mpi_io( 'bc_e_b', bc_e_b )
       CALL wrd_mpi_io( 'bc_lr', bc_lr )
       CALL wrd_mpi_io( 'bc_ns', bc_ns )
       CALL wrd_mpi_io( 'bc_p_b', bc_p_b )
       CALL wrd_mpi_io( 'bc_p_t', bc_p_t )
       CALL wrd_mpi_io( 'bc_pt_b', bc_pt_b )
       CALL wrd_mpi_io( 'bc_pt_t', bc_pt_t )
       CALL wrd_mpi_io( 'bc_pt_t_val', bc_pt_t_val )
       CALL wrd_mpi_io( 'bc_q_b', bc_q_b )
       CALL wrd_mpi_io( 'bc_q_t', bc_q_t )
       CALL wrd_mpi_io( 'bc_q_t_val', bc_q_t_val )
       CALL wrd_mpi_io( 'bc_s_b', bc_s_b )
       CALL wrd_mpi_io( 'bc_s_t', bc_s_t )
       CALL wrd_mpi_io( 'bc_uv_b', bc_uv_b )
       CALL wrd_mpi_io( 'bc_uv_t', bc_uv_t )
       CALL wrd_mpi_io( 'biometeorology', biometeorology )
       CALL wrd_mpi_io( 'building_height', building_height )
       CALL wrd_mpi_io( 'building_length_x', building_length_x )
       CALL wrd_mpi_io( 'building_length_y', building_length_y )
       CALL wrd_mpi_io( 'building_wall_left', building_wall_left )
       CALL wrd_mpi_io( 'building_wall_south', building_wall_south )
       CALL wrd_mpi_io( 'bulk_cloud_model', bulk_cloud_model )
       CALL wrd_mpi_io( 'call_psolver_at_all_substeps', call_psolver_at_all_substeps )
       CALL wrd_mpi_io( 'canyon_height', canyon_height )
       CALL wrd_mpi_io( 'canyon_wall_left', canyon_wall_left )
       CALL wrd_mpi_io( 'canyon_wall_south', canyon_wall_south )
       CALL wrd_mpi_io( 'canyon_width_x',  canyon_width_x )
       CALL wrd_mpi_io( 'canyon_width_y', canyon_width_y )
       CALL wrd_mpi_io( 'cfl_factor', cfl_factor )
       CALL wrd_mpi_io( 'cloud_droplets',  cloud_droplets )
       CALL wrd_mpi_io( 'collective_wait', collective_wait )
       CALL wrd_mpi_io( 'conserve_volume_flow', conserve_volume_flow )
       CALL wrd_mpi_io( 'conserve_volume_flow_mode', conserve_volume_flow_mode )
       CALL wrd_mpi_io( 'constant_flux_layer', constant_flux_layer )
       CALL wrd_mpi_io( 'coupling_start_time', coupling_start_time )
       CALL wrd_mpi_io( 'current_timestep_number', current_timestep_number )
       CALL wrd_mpi_io( 'cycle_mg', cycle_mg )
       CALL wrd_mpi_io( 'damp_level_1d', damp_level_1d )
       CALL wrd_mpi_io( 'dissipation_1d', dissipation_1d )
       CALL wrd_mpi_io( 'dp_external', dp_external )
       CALL wrd_mpi_io( 'dp_level_b', dp_level_b )
       CALL wrd_mpi_io( 'dp_smooth', dp_smooth )
       CALL wrd_mpi_io_global_array( 'dpdxy', dpdxy )
       CALL wrd_mpi_io( 'dt_3d', dt_3d )
       CALL wrd_mpi_io( 'dt_pr_1d', dt_pr_1d )
       CALL wrd_mpi_io( 'dt_run_control_1d', dt_run_control_1d )
       CALL wrd_mpi_io( 'dx', dx )
       CALL wrd_mpi_io( 'dy', dy )
       CALL wrd_mpi_io_global_array( 'dz', dz )
       CALL wrd_mpi_io( 'dz_max', dz_max )
       CALL wrd_mpi_io( 'dz_stretch_factor', dz_stretch_factor )
       CALL wrd_mpi_io_global_array( 'dz_stretch_factor_array', dz_stretch_factor_array )
       CALL wrd_mpi_io( 'dz_stretch_level', dz_stretch_level )
       CALL wrd_mpi_io_global_array( 'dz_stretch_level_end', dz_stretch_level_end )
       CALL wrd_mpi_io_global_array( 'dz_stretch_level_start', dz_stretch_level_start )
       CALL wrd_mpi_io( 'e_min', e_min )
       CALL wrd_mpi_io( 'end_time_1d', end_time_1d )
       CALL wrd_mpi_io( 'fft_method', fft_method )
       CALL wrd_mpi_io( 'first_call_lpm', first_call_lpm )
       CALL wrd_mpi_io( 'galilei_transformation', galilei_transformation )
       CALL wrd_mpi_io( 'gust_module_enabled', gust_module_enabled )
       CALL wrd_mpi_io_global_array( 'hom', hom )
       CALL wrd_mpi_io_global_array( 'hom_sum', hom_sum )
       CALL wrd_mpi_io( 'homogenize_surface_temperature', homogenize_surface_temperature )
       CALL wrd_mpi_io( 'humidity', humidity )
       CALL wrd_mpi_io( 'inflow_disturbance_begin', inflow_disturbance_begin )
       CALL wrd_mpi_io( 'inflow_disturbance_end', inflow_disturbance_end )
       CALL wrd_mpi_io( 'km_constant', km_constant )
       CALL wrd_mpi_io( 'large_scale_forcing', large_scale_forcing )
       CALL wrd_mpi_io( 'large_scale_subsidence', large_scale_subsidence )
       CALL wrd_mpi_io( 'latitude', latitude )
       CALL wrd_mpi_io( 'longitude', longitude )
       CALL wrd_mpi_io( 'loop_optimization', loop_optimization )
       CALL wrd_mpi_io( 'masking_method', masking_method )
       IF ( ALLOCATED( mean_inflow_profiles ) )  THEN
          CALL wrd_mpi_io_global_array( 'mean_inflow_profiles', mean_inflow_profiles )
       ENDIF
       CALL wrd_mpi_io( 'mg_cycles', mg_cycles )
       CALL wrd_mpi_io( 'mg_switch_to_pe0_level', mg_switch_to_pe0_level )
       CALL wrd_mpi_io( 'mixing_length_1d', mixing_length_1d )
       CALL wrd_mpi_io( 'momentum_advec', momentum_advec )
!
!--    There is no module procedure for CHARACTER arrays
       DO  i = 1, SIZE( netcdf_precision, 1 )
          WRITE ( tmp_name, '(A,I2.2)' )  'netcdf_precision_', i
          CALL wrd_mpi_io( TRIM( tmp_name ), netcdf_precision(i) )
       ENDDO
       CALL wrd_mpi_io( 'neutral', neutral )
       CALL wrd_mpi_io( 'ngsrb', ngsrb )
       CALL wrd_mpi_io( 'nsor', nsor )
       CALL wrd_mpi_io( 'nsor_ini', nsor_ini )
       CALL wrd_mpi_io( 'nudging', nudging )
       CALL wrd_mpi_io( 'num_leg', num_leg )
       CALL wrd_mpi_io( 'nx', nx )
       CALL wrd_mpi_io( 'ny', ny )
       CALL wrd_mpi_io( 'ocean_mode', ocean_mode )
       CALL wrd_mpi_io( 'omega', omega )
       CALL wrd_mpi_io( 'omega_sor', omega_sor )
       CALL wrd_mpi_io( 'origin_date_time', origin_date_time )
       CALL wrd_mpi_io( 'outflow_damping_factor', outflow_damping_factor )
       CALL wrd_mpi_io( 'outflow_damping_width', outflow_damping_width )
       CALL wrd_mpi_io( 'output_for_t0', output_for_t0 )
       CALL wrd_mpi_io( 'particle_advection', particle_advection )
       CALL wrd_mpi_io( 'passive_scalar', passive_scalar )
       CALL wrd_mpi_io( 'prandtl_number', prandtl_number )
       CALL wrd_mpi_io( 'psolver', psolver )
       CALL wrd_mpi_io( 'pt_damping_factor', pt_damping_factor )
       CALL wrd_mpi_io( 'pt_damping_width', pt_damping_width )
       CALL wrd_mpi_io_global_array( 'pt_init', pt_init )
       CALL wrd_mpi_io( 'pt_reference', pt_reference )
       CALL wrd_mpi_io( 'pt_surface', pt_surface )
       CALL wrd_mpi_io( 'pt_surface_heating_rate', pt_surface_heating_rate )
       CALL wrd_mpi_io_global_array( 'pt_vertical_gradient', pt_vertical_gradient )
       CALL wrd_mpi_io_global_array( 'pt_vertical_gradient_level', pt_vertical_gradient_level )
       CALL wrd_mpi_io_global_array( 'pt_vertical_gradient_level_ind', pt_vertical_gradient_level_ind )
       CALL wrd_mpi_io_global_array( 'q_init', q_init )
       CALL wrd_mpi_io( 'q_surface', q_surface )
       CALL wrd_mpi_io_global_array( 'q_vertical_gradient', q_vertical_gradient )
       CALL wrd_mpi_io_global_array( 'q_vertical_gradient_level', q_vertical_gradient_level )
       CALL wrd_mpi_io_global_array( 'q_vertical_gradient_level_ind', q_vertical_gradient_level_ind )
       CALL wrd_mpi_io( 'random_generator', random_generator )
       CALL wrd_mpi_io( 'random_heatflux', random_heatflux )
       CALL wrd_mpi_io( 'rans_mode', rans_mode )
       CALL wrd_mpi_io( 'rayleigh_damping_factor', rayleigh_damping_factor )
       CALL wrd_mpi_io( 'rayleigh_damping_height', rayleigh_damping_height )
       CALL wrd_mpi_io_global_array( 'ref_state', ref_state )
       CALL wrd_mpi_io( 'reference_state', reference_state )
       CALL wrd_mpi_io( 'residual_limit', residual_limit )
       CALL wrd_mpi_io( 'roughness_length', roughness_length )
       CALL wrd_mpi_io( 'runnr', runnr )
       CALL wrd_mpi_io_global_array( 's_init', s_init )
       CALL wrd_mpi_io( 's_surface', s_surface )
       CALL wrd_mpi_io_global_array( 's_vertical_gradient', s_vertical_gradient )
       CALL wrd_mpi_io_global_array( 's_vertical_gradient_level', s_vertical_gradient_level )
       CALL wrd_mpi_io_global_array( 's_vertical_gradient_level_ind', s_vertical_gradient_level_ind )
       CALL wrd_mpi_io( 'scalar_advec', scalar_advec )
       CALL wrd_mpi_io( 'simulated_time', simulated_time )
       IF ( ALLOCATED( spectrum_x ) )  THEN
          CALL wrd_mpi_io_global_array( 'spectrum_x', spectrum_x )
          CALL wrd_mpi_io_global_array( 'spectrum_y', spectrum_y )
       ENDIF
       CALL wrd_mpi_io( 'spinup_time ', spinup_time )
       CALL wrd_mpi_io_global_array( 'subs_vertical_gradient', subs_vertical_gradient )
       CALL wrd_mpi_io_global_array( 'subs_vertical_gradient_level', subs_vertical_gradient_level )
       CALL wrd_mpi_io_global_array( 'subs_vertical_gradient_level_i', subs_vertical_gradient_level_i )
       CALL wrd_mpi_io( 'surface_heatflux', surface_heatflux )
       CALL wrd_mpi_io( 'surface_output', surface_output )
       CALL wrd_mpi_io( 'surface_pressure', surface_pressure )
       CALL wrd_mpi_io( 'surface_scalarflux', surface_scalarflux )
       CALL wrd_mpi_io( 'surface_waterflux', surface_waterflux )
       CALL wrd_mpi_io( 'time_coupling', time_coupling )
       CALL wrd_mpi_io( 'time_disturb', time_disturb )
       CALL wrd_mpi_io( 'time_do2d_xy', time_do2d_xy )
       CALL wrd_mpi_io( 'time_do2d_xz', time_do2d_xz )
       CALL wrd_mpi_io( 'time_do2d_yz', time_do2d_yz )
       CALL wrd_mpi_io( 'time_do3d', time_do3d )
       CALL wrd_mpi_io( 'time_do_av', time_do_av )
       CALL wrd_mpi_io( 'time_do_sla', time_do_sla )
       CALL wrd_mpi_io_global_array( 'time_domask', time_domask )
       CALL wrd_mpi_io( 'time_dopr', time_dopr )
       CALL wrd_mpi_io( 'time_dopr_av', time_dopr_av )
       CALL wrd_mpi_io( 'time_dopr_listing', time_dopr_listing )
       CALL wrd_mpi_io( 'time_dopts', time_dopts )
       CALL wrd_mpi_io( 'time_dosp', time_dosp )
       CALL wrd_mpi_io( 'time_dots', time_dots )
       CALL wrd_mpi_io( 'time_indoor', time_indoor )
       CALL wrd_mpi_io( 'time_radiation', time_radiation )
       CALL wrd_mpi_io( 'time_restart', time_restart )
       CALL wrd_mpi_io( 'time_run_control', time_run_control )
       CALL wrd_mpi_io( 'time_since_reference_point', time_since_reference_point )
       CALL wrd_mpi_io( 'time_virtual_measurement_pr', time_virtual_measurement_pr )
       CALL wrd_mpi_io( 'time_virtual_measurement_ts', time_virtual_measurement_ts )
       CALL wrd_mpi_io( 'time_virtual_measurement_tr', time_virtual_measurement_tr )
       CALL wrd_mpi_io( 'timestep_scheme', timestep_scheme )
       CALL wrd_mpi_io( 'top_heatflux', top_heatflux )
       CALL wrd_mpi_io( 'top_momentumflux_u', top_momentumflux_u )
       CALL wrd_mpi_io( 'top_momentumflux_v', top_momentumflux_v )
       CALL wrd_mpi_io( 'top_scalarflux', top_scalarflux )
       CALL wrd_mpi_io( 'topography', topography )
       CALL wrd_mpi_io( 'topography_grid_convention', topography_grid_convention )
       CALL wrd_mpi_io_global_array( 'tsc', tsc )
       CALL wrd_mpi_io( 'tunnel_height', tunnel_height )
       CALL wrd_mpi_io( 'tunnel_length', tunnel_length )
       CALL wrd_mpi_io( 'tunnel_wall_depth', tunnel_wall_depth )
       CALL wrd_mpi_io( 'tunnel_width_x', tunnel_width_x )
       CALL wrd_mpi_io( 'tunnel_width_y', tunnel_width_y )
       CALL wrd_mpi_io( 'turbulence_closure', turbulence_closure )
       CALL wrd_mpi_io( 'u_bulk', u_bulk )
       CALL wrd_mpi_io_global_array( 'u_init', u_init )
       CALL wrd_mpi_io( 'u_max', u_max )
       CALL wrd_mpi_io_global_array( 'u_max_ijk', u_max_ijk )
       CALL wrd_mpi_io_global_array( 'ug', ug )
       CALL wrd_mpi_io( 'ug_surface', ug_surface )
       CALL wrd_mpi_io_global_array( 'ug_vertical_gradient', ug_vertical_gradient )
       CALL wrd_mpi_io_global_array( 'ug_vertical_gradient_level', ug_vertical_gradient_level )
       CALL wrd_mpi_io_global_array( 'ug_vertical_gradient_level_ind', ug_vertical_gradient_level_ind )
       CALL wrd_mpi_io( 'use_surface_fluxes', use_surface_fluxes )
       CALL wrd_mpi_io( 'use_top_fluxes', use_top_fluxes )
       CALL wrd_mpi_io( 'use_ug_for_galilei_tr', use_ug_for_galilei_tr )
       CALL wrd_mpi_io( 'use_upstream_for_tke', use_upstream_for_tke )
       CALL wrd_mpi_io( 'user_module_enabled', user_module_enabled )
       CALL wrd_mpi_io( 'v_bulk', v_bulk )
       CALL wrd_mpi_io_global_array( 'v_init', v_init )
       CALL wrd_mpi_io( 'v_max', v_max )
       CALL wrd_mpi_io_global_array( 'v_max_ijk', v_max_ijk )
       CALL wrd_mpi_io_global_array( 'vg', vg )
       CALL wrd_mpi_io( 'vg_surface', vg_surface )
       CALL wrd_mpi_io_global_array( 'vg_vertical_gradient', vg_vertical_gradient )
       CALL wrd_mpi_io_global_array( 'vg_vertical_gradient_level', vg_vertical_gradient_level )
       CALL wrd_mpi_io_global_array( 'vg_vertical_gradient_level_ind', vg_vertical_gradient_level_ind )
       CALL wrd_mpi_io( 'virtual_flight', virtual_flight )
       CALL wrd_mpi_io_global_array( 'volume_flow_area', volume_flow_area )
       CALL wrd_mpi_io_global_array( 'volume_flow_initial', volume_flow_initial )
       CALL wrd_mpi_io( 'w_max', w_max )
       CALL wrd_mpi_io_global_array( 'w_max_ijk', w_max_ijk )
       CALL wrd_mpi_io( 'wall_adjustment', wall_adjustment )
       CALL wrd_mpi_io_global_array( 'wall_heatflux', wall_heatflux )
       CALL wrd_mpi_io_global_array( 'wall_humidityflux', wall_humidityflux )
       CALL wrd_mpi_io_global_array( 'wall_scalarflux', wall_scalarflux )
       CALL wrd_mpi_io( 'y_shift', y_shift )
       CALL wrd_mpi_io( 'z0h_factor', z0h_factor )
       CALL wrd_mpi_io( 'zeta_max', zeta_max )
       CALL wrd_mpi_io( 'zeta_min', zeta_min )
       CALL wrd_mpi_io_global_array( 'z_i', z_i )

    ENDIF
!
!-- Write restart data of surface_data_output_mod
    IF ( surface_output )  CALL surface_data_output_wrd_global
!
!-- Write restart data of the other modules
    CALL module_interface_wrd_global


 END SUBROUTINE wrd_global

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write global parameters to spinup restart file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_global_spinup
!
!-- Write nx and ny. This is actually only required to check if the spinup-surface data matches
!-- the dimensions in the main run, in order to enable the cyclic-fill control flag.
    CALL wrd_mpi_io( 'nx', nx )
    CALL wrd_mpi_io( 'ny', ny )

 END SUBROUTINE wrd_global_spinup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Processor specific data of variables and arrays is written out for restarts (binary format).
!> This information is written to the file opened by each PE.
!--------------------------------------------------------------------------------------------------!


 SUBROUTINE wrd_local


    CHARACTER(LEN=10) ::  binary_version_local  !<
    CHARACTER(LEN=20) ::  tmp_name              !<  temporary variable

    INTEGER ::  i                               !<  loop index

!
!-- Write arrays.
    binary_version_local = '22.04'

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
!
!--    Output in Fortran binary format

       CALL wrd_write_string( 'binary_version_local' )
       WRITE ( 14 )  binary_version_local

       WRITE ( 14 )  numprocs, myid, nxl, nxr, nys, nyn, nzb, nzt

!
!--    Attention: After changes to the following output commands the version number of the variable
!--    ---------- binary_version_local must be changed!
!--               Also, the list of arrays to be read in rrd_local must be adjusted accordingly.
       CALL wrd_write_string( 'e' )
       WRITE ( 14 )  e

       IF ( ALLOCATED( e_av ) )  THEN
          CALL wrd_write_string( 'e_av' )
          WRITE ( 14 )  e_av
       ENDIF

       CALL wrd_write_string( 'kh' )
       WRITE ( 14 )  kh


       IF ( ALLOCATED( kh_av ) )  THEN
          CALL wrd_write_string( 'kh_av' )
          WRITE ( 14 )  kh_av
       ENDIF

       CALL wrd_write_string( 'km' )
       WRITE ( 14 )  km

       IF ( ALLOCATED( km_av ) )  THEN
          CALL wrd_write_string( 'km_av' )
          WRITE ( 14 )  km_av
       ENDIF

       IF ( ALLOCATED( lpt_av ) )  THEN
          CALL wrd_write_string( 'lpt_av' )
          WRITE ( 14 )  lpt_av
       ENDIF

       IF ( ALLOCATED( lwp_av ) )  THEN
          CALL wrd_write_string( 'lwp_av' )
          WRITE ( 14 )  lwp_av
       ENDIF

       IF ( ALLOCATED( ol_av ) )  THEN
          CALL wrd_write_string( 'ol_av' )
          WRITE ( 14 )  ol_av
       ENDIF

       CALL wrd_write_string( 'p' )
       WRITE ( 14 )  p

       IF ( ALLOCATED( p_av ) )  THEN
          CALL wrd_write_string( 'p_av' )
          WRITE ( 14 )  p_av
       ENDIF

       IF ( ALLOCATED( pc_av ) )  THEN
          CALL wrd_write_string( 'pc_av' )
          WRITE ( 14 )  pc_av
       ENDIF

       IF ( ALLOCATED( pr_av ) )  THEN
          CALL wrd_write_string( 'pr_av' )
          WRITE ( 14 )  pr_av
       ENDIF

       IF ( ALLOCATED( pres_drag_x_av ) )  THEN
          CALL wrd_write_string( 'pres_drag_x_av' )
          WRITE ( 14 )  pres_drag_x_av
       ENDIF

       IF ( ALLOCATED( pres_drag_y_av ) )  THEN
          CALL wrd_write_string( 'pres_drag_y_av' )
          WRITE ( 14 )  pres_drag_y_av
       ENDIF

       CALL wrd_write_string( 'pt' )
       WRITE ( 14 )  pt

       IF ( ALLOCATED( pt_av ) )  THEN
          CALL wrd_write_string( 'pt_av' )
          WRITE ( 14 )  pt_av
       ENDIF

       IF ( humidity )  THEN

          CALL wrd_write_string( 'q' )
          WRITE ( 14 )  q

          IF ( ALLOCATED( q_av ) )  THEN
             CALL wrd_write_string( 'q_av' )
             WRITE ( 14 )  q_av
          ENDIF

          IF ( cloud_droplets )  THEN

             CALL wrd_write_string( 'ql' )
             WRITE ( 14 )  ql

             IF ( ALLOCATED( ql_av ) )  THEN
                CALL wrd_write_string( 'ql_av' )
                WRITE ( 14 )  ql_av
             ENDIF

          ENDIF

          IF ( ALLOCATED( qsurf_av ) )  THEN
             CALL wrd_write_string( 'qsurf_av' )
             WRITE ( 14 )  qsurf_av
          ENDIF

          IF ( ALLOCATED( qsws_av ) )  THEN
             CALL wrd_write_string( 'qsws_av' )
             WRITE ( 14 )  qsws_av
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN

          CALL wrd_write_string( 's' )
          WRITE ( 14 )  s

          IF ( ALLOCATED( s_av ) )  THEN
             CALL wrd_write_string( 's_av' )
             WRITE ( 14 )  s_av
          ENDIF

          IF ( ALLOCATED( ssurf_av ) )  THEN
             CALL wrd_write_string( 'ssurf_av' )
             WRITE ( 14 )  ssurf_av
          ENDIF

          IF ( ALLOCATED( ssws_av ) )  THEN
             CALL wrd_write_string( 'ssws_av' )
             WRITE ( 14 )  ssws_av
          ENDIF

       ENDIF

       IF ( ALLOCATED( ql_c_av ) )  THEN
          CALL wrd_write_string( 'ql_c_av' )
          WRITE ( 14 )  ql_c_av
       ENDIF

       IF ( ALLOCATED( ql_v_av ) )  THEN
          CALL wrd_write_string( 'ql_v_av' )
          WRITE ( 14 )  ql_v_av
       ENDIF

       IF ( ALLOCATED( ql_vp_av ) )  THEN
          CALL wrd_write_string( 'ql_vp_av' )
          WRITE ( 14 )  ql_vp_av
       ENDIF

       IF ( ALLOCATED( qv_av ) )  THEN
          CALL wrd_write_string( 'qv_av' )
          WRITE ( 14 )  qv_av
       ENDIF

       CALL wrd_write_string( 'random_iv' )
       WRITE ( 14 )  random_iv
       WRITE ( 14 )  random_iy

       IF ( ALLOCATED( seq_random_array ) )  THEN
          CALL wrd_write_string( 'seq_random_array' )
          WRITE ( 14 )  id_random_array
          WRITE ( 14 )  seq_random_array
       ENDIF

       IF ( ALLOCATED( s_av ) )  THEN
          CALL wrd_write_string( 's_av' )
          WRITE ( 14 )  s_av
       ENDIF

       IF ( ALLOCATED( shf_av ) )  THEN
          CALL wrd_write_string( 'shf_av' )
          WRITE ( 14 )  shf_av
       ENDIF

       IF ( ALLOCATED( ts_av ) )  THEN
          CALL wrd_write_string( 'ts_av' )
          WRITE ( 14 ) ts_av
       ENDIF

       CALL wrd_write_string( 'u' )
       WRITE ( 14 )  u

       IF ( ALLOCATED( u_av ) )  THEN
          CALL wrd_write_string( 'u_av' )
          WRITE ( 14 )  u_av
       ENDIF

       IF ( ALLOCATED( us_av ) )  THEN
          CALL wrd_write_string( 'us_av' )
          WRITE ( 14 )  us_av
       ENDIF

       CALL wrd_write_string( 'v' )
       WRITE ( 14 )  v

       IF ( ALLOCATED( v_av ) )  THEN
          CALL wrd_write_string( 'v_av' )
          WRITE ( 14 )  v_av
       ENDIF

       IF ( humidity )  THEN

          CALL wrd_write_string( 'vpt' )
          WRITE ( 14 )  vpt

          IF ( ALLOCATED( vpt_av ) )  THEN
             CALL wrd_write_string( 'vpt_av' )
             WRITE ( 14 ) vpt_av
          ENDIF

       ENDIF

       CALL wrd_write_string( 'w' )
       WRITE ( 14 )  w

       IF ( ALLOCATED( w_av ) )  THEN
          CALL wrd_write_string( 'w_av' )
          WRITE ( 14 )  w_av
       ENDIF

       IF ( ALLOCATED( z0_av ) )  THEN
          CALL wrd_write_string( 'z0_av' )
          WRITE ( 14 )  z0_av
       ENDIF

      IF ( ALLOCATED( z0h_av ) )  THEN
          CALL wrd_write_string( 'z0h_av' )
          WRITE ( 14 )  z0h_av
       ENDIF

       IF ( ALLOCATED( z0q_av ) )  THEN
          CALL wrd_write_string( 'z0q_av' )
          WRITE ( 14 )  z0q_av
       ENDIF

       IF ( land_surface .OR. urban_surface )  THEN

          IF ( ALLOCATED( ghf_av ) )  THEN
             CALL wrd_write_string( 'ghf_av' )
             WRITE ( 14 )  ghf_av
          ENDIF

          IF ( ALLOCATED( r_a_av ) )  THEN
             CALL wrd_write_string( 'r_a_av' )
             WRITE ( 14 )  r_a_av
          ENDIF

       ENDIF

       IF ( ALLOCATED( tsurf_av ) )  THEN
          CALL wrd_write_string( 'tsurf_av' )
          WRITE ( 14 )  tsurf_av
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--    Write local restart data using MPI-IO
!--    Attention: After changes to the following output commands the version number
!--    ---------  of the variable binary_version must be changed!
!--               Also, the list of arrays to be read in rrd_local must be adjusted accordingly.
       CALL wrd_mpi_io( 'e',  e)
       IF ( ALLOCATED( e_av ) )  CALL wrd_mpi_io( 'e_av',  e_av)
       CALL wrd_mpi_io( 'kh', kh )
       IF ( ALLOCATED( kh_av ) )  CALL wrd_mpi_io( 'kh_av', kh_av )
       CALL wrd_mpi_io( 'km' , km )
       IF ( ALLOCATED( km_av ) )  CALL wrd_mpi_io( 'km_av', km_av )
       IF ( ALLOCATED( lpt_av ) )  CALL wrd_mpi_io( 'lpt_av', lpt_av )
       IF ( ALLOCATED( lwp_av ) )  CALL wrd_mpi_io( 'lwp_av', lwp_av )
       IF ( ALLOCATED( ol_av ) )  CALL wrd_mpi_io( 'ol_av', ol_av )
       CALL wrd_mpi_io( 'p', p )
       IF ( ALLOCATED( p_av ) )  CALL wrd_mpi_io( 'p_av', p_av )
       IF ( ALLOCATED( pc_av ) )  CALL wrd_mpi_io( 'pc_av', pc_av )
       IF ( ALLOCATED( pr_av ) )  CALL wrd_mpi_io( 'pr_av', pr_av )
       IF ( ALLOCATED( pres_drag_x_av ) )  CALL wrd_mpi_io( 'pres_drag_x_av', pres_drag_x_av )
       IF ( ALLOCATED( pres_drag_y_av ) )  CALL wrd_mpi_io( 'pres_drag_y_av', pres_drag_y_av )
       CALL wrd_mpi_io( 'pt', pt )
       IF ( ALLOCATED( pt_av ) )  CALL wrd_mpi_io( 'pt_av', pt_av )

       IF ( humidity )  THEN

          CALL wrd_mpi_io( 'q', q )
          IF ( ALLOCATED( q_av ) )  CALL wrd_mpi_io( 'q_av', q_av )

          IF ( cloud_droplets )  THEN

             CALL wrd_mpi_io( 'ql', ql )
             IF ( ALLOCATED( ql_av ) )  CALL wrd_mpi_io( 'ql_av', ql_av )

          ENDIF

          IF ( ALLOCATED( qsurf_av ) )  CALL wrd_mpi_io( 'qsurf_av', qsurf_av )
          IF ( ALLOCATED( qsws_av ) )  CALL wrd_mpi_io( 'qsws_av', qsws_av )

       ENDIF

       IF ( passive_scalar )  THEN

          CALL wrd_mpi_io( 's', s )
          IF ( ALLOCATED( s_av ) )  CALL wrd_mpi_io( 's_av',  s_av )
          IF ( ALLOCATED( ssurf_av ) )  CALL wrd_mpi_io( 'ssurf_av', ssurf_av )
          IF ( ALLOCATED( ssws_av ) )  CALL wrd_mpi_io( 'ssws_av', ssws_av )

       ENDIF

       IF ( ALLOCATED( ql_c_av ) )  CALL wrd_mpi_io( 'ql_c_av', ql_c_av )
       IF ( ALLOCATED( ql_v_av ) )  CALL wrd_mpi_io( 'ql_v_av', ql_v_av )
       IF ( ALLOCATED( ql_vp_av ) )  CALL wrd_mpi_io( 'ql_vp_av', ql_vp_av )
       IF ( ALLOCATED( qv_av ) )  CALL wrd_mpi_io( 'qv_av', qv_av )
!
!--    Only PE0 writes random_iv and random_iy to restart file.
!--    ATTENTION: If one value for every PE is required, the general approach of PE indendent
!--    restart will be lost. That means that in general the parallel random number generator
!--    in random_generatot_parallel_mod should be used!
       CALL wrd_mpi_io_global_array( 'random_iv', random_iv )
       CALL wrd_mpi_io( 'random_iy', random_iy )

       IF ( ALLOCATED( seq_random_array ) )  THEN
          CALL wrd_mpi_io( 'id_random_array', id_random_array )
          DO  i = 1, SIZE( seq_random_array, 1 )
             WRITE( tmp_name, '(A,I2.2)' )  'seq_random_array', i
             CALL wrd_mpi_io( TRIM( tmp_name ), seq_random_array(i,:,:) )
          ENDDO
       ENDIF
       IF ( ALLOCATED( s_av ) )  CALL wrd_mpi_io( 's_av', s_av )
       IF ( ALLOCATED( shf_av ) )  CALL wrd_mpi_io( 'shf_av', shf_av )
       IF ( ALLOCATED( ts_av ) )  CALL wrd_mpi_io( 'ts_av', ts_av )
       CALL wrd_mpi_io( 'u', u)
       IF ( ALLOCATED( u_av ) )  CALL wrd_mpi_io( 'u_av', u_av )
       IF ( ALLOCATED( us_av ) )  CALL wrd_mpi_io( 'us_av', us_av )
       CALL wrd_mpi_io( 'v', v )
       IF ( ALLOCATED( v_av ) )  CALL wrd_mpi_io( 'v_av', v_av )
       IF ( humidity )  THEN
          CALL wrd_mpi_io( 'vpt',  vpt )
          IF ( ALLOCATED( vpt_av ) )  CALL wrd_mpi_io( 'vpt_av', vpt_av )
       ENDIF
       CALL wrd_mpi_io( 'w', w)
       IF ( ALLOCATED( w_av ) )  CALL wrd_mpi_io( 'w_av', w_av )
       IF ( ALLOCATED( z0_av ) )  CALL wrd_mpi_io( 'z0_av', z0_av )
       IF ( ALLOCATED( z0h_av ) )  CALL wrd_mpi_io( 'z0h_av', z0h_av )
       IF ( ALLOCATED( z0q_av ) )  CALL wrd_mpi_io( 'z0q_av', z0q_av )
       IF ( land_surface  .OR.  urban_surface )  THEN
          IF ( ALLOCATED( ghf_av ) )  CALL wrd_mpi_io( 'ghf_av', ghf_av )
          IF ( ALLOCATED( r_a_av ) )  CALL wrd_mpi_io( 'r_a_av', r_a_av )
       ENDIF
       IF ( ALLOCATED( tsurf_av ) )  CALL wrd_mpi_io( 'tsurf_av', tsurf_av )

    ENDIF

    CALL surface_wrd_local
!
!-- Write restart data of surface_data_output_mod
    IF ( surface_output )  CALL surface_data_output_wrd_local
!
!-- Write restart data of other modules
    CALL module_interface_wrd_local

!
!-- Write end label
    IF (  restart_data_format_output == 'fortran_binary' )  THEN
       CALL wrd_write_string( '*** end ***' )
    ENDIF

 END SUBROUTINE wrd_local

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write local variables to spinup file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_local_spinup

    CALL surface_wrd_local
!
!-- Write restart data of other modules
    CALL module_interface_wrd_local_spinup

 END SUBROUTINE wrd_local_spinup


 END MODULE write_restart_data_mod
