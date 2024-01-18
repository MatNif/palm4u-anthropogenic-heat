!> @file parin.f90
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
!> This subroutine reads variables controling the run from the NAMELIST files
!>
!> @todo: Revise max_pr_cs (profiles for chemistry)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE parin

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  pt_init,                                                                            &
               q_init,                                                                             &
               ref_state,                                                                          &
               s_init,                                                                             &
               sa_init,                                                                            &
               ug,                                                                                 &
               u_init,                                                                             &
               v_init,                                                                             &
               vg

    USE control_parameters

    USE cpulog,                                                                                    &
        ONLY:  cpu_log_barrierwait

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               ny,                                                                                 &
               nz

    USE kinds

    USE model_1d_mod,                                                                              &
        ONLY:  damp_level_1d,                                                                      &
               dt_pr_1d,                                                                           &
               dt_run_control_1d,                                                                  &
               end_time_1d

    USE module_interface,                                                                          &
        ONLY:  module_interface_parin

    USE netcdf_interface,                                                                          &
        ONLY:  netcdf_data_format,                                                                 &
               netcdf_deflate,                                                                     &
               netcdf_precision

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  nesting_bounds

    USE profil_parameter,                                                                          &
        ONLY:  cross_profiles,                                                                     &
               profile_columns,                                                                    &
               profile_rows

    USE progress_bar,                                                                              &
        ONLY:  progress_bar_disabled

    USE read_restart_data_mod,                                                                     &
        ONLY:  rrd_global

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               hom_sum,                                                                            &
               pr_max,                                                                             &
               statistic_regions

    USE surface_data_output_mod,                                                                   &
        ONLY: surface_data_output_parin

    USE turbulence_closure_mod,                                                                    &
        ONLY:  rans_const_c,                                                                       &
               rans_const_sigma


    IMPLICIT NONE

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter
                                 !< file
#if defined( __parallel )
    CHARACTER(LEN=100) ::  root_initializing_actions = ' '  !< local copy of the root value of
                                                            !< initializing_actions to be
                                                            !< broadcast to all children
#endif

    INTEGER(iwp) ::  global_id      !< process id with respect to MPI_COMM_WORLD
    INTEGER(iwp) ::  global_procs   !< # of procs with respect to MPI_COMM_WORLD
    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  io_status      !< status after reading the namelist files


    NAMELIST /initialization_parameters/  advanced_div_correction,                                 &
                                          allow_negative_scalar_values,                            &
                                          allow_roughness_limitation,                              &
                                          alpha_surface,                                           &
                                          approximation,                                           &
                                          bc_e_b,                                                  &
                                          bc_lr,                                                   &
                                          bc_ns,                                                   &
                                          bc_p_b,                                                  &
                                          bc_p_t,                                                  &
                                          bc_pt_b,                                                 &
                                          bc_pt_t,                                                 &
                                          bc_q_b,                                                  &
                                          bc_q_t,                                                  &
                                          bc_s_b,                                                  &
                                          bc_s_t,                                                  &
                                          bc_uv_b,                                                 &
                                          bc_uv_t,                                                 &
                                          building_height,                                         &
                                          building_length_x,                                       &
                                          building_length_y,                                       &
                                          building_wall_left,                                      &
                                          building_wall_south,                                     &
                                          calc_soil_moisture_during_spinup,                        &
                                          call_psolver_at_all_substeps,                            &
                                          canyon_height,                                           &
                                          canyon_wall_left,                                        &
                                          canyon_wall_south,                                       &
                                          canyon_width_x,                                          &
                                          canyon_width_y,                                          &
                                          cfl_factor,                                              &
                                          check_realistic_q,                                       &
                                          cloud_droplets,                                          &
                                          collective_wait,                                         &
                                          conserve_volume_flow,                                    &
                                          conserve_volume_flow_mode,                               &
                                          constant_flux_layer,                                     &
                                          coupling_start_time,                                     &
                                          cycle_mg,                                                &
                                          damp_level_1d,                                           &
                                          data_output_during_spinup,                               &
                                          dissipation_1d,                                          &
                                          dp_external,                                             &
                                          dp_level_b,                                              &
                                          dp_smooth, dpdxy,                                        &
                                          dt,                                                      &
                                          dt_pr_1d,                                                &
                                          dt_run_control_1d,                                       &
                                          dt_spinup,                                               &
                                          dx,                                                      &
                                          dy,                                                      &
                                          dz,                                                      &
                                          dz_max,                                                  &
                                          dz_stretch_factor,                                       &
                                          dz_stretch_level,                                        &
                                          dz_stretch_level_start,                                  &
                                          dz_stretch_level_end,                                    &
                                          e_init,                                                  &
                                          e_min,                                                   &
                                          end_time_1d,                                             &
                                          ensemble_member_nr,                                      &
                                          fft_method,                                              &
                                          flux_input_mode,                                         &
                                          flux_output_mode,                                        &
                                          galilei_transformation,                                  &
                                          homogenize_surface_temperature,                          &
                                          humidity,                                                &
                                          inflow_disturbance_begin,                                &
                                          inflow_disturbance_end,                                  &
                                          initializing_actions,                                    &
                                          km_constant,                                             &
                                          large_scale_forcing,                                     &
                                          large_scale_subsidence,                                  &
                                          latitude,                                                &
                                          longitude,                                               &
                                          loop_optimization,                                       &
                                          lsf_exception,                                           &
                                          masking_method,                                          &
                                          mg_cycles,                                               &
                                          mg_switch_to_pe0_level,                                  &
                                          mixing_length_1d,                                        &
                                          momentum_advec,                                          &
                                          monotonic_limiter_z,                                     &
                                          netcdf_precision,                                        &
                                          neutral,                                                 &
                                          ngsrb,                                                   &
                                          nsor,                                                    &
                                          nsor_ini,                                                &
                                          nudging,                                                 &
                                          nx,                                                      &
                                          ny,                                                      &
                                          nz,                                                      &
                                          omega,                                                   &
                                          omega_sor,                                               &
                                          origin_date_time,                                        &
                                          outflow_damping_factor,                                  &
                                          outflow_damping_width,                                   &
                                          outflow_source_plane,                                    &
                                          passive_scalar,                                          &
                                          prandtl_number,                                          &
                                          psolver,                                                 &
                                          pt_damping_factor,                                       &
                                          pt_damping_width,                                        &
                                          pt_reference,                                            &
                                          pt_surface,                                              &
                                          pt_surface_heating_rate,                                 &
                                          pt_surface_initial_change,                               &
                                          pt_vertical_gradient,                                    &
                                          pt_vertical_gradient_level,                              &
                                          q_surface,                                               &
                                          q_surface_initial_change,                                &
                                          q_vertical_gradient,                                     &
                                          q_vertical_gradient_level,                               &
                                          random_generator,                                        &
                                          random_heatflux,                                         &
                                          rans_const_c,                                            &
                                          rans_const_sigma,                                        &
                                          rayleigh_damping_factor,                                 &
                                          rayleigh_damping_height,                                 &
                                          reference_state,                                         &
                                          residual_limit,                                          &
                                          restart_data_format,                                     &
                                          restart_data_format_input,                               &
                                          restart_data_format_output,                              &
                                          rotation_angle,                                          &
                                          roughness_length,                                        &
                                          s_surface,                                               &
                                          s_surface_initial_change,                                &
                                          s_vertical_gradient,                                     &
                                          s_vertical_gradient_level,                               &
                                          scalar_advec,                                            &
                                          scalar_rayleigh_damping,                                 &
                                          spinup_time,                                             &
                                          spinup_pt_amplitude,                                     &
                                          spinup_pt_mean,                                          &
                                          statistic_regions,                                       &
                                          subs_vertical_gradient,                                  &
                                          subs_vertical_gradient_level,                            &
                                          surface_heatflux,                                        &
                                          surface_pressure,                                        &
                                          surface_scalarflux,                                      &
                                          surface_waterflux,                                       &
                                          terrain_following_mapping,                               &
                                          timestep_scheme,                                         &
                                          topography,                                              &
                                          topography_grid_convention,                              &
                                          top_heatflux,                                            &
                                          top_momentumflux_u,                                      &
                                          top_momentumflux_v,                                      &
                                          top_scalarflux,                                          &
                                          tunnel_height,                                           &
                                          tunnel_length,                                           &
                                          tunnel_wall_depth,                                       &
                                          tunnel_width_x,                                          &
                                          tunnel_width_y,                                          &
                                          turbulence_closure,                                      &
                                          turbulent_outflow,                                       &
                                          u_bulk,                                                  &
                                          u_profile,                                               &
                                          ug_surface,                                              &
                                          ug_vertical_gradient,                                    &
                                          ug_vertical_gradient_level,                              &
                                          use_fixed_date,                                          &
                                          use_fixed_time,                                          &
                                          use_free_convection_scaling,                             &
                                          use_ug_for_galilei_tr,                                   &
                                          use_subsidence_tendencies,                               &
                                          use_surface_fluxes,                                      &
                                          use_top_fluxes,                                          &
                                          use_upstream_for_tke,                                    &
                                          uv_heights,                                              &
                                          v_bulk,                                                  &
                                          v_profile,                                               &
                                          vdi_checks,                                              &
                                          vg_surface,                                              &
                                          vg_vertical_gradient,                                    &
                                          vg_vertical_gradient_level,                              &
                                          wall_adjustment,                                         &
                                          wall_heatflux,                                           &
                                          wall_humidityflux,                                       &
                                          wall_scalarflux,                                         &
                                          y_shift,                                                 &
                                          zeta_max,                                                &
                                          zeta_min,                                                &
                                          z0h_factor

    NAMELIST /runtime_parameters/  averaging_interval,                                             &
                                   averaging_interval_pr,                                          &
                                   bc_pt_b,                                                        &
                                   cpu_log_barrierwait,                                            &
                                   create_disturbances,                                            &
                                   cross_profiles,                                                 &
                                   data_output,                                                    &
                                   data_output_2d_on_each_pe,                                      &
                                   data_output_masks,                                              &
                                   data_output_pr,                                                 &
                                   debug_output,                                                   &
                                   debug_output_timestep,                                          &
                                   disturbance_amplitude,                                          &
                                   disturbance_energy_limit,                                       &
                                   disturbance_level_b,                                            &
                                   disturbance_level_t,                                            &
                                   do2d_at_begin,                                                  &
                                   do3d_at_begin,                                                  &
                                   dt,                                                             &
                                   dt_averaging_input,                                             &
                                   dt_averaging_input_pr,                                          &
                                   dt_coupling,                                                    &
                                   dt_data_output,                                                 &
                                   dt_data_output_av,                                              &
                                   dt_disturb,                                                     &
                                   dt_domask,                                                      &
                                   dt_dopr,                                                        &
                                   dt_dopr_listing,                                                &
                                   dt_dots,                                                        &
                                   dt_do2d_xy,                                                     &
                                   dt_do2d_xz,                                                     &
                                   dt_do2d_yz,                                                     &
                                   dt_do3d,                                                        &
                                   dt_max,                                                         &
                                   dt_restart,                                                     &
                                   dt_run_control,                                                 &
                                   end_time,                                                       &
                                   force_print_header,                                             &
                                   homogenize_surface_temperature,                                 &
                                   interpolate_to_grid_center,                                     &
                                   mask_k_over_surface,                                            &
                                   mask_scale_x,                                                   &
                                   mask_scale_y,                                                   &
                                   mask_scale_z,                                                   &
                                   mask_x,                                                         &
                                   mask_y,                                                         &
                                   mask_z,                                                         &
                                   mask_x_loop,                                                    &
                                   mask_y_loop,                                                    &
                                   mask_z_loop,                                                    &
                                   netcdf_data_format,                                             &
                                   netcdf_deflate,                                                 &
                                   normalizing_region,                                             &
                                   npex,                                                           &
                                   npey,                                                           &
                                   nz_do3d,                                                        &
                                   open_debug_files,                                               &
                                   profile_columns,                                                &
                                   profile_rows,                                                   &
                                   pt_surface_initial_change,                                      &
                                   pt_surface_heating_rate,                                        &
                                   restart_time,                                                   &
                                   restart_data_format,                                            &
                                   restart_data_format_input,                                      &
                                   restart_data_format_output,                                     &
                                   section_xy,                                                     &
                                   section_xy_m,                                                   &
                                   section_xz,                                                     &
                                   section_xz_m,                                                   &
                                   section_yz,                                                     &
                                   section_yz_m,                                                   &
                                   skip_time_data_output,                                          &
                                   skip_time_data_output_av,                                       &
                                   skip_time_dopr,                                                 &
                                   skip_time_do2d_xy,                                              &
                                   skip_time_do2d_xz,                                              &
                                   skip_time_do2d_yz,                                              &
                                   skip_time_do3d,                                                 &
                                   skip_time_domask,                                               &
                                   surface_heatflux,                                               &
                                   synchronous_exchange,                                           &
                                   termination_time_needed,                                        &
                                   wall_heatflux

    NAMELIST /envpar/  host,                                                                       &
                       maximum_cpu_time_allowed,                                                   &
                       maximum_parallel_io_streams,                                                &
                       progress_bar_disabled,                                                      &
                       read_svf,                                                                   &
                       version_string,                                                             &
                       run_identifier,                                                             &
                       tasks_per_node,                                                             &
                       write_binary,                                                               &
                       write_spinup_data,                                                          &
                       write_svf
!
!-- Set default value for initializing_actions in case of child domain in a nested run
!-- before the namelist initialization_parameters is read.
#if defined( __parallel )
    IF ( nested_run )  THEN
       IF ( child_domain )  THEN
          initializing_actions = 'interpolate_from_parent'
       ENDIF
    ENDIF
#endif
!
!-- First read values of environment variables (this NAMELIST file is generated by palmrun)
    CALL location_message( 'reading environment parameters from ENVPAR', 'start' )

    OPEN( 90, FILE='ENVPAR', STATUS='OLD', FORM='FORMATTED', IOSTAT=io_status )

    IF ( io_status /= 0 )  THEN
       message_string = 'Local file ENVPAR not found.' //                                          &
                        '&Some variables for steering may not be properly set.'
       CALL message( 'parin', 'PAC0253', 0, 1, 0, 6, 0 )
    ELSE
       READ( 90, envpar, IOSTAT=io_status )
       IF ( io_status < 0 )  THEN
          message_string = 'No envpar-NAMELIST found in local file ENVPAR.&'  //                   &
                           'Some variables for steering may not be properly set.'
          CALL message( 'parin', 'PAC0255', 0, 1, 0, 6, 0 )
       ELSEIF ( io_status > 0 )  THEN
          BACKSPACE(90 )
          READ( 90 , '(A)' ) line
          message_string = 'Errors in local file ENVPAR within line:&' // line //                  &
                           '&Some variables for steering may not be properly set.'
          CALL message( 'parin', 'PAC0254', 0, 1, 0, 6, 0 )
       ENDIF
       CLOSE( 90 )
    ENDIF

    CALL location_message( 'reading environment parameters from ENVPAR', 'finished' )

    CALL location_message( 'reading NAMELIST parameters from PARIN', 'start' )

!
!-- Open the NAMELIST-file which is send with this job
    CALL check_open( 11 )

!
!-- Read the control parameters for initialization.
!-- The namelist "initialisation_parameters" must be provided in the NAMELIST-file.
    READ( 11, initialization_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status > 0 )  THEN
!
!--    initialisation_parameters namelist was found but countained errors. Print an error message
!-     including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'initialization_parameters', line )

    ELSEIF ( io_status < 0 )  THEN
!
!--    initialisation_parametes namelist was not found. Return a message.
       message_string = 'no initialization_parameters namelist found'
       CALL message( 'parin', 'PAC0256', 1, 2, 0, 6, 0 )

    ENDIF
!
!-- Read surface data output namelist
    CALL surface_data_output_parin
!
!-- Check for module namelists and read them. Please note, reading of module-specific namelists
!-- must be performed before module specific global restart data is read. This is because
!-- the control flags switching the modules on/off are set here.
    CALL module_interface_parin
!
!-- Set internal switches based on initialization parameter initializing_actions.
    IF ( INDEX( initializing_actions, 'cyclic_fill' ) /= 0 )  cyclic_fill_initialization = .TRUE.

!
!-- Calculate the number of groups into which parallel I/O is split.
!-- The default for files which are opened by all PEs (or where each PE opens its own independent
!-- file) is, that all PEs are doing input/output in parallel at the same time. This might cause
!-- performance or even more severe problems depending on the configuration of the underlying file
!-- system.
!-- Calculation of the number of blocks and the I/O group must be based on all PEs involved in this
!-- run. Since myid and numprocs are related to the comm2d communicator, which gives only a subset
!-- of all PEs in case of nested runs, that information must be inquired again from the global
!-- communicator.
!-- First, set the default:
#if defined( __parallel )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, global_id, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, global_procs, ierr )
#else
    global_id    = 0
    global_procs = 1
#endif
    IF ( maximum_parallel_io_streams == -1  .OR.  maximum_parallel_io_streams > global_procs )  THEN
       maximum_parallel_io_streams = global_procs
    ENDIF
!
!-- Now calculate the number of io_blocks and the io_group to which the respective PE belongs. I/O
!-- of the groups is done in serial, but in parallel for all PEs belonging to the same group.
    io_blocks = global_procs / maximum_parallel_io_streams
    io_group  = MOD( global_id+1, io_blocks )
!
!-- In case of nested run, make sure that initializing_actions = 'interpolate_from_parent' which
!-- leads to the same branch in init_3d_model as 'set_constant_profiles' for all the children. The
!-- set_constant_profiles-branch in the init_3d_model is executed even though the constant-profiles
!-- initializations for the prognostic variables will be overwritten by pmci_child_initialize
!-- and pmci_parent_initialize. This is, however, important to make sure that e.g. diagnostic
!-- variables are set properly. An exception is made in case of restart runs, if cyclic-fill is
!-- used, and if user decides to initialize everything. A copy of the root value is first broadcast
!-- to all children so that if the root has 'read_restart_data', it will be used also for all the
!-- children. This arrangement means that initializing_actions never can't be 'by user' for any
!-- child domain.
#if defined( __parallel )
    IF ( nested_run )  THEN

       IF ( .NOT. child_domain )  THEN
          root_initializing_actions = initializing_actions
       ENDIF

       CALL MPI_BCAST( root_initializing_actions, LEN( root_initializing_actions ),                &
                       MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

       IF ( child_domain )  THEN

          IF ( INDEX( root_initializing_actions, 'read_restart_data' ) /= 0 )  THEN
             IF ( INDEX( initializing_actions, 'read_restart_data') == 0 )  THEN
                message_string = 'initializing_actions = ' // TRIM( initializing_actions ) //      &
                                 ' has been changed to "read_restart_data" in child domain.'
                CALL message( 'parin', 'PAC0257', 0, 0, 0, 6, 0 )
                initializing_actions = root_initializing_actions
             ENDIF
          ENDIF

          IF ( .NOT. ( INDEX( initializing_actions, 'read_restart_data' ) /= 0  .OR.               &
                       INDEX( initializing_actions, 'by_user' ) /= 0            .OR.               &
                       cyclic_fill_initialization ) )                                              &
          THEN
             IF ( INDEX( initializing_actions, 'interpolate_from_parent') == 0 )  THEN
                message_string = 'initializing_actions = ' // TRIM( initializing_actions ) //      &
                                 ' has been changed to "interpolate_from_parent" in child domain.'
                CALL message( 'parin', 'PAC0257', 0, 0, 0, 6, 0 )
                initializing_actions = 'interpolate_from_parent'
             ENDIF
          ENDIF

       ENDIF

    ENDIF
#endif
!
!-- Try to read runtime parameters given by the user for this run (namelist "&runtime_parameters").
!-- This namelist can be omitted. In that case default values are used for the parameters.
    REWIND( 11 )
    READ( 11, runtime_parameters, IOSTAT=io_status )

    IF ( io_status > 0 )  THEN
!
!--    Namelist runtime_parameters was found but contained errors. Print an error message including
!--    the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'runtime_parameters', line )

    ENDIF

!
!-- If required, read control parameters from restart file (produced by a prior run). All PEs are
!-- reading from file created by PE0 (see check_open)
    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN

!
!--    If not set by the user in the runtime parameters, the data format for restart input needs to
!--    be set now! This is normally done later in check parameters.
       IF ( TRIM( restart_data_format ) /= 'fortran_binary'  .AND.                                  &
            TRIM( restart_data_format ) /= 'mpi'             .AND.                                  &
            TRIM( restart_data_format ) /= 'mpi_shared_memory' )  THEN
          message_string = 'illegal restart data format "' // TRIM( restart_data_format ) // '"'
          CALL message( 'parin', 'PAC0013', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( TRIM( restart_data_format_input ) == 'undefined' )  THEN
          restart_data_format_input = restart_data_format
       ENDIF

!
!--    Blockwise I/O does not work together with MPI-I/O
       IF ( restart_data_format_input(1:3) == 'mpi'  .AND.  io_blocks /= 1 )  THEN
          CALL rrd_global
       ELSE
!
!--       Data is read in parallel by groups of PEs
          DO  i = 0, io_blocks-1
             IF ( i == io_group )  THEN
                CALL rrd_global
             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO
       ENDIF

!
!--    Increment the run count
       runnr = runnr + 1
!
!--    In case of a restart run, the number of user-defined profiles on the restart file (already
!--    stored in max_pr_user) has to match the one given for the current run. max_pr_user_tmp is
!--    calculated in user_parin and max_pr_user is read in via rrd_global.
       IF ( max_pr_user /= max_pr_user_tmp )  THEN
          WRITE( message_string, * ) 'the number of user-defined ',                                &
                                     'profiles given in data_output_pr (', max_pr_user_tmp,        &
                                     ') does not match the one ',                                  &
                                     'found in the restart file (', max_pr_user, ')'
          CALL message( 'user_parin', 'USR0008', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read runtime parameters again, in order to follow the rule that they have to overwrite the
!--    parameter values that have been read from the restart file.
       REWIND( 11 )
       READ( 11, runtime_parameters )

    ELSE

       max_pr_user = max_pr_user_tmp

    ENDIF

!
!-- Activate spinup
    IF ( land_surface  .OR.  urban_surface )  THEN
       IF ( spinup_time > 0.0_wp )  THEN
          coupling_start_time = spinup_time
          time_since_reference_point = simulated_time - coupling_start_time
          IF ( spinup_pt_mean == 9999999.9_wp )  THEN
             spinup_pt_mean = pt_surface
          ENDIF
          end_time = end_time + spinup_time
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
             spinup = .TRUE.
          ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.                      &
                   time_since_reference_point > 0.0_wp )                                           &
          THEN
             data_output_during_spinup = .FALSE.  !< required for correct ntdim calculation
                                                  !< in check_parameters for restart run
          ENDIF
       ENDIF
    ENDIF

!
!-- In case of nested runs, explicitly set nesting boundary conditions along x and y for childs
!-- in case of 3d-nesting. Top boundary conditions for childs need to be set 'nested', except for
!-- pressure, which requires a Neumann condition.
!-- This will overwrite any user settings and basic defaults for childs.
    IF ( nested_run  .AND.  child_domain )  THEN
!
!--    Set conditions along x.
       IF ( .NOT. ( nesting_bounds == 'vertical_only'  .OR.  nesting_bounds == 'cyclic_along_x' ) )&
       THEN
          bc_lr = 'nested'
       ENDIF
!
!--    Set conditions along y.
       IF ( .NOT. ( nesting_bounds == 'vertical_only'  .OR.  nesting_bounds == 'cyclic_along_y' ) )&
       THEN
          bc_ns = 'nested'
       ENDIF

    ENDIF
!
!-- Set boundary conditions also in case the model is offline-nested in larger-scale models.
!>  @todo boundary conditions for scalars, chemical-, and aerosol-species need to be added here.
    IF ( nesting_offline )  THEN
       bc_lr   = 'nesting_offline'
       bc_ns   = 'nesting_offline'
       bc_uv_t = 'nesting_offline'
       bc_pt_t = 'nesting_offline'
       bc_q_t  = 'nesting_offline'
       bc_p_t  = 'neumann'
    ENDIF
!
!-- Check validity of lateral boundary conditions. This has to be done here because they are already
!-- used in init_pegrid and init_grid and therefore cannot be checked in check_parameters.
    IF ( bc_lr /= 'cyclic'               .AND.  bc_lr /= 'dirichlet/radiation'  .AND.              &
         bc_lr /= 'radiation/dirichlet'  .AND.  bc_lr /= 'nested'               .AND.              &
         bc_lr /= 'nesting_offline' )  THEN
       message_string = 'unknown boundary condition: bc_lr = "' // TRIM( bc_lr ) // '"'
       CALL message( 'parin', 'PAC0258', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( bc_ns /= 'cyclic'               .AND.  bc_ns /= 'dirichlet/radiation'  .AND.              &
         bc_ns /= 'radiation/dirichlet'  .AND.  bc_ns /= 'nested'               .AND.              &
         bc_ns /= 'nesting_offline' )  THEN
       message_string = 'unknown boundary condition: bc_ns = "' // TRIM( bc_ns ) // '"'
       CALL message( 'parin', 'PAC0259', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set internal variables used for speed optimization in if clauses
    IF ( bc_lr /= 'cyclic' )               bc_lr_cyc    = .FALSE.
    IF ( bc_lr == 'dirichlet/radiation' )  bc_lr_dirrad = .TRUE.
    IF ( bc_lr == 'radiation/dirichlet' )  bc_lr_raddir = .TRUE.
    IF ( bc_ns /= 'cyclic' )               bc_ns_cyc    = .FALSE.
    IF ( bc_ns == 'dirichlet/radiation' )  bc_ns_dirrad = .TRUE.
    IF ( bc_ns == 'radiation/dirichlet' )  bc_ns_raddir = .TRUE.
!
!-- Radiation-Dirichlet conditions are allowed along one of the horizontal directions only.
!-- In general, such conditions along x and y may work, but require a) some code changes (e.g.
!-- concerning allocation of c_u, c_v ... arrays), and b) a careful model setup by the user, in
!-- order to guarantee that there is a clearly defined outflow at two sides.
!-- Otherwise, the radiation condition may produce wrong results.
    IF ( ( bc_lr_dirrad .OR. bc_lr_raddir )  .AND.  ( bc_ns_dirrad .OR. bc_ns_raddir ) )  THEN
       message_string = 'bc_lr = "' // TRIM( bc_lr ) // '" and bc_ns = "' // TRIM( bc_ns ) //      &
                        '" are not allowed to be set at the same time'
       CALL message( 'parin', 'PAC0260', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check in case of initial run, if the grid point numbers are well defined and allocate some
!-- arrays which are already needed in init_pegrid or check_parameters. During restart jobs, these
!-- arrays will be allocated in rrd_global. All other arrays are allocated in init_3d_model.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

       IF ( nx <= 0 )  THEN
          WRITE( message_string, * ) 'no value or wrong value given', ' for nx: nx=', nx
          CALL message( 'parin', 'PAC0261', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( ny <= 0 )  THEN
          WRITE( message_string, * ) 'no value or wrong value given', ' for ny: ny=', ny
          CALL message( 'parin', 'PAC0262', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( nz <= 0 )  THEN
          WRITE( message_string, * ) 'no value or wrong value given', ' for nz: nz=', nz
          CALL message( 'parin', 'PAC0263', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    As a side condition, routine flow_statistics require at least 14 vertical grid levels (they
!--    are used to store time-series data)
!>     @todo   Remove this restriction
       IF ( nz < 14 )  THEN
          WRITE( message_string, * ) 'nz < 14'
          CALL message( 'parin', 'PAC0264', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    ATTENTION: in case of changes to the following statement please also check the allocate
!--               statement in routine rrd_global
       ALLOCATE( pt_init(0:nz+1), q_init(0:nz+1), s_init(0:nz+1),                                  &
                 ref_state(0:nz+1), sa_init(0:nz+1), ug(0:nz+1),                                   &
                 u_init(0:nz+1), v_init(0:nz+1), vg(0:nz+1),                                       &
                 hom(0:nz+1,2,pr_max,0:statistic_regions),                                         &
                 hom_sum(0:nz+1,pr_max,0:statistic_regions) )

       hom = 0.0_wp
       hom_sum = 0.0_wp

    ENDIF

!
!-- NAMELIST-file is not needed anymore
    CALL close_file( 11 )

    CALL location_message( 'reading NAMELIST parameters from PARIN', 'finished' )

 END SUBROUTINE parin
