!> @file check_parameters.f90
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
!> Check control parameters and deduce further quantities.
!
!> @todo Increase character length of unit and corresponding characters to LEN>=8 in order to allow
!>       units like degree_C (05.08.2020)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE check_parameters

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d

    USE basic_constants_and_equations_mod

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters

    USE grid_variables

    USE kinds

    USE indices

    USE model_1d_mod,                                                                              &
        ONLY:  damp_level_1d,                                                                      &
               damp_level_ind_1d

    USE module_interface,                                                                          &
        ONLY:  module_interface_check_data_output,                                                 &
               module_interface_check_data_output_pr,                                              &
               module_interface_check_data_output_ts,                                              &
               module_interface_check_parameters

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  init_model, input_pids_static,                                                      &
               netcdf_data_input_check_dynamic,                                                    &
               netcdf_data_input_check_static

    USE netcdf_interface,                                                                          &
        ONLY:  do2d_unit,                                                                          &
               do3d_unit,                                                                          &
               dopr_unit,                                                                          &
               dots_max,                                                                           &
               dots_num,                                                                           &
               dots_unit,                                                                          &
               heatflux_output_unit,                                                               &
               momentumflux_output_unit,                                                           &
               netcdf_data_format,                                                                 &
               netcdf_data_format_string,                                                          &
               waterflux_output_unit

    USE particle_attributes,                                                                       &
        ONLY:  particle_advection,                                                                 &
               use_sgs_for_particles

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run,                                                       &
               cpl_id,                                                                             &
               lower_left_coord_x,                                                                 &
               lower_left_coord_y

    USE profil_parameter

    USE statistics

    USE subsidence_mod

    USE surface_data_output_mod,                                                                   &
        ONLY: surface_data_output_check_parameters

    IMPLICIT NONE

    CHARACTER (LEN=varnamelength)  ::  var           !< variable name
    CHARACTER (LEN=7)   ::  unit                     !< unit of variable
    CHARACTER (LEN=8)   ::  date                     !< current date string
    CHARACTER (LEN=10)  ::  time                     !< current time string
    CHARACTER (LEN=20)  ::  ensemble_string          !< string containing number of ensemble member
    CHARACTER (LEN=15)  ::  nest_string              !< string containing id of nested domain
    CHARACTER (LEN=40)  ::  coupling_string          !< string containing type of coupling
    CHARACTER (LEN=100) ::  action                   !< flag string

    CHARACTER (LEN=varnamelength), DIMENSION(:), ALLOCATABLE ::  data_output_filtered  !< filtered list of output variables

    INTEGER(iwp) ::  i                               !< loop index
    INTEGER(iwp) ::  ilen                            !< string length
    INTEGER(iwp) ::  j                               !< loop index
    INTEGER(iwp) ::  k                               !< loop index
    INTEGER(iwp) ::  kk                              !< loop index
    INTEGER(iwp) ::  mid                             !< masked output running index
    INTEGER(iwp) ::  netcdf_data_format_save         !< initial value of netcdf_data_format
    INTEGER(iwp) ::  next_index                      !< index of next free section array element
    INTEGER(iwp) ::  position                        !< index position of string
#if defined( __parallel )
    INTEGER(iwp) ::  comm_node                       !< local communicator required to determine number of PEs per node
    INTEGER(iwp) ::  max_npes_per_node               !< maximum number pf PEs per node
    INTEGER(iwp) ::  npes_per_node                   !< local number of PEs per node
#endif
    INTEGER(iwp) ::  nr_unique_sections              !< number of unique layers in section output

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  section_tmp !< dummy array used to sort out double occurrences in output sections

    LOGICAL     ::  file_exists                      !< flag checking if a file exists
    LOGICAL     ::  found                            !< flag, true if output variable is already marked for averaging
    
    LOGICAL     ::  found_err                        !< global flag, true if error is found on one core
    LOGICAL     ::  found_err_l = .FALSE.            !< local flag, true if error is found on one core

    REAL(wp)    ::  gradient                         !< local gradient
    REAL(wp)    ::  section_position                 !< relative position of cross sections in m
#if defined( __parallel )
    REAL(wp)    ::  dt_spinup_max                    !< maximum spinup timestep in nested domains
    REAL(wp)    ::  spinup_time_max                  !< maximum spinup time in nested domains
#endif


    CALL location_message( 'checking parameters', 'start' )
!
!-- At first, check static and dynamic input for consistency.
    CALL netcdf_data_input_check_dynamic
    CALL netcdf_data_input_check_static


!
!-- Check if humidity is set to .TRUE. in case of atmospheric run that is coupled to ocean.
    IF ( atmosphere_run_coupled_to_ocean  .AND.  .NOT. humidity ) THEN
       message_string = 'humidity has to be set to .T. in the atmosphere _p3d file ' //            &
                        'for coupled runs between ocean and atmosphere.'
       CALL message( 'check_parameters', 'PAC0012', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check and set the restart data format variables
    IF ( TRIM( restart_data_format ) /= 'fortran_binary'  .AND.                                    &
         TRIM( restart_data_format ) /= 'mpi'             .AND.                                    &
         TRIM( restart_data_format ) /= 'mpi_shared_memory' )  THEN
       message_string = 'illegal restart data format "' // TRIM( restart_data_format ) // '"'
       CALL message( 'check_parameters', 'PAC0013', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( restart_data_format_input ) == 'undefined' )  THEN
       restart_data_format_input = restart_data_format
    ENDIF
    IF ( TRIM( restart_data_format_output ) == 'undefined' )  THEN
       restart_data_format_output = restart_data_format
    ENDIF

    IF ( TRIM( restart_data_format_input ) /= 'fortran_binary'  .AND.                              &
         TRIM( restart_data_format_input ) /= 'mpi'             .AND.                              &
         TRIM( restart_data_format_input ) /= 'mpi_shared_memory' )  THEN
       message_string = 'illegal restart input data format "' //                                   &
                        TRIM( restart_data_format_input ) // '"'
       CALL message( 'check_parameters', 'PAC0014', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( restart_data_format_output ) /= 'fortran_binary'  .AND.                             &
         TRIM( restart_data_format_output ) /= 'mpi'             .AND.                             &
         TRIM( restart_data_format_output ) /= 'mpi_shared_memory' )  THEN
       message_string = 'illegal restart output data format "' //                                  &
                        TRIM( restart_data_format_output ) // '"'
       CALL message( 'check_parameters', 'PAC0015', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( restart_data_format_output ) == 'mpi_shared_memory'  .AND.  write_binary  .AND.     &
         particle_advection )                                                                      &
    THEN
       message_string = 'mpi_shared_memory is not implemented for particle I/O'
       CALL message( 'check_parameters', 'PAC0348', 1, 2, 0, 6, 0 )
    ENDIF

#if defined( __parallel )
!
!-- Restart I/O with shared memory MPI is not efficient for small setups and processors with
!-- less or equal 8 cores/node.
    IF ( TRIM( restart_data_format_input ) == 'mpi_shared_memory'  .OR.                            &
         TRIM( restart_data_format_output ) == 'mpi_shared_memory' )  THEN
!
!--    Check the number of available cores per node.
       CALL MPI_COMM_SPLIT_TYPE( comm2d, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_node, ierr )
       CALL MPI_COMM_SIZE( comm_node, npes_per_node, ierr )
       CALL MPI_ALLREDUCE( npes_per_node, max_npes_per_node, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
       CALL MPI_COMM_FREE( comm_node, ierr )

       IF ( max_npes_per_node <= 8  .OR.  nx < 64  .OR.  ny < 64 )  THEN
          IF ( TRIM( restart_data_format_input ) == 'mpi_shared_memory' )  THEN
             restart_data_format_input = 'mpi'
          ENDIF
          IF ( TRIM( restart_data_format_output ) == 'mpi_shared_memory' )  THEN
             restart_data_format_output = 'mpi'
          ENDIF
       ENDIF
    ENDIF
#endif

!
!-- User settings for restart times requires that "restart" has been given as file activation
!-- string. Otherwise, binary output would not be saved by palmrun.
    IF ( ( restart_time /= 9999999.9_wp  .OR.  dt_restart /= 9999999.9_wp )                        &
         .AND.  .NOT. write_binary )  THEN
       WRITE( message_string, * ) 'manual restart settings requires file ',                        &
                                  'activation string "restart"'
       CALL message( 'check_parameters', 'PAC0016', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set time for the next user defined restart (time_restart is the internal parameter for steering
!-- restart events)
    IF ( restart_time /= 9999999.9_wp )  THEN
!
!--    Set the internal parameter only if the given restart time lies after the start time
!--    of the curent run. If it lies in the past, the value of the internal parameter is
!--    determined by the restart data from the previous run (if it is a restart run).
       IF ( restart_time > time_since_reference_point )  THEN
          time_restart = restart_time
       ENDIF
    ELSE
!
!--    In case of a restart run, set internal parameter to default (no restart) if the
!--    NAMELIST-parameter restart_time is omitted
       time_restart = 9999999.9_wp
    ENDIF

!
!-- Generate the file header which is used as a header for most of PALM's output files
    CALL DATE_AND_TIME( date, time, run_zone )
    run_date = date(1:4) // '-' // date(5:6) // '-' // date(7:8)
    run_time = time(1:2) // ':' // time(3:4) // ':' // time(5:6)

    coupling_string = ''
    IF ( ensemble_member_nr /= 0 )  THEN
       WRITE( ensemble_string, '(2X,A,I2.2)' )  'en-no: ', ensemble_member_nr
    ELSE
       ensemble_string = ''
    ENDIF
    IF ( nested_run )  THEN
       WRITE( nest_string, '(2X,A,I2.2)' )  'nest-id: ', cpl_id
    ELSEIF ( atmosphere_ocean_coupled_run )  THEN
       coupling_string = 'atmos_ocean'
       nest_string     = '_coupling'
    ELSE
       nest_string = ''
    ENDIF

    WRITE ( run_description_header, '(A,2X,A,A,A,I2.2,A,A,A,2X,A,A,2X,A,1X,A)' )                   &
          TRIM( version_string ), 'run: ', TRIM( run_identifier ), '.', runnr,                     &
          TRIM( coupling_string ), TRIM( nest_string ), TRIM( ensemble_string), 'host: ',          &
          TRIM( host ), run_date, run_time

!
!-- Check the general loop optimization method
    SELECT CASE ( TRIM( loop_optimization ) )

       CASE ( 'cache', 'vector' )
          CONTINUE

       CASE DEFAULT
          message_string = 'illegal value given for loop_optimization: "' //                       &
                           TRIM( loop_optimization ) // '"'
          CALL message( 'check_parameters', 'PAC0017', 1, 2, 0, 6, 0 )

    END SELECT

!
!-- Check topography setting (check for illegal parameter combinations)
    IF ( topography /= 'flat' )  THEN
       action = ' '
       IF ( scalar_advec /= 'pw-scheme'  .AND.  scalar_advec /= 'ws-scheme' )  THEN
          WRITE( action, '(A,A)' )  'scalar_advec = ', scalar_advec
       ENDIF
       IF ( momentum_advec /= 'pw-scheme'  .AND.  momentum_advec /= 'ws-scheme' )  THEN
          WRITE( action, '(A,A)' )  'momentum_advec = ', momentum_advec
       ENDIF
       IF ( psolver == 'sor' )  THEN
          WRITE( action, '(A,A)' )  'psolver = ', psolver
       ENDIF
       IF ( sloping_surface )  THEN
          WRITE( action, '(A)' )  'sloping surface = .TRUE.'
       ENDIF
       IF ( galilei_transformation )  THEN
          WRITE( action, '(A)' )  'galilei_transformation = .TRUE.'
       ENDIF
       IF ( cloud_droplets )  THEN
          WRITE( action, '(A)' )  'cloud_droplets = .TRUE.'
       ENDIF
       IF ( .NOT. constant_flux_layer  .AND.  topography /= 'closed_channel' )  THEN
          WRITE( action, '(A)' )  'constant_flux_layer = .FALSE.'
       ENDIF
       IF ( action /= ' ' )  THEN
          message_string = 'the specified topography does not allow ' // TRIM( action )
          CALL message( 'check_parameters', 'PAC0018', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check illegal/untested parameter combinations for closed channel
       If ( topography == 'closed_channel' )  THEN
          symmetry_flag = 1
          message_string = 'bottom and top boundary are treated equal'
          CALL message( 'check_parameters', 'PAC0019', 0, 0, 0, 6, 0 )

          IF ( dz(1) /= dz(COUNT( dz /= -1.0_wp ))  .OR.  dz_stretch_level /= -9999999.9_wp)  THEN
             WRITE( message_string, * )  'dz should be equal close to the ' //                     &
                                         'boundaries due to symmetrical problem'
             CALL message( 'check_parameters', 'PAC0020', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( constant_flux_layer )  THEN
             WRITE( message_string, * )  'a constant flux layer is not ' //                        &
                                         'allowed if a closed channel shall be used'
             CALL message( 'check_parameters', 'PAC0021', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( ocean_mode )  THEN
             WRITE( message_string, * )  'ocean mode is not allowed if ' //                        &
                                         'a closed channel shall be used'
             CALL message( 'check_parameters', 'PAC0022', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( momentum_advec /= 'ws-scheme'  .OR.                                                 &
               scalar_advec /= 'ws-scheme' )  THEN
             WRITE( message_string, * )  'closed channel requires Wicker-Skamarock as advection ', &
                                         'scheme'
             CALL message( 'check_parameters', 'PAC0023', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Check approximation
    IF ( TRIM( approximation ) /= 'boussinesq'  .AND.  TRIM( approximation ) /= 'anelastic' )  THEN
       message_string = 'unknown approximation: approximation = "' // TRIM( approximation ) // '"'
       CALL message( 'check_parameters', 'PAC0024', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check approximation requirements
    IF ( TRIM( approximation ) == 'anelastic'  .AND.  TRIM( momentum_advec ) /= 'ws-scheme' )  THEN
       message_string = 'anelastic approximation requires momentum_advec = "ws-scheme"'
       CALL message( 'check_parameters', 'PAC0025', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( approximation ) == 'anelastic'  .AND.  conserve_volume_flow )  THEN
       message_string = 'anelastic approximation is not allowed with conserve_volume_flow = .TRUE.'
       CALL message( 'check_parameters', 'PAC0026', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check flux input mode
    IF ( TRIM( flux_input_mode ) /= 'dynamic'  .AND.  TRIM( flux_input_mode ) /= 'kinematic'       &
         .AND.  TRIM( flux_input_mode ) /= 'application-specific' )  THEN
       message_string = 'unknown flux input mode: flux_input_mode = "' //                          &
                        TRIM( flux_input_mode ) // '"'
       CALL message( 'check_parameters', 'PAC0027', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set meaningful default flux input mode if nothing else is prescribed
    IF ( TRIM( flux_input_mode ) == 'application-specific' )  THEN
!
!--    Set flux input mode according to approximation if applicable
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_input_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_input_mode = 'kinematic'
       ENDIF
!
!--    When the land- or urban-surface model is used, the flux input must be dynamic
       IF ( land_surface  .OR.  urban_surface )  THEN
          flux_input_mode = 'dynamic'
       ENDIF
    ENDIF

!
!-- Check flux output mode
    IF ( TRIM( flux_output_mode ) /= 'dynamic'  .AND.  TRIM( flux_output_mode ) /= 'kinematic'     &
         .AND.  TRIM( flux_output_mode ) /= 'application-specific' )  THEN
       message_string = 'unknown flux output mode: flux_output_mode = "' //                        &
                        TRIM( flux_output_mode ) // '"'
       CALL message( 'check_parameters', 'PAC0028', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if maximum number of allowed timeseries is exceeded
    IF ( dots_num > dots_max )  THEN
       WRITE( message_string, * ) 'number of time series quantities exceeds',                      &
                                  ' its maximum of dots_max = ', dots_max
       CALL message( 'check_parameters', 'PAC0029', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check whether there are any illegal values
!-- Pressure solver:
    IF ( psolver /= 'poisfft'  .AND.  psolver /= 'poisfft_sm' .AND.  psolver /= 'sor'  .AND.       &
         psolver /= 'multigrid'  .AND.  psolver /= 'multigrid_noopt' )                             &
    THEN
       message_string = 'unknown solver for perturbation pressure: psolver = "' //                 &
                        TRIM( psolver ) // '"'
       CALL message( 'check_parameters', 'PAC0030', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( psolver(1:9) == 'multigrid' )  THEN
       IF ( cycle_mg == 'w' )  THEN
          gamma_mg = 2
       ELSEIF ( cycle_mg == 'v' )  THEN
          gamma_mg = 1
       ELSE
          message_string = 'unknown multigrid cycle: cycle_mg = "' //  TRIM( cycle_mg ) // '"'
          CALL message( 'check_parameters', 'PAC0031', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( fft_method /= 'singleton-algorithm'  .AND.  fft_method /= 'temperton-algorithm'  .AND.    &
         fft_method /= 'fftw'                 .AND.  fft_method /= 'system-specific' )  THEN
       message_string = 'unknown fft-algorithm: fft_method = "' // TRIM( fft_method ) // '"'
       CALL message( 'check_parameters', 'PAC0032', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( ( fft_method == 'fftw'  .OR.  fft_method ==  'singleton-algorithm' )  .AND.               &
         ( ( npey == 1  .AND.  npex > 1 )  .OR.  ( npex == 1 .AND. npey > 1 ) )  .AND.             &
         loop_optimization == 'vector' )                                                           &
    THEN
       message_string = 'fft_method = "' // TRIM( fft_method ) //                                  &
                        '" is not available for loop_optimization =& "vector" with ' //            &
                        '1d-domain-decomposition'
       CALL message( 'check_parameters', 'PAC0033', 1, 2, 0, 6, 0 )
    ENDIF

    IF( momentum_advec == 'ws-scheme' .AND.  .NOT. call_psolver_at_all_substeps  ) THEN
        message_string = 'psolver must be called at each RK3 substep when "'//                     &
                         TRIM(momentum_advec) // ' "is used for momentum_advec'
        CALL message( 'check_parameters', 'PAC0034', 1, 2, 0, 6, 0 )
    END IF
!
!-- Advection schemes:
    IF ( momentum_advec /= 'pw-scheme'  .AND.  momentum_advec /= 'ws-scheme'  .AND.                &
         momentum_advec /= 'up-scheme' )  THEN
       message_string = 'unknown momentum advection scheme = "' // TRIM( momentum_advec ) // '"'
       CALL message( 'check_parameters', 'PAC0035', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ( momentum_advec == 'ws-scheme' .OR. scalar_advec == 'ws-scheme' )                        &
         .AND. ( timestep_scheme == 'euler' .OR.  timestep_scheme == 'runge-kutta-2' ) )  THEN
       message_string = 'momentum_advec or scalar_advec = "' // TRIM( momentum_advec ) //          &
                        '" is not allowed with timestep_scheme = "' //                             &
                        TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PAC0036', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec /= 'pw-scheme'  .AND.  scalar_advec /= 'ws-scheme'  .AND.                    &
         scalar_advec /= 'bc-scheme'  .AND.  scalar_advec /= 'up-scheme' )  THEN
       message_string = 'unknown scalar advection scheme = "' // TRIM( scalar_advec ) // '"'
       CALL message( 'check_parameters', 'PAC0037', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec == 'bc-scheme'  .AND.  loop_optimization == 'cache' )  THEN
       message_string = 'advection_scheme scalar_advec = "' // TRIM( scalar_advec ) //             &
                        '" not implemented for loop_optimization = "' //                           &
                        TRIM( loop_optimization ) // '"'
       CALL message( 'check_parameters', 'PAC0038', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets  .AND.  .NOT. use_upstream_for_tke     &
         .AND.  scalar_advec /= 'ws-scheme' )  THEN
       use_upstream_for_tke = .TRUE.
       message_string = 'use_upstream_for_tke is set to .TRUE. because ' //                        &
                        'use_sgs_for_particles = .TRUE. and scalar_advec /= ws-scheme'
       CALL message( 'check_parameters', 'PAC0039', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Perform checks for the advanced divergence correction method. This is only implemented in the
!-- vector-machine-optimized branch. Moreover, with this method, the Courant number should not
!-- exceed 0.5 (empirically derived value).
    IF ( scalar_advec == 'ws-scheme'  .AND.  advanced_div_correction  .AND.                        &
         TRIM( loop_optimization ) /= 'vector' )  THEN
       message_string = 'advanced flow-divergence correction is only possible with &'    //        &
                        'loop_optimization = vector'
       CALL message( 'check_parameters', 'PAC0040', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( scalar_advec == 'ws-scheme'  .AND.  advanced_div_correction  .AND.                        &
         cfl_factor >= 0.5_wp )  THEN
       message_string = 'Advanced flow-divergence correction requires a smaller timestep. &' //    &
                        'Please set cfl_factor to a value < 0.5.'
       CALL message( 'check_parameters', 'PAC0041', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set LOGICAL switches to enhance performance
    IF ( momentum_advec == 'ws-scheme' )  ws_scheme_mom = .TRUE.
    IF ( scalar_advec   == 'ws-scheme' )  ws_scheme_sca = .TRUE.


!
!-- Timestep schemes:
    SELECT CASE ( TRIM( timestep_scheme ) )

       CASE ( 'euler' )
          intermediate_timestep_count_max = 1

       CASE ( 'runge-kutta-2' )
          intermediate_timestep_count_max = 2

       CASE ( 'runge-kutta-3' )
          intermediate_timestep_count_max = 3

       CASE DEFAULT
          message_string = 'unknown timestep scheme: timestep_scheme = "' //                       &
                           TRIM( timestep_scheme ) // '"'
          CALL message( 'check_parameters', 'PAC0042', 1, 2, 0, 6, 0 )

    END SELECT

    IF ( ( momentum_advec /= 'pw-scheme' .AND. momentum_advec /= 'ws-scheme' )                     &
         .AND.  timestep_scheme(1:5) == 'runge' ) THEN
       message_string = 'momentum advection scheme "' // TRIM( momentum_advec ) //                 &
                        '" & does not work with timestep_scheme "' // TRIM( timestep_scheme )      &
                        // '"'
       CALL message( 'check_parameters', 'PAC0043', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Initializing actions must have been set by the user
    IF ( TRIM( initializing_actions ) == '' )  THEN
       message_string = 'no value specified for initializing_actions'
       CALL message( 'check_parameters', 'PAC0044', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.                                &
         .NOT. cyclic_fill_initialization )                                                        &
    THEN
!
!--    No restart run or main run following a pre-run: several initialising actions are possible
       action = initializing_actions
       DO  WHILE ( TRIM( action ) /= '' )
          position = INDEX( action, ' ' )
          SELECT CASE ( action(1:position-1) )

             CASE ( 'set_constant_profiles', 'set_1d-model_profiles', 'by_user',                   &
                    'initialize_vortex', 'initialize_ptanom', 'initialize_bubble',                 &
                    'interpolate_from_parent', 'read_from_file', 'read_spinup_data' )
                action = action(position+1:)

             CASE DEFAULT
                message_string = 'initializing_action = "' //                                      &
                                 TRIM( action ) // '" unknown or not allowed'
                CALL message( 'check_parameters', 'PAC0046', 1, 2, 0, 6, 0 )

          END SELECT
       ENDDO
    ENDIF

    IF ( TRIM( initializing_actions ) == 'initialize_vortex'  .AND.  conserve_volume_flow ) THEN
         message_string = 'initializing_actions = "initialize_vortex"' //                          &
                          ' is not allowed with conserve_volume_flow = .T.'
       CALL message( 'check_parameters', 'PAC0047', 1, 2, 0, 6, 0 )
    ENDIF


    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.                        &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //                        &
                        ' and "set_1d-model_profiles" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PAC0048', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.                        &
         INDEX( initializing_actions, 'by_user' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //                        &
                        ' and "by_user" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PAC0049', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'by_user' ) /= 0  .AND.                                      &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "by_user" and ' //                                 &
                        '"set_1d-model_profiles" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PAC0050', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Only allow setting of spinup_time if land- and/or urban-surface model is switched on.
    IF ( spinup_time > 0.0_wp  .AND.  .NOT. ( land_surface  .OR.  urban_surface ) )  THEN
       message_string = 'spinup_time > 0.0 is not allowed when both the land- and urban-' //       &
                        'surface model are switched-off'
       CALL message( 'check_parameters', 'PAC0349', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set flag for reading surface spinup data.
    read_spinup_data = INDEX( initializing_actions, 'read_spinup_data' ) /= 0
!
!-- Check if the SPINUPIN file is present within the working directory.
    INQUIRE( FILE = 'SPINUPIN' // TRIM( coupling_char ), EXIST = file_exists )
    IF ( read_spinup_data  .AND.  .NOT. file_exists )  THEN
       message_string = 'initializing_actions = "read_spinup_data" requires a spinup file ' //     &
                        'within the working directory'
       CALL message( 'check_parameters', 'PAC0051', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Reading spinup data together with non-zero spinup time is not allowed
    IF ( read_spinup_data  .AND.  spinup )  THEN
       message_string = 'initializing_actions = "read_spinup_data" requires spinup_time = 0.0'
       CALL message( 'check_parameters', 'PAC0052', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Reading spinup data is only allowed in initial runs, not in restart runs.
    IF ( read_spinup_data  .AND.  INDEX( initializing_actions, 'read_restart_data' ) /= 0 )  THEN
       message_string = 'initializing_actions = "read_spinup_data read_restart_data"' //           &
                        'is not possible'
       CALL message( 'check_parameters', 'PAC0053', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Reading spinup data with Fortran IO is not possible
    IF ( read_spinup_data  .AND.  TRIM( restart_data_format_input ) == 'fortran_binary' )  THEN
       message_string = 'reading spinup data requires restart_data_format_input = ' //             &
                        '"mpi" or "mpi_shared_memory"'
       CALL message( 'check_parameters', 'PAC0054', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Writing spinup data with Fortran IO is not possible
    IF ( write_spinup_data  .AND.  TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
       message_string = 'writing spinup data requires restart_data_format_input = ' //             &
                        '"mpi" or "mpi_shared_memory"'
       CALL message( 'check_parameters', 'PAC0055', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Give a warning if surface spinup data shall be written out but no spinup has been enabled.
    IF ( write_spinup_data  .AND.  spinup_time == 0.0_wp )  THEN
       message_string = 'writing spinup data is enabled but surface_spinup is switched-off'
       CALL message( 'check_parameters', 'PAC0056', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- In case of spinup and nested run, spinup end time must be identical in order to have
!-- synchronously running simulations.
    IF ( nested_run )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( spinup_time, spinup_time_max, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD,     &
                           ierr )
       CALL MPI_ALLREDUCE( dt_spinup,   dt_spinup_max,   1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD,     &
                           ierr )

       IF ( spinup_time /= spinup_time_max  .OR.  dt_spinup /= dt_spinup_max )  THEN
          message_string = 'nesting, spinup_time, and dt_spinup are not identical in all parent' //&
                           ' andd child domains'
          CALL message( 'check_parameters', 'PAC0057', 3, 2, 0, 6, 0 )
       ENDIF
#endif
    ENDIF

    IF ( humidity  .AND.  sloping_surface )  THEN
       message_string = 'humidity = .TRUE. and sloping_surface = .TRUE. ' //                       &
                        'are not allowed simultaneously'
       CALL message( 'check_parameters', 'PAC0058', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check parameters of surface data output in surface_data_output_mod and in module interface
    IF ( surface_output )  THEN
       CALL surface_data_output_check_parameters
    ENDIF

!-- Check the module settings
    CALL module_interface_check_parameters

!
!-- In case of no restart run, check initialising parameters and calculate further quantities
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

!
!--    Initial profiles for 1D and 3D model, respectively (u,v further below)
       pt_init = pt_surface
       IF ( humidity       )  q_init  = q_surface
       IF ( passive_scalar )  s_init  = s_surface

!--
!--    If required, compute initial profile of the geostrophic wind (component ug)
       i = 1
       gradient = 0.0_wp

       IF ( .NOT. ocean_mode )  THEN

          ug_vertical_gradient_level_ind(1) = 0
          ug(0) = ug_surface
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) < zu(k)  .AND.                                  &
                     ug_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   ug(k) = ug(k-1) + dzu(k) * gradient
                ELSE
                   ug(k) = ug_surface + dzu(k) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k-1)
             ENDIF
          ENDDO

       ELSE

          ug_vertical_gradient_level_ind(1) = nzt+1
          ug(nzt+1) = ug_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) > zu(k)  .AND.                                  &
                     ug_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   ug(k) = ug(k+1) - dzu(k+1) * gradient
                ELSE
                   ug(k)   = ug_surface - 0.5_wp * dzu(k+1) * gradient
                   ug(k+1) = ug_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k+1)
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients for ug, choose a zero gradient
       IF ( ug_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          ug_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--
!--    If required, compute initial profile of the geostrophic wind (component vg)
       i = 1
       gradient = 0.0_wp

       IF ( .NOT. ocean_mode )  THEN

          vg_vertical_gradient_level_ind(1) = 0
          vg(0) = vg_surface
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) < zu(k)  .AND.                                  &
                     vg_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   vg(k) = vg(k-1) + dzu(k) * gradient
                ELSE
                   vg(k) = vg_surface + dzu(k) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k-1)
             ENDIF
          ENDDO

       ELSE

          vg_vertical_gradient_level_ind(1) = nzt+1
          vg(nzt+1) = vg_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) > zu(k)  .AND.                                  &
                     vg_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   vg(k) = vg(k+1) - dzu(k+1) * gradient
                ELSE
                   vg(k)   = vg_surface - 0.5_wp * dzu(k+1) * gradient
                   vg(k+1) = vg_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k+1)
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients for vg, choose a zero gradient
       IF ( vg_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          vg_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--    Let the initial wind profiles be the calculated ug/vg profiles or interpolate them from wind
!--    profile data (if given)
       IF ( u_profile(1) == 9999999.9_wp  .AND.  v_profile(1) == 9999999.9_wp )  THEN

          u_init = ug
          v_init = vg

       ELSEIF ( u_profile(1) == 0.0_wp  .AND.  v_profile(1) == 0.0_wp )  THEN

          IF ( uv_heights(1) /= 0.0_wp )  THEN
             message_string = 'uv_heights(1) must be 0.0'
             CALL message( 'check_parameters', 'PAC0059', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( omega /= 0.0_wp )  THEN
             message_string = 'Coriolis force is switched on (omega /= 0.0) and' //                &
                              ' initial profiles u_profile and v_profile are prescribed'
             CALL message( 'check_parameters', 'PAC0060', 0, 0, 0, 6, 0 )
          ENDIF

          use_prescribed_profile_data = .TRUE.

          kk = 1
          u_init(0) = 0.0_wp
          v_init(0) = 0.0_wp

          DO  k = 1, nz+1

             IF ( kk < 200 )  THEN
                DO  WHILE ( uv_heights(kk+1) <= zu(k) )
                   kk = kk + 1
                   IF ( kk == 200 )  EXIT
                ENDDO
             ENDIF

             IF ( kk < 200  .AND.  uv_heights(kk+1) /= 9999999.9_wp )  THEN
                u_init(k) = u_profile(kk) + ( zu(k) - uv_heights(kk) ) /                           &
                                            ( uv_heights(kk+1) - uv_heights(kk) ) *                &
                                            ( u_profile(kk+1) - u_profile(kk) )
                v_init(k) = v_profile(kk) + ( zu(k) - uv_heights(kk) ) /                           &
                                            ( uv_heights(kk+1) - uv_heights(kk) ) *                &
                                            ( v_profile(kk+1) - v_profile(kk) )
             ELSE
                u_init(k) = u_profile(kk)
                v_init(k) = v_profile(kk)
             ENDIF

          ENDDO

       ELSE

          message_string = 'u_profile(1) and v_profile(1) must be 0.0'
          CALL message( 'check_parameters', 'PAC0061', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Compute initial temperature profile using the given temperature gradients
       IF (  .NOT.  neutral )  THEN
          CALL init_vertical_profiles( pt_vertical_gradient_level_ind, pt_vertical_gradient_level, &
                                       pt_vertical_gradient, pt_init, pt_surface, bc_pt_t_val )
       ENDIF
!
!--    Compute initial humidity profile using the given humidity gradients
       IF ( humidity )  THEN
          CALL init_vertical_profiles( q_vertical_gradient_level_ind, q_vertical_gradient_level,   &
                                       q_vertical_gradient, q_init, q_surface, bc_q_t_val )
       ENDIF
!
!--    Compute initial scalar profile using the given scalar gradients
       IF ( passive_scalar )  THEN
          CALL init_vertical_profiles( s_vertical_gradient_level_ind, s_vertical_gradient_level,   &
                                       s_vertical_gradient, s_init, s_surface, bc_s_t_val )
       ENDIF
!
!--    TODO
!--    Compute initial chemistry profile using the given chemical species gradients
!--    Russo: Is done in chem_init --> kanani: Revise

    ENDIF

!
!-- Check if the control parameter use_subsidence_tendencies is used correctly
    IF ( use_subsidence_tendencies  .AND.  .NOT.  large_scale_subsidence )  THEN
       message_string = 'usage of use_subsidence_tendencies ' //                                   &
                        'requires large_scale_subsidence = .TRUE.'
       CALL message( 'check_parameters', 'PAC0062', 1, 2, 0, 6, 0 )
    ELSEIF ( use_subsidence_tendencies  .AND.  .NOT. large_scale_forcing )  THEN
       message_string = 'usage of use_subsidence_tendencies ' //                                   &
                        'requires large_scale_forcing = .TRUE.'
       CALL message( 'check_parameters', 'PAC0062', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize large scale subsidence if required
    If ( large_scale_subsidence )  THEN
       IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp  .AND. .NOT. large_scale_forcing )    &
       THEN
          CALL init_w_subsidence
       ENDIF
!
!--    In case large_scale_forcing is used, profiles for subsidence velocity are read in from file
!--    LSF_DATA.
       IF ( subs_vertical_gradient_level(1) == -9999999.9_wp  .AND. .NOT. large_scale_forcing )    &
       THEN
          message_string = 'missing large scale vertical velocity profile'
          CALL message( 'check_parameters', 'PAC0063', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
        IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp )  THEN
           message_string = 'large scale subsidence is not activated'
          CALL message( 'check_parameters', 'PAC0064', 1, 2, 0, 6, 0 )
        ENDIF
    ENDIF

!
!-- Overwrite parameters from namelist if necessary and compute Coriolis parameter.
!-- @todo - move initialization of f and fs to coriolis_mod.
    IF ( input_pids_static )  THEN
       latitude       = init_model%latitude
       longitude      = init_model%longitude
       rotation_angle = init_model%rotation_angle
    ENDIF

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 2.0_wp * omega * COS( latitude / 180.0_wp * pi )

!
!-- Check and set buoyancy related parameters and switches
    IF ( reference_state == 'horizontal_average' )  THEN
       CONTINUE
    ELSEIF ( reference_state == 'initial_profile' )  THEN
       use_initial_profile_as_reference = .TRUE.
    ELSEIF ( reference_state == 'single_value' )  THEN
       use_single_reference_value = .TRUE.
       IF ( pt_reference == 9999999.9_wp )  pt_reference = pt_surface
       vpt_reference = pt_reference * ( 1.0_wp + 0.61_wp * q_surface )
    ELSE
       message_string = 'illegal value for reference_state: "' // TRIM( reference_state ) // '"'
       CALL message( 'check_parameters', 'PAC0065', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- In case of a given slope, compute the relevant quantities
    IF ( alpha_surface /= 0.0_wp )  THEN
       IF ( ABS( alpha_surface ) > 90.0_wp )  THEN
          WRITE( message_string, * ) 'ABS( alpha_surface = ', alpha_surface, ' ) must be < 90.0'
          CALL message( 'check_parameters', 'PAC0066', 1, 2, 0, 6, 0 )
       ENDIF
       sloping_surface = .TRUE.
       cos_alpha_surface = COS( alpha_surface / 180.0_wp * pi )
       sin_alpha_surface = SIN( alpha_surface / 180.0_wp * pi )
    ENDIF

!
!-- Check time step and cfl_factor
    IF ( dt /= -1.0_wp )  THEN
       IF ( dt <= 0.0_wp )  THEN
          WRITE( message_string, * ) 'dt = ', dt , ' <= 0.0'
          CALL message( 'check_parameters', 'PAC0067', 1, 2, 0, 6, 0 )
       ENDIF
       dt_3d = dt
       dt_fixed = .TRUE.
    ENDIF

    IF ( cfl_factor <= 0.0_wp  .OR.  cfl_factor > 1.0_wp )  THEN
       IF ( cfl_factor == -1.0_wp )  THEN
          IF ( timestep_scheme == 'runge-kutta-2' )  THEN
             cfl_factor = 0.8_wp
          ELSEIF ( timestep_scheme == 'runge-kutta-3' )  THEN
             cfl_factor = 0.9_wp
          ELSE
             cfl_factor = 0.9_wp
          ENDIF
       ELSE
          WRITE( message_string, * ) 'cfl_factor = ', cfl_factor, ' out of range'
          CALL message( 'check_parameters', 'PAC0068', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Store simulated time at begin
    simulated_time_at_begin = simulated_time

!
!-- Store reference time for coupled runs
    IF ( simulated_time == 0.0_wp  .AND.  coupling_start_time == 0.0_wp )  THEN
       time_since_reference_point = 0.0_wp
    ENDIF

!
!-- Set wind speed in the Galilei-transformed system
    IF ( galilei_transformation )  THEN
       IF ( use_ug_for_galilei_tr                    .AND.                                         &
            ug_vertical_gradient_level(1) == 0.0_wp  .AND.                                         &
            ug_vertical_gradient(1)       == 0.0_wp  .AND.                                         &
            vg_vertical_gradient_level(1) == 0.0_wp  .AND.                                         &
            vg_vertical_gradient(1)       == 0.0_wp )  THEN
          u_gtrans = ug_surface * 0.6_wp
          v_gtrans = vg_surface * 0.6_wp
       ELSEIF ( use_ug_for_galilei_tr  .AND.  ( ug_vertical_gradient_level(1) /= 0.0_wp .OR.       &
                                                ug_vertical_gradient(1) /= 0.0_wp ) )  THEN
          message_string = 'galilei transformation does not allow height depending ug'
          CALL message( 'check_parameters', 'PAC0069', 1, 2, 0, 6, 0 )
       ELSEIF ( use_ug_for_galilei_tr  .AND.  ( vg_vertical_gradient_level(1) /= 0.0_wp  .OR.      &
                                                vg_vertical_gradient(1) /= 0.0_wp ) )  THEN
          message_string = 'galilei transformation does not allow height depending vg'
          CALL message( 'check_parameters', 'PAC0069', 1, 2, 0, 6, 0 )
       ELSE
          message_string = 'variable translation speed used for Galilei-transformation, which ' // &
                           'may cause & instabilities in stably stratified regions'
          CALL message( 'check_parameters', 'PAC0070', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- In case of using a constant flux layer, calculated (or prescribed) surface fluxes have to be
!-- used in the diffusion-terms
    IF ( constant_flux_layer )  use_surface_fluxes = .TRUE.

!
!-- The Poisson-FFT-solver requires even number of grid points along x and y in case of cylic
!-- boundary conditions.
    IF ( psolver(1:7) == 'poisfft'  .AND.  bc_lr_cyc )  THEN
       IF ( MOD( nx+1, 2 ) == 1 )  THEN
          message_string = 'psolver = "poisfft" or "poisfft_sm" requires the number of grid ' //   &
                           'points& along x-direction (nx+1) to be even '
          CALL message( 'check_parameters', 'PAC0071', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( psolver(1:7) == 'poisfft'  .AND.  bc_ns_cyc )  THEN
       IF ( MOD( ny+1, 2 ) == 1 )  THEN
          message_string = 'psolver = "poisfft" or "poisfft_sm" requires the number of grid ' //  &
                           'points & along y-direction (ny+1) to be even '
          CALL message( 'check_parameters', 'PAC0072', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Check boundary conditions and set internal variables:
!-- Attention: the lateral boundary conditions have been already checked in parin
!
!-- Non-cyclic lateral boundaries require specific FFT methods and Piascek-Williams or
!-- Wicker-Skamarock advection scheme. Several schemes and tools do not work with non-cyclic
!-- boundary conditions.
    IF ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )  THEN
       IF ( psolver(1:9) /= 'multigrid' )  THEN
!
!--       Non-cyclic boundary conditions are realized in the pressure solver only for specific FFTs
          IF ( fft_method /= 'fftw'  .AND.  fft_method /= 'temperton-algorithm'  .AND.             &
               fft_method /= 'singleton-algorithm' )  THEN
             message_string = 'non-cyclic lateral boundaries do not allow ' // 'psolver = "' //    &
                              TRIM( psolver ) // '" with FFT method = "' //                        &
                              TRIM( fft_method ) // '"'
             CALL message( 'check_parameters', 'PAC0073', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( momentum_advec /= 'pw-scheme'  .AND.  momentum_advec /= 'ws-scheme' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow momentum_advec = "' //      &
                           TRIM( momentum_advec ) // '"'
          CALL message( 'check_parameters', 'PAC0074', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( scalar_advec /= 'pw-scheme'  .AND.  scalar_advec /= 'ws-scheme' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow scalar_advec = "' //        &
                           TRIM( scalar_advec ) // '"'
          CALL message( 'check_parameters', 'PAC0075', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( galilei_transformation )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow galilei_transformation = .T.'
          CALL message( 'check_parameters', 'PAC0076', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( ( ( bc_lr_cyc  .AND.  .NOT.  bc_ns_cyc )  .OR.  ( .NOT. bc_lr_cyc  .AND.  bc_ns_cyc ) )   &
         .AND.  psolver(1:9) /= 'multigrid' )  THEN
       message_string = 'psolver = "poisfft" requires that boundary conditions along x and y ' //  &
                        'are both either cyclic or non-cyclic'
       CALL message( 'check_parameters', 'PAC0077', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Bottom boundary condition for the turbulent Kinetic energy
    IF ( bc_e_b == 'neumann' )  THEN
       ibc_e_b = 1
    ELSEIF ( bc_e_b == '(u*)**2+neumann' )  THEN
       ibc_e_b = 2
       IF ( .NOT. constant_flux_layer )  THEN
          bc_e_b = 'neumann'
          ibc_e_b = 1
          message_string = 'boundary condition bc_e_b changed to "' // TRIM( bc_e_b ) // '"'
          CALL message( 'check_parameters', 'PAC0078', 0, 1, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_e_b = "' // TRIM( bc_e_b ) // '"'
       CALL message( 'check_parameters', 'PAC0079', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Bottom boundary condition for perturbation pressure.
    IF ( bc_p_b == 'dirichlet' )  THEN
       ibc_p_b = 0
    ELSEIF ( bc_p_b == 'neumann' )  THEN
       ibc_p_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_b = "' // TRIM( bc_p_b ) // '"'
       CALL message( 'check_parameters', 'PAC0080', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Top boundary condition for perturbation pressure. The default value depends on the selected
!-- pressure solver. A Neumann-condition is avoided for the multigrid solver since it may cause
!-- convergence problems in some cases.
    IF ( bc_p_t == 'default' )  THEN
       IF ( psolver(1:7) == 'poisfft' )  THEN
          bc_p_t = 'neumann'
       ELSE
          bc_p_t = 'dirichlet'
       ENDIF
    ENDIF

    IF ( bc_p_t == 'dirichlet' )  THEN
       ibc_p_t = 0
    ELSEIF ( bc_p_t == 'neumann' .OR. bc_p_t == 'nested' )  THEN
       ibc_p_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_t = "' // TRIM( bc_p_t ) // '"'
       CALL message( 'check_parameters', 'PAC0081', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for potential temperature
    IF ( atmosphere_run_coupled_to_ocean )  THEN
       ibc_pt_b = 2
    ELSE
       IF ( bc_pt_b == 'dirichlet' )  THEN
          ibc_pt_b = 0
       ELSEIF ( bc_pt_b == 'neumann' )  THEN
          ibc_pt_b = 1
       ELSEIF ( bc_pt_b == 'nested' )  THEN
          ibc_pt_b = 3
       ELSE
          message_string = 'unknown boundary condition: bc_pt_b = "' // TRIM( bc_pt_b ) // '"'
          CALL message( 'check_parameters', 'PAC0082', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_pt_t == 'dirichlet' )  THEN
       ibc_pt_t = 0
    ELSEIF ( bc_pt_t == 'neumann' )  THEN
       ibc_pt_t = 1
    ELSEIF ( bc_pt_t == 'initial_gradient' )  THEN
       ibc_pt_t = 2
    ELSEIF ( bc_pt_t == 'nested'  .OR.  bc_pt_t == 'nesting_offline' )  THEN
       ibc_pt_t = 3
    ELSE
       message_string = 'unknown boundary condition: bc_pt_t = "' // TRIM( bc_pt_t ) // '"'
       CALL message( 'check_parameters', 'PAC0083', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( ANY( wall_heatflux /= 0.0_wp )  .AND.  surface_heatflux == 9999999.9_wp )  THEN
       message_string = 'wall_heatflux additionally requires setting of surface_heatflux'
       CALL message( 'check_parameters', 'PAC0084', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- This IF clause needs revision, got too complex!!
    IF ( surface_heatflux == 9999999.9_wp  )  THEN
       constant_heatflux = .FALSE.
       IF ( large_scale_forcing  .OR.  land_surface  .OR.  urban_surface )  THEN
          IF ( ibc_pt_b == 0 )  THEN
             constant_heatflux = .FALSE.
          ELSEIF ( ibc_pt_b == 1 )  THEN
             constant_heatflux = .TRUE.
             surface_heatflux = 0.0_wp
          ENDIF
       ENDIF
    ELSE
       constant_heatflux = .TRUE.
    ENDIF

    IF ( top_heatflux == 9999999.9_wp )  constant_top_heatflux = .FALSE.

    IF ( neutral )  THEN

       IF ( surface_heatflux /= 0.0_wp  .AND.  surface_heatflux /= 9999999.9_wp )  THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PAC0085', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( top_heatflux /= 0.0_wp  .AND.  top_heatflux /= 9999999.9_wp )  THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PAC0085', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

    IF ( top_momentumflux_u /= 9999999.9_wp  .AND.  top_momentumflux_v /= 9999999.9_wp )  THEN
       constant_top_momentumflux = .TRUE.
    ELSEIF ( .NOT. ( top_momentumflux_u == 9999999.9_wp  .AND.                                     &
           top_momentumflux_v == 9999999.9_wp ) )  THEN
       message_string = 'both, top_momentumflux_u AND top_momentumflux_v must be set'
       CALL message( 'check_parameters', 'PAC0086', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given surface temperature implies Dirichlet boundary condition for temperature. In this case
!-- specification of a constant heat flux is forbidden.
!-- ATTENTION: the additional condition ibc_pt_b == 3 may be removed lateron, if it turns out that
!--            surface layer fluxes are not calculated for bottom boundary grid points that are
!--            located in the atmosphere!!!
    IF ( ( ibc_pt_b == 0  .OR.  ibc_pt_b == 3 )  .AND.                                             &
         constant_heatflux  .AND.  surface_heatflux /= 0.0_wp )                                    &
    THEN
       message_string = 'boundary_condition: bc_pt_b = "' // TRIM( bc_pt_b ) //                    &
                        '& is not allowed with constant_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PAC0087', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( constant_heatflux  .AND.  pt_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( message_string, * )  'constant_heatflux = .TRUE. is not allo',                      &
               'wed with pt_surface_initial_change (/=0) = ', pt_surface_initial_change
       CALL message( 'check_parameters', 'PAC0088', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( constant_heatflux  .AND.  pt_surface_heating_rate /= 0.0_wp )  THEN
       WRITE ( message_string, * )  'constant_heatflux = .TRUE. is not allo',                      &
               'wed with pt_surface_heating_rate (/=0) = ', pt_surface_heating_rate
       CALL message( 'check_parameters', 'PAC0089', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given temperature at the top implies Dirichlet boundary condition for temperature. In this
!-- case specification of a constant heat flux is forbidden.
    IF ( ibc_pt_t == 0  .AND.  constant_top_heatflux  .AND.  top_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_t = "' // TRIM( bc_pt_t ) //                    &
                        '" is not allowed with constant_top_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PAC0090', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set boundary conditions for total water content
    IF ( humidity )  THEN

       IF ( ANY( wall_humidityflux /= 0.0_wp )  .AND.  surface_waterflux == 9999999.9_wp )  THEN
          message_string = 'wall_humidityflux additionally requires setting of surface_waterflux'
          CALL message( 'check_parameters', 'PAC0091', 1, 2, 0, 6, 0 )
       ENDIF

       CALL set_bc_scalars( 'q', bc_q_b, bc_q_t, ibc_q_b, ibc_q_t, 'PAC0092', 'PAC0093' )

       IF ( surface_waterflux == 9999999.9_wp  )  THEN
          constant_waterflux = .FALSE.
          IF ( large_scale_forcing .OR. land_surface )  THEN
             IF ( ibc_q_b == 0 )  THEN
                constant_waterflux = .FALSE.
             ELSEIF ( ibc_q_b == 1 )  THEN
                constant_waterflux = .TRUE.
             ENDIF
          ENDIF
       ELSE
          constant_waterflux = .TRUE.
       ENDIF

       CALL check_bc_scalars( 'q', bc_q_b, ibc_q_b, 'PAC0094', 'PAC0095', constant_waterflux,      &
                              q_surface_initial_change )

    ENDIF

    IF ( passive_scalar )  THEN

       IF ( ANY( wall_scalarflux /= 0.0_wp )  .AND.  surface_scalarflux == 9999999.9_wp )  THEN
          message_string = 'wall_scalarflux additionally requires setting of surface_scalarflux'
          CALL message( 'check_parameters', 'PAC0096', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( surface_scalarflux == 9999999.9_wp )  constant_scalarflux = .FALSE.

       CALL set_bc_scalars( 's', bc_s_b, bc_s_t, ibc_s_b, ibc_s_t, 'PAC0092', 'PAC0093' )

       CALL check_bc_scalars( 's', bc_s_b, ibc_s_b, 'PAC0094', 'PAC0095', constant_scalarflux,     &
                              s_surface_initial_change )

       IF ( top_scalarflux == 9999999.9_wp )  constant_top_scalarflux = .FALSE.
!
!--    A fixed scalar concentration at the top implies Dirichlet boundary condition for scalar.
!--    Hence, in this case specification of a constant scalar flux is forbidden.
       IF ( ( ibc_s_t == 0 .OR. ibc_s_t == 2 )  .AND.  constant_top_scalarflux  .AND.              &
              top_scalarflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_s_t = "' // TRIM( bc_s_t ) //                   &
                           '" is not allowed with top_scalarflux /= 0.0'
          CALL message( 'check_parameters', 'PAC0097', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Boundary conditions for horizontal components of wind speed
    IF ( bc_uv_b == 'dirichlet' )  THEN
       ibc_uv_b = 0
    ELSEIF ( bc_uv_b == 'neumann' )  THEN
       ibc_uv_b = 1
       IF ( constant_flux_layer )  THEN
          message_string = 'boundary condition: bc_uv_b = "' // TRIM( bc_uv_b ) //                 &
                           '" is not allowed with constant_flux_layer = .TRUE.'
          CALL message( 'check_parameters', 'PAC0098', 1, 2, 0, 6, 0 )
       ENDIF
    ELSEIF ( bc_uv_b == 'nested' )  THEN
       ibc_uv_b = 3
    ELSE
       message_string = 'unknown boundary condition: bc_uv_b = "' // TRIM( bc_uv_b ) // '"'
       CALL message( 'check_parameters', 'PAC0099', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- In case of coupled atmosphere-ocean simulations u and v at the atmosphere's surface will be
!-- set to the u and v values of the ocean surface.
    IF ( atmosphere_run_coupled_to_ocean )  THEN
       ibc_uv_b = 2
    ENDIF

    IF ( ocean_run_coupled_to_atmosphere )  THEN
       bc_uv_t = 'neumann'
       ibc_uv_t = 1
    ELSE
       IF ( bc_uv_t == 'dirichlet' .OR. bc_uv_t == 'dirichlet_0' )  THEN
          ibc_uv_t = 0
          IF ( bc_uv_t == 'dirichlet_0' )  THEN
!
!--          Velocities for the initial u,v-profiles are set zero at the top in case of dirichlet_0
!--          conditions
             u_init(nzt+1)    = 0.0_wp
             v_init(nzt+1)    = 0.0_wp
          ENDIF
       ELSEIF ( bc_uv_t == 'neumann' )  THEN
          ibc_uv_t = 1
       ELSEIF ( bc_uv_t == 'nested'  .OR.  bc_uv_t == 'nesting_offline' )  THEN
          ibc_uv_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_uv_t = "' // TRIM( bc_uv_t ) // '"'
          CALL message( 'check_parameters', 'PAC0100', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Check the Rayleigh damping factor.
    IF ( rayleigh_damping_factor < 0.0_wp  .OR.  rayleigh_damping_factor > 1.0_wp )  THEN
       WRITE( message_string, * )  'rayleigh_damping_factor = ', rayleigh_damping_factor,       &
              ' out of range'
       CALL message( 'check_parameters', 'PAC0101', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( rayleigh_damping_height == -1.0_wp )  THEN
       IF ( .NOT. ocean_mode )  THEN
          rayleigh_damping_height = 0.66666666666_wp * zu(nzt)
       ELSE
          rayleigh_damping_height = 0.66666666666_wp * zu(nzb)
       ENDIF
    ELSE
       IF ( .NOT. ocean_mode )  THEN
          IF ( rayleigh_damping_height < 0.0_wp  .OR.  rayleigh_damping_height > zu(nzt) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ',  rayleigh_damping_height,   &
                    ' out of range [0.0,', zu(nzt), ']'
             CALL message( 'check_parameters', 'PAC0102', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          IF ( rayleigh_damping_height > 0.0_wp  .OR.  rayleigh_damping_height < zu(nzb) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ', rayleigh_damping_height,    &
                    ' out of range [0.0,', zu(nzb), ']'
             CALL message( 'check_parameters', 'PAC0102', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Check number of chosen statistic regions
    IF ( statistic_regions < 0 )  THEN
       WRITE ( message_string, * ) 'number of statistic_regions = ', statistic_regions+1,          &
               ' is not allowed'
       CALL message( 'check_parameters', 'PAC0103', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( normalizing_region > statistic_regions  .OR.  normalizing_region < 0)  THEN
       WRITE ( message_string, * ) 'normalizing_region = ', normalizing_region,                    &
               ' must be >= 0 and <= ',statistic_regions, ' (value of statistic_regions)'
       CALL message( 'check_parameters', 'PAC0104', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default intervals for data output, if necessary
!-- NOTE: dt_dosp has already been set in spectra_parin
    IF ( dt_data_output /= 9999999.9_wp )  THEN
       IF ( dt_dopr           == 9999999.9_wp )  dt_dopr           = dt_data_output
       IF ( dt_do2d_xy        == 9999999.9_wp )  dt_do2d_xy        = dt_data_output
       IF ( dt_do2d_xz        == 9999999.9_wp )  dt_do2d_xz        = dt_data_output
       IF ( dt_do2d_yz        == 9999999.9_wp )  dt_do2d_yz        = dt_data_output
       IF ( dt_do3d           == 9999999.9_wp )  dt_do3d           = dt_data_output
       IF ( dt_data_output_av == 9999999.9_wp )  dt_data_output_av = dt_data_output
       DO  mid = 1, max_masks
          IF ( dt_domask(mid) == 9999999.9_wp )  dt_domask(mid)    = dt_data_output
       ENDDO
    ENDIF

!
!-- Set the default skip time intervals for data output, if necessary
    IF ( skip_time_dopr    == 9999999.9_wp )  skip_time_dopr    = skip_time_data_output
    IF ( skip_time_do2d_xy == 9999999.9_wp )  skip_time_do2d_xy = skip_time_data_output
    IF ( skip_time_do2d_xz == 9999999.9_wp )  skip_time_do2d_xz = skip_time_data_output
    IF ( skip_time_do2d_yz == 9999999.9_wp )  skip_time_do2d_yz = skip_time_data_output
    IF ( skip_time_do3d    == 9999999.9_wp )  skip_time_do3d    = skip_time_data_output
    IF ( skip_time_data_output_av == 9999999.9_wp )                                                &
                                       skip_time_data_output_av = skip_time_data_output
    DO  mid = 1, max_masks
       IF ( skip_time_domask(mid) == 9999999.9_wp )                                                &
                                       skip_time_domask(mid)    = skip_time_data_output
    ENDDO

!
!-- Check the average intervals (first for 3d-data, then for profiles)
    IF ( averaging_interval > dt_data_output_av )  THEN
       WRITE( message_string, * )  'averaging_interval = ', averaging_interval,                    &
              ' must be <= dt_data_output_av = ', dt_data_output_av
       CALL message( 'check_parameters', 'PAC0105', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( averaging_interval_pr == 9999999.9_wp )  THEN
       averaging_interval_pr = averaging_interval
    ENDIF

    IF ( averaging_interval_pr > dt_dopr )  THEN
       WRITE( message_string, * )  'averaging_interval_pr = ', averaging_interval_pr,              &
              ' must be <= dt_dopr = ', dt_dopr
       CALL message( 'check_parameters', 'PAC0106', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default interval for profiles entering the temporal average
    IF ( dt_averaging_input_pr == 9999999.9_wp )  THEN
       dt_averaging_input_pr = dt_averaging_input
    ENDIF

!
!-- Set meaningful default flux output mode if nothing else is prescribed
    IF ( TRIM( flux_output_mode ) == 'application-specific' )  THEN
!
!--    Set flux output mode according to approximation if applicable
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_output_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_output_mode = 'kinematic'
       ENDIF
!
!--    When the land- or urban-surface model is used, the flux output must be dynamic
       IF ( land_surface  .OR.  urban_surface )  THEN
          flux_output_mode = 'dynamic'
       ENDIF
    ENDIF

!
!-- Set the flux output units according to flux_output_mode
    IF ( TRIM( flux_output_mode ) == 'kinematic' )  THEN
        heatflux_output_unit              = 'K m/s'
        waterflux_output_unit             = 'kg/kg m/s'
        momentumflux_output_unit          = 'm2/s2'
    ELSEIF ( TRIM( flux_output_mode ) == 'dynamic' )  THEN
        heatflux_output_unit              = 'W/m2'
        waterflux_output_unit             = 'W/m2'
        momentumflux_output_unit          = 'N/m2'
    ENDIF

!
!-- set time series output units for fluxes
    dots_unit(14:16) = TRIM( heatflux_output_unit )
    dots_unit(21)    = TRIM( waterflux_output_unit )
    dots_unit(19:20) = TRIM( momentumflux_output_unit )

!
!-- Set the default interval for the output of timeseries to a reasonable value (tries to minimize
!-- the number of calls of flow_statistics)
    IF ( dt_dots == 9999999.9_wp )  THEN
       IF ( averaging_interval_pr == 0.0_wp )  THEN
          dt_dots = MIN( dt_run_control, dt_dopr )
       ELSE
          dt_dots = MIN( dt_run_control, dt_averaging_input_pr )
       ENDIF
    ENDIF

!
!-- Check the sample rate for averaging (first for 3d-data, then for profiles).
!-- Default values for dt_averaging_input and averaging interval are set to zero. Hence, this check
!-- should be only applied for values of dt_averaging_input > 0.0. Else, it is not relevant.
    IF ( dt_averaging_input >= averaging_interval  .AND.  dt_averaging_input > 0.0_wp )  THEN
       WRITE( message_string, * )  'dt_averaging_input = ', dt_averaging_input,                    &
              ' must be < averaging_interval = ', averaging_interval
       CALL message( 'check_parameters', 'PAC0107', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( dt_averaging_input_pr >= averaging_interval_pr  .AND.  dt_averaging_input_pr > 0.0_wp )   &
    THEN
       WRITE( message_string, * )  'dt_averaging_input_pr = ', dt_averaging_input_pr,              &
              ' must be < averaging_interval_pr = ', averaging_interval_pr
       CALL message( 'check_parameters', 'PAC0108', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check data-output variable list for any duplicates and remove them
    ALLOCATE( data_output_filtered(1:UBOUND( data_output_pr, DIM=1 )) )
    CALL filter_duplicate_strings( varnamelength, data_output_pr, data_output_filtered )
    data_output_pr = data_output_filtered
    DEALLOCATE( data_output_filtered )
!
!-- Determine the number of output profiles and check whether they are permissible
    DO  WHILE ( data_output_pr(dopr_n+1) /= '          ' )

       dopr_n = dopr_n + 1
       i = dopr_n

!
!--    Size check. Some arrays in the netcdf interface module are internally limited to 500 elements
!--    and would need to be increased for output of more than 500 profiles.
       IF ( dopr_n > 500 )  THEN
          message_string = 'output of more than 500 profiles requested'
          CALL message( 'check_parameters', 'PAC0109', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Determine internal profile number (for hom, homs) and store height levels
       SELECT CASE ( TRIM( data_output_pr(i) ) )

          CASE ( 'u', '#u' )
             dopr_index(i) = 1
             dopr_unit(i)  = 'm/s'
             hom(:,2,1,:)  = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 5
                hom(:,2,5,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'v', '#v' )
             dopr_index(i) = 2
             dopr_unit(i)  = 'm/s'
             hom(:,2,2,:)  = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 6
                hom(:,2,6,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w' )
             dopr_index(i) = 3
             dopr_unit(i)  = 'm/s'
             hom(:,2,3,:)  = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'theta', '#theta' )
             IF ( .NOT. bulk_cloud_model )  THEN
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 7
                   hom(:,2,7,:)          = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ELSE
                dopr_index(i) = 43
                dopr_unit(i)  = 'K'
                hom(:,2,43,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 28
                   hom(:,2,28,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'e', '#e' )
             dopr_index(i)  = 8
             dopr_unit(i)   = 'm2/s2'
             hom(:,2,8,:)   = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 8
                hom(:,2,8,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'km', '#km' )
             dopr_index(i)  = 9
             dopr_unit(i)   = 'm2/s'
             hom(:,2,9,:)   = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 23
                hom(:,2,23,:)         = hom(:,2,9,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'kh', '#kh' )
             dopr_index(i)   = 10
             dopr_unit(i)    = 'm2/s'
             hom(:,2,10,:)   = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 24
                hom(:,2,24,:)         = hom(:,2,10,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'l', '#l' )
             dopr_index(i)   = 11
             dopr_unit(i)    = 'm'
             hom(:,2,11,:)   = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 25
                hom(:,2,25,:)         = hom(:,2,11,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w"u"' )
             dopr_index(i) = 12
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,12,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*u*' )
             dopr_index(i) = 13
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,13,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"v"' )
             dopr_index(i) = 14
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,14,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*v*' )
             dopr_index(i) = 15
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,15,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"theta"' )
             dopr_index(i) = 16
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,16,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*theta*' )
             dopr_index(i) = 17
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,17,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wtheta' )
             dopr_index(i) = 18
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,18,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wu' )
             dopr_index(i) = 19
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,19,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wv' )
             dopr_index(i) = 20
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,20,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*theta*BC' )
             dopr_index(i) = 21
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,21,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wthetaBC' )
             dopr_index(i) = 22
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,22,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'u*2' )
             dopr_index(i) = 30
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,30,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v*2' )
             dopr_index(i) = 31
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,31,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*2' )
             dopr_index(i) = 32
             dopr_unit(i)  = 'm2/s2'
             IF ( .NOT. ws_scheme_mom  )  THEN
                hom(:,2,32,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                hom(:,2,32,0) = zu
                If ( statistic_regions > 0 )  THEN
                   hom(:,2,32,1:statistic_regions) = SPREAD( zw, 2, statistic_regions )
                ENDIF
             ENDIF

          CASE ( 'theta*2' )
             dopr_index(i) = 33
             dopr_unit(i)  = 'K2'
             hom(:,2,33,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'e*' )
             dopr_index(i) = 34
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,34,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*2theta*' )
             dopr_index(i) = 35
             dopr_unit(i)  = 'K m2/s2'
             hom(:,2,35,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*theta*2' )
             dopr_index(i) = 36
             dopr_unit(i)  = 'K2 m/s'
             hom(:,2,36,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*e*' )
             dopr_index(i) = 37
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,37,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*3' )
             dopr_index(i) = 38
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,38,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'Sw' )
             dopr_index(i) = 39
             dopr_unit(i)  = 'none'
             hom(:,2,39,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'p' )
             dopr_index(i) = 40
             dopr_unit(i)  = 'Pa'
             hom(:,2,40,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'q', '#q' )
             IF ( .NOT. humidity )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0110', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 41
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,41,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 's', '#s' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PAC0111', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 115
                dopr_unit(i)  = 'kg/m3'
                hom(:,2,115,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 121
                   hom(:,2,121,:)        = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'qv', '#qv' )
             IF ( .NOT. bulk_cloud_model ) THEN
                dopr_index(i) = 41
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,41,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ELSE
                dopr_index(i) = 42
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,42,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 27
                   hom(:,2,27,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'thetal', '#thetal' )
             IF ( .NOT. bulk_cloud_model ) THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for bulk_cloud_model = .FALSE.'
                CALL message( 'check_parameters', 'PAC0112', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 7
                   hom(:,2,7,:)          = SPREAD( zu, 2, statistic_regions+1 )
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'thetav', '#thetav' )
             dopr_index(i) = 44
             dopr_unit(i)  = 'K'
             hom(:,2,44,:) = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 29
                hom(:,2,29,:)         = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w"thetav"' )
             dopr_index(i) = 45
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,45,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*thetav*' )
             dopr_index(i) = 46
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,46,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wthetav' )
             dopr_index(i) = 47
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,47,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"q"' )
             IF ( .NOT. humidity )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0110', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 48
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,48,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*q*' )
             IF ( .NOT. humidity )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0110', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 49
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,49,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'wq' )
             IF ( .NOT. humidity )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0110', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 50
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,50,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"s"' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PAC0111', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 117
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,117,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*s*' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PAC0111', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 114
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,114,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'ws' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PAC0111', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 118
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,118,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"qv"' )
             IF ( humidity  .AND.  .NOT. bulk_cloud_model )  THEN
                dopr_index(i) = 48
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,48,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF ( humidity  .AND.  bulk_cloud_model )  THEN
                dopr_index(i) = 51
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,51,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for bulk_cloud_model = .FALSE. ' //          &
                                 'and humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0113', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'w*qv*' )
             IF ( humidity  .AND.  .NOT. bulk_cloud_model )  THEN
                dopr_index(i) = 49
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,49,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF( humidity  .AND.  bulk_cloud_model )  THEN
                dopr_index(i) = 52
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,52,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for bulk_cloud_model = .FALSE. ' //          &
                                 'and humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0113', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'wqv' )
             IF ( humidity  .AND.  .NOT. bulk_cloud_model )  THEN
                dopr_index(i) = 50
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,50,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF ( humidity  .AND.  bulk_cloud_model )  THEN
                dopr_index(i) = 53
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,53,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for bulk_cloud_model = .FALSE. ' //          &
                                 'and humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0113', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'ql' )
             IF ( .NOT. bulk_cloud_model  .AND.  .NOT. cloud_droplets )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for bulk_cloud_model = .FALSE. ' //          &
                                 'and cloud_droplets = .FALSE.'
                CALL message( 'check_parameters', 'PAC0114', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 54
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,54,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*u*u*ddz' )
             dopr_index(i) = 55
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,55,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*p*ddz' )
             dopr_index(i) = 56
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,56,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"eddz' )
             dopr_index(i) = 57
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,57,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'u"theta"' )
             dopr_index(i) = 58
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,58,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'u*theta*' )
             dopr_index(i) = 59
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,59,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'utheta_t' )
             dopr_index(i) = 60
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,60,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v"theta"' )
             dopr_index(i) = 61
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,61,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v*theta*' )
             dopr_index(i) = 62
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,62,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'vtheta_t' )
             dopr_index(i) = 63
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,63,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*p*' )
             dopr_index(i) = 68
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,68,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w"e' )
             dopr_index(i) = 69
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,69,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'q*2' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PAC0110', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 70
                dopr_unit(i)  = 'kg2/kg2'
                hom(:,2,70,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'hyp' )
             dopr_index(i) = 72
             dopr_unit(i)  = 'hPa'
             hom(:,2,72,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho' )
             dopr_index(i)  = 119
             dopr_unit(i)   = 'kg/m3'
             hom(:,2,119,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho_zw' )
             dopr_index(i)  = 120
             dopr_unit(i)   = 'kg/m3'
             hom(:,2,120,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'ug' )
             dopr_index(i) = 78
             dopr_unit(i)  = 'm/s'
             hom(:,2,78,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'vg' )
             dopr_index(i) = 79
             dopr_unit(i)  = 'm/s'
             hom(:,2,79,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w_subs' )
             IF (  .NOT.  large_scale_subsidence )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for large_scale_subsidence = .FALSE.'
                CALL message( 'check_parameters', 'PAC0115', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 80
                dopr_unit(i)  = 'm/s'
                hom(:,2,80,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 's*2' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(i) ) //               &
                                 ' is not implemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PAC0116', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 116
                dopr_unit(i)  = 'kg2/m6'
                hom(:,2,116,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

           CASE ( 'eta' )
              dopr_index(i) = 121
              dopr_unit(i)  = 'm'
              hom(:,2,121,:) = SPREAD( zu, 2, statistic_regions+1 )

              kolmogorov_length_scale = .TRUE.

          CASE DEFAULT
             unit = 'illegal'
!
!--          Check for other modules
             CALL module_interface_check_data_output_pr( data_output_pr(i), i, unit, dopr_unit(i) )

!
!--          No valid quantity found
             IF ( unit == 'illegal' )  THEN
                IF ( data_output_pr_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output_pr or ' //                      &
                                    'data_output_pr_user = "' // TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PAC0117', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output_pr = "' //                      &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PAC0118', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT

    ENDDO


!
!-- Append user-defined data output variables to the standard data output
    IF ( data_output_user(1) /= ' ' )  THEN
       i = 1
       DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
          i = i + 1
       ENDDO
       j = 1
       DO  WHILE ( data_output_user(j) /= ' '  .AND.  j <= 500 )
          IF ( i > 500 )  THEN
             message_string = 'number of output quantitities given by data' //                     &
                              '_output and data_output_user exceeds the limit of 500'
             CALL message( 'check_parameters', 'PAC0119', 1, 2, 0, 6, 0 )
          ENDIF
          data_output(i) = data_output_user(j)
          i = i + 1
          j = j + 1
       ENDDO
    ENDIF
!
!-- Check data-output variable list for any duplicates and remove them
    ALLOCATE( data_output_filtered(1:UBOUND( data_output, DIM=1 )) )
    CALL filter_duplicate_strings( varnamelength, data_output, data_output_filtered )
    data_output = data_output_filtered
    DEALLOCATE( data_output_filtered )
!
!-- Check and set steering parameters for 2d/3d data output and averaging
    i   = 1
    DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
!
!--    Check for data averaging
       ilen = LEN_TRIM( data_output(i) )
       j = 0                                                 ! no data averaging
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_av' )  THEN
             j = 1                                           ! data averaging
             data_output(i) = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF
!
!--    Check for cross section or volume data
       ilen = LEN_TRIM( data_output(i) )
       k = 0                                                   ! 3d data
       var = data_output(i)(1:ilen)
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_xy'  .OR.  data_output(i)(ilen-2:ilen) == '_xz'    &
               .OR.  data_output(i)(ilen-2:ilen) == '_yz' )  THEN
             k = 1                                             ! 2d data
             var = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF

!
!--    Check for allowed value and set units
       SELECT CASE ( TRIM( var ) )

          CASE ( 'e' )
             IF ( constant_diffusion )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'constant_diffusion = .FALSE.'
                CALL message( 'check_parameters', 'PAC0120', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'm2/s2'

          CASE ( 'thetal' )
             IF ( .NOT. bulk_cloud_model )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'bulk_cloud_model = .TRUE.'
                CALL message( 'check_parameters', 'PAC0121', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'K'

          CASE ( 'pc', 'pr' )
             IF ( .NOT. particle_advection )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'a "particle_parameters" namelist in the parameter file (PARIN)'
                CALL message( 'check_parameters', 'PAC0122', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'pc' )  unit = 'number'
             IF ( TRIM( var ) == 'pr' )  unit = 'm'

          CASE ( 'q', 'thetav' )
             IF ( .NOT. humidity )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires humidity = .TRUE.'
                CALL message( 'check_parameters', 'PAC0123', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'q'   )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'thetav' )  unit = 'K'

          CASE ( 'ql' )
             IF ( .NOT.  ( bulk_cloud_model  .OR.  cloud_droplets ) )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'bulk_cloud_model = .TRUE. or cloud_droplets = .TRUE.'
                CALL message( 'check_parameters', 'PAC0124', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'ql_c', 'ql_v', 'ql_vp' )
             IF ( .NOT. cloud_droplets )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'cloud_droplets = .TRUE.'
                CALL message( 'check_parameters', 'PAC0125', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'ql_c'  )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'ql_v'  )  unit = 'm3'
             IF ( TRIM( var ) == 'ql_vp' )  unit = 'none'

          CASE ( 'qv' )
             IF ( .NOT. bulk_cloud_model )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'bulk_cloud_model = .TRUE.'
                CALL message( 'check_parameters', 'PAC0121', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 's' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PAC0126', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/m3'

          CASE ( 'p', 'theta', 'u', 'v', 'w' )
             IF ( TRIM( var ) == 'p'  )  unit = 'Pa'
             IF ( TRIM( var ) == 'theta' )  unit = 'K'
             IF ( TRIM( var ) == 'u'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'v'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'w'  )  unit = 'm/s'
             CONTINUE

          CASE ( 'ghf*', 'lwp*', 'ol*', 'pres_drag_x*', 'pres_drag_y*', 'qsurf*', 'qsws*', 'r_a*', &
                 'shf*', 'ssurf*', 'ssws*', 't*', 'tsurf*', 'us*', 'z0*', 'z0h*', 'z0q*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' // TRIM( var ) //              &
                                 '" & only 2d-horizontal cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PAC0127', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'lwp*'  .AND.  .NOT. bulk_cloud_model )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'bulk_cloud_model = .TRUE.'
                CALL message( 'check_parameters', 'PAC0121', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'qsws*'  .AND.  .NOT.  humidity )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires humidity = .TRUE.'
                CALL message( 'check_parameters', 'PAC0128', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'ghf*'  .AND.  .NOT.  land_surface )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requires land_surface = .TRUE.'
                CALL message( 'check_parameters', 'PAC0129', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( (  TRIM( var ) == 'r_a*' .OR. TRIM( var ) == 'ghf*' )  .AND.  .NOT. land_surface &
                   .AND.  .NOT. urban_surface )                                                    &
             THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'land_surface = .TRUE. or ' // 'urban_surface = .TRUE.'
                CALL message( 'check_parameters', 'PAC0130', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( ( TRIM( var ) == 'ssws*' .OR.  TRIM( var ) == 'ssurf*' ) .AND.                   &
                   .NOT. passive_scalar )                                                          &
             THEN
                message_string = 'output of "' // TRIM( var ) // '" requires ' //                  &
                                 'passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PAC0131', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'ghf*'         )  unit = 'W/m2'
             IF ( TRIM( var ) == 'lwp*'         )  unit = 'kg/m2'
             IF ( TRIM( var ) == 'ol*'          )  unit = 'm'
             IF ( TRIM( var ) == 'pres_drag_x*' )  unit = 'kgm/s2'
             IF ( TRIM( var ) == 'pres_drag_y*' )  unit = 'kgm/s2'
             IF ( TRIM( var ) == 'qsurf*'       )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'qsws*'        )  unit = TRIM ( waterflux_output_unit )
             IF ( TRIM( var ) == 'r_a*'         )  unit = 's/m'
             IF ( TRIM( var ) == 'shf*'         )  unit = TRIM ( heatflux_output_unit )
             IF ( TRIM( var ) == 'ssurf*'       )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'ssws*'        )  unit = 'kg/m2*s'
             IF ( TRIM( var ) == 't*'           )  unit = 'K'
             IF ( TRIM( var ) == 'tsurf*'       )  unit = 'K'
             IF ( TRIM( var ) == 'us*'          )  unit = 'm/s'
             IF ( TRIM( var ) == 'z0*'          )  unit = 'm'
             IF ( TRIM( var ) == 'z0h*'         )  unit = 'm'

          CASE DEFAULT
!
!--          Check for other modules
             CALL module_interface_check_data_output( var, unit, i, j, ilen, k )

             IF ( unit == 'illegal' )  THEN
                IF ( data_output_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output or ' //                         &
                                    'data_output_user = "' // TRIM( data_output(i) ) // '"'
                   CALL message( 'check_parameters', 'PAC0132', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output = "' //                         &
                                    TRIM( data_output(i) ) // '"'
                   CALL message( 'check_parameters', 'PAC0133', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT
!
!--    Set the internal steering parameters appropriately
       IF ( k == 0 )  THEN
          do3d_no(j)              = do3d_no(j) + 1
          do3d(j,do3d_no(j))      = data_output(i)
          do3d_unit(j,do3d_no(j)) = unit
       ELSE
          do2d_no(j)              = do2d_no(j) + 1
          do2d(j,do2d_no(j))      = data_output(i)
          do2d_unit(j,do2d_no(j)) = unit
          IF ( data_output(i)(ilen-2:ilen) == '_xy' )  THEN
             data_output_xy(j) = .TRUE.
          ENDIF
          IF ( data_output(i)(ilen-2:ilen) == '_xz' )  THEN
             data_output_xz(j) = .TRUE.
          ENDIF
          IF ( data_output(i)(ilen-2:ilen) == '_yz' )  THEN
             data_output_yz(j) = .TRUE.
          ENDIF
       ENDIF

       IF ( j == 1 )  THEN
!
!--       Check, if variable is already subject to averaging
          found = .FALSE.
          DO  k = 1, doav_n
             IF ( TRIM( doav(k) ) == TRIM( var ) )  found = .TRUE.
          ENDDO

          IF ( .NOT. found )  THEN
             doav_n = doav_n + 1
             doav(doav_n) = var
          ENDIF
       ENDIF

       i = i + 1
    ENDDO

!
!-- Averaged 2d or 3d output requires that an averaging interval has been set
    IF ( doav_n > 0  .AND.  averaging_interval == 0.0_wp )  THEN
       WRITE( message_string, * )  'output of averaged quantity "', TRIM( doav(1) ),               &
                                   '_av" requires to set a ', 'non-zero averaging interval'
       CALL message( 'check_parameters', 'PAC0134', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check sectional planes given in physical coordinates, calculate the respective grid index,
!-- and add the index to the sectional planes that are given as indices.
    IF ( section_xy_m(1) /= -9999999.9_wp )  THEN
!
!--    First determine, how many sections are given as indices.
       next_index = COUNT( section_xy /= -9999 ) + 1

       k = 1
       DO  WHILE ( section_xy_m(k) /= -9999999.9_wp  .AND.  k < 100 )
!
!--       Check, if the section lies within the domain limits:
!--       Calculate the relative height with respect to the bottom domain boundary.
          section_position = section_xy_m(k) - init_model%origin_z
!
!--       Check position and abort, if it lies outside the domain.
          IF ( section_position < 0.0_wp  .OR.  section_position > zu(nz+1) )  THEN
             WRITE( message_string, '(A,F8.2,A,F8.2,A)' )  'section_xy_m must be >= ',             &
                                                           init_model%origin_z, 'm and <= ',       &
                                                           init_model%origin_z + zu(nz+1), ' m'
             CALL message( 'check_parameters', 'PAC0135', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       For the given position, calculate the respective index and add it to the section index
!--       array.
          DO  kk = 0, nz+1
             IF ( section_position <= zu(kk) )  THEN
                IF ( next_index <= 100 )  THEN
                   section_xy(next_index) = kk
                   next_index = next_index + 1
                ENDIF
                EXIT
             ENDIF
          ENDDO

          k = k + 1

       ENDDO
!
!--    Check for multiple occurrences of output layers and remove them if required.
       ALLOCATE( section_tmp(1:COUNT( section_xy /= -9999 )) )
       section_tmp = HUGE( 1 )

       kk = 1
       DO  k = 1, COUNT( section_xy /= -9999 )
          IF ( .NOT. ANY( section_xy(k) == section_tmp ) )  THEN
             section_tmp(kk) = section_xy(k)
             kk = kk + 1
          ENDIF
       ENDDO
       nr_unique_sections = COUNT( section_tmp /= HUGE( 1 ) )

       section_xy = -9999
       section_xy(1:nr_unique_sections) = section_tmp(1:nr_unique_sections)
       DEALLOCATE( section_tmp )

    ENDIF

!
!-- Same for xz-cross-sections given in physical coordinates.
    IF ( section_xz_m(1) /= -9999999.9_wp )  THEN
!
!--    First determine, how many sections are given as indices.
       next_index = COUNT( section_xz /= -9999 ) + 1

       j = 1
       DO  WHILE ( section_xz_m(j) /= -9999999.9_wp  .AND.  j < 100 )
!
!--       Check, if the section lies within the domain limits:
!--       Calculate the relative distance with respect to the left domain boundary.
          section_position = section_xz_m(j) - lower_left_coord_y
!
!--       Check position and abort, if it lies outside the domain.
          IF ( section_position < 0.0_wp  .OR.  section_position > ( ny * dy ) )  THEN
             WRITE( message_string, '(A,F8.2,A,F8.2,A)' )  'section_xz_m must be >= ',             &
                                                           lower_left_coord_y, 'm and <= ',        &
                                                           lower_left_coord_y + ( ny * dy ), ' m'
             CALL message( 'check_parameters', 'PAC0136', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       For the given position, calculate the respective index and add it to the section index
!--       array.
          section_xz(next_index) = NINT( section_position / dy )
          next_index = next_index + 1

          j = j + 1

       ENDDO
!
!--    Check for multiple occurrences of output layers and remove them if required.
       ALLOCATE( section_tmp(1:COUNT( section_xz /= -9999 )) )
       section_tmp = HUGE( 1 )

       kk = 1
       DO  k = 1, COUNT( section_xz /= -9999 )
          IF ( .NOT. ANY( section_xz(k) == section_tmp ) )  THEN
             section_tmp(kk) = section_xz(k)
             kk = kk + 1
          ENDIF
       ENDDO
       nr_unique_sections = COUNT( section_tmp /= HUGE( 1 ) )

       section_xz = -9999
       section_xz(1:nr_unique_sections) = section_tmp(1:nr_unique_sections)
       DEALLOCATE( section_tmp )

    ENDIF

!
!-- Same for yz-cross-sections given in physical coordinates.
    IF ( section_yz_m(1) /= -9999999.9_wp )  THEN
!
!--    First determine, how many sections are given as indices.
       next_index = COUNT( section_yz /= -9999 ) + 1

       i = 1
       DO  WHILE ( section_yz_m(i) /= -9999999.9_wp  .AND.  i < 100 )
!
!--       Check, if the section lies within the domain limits:
!--       Calculate the relative distance with respect to the left domain boundary.
          section_position = section_yz_m(i) - lower_left_coord_x
!
!--       Check position and abort, if it lies outside the domain.
          IF ( section_position < 0.0_wp  .OR.  section_position > ( nx * dx ) )  THEN
             WRITE( message_string, '(A,F8.2,A,F8.2,A)' )  'section_yz_m must be >= ',             &
                                                           lower_left_coord_x, 'm and <= ',        &
                                                           lower_left_coord_x + ( nx * dx ), ' m'
             CALL message( 'check_parameters', 'PAC0137', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       For the given position, calculate the respective index and add it to the section index
!--       array.
          section_yz(next_index) = NINT( section_position / dx )
          next_index = next_index + 1

          i = i + 1

       ENDDO
!
!--    Check for multiple occurrences of output layers and remove them if required.
       ALLOCATE( section_tmp(1:COUNT( section_yz /= -9999 )) )
       section_tmp = HUGE( 1 )

       kk = 1
       DO  k = 1, COUNT( section_yz /= -9999 )
          IF ( .NOT. ANY( section_yz(k) == section_tmp ) )  THEN
             section_tmp(kk) = section_yz(k)
             kk = kk + 1
          ENDIF
       ENDDO
       nr_unique_sections = COUNT( section_tmp /= HUGE( 1 ) )

       section_yz = -9999
       section_yz(1:nr_unique_sections) = section_tmp(1:nr_unique_sections)
       DEALLOCATE( section_tmp )

    ENDIF

!
!-- Check sectional planes indices and store them in one shared array.
    IF ( ANY( section_xy > nz + 1 ) )  THEN
       WRITE( message_string, * )  'section_xy must be <= nz + 1 = ', nz + 1
       CALL message( 'check_parameters', 'PAC0138', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_xz > ny ) )  THEN
       WRITE( message_string, * )  'section_xz must be <= ny = ', ny
       CALL message( 'check_parameters', 'PAC0139', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_yz > nx ) )  THEN
       WRITE( message_string, * )  'section_yz must be <= nx = ', nx
       CALL message( 'check_parameters', 'PAC0140', 1, 2, 0, 6, 0 )
    ENDIF
    section(:,1) = section_xy
    section(:,2) = section_xz
    section(:,3) = section_yz

    IF ( ANY( data_output_xy )  .AND.  .NOT. ANY( section(:,1) /= -9999 ) )  THEN
       WRITE( message_string, * )  'section_xy not defined for requested xy-cross section output'
       CALL message( 'check_parameters', 'PAC0141', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( data_output_xz )  .AND.  .NOT. ANY( section(:,2) /= -9999 ) )  THEN
       WRITE( message_string, * )  'section_xz not defined for requested xz-cross section output'
       CALL message( 'check_parameters', 'PAC0141', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( data_output_yz )  .AND.  .NOT. ANY( section(:,3) /= -9999 ) )  THEN
       WRITE( message_string, * )  'section_yz not defined for requested yz-cross section output'
       CALL message( 'check_parameters', 'PAC0141', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Upper plot limit for 3D arrays
    IF ( nz_do3d == -9999 )  nz_do3d = nzt + 1

!
!-- Set output format string (used in header)
    SELECT CASE ( netcdf_data_format )
       CASE ( 1 )
          netcdf_data_format_string = 'netCDF classic'
       CASE ( 2 )
          netcdf_data_format_string = 'netCDF 64bit offset'
       CASE ( 3 )
          netcdf_data_format_string = 'netCDF4/HDF5'
       CASE ( 4 )
          netcdf_data_format_string = 'netCDF4/HDF5 classic'
       CASE ( 5 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5'
       CASE ( 6 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5 classic'

    END SELECT

!
!-- Check mask conditions
    DO mid = 1, max_masks
       IF ( data_output_masks(mid,1) /= ' '  .OR.  data_output_masks_user(mid,1) /= ' ' )  THEN
          masks = masks + 1
       ENDIF
    ENDDO

    IF ( masks < 0  .OR.  masks > max_masks )  THEN
       WRITE( message_string, * )  'illegal value: masks must be >= 0 and ', '<= ', max_masks,     &
              ' (=max_masks)'
       CALL message( 'check_parameters', 'PAC0142', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( masks > 0 )  THEN
       mask_scale(1) = mask_scale_x
       mask_scale(2) = mask_scale_y
       mask_scale(3) = mask_scale_z
       IF ( ANY( mask_scale <= 0.0_wp ) )  THEN
          WRITE( message_string, * )  'illegal value: mask_scale_x, mask_scale_y and mask_scale_z',&
                 'must be > 0.0'
          CALL message( 'check_parameters', 'PAC0143', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Generate masks for masked data output.
!--    Parallel netcdf output is not tested so far for masked data, hence netcdf_data_format is
!--    switched back to non-parallel output.
       netcdf_data_format_save = netcdf_data_format
       IF ( netcdf_data_format > 4 )  THEN
          IF ( netcdf_data_format == 5 ) netcdf_data_format = 3
          IF ( netcdf_data_format == 6 ) netcdf_data_format = 4
          message_string = 'netCDF file formats '// '5 (parallel netCDF 4) and 6 (parallel ' //    &
                           'netCDF 4 Classic model) & are currently not supported (not yet ' //    &
                           'tested) for masked data. &Using respective non-parallel' //            &
                           ' output for masked data.'
          CALL message( 'check_parameters', 'PAC0144', 0, 0, 0, 6, 0 )
       ENDIF
       CALL init_masks
       netcdf_data_format = netcdf_data_format_save
    ENDIF

!
!-- Check the NetCDF data format
    IF ( netcdf_data_format > 2 )  THEN
#if defined( __netcdf4 )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 format requested but no ' //                              &
                        'cpp-directive __netcdf4 given & switch back to 64-bit offset format'
       CALL message( 'check_parameters', 'PAC0145', 0, 1, 0, 6, 0 )
       netcdf_data_format = 2
#endif
    ENDIF
    IF ( netcdf_data_format > 4 )  THEN
#if defined( __netcdf4 ) && defined( __netcdf4_parallel )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 parallel output requested but no ' //                     &
                        'cpp-directive __netcdf4_parallel given & switch ' //                      &
                        'back to netCDF4 non-parallel output'
       CALL message( 'check_parameters', 'PAC0146', 0, 1, 0, 6, 0 )
       netcdf_data_format = netcdf_data_format - 2
#endif
    ENDIF

!
!-- Calculate fixed number of output time levels for parallel netcdf output.
!-- The time dimension has to be defined as limited for parallel output, because otherwise the I/O
!-- performance drops significantly.
    IF ( netcdf_data_format > 4 )  THEN

!
!--    Check if any of the follwoing data output interval is 0.0s, which is not allowed for parallel
!--    output.
       CALL check_dt_do( dt_do3d,           'dt_do3d'           )
       CALL check_dt_do( dt_do2d_xy,        'dt_do2d_xy'        )
       CALL check_dt_do( dt_do2d_xz,        'dt_do2d_xz'        )
       CALL check_dt_do( dt_do2d_yz,        'dt_do2d_yz'        )
       CALL check_dt_do( dt_data_output_av, 'dt_data_output_av' )

!
!--    Set needed time levels (ntdim) to saved time levels + to be saved time levels.
       ntdim_3d(0) = CEILING(                                                                      &
                     ( MIN( end_time, time_restart ) - MAX(                                        &
                        MERGE( skip_time_do3d, skip_time_do3d + spinup_time,                       &
                               data_output_during_spinup ),                                        &
                                                            simulated_time_at_begin )              &
                     ) / dt_do3d )
       IF ( do3d_at_begin ) ntdim_3d(0) = ntdim_3d(0) + 1

       ntdim_3d(1) = CEILING(                                                                      &
                     ( MIN( end_time, time_restart ) - MAX(                                        &
                        MERGE( skip_time_data_output_av, skip_time_data_output_av + spinup_time,   &
                               data_output_during_spinup ),                                        &
                                                            simulated_time_at_begin )              &
                     ) / dt_data_output_av )

       ntdim_2d_xy(0) = CEILING(                                                                   &
                        ( MIN( end_time, time_restart ) - MAX(                                     &
                           MERGE( skip_time_do2d_xy, skip_time_do2d_xy + spinup_time,              &
                                  data_output_during_spinup ),                                     &
                                                               simulated_time_at_begin )           &
                        ) / dt_do2d_xy )

       ntdim_2d_xz(0) = CEILING(                                                                   &
                        ( MIN( end_time, time_restart ) - MAX(                                     &
                           MERGE( skip_time_do2d_xz, skip_time_do2d_xz + spinup_time,              &
                                  data_output_during_spinup ),                                     &
                                                               simulated_time_at_begin )           &
                        ) / dt_do2d_xz )

       ntdim_2d_yz(0) = CEILING(                                                                   &
                        ( MIN( end_time, time_restart ) - MAX(                                     &
                           MERGE( skip_time_do2d_yz, skip_time_do2d_yz + spinup_time,              &
                                  data_output_during_spinup ),                                     &
                                                               simulated_time_at_begin )           &
                        ) / dt_do2d_yz )

       IF ( do2d_at_begin )  THEN
          ntdim_2d_xy(0) = ntdim_2d_xy(0) + 1
          ntdim_2d_xz(0) = ntdim_2d_xz(0) + 1
          ntdim_2d_yz(0) = ntdim_2d_yz(0) + 1
       ENDIF
!
!--    Please note, for averaged 2D data skip_time_data_output_av is the relavant output control
!--    parameter.
       ntdim_2d_xy(1) = CEILING(                                                                   &
                        ( MIN( end_time, time_restart ) - MAX(                                     &
                           MERGE( skip_time_data_output_av, skip_time_data_output_av + spinup_time,&
                                  data_output_during_spinup ),                                     &
                                                               simulated_time_at_begin )           &
                        ) / dt_data_output_av )

       ntdim_2d_xz(1) = CEILING(                                                                   &
                        ( MIN( end_time, time_restart ) - MAX(                                     &
                           MERGE( skip_time_data_output_av, skip_time_data_output_av + spinup_time,&
                                  data_output_during_spinup ),                                     &
                                                               simulated_time_at_begin )           &
                        ) / dt_data_output_av )

       ntdim_2d_yz(1) = CEILING(                                                                   &
                        ( MIN( end_time, time_restart ) - MAX(                                     &
                           MERGE( skip_time_data_output_av, skip_time_data_output_av + spinup_time,&
                                  data_output_during_spinup ),                                     &
                                                               simulated_time_at_begin )           &
                        ) / dt_data_output_av )

    ENDIF

!
!-- Check, whether a constant diffusion coefficient shall be used
    IF ( km_constant /= -1.0_wp )  THEN
       IF ( km_constant < 0.0_wp )  THEN
          WRITE( message_string, * )  'km_constant = ', km_constant, ' < 0.0'
          CALL message( 'check_parameters', 'PAC0147', 1, 2, 0, 6, 0 )
       ELSE
          IF ( prandtl_number < 0.0_wp )  THEN
             WRITE( message_string, * )  'prandtl_number = ', prandtl_number, ' < 0.0'
             CALL message( 'check_parameters', 'PAC0148', 1, 2, 0, 6, 0 )
          ENDIF
          constant_diffusion = .TRUE.

          IF ( constant_flux_layer )  THEN
             message_string = 'constant_flux_layer is not allowed with fixed value of km'
             CALL message( 'check_parameters', 'PAC0149', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- In case of non-cyclic lateral boundaries and a damping layer for the potential temperature,
!-- check the width of the damping layer
    IF ( bc_lr /= 'cyclic' ) THEN
       IF ( pt_damping_width < 0.0_wp  .OR.  pt_damping_width > REAL( (nx+1) * dx ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PAC0150', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_ns /= 'cyclic' )  THEN
       IF ( pt_damping_width < 0.0_wp  .OR.  pt_damping_width > REAL( (ny+1) * dy ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PAC0150', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- In case of non-cyclic lateral boundaries and an outflow damping layer,
!-- check the width of the damping layer
    IF ( bc_lr /= 'cyclic' )  THEN
       IF ( outflow_damping_width < 0.0_wp  .OR.  outflow_damping_width > REAL( (nx+1) * dx ) )    &
       THEN
          message_string = 'outflow_damping_width out of range'
          CALL message( 'check_parameters', 'PAC0351', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_ns /= 'cyclic' )  THEN
       IF ( outflow_damping_width < 0.0_wp  .OR.  outflow_damping_width > REAL( (ny+1) * dy ) )    &
       THEN
          message_string = 'outflow_damping_width out of range'
          CALL message( 'check_parameters', 'PAC0351', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- In case of an existing outflow damping layer, check whether bc_lr or bc_ns are non-cyclic and 
!-- whether topography height is zero at all inflow and outflow boundary points.
    IF ( outflow_damping_factor > 0.0_wp )  THEN
       IF ( bc_lr == 'dirichlet/radiation'  .OR.  bc_lr == 'radiation/dirichlet' )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
             DO  j = nys, nyn
                IF ( .NOT. BTEST( topo_flags(nzb + 1,j,nxl), 0) )  THEN
                   found_err_l = .TRUE.
                   EXIT
                ENDIF
             ENDDO
          ELSEIF ( bc_radiation_r  .OR.  bc_dirichlet_r )  THEN
             DO  j = nys, nyn
                IF ( .NOT. BTEST( topo_flags(nzb + 1,j,nxr), 0) )  THEN
                   found_err_l = .TRUE.
                   EXIT
                ENDIF
             ENDDO
          ENDIF
       ELSEIF ( bc_ns == 'dirichlet/radiation'  .OR.  bc_ns == 'radiation/dirichlet' )  THEN
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
             DO  i = nxl, nxr
                IF ( .NOT. BTEST( topo_flags(nzb + 1,nys,i), 0) )  THEN
                   found_err_l = .TRUE.
                   EXIT
                ENDIF
             ENDDO
          ELSEIF ( bc_radiation_n  .OR.  bc_dirichlet_n )  THEN
             DO  i = nxl, nxr
                IF ( .NOT. BTEST( topo_flags(nzb + 1,nyn,i), 0) )  THEN
                   found_err_l = .TRUE.
                   EXIT
                ENDIF
             ENDDO
          ENDIF
      ELSE
         message_string = 'given lateral boundary condition(s) not allowed for outflow damping'
         CALL message( 'check_parameters', 'PAC0352', 1, 2, 0, 6, 0 ) 
       ENDIF
    ENDIF

#if defined( __parallel )
    CALL MPI_ALLREDUCE( found_err_l, found_err, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr  )
#else
    found_err = found_err_l
#endif
    IF ( found_err )  THEN
       message_string = 'topography height is not zero at the inflow and/or outflow boundary'
       CALL message( 'check_parameters', 'PAC0353', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check value range for zeta = z/L
    IF ( zeta_min >= zeta_max )  THEN
       WRITE( message_string, * )  'zeta_min = ', zeta_min, ' must be less ', 'than zeta_max = ',  &
              zeta_max
       CALL message( 'check_parameters', 'PAC0151', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check random generator
    IF ( (random_generator /= 'system-specific'      .AND.                                         &
          random_generator /= 'random-parallel'   )  .AND.                                         &
          random_generator /= 'numerical-recipes' )  THEN
       message_string = 'unknown random generator: random_generator = "' //                        &
                        TRIM( random_generator ) // '"'
       CALL message( 'check_parameters', 'PAC0152', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine upper and lower hight level indices for random perturbations
    IF ( disturbance_level_b == -9999999.9_wp )  THEN
       IF ( ocean_mode )  THEN
          disturbance_level_b     = zu((nzt*2)/3)
          disturbance_level_ind_b = ( nzt * 2 ) / 3
       ELSE
          disturbance_level_b     = zu(nzb+3)
          disturbance_level_ind_b = nzb + 3
       ENDIF
    ELSEIF ( disturbance_level_b < zu(3) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ', disturbance_level_b, ' must be >= ',  &
              zu(3), '(zu(3))'
       CALL message( 'check_parameters', 'PAC0153', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_b > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ', disturbance_level_b, ' must be <= ',  &
              zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PAC0154', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_b <= zu(k) )  THEN
             disturbance_level_ind_b = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( disturbance_level_t == -9999999.9_wp )  THEN
       IF ( ocean_mode )  THEN
          disturbance_level_t     = zu(nzt-3)
          disturbance_level_ind_t = nzt - 3
       ELSE
          disturbance_level_t     = zu(nzt/3)
          disturbance_level_ind_t = nzt / 3
       ENDIF
    ELSEIF ( disturbance_level_t > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ', disturbance_level_t, ' must be <= ',  &
              zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PAC0155', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_t < disturbance_level_b )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ', disturbance_level_t,                  &
             ' must be >= disturbance_level_b = ', disturbance_level_b
       CALL message( 'check_parameters', 'PAC0156', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_t <= zu(k) )  THEN
             disturbance_level_ind_t = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

!
!-- Check again whether the levels determined this way are ok.
!-- Error may occur at automatic determination and too few grid points in z-direction.
    IF ( disturbance_level_ind_t < disturbance_level_ind_b )  THEN
       WRITE( message_string, * )  'disturbance_level_ind_t = ', disturbance_level_ind_t,          &
              ' must be >= ', 'disturbance_level_ind_b = ', disturbance_level_ind_b
       CALL message( 'check_parameters', 'PAC0157', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine the horizontal index range for random perturbations.
!-- In case of non-cyclic horizontal boundaries, no perturbations are imposed near the inflow and
!-- the perturbation area is further limited to ...(1) after the initial phase of the flow.

    IF ( bc_lr /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, nx/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > nx )  THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PAC0158', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*nx/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > nx )  THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PAC0159', 1, 2, 0, 6, 0 )
       ENDIF
    ELSEIF ( bc_ns /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, ny/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > ny )  THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PAC0158', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*ny/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > ny )  THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PAC0159', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( random_generator == 'random-parallel' )  THEN
       dist_nxl = nxl;  dist_nxr = nxr
       dist_nys = nys;  dist_nyn = nyn
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
          dist_nxl(1) = MAX( nx - inflow_disturbance_end, nxl )
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr(1) = MIN( inflow_disturbance_end, nxr )
       ELSEIF ( bc_lr == 'nested'  .OR.  bc_lr == 'nesting_offline' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
          dist_nys(1) = MAX( ny - inflow_disturbance_end, nys )
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn(1) = MIN( inflow_disturbance_end, nyn )
       ELSEIF ( bc_ns == 'nested'  .OR.  bc_ns == 'nesting_offline' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
       ENDIF
    ELSE
       dist_nxl = 0;  dist_nxr = nx
       dist_nys = 0;  dist_nyn = ny
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl(1) = nx - inflow_disturbance_end
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = inflow_disturbance_begin
          dist_nxr(1) = inflow_disturbance_end
       ELSEIF ( bc_lr == 'nested'  .OR.  bc_lr == 'nesting_offline' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl    = inflow_disturbance_begin
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys(1) = ny - inflow_disturbance_end
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = inflow_disturbance_begin
          dist_nyn(1) = inflow_disturbance_end
       ELSEIF ( bc_ns == 'nested'  .OR.  bc_ns == 'nesting_offline' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys    = inflow_disturbance_begin
       ENDIF
    ENDIF

    IF ( turbulent_outflow )  THEN
!
!--    Turbulent outflow requires Dirichlet conditions at the respective inflow boundary (so far, a
!--    turbulent outflow is realized at the right side only).
       IF ( bc_lr /= 'dirichlet/radiation' )  THEN
          message_string = 'turbulent_outflow = .T. requires bc_lr = "dirichlet/radiation"'
          CALL message( 'check_parameters', 'PAC0160', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    The ouflow-source plane must lay inside the model domain
       IF ( outflow_source_plane < dx  .OR.  outflow_source_plane > nx * dx )  THEN
          WRITE( message_string, * )  'illegal value for outflow_source_plane: ',                  &
                                      outflow_source_plane
          CALL message( 'check_parameters', 'PAC0161', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Determine damping level index for 1D model
    IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       IF ( damp_level_1d == -1.0_wp )  THEN
          damp_level_1d     = zu(nzt+1)
          damp_level_ind_1d = nzt + 1
       ELSEIF ( damp_level_1d < 0.0_wp  .OR.  damp_level_1d > zu(nzt+1) )  THEN
          WRITE( message_string, * )  'damp_level_1d = ', damp_level_1d, ' must be >= 0.0 and <= ',&
                 zu(nzt+1), ' (zu(nzt+1))'
          CALL message( 'check_parameters', 'PAC0162', 1, 2, 0, 6, 0 )
       ELSE
          DO  k = 1, nzt+1
             IF ( damp_level_1d <= zu(k) )  THEN
                damp_level_ind_1d = k
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDIF

!
!-- Check some other 1d-model parameters
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                                       &
         TRIM( mixing_length_1d ) /= 'blackadar' )  THEN
       message_string = 'mixing_length_1d = "' // TRIM( mixing_length_1d ) // '" is unknown'
       CALL message( 'check_parameters', 'PAC0163', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( dissipation_1d ) /= 'as_in_3d_model'  .AND.                                         &
         TRIM( dissipation_1d ) /= 'detering'        .AND.                                         &
         TRIM( dissipation_1d ) /= 'prognostic' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) // '" is unknown'
       CALL message( 'check_parameters', 'PAC0164', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                                       &
         TRIM( dissipation_1d ) == 'as_in_3d_model' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //                          &
                        '" requires mixing_length_1d = "as_in_3d_model"'
       CALL message( 'check_parameters', 'PAC0165', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check pressure gradient conditions
    IF ( dp_external  .AND.  conserve_volume_flow )  THEN
       WRITE( message_string, * )  'both dp_external and conserve_volume_flo',                     &
              'w are .TRUE. but one of them must be .FALSE.'
       CALL message( 'check_parameters', 'PAC0166', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( dp_external )  THEN
       IF ( dp_level_b < zu(nzb)  .OR.  dp_level_b > zu(nzt) )  THEN
          WRITE( message_string, * )  'dp_level_b = ', dp_level_b, ' is out ',                     &
                 ' of range [zu(nzb), zu(nzt)]'
          CALL message( 'check_parameters', 'PAC0167', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( .NOT. ANY( dpdxy /= 0.0_wp ) )  THEN
          WRITE( message_string, * )  'dp_external is .TRUE. but dpdxy is ze',                     &
                 'ro, i.e. the external pressure gradient will not be applied'
          CALL message( 'check_parameters', 'PAC0168', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ANY( dpdxy /= 0.0_wp )  .AND.  .NOT.  dp_external )  THEN
       WRITE( message_string, * )  'dpdxy is nonzero but dp_external is ',                         &
              '.FALSE., i.e. the external pressure gradient & will not be applied'
       CALL message( 'check_parameters', 'PAC0169', 0, 1, 0, 6, 0 )
    ENDIF
    IF ( conserve_volume_flow )  THEN
       IF ( TRIM( conserve_volume_flow_mode ) /= 'initial_profiles'  .AND.                         &
            TRIM( conserve_volume_flow_mode ) /= 'bulk_velocity' )                                 &
       THEN
          WRITE( message_string, * )  'unknown conserve_volume_flow_mode: "',                      &
                 TRIM( conserve_volume_flow_mode ), '"'
          CALL message( 'check_parameters', 'PAC0170', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic')  .AND.                                   &
            TRIM( conserve_volume_flow_mode ) == 'bulk_velocity' )                                 &
       THEN
          WRITE( message_string, * )  'non-cyclic boundary conditions ',                           &
                 'require  conserve_volume_flow_mode = "initial_profiles"'
          CALL message( 'check_parameters', 'PAC0171', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ( u_bulk /= 0.0_wp .OR. v_bulk /= 0.0_wp )  .AND.                                         &
         ( .NOT. conserve_volume_flow .OR. TRIM( conserve_volume_flow_mode ) /= 'bulk_velocity' ) )&
    THEN
       WRITE( message_string, * )  'nonzero bulk velocity requires ',                              &
              'conserve_volume_flow = .T. and ', 'conserve_volume_flow_mode = "bulk_velocity"'
       CALL message( 'check_parameters', 'PAC0172', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check roughness length, which has to be smaller than dz/2
    IF ( ( constant_flux_layer .OR. INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  &
         .AND.  roughness_length >= 0.5 * dz(1) )  THEN
       message_string = 'roughness_length must be smaller than dz/2'
       CALL message( 'check_parameters', 'PAC0173', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if topography is read from file in case of complex terrain simulations.
    IF ( terrain_following_mapping  .AND.  TRIM( topography ) /= 'read_from_file' )  THEN
       message_string = 'terrain_following_mapping requires topography = "read_from_file"'
       CALL message( 'check_parameters', 'PAC0174', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if vertical grid stretching is switched off in case of complex terrain simulations
    IF ( terrain_following_mapping  .AND.  ANY( dz_stretch_level_start /= -9999999.9_wp ) )  THEN
       message_string = 'vertical grid stretching is not allowed for ' //                          &
                        'terrain_following_mapping = .TRUE.'
       CALL message( 'check_parameters', 'PAC0175', 1, 2, 0, 6, 0 )
    ENDIF

    CALL location_message( 'checking parameters', 'finished' )

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the length of data output intervals. In case of parallel NetCDF output the time levels of
!> the output files need to be fixed. Therefore setting the output interval to 0.0s (usually used to
!> output each timestep) is not possible as long as a non-fixed timestep is used.
!--------------------------------------------------------------------------------------------------!

    SUBROUTINE check_dt_do( dt_do, dt_do_name )

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT (IN) :: dt_do_name !< parin variable name

       REAL(wp), INTENT (INOUT)       :: dt_do      !< data output interval

       IF ( dt_do == 0.0_wp )  THEN
          IF ( dt_fixed )  THEN
             WRITE( message_string, '(A,F9.4,A)' )  'Output at every timestep is wanted (' //      &
                    dt_do_name // ' = 0.0).&'//                                                    &
                    'The output interval is set to the fixed timestep dt = ', dt, ' s.'
             CALL message( 'check_parameters', 'PAC0176', 0, 0, 0, 6, 0 )
             dt_do = dt
          ELSE
             message_string = dt_do_name // ' = 0.0 while using a ' //                             &
                              'variable timestep and parallel netCDF4 is not allowed.'
             CALL message( 'check_parameters', 'PAC0177', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

    END SUBROUTINE check_dt_do



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the bottom and top boundary conditions for humidity and scalars.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE set_bc_scalars( sq, bc_b, bc_t, ibc_b, ibc_t, err_nr_b, err_nr_t )

       IMPLICIT NONE

       CHARACTER (LEN=*)   ::  bc_b       !< bottom boundary condition
       CHARACTER (LEN=*)   ::  bc_t       !< top boundary condition
       CHARACTER (LEN=*)   ::  err_nr_b   !< error number if bottom bc is unknown
       CHARACTER (LEN=*)   ::  err_nr_t   !< error number if top bc is unknown
       CHARACTER (LEN=1)   ::  sq         !< name of scalar quantity


       INTEGER(iwp)        ::  ibc_b      !< index for bottom boundary condition
       INTEGER(iwp)        ::  ibc_t      !< index for top boundary condition

!
!--    Set Integer flags and check for possilbe errorneous settings for bottom boundary condition
       IF ( bc_b == 'dirichlet' )  THEN
          ibc_b = 0
       ELSEIF ( bc_b == 'neumann' )  THEN
          ibc_b = 1
       ELSEIF ( bc_b == 'nested' )  THEN
          ibc_b = 3
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) // '_b ="' //           &
                           TRIM( bc_b ) // '"'
          CALL message( 'check_parameters', err_nr_b, 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Set Integer flags and check for possilbe errorneous settings for top boundary condition
       IF ( bc_t == 'dirichlet' )  THEN
          ibc_t = 0
       ELSEIF ( bc_t == 'neumann' )  THEN
          ibc_t = 1
       ELSEIF ( bc_t == 'initial_gradient' )  THEN
          ibc_t = 2
       ELSEIF ( bc_t == 'nested'  .OR.  bc_t == 'nesting_offline' )  THEN
          ibc_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) // '_t ="' //           &
                           TRIM( bc_t ) // '"'
          CALL message( 'check_parameters', err_nr_t, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE set_bc_scalars


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check for consistent settings of bottom boundary conditions for humidity and scalars.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE check_bc_scalars( sq, bc_b, ibc_b, err_nr_1, err_nr_2, constant_flux,               &
                                 surface_initial_change )

       IMPLICIT NONE

       CHARACTER (LEN=*)   ::  bc_b                     !< bottom boundary condition
       CHARACTER (LEN=*)   ::  err_nr_1                 !< error number of first error
       CHARACTER (LEN=*)   ::  err_nr_2                 !< error number of second error
       CHARACTER (LEN=1)   ::  sq                       !< name of scalar quantity


       INTEGER(iwp)        ::  ibc_b                    !< index of bottom boundary condition

       LOGICAL             ::  constant_flux            !< flag for constant-flux layer

       REAL(wp)            ::  surface_initial_change   !< value of initial change at the surface

!
!--    A given surface value implies Dirichlet boundary condition for the respective quantity. In
!--    this case specification of a constant flux is forbidden. However, an exception is made for
!--    large-scale forcing as well as land-surface model.
       IF ( .NOT. land_surface  .AND.  .NOT. large_scale_forcing )  THEN
          IF ( ibc_b == 0  .AND.  constant_flux )  THEN
             message_string = 'boundary condition: bc_' // TRIM( sq ) // '_b = "' //               &
                              TRIM( bc_b ) // '" is not allowed with prescribed surface flux'
             CALL message( 'check_parameters', err_nr_1, 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( constant_flux  .AND.  surface_initial_change /= 0.0_wp )  THEN
          WRITE( message_string, * )  'a prescribed surface flux is not allowed with ', sq,        &
                '_surface_initial_change (/=0) = ', surface_initial_change
          CALL message( 'check_parameters', err_nr_2, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE check_bc_scalars



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Filter a list of strings for duplicates. The output list is of same length but does not contain
!> any duplicates.
!--------------------------------------------------------------------------------------------------!

    SUBROUTINE filter_duplicate_strings( string_length, string_list, string_list_filtered )

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: string_length  !< length of a single string in list

       CHARACTER (LEN=*),             DIMENSION(:), INTENT(IN)               ::  string_list           !< non-filtered list
       CHARACTER (LEN=string_length), DIMENSION(:), INTENT(OUT), ALLOCATABLE ::  string_list_filtered  !< filtered list

       CHARACTER (LEN=string_length) ::  string_to_check  !< string to be check for duplicates

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index
       INTEGER(iwp) ::  nstrings  !< number of strings in list

       LOGICAL ::  is_duplicate  !< true if string has duplicates in list


       nstrings = UBOUND( string_list, DIM=1 )
       ALLOCATE( string_list_filtered(1:nstrings) )
       string_list_filtered(:) = ' '

       k = 0
       DO  i = 1, nstrings
          string_to_check = string_list(i)

          IF ( string_to_check == ' ' )  EXIT

          is_duplicate = .FALSE.
          DO  j = 1, nstrings
             IF ( string_list_filtered(j) == ' ' )  THEN
                EXIT
             ELSEIF ( string_to_check == string_list_filtered(j) )  THEN
                is_duplicate = .TRUE.
                EXIT
             ENDIF
          ENDDO

          IF ( .NOT. is_duplicate )  THEN
             k = k + 1
             string_list_filtered(k) = string_to_check
          ENDIF

       ENDDO

    END SUBROUTINE filter_duplicate_strings


 END SUBROUTINE check_parameters
