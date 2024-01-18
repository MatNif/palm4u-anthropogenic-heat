! !> @file palm.f90
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
!> Large-Eddy Simulation (LES) model for atmospheric and oceanic boundary-layer flows,
!> see the PALM homepage https://palm-model.org for further information
!--------------------------------------------------------------------------------------------------!
 PROGRAM palm

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d

    USE control_parameters,                                                                        &
        ONLY:  coupling_char,                                                                      &
               debug_output,                                                                       &
               do2d_at_begin,                                                                      &
               do3d_at_begin,                                                                      &
               io_blocks,                                                                          &
               io_group,                                                                           &
               open_debug_files,                                                                   &
               restart_data_format_output,                                                         &
               runnr,                                                                              &
               simulated_time_chr,                                                                 &
               spinup,                                                                             &
               time_since_reference_point,                                                         &
               write_binary,                                                                       &
               write_spinup_data

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  initializing_actions,                                                               &
               nested_run
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_statistics,                                                                     &
               log_point

#if defined( __parallel )
    USE cpulog,                                                                                    &
        ONLY:  log_point_s
#endif

    USE diagnostic_output_quantities_mod,                                                          &
        ONLY:  doq_calculate

    USE kinds

    USE module_interface,                                                                          &
        ONLY:  module_interface_init_numerics,                                                     &
               module_interface_init_output,                                                       &
               module_interface_last_actions


    USE multi_agent_system_mod,                                                                    &
        ONLY:  agents_active,                                                                      &
               mas_last_actions

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  netcdf_data_input_inquire_file,                                                     &
               netcdf_data_input_init,                                                             &
               netcdf_data_input_surface_data

    USE pegrid

#if defined( __parallel )
    USE particle_attributes,                                                                       &
        ONLY:  particle_advection

    USE pmc_particle_interface,                                                                    &
        ONLY:  pmcp_g_alloc_win

    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run,                                                       &
               pmci_child_initialize,                                                              &
               pmci_finalize,                                                                      &
               pmci_init,                                                                          &
               pmci_modelconfiguration,                                                            &
               pmci_parent_initialize
#endif

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_close,                                                                    &
               rd_mpi_io_open

    USE surface_data_output_mod,                                                                   &
        ONLY:  surface_data_output_last_action

    USE topography_mod,                                                                           &
        ONLY:  init_topography

    USE write_restart_data_mod,                                                                    &
        ONLY:  wrd_global,                                                                         &
               wrd_global_spinup,                                                                  &
               wrd_local,                                                                          &
               wrd_local_spinup

#if defined( __parallel )  &&  defined( _OPENACC )
    USE openacc
#endif


    IMPLICIT NONE

!
!-- Local variables
    CHARACTER(LEN=9) ::  time_to_string  !<

    INTEGER(iwp)             ::  i               !< loop counter for blocked I/O
#if defined( __parallel) && defined( _OPENACC )
    INTEGER(acc_device_kind) ::  device_type      !< device type for OpenACC
    INTEGER(iwp)             ::  local_comm       !< local communicator (shared memory)
    INTEGER(iwp)             ::  local_num_procs  !< local number of processes
    INTEGER(iwp)             ::  local_id         !< local id
    INTEGER(iwp)             ::  num_devices      !< number of devices visible to OpenACC
    INTEGER(iwp)             ::  my_device        !< device used by this process
#endif

#if defined( __parallel )
!
!-- MPI initialisation. comm2d is preliminary set, because it will be defined in init_pegrid but is
!-- used before in cpu_log.
    CALL MPI_INIT( ierr )

!
!-- Initialize the coupling for nested-domain runs comm_palm is the communicator which includes all
!-- PEs (MPI processes) available for this (nested) model. If it is not a nested run, comm_palm is
!-- returned as MPI_COMM_WORLD.
    CALL cpu_log( log_point_s(70), 'pmci_init', 'start' )
    CALL pmci_init( comm_palm )
    CALL cpu_log( log_point_s(70), 'pmci_init', 'stop' )
    comm2d = comm_palm
!
!-- Get the (preliminary) number of MPI processes and the local PE-id (in case of a further
!-- communicator splitting in init_coupling, these numbers will be changed in init_pegrid).
    IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN

       CALL MPI_COMM_SIZE( comm_palm, numprocs, ierr )
       CALL MPI_COMM_RANK( comm_palm, myid, ierr )

    ELSE

       CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    ENDIF

!
!-- Set the suffixes for files to be read or written by cores separately, based on the above
!-- preliminary communicators. myid_char will be set again in routine init_pegrid, based on the
!-- final communicator comm2d, however it is required already when the global restart data are read
!-- in Fortran binary format from routine parin via rrd_global.
!-- ATTENTION: It is still unclear, if the ID numbering for comm2d is the same as for the
!-- preliminary communicators. If the order is the same (and it should be, otherwise restart files
!-- may be out of order) the setting of myid_char within init_pegrid could be removed.
    WRITE( myid_char,'(''_'',I6.6)' )  myid

#ifdef _OPENACC
!
!-- Select OpenACC device to use in this process. For this find out how many neighbors there are
!-- running on the same node and which id this process is.
    IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN
       CALL MPI_COMM_SPLIT_TYPE( comm_palm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, local_comm,    &
                                 ierr )
    ELSE
       CALL MPI_COMM_SPLIT_TYPE( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,           &
                                 local_comm, ierr )
    ENDIF
    CALL MPI_COMM_SIZE( local_comm, local_num_procs, ierr )
    CALL MPI_COMM_RANK( local_comm, local_id, ierr )

!
!-- This loop including the barrier is a workaround for PGI compiler versions up to and including
!-- 18.4. Later releases are able to select their GPUs in parallel, without running into spurious
!-- errors.
    DO i = 0, local_num_procs-1
       CALL MPI_BARRIER( local_comm, ierr )

       IF ( i == local_id )  THEN
          device_type = acc_get_device_type()
          num_devices = acc_get_num_devices( device_type )
          my_device = MOD( local_id, num_devices )
          CALL acc_set_device_num( my_device, device_type )
       ENDIF
    ENDDO

    CALL MPI_COMM_FREE( local_comm, ierr )
#endif
#endif

!
!-- Initialize measuring of the CPU-time remaining to the run.
    CALL local_tremain_ini

!
!-- Start of total CPU time measuring.
    CALL cpu_log( log_point(1), 'total', 'start' )
    CALL cpu_log( log_point(2), 'initialisation', 'start' )

!
!-- Read control parameters from NAMELIST files and read environment-variables.
    CALL parin

!
!-- Determine processor topology and local array indices.
    CALL init_pegrid
!
!-- Open a file for debug output.
    IF ( open_debug_files .OR. debug_output )  THEN
       OPEN( 9, FILE='DEBUG'//TRIM( coupling_char )//myid_char, FORM='FORMATTED' )
    ENDIF
!
!-- Check if input file according to input-data standard exists.
    CALL netcdf_data_input_inquire_file
!
!-- Generate grid parameters.
    CALL init_grid
!
!-- Initialize topography.
    CALL init_topography
!
!-- Initialize boundary conditions and numerics such as the multigrid solver or the advection
!-- routine.
    CALL module_interface_init_numerics
!
!-- Read global attributes if available.
    CALL netcdf_data_input_init
!
!-- Read surface classification data, e.g. vegetation and soil types, water surfaces, etc., if
!-- available. Some of these data is required before check parameters is invoked.
    CALL netcdf_data_input_surface_data
!
!-- Check control parameters and deduce further quantities.
    CALL check_parameters

    CALL init_3d_model

    CALL module_interface_init_output

#if defined( __parallel )
!
!-- Coupling protocol setup for nested-domain runs
    IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN
       CALL pmci_modelconfiguration
!
!--    Receive and interpolate initial data on children.
!--    Child initialization must be made first if the model is both child and parent if necessary.
       IF ( TRIM( initializing_actions )  /= 'read_restart_data' )  THEN
          CALL pmci_child_initialize
!
!--       Send initial condition data from parent to children
          CALL pmci_parent_initialize
       ENDIF

!
!--    Allocation of MPI windows for particle transfer
!>     TODO: This needs to be moved to the place where the MPI windows for the Eulerian quantities
!>           are allocated. It must be located after pmci_child/parent_initialize
       IF ( particle_advection  .AND.  nested_run )  THEN
          CALL pmcp_g_alloc_win
       ENDIF
    ENDIF
#endif

!
!-- Output of program header
    IF ( myid == 0 )  CALL header

    CALL cpu_log( log_point(2), 'initialisation', 'stop' )
!
!-- Integration of the non-atmospheric equations (land surface model, urban surface model)
    IF ( spinup )  THEN
       CALL cpu_log( log_point(41), 'wall/soil spinup', 'start' )
       CALL time_integration_spinup
       CALL cpu_log( log_point(41), 'wall/soil spinup', 'stop' )
!
!--    Write spinup data to file
       IF ( write_spinup_data )  THEN
          CALL cpu_log( log_point(27), 'write-spinup-data', 'start' )
          CALL location_message( 'writing spinup data', 'start' )
!
!--       Open MPI-IO spinup file
          CALL rd_mpi_io_open( 'write', 'SPINUPOUT' // TRIM( coupling_char ) )
!
!--       Write control parameters and other global variables for spinup file.
          CALL wrd_global_spinup
!
!--       Write processor specific surface data for spinup runs
          CALL wrd_local_spinup
!
!--       Close spinup File
          CALL rd_mpi_io_close

          CALL location_message( 'writing spinup data', 'finished' )
          CALL cpu_log( log_point(27), 'write-spinup-data', 'stop' )
       ENDIF
    ENDIF
!
!-- Set start time in format hh:mm:ss
    simulated_time_chr = time_to_string( time_since_reference_point )
!
!-- If required, output of initial arrays
    IF ( do2d_at_begin )  THEN
       CALL doq_calculate    !TODO, will be called twice

       CALL data_output_2d( 'xy', 0 )
       CALL data_output_2d( 'xz', 0 )
       CALL data_output_2d( 'yz', 0 )
    ENDIF

    IF ( do3d_at_begin )  THEN
       CALL doq_calculate    !TODO, will be called twice

       CALL data_output_3d( 0 )
    ENDIF
!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

!
!-- If required, write binary data for restart runs
    IF ( write_binary )  THEN

       CALL cpu_log( log_point(22), 'write-restart-data', 'start' )

       CALL location_message( 'writing restart data', 'start' )

       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          DO  i = 0, io_blocks-1
             IF ( i == io_group )  THEN

!
!--             Open binary file
                CALL check_open( 14 )
!
!--             Write control parameters and other global variables for restart.
                IF ( myid == 0 )  CALL wrd_global
!
!--             Write processor specific flow field data for restart runs
                CALL wrd_local
!
!--             Close binary file
                CALL close_file( 14 )

             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO

       ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--       Open MPI-IO restart file
          CALL rd_mpi_io_open( 'write', 'BINOUT' // TRIM( coupling_char ) )
!
!--       Write control parameters and other global variables for restart.
          CALL wrd_global
!
!--       Write processor specific flow field data for restart runs
          CALL wrd_local
!
!--       Close restart File
          CALL rd_mpi_io_close

       ENDIF

       CALL location_message( 'writing restart data', 'finished' )

       CALL cpu_log( log_point(22), 'write-restart-data', 'stop' )

    ENDIF
!
!-- Last actions for surface output, for instantaneous and time-averaged data
    CALL surface_data_output_last_action( 0 )
    CALL surface_data_output_last_action( 1 )

!
!-- If required, repeat output of header including the required CPU-time
    IF ( myid == 0 )  CALL header
!
!-- Perform module specific last actions
    CALL cpu_log( log_point(4), 'last actions', 'start' )

    IF ( myid == 0 .AND. agents_active )  CALL mas_last_actions ! ToDo: move to module_interface

    CALL module_interface_last_actions

    CALL cpu_log( log_point(4), 'last actions', 'stop' )

!
!-- Close files
    CALL close_file( 0 )

!
!-- Write run number to file (used by palmrun to create unified cycle numbers for output files).
    IF ( myid == 0  .AND.  runnr > 0 )  THEN
       OPEN( 90, FILE='RUN_NUMBER', FORM='FORMATTED' )
       WRITE( 90, '(I4)' )  runnr
       CLOSE( 90 )
    ENDIF

#if defined( __parallel )
!
!-- Free all MPI windows used by the pmc. This is done to avoid aborts or hanging situations
!-- in MPI_FINALIZE which sometimes appear when Intel MPI is used.
    IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN
       CALL pmci_finalize
    ENDIF
#endif

!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
    CALL cpu_statistics

#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

 END PROGRAM palm
