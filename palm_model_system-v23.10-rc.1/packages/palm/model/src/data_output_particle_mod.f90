!> @file data_output_particle_mod.f90
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
!> Output of particle time series
!>
!> @bug Calling netcdf4_init_module will break the setup of the module if used for output by other
!>      modules. This call must be removed.
!> @todo Remove any calls to data_output_netcdf4_module and replace by calls to data_output_module.
!--------------------------------------------------------------------------------------------------!
 MODULE data_output_particle_mod

#if defined( __parallel )
    USE MPI
#endif

#if defined( __netcdf4 )
    USE NETCDF
#endif

    USE control_parameters,                                                                        &
        ONLY:  coupling_char,                                                                      &
               end_time,                                                                           &
               message_string,                                                                     &
               nested_run,                                                                         &
               run_description_header,                                                             &
               simulated_time,                                                                     &
               simulated_time_at_begin

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE data_output_netcdf4_module,                                                                &
        ONLY:  netcdf4_finalize,                                                                   &
               netcdf4_get_error_message,                                                          &
               netcdf4_init_dimension,                                                             &
               netcdf4_init_module,                                                                &
               netcdf4_init_variable,                                                              &
               netcdf4_inquire_dimension,                                                          &
               netcdf4_inq_varid,                                                                  &
               netcdf4_open_file,                                                                  &
               netcdf4_stop_file_header_definition,                                                &
               netcdf4_write_attribute,                                                            &
               netcdf4_write_variable

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nnx,                                                                                &
               nny,                                                                                &
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

    USE particle_attributes,                                                                       &
        ONLY:  data_output_pts,                                                                    &
               dt_dopts,                                                                           &
               extend_prts_filesize,                                                               &
               grid_particles,                                                                     &
               grid_particle_def,                                                                  &
               maximum_number_of_output_particles,                                                 &
               number_of_particles,                                                                &
               number_of_particle_groups,                                                          &
               particles,                                                                          &
               particle_advection_start,                                                           &
               particle_id_file_found,                                                             &
               particle_type,                                                                      &
               prt_count,                                                                          &
               pts_increment,                                                                      &
               pts_percentage,                                                                     &
               unlimited_dimension,                                                                &
               zero_particle

    USE pegrid,                                                                                    &
        ONLY:  collective_wait,                                                                    &
               comm1dx,                                                                            &
               comm1dy,                                                                            &
               comm2d,                                                                             &
               myid,                                                                               &
               myidx,                                                                              &
               myidy,                                                                              &
               npex,                                                                               &
               npey,                                                                               &
               numprocs

#if defined( __parallel )
    USE pmc_interface,                                                                             &
        ONLY:  comm_world_nesting,                                                                 &
               particle_coupling,                                                                  &
               pmc_get_model_info,                                                                 &
               root_model

    USE pmc_particle_interface,                                                                    &
        ONLY:  cnpo,                                                                               &
               pmcp_g_nested_init
#else
    USE pmc_interface,                                                                             &
        ONLY:  particle_coupling,                                                                  &
               root_model
#endif
    USE shared_memory_io_mod,                                                                      &
        ONLY:  sm_class

    IMPLICIT NONE

    PRIVATE

    SAVE

    CHARACTER(LEN=32)   ::  file_name        !< name of particle data NetCDF file
    CHARACTER(LEN=4000) ::  var_list  = ' '  !< variable list for global attribute

    INTEGER, PARAMETER ::  max_nr_prt_quantities         = 128  !< maximum number of allowed particle quantities for output
    INTEGER, PARAMETER ::  max_nr_prt_quantities_const   =  16  !< maximum number of allowed constant particle quantitities for output
    INTEGER, PARAMETER, PUBLIC ::  nr_prt_quantities_statistical =  26  !< number of statistical particle quantities for output

    INTEGER(iwp) ::  io_end_index            !<
    INTEGER(iwp) ::  io_start_index          !< start index of the output area on IO thread
    INTEGER(iwp) ::  nr_out_prts_on_this_pe  !< number of particles assigned for output on this PE
    INTEGER(iwp) ::  nr_out_prts_rest        !< number of rest particles in case of irregular distribution of particles on PEs
    INTEGER(iwp) ::  nr_time_values          !< number of values on time axis
    INTEGER(iwp) ::  pe_end_index            !< start index of the output area on this PE
    INTEGER(iwp) ::  pe_start_index          !< start index of the output area on this PE
    INTEGER(iwp) ::  pts_time_index = 0      !< index of time axis

    INTEGER(iwp), PUBLIC ::  max_nr_out_prts  !< maximum number of particles allowed for output; public because contained in restart file
    INTEGER(iwp), PUBLIC ::  nr_out_prts      !< number of particles currently scheduled for output; may increase during run in case of additionally released particles; public because contained in restart file

    INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  nr_prts_prev_pes_ycolumn  !< total number of particles in previous PE columns along y
    INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  start_prt_count_yz        !< value of particle counter at (nys,nzb) of a yz-cross-section at i

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  io_indices              !< indices on IO processes
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  mo_indices              !< indices for model communicator
    INTEGER(idp), ALLOCATABLE, DIMENSION(:,:) ::  nr_prts_along_yz_on_pe  !< number of particles for each yz-plane on PE
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  nr_out_prts_remote      !< number of particles scheduled for output on all remote PEs
                                                                          !< first index 1: accumulated value up to value of second index
                                                                          !< first index 2: value on second index PE
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  rma_particles           !< itart address and number of remote particles
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  sh_indices              !< indices in shared memory group
    INTEGER(idp), ALLOCATABLE, DIMENSION(:,:) ::  nr_prts_on_pe           !< total number of particles on PE

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:,:) ::  new_prt_count  !< prt_count for newly released particles

    LOGICAL ::  irregular_distribution  !< irregular distribution of output particles
    LOGICAL, PUBLIC ::  dop_individual_ts = .FALSE.  !< switch to activate output of individual particle time series data

    REAl(wp), PARAMETER ::  eps = 0.00001  !<

    REAL(sp), POINTER, DIMENSION(:), CONTIGUOUS ::  time_axis_values  !< time axis Values

    TYPE particle_feature
       INTEGER(iwp)      ::  var_id                !<
       CHARACTER(LEN=32) ::  name                  !<
       CHARACTER(LEN=32) ::  units                 !<
       LOGICAL           ::  is_integer = .FALSE.  !<
    END TYPE particle_feature

    TYPE, EXTENDS(particle_feature) :: statistical_particle_feature
       REAL(sp), DIMENSION(:), POINTER, CONTIGUOUS ::  value  !<
    END TYPE statistical_particle_feature

    TYPE(sm_class) ::  prt_mio  !< manage communicator for particle IO

    TYPE(particle_feature), DIMENSION(max_nr_prt_quantities_const) ::  prt_quantity_const
    TYPE(particle_feature), DIMENSION(max_nr_prt_quantities) ::  prt_quantity
    TYPE(statistical_particle_feature), DIMENSION(:,:), ALLOCATABLE ::  prt_quantity_statistical

    INTEGER, DIMENSION(2) ::  dimension_ids  !<


!
!-- NetCDF.
    INTEGER(iwp) ::  file_id = -1                       !< id of Netcdf file
    INTEGER(iwp) ::  nr_prt_quantities_const_scheduled  !< number of constant particle quantities scheduled for output
    INTEGER(iwp) ::  nr_prt_quantities_scheduled        !< number of particle quantities scheduled for output

    TYPE dimension_id
       INTEGER(iwp) ::  prt
       INTEGER(iwp) ::  time
    END TYPE dimension_id

    TYPE variable_id
       INTEGER(iwp) ::  prt
       INTEGER(iwp) ::  time
    END TYPE variable_id

    TYPE(dimension_id) ::  did
    TYPE(variable_id)  ::  var_id

!
!-- Shared memory buffer.
    INTEGER(iwp), POINTER, CONTIGUOUS, DIMENSION(:) ::  out_buf_i  !< integer output buffer
    REAL(sp), POINTER, CONTIGUOUS, DIMENSION(:)     ::  out_buf_r  !< real output buffer

!
!-- Particle list in file.
    INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  particle_id_scheduled_for_output

!
!-- RMA window.
#if defined( __parallel )
    INTEGER(iwp) ::  win_rma_buf_i   !<
    INTEGER(iwp) ::  win_rma_buf_r   !<
    INTEGER(iwp) ::  win_prt_i = -1  !< integer MPI shared memory window
    INTEGER(iwp) ::  win_prt_r = -1  !< real MPI shared memory window
#endif

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  remote_indices  !< particle number of the remote particles,
                                                                !< used as indices in the output array
    INTEGER(iwp), POINTER, DIMENSION(:) ::  transfer_buffer_i  !< rma window to provide data, which can be
                                                               !< fetch via MPI_Get
    REAL(sp), POINTER, DIMENSION(:) ::  transfer_buffer_r  !< same for REAL

!
!-- Public subroutine interface.
#if defined( __parallel )
    INTERFACE dop_alloc_rma_mem
       MODULE PROCEDURE dop_alloc_rma_mem_i1
       MODULE PROCEDURE dop_alloc_rma_mem_r1
    END INTERFACE dop_alloc_rma_mem
#endif

    INTERFACE dop_collect_statistics
       MODULE PROCEDURE dop_collect_statistics
    END INTERFACE dop_collect_statistics

    INTERFACE dop_data_output_ptseries
       MODULE PROCEDURE dop_data_output_ptseries
    END INTERFACE dop_data_output_ptseries

    INTERFACE dop_finalize
       MODULE PROCEDURE dop_finalize
    END INTERFACE dop_finalize

    INTERFACE dop_init
       MODULE PROCEDURE dop_init
    END INTERFACE dop_init

    INTERFACE dop_write_tseries_data
       MODULE PROCEDURE dop_write_tseries_data
    END INTERFACE dop_write_tseries_data

#if defined( __parallel )
    PUBLIC dop_alloc_rma_mem  ! must be PUBLIC on NEC, although if it is only used in submodule
#endif
    PUBLIC dop_collect_statistics,                                                                 &
           dop_data_output_ptseries,                                                               &
           dop_finalize,                                                                           &
           dop_init,                                                                               &
           dop_write_tseries_data


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_init( read_restart )

    IMPLICIT NONE

#if defined( __parallel )
    INTEGER(iwp) ::  ierr  !< MPI error code
#endif
    INTEGER(iwp) ::  nr_local_last_pe  !< number of output particles on myid == numprocs-2
    INTEGER(idp) ::  nr_out_prts_8     !< total number of output particles in 64 bit

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  nr_out_prts_on_pe    !< number of output particles on the PEs
    INTEGER(idp), DIMENSION(0:numprocs-1) ::  nr_out_prts_on_pe_8  !< number of output particles on the PEs in 64bit

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  io_indices_s  !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  mo_indices_s  !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  sh_indices_s  !<

    LOGICAL, INTENT(IN) ::  read_restart  !<

    REAL(wp) ::  rest    !< temporary used for calculating the number of output time levels
    REAL(wp) ::  rvalue  !< temporary used for calculating the number of output time levels


    IF ( dop_individual_ts )  THEN

       IF ( root_model )  THEN

          ALLOCATE( nr_prts_on_pe(0:npey-1,0:npex-1) )
          ALLOCATE( nr_prts_along_yz_on_pe(nxl:nxr,0:npey-1) )
          ALLOCATE( start_prt_count_yz(nxl:nxr) )
          ALLOCATE( nr_prts_prev_pes_ycolumn(0:npex-1) )
          ALLOCATE( new_prt_count(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

          IF ( .NOT. particle_id_file_found )  THEN

             IF ( .NOT. read_restart)  THEN

                CALL set_output_particle_ids_undefined

                nr_out_prts_on_pe(:) = 0
                CALL count_output_particles( nr_out_prts_on_pe(myid), prt_count )

#if defined( __parallel )
                CALL MPI_ALLREDUCE( MPI_IN_PLACE, nr_out_prts_on_pe, SIZE( nr_out_prts_on_pe ),    &
                                    MPI_INTEGER, MPI_SUM, comm2d, ierr )
#endif

                nr_out_prts = SUM( nr_out_prts_on_pe )
!
!--             Use 64-bit INTEGER, because number of particles may exceed 2 GByte.
!--             Anyhow, output of particle numbers > 2**32 is not allowed due to MPI array size
!--             restrictions.
                nr_out_prts_on_pe_8 = nr_out_prts_on_pe
                nr_out_prts_8 = SUM( nr_out_prts_on_pe_8 )
                IF ( nr_out_prts_8 > HUGE( nr_out_prts_on_pe ) )  THEN
                   WRITE( message_string, '(A,I20)' ) 'number of particles for output too large:', &
                                                      nr_out_prts_8
                   CALL message( 'dop_init', 'LPM0029', 3, 2, 0, 6, 0 )
                ENDIF
!
!--             Calculate/set the maximum number of particles allowed for output. Maybe reserve
!--             space for particles that are released at later times. Reservation is steered by
!--             user via namelist parameters maximum_number_of_output_particles or
!--             extend_prts_filesize.
!--             The variable max_nr_out_prts is used to defined the particle axis (dimension)
!--             size of the NetCDF file. maximum_number_of_output_particles is given by user
!--             via namelist.
                max_nr_out_prts = nr_out_prts_8
                IF ( maximum_number_of_output_particles > 0 )  THEN
                   max_nr_out_prts = MAX( maximum_number_of_output_particles, max_nr_out_prts )
                ELSEIF ( extend_prts_filesize > 0.0_wp )  THEN
!
!--                Double precision is required to avoid round-off errors.
                   max_nr_out_prts = max_nr_out_prts *                                             &
                                     ( 100.0_dp + extend_prts_filesize ) / 100.0_dp
                ENDIF

             ENDIF

          ELSE

             CALL set_output_particle_ids_undefined
             CALL dop_read_ids_scheduled_for_ts_output( nr_out_prts )
             max_nr_out_prts = nr_out_prts

          ENDIF

!
!--       The maximum number of particles must be at least the number of MPI processes.
          max_nr_out_prts = MAX( max_nr_out_prts, numprocs )
!
!--       Number of particles scheduled for output on this PE.
          nr_out_prts_on_this_pe = ( max_nr_out_prts + numprocs - 1 ) / numprocs
!
!--       Output numbering on this PE.
          pe_start_index = myid * nr_out_prts_on_this_pe + 1
          pe_end_index   = MIN( ( myid + 1 ) * nr_out_prts_on_this_pe, max_nr_out_prts )

          irregular_distribution = .FALSE.

#if defined( __parallel )
!
!--       In case of few particles, it may happen that not only the last PE gets fewer output
!--       particles. In this case, the local number of particles on PE numprocs-1 will be < 1.
!--       If this happens, an irregular distribution of output particles will be used, where also
!--       other PEs than the last one will be assigned fewer particles.
          IF ( myid == numprocs-1 )  THEN
             nr_local_last_pe = pe_end_index - pe_start_index + 1
          ELSE
             nr_local_last_pe = 0
          ENDIF
          CALL MPI_BCAST( nr_local_last_pe, 1, MPI_INTEGER, numprocs-1, comm2d, ierr )
#else
          nr_local_last_pe = nr_out_prts_on_this_pe
#endif
          IF ( nr_local_last_pe < 1 )  THEN
             irregular_distribution = .TRUE.
             CALL dop_setup_ireg_distribution
          ENDIF
!
!--       Set contiguous particle ID for the output particles.
          IF ( .NOT. read_restart  .AND.  .NOT. particle_id_file_found )  THEN
             CALL dop_set_io_ids
          ENDIF
!
!--       Set particle IDs for particles that have been selected by the user.
          IF ( particle_id_file_found )  THEN
             CALL dop_mark_ids_scheduled_for_output_via_file
          ENDIF
!
!--       Prepare/provide the shared memory arrays for particle output.
          CALL prt_mio%sm_init_data_output_particles( )

          ALLOCATE( sh_indices_s(2,0:prt_mio%sh_npes-1) )
          ALLOCATE( sh_indices(2,0:prt_mio%sh_npes-1) )

#if defined( __parallel )
          sh_indices_s = 0
          sh_indices_s(1,prt_mio%sh_rank) = pe_start_index
          sh_indices_s(2,prt_mio%sh_rank) = pe_end_index

          CALL MPI_ALLREDUCE( sh_indices_s, sh_indices, 2*prt_mio%sh_npes, MPI_INTEGER, MPI_SUM,   &
                              prt_mio%comm_shared, ierr )
!
!--       Output numbering on actual IO PE.
          io_start_index = sh_indices(1,0)
          io_end_index   = sh_indices(2,prt_mio%sh_npes-1)
#else
!
!--       Output numbering.
          io_start_index = pe_start_index
          io_end_index   = pe_end_index
#endif


#if defined( __parallel )
          CALL MPI_BCAST( prt_mio%io_npes, 1, MPI_INTEGER, 0,  prt_mio%comm_shared, ierr )
#endif
          ALLOCATE( io_indices(2,0:prt_mio%io_npes-1) )
          IF ( prt_mio%iam_io_pe )  THEN
             ALLOCATE( io_indices_s(2,0:prt_mio%io_npes-1) )

             io_indices_s = 0
             io_indices_s(1,prt_mio%io_rank) = io_start_index
             io_indices_s(2,prt_mio%io_rank) = io_end_index

#if defined( __parallel )
             CALL MPI_ALLREDUCE( io_indices_s, io_indices, 2 * prt_mio%io_npes, MPI_INTEGER,       &
                                 MPI_SUM, prt_mio%comm_io, ierr )
#else
             io_indices = io_indices_s
#endif
          ENDIF

#if defined( __parallel )
          CALL MPI_BCAST( io_indices, SIZE( io_indices ), MPI_INTEGER, 0,  prt_mio%comm_shared,    &
                          ierr )
#endif

          ALLOCATE( nr_out_prts_remote(2,0:numprocs-1) )
          ALLOCATE( rma_particles(2,0:numprocs-1) )

          ALLOCATE( mo_indices(2,0:numprocs-1) )
          ALLOCATE( mo_indices_s(2,0:numprocs-1) )

          mo_indices_s = 0
          mo_indices_s(1,myid) = pe_start_index
          mo_indices_s(2,myid) = pe_end_index

#if defined( __parallel )
          CALL MPI_ALLREDUCE( mo_indices_s, mo_indices, 2 * numprocs, MPI_INTEGER, MPI_SUM,        &
                              comm2d, ierr )
#else
          mo_indices = mo_indices_s
#endif
!
!--       Allocate output buffer
#if defined( __parallel )
          CALL prt_mio%sm_allocate_shared( out_buf_r, io_start_index, io_end_index, win_prt_r )
          CALL prt_mio%sm_allocate_shared( out_buf_i, io_start_index, io_end_index, win_prt_i )
#else
          ALLOCATE( out_buf_r(io_start_index:io_end_index) )
          ALLOCATE( out_buf_i(io_start_index:io_end_index) )
#endif
       ENDIF

!
!--    Prepare nested particle output.
       IF ( nested_run  .AND.  particle_coupling )  THEN

#if defined( __parallel )
          CALL pmcp_g_nested_init

          IF ( root_model )  THEN

             CALL MPI_BCAST( mo_indices, SIZE(mo_indices), MPI_INTEGER, 0, comm_world_nesting,     &
                             ierr )
!
!--          Prepare the NetCDF output.
             CALL dop_netcdf_setup( )

             CALL MPI_BCAST( nr_prt_quantities_const_scheduled, 1, MPI_INTEGER, 0,                 &
                             prt_mio%comm_shared, ierr )
             CALL MPI_BCAST( nr_prt_quantities_scheduled, 1, MPI_INTEGER, 0, prt_mio%comm_shared,  &
                             ierr )

             CALL dop_count_particles_to_be_transfered_to_root
             CALL dop_write_constant_quantities

             CALL deallocate_and_free

          ELSE

             ALLOCATE( mo_indices(2,0:cnpo%rem_size-1) )
             ALLOCATE( nr_out_prts_remote(2,0:cnpo%rem_size-1) )
             nr_out_prts_remote = 0

             CALL MPI_BCAST( mo_indices, SIZE(mo_indices), MPI_INTEGER, 0, comm_world_nesting,     &
                             ierr )
          ENDIF
#endif
       ELSE
!
!--       Prepare the NetCDF output.
          CALL dop_netcdf_setup( )

#if defined( __parallel )
          CALL MPI_BCAST( nr_prt_quantities_const_scheduled, 1, MPI_INTEGER, 0,                    &
                          prt_mio%comm_shared, ierr )
          CALL MPI_BCAST( nr_prt_quantities_scheduled, 1, MPI_INTEGER, 0, prt_mio%comm_shared,     &
                          ierr )
#endif

          CALL dop_count_particles_to_be_transfered_to_root
          CALL dop_write_constant_quantities

          CALL deallocate_and_free

       ENDIF

    ELSE

       IF ( root_model )  THEN
          CALL dop_netcdf_setup_statistic_only
       ENDIF

    ENDIF

!
!-- Calculate the number of output intervals for the case of a limited time dimension.
!-- Be aware that the CEILING function returns a value of 1, even if dt_dopts is set larger than
!-- end_time.
    rvalue = ( end_time - MAX( particle_advection_start, simulated_time_at_begin ) ) / dt_dopts
    rest = rvalue - AINT( rvalue )
    nr_time_values = CEILING( rvalue )
!
!-- Take into account, that data at the particle start time are always output.
    IF ( rest == 0.0_wp )  nr_time_values = nr_time_values + 1

!
!-- This array is used to store the time axis values, used for both the individual time series and
!-- the statistics.
    ALLOCATE( time_axis_values(nr_time_values) )

    IF ( root_model )  THEN
       CALL dop_statistics_init
    ENDIF


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> TODO: add explanation about how the particles are distributed in this iregular case,
!>       maybe by giving an example.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_setup_ireg_distribution

    IMPLICIT NONE

!
!-- Number of particles scheduled for output on this PE.
    nr_out_prts_on_this_pe = ( max_nr_out_prts ) / numprocs
    nr_out_prts_rest       = max_nr_out_prts - numprocs * nr_out_prts_on_this_pe
    nr_out_prts_on_this_pe = nr_out_prts_on_this_pe + 1
!
!-- Output numbering on this PE.
    IF ( myid < nr_out_prts_rest )  THEN
       pe_start_index = myid * nr_out_prts_on_this_pe + 1
       pe_end_index   = MIN( ( myid + 1 ) * nr_out_prts_on_this_pe, max_nr_out_prts)
    ELSE
       pe_start_index = nr_out_prts_rest * ( nr_out_prts_on_this_pe )                              &
                        + ( myid - nr_out_prts_rest ) * ( nr_out_prts_on_this_pe - 1 ) + 1
       pe_end_index   = MIN( pe_start_index + nr_out_prts_on_this_pe - 2, max_nr_out_prts )
    ENDIF

 END SUBROUTINE dop_setup_ireg_distribution


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define statistical particle NetCDF quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_statistics_init

    IMPLICIT NONE

    CHARACTER(LEN=3)  ::  suffix     !<
    CHARACTER(LEN=8)  ::  mode       !<
    CHARACTER(LEN=16) ::  long_name  !<

    CHARACTER(LEN=7), DIMENSION(nr_prt_quantities_statistical) :: dopts_label =                    &
          (/ 'tnpt   ', 'x_     ', 'y_     ', 'z_     ', 'z_abs  ', 'u      ',                     &
             'v      ', 'w      ', 'u"     ', 'v"     ', 'w"     ', 'npt_up ',                     &
             'w_up   ', 'w_down ', 'radius_', 'r_min  ', 'r_max  ', 'x*2    ',                     &
             'y*2    ', 'z*2    ', 'u*2    ', 'v*2    ', 'w*2    ', 'u"2    ',                     &
             'v"2    ', 'w"2    ' /)

    CHARACTER(LEN=7), DIMENSION(nr_prt_quantities_statistical) :: dopts_unit =                     &
          (/ '       ', 'm      ', 'm      ', 'm      ', 'm      ', 'm/s    ',                     &
             'm/s    ', 'm/s    ', 'm/s    ', 'm/s    ', 'm/s    ', '       ',                     &
             'm/s    ', 'm/s    ', 'm      ', 'm      ', 'm      ', 'm2     ',                     &
             'm2     ', 'm2     ', 'm2/s2  ', 'm2/s2  ', 'm2/s2  ', 'm2/s2  ',                     &
             'm2/s2  ', 'm2/s2  ' /)

    INTEGER, PARAMETER ::  global_id_in_file = -1  !<

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  j             !<
    INTEGER(iwp) ::  return_value  !<


    ALLOCATE( prt_quantity_statistical(0:number_of_particle_groups,nr_prt_quantities_statistical) )

    IF ( dop_individual_ts )  THEN
       mode = 'parallel'
    ELSE
       mode = 'serial'
    ENDIF

    DO  i = 1, nr_prt_quantities_statistical
       DO  j = 0, number_of_particle_groups

          IF ( j == 0 )  THEN
             suffix = ''
          ELSE
             WRITE ( suffix, '(''_'',I2.2)' )  j
          ENDIF

          prt_quantity_statistical(j,i)%name  = TRIM(dopts_label(i)) // TRIM(suffix)
          prt_quantity_statistical(j,i)%units = dopts_unit(i)

          CALL netcdf4_init_variable( TRIM( mode ), file_id,                                       &
                                      prt_quantity_statistical(j,i)%var_id,                        &
                                      prt_quantity_statistical(j,i)%name, 'real32',                &
                                      dimension_ids(2:2), 'master-only', return_value )

          CALL netcdf4_write_attribute( TRIM( mode ), file_id,                                     &
                                        prt_quantity_statistical(j,i)%var_id, 'units',             &
                                        value_char=TRIM( prt_quantity_statistical(j,i)%units ),    &
                                        return_value=return_value )

          IF ( j == 0 )  THEN
             long_name = TRIM( prt_quantity_statistical(j,i)%name )
          ELSE
             WRITE( long_name, '(A,1X,''PG'',I2.2)' )  TRIM( prt_quantity_statistical(j,i)%name ), j
          ENDIF
          CALL netcdf4_write_attribute( TRIM( mode ), file_id,                                     &
                                        prt_quantity_statistical(j,i)%var_id,                      &
                                        'long_name', value_char=TRIM( long_name ),                 &
                                        return_value=return_value )

          var_list = TRIM( var_list ) // TRIM( prt_quantity_statistical(j,i)%name ) // '; '

          ALLOCATE( prt_quantity_statistical(j,i)%value(nr_time_values) )

          IF ( number_of_particle_groups == 1 )  EXIT

       ENDDO
    ENDDO

!
!-- Global attributes.
    CALL netcdf4_write_attribute( TRIM( mode ), file_id, global_id_in_file, 'title',               &
                                  TRIM( run_description_header ), return_value=return_value )

    CALL netcdf4_write_attribute( TRIM( mode ), file_id, global_id_in_file, 'VAR_LIST',            &
                                  TRIM( var_list ), return_value=return_value )
!
!-- Time attributes.
    CALL netcdf4_write_attribute( TRIM( mode ), file_id, var_id%time, 'unit', 'seconds',           &
                                  return_value=return_value )

    CALL netcdf4_write_attribute( TRIM( mode ), file_id, var_id%time, 'standard_name', 'time',     &
                                  return_value=return_value )

    CALL netcdf4_write_attribute( TRIM( mode ), file_id, var_id%time, 'long_name', 'time',         &
                                  return_value=return_value )

    CALL netcdf4_write_attribute( TRIM( mode ), file_id, var_id%time, 'axis', 'T',                 &
                                  return_value=return_value )

 END SUBROUTINE dop_statistics_init

 END SUBROUTINE dop_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine collects and calculates particle data for timeseries output and initiates the
!> output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_data_output_ptseries

    INTEGER(iwp) ::  i     !<
#if defined( __parallel )
    INTEGER      ::  ierr  !<
#endif
    INTEGER(iwp) ::  inum  !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  jg    !<
    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  n     !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pts_value    !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pts_value_l  !<


    CALL cpu_log( log_point(36), 'data_output_ptseries', 'start' )

!
!-- Add time level to the time axis data.
    pts_time_index = pts_time_index + 1
    IF ( myid == 0  .AND.  root_model )  THEN
       time_axis_values(pts_time_index) = simulated_time
    ENDIF

    IF ( dop_individual_ts )  THEN
       CALL dop_write_tseries_data
    ENDIF

    ALLOCATE( pts_value(0:number_of_particle_groups,nr_prt_quantities_statistical),                &
              pts_value_l(0:number_of_particle_groups,nr_prt_quantities_statistical) )

    pts_value_l = 0.0_wp
    pts_value_l(:,16) = 9999999.9_wp    ! for calculation of minimum radius

!
!-- Calculate or collect the particle time series quantities for all particles and seperately for
!-- each particle group (if there is more than one group)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             number_of_particles = prt_count(k,j,i)
             IF (number_of_particles <= 0)  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                IF ( particles(n)%particle_mask )  THEN  ! Restrict analysis to active particles

                   pts_value_l(0,1)  = pts_value_l(0,1) + 1.0_wp                   ! total # of particles
                   pts_value_l(0,2)  = pts_value_l(0,2) +                                          &
                                       ( particles(n)%x - particles(n)%origin_x )  ! mean x
                   pts_value_l(0,3)  = pts_value_l(0,3) +                                          &
                                       ( particles(n)%y - particles(n)%origin_y )  ! mean y
                   pts_value_l(0,4)  = pts_value_l(0,4) +                                          &
                                       ( particles(n)%z - particles(n)%origin_z )  ! mean z
                   pts_value_l(0,5)  = pts_value_l(0,5) + particles(n)%z           ! mean z (absolute)
                   pts_value_l(0,6)  = pts_value_l(0,6) + particles(n)%speed_x     ! mean u
                   pts_value_l(0,7)  = pts_value_l(0,7) + particles(n)%speed_y     ! mean v
                   pts_value_l(0,8)  = pts_value_l(0,8) + particles(n)%speed_z     ! mean w
                   pts_value_l(0,9)  = pts_value_l(0,9)  + particles(n)%rvar1      ! mean sgsu
                   pts_value_l(0,10) = pts_value_l(0,10) + particles(n)%rvar2      ! mean sgsv
                   pts_value_l(0,11) = pts_value_l(0,11) + particles(n)%rvar3      ! mean sgsw
                   IF ( particles(n)%speed_z > 0.0_wp )  THEN
                      pts_value_l(0,12) = pts_value_l(0,12) + 1.0_wp                ! # of upward moving prts
                      pts_value_l(0,13) = pts_value_l(0,13) +  particles(n)%speed_z ! mean w upw.
                   ELSE
                      pts_value_l(0,14) = pts_value_l(0,14) + particles(n)%speed_z  ! mean w down
                   ENDIF
                   pts_value_l(0,15) = pts_value_l(0,15) + particles(n)%radius       ! mean rad
                   pts_value_l(0,16) = MIN( pts_value_l(0,16), particles(n)%radius ) ! minrad
                   pts_value_l(0,17) = MAX( pts_value_l(0,17), particles(n)%radius ) ! maxrad
!
!--                Repeat the same for the respective particle group
                   IF ( number_of_particle_groups > 1 )  THEN
                      jg = particles(n)%group

                      pts_value_l(jg,1)  = pts_value_l(jg,1) + 1.0_wp
                      pts_value_l(jg,2)  = pts_value_l(jg,2) +                                     &
                                           ( particles(n)%x  - particles(n)%origin_x )
                      pts_value_l(jg,3)  = pts_value_l(jg,3) +                                     &
                                           ( particles(n)%y  - particles(n)%origin_y )
                      pts_value_l(jg,4)  = pts_value_l(jg,4) +                                     &
                                           ( particles(n)%z  - particles(n)%origin_z )
                      pts_value_l(jg,5)  = pts_value_l(jg,5) + particles(n)%z
                      pts_value_l(jg,6)  = pts_value_l(jg,6) + particles(n)%speed_x
                      pts_value_l(jg,7)  = pts_value_l(jg,7) + particles(n)%speed_y
                      pts_value_l(jg,8)  = pts_value_l(jg,8) + particles(n)%speed_z
                      pts_value_l(jg,9)  = pts_value_l(jg,9)  + particles(n)%rvar1
                      pts_value_l(jg,10) = pts_value_l(jg,10) + particles(n)%rvar2
                      pts_value_l(jg,11) = pts_value_l(jg,11) + particles(n)%rvar3
                      IF ( particles(n)%speed_z > 0.0_wp )  THEN
                         pts_value_l(jg,12) = pts_value_l(jg,12) + 1.0_wp
                         pts_value_l(jg,13) = pts_value_l(jg,13) + particles(n)%speed_z
                      ELSE
                         pts_value_l(jg,14) = pts_value_l(jg,14) + particles(n)%speed_z
                      ENDIF
                      pts_value_l(jg,15) = pts_value_l(jg,15) + particles(n)%radius
                      pts_value_l(jg,16) = MIN( pts_value_l(jg,16), particles(n)%radius )
                      pts_value_l(jg,17) = MAX( pts_value_l(jg,17), particles(n)%radius )
                   ENDIF

                ENDIF

             ENDDO

          ENDDO
       ENDDO
    ENDDO


#if defined( __parallel )
!
!-- Sum values of the subdomains
!-- Particle output with dop is done only on the root model. Therefore MPI_COMM_WORLD is used to
!-  compute the sum across all nested models.
    inum = number_of_particle_groups + 1

    IF ( collective_wait )  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,1), pts_value(0,1), 15*inum, MPI_REAL, MPI_SUM,              &
                        MPI_COMM_WORLD, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,16), pts_value(0,16), inum, MPI_REAL, MPI_MIN,               &
                        MPI_COMM_WORLD, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,17), pts_value(0,17), inum, MPI_REAL, MPI_MAX,               &
                        MPI_COMM_WORLD, ierr )
#else
    pts_value(:,1:17) = pts_value_l(:,1:17)
#endif

!
!-- Normalize the above calculated quantities (except min/max values) with the total number of
!-- particles
    IF ( number_of_particle_groups > 1 )  THEN
       inum = number_of_particle_groups
    ELSE
       inum = 0
    ENDIF

    DO  j = 0, inum

       IF ( pts_value(j,1) > 0.0_wp )  THEN

          pts_value(j,2:15) = pts_value(j,2:15) / pts_value(j,1)
          IF ( pts_value(j,12) > 0.0_wp  .AND.  pts_value(j,12) < 1.0_wp )  THEN
             pts_value(j,13) = pts_value(j,13) / pts_value(j,12)
             pts_value(j,14) = pts_value(j,14) / ( 1.0_wp - pts_value(j,12) )
          ELSEIF ( pts_value(j,12) == 0.0_wp )  THEN
             pts_value(j,13) = -1.0_wp
          ELSE
             pts_value(j,14) = -1.0_wp
          ENDIF

       ENDIF

    ENDDO

!
!-- Calculate higher order moments of particle time series quantities, seperately for each particle
!-- group (if there is more than one group)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             number_of_particles = prt_count(k,j,i)
             IF (number_of_particles <= 0)  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                pts_value_l(0,18) = pts_value_l(0,18) + ( particles(n)%x -                         &
                                       particles(n)%origin_x - pts_value(0,2) )**2 ! x*2
                pts_value_l(0,19) = pts_value_l(0,19) + ( particles(n)%y -                         &
                                       particles(n)%origin_y - pts_value(0,3) )**2 ! y*2
                pts_value_l(0,20) = pts_value_l(0,20) + ( particles(n)%z -                         &
                                       particles(n)%origin_z - pts_value(0,4) )**2 ! z*2
                pts_value_l(0,21) = pts_value_l(0,21) + ( particles(n)%speed_x -                   &
                                                          pts_value(0,6) )**2      ! u*2
                pts_value_l(0,22) = pts_value_l(0,22) + ( particles(n)%speed_y -                   &
                                                          pts_value(0,7) )**2      ! v*2
                pts_value_l(0,23) = pts_value_l(0,23) + ( particles(n)%speed_z -                   &
                                                          pts_value(0,8) )**2      ! w*2
                pts_value_l(0,24) = pts_value_l(0,24) + ( particles(n)%rvar1 -                     &
                                                          pts_value(0,9) )**2      ! u"2
                pts_value_l(0,25) = pts_value_l(0,25) + ( particles(n)%rvar2 -                     &
                                                          pts_value(0,10) )**2     ! v"2
                pts_value_l(0,26) = pts_value_l(0,26) + ( particles(n)%rvar3 -                     &
                                                          pts_value(0,11) )**2  ! w"2
!
!--             Repeat the same for the respective particle group
                IF ( number_of_particle_groups > 1 )  THEN
                   jg = particles(n)%group

                   pts_value_l(jg,18) = pts_value_l(jg,18) + ( particles(n)%x -                    &
                                           particles(n)%origin_x - pts_value(jg,2) )**2
                   pts_value_l(jg,19) = pts_value_l(jg,19) + ( particles(n)%y -                    &
                                           particles(n)%origin_y - pts_value(jg,3) )**2
                   pts_value_l(jg,20) = pts_value_l(jg,20) + ( particles(n)%z -                    &
                                           particles(n)%origin_z - pts_value(jg,4) )**2
                   pts_value_l(jg,21) = pts_value_l(jg,21) + ( particles(n)%speed_x -              &
                                                               pts_value(jg,6) )**2
                   pts_value_l(jg,22) = pts_value_l(jg,22) + ( particles(n)%speed_y -              &
                                                               pts_value(jg,7) )**2
                   pts_value_l(jg,23) = pts_value_l(jg,23) + ( particles(n)%speed_z -              &
                                                               pts_value(jg,8) )**2
                   pts_value_l(jg,24) = pts_value_l(jg,24) + ( particles(n)%rvar1 -                &
                                                               pts_value(jg,9) )**2
                   pts_value_l(jg,25) = pts_value_l(jg,25) + ( particles(n)%rvar2 -                &
                                                               pts_value(jg,10) )**2
                   pts_value_l(jg,26) = pts_value_l(jg,26) + ( particles(n)%rvar3 -                &
                                                               pts_value(jg,11) )**2
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

#if defined( __parallel )
!
!-- Sum values of the subdomains
    inum = number_of_particle_groups + 1

    IF ( collective_wait )  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,18), pts_value(0,18), inum*9, MPI_REAL, MPI_SUM, MPI_COMM_WORLD,     &
                        ierr )
#else
    pts_value(:,18:26) = pts_value_l(:,18:26)
#endif

!
!-- Normalize the above calculated quantities with the total number of particles
    IF ( number_of_particle_groups > 1 )  THEN
       inum = number_of_particle_groups
    ELSE
       inum = 0
    ENDIF

    DO  j = 0, inum

       IF ( pts_value(j,1) > 0.0_wp )  THEN
          pts_value(j,18:26) = pts_value(j,18:26) / pts_value(j,1)
       ENDIF

    ENDDO

#if defined( __netcdf )
!
!-- Just collect statistical time series quantities to be output in NetCDF format at the end of the
!-- run.
    CALL dop_collect_statistics( pts_value )
#endif

    DEALLOCATE( pts_value, pts_value_l )

    CALL cpu_log( log_point(36), 'data_output_ptseries', 'stop' )

 END SUBROUTINE dop_data_output_ptseries


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Output of particle time series data at the selected time steps.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_write_tseries_data

    IMPLICIT NONE

    INTEGER(iwp)        ::  i              !<
    INTEGER(iwp)        ::  return_value   !< Return value data_output_netcdf4 .. routines
#if defined( __parallel )
    INTEGER(iwp)        ::  ierr           !< MPI error code
#endif
    INTEGER(iwp), DIMENSION(2) ::  bounds_origin  !<
    INTEGER(iwp), DIMENSION(2) ::  bounds_start   !<
    INTEGER(iwp), DIMENSION(2) ::  value_counts   !<

!
!-- Particles are output by root PEs only. Collect particle data from particles located in childs.
    IF ( nested_run )  THEN
       CALL dop_transfer_from_child_to_root
    ENDIF

!
!-- Now output by root model only.
    IF ( root_model )  THEN

!
!--    Output particle IDs for masked particles are set to -1 so that no time series data will be
!--    written.
       CALL dop_mark_masked_particles

       IF ( particle_id_file_found )  THEN
!
!--       Mark those particles scheduled for output via a particle ID file by setting their
!--       output ID to 1.
          CALL dop_mark_ids_scheduled_for_output_via_file
       ELSE
!
!--       Mark those particles scheduled for output via parameters pts_increment or pts_percentage.
          CALL dop_mark_ids_scheduled_for_output
       ENDIF

       CALL dop_count_particles_to_be_transfered_to_root

       bounds_origin    = 1
       bounds_start(1)  = io_start_index
       bounds_start(2)  = pts_time_index

       value_counts(1)  = io_end_index-io_start_index + 1
       value_counts(2)  = 1

       DO  i = 1, nr_prt_quantities_scheduled
#if defined( __netcdf4 )
          IF ( prt_quantity(i)%is_integer )  THEN
             out_buf_i(pe_start_index:pe_end_index) = NF90_FILL_INT
          ELSE
             out_buf_r(pe_start_index:pe_end_index) = NF90_FILL_REAL
          ENDIF
#endif
          CALL prt_mio%sm_node_barrier( )

          CALL cpu_log( log_point_s(99), 'dop_fill_buffers', 'start' )
!
!--       Copy the output particles into buffers (output buffer on root, transfer buffer on childs).
          CALL dop_fill_buffers( prt_quantity(i) )
#if defined( __parallel )
          IF ( nested_run  .AND.  particle_coupling )  THEN
             CALL prt_mio%sm_node_barrier( )
             CALL dop_child_var_to_out_buf( prt_quantity(i) )
          ENDIF
#endif
          CALL cpu_log( log_point_s(99), 'dop_fill_buffers', 'stop' )

          CALL prt_mio%sm_node_barrier( )

          CALL cpu_log( log_point_s(88), 'dop_get_remote_particles', 'start' )
#if defined( __parallel )
          CALL  dop_get_remote_particle( prt_quantity(i)%is_integer )
#endif
          CALL cpu_log( log_point_s(88), 'dop_get_remote_particles', 'stop' )

          CALL cpu_log( log_point_s(89), 'particle NetCDF output', 'start' )
          CALL prt_mio%sm_node_barrier( )
          IF ( prt_mio%iam_io_pe )  THEN
             IF ( prt_quantity(i)%is_integer )  THEN
                CALL netcdf4_write_variable( 'parallel', file_id, prt_quantity(i)%var_id,          &
                                             bounds_start, value_counts, bounds_origin,            &
                                             'collective', values_int32_1d=out_buf_i,              &
                                             return_value=return_value )
             ELSE
                CALL netcdf4_write_variable( 'parallel', file_id, prt_quantity(i)%var_id,          &
                                             bounds_start, value_counts, bounds_origin,            &
                                             'collective', values_real32_1d=out_buf_r,             &
                                             return_value=return_value )
             ENDIF
          ENDIF
          CALL prt_mio%sm_node_barrier( )
          CALL cpu_log( log_point_s(89), 'particle NetCDF output', 'stop' )

#if defined( __parallel )
!
!--       This barrier is required, although the reason for it is unclear.
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       CALL deallocate_and_free

    ENDIF

 END SUBROUTINE dop_write_tseries_data


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collect statistical particle data at the selected output time steps.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_collect_statistics( pts_value )

    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(0:,:) ::  pts_value  !<

    INTEGER(iwp)        ::  i           !<
    INTEGER(iwp)        ::  j           !<


    IF ( myid == 0  .AND.  root_model )  THEN

       DO  i = 1, nr_prt_quantities_statistical
          DO  j = 0, number_of_particle_groups
             prt_quantity_statistical(j,i)%value(pts_time_index) = pts_value(j,i)
             IF ( number_of_particle_groups == 1 )  EXIT
          ENDDO
       ENDDO

    ENDIF

 END SUBROUTINE dop_collect_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Close all MPI windows for one-sided particle data exchange.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_finalize

    IMPLICIT NONE

#if defined( __parallel )
    INTEGER(iwp) :: ierr          !< MPI error code
#endif
    INTEGER(iwp) :: return_value  !< return value for data_output_netcdf4 routines

#if defined( __parallel )
    IF ( nested_run  .AND.  particle_coupling  .AND.  dop_individual_ts )  THEN
       CALL MPI_WIN_FREE (cnpo%rem_out_win, ierr)
    ENDIF
#endif

    IF ( root_model )  THEN

       IF ( dop_individual_ts )  THEN

#if defined( __parallel )
          IF ( win_prt_i /= -1 )  THEN
             CALL prt_mio%sm_free_shared( win_prt_i )
          ENDIF
          IF ( win_prt_r /= -1 )  THEN
             CALL prt_mio%sm_free_shared( win_prt_r )
          ENDIF
#endif
          IF ( file_id /= -1  .AND.  prt_mio%iam_io_pe )  THEN
             CALL netcdf4_finalize( 'parallel', file_id, return_value )
             file_id = -1
          ENDIF
       ENDIF

       CALL dop_statistics_and_finalize
    ENDIF

 END SUBROUTINE dop_finalize


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> output of earlier colected statistical particle data and finalize statistical output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_statistics_and_finalize

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  j             !<
    INTEGER(iwp) ::  lo_var_id     !<
    INTEGER(iwp) ::  return_value  !< return value data_output_netcdf4 routines
    INTEGER(iwp) ::  var_len       !<

    INTEGER(iwp), DIMENSION(1) ::  bounds_origin  !<
    INTEGER(iwp), DIMENSION(1) ::  bounds_start   !<
    INTEGER(iwp), DIMENSION(1) ::  value_counts   !<

    REAL(sp), POINTER, DIMENSION(:), CONTIGUOUS ::  outbuf  !<


    IF ( myid == 0 )  THEN

       bounds_start  = 1
       bounds_origin = 1

       IF ( file_id == -1 )  THEN
          CALL netcdf4_open_file( 'serial', TRIM( file_name ), file_id, return_value,              &
                                  reopen=.TRUE. )
          CALL netcdf4_inquire_dimension( file_id, did%time, return_value,                         &
                                  dimension_length = var_len)
       ELSE
          var_len = pts_time_index
       ENDIF

       value_counts(1) = var_len
       DO  i = 1, nr_prt_quantities_statistical
          DO  j = 0, number_of_particle_groups

             CALL netcdf4_inq_varid( file_id, prt_quantity_statistical(j,i)%name, lo_var_id,       &
                                     return_value)
             outbuf => prt_quantity_statistical(j,i)%value(1:var_len)
             CALL netcdf4_write_variable( 'serial', file_id, lo_var_id, bounds_start,              &
                                          value_counts, bounds_origin, 'master-only',              &
                                          values_real32_1d=outbuf, return_value=return_value )

             IF ( number_of_particle_groups == 1 )  EXIT

          ENDDO
       ENDDO

       CALL netcdf4_inq_varid( file_id, 'time', lo_var_id, return_value )
!
!--    ATTENTION:The following statement works with gfortran, but not with the Intel 19.0.5
!       CALL  netcdf4_write_variable( 'serial', file_id, lo_var_id,           &
!          bounds_start, value_counts, bounds_origin, 'master-only',          &
!          values_real32_1d=time_axis_values(1:var_len), return_value=return_value )
!
!--    Therefore, the following workaround:
       outbuf => time_axis_values(1:var_len)
       CALL netcdf4_write_variable( 'serial', file_id, lo_var_id, bounds_start, value_counts,      &
                                    bounds_origin, 'master-only', values_real32_1d=outbuf,         &
                                    return_value=return_value )

       CALL netcdf4_finalize( 'serial', file_id, return_value )

    ENDIF

 END SUBROUTINE dop_statistics_and_finalize


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set all output particle IDs to undefined (-1).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE set_output_particle_ids_undefined

    IMPLICIT NONE

    INTEGER(iwp) :: i  !<
    INTEGER(iwp) :: j  !<
    INTEGER(iwp) :: k  !<
    INTEGER(iwp) :: n  !<

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             DO  n = 1, SIZE( grid_particles(k,j,i)%particles )
                grid_particles(k,j,i)%particles(n)%io_id = -1
             ENDDO
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE set_output_particle_ids_undefined


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Count particles scheduled for output.
!> Here pts_increment and pts_percentage are used to select output particles.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE count_output_particles( pcount, nr_prts_in_gridbox )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(OUT) ::  pcount  !<

    INTEGER(iwp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  nr_prts_in_gridbox  !< number of particles in grid box (k,j,i)

    INTEGER(iwp) ::  bcount  !< basis for particle count in case of pts_percentage
    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  k       !<
    INTEGER(iwp) ::  n       !<
    INTEGER(iwp) ::  n_all   !< count all particles for MOD function

#if defined( __parallel )
    INTEGER(iwp) ::  ierr    !< MPI error code
#endif

    REAL(dp) ::  fcount  !<
    REAL(dp) ::  finc    !<


!
!-- Compute number of particles on local PE and broadcast it to all PEs.
    nr_prts_on_pe = 0_idp
    nr_prts_on_pe(myidy,myidx) = SUM( nr_prts_in_gridbox )
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nr_prts_on_pe, SIZE( nr_prts_on_pe ), MPI_INTEGER8, MPI_SUM, &
                        comm2d, ierr )
#endif
!
!-  Compute the number of particles of previous PE columns along y.
    nr_prts_prev_pes_ycolumn(0) = 0
    DO  i = 1, npex-1
       nr_prts_prev_pes_ycolumn(i) = nr_prts_prev_pes_ycolumn(i-1) + SUM( nr_prts_on_pe(:,i-1) )
    ENDDO
!
!-- For each grid point along x (Eulerian coordinate), calculate the sum of particles over the
!-- respective yz-cross-section of all PE subdomains with the same virtual PE index y.
    nr_prts_along_yz_on_pe = 0_idp
    DO  i = nxl, nxr
       nr_prts_along_yz_on_pe(i,myidy) = SUM( nr_prts_in_gridbox(:,:,i) )
    ENDDO
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nr_prts_along_yz_on_pe, SIZE( nr_prts_along_yz_on_pe ),      &
                        MPI_INTEGER8, MPI_SUM, comm1dy, ierr )
#endif
!
!-- Compute global particle counter at begin of every yz-cross-section.
    start_prt_count_yz = 0_dp
!
!-- First: add particles for i < given i and summarized over all PEs along y (which gives the
!-- number of particles of yz-cross-sections of the total domain with i < given i.
    DO  i = nxl+1, nxr
       start_prt_count_yz(i) = SUM( nr_prts_along_yz_on_pe(nxl:i-1,:) )
    ENDDO
!
!-- Second: for given i, add particles of yz-cross sections of all PEs smaller than myidy.
    DO  i = nxl, nxr
       DO  j = 1, myidy
          start_prt_count_yz(i) = start_prt_count_yz(i) + nr_prts_along_yz_on_pe(i,j-1)
       ENDDO
    ENDDO
!
!-- Third: add all particles from PEs with an id < myidx.
    DO  i = nxl, nxr
       start_prt_count_yz(i) = start_prt_count_yz(i) + nr_prts_prev_pes_ycolumn(myidx)
    ENDDO
!
!-- Compute number of output particles, depending on if they are selected via an increment or via
!-- percentage.
    pcount = 0
    IF ( pts_increment > 0 )  THEN

       DO  i = nxl, nxr
          n_all = start_prt_count_yz(i)
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, nr_prts_in_gridbox(k,j,i)
                   IF ( MOD( n_all, pts_increment ) == 0 )  THEN
                      pcount = pcount + 1
                   ENDIF
                   n_all = n_all + 1
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( pts_percentage > 0.0_wp )  THEN

       finc   = pts_percentage / 100.0_idp
       DO  i = nxl, nxr

          n_all  = start_prt_count_yz(i)
          fcount = REAL( n_all, dp ) * finc + eps
          bcount = INT( REAL( n_all-1, dp ) * finc ) + 1
          IF ( n_all == 0 )  bcount = 0
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, nr_prts_in_gridbox(k,j,i)
                   IF ( fcount >= REAL( bcount, dp ) )  THEN
                      pcount = pcount + 1
                      bcount = bcount + 1
                   ENDIF
                   fcount = fcount + finc
                   n_all = n_all + 1
                ENDDO
             ENDDO
          ENDDO

       ENDDO

    ENDIF

 END SUBROUTINE count_output_particles


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read particle list from the particle-ID file. These particles are scheduled for time series
!> output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_read_ids_scheduled_for_ts_output( nr_out_prts )

    IMPLICIT NONE

    INTEGER(idp) ::  dummy      !<
    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  io_status  !<

    INTEGER(iwp), INTENT(OUT) ::  nr_out_prts  !< number of particles scheduled for time series output


    OPEN( 90, FILE='PARTICLE_IDS' )

!
!-- First stride: count number of particle-IDs given in file.
    nr_out_prts = 0
    io_status = 0
    DO  WHILE( io_status == 0 )
       READ( 90, *, IOSTAT=io_status )  dummy
       nr_out_prts = nr_out_prts + 1
       IF ( io_status > 0 )  THEN
          WRITE( message_string, '(A,1X,I5)' )  'error while reading file PARTICLE_IDS in line',   &
                                                nr_out_prts
          CALL message( 'dop_read_ids_scheduled_for_ts_output', 'LPM0030', 3, 2, 0, 6, 1 )
       ENDIF
    ENDDO
!
!-- Subtract 1 for end of file read.
    nr_out_prts = nr_out_prts - 1

    ALLOCATE( particle_id_scheduled_for_output(nr_out_prts) )

    REWIND( 90 )
!
!-- Second stride, read particle IDs that are scheduled for time series output.
    DO  i = 1, nr_out_prts
       READ( 90, * )  particle_id_scheduled_for_output(i)
    ENDDO

    CLOSE( 90 )

 END SUBROUTINE dop_read_ids_scheduled_for_ts_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set output particle number for particles that are scheduled for output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_set_io_ids

    IMPLICIT NONE

    INTEGER(iwp) ::  bcount  !< basis for particle count in case of pts_percentage
    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  k       !<
    INTEGER(iwp) ::  n       !<
    INTEGER(iwp) ::  n_all   !< count all particles for MOD function
    INTEGER(iwp) ::  pcount  !< local particle count in case of pts_percentage

    REAL(dp) ::  fcount  !< partical progress in %/100
    REAL(dp) ::  finc    !< increment of particle


!
!-- Set particle output IDs depending on if they are selected via an increment or via percentage.
    pcount = 0
    IF ( pts_increment > 0 )  THEN

       DO  i = nxl, nxr
          n_all  = start_prt_count_yz(i)
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                DO  n = 1, prt_count(k,j,i)
                   IF ( MOD( n_all, pts_increment ) == 0 )  THEN
                      grid_particles(k,j,i)%particles(n)%io_id = (n_all/pts_increment) + 1
                   ELSE
                      grid_particles(k,j,i)%particles(n)%io_id = -2
                   ENDIF
                   n_all = n_all + 1
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( pts_percentage > 0.0_wp )  THEN

       DO  i = nxl, nxr

          n_all  = start_prt_count_yz(i)
          finc   = pts_percentage / 100.0_wp
          fcount = REAL( n_all, dp ) * finc + eps
          bcount = INT( REAL( n_all-1, dp ) * finc ) + 1
          IF ( n_all == 0 )  bcount = 0

          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, prt_count(k,j,i)
!
!--                Particles to be output defined as percentage of all particles.
!--                If e.g. finc = 0.2, then every 5th particle will be output.
                   IF ( fcount >= REAL( bcount, dp ) )  THEN
                      grid_particles(k,j,i)%particles(n)%io_id = INT( fcount ) + 1
                      bcount = bcount + 1
                   ELSE
                      grid_particles(k,j,i)%particles(n)%io_id = -2
                   ENDIF
                   fcount = fcount + finc
                   n_all = n_all + 1
                ENDDO
             ENDDO
          ENDDO

       ENDDO

    ENDIF

 END SUBROUTINE dop_set_io_ids


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Find/identify those particles that are defined in the given particle-ID list, and set their
!> output ID to 1.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_mark_ids_scheduled_for_output_via_file

    IMPLICIT NONE

#if defined( __parallel )
    INTEGER(iwp) :: ierr  !< MPI error code
#endif
    INTEGER(iwp) :: i                  !< loop index
    INTEGER(iwp) :: ioid               !< loop index
    INTEGER(iwp) :: j                  !< loop index
    INTEGER(iwp) :: k                  !< loop index
    INTEGER(iwp) :: n                  !< loop index
    INTEGER(iwp) :: nr_out_prts_found  !< counter for number of particle-IDs found in file


    nr_out_prts_found = 0
!
!>  TODO: If there is a long particle output list, it may become necessary to optimize the following
!>        loop for performance reasons.
!
!-- Decode the particle id in i,j,k,n.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             DO  n = 1, prt_count(k,j,i)
                DO  ioid = 1, SIZE( particle_id_scheduled_for_output )
                   IF ( grid_particles(k,j,i)%particles(n)%id ==                                   &
                        particle_id_scheduled_for_output(ioid) )                                   &
                   THEN
                      grid_particles(k,j,i)%particles(n)%io_id = ioid
                      nr_out_prts_found = nr_out_prts_found + 1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

#if defined( __parallel )
    CALL MPI_ALLREDUCE( nr_out_prts_found, nr_out_prts, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    nr_out_prts = nr_out_prts_found
#endif

 END SUBROUTINE dop_mark_ids_scheduled_for_output_via_file


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Open the NetCDF File DATA_1D_PTS_NETCDF. Define dimensions and output quantitiess,
!> and write the constant quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_netcdf_setup

    IMPLICIT NONE

    INTEGER, PARAMETER ::  global_id_in_file = -1  !<

    INTEGER(iwp) ::  ind  !< quantity index
    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  return_value  !<

    LOGICAL ::  const_flag  !<


    IF ( prt_mio%iam_io_pe )  THEN
       CALL netcdf4_init_module( "", prt_mio%comm_io, 0, 9, .TRUE., -1 )

       file_name = 'DATA_1D_PTS_NETCDF' // TRIM( coupling_char )
#if defined( __parallel )
       CALL netcdf4_open_file( 'parallel', TRIM( file_name ), file_id, return_value )
#else
       CALL netcdf4_open_file( 'serial', TRIM( file_name ), file_id, return_value )
#endif
!
!--    Global attributes.
       CALL netcdf4_write_attribute( 'parallel', file_id, global_id_in_file, 'comment',            &
                                     'Particle output created by PALM module data_output_particle',&
                                     return_value=return_value )

       CALL netcdf4_write_attribute( 'parallel', file_id, global_id_in_file,                       &
                                     'initial_nr_particles',                                       &
                                     value_int32=nr_out_prts, return_value=return_value )
!
!--    Define dimensions.
       CALL netcdf4_init_dimension( 'parallel', file_id, did%prt, var_id%prt, 'prt', 'int32' ,     &
                                     max_nr_out_prts, 'master-only', return_value )


       CALL netcdf4_write_attribute( 'parallel', file_id, var_id%prt, 'units',                     &
                                     '      ', return_value=return_value )

       CALL netcdf4_write_attribute( 'parallel', file_id, var_id%prt, 'long_name',                 &
                                     'particle_number', return_value=return_value )

       CALL netcdf4_write_attribute( 'parallel', file_id, var_id%prt, 'axis',                      &
                                     'N', return_value=return_value )

       IF ( unlimited_dimension )  THEN
          CALL netcdf4_init_dimension( 'parallel', file_id, did%time, var_id%time, 'time','real32',&
                                       -1, 'master-only', return_value )
       ELSE
          CALL netcdf4_init_dimension( 'parallel', file_id, did%time, var_id%time, 'time','real32',&
                                       nr_time_values, 'master-only', return_value )
       ENDIF
    ENDIF

!
!-- Quantities which do not depend on time (i.e. without time axis).
!-- These quantities will always but only be written once at the beginning of the output.
    dimension_ids(1) = did%prt
    ind = 1
    prt_quantity_const(ind)%name  = 'origin_x'
    prt_quantity_const(ind)%units = 'meter'

    IF ( prt_mio%iam_io_pe )  THEN
       CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity_const(ind)%var_id,            &
                                   prt_quantity_const(ind)%name, 'real32', dimension_ids(1:1),     &
                                   'collective', return_value )

       CALL netcdf4_write_attribute( 'parallel', file_id, prt_quantity_const(ind)%var_id, 'units', &
                                     value_char=TRIM( prt_quantity_const(ind)%units ),             &
                                     return_value=return_value )

       var_list = TRIM( var_list ) // TRIM( prt_quantity_const(ind)%name ) // '; '
    ENDIF

    ind = ind + 1
    prt_quantity_const(ind)%name  = 'origin_y'
    prt_quantity_const(ind)%units = 'meter'

    IF ( prt_mio%iam_io_pe )  THEN
       CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity_const(ind)%var_id,            &
                                   prt_quantity_const(ind)%name, 'real32', dimension_ids(1:1),     &
                                   'collective', return_value )

       CALL netcdf4_write_attribute( 'parallel', file_id, prt_quantity_const(ind)%var_id, 'units', &
                                     value_char=TRIM( prt_quantity_const(ind)%units ),             &
                                     return_value=return_value )

       var_list = TRIM( var_list ) // TRIM( prt_quantity_const(ind)%name ) // '; '
    ENDIF

    ind = ind + 1
    prt_quantity_const(ind)%name  = 'origin_z'
    prt_quantity_const(ind)%units = 'meter'

    IF ( prt_mio%iam_io_pe )  THEN
       CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity_const(ind)%var_id,            &
                                   prt_quantity_const(ind)%name, 'real32', dimension_ids(1:1),     &
                                   'collective', return_value )

       CALL netcdf4_write_attribute( 'parallel', file_id, prt_quantity_const(ind)%var_id, 'units', &
                                     value_char=TRIM( prt_quantity_const(ind)%units ),             &
                                     return_value=return_value )

       var_list = TRIM( var_list ) // TRIM( prt_quantity_const(ind)%name ) // '; '
    ENDIF

!
!-- These quantities are written if name ends with '_const', which means that they do not depend
!-- on time.
    DO  i = 1, SIZE( data_output_pts )
       const_flag = ( INDEX( TRIM( data_output_pts(i) ), '_const' ) > 0 )
       IF ( LEN( TRIM( data_output_pts(i) ) ) > 0  .AND.  const_flag )  THEN
          ind = ind + 1
          prt_quantity_const(ind)%name  = TRIM( data_output_pts(i) )

          SELECT CASE ( TRIM( prt_quantity_const(ind)%name ) )

             CASE ( 'radius_const' )
                prt_quantity_const(ind)%units = 'meter'
                prt_quantity_const(ind)%is_integer = .FALSE.

             CASE ( 'aux1_const' )
                prt_quantity_const(ind)%units = 'depend_on_setup'
                prt_quantity_const(ind)%is_integer = .FALSE.

             CASE ( 'aux2_const' )
                prt_quantity_const(ind)%units = 'depend_on_setup'
                prt_quantity_const(ind)%is_integer = .FALSE.

             CASE ( 'rvar1_const' )
                prt_quantity_const(ind)%units = 'depend_on_setup'
                prt_quantity_const(ind)%is_integer = .FALSE.

             CASE ( 'rvar2_const' )
                prt_quantity_const(ind)%units = 'depend_on_setup'
                prt_quantity_const(ind)%is_integer = .FALSE.

             CASE ( 'rvar3_const' )
                prt_quantity_const(ind)%units = 'depend_on_setup'
                prt_quantity_const(ind)%is_integer = .FALSE.

          END SELECT

          IF ( prt_mio%iam_io_pe )  THEN
             IF ( prt_quantity_const(ind)%is_integer )  THEN
                CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity_const(ind)%var_id,   &
                                            prt_quantity_const(ind)%name, 'int32',                 &
                                            dimension_ids(1:1), 'collective', return_value )
             ELSE
                CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity_const(ind)%var_id,   &
                                             prt_quantity_const(ind)%name, 'real32',               &
                                             dimension_ids(1:1), 'collective', return_value )
             ENDIF

             CALL netcdf4_write_attribute( 'parallel', file_id, prt_quantity_const(ind)%var_id,    &
                                           'units',                                                &
                                           value_char=TRIM( prt_quantity_const(ind)%units ),       &
                                           return_value=return_value )

             var_list = TRIM( var_list ) // TRIM( prt_quantity_const(ind)%name ) // '; '
          ENDIF

       ENDIF
    ENDDO
!
!-- Save the number of constant particle quantities scheduled for output.
    nr_prt_quantities_const_scheduled = ind

!
!-- Variables time axis.
!-- Following quantities will always be written at each output timestep.
    dimension_ids(1) = did%prt
    dimension_ids(2) = did%time

    ind = 0

    DO  i = 1, SIZE( data_output_pts )
       const_flag = ( INDEX( TRIM( data_output_pts(i) ), '_const' ) > 0 )
       IF ( LEN( TRIM( data_output_pts(i) ) ) > 0  .AND.  .NOT. const_flag )  THEN
          ind = ind + 1
          prt_quantity(ind)%name  = TRIM( data_output_pts(i) )

          SELECT CASE ( TRIM( prt_quantity(ind)%name) )

             CASE ( 'id' )
                prt_quantity(ind)%name  = TRIM( data_output_pts(i) ) // '_low'
                prt_quantity(ind)%units = '      '
                prt_quantity(ind)%is_integer = .TRUE.
                IF ( prt_mio%iam_io_pe )  THEN
                   CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity(ind)%var_id,      &
                                               prt_quantity(ind)%name, 'int32',                    &
                                               dimension_ids(1:2), 'collective', return_value )
                   CALL netcdf4_write_attribute( 'parallel', file_id, prt_quantity(ind)%var_id,    &
                                                 'units',                                          &
                                                 value_char=TRIM( prt_quantity(ind)%units ),       &
                                                 return_value=return_value )
                ENDIF

                ind = ind + 1
                prt_quantity(ind)%name  = TRIM( data_output_pts(i) ) // '_high'
                prt_quantity(ind)%units = '      '
                prt_quantity(ind)%is_integer = .TRUE.

             CASE ( 'particle_io_id' )
                prt_quantity(ind)%units = '      '
                prt_quantity(ind)%is_integer = .TRUE.

             CASE ( 'class' )
                prt_quantity(ind)%units = '      '
                prt_quantity(ind)%is_integer = .TRUE.

             CASE ( 'group' )
                prt_quantity(ind)%units = '      '
                prt_quantity(ind)%is_integer = .TRUE.

             CASE ( 'x' )
                prt_quantity(ind)%units = 'meter'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'y' )
                prt_quantity(ind)%units = 'meter'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'z' )
                prt_quantity(ind)%units = 'meter'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'speed_x' )
                prt_quantity(ind)%units = 'm/s'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'speed_y' )
                prt_quantity(ind)%units = 'm/s'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'speed_z' )
                prt_quantity(ind)%units = 'm/s'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'radius' )
                prt_quantity(ind)%units = 'meter'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'age' )
                prt_quantity(ind)%units = 'sec'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'age_m' )
                prt_quantity(ind)%units = 'sec'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'dt_sum' )
                prt_quantity(ind)%units = 'sec'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'e_m' )
                prt_quantity(ind)%units = 'Ws'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE( 'weight_factor' )
                prt_quantity(ind)%units = 'factor'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'aux1' )
                prt_quantity(ind)%units = 'depend_on_setup'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'aux2' )
                prt_quantity(ind)%units = 'depend_on_setup'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'rvar1' )
                prt_quantity(ind)%units = 'depend_on_setup'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'rvar2' )
                prt_quantity(ind)%units = 'depend_on_setup'
                prt_quantity(ind)%is_integer = .FALSE.

             CASE ( 'rvar3' )
                prt_quantity(ind)%units = 'depend_on_setup'
                prt_quantity(ind)%is_integer = .FALSE.

          END SELECT

          IF ( prt_mio%iam_io_pe )  THEN
             IF ( prt_quantity(ind)%is_integer )  THEN
                CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity(ind)%var_id,         &
                                            prt_quantity(ind)%name, 'int32', dimension_ids(1:2),   &
                                            'collective', return_value )
             ELSE
                CALL netcdf4_init_variable( 'parallel', file_id, prt_quantity(ind)%var_id,         &
                                            prt_quantity(ind)%name, 'real32', dimension_ids(1:2),  &
                                            'collective', return_value )
             ENDIF

             CALL netcdf4_write_attribute( 'parallel', file_id, prt_quantity(ind)%var_id, 'units', &
                                           value_char=TRIM( prt_quantity(ind)%units ),             &
                                           return_value=return_value )


             var_list = TRIM( var_list ) // TRIM( prt_quantity(ind)%name ) // '; '
          ENDIF

       ENDIF
    ENDDO
!
!-- Store the number of particle quantities scheduled for output.
    nr_prt_quantities_scheduled = ind

    IF ( prt_mio%iam_io_pe )  THEN
       CALL netcdf4_stop_file_header_definition( 'parallel', file_id, return_value )
    ENDIF

    CALL dop_write_particle_axis

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the particle coordinate axis (numbered ny particle I/O ID).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_write_particle_axis

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  return_value  !< return value of data_output_netcdf4... routines

    INTEGER,DIMENSION(1) ::  bounds_origin  !<
    INTEGER,DIMENSION(1) ::  bounds_start   !<
    INTEGER,DIMENSION(1) ::  value_counts   !<

    INTEGER, POINTER, CONTIGUOUS, DIMENSION(:) ::  prt_val  !<


    bounds_origin = 1
    bounds_start(1) = 1

    IF ( myid == 0 )  THEN

       ALLOCATE( prt_val(max_nr_out_prts) )
       DO  i = 1, max_nr_out_prts
          prt_val(i) = i
       ENDDO
       value_counts(1) = max_nr_out_prts

       CALL netcdf4_write_variable( 'parallel', file_id, var_id%prt, bounds_start, value_counts,   &
                                    bounds_origin, 'master-only', values_int32_1d=prt_val,         &
                                    return_value=return_value )
       DEALLOCATE( prt_val )

    ENDIF

 END SUBROUTINE dop_write_particle_axis

 END SUBROUTINE dop_netcdf_setup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Open the NetCDF File DATA_1D_PTS_NETCDF for output of particle statistics only.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_netcdf_setup_statistic_only

    IMPLICIT NONE

    INTEGER(iwp) ::  return_value  !<


    IF ( myid == 0 )  THEN

       CALL netcdf4_init_module( "", comm2d, 0, 9, .TRUE., -1 )

       file_name = 'DATA_1D_PTS_NETCDF' // TRIM( coupling_char )

       CALL netcdf4_open_file( 'serial', TRIM( file_name ), file_id, return_value )

       CALL netcdf4_init_dimension( 'serial', file_id, did%time, var_id%time, 'time','real32',     &
                                    -1, 'master-only', return_value )
       dimension_ids(2) = did%time

       CALL netcdf4_stop_file_header_definition( 'serial', file_id, return_value )

    ENDIF

 END SUBROUTINE dop_netcdf_setup_statistic_only


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write those particle quantities that do not depend on time.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_write_constant_quantities

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  return_value  !<

    INTEGER(iwp),DIMENSION(1) ::  bounds_origin  !<
    INTEGER(iwp),DIMENSION(1) ::  bounds_start   !<
    INTEGER(iwp),DIMENSION(1) ::  value_counts   !<


    bounds_origin(1) = 1
    bounds_start(1)  = io_start_index
    value_counts(1)  = io_end_index - io_start_index + 1


    DO  i = 1, nr_prt_quantities_const_scheduled

#if defined( __netcdf4 )
       IF ( prt_quantity_const(i)%is_integer )  THEN
          out_buf_i(pe_start_index:pe_end_index) = NF90_FILL_INT
       ELSE
          out_buf_r(pe_start_index:pe_end_index) = NF90_FILL_REAL
       ENDIF
#endif
       CALL prt_mio%sm_node_barrier( )
!
!--    Copy the output particles into buffers (output buffer on root, transfer buffer on childs).
       CALL dop_fill_buffers( prt_quantity_const(i) )

#if defined( __parallel )
       IF ( nested_run  .AND.  particle_coupling )  THEN
          CALL prt_mio%sm_node_barrier( )
          CALL dop_child_var_to_out_buf( prt_quantity_const(i) )
       ENDIF
#endif

       CALL prt_mio%sm_node_barrier( )
       CALL dop_get_remote_particle( prt_quantity_const(i)%is_integer )

       CALL prt_mio%sm_node_barrier( )
       IF ( prt_mio%iam_io_pe )  THEN
          IF ( prt_quantity_const(i)%is_integer )  THEN
             CALL netcdf4_write_variable( 'parallel', file_id, prt_quantity_const(i)%var_id,       &
                                          bounds_start, value_counts, bounds_origin,               &
                                          'collective', values_int32_1d=out_buf_i,                 &
                                          return_value=return_value )
          ELSE
             CALL  netcdf4_write_variable( 'parallel', file_id, prt_quantity_const(i)%var_id,      &
                                           bounds_start, value_counts, bounds_origin,              &
                                           'collective', values_real32_1d=out_buf_r,               &
                                           return_value=return_value )
          ENDIF
       ENDIF
       CALL prt_mio%sm_node_barrier( )
    ENDDO

 END SUBROUTINE dop_write_constant_quantities


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prevent time series output of particles that have been masked by setting their output ID to -1.
!> Anyhow, output is made for those particles, but only fill values are written to file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_mark_masked_particles

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  n  !<


!
!-- Mark all particles with particle_mask == .FALSE. by giving them an output ID of -1.
!-- Fill values will be written for these particles later.
    nr_out_prts_remote = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             DO  n = 1, prt_count(k,j,i)
                IF ( .NOT. grid_particles(k,j,i)%particles(n)%particle_mask )  THEN
                   grid_particles(k,j,i)%particles(n)%io_id = -1
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE dop_mark_masked_particles


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Count number of newly generated particles that are scheduled for time series output and give
!> them a contiguous I/O-ID.
!> The condition for a new particle is particle_mask = .TRUE. and a particle I/O ID of -1.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_mark_ids_scheduled_for_output

    IMPLICIT NONE

    INTEGER(iwp) ::  bcount               !< basis for particle count in case of pts_percentage
    INTEGER(iwp) ::  i                    !<
#if defined( __parallel )
    INTEGER(iwp) ::  ierr                 !< MPI error code
#endif
    INTEGER(iwp) ::  j                    !<
    INTEGER(iwp) ::  k                    !<
    INTEGER(iwp) ::  n                    !<
    INTEGER(iwp) ::  n_all                !< particle counter for MOD function
    INTEGER(iwp) ::  nr_out_prts_new      !< total number of new output particles over all PEs
    INTEGER(iwp) ::  particle_io_id       !<
    INTEGER(iwp) ::  start_new_numbering  !<

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  nr_out_prts_new_pe  !< number of new output particles per PE

    LOGICAL ::  new_particles_found  !< switch to decide, if new particles are found

    REAL(dp) ::  fcount  !<
    REAL(dp) ::  finc    !<


!
!>  TODO:
!>  For performance reasons, this subroutine may be combined later with dop_mark_masked_particles.
    new_prt_count     = 0
    nr_out_prts_new_pe(:) = 0
!
!-- Count particles depending on if they are selected via an increment or via percentage.
    IF ( pts_increment > 0 )  THEN
!
!--    Count newly generated particles (all, not only output particles).
       new_particles_found = .FALSE.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, prt_count(k,j,i)
                   IF ( grid_particles(k,j,i)%particles(n)%particle_mask )  THEN
                      IF ( grid_particles(k,j,i)%particles(n)%io_id == -1 )  THEN
                         new_prt_count(k,j,i) = new_prt_count(k,j,i) + 1
                         new_particles_found = .TRUE.
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, new_particles_found, 1, MPI_INTEGER8, MPI_LOR, comm2d,    &
                           ierr )
#endif

       IF ( new_particles_found )  THEN
          CALL count_output_particles( nr_out_prts_new_pe(myid), new_prt_count )
       ENDIF

    ELSEIF ( pts_percentage > 0.0_wp )  THEN

       new_particles_found = .FALSE.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, prt_count(k,j,i)
                   IF ( grid_particles(k,j,i)%particles(n)%particle_mask )  THEN
                      IF ( grid_particles(k,j,i)%particles(n)%io_id == -1 )  THEN
                         new_prt_count(k,j,i) = new_prt_count(k,j,i) + 1
                         new_particles_found = .TRUE.
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, new_particles_found, 1, MPI_INTEGER8, MPI_LOR, comm2d,    &
                           ierr )
#endif

       IF ( new_particles_found )  THEN
          CALL count_output_particles( nr_out_prts_new_pe(myid), new_prt_count )
       ENDIF

    ENDIF
!
!-- Determine the total number of new particles scheduled for output.
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nr_out_prts_new_pe, numprocs, MPI_INTEGER, MPI_SUM, comm2d,  &
                        ierr )
    nr_out_prts_new = SUM( nr_out_prts_new_pe )
#else
    nr_out_prts_new = nr_out_prts_new_pe(0)
#endif

!
!-- Abort, if selected particles from new particle set would exceed the size of the particle axis
!-- dimension.
    IF ( ( nr_out_prts_new + nr_out_prts ) > max_nr_out_prts )  THEN
       WRITE( message_string, '(A,I10,A,I10,A,I10 )' )                                             &
              'newly released particles exceed maximum number of particles allocated for time ' // &
              'series&allocated particles: ', max_nr_out_prts, '&current output:      ',           &
              nr_out_prts, '&newly added output:  ', nr_out_prts_new
       CALL message( 'dop_mark_ids_scheduled_for_output', 'LPM0031', 3, 2, 0, 6, 0 )
    ENDIF

    start_new_numbering = nr_out_prts + 1
    nr_out_prts = nr_out_prts + nr_out_prts_new

!
!-- Set output IDs of new particles depending on if particles are selected via an increment or via
!-- percentage.
    particle_io_id = start_new_numbering

    IF ( pts_increment > 0 )  THEN

       DO  i = nxl, nxr
          n_all = start_prt_count_yz(i)
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, prt_count(k,j,i)
                   IF ( grid_particles(k,j,i)%particles(n)%particle_mask)  THEN
                      IF ( grid_particles(k,j,i)%particles(n)%io_id == -1 )  THEN
                         IF ( MOD( n_all, pts_increment ) == 0 )  THEN
                            grid_particles(k,j,i)%particles(n)%io_id =                    &
                                                   ( n_all / pts_increment ) + start_new_numbering
                         ELSE
                            grid_particles(k,j,i)%particles(n)%io_id = -2
                         ENDIF
                         n_all = n_all + 1
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( pts_percentage > 0.0_wp )  THEN

       finc   = pts_percentage / 100.0_wp

       DO  i = nxl, nxr
          n_all  = start_prt_count_yz(i)
          finc   = pts_percentage / 100.0_wp
          fcount = REAL( n_all, dp ) * finc + eps
          bcount = INT( REAL( n_all-1, dp ) * finc ) + 1
          IF ( n_all == 0 )  bcount = 0
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                DO  n = 1, prt_count(k,j,i)
                   IF ( grid_particles(k,j,i)%particles(n)%particle_mask )  THEN
                      IF ( grid_particles(k,j,i)%particles(n)%io_id == -1 )  THEN
                         IF ( fcount >= REAL( bcount, dp ) )  THEN
                            particle_io_id = INT( fcount + 0.00001_wp) + start_new_numbering
                            grid_particles(k,j,i)%particles(n)%io_id = particle_io_id
                            bcount = bcount + 1
                         ELSE
                            grid_particles(k,j,i)%particles(n)%io_id = -2
                         ENDIF
                         fcount = fcount + finc
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDIF

 END SUBROUTINE dop_mark_ids_scheduled_for_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Count particles to be transfered to the root model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_count_particles_to_be_transfered_to_root

    IMPLICIT NONE

#if defined( __parallel )
    INTEGER(iwp) :: i               !<
    INTEGER(iwp) :: ierr            !< MPI error code
    INTEGER(iwp) :: j               !<
    INTEGER(iwp) :: k               !<
    INTEGER(iwp) :: n               !<
    INTEGER(iwp) :: iop             !<
    INTEGER(iwp) :: particle_io_id  !<
    INTEGER(iwp) :: pe_nr           !<
    INTEGER(iwp) :: win_size        !<

    INTEGER(iwp), DIMENSION(0:numprocs-1) :: part_ind  !<


!
!-- Count remote particles.
    nr_out_prts_remote = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             DO  n = 1, prt_count(k,j,i)
                particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                IF ( particle_io_id > 0 )  THEN
!>                 TODO: kk this loop has to be optimized
                   DO  iop = 0, numprocs-1
!
!--                   Although the counting is local PE based, the following IF is I/O-PE based,
!--                   because particles in MPI shared memory do not have to be transfered.
                      IF ( particle_io_id < io_start_index  .OR.  particle_io_id > io_end_index )  &
                      THEN
                         IF ( particle_io_id >= mo_indices(1,iop)  .AND.                           &
                              particle_io_id <= mo_indices(2,iop) )                                &
                         THEN
                            nr_out_prts_remote(2,iop) = nr_out_prts_remote(2,iop) + 1
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    nr_out_prts_remote(1,0) = 0
    DO  i = 1, numprocs-1
       nr_out_prts_remote(1,i) = nr_out_prts_remote(1,i-1) + nr_out_prts_remote(2,i-1)
    ENDDO

    win_size = SUM( nr_out_prts_remote(2,:) )
    CALL dop_alloc_rma_mem( transfer_buffer_i, win_size, win_rma_buf_i )
    CALL dop_alloc_rma_mem( transfer_buffer_r, win_size, win_rma_buf_r )

    CALL MPI_ALLTOALL( nr_out_prts_remote, 2, MPI_INTEGER, rma_particles, 2, MPI_INTEGER, comm2d, &
                       ierr)
!
!-- The particle indices are the same for all output variables during one time step.
!-- Therefore, the indices are transfered here only once.
    part_ind = nr_out_prts_remote(1,:)
    transfer_buffer_i = -9999

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

             DO  n = 1, prt_count(k,j,i)
                particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                IF ( particle_io_id < io_start_index  .OR.  particle_io_id > io_end_index )  THEN
                   IF ( particle_io_id > 0 )  THEN
                      pe_nr = get_pe_of_output_particle( particle_io_id )
                      transfer_buffer_i(part_ind(pe_nr)) = particle_io_id
                      part_ind(pe_nr) = part_ind(pe_nr) + 1
                   ENDIF
                ENDIF
             ENDDO

          ENDDO
       ENDDO
    ENDDO

    CALL MPI_BARRIER( comm2d, ierr )

    CALL dop_get_remote_indices

#endif

 END SUBROUTINE dop_count_particles_to_be_transfered_to_root


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In case of nested runs, transfer output particles from all childs to the root model.
!> For the transfer of output particles, there is no parent-child relationship. For this reason,
!> there would be no typical pmcp_c- or pmcp_p-subroutines in pmc_particle_interface.f90.
!> The transfer is done from all childs directly to the root model, instead.
!> All MPI data transfer for output particles is done inside this module, so it does not
!> disturb the general concept, if the transfer of output particles between the models is done also
!> here.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_transfer_from_child_to_root

    IMPLICIT NONE

#if defined( __parallel )
    INTEGER(KIND=MPI_OFFSET_KIND)  ::  disp           !< displacement of actual indices
    INTEGER(iwp)                   ::  ierr           !< MPI error code
    INTEGER(iwp)                   ::  iop            !<
    INTEGER(iwp)                   ::  n              !<
    INTEGER(iwp)                   ::  particle_size  !< particle size in byte
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize        !< size of RMA window

    INTEGER(iwp),DIMENSION(0:cnpo%rem_size-1) ::  local_index  !<

    INTEGER(iwp),DIMENSION(2,cnpo%rem_size) ::  rem_part  !<

    TYPE(particle_type),ALLOCATABLE,DIMENSION(:) ::  local_particle_buffer  !<


    particle_size = STORAGE_SIZE( zero_particle ) / 8

    IF ( root_model )  THEN

       cnpo%index_and_size_buffer = 0

       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )
       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )

       cnpo%nr_rem_out_particle = SUM( cnpo%index_and_size_buffer )

       rem_part(2,:) = cnpo%index_and_size_buffer
       rem_part(1,1) = 0
       DO  n = 2, cnpo%rem_size
          rem_part(1,n) = rem_part(1,n-1) + rem_part(2,n-1)
       ENDDO
       cnpo%index_and_size_buffer = rem_part(1,:)

       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )
       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )

       IF ( ALLOCATED(cnpo%transfer_buffer ) )  DEALLOCATE( cnpo%transfer_buffer )
       ALLOCATE( cnpo%transfer_buffer(cnpo%nr_rem_out_particle) )

       winsize = cnpo%nr_rem_out_particle * particle_size
       CALL MPI_WIN_CREATE( cnpo%transfer_buffer, winsize, particle_size, MPI_INFO_NULL,           &
                            comm_world_nesting, cnpo%buf_win, ierr )

       CALL MPI_WIN_FENCE( 0, cnpo%buf_win, ierr )
       CALL MPI_WIN_FENCE( 0, cnpo%buf_win, ierr )

       CALL MPI_WIN_FREE( cnpo%buf_win, ierr )

    ELSE

       CALL dop_count_child_particles

       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )
       DO  iop = 0, cnpo%rem_size-1
          disp = cnpo%my_rank
          CALL MPI_PUT( nr_out_prts_remote(2,iop), 1, MPI_INTEGER, iop, disp, 1, MPI_INTEGER,      &
                        cnpo%rem_out_win, ierr )
       ENDDO

       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )
       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )

       DO  iop = 0, cnpo%rem_size-1
          disp = cnpo%my_rank
          CALL MPI_GET( nr_out_prts_remote(1,iop), 1, MPI_INTEGER, iop, disp, 1, MPI_INTEGER,      &
                        cnpo%rem_out_win, ierr )
       ENDDO
       CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )

       ALLOCATE( local_particle_buffer(cnpo%nr_rem_out_particle) )
       CALL copy_child_particles( local_particle_buffer )

       ALLOCATE( cnpo%transfer_buffer(1) )
       winsize = 1 * particle_size
       CALL MPI_WIN_CREATE( cnpo%transfer_buffer, winsize, particle_size, MPI_INFO_NULL,           &
                            comm_world_nesting, cnpo%buf_win, ierr )
       CALL MPI_WIN_FENCE( 0, cnpo%buf_win, ierr )

       local_index(0) = 1
       DO  iop = 1, cnpo%rem_size-1
          local_index(iop) = local_index(iop-1) + nr_out_prts_remote(2,iop-1)
       ENDDO
       DO  iop = 0, cnpo%rem_size-1
          n    = nr_out_prts_remote(2,iop) * particle_size
          disp = nr_out_prts_remote(1,iop)
          IF ( n > 0 )  THEN
             CALL MPI_PUT( local_particle_buffer(local_index(iop)), n, MPI_BYTE, iop, disp, n,     &
                           MPI_BYTE, cnpo%buf_win, ierr )
          ENDIF
       ENDDO
       CALL MPI_WIN_FENCE( 0, cnpo%buf_win, ierr )

       CALL MPI_WIN_FREE( cnpo%buf_win, ierr )
       DEALLOCATE( cnpo%transfer_buffer )
       DEALLOCATE( local_particle_buffer )

    ENDIF
#endif

 END SUBROUTINE dop_transfer_from_child_to_root


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In case of nested runs, count particles on child processes to prepare sending to root model.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE dop_count_child_particles

    IMPLICIT NONE

    INTEGER(iwp) :: i               !<
    INTEGER(iwp) :: j               !<
    INTEGER(iwp) :: k               !<
    INTEGER(iwp) :: n               !<
    INTEGER(iwp) :: iop             !<
    INTEGER(iwp) :: particle_io_id  !<


!
!-- Count remote particles.
    nr_out_prts_remote = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             DO  n = 1, prt_count(k,j,i)
                particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                IF ( particle_io_id > 0 )  THEN
!>                 TODO: kk this loop has to be optimized
                   DO  iop = 0, cnpo%rem_size-1
                      IF ( particle_io_id >= mo_indices(1,iop)  .AND.                              &
                           particle_io_id <= mo_indices(2,iop) )                                   &
                      THEN
                         nr_out_prts_remote(2,iop) = nr_out_prts_remote(2,iop) + 1
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    cnpo%nr_rem_out_particle = SUM( nr_out_prts_remote(2,:) )

END SUBROUTINE dop_count_child_particles
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In case of nested runs, copy particles on child processes to prepare sending to root model.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE copy_child_particles( local_particle_buffer )

    IMPLICIT NONE

    TYPE(particle_type), DIMENSION(:) ::  local_particle_buffer  !<

    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  iop             !<
    INTEGER(iwp) ::  j               !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  n               !<
    INTEGER(iwp) ::  particle_io_id  !<

    INTEGER(iwp), DIMENSION(0:cnpo%rem_size-1) ::  index_count  !<

    REAL(wp) ::  lower_left_x  !<
    REAL(wp) ::  lower_left_y  !<


!
!-- Count remote particles.
    index_count(0) = 1
    DO  iop = 1, cnpo%rem_size-1
       index_count(iop) = index_count(iop-1) + nr_out_prts_remote(2,iop-1)
    ENDDO

    CALL pmc_get_model_info( lower_left_x = lower_left_x, lower_left_y = lower_left_y )
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             DO  n = 1, prt_count(k,j,i)
                particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                IF ( particle_io_id > 0 )  THEN
!>                 TODO: kk this loop has to be optimized
                   DO  iop = 0, cnpo%rem_size-1
                      IF ( particle_io_id >= mo_indices(1,iop)  .AND.                              &
                           particle_io_id <= mo_indices(2,iop) )                                   &
                      THEN
                         local_particle_buffer(index_count(iop)) =                                 &
                                           grid_particles(k,j,i)%particles(n)
                         local_particle_buffer(index_count(iop))%x =                               &
                                           local_particle_buffer(index_count(iop))%x + lower_left_x
                         local_particle_buffer(index_count(iop))%y =                               &
                                           local_particle_buffer(index_count(iop))%y + lower_left_y
                         index_count(iop) = index_count(iop)+1
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

END SUBROUTINE copy_child_particles
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fill output buffers.
!> Local variable values are copied into output buffer. Here, local means that they belong to the
!> same shared memory group.
!> The child variable values are copied into the transfer buffer.
!> This routine is called by all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_fill_buffers( var )

    IMPLICIT NONE

    CHARACTER(LEN=32) :: local_name

    INTEGER(iwp) :: i               !<
    INTEGER(iwp) :: j               !<
    INTEGER(iwp) :: k               !<
    INTEGER(iwp) :: local_len       !<
    INTEGER(iwp) :: n               !<
    INTEGER(iwp) :: particle_io_id  !<
    INTEGER(iwp) :: pe_nr           !<
    INTEGER(idp) :: pval            !<

    INTEGER(iwp), DIMENSION(0:numprocs-1) :: part_ind  !<

    TYPE(particle_feature), INTENT(IN) :: var  !<


    part_ind = nr_out_prts_remote(1,:)
    transfer_buffer_i = -9998

!
!-- Filling output buffer is the same for variable name and variable name_const, therefore set
!-- local_name without _const.
    local_len = INDEX( TRIM( var%name ), '_const' )
    IF ( local_len == 0 )  THEN
       local_name = var%name
    ELSE
       local_name = var%name(1:local_len-1)
    ENDIF

!
!-- All particles which are located in the share memory area of the respective IO thread are copied
!-- into the output buffer. The other output particle are copied into the transfer buffer.
    SELECT CASE ( TRIM( local_name ) )

       CASE ( 'origin_x' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n =1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%origin_x
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                        grid_particles(k,j,i)%particles(n)%origin_x
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'origin_y' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%origin_y
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                        grid_particles(k,j,i)%particles(n)%origin_y
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'origin_z' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%origin_z
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                        grid_particles(k,j,i)%particles(n)%origin_z
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'id_low' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      pval = IBITS( grid_particles(k,j,i)%particles(n)%id,0,32 )
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_i(particle_io_id) = INT( pval, 4 )
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_i(part_ind(pe_nr)) = INT( pval, 4 )
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'id_high' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      pval = IBITS( grid_particles(k,j,i)%particles(n)%id,32,32 )
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_i(particle_io_id) = INT( pval, 4 )
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_i(part_ind(pe_nr)) = INT( pval, 4 )
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'particle_io_id' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_i(particle_io_id) = grid_particles(k,j,i)%particles(n)%io_id
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_i(part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%io_id
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'class' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_i(particle_io_id) = grid_particles(k,j,i)%particles(n)%class
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_i(part_ind(pe_nr)) =                                      &
                                                           grid_particles(k,j,i)%particles(n)%class
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'group' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_i(particle_io_id) = grid_particles(k,j,i)%particles(n)%group
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_i(part_ind(pe_nr)) =                                      &
                                                           grid_particles(k,j,i)%particles(n)%group
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'x' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%x
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%x
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'y' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%y
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%y
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'z' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%z
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%z
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'speed_x' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%speed_x
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                         grid_particles(k,j,i)%particles(n)%speed_x
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'speed_y' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%speed_y
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                         grid_particles(k,j,i)%particles(n)%speed_y
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'speed_z' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%speed_z
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r (part_ind(pe_nr)) =                                     &
                                                         grid_particles(k,j,i)%particles(n)%speed_z
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'radius' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k= nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%radius
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                          grid_particles(k,j,i)%particles(n)%radius
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'age' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%age
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%age
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'age_m' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%age_m
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                           grid_particles(k,j,i)%particles(n)%age_m
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dt_sum' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%dt_sum
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                          grid_particles(k,j,i)%particles(n)%dt_sum
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE( 'e_m' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%e_m
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%e_m
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'weight_factor' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%weight_factor
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                   grid_particles(k,j,i)%particles(n)%weight_factor
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'aux1' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%aux1
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                            grid_particles(k,j,i)%particles(n)%aux1
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'aux2' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%aux2
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                            grid_particles(k,j,i)%particles(n)%aux2
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'rvar1' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%rvar1
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                           grid_particles(k,j,i)%particles(n)%rvar1
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'rvar2' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%rvar2
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r (part_ind(pe_nr)) =                                     &
                                                           grid_particles(k,j,i)%particles(n)%rvar2
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'rvar3' )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   DO  n = 1, prt_count(k,j,i)
                      particle_io_id = grid_particles(k,j,i)%particles(n)%io_id
                      IF ( particle_io_id >= io_start_index .AND. particle_io_id <= io_end_index ) &
                      THEN
                         out_buf_r(particle_io_id) = grid_particles(k,j,i)%particles(n)%rvar3
                      ELSEIF ( particle_io_id > 0 )  THEN
                         pe_nr = get_pe_of_output_particle( particle_io_id )
                         transfer_buffer_r(part_ind(pe_nr)) =                                      &
                                                           grid_particles(k,j,i)%particles(n)%rvar3
                         part_ind(pe_nr) = part_ind(pe_nr) + 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
    END SELECT

 END SUBROUTINE dop_fill_buffers


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fill ouput buffer.
!> Copy child particles into output buffer.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE dop_child_var_to_out_buf( var )

    IMPLICIT NONE

    CHARACTER(LEN=32) ::  local_name  !<

    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  local_len       !<
    INTEGER(iwp) ::  particle_io_id  !<
    INTEGER(idp) ::  pval            !<

    TYPE(particle_feature), INTENT(IN) ::  var  !<


!
!-- Filling output buffer is the same for variable name and variable name_const, therefore set
!-- local_name without _const.
    local_len = INDEX( TRIM( var%name ), '_const' )
    IF ( local_len == 0 )  THEN
       local_name = var%name
    ELSE
       local_name = var%name(1:local_len-1)
    ENDIF

!
!-- Return, if no particles are in the respective buffer.
    IF ( .NOT. ALLOCATED( cnpo%transfer_buffer ) )  THEN
       RETURN
    ENDIF

!
!-- All particles which are located in the share memory area of the respective IO thread are copied
!-- into the output buffer. The other output particle are copied into the transfer buffer.
    SELECT CASE ( TRIM( local_name ) )

       CASE ( 'origin_x' )

          DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%origin_x
              ENDIF
           ENDDO

       CASE ( 'origin_y' )

          DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%origin_y
              ENDIF
           ENDDO

       CASE ( 'origin_z' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%origin_z
              ENDIF
           ENDDO

       CASE ( 'id_low' )

          DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              pval = IBITS( cnpo%transfer_buffer(i)%id, 0, 32 )
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_i(particle_io_id) = pval
              ENDIF
           ENDDO

       CASE ( 'id_high' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              pval = IBITS( cnpo%transfer_buffer(i)%id, 32, 32 )
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_i(particle_io_id) = pval
              ENDIF
           ENDDO

       CASE ( 'particle_io_id' )

          DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%io_id
              ELSE
!
!--              Particle I/O ID is not in the define range. This should never happen, but for
!--              safety reasons it is checked here.
                 WRITE( message_string, '(A,3(1X,I10))' )  'particle I/O-id not in defined range', &
                                                           particle_io_id, pe_start_index,         &
                                                           pe_end_index
                 CALL message( 'dop_child_var_to_out_buf', 'LPM0032', 3, 2, 0, 6, 0 )
              ENDIF
           ENDDO

       CASE ( 'class' )

          DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_i(particle_io_id) = cnpo%transfer_buffer(i)%class
              ENDIF
           ENDDO

       CASE ( 'group' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_i(particle_io_id) = cnpo%transfer_buffer(i)%group
              ENDIF
           ENDDO

       CASE ( 'x' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%x
              ENDIF
           ENDDO

       CASE ( 'y' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%y
              ENDIF
           ENDDO

       CASE ( 'z' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%z
              ENDIF
           ENDDO

       CASE ( 'speed_x' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%speed_x
              ENDIF
           ENDDO

       CASE ( 'speed_y' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%speed_y
              ENDIF
           ENDDO

       CASE ( 'speed_z' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%speed_z
              ENDIF
           ENDDO

       CASE ( 'radius' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%radius
              ENDIF
           ENDDO

       CASE ( 'age' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%age
              ENDIF
           ENDDO

       CASE ( 'age_m' )
           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%age_m
              ENDIF
           ENDDO

       CASE ( 'dt_sum' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%dt_sum
              ENDIF
           ENDDO

       CASE( 'e_m' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%e_m
              ENDIF
           ENDDO

       CASE ( 'weight_factor' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%weight_factor
              ENDIF
           ENDDO

       CASE ( 'aux1' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%aux1
              ENDIF
           ENDDO

       CASE ( 'aux2' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%aux2
              ENDIF
           ENDDO

       CASE ( 'rvar1' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%rvar1
              ENDIF
           ENDDO

       CASE ( 'rvar2' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%rvar2
              ENDIF
           ENDDO

       CASE ( 'rvar3' )

           DO  i = 1, SIZE( cnpo%transfer_buffer )
              particle_io_id = cnpo%transfer_buffer(i)%io_id
              IF ( particle_io_id >= pe_start_index  .AND.  particle_io_id <= pe_end_index )  THEN
                 out_buf_r(particle_io_id) = cnpo%transfer_buffer(i)%rvar3
              ENDIF
           ENDDO

    END SELECT

 END SUBROUTINE dop_child_var_to_out_buf
#endif


#if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Get indices (displacement) of remote particles.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_get_remote_indices

    IMPLICIT NONE

    INTEGER(iwp) ::  bufsize    !< size of remote indices array
    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  ierr       !< MPI error code
    INTEGER(iwp) ::  ind_local  !< index in remore indices array

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  disp  !< displacement in RMA window


    bufsize = SUM( rma_particles(2,:) )
    ALLOCATE( remote_indices(0:bufsize) )
    remote_indices = -1

    ind_local = 0
    CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )
    DO  i = 0, numprocs-1
       IF ( rma_particles(2,i) > 0 )  THEN
          disp = rma_particles(1,i)
          IF ( rma_particles(2,i) > 0 )  THEN
             CALL MPI_GET( remote_indices(ind_local), rma_particles(2,i), MPI_INTEGER, i, disp,    &
                           rma_particles(2,i), MPI_INTEGER, win_rma_buf_i, ierr )
             ind_local = ind_local + rma_particles(2,i)
          ENDIF
       ENDIF
    ENDDO
    CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )

 END SUBROUTINE dop_get_remote_indices
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transfer I/O particles from childs to root.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_get_remote_particle( is_integer )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) ::  is_integer  !<

#if defined( __parallel )
    INTEGER(iwp) ::  bufsize         !< size of remote data array
    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  ierr            !< MPI error code
    INTEGER(iwp) ::  ind_local       !< index in remore indices array
    INTEGER(iwp) ::  j               !<
    INTEGER(iwp) ::  particle_io_id  !< I/O ID of particle

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  disp  !< displacement in RMA window

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) :: rma_buf_i  !< buffer to receive remote data (INTEGER)

    REAL(sp), ALLOCATABLE, DIMENSION(:) ::  rma_buf_r  !< buffer to receive remote data (REAL)


    bufsize   = SUM( rma_particles(2,:) )
    ind_local = 0
    ALLOCATE( rma_buf_r(0:bufsize-1) )
    ALLOCATE( rma_buf_i(0:bufsize-1) )

    IF ( is_integer )  THEN
       CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )
    ELSE
       CALL MPI_WIN_FENCE( 0, win_rma_buf_r, ierr )
    ENDIF
    DO  i = 0, numprocs-1
       IF ( rma_particles(2,i) > 0 )  THEN
          IF ( is_integer )  THEN
             disp = rma_particles(1,i)
             CALL MPI_GET( rma_buf_i(ind_local), rma_particles(2,i), MPI_INTEGER, i, disp,         &
                           rma_particles(2,i), MPI_INTEGER,  win_rma_buf_i, ierr )
             ind_local = ind_local + rma_particles(2,i)
          ELSE
             disp = rma_particles(1,i)
             CALL MPI_GET( rma_buf_r(ind_local), rma_particles(2,i), MPI_real, i, disp,            &
                           rma_particles(2,i), MPI_real,  win_rma_buf_r, ierr )
             ind_local = ind_local + rma_particles(2,i)
          ENDIF
       ENDIF
    ENDDO
    IF ( is_integer )  THEN
       CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )
    ELSE
       CALL MPI_WIN_FENCE( 0, win_rma_buf_r, ierr )
    ENDIF

    ind_local = 0

    DO  i = 0, numprocs-1
       IF ( rma_particles(2,i) > 0 )  THEN
          IF ( is_integer )  THEN
!
!--          Copy data from remote PEs into output array.
             DO  j = 0, rma_particles(2,i)-1
                particle_io_id = remote_indices(ind_local)
                out_buf_i(particle_io_id) = rma_buf_i(ind_local)
                ind_local = ind_local + 1
             ENDDO
          ELSE
!
!--          Copy data from remote PEs into output array.
             DO   j = 0, rma_particles(2,i)-1
                particle_io_id = remote_indices(ind_local)
                out_buf_r(particle_io_id) = rma_buf_r(ind_local)
                ind_local = ind_local + 1
             ENDDO
          ENDIF
       ENDIF
    ENDDO

    IF ( ALLOCATED( rma_buf_r) )  DEALLOCATE( rma_buf_r )
    IF ( ALLOCATED( rma_buf_i) )  DEALLOCATE( rma_buf_i )
#else
   IF ( is_integer )  THEN
   ENDIF
#endif

 END SUBROUTINE dop_get_remote_particle


#if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate memory and create window for one-sided communication (INTEGER 1-D array).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_alloc_rma_mem_i1( array, idim1, win )

    IMPLICIT NONE

    INTEGER(iwp) ::  ierr     !< MPI error code

    INTEGER(iwp), INTENT(IN)  ::  idim1  !<
    INTEGER(iwp), INTENT(OUT) ::  win    !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window

    INTEGER(isp), DIMENSION(:), POINTER, INTENT(INOUT) ::  array  !<


    winsize = MAX( idim1, 2 )

    ALLOCATE( array(0:winsize-1) )

    winsize = winsize * isp

    CALL MPI_WIN_CREATE( array, winsize, isp, MPI_INFO_NULL, comm2d, win, ierr )

    array = -1

    CALL MPI_WIN_FENCE( 0, win, ierr )

 END SUBROUTINE dop_alloc_rma_mem_i1
#endif



#if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate memory and create window for one-sided communication (REAL 1-D array).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dop_alloc_rma_mem_r1( array, idim1, win )

    IMPLICIT NONE

    INTEGER(iwp) ::  ierr  !< MPI error code

    INTEGER(iwp), INTENT(IN)  ::  idim1  !<
    INTEGER(iwp), INTENT(OUT) ::  win    !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window

    REAL(sp), DIMENSION(:), POINTER, INTENT(INOUT) ::  array   !<


    winsize = MAX( idim1, 2 )

    ALLOCATE( array(0:winsize-1) )

    winsize = winsize * sp

    CALL MPI_WIN_CREATE( array, winsize, sp, MPI_INFO_NULL, comm2d, win, ierr )

    array = -1.0_wp

    CALL MPI_WIN_FENCE( 0, win, ierr )

 END SUBROUTINE dop_alloc_rma_mem_r1
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Free all MPI windows. Required because of unclear MPI_FINALIZE problems.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE deallocate_and_free

    IMPLICIT NONE

#if defined( __parallel )
    INTEGER(iwp) ::  ierr  !< MPI error code
#endif

#if defined( __parallel )
    CALL MPI_WIN_FREE( win_rma_buf_i, ierr )
    CALL MPI_WIN_FREE( win_rma_buf_r, ierr )
#endif
    IF ( ALLOCATED( remote_indices ) )  DEALLOCATE( remote_indices )

    DEALLOCATE( transfer_buffer_i )
    DEALLOCATE( transfer_buffer_r )

 END SUBROUTINE deallocate_and_free


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Find the number of the PE from which the respective output particle will be transferred to the
!> root model.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_pe_of_output_particle( particle_io_id )  RESULT( pe_nr )

    IMPLICIT NONE

    INTEGER(iwp) :: base   !<
    INTEGER(iwp) :: pe_nr  !<
    INTEGER(iwp) :: pnr    !<

    INTEGER(iwp), INTENT(IN) ::  particle_io_id  !<


    IF ( irregular_distribution )  THEN
       IF ( particle_io_id <= nr_out_prts_rest * nr_out_prts_on_this_pe )  THEN
          pe_nr = ( particle_io_id - 1 ) / nr_out_prts_on_this_pe
       ELSE
          base  = nr_out_prts_rest * nr_out_prts_on_this_pe
          pnr   = particle_io_id - base
          pe_nr = ( pnr - 1 ) / ( nr_out_prts_on_this_pe - 1 )
          pe_nr = pe_nr + nr_out_prts_rest
       ENDIF
    ELSE
       pe_nr = ( particle_io_id - 1 ) / nr_out_prts_on_this_pe
    ENDIF

 END FUNCTION get_pe_of_output_particle

 END MODULE data_output_particle_mod
