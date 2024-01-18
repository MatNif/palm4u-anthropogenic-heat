!> @file restart_data_mpi_io_mod.f90
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
!
! Description:
! ------------
!> Routines for restart data handling using MPI-IO.
!--------------------------------------------------------------------------------------------------!
 MODULE restart_data_mpi_io_mod

#if defined( __parallel )
    USE MPI
#else
    USE posix_interface,                                                                           &
        ONLY:  posix_close,                                                                        &
               posix_lseek,                                                                        &
               posix_open,                                                                         &
               posix_read,                                                                         &
               posix_write
#endif

    USE, INTRINSIC ::  ISO_C_BINDING

    USE boundary_settings_mod,                                                                     &
        ONLY:  set_lateral_neumann_bc

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               cyclic_fill_initialization,                                                         &
               debug_output,                                                                       &
               debug_string,                                                                       &
               message_string,                                                                     &
               restart_data_format_input,                                                          &
               restart_data_format_output,                                                         &
               restart_file_size

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz,                                                                     &
               exchange_horiz_2d

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nnx,                                                                                &
               nny,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nx_on_file,                                                                         &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               ny_on_file,                                                                         &
               nys,                                                                                &
               nysg,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE kinds

    USE particle_attributes,                                                                       &
        ONLY:  grid_particles,                                                                     &
               particles,                                                                          &
               particle_type,                                                                      &
               prt_count,                                                                          &
               zero_particle

    USE pegrid,                                                                                    &
        ONLY:  comm1dx,                                                                            &
               comm1dy,                                                                            &
               comm2d,                                                                             &
               communicator_configurations,                                                        &
               myid,                                                                               &
               myidx,                                                                              &
               myidy,                                                                              &
               npex,                                                                               &
               npey,                                                                               &
               numprocs

    USE shared_memory_io_mod,                                                                      &
        ONLY:  domain_decomposition_grid_features,                                                 &
               sm_class

    IMPLICIT NONE

    CHARACTER(LEN=128) ::  io_file_name  !> internal variable to communicate filename between different subroutines

#if defined( __parallel )
    INTEGER(iwp)            ::  ierr                              !< error status of MPI-calls
    INTEGER(iwp), PARAMETER ::  rd_offset_kind = MPI_OFFSET_KIND  !< Adress or Offset kind
    INTEGER(iwp), PARAMETER ::  rd_status_size = MPI_STATUS_SIZE  !<
#else
    INTEGER(iwp), PARAMETER ::  rd_offset_kind = C_SIZE_T         !<
    INTEGER(iwp), PARAMETER ::  rd_status_size = 1                !< Not required in sequential mode
#endif

    INTEGER(iwp)            ::  comm_io          !< Communicator for MPI-IO
    INTEGER(iwp)            ::  fh = -1          !< MPI-IO file handle
#if defined( __parallel )
    INTEGER(iwp)            ::  ft_surf = -1     !< MPI filetype surface data
    INTEGER(iwp)            ::  ft_2di_nb        !< MPI filetype 2D array INTEGER no outer boundary
    INTEGER(iwp)            ::  ft_2di8_nb       !< MPI filetype 2D array INTEGER no outer boundary
    INTEGER(iwp)            ::  ft_2d            !< MPI filetype 2D array REAL with outer boundaries
    INTEGER(iwp)            ::  ft_3d            !< MPI filetype 3D array REAL with outer boundaries
    INTEGER(iwp)            ::  ft_3di4 = -1     !< MPI filetype 3D array INTEGER*4
    INTEGER(iwp)            ::  ft_3di8 = -1     !< MPI filetype 3D array INTEGER*8
    INTEGER(iwp)            ::  ft_3dsoil        !< MPI filetype for 3d-soil array
#endif
    INTEGER(iwp)            ::  glo_start        !< global start index on this PE
    INTEGER(isp)            ::  i4_rep           !<  representive variable for INTEGER*4, can be used with HUGE
    INTEGER(idp)            ::  total_number_of_surface_elements  !< total number of surface elements for one variable
#if defined( __parallel )
    INTEGER(iwp)            ::  win_2di          !< MPI window variables
    INTEGER(iwp)            ::  win_2di8         !<
    INTEGER(iwp)            ::  win_2dr          !<
    INTEGER(iwp)            ::  win_3di4 = -1    !<
    INTEGER(iwp)            ::  win_3di8 = -1    !<
    INTEGER(iwp)            ::  win_3dr          !<
    INTEGER(iwp)            ::  win_3ds          !<
    INTEGER(iwp)            ::  win_end   = -1   !<
    INTEGER(iwp)            ::  win_glost = -1   !<
    INTEGER(iwp)            ::  win_out   = -1   !<
    INTEGER(iwp)            ::  win_start = -1   !<
    INTEGER(iwp)            ::  win_surf = -1    !<
#endif

    INTEGER(KIND=rd_offset_kind) ::  array_position   !<
    INTEGER(KIND=rd_offset_kind) ::  header_position  !<
#if defined( __parallel )
    INTEGER(KIND=rd_offset_kind) ::  local_offset     !<
#endif

    INTEGER(iwp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  array_2di   !<
    INTEGER(idp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  array_2di8  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  e_end_index       !< extended end index, every grid cell has at least one value
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  e_start_index     !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  m_end_index       !< module copy of end_index
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  m_start_index     !<
    INTEGER(idp), DIMENSION(:),   ALLOCATABLE ::  nr_surfaces_in_tb !< number of surface elements in transfer buffer
    INTEGER(idp), DIMENSION(:),   ALLOCATABLE ::  s_index_in_tb     !< start index of surface elements in transfer buffer
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  s_index_in_window !< start index for the pe in the rma window
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  transfer_index    !<
!
!-- Indices for cyclic fill
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  c_end_index     !< extended end index cyclic fill
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  c_global_end    !< global end index cyclic-fill
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  c_global_start  !< global start index cyclic-fill
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  c_start_index   !< extended start index cyclic fill
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  o_end_index     !< end index cyclic-fill
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  o_start_index   !< start index cyclic-fill

    INTEGER(isp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  array_3di4  !<
    INTEGER(idp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  array_3di8  !<

    LOGICAL ::  filetypes_created                !<
    LOGICAL ::  io_on_limited_cores_per_node     !< switch to shared memory MPI-IO
    LOGICAL ::  less_than_2g_elements            !< TRUE if less than 2**31-1 elements
    LOGICAL ::  rd_flag                          !< file is opened for read
    LOGICAL ::  total_domain_boundaries_on_file  !< total domain boundaries are read from or written to file
    LOGICAL ::  wr_flag                          !< file is opened for write

#if defined( __parallel )
    INTEGER(idp) ::  array_out_reduction  !< reduction of array_out to allow start with index 1

    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE  ::  local_indices
#endif

    REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS     ::  array_out      !<
#if defined( __parallel )
    REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS     ::  array_1d       !<
#endif
    REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS   ::  array_2d       !<
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  array_3d       !<
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  array_3d_soil  !<

!
!-- Variable to store the grid features (index bounds) of the temporary arrays that are used
!-- to read and write the restart data. They differ depending on if the outer boundary of the
!-- total domain is contained in the restart data or not. iog stands for IO-grid.
    TYPE(domain_decomposition_grid_features) ::  iog  !<

!
!-- General Header (first 32 byte in restart file)
    TYPE general_header
       INTEGER(iwp) :: endian       !< little endian (1) or big endian (2) internal format
       INTEGER(iwp) :: outer_ghost_layers_included  !< if 1, file contains the outer ghost layer data
       INTEGER(iwp) :: nr_arrays    !< number of arrays in restart files
       INTEGER(iwp) :: nr_char      !< number of Text strings entries in header
       INTEGER(iwp) :: nr_int       !< number of INTEGER entries in header
       INTEGER(iwp) :: nr_real      !< number of REAL entries in header
       INTEGER(iwp) :: pes_along_x  !< number of PEs along x-direction during writing restart file
       INTEGER(iwp) :: pes_along_y  !< number of PEs along y-direction during writing restart file
       INTEGER(iwp) :: total_nx     !< total number of points in x-direction
       INTEGER(iwp) :: total_ny     !< total number of points in y-direction
    END TYPE general_header

    TYPE(general_header), TARGET, PUBLIC ::  tgh    !<

    TYPE(sm_class)               ::  sm_io  !<

!
!-- Declaration of varibales for file header section
    INTEGER(KIND=rd_offset_kind)                ::  header_int_index
    INTEGER, PARAMETER                          ::  max_nr_int=256
    CHARACTER(LEN=32), DIMENSION(max_nr_int)    ::  int_names
    INTEGER(KIND=iwp), DIMENSION(max_nr_int)    ::  int_values

    INTEGER(KIND=rd_offset_kind)                ::  header_char_index
    INTEGER, PARAMETER                          ::  max_nr_char=128
    CHARACTER(LEN=128), DIMENSION(max_nr_char)  ::  text_lines

    INTEGER(KIND=rd_offset_kind)                ::  header_real_index
    INTEGER, PARAMETER                          ::  max_nr_real=256
    CHARACTER(LEN=32), DIMENSION(max_nr_real)   ::  real_names
    REAL(KIND=wp), DIMENSION(max_nr_real)       ::  real_values

    INTEGER(KIND=rd_offset_kind)                ::  header_array_index
    INTEGER, PARAMETER                          ::  max_nr_arrays=600
    CHARACTER(LEN=32), DIMENSION(max_nr_arrays) ::  array_names
    INTEGER(KIND=rd_offset_kind), DIMENSION(max_nr_arrays) :: array_offset

!
!-- Variables to handle the cyclic fill initialization mode
    INTEGER ::  comm_cyclic_fill  !< communicator for cyclic fill PEs
#if defined( __parallel )
    INTEGER ::  rmawin_2di        !< RMA window 2d INTEGER
    INTEGER ::  rmawin_2di8       !< RMA window 2d INTEGER*8
    INTEGER ::  rmawin_2d         !< RMA window 2d REAL
    INTEGER ::  rmawin_3d         !< RMA window 3d
#endif

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  remote_pe      !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  remote_pe_s    !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  rma_offset     !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  rma_offset_s   !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  rmabuf_2di     !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  rmabuf_2di8    !<

    LOGICAL ::  cyclic_fill_mode            !< arrays are filled cyclically with data from prerun
    LOGICAL ::  pe_active_for_read = .TRUE. !< this PE is active for reading data from prerun or
                                            !< restart run. For restarts all PEs are active.

    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: rmabuf_2d
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: rmabuf_3d

    TYPE(domain_decomposition_grid_features) ::  mainrun_grid  !< grid variables for the main run
    TYPE(domain_decomposition_grid_features) ::  prerun_grid   !< grid variables for the prerun

!
!-- MPI_INTEGER8 is not standard MPI, but is supported on most MPI distibutions
!-- If not suppported, a workaround could be enabled with the following preprocessor directive
!#if defined( __NO_INTEGER8)
!    INTEGER ::  MPI_INTEGER8  !< MPI data type INTEGER8
!#endif

    SAVE

    PRIVATE

    PUBLIC  restart_file_size, total_number_of_surface_elements

!
!-- PALM interfaces
    INTERFACE rd_mpi_io_check_array
       MODULE PROCEDURE rd_mpi_io_check_array
    END INTERFACE rd_mpi_io_check_array

    INTERFACE rd_mpi_io_close
       MODULE PROCEDURE rd_mpi_io_close
    END INTERFACE rd_mpi_io_close

    INTERFACE rd_mpi_io_check_open
       MODULE PROCEDURE rd_mpi_io_check_open
    END INTERFACE rd_mpi_io_check_open

    INTERFACE rd_mpi_io_open
       MODULE PROCEDURE rd_mpi_io_open
    END INTERFACE rd_mpi_io_open

    INTERFACE rrd_mpi_io
       MODULE PROCEDURE rrd_mpi_io_char
       MODULE PROCEDURE rrd_mpi_io_int
       MODULE PROCEDURE rrd_mpi_io_int_2d
       MODULE PROCEDURE rrd_mpi_io_int8
       MODULE PROCEDURE rrd_mpi_io_int8_2d
       MODULE PROCEDURE rrd_mpi_io_int4_3d
       MODULE PROCEDURE rrd_mpi_io_int4_3d_special
       MODULE PROCEDURE rrd_mpi_io_int8_3d
       MODULE PROCEDURE rrd_mpi_io_logical
       MODULE PROCEDURE rrd_mpi_io_logical_3d
       MODULE PROCEDURE rrd_mpi_io_real
       MODULE PROCEDURE rrd_mpi_io_real_2d
       MODULE PROCEDURE rrd_mpi_io_real_3d
       MODULE PROCEDURE rrd_mpi_io_real_3d_soil
    END INTERFACE rrd_mpi_io

    INTERFACE rrd_mpi_io_global_array
       MODULE PROCEDURE rrd_mpi_io_global_array_int_1d
       MODULE PROCEDURE rrd_mpi_io_global_array_real_1d
       MODULE PROCEDURE rrd_mpi_io_global_array_real_2d
       MODULE PROCEDURE rrd_mpi_io_global_array_real_3d
       MODULE PROCEDURE rrd_mpi_io_global_array_real_4d
    END INTERFACE rrd_mpi_io_global_array

    INTERFACE rrd_mpi_io_surface
       MODULE PROCEDURE rrd_mpi_io_surface
       MODULE PROCEDURE rrd_mpi_io_surface_2d
    END INTERFACE rrd_mpi_io_surface

    INTERFACE rrd_mpi_io_particles
       MODULE PROCEDURE rrd_mpi_io_particles
    END INTERFACE rrd_mpi_io_particles

    INTERFACE rd_mpi_io_particle_filetypes
       MODULE PROCEDURE rd_mpi_io_particle_filetypes
    END INTERFACE rd_mpi_io_particle_filetypes

    INTERFACE rd_mpi_io_surface_filetypes
       MODULE PROCEDURE rd_mpi_io_surface_filetypes
    END INTERFACE rd_mpi_io_surface_filetypes

    INTERFACE wrd_mpi_io
       MODULE PROCEDURE wrd_mpi_io_char
       MODULE PROCEDURE wrd_mpi_io_int
       MODULE PROCEDURE wrd_mpi_io_int_2d
       MODULE PROCEDURE wrd_mpi_io_int8
       MODULE PROCEDURE wrd_mpi_io_int8_2d
       MODULE PROCEDURE wrd_mpi_io_int4_3d
       MODULE PROCEDURE wrd_mpi_io_int4_3d_special
       MODULE PROCEDURE wrd_mpi_io_int8_3d
       MODULE PROCEDURE wrd_mpi_io_logical
       MODULE PROCEDURE wrd_mpi_io_logical_3d
       MODULE PROCEDURE wrd_mpi_io_real
       MODULE PROCEDURE wrd_mpi_io_real_2d
       MODULE PROCEDURE wrd_mpi_io_real_3d
       MODULE PROCEDURE wrd_mpi_io_real_3d_soil
    END INTERFACE wrd_mpi_io

    INTERFACE wrd_mpi_io_global_array
       MODULE PROCEDURE wrd_mpi_io_global_array_int_1d
       MODULE PROCEDURE wrd_mpi_io_global_array_real_1d
       MODULE PROCEDURE wrd_mpi_io_global_array_real_2d
       MODULE PROCEDURE wrd_mpi_io_global_array_real_3d
       MODULE PROCEDURE wrd_mpi_io_global_array_real_4d
    END INTERFACE wrd_mpi_io_global_array

    INTERFACE wrd_mpi_io_particles
       MODULE PROCEDURE wrd_mpi_io_particles
    END INTERFACE wrd_mpi_io_particles

    INTERFACE wrd_mpi_io_surface
       MODULE PROCEDURE wrd_mpi_io_surface
       MODULE PROCEDURE wrd_mpi_io_surface_2d
    END INTERFACE wrd_mpi_io_surface

    PUBLIC  rd_mpi_io_check_array,                                                                 &
            rd_mpi_io_check_open,                                                                  &
            rd_mpi_io_close,                                                                       &
            rd_mpi_io_open,                                                                        &
            rd_mpi_io_particle_filetypes,                                                          &
            rd_mpi_io_surface_filetypes,                                                           &
            rrd_mpi_io,                                                                            &
            rrd_mpi_io_global_array,                                                               &
            rrd_mpi_io_particles,                                                                  &
            rrd_mpi_io_surface,                                                                    &
            wrd_mpi_io,                                                                            &
            wrd_mpi_io_global_array,                                                               &
            wrd_mpi_io_particles,                                                                  &
            wrd_mpi_io_surface

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Open restart file for read or write with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_open( action, file_name, open_for_global_io_only,                            &
                            write_total_domain_boundaries )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  action     !<
    CHARACTER(LEN=*), INTENT(IN) ::  file_name  !<

    INTEGER(iwp)                 ::  i          !<
    INTEGER(iwp)                 ::  gh_size    !<

    INTEGER(KIND=rd_offset_kind) ::  offset     !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    LOGICAL, INTENT(IN), OPTIONAL ::  open_for_global_io_only        !<
    LOGICAL                       ::  set_filetype                   !<
    LOGICAL, INTENT(IN), OPTIONAL ::  write_total_domain_boundaries  !<

#if ! defined( __parallel )
    TYPE(C_PTR)                   ::  buf_ptr  !<
#endif


    offset = 0
    io_on_limited_cores_per_node = .FALSE.

    rd_flag = ( TRIM( action ) == 'READ'  .OR. TRIM( action ) == 'read'  )
    wr_flag = ( TRIM( action ) == 'WRITE' .OR. TRIM( action ) == 'write' )

    IF ( .NOT. ( rd_flag .OR. wr_flag ) )  THEN
       message_string = 'illegal action "' // TRIM( action ) // '" for opening restart files'
       CALL message( 'restart_data_mpi_io_mod', 'PAC0294', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Store name of I/O file to communicate it internally within this module.
    io_file_name = file_name

!
!-- Setup for IO on a limited number of PEs per node (using shared memory MPI)
    IF ( rd_flag )  THEN
       set_filetype = .TRUE.
       IF ( TRIM( restart_data_format_input ) == 'mpi_shared_memory' )  THEN
          io_on_limited_cores_per_node = .TRUE.
       ENDIF
    ENDIF

    IF ( TRIM( restart_data_format_output ) == 'mpi_shared_memory' .AND.  wr_flag )  THEN
       io_on_limited_cores_per_node = .TRUE.
    ENDIF

!
!-- Shared memory MPI is not used for reading of global data
    IF ( PRESENT( open_for_global_io_only )  .AND.  rd_flag )  THEN
       IF ( open_for_global_io_only )  THEN
          io_on_limited_cores_per_node = .FALSE.
          set_filetype                 = .FALSE.
       ENDIF
    ENDIF

!
!-- Determine, if prerun data shall be read and mapped cyclically to the mainrun arrays.
!-- In cyclic fill initialization mode only a subset of the PEs will read.
    cyclic_fill_mode   = .FALSE.
    pe_active_for_read = .TRUE.

    IF ( rd_flag  .AND.  .NOT. PRESENT( open_for_global_io_only )  .AND.                           &
         cyclic_fill_initialization )                                                              &
    THEN
       cyclic_fill_mode = .TRUE.
       CALL setup_cyclic_fill
!
!--    Shared memory IO on limited cores is not allowed for cyclic fill mode.
       CALL sm_io%sm_init_comm( .FALSE. )
    ELSE
       CALL sm_io%sm_init_comm( io_on_limited_cores_per_node )
    ENDIF

!
!-- TODO: add a more detailed meaningful comment about what is happening here
!-- activate model grid
    IF ( cyclic_fill_mode  .AND.  .NOT. pe_active_for_read )  THEN
      CALL mainrun_grid%activate_grid_from_this_class()
      RETURN
    ENDIF

!
!-- Set communicator to be used. If all cores are doing I/O, comm2d is used as usual.
    IF( sm_io%is_sm_active() )  THEN
       comm_io = sm_io%comm_io
    ELSEIF ( cyclic_fill_mode )  THEN
       comm_io = comm_cyclic_fill
    ELSE
       comm_io = comm2d
    ENDIF

!
!-- Decide if the total domain boundaries are written to file.
    IF ( wr_flag )  THEN
       IF ( PRESENT( write_total_domain_boundaries ) )  THEN
          total_domain_boundaries_on_file = write_total_domain_boundaries
       ELSE
!
!--       Boundaries will and need only be written if non-cyclic boundary conditions are used at
!--       least along one direction.
          total_domain_boundaries_on_file = .NOT. ( bc_lr_cyc  .AND.  bc_ns_cyc )
       ENDIF
    ENDIF

!
!-- Create subarrays and file types. In case of read it is not known yet if data include total
!-- domain boundaries. Therefore, filetypes will be created further below.
    IF ( wr_flag )  THEN
       CALL rd_mpi_io_create_filetypes
       filetypes_created = .TRUE.
    ELSE
       filetypes_created = .FALSE.
    ENDIF

!
!-- Open file for MPI-IO
#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN

       IF ( rd_flag )  THEN

          IF ( debug_output )  THEN
             WRITE( debug_string, * )  'open joint restart file "' // TRIM( io_file_name ) //      &
                                       '" for read with MPI-IO'
             CALL debug_message( debug_string, 'start' )
          ENDIF

          CALL MPI_FILE_OPEN( comm_io, TRIM( io_file_name ), MPI_MODE_RDONLY, MPI_INFO_NULL, fh,   &
                              ierr )

          IF ( ierr /= 0 )  THEN
             message_string = 'error opening restart file "' // TRIM( io_file_name ) //            &
                              '" for reading with MPI-IO'
             CALL message( 'rrd_mpi_io_open', 'PAC0295', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( debug_output )  THEN
             WRITE( debug_string, * )  'open joint restart file "' // TRIM( io_file_name ) //      &
                                       '" for read with MPI-IO'
             CALL debug_message( debug_string, 'end' )
          ENDIF

       ELSEIF ( wr_flag )  THEN

          IF ( debug_output )  THEN
             WRITE( debug_string, * )  'open joint restart file "' // TRIM( io_file_name ) //      &
                                       '" for write with MPI-IO'
             CALL debug_message( debug_string, 'start' )
          ENDIF

          CALL MPI_FILE_OPEN( comm_io, TRIM( io_file_name ), MPI_MODE_CREATE+MPI_MODE_WRONLY,      &
                              MPI_INFO_NULL, fh, ierr )

          IF ( ierr /= 0 )  THEN
             message_string = 'error opening restart file "' // TRIM( io_file_name ) //            &
                              '" for writing with MPI-IO'
             CALL message( 'rrd_mpi_io_open', 'PAC0296', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( debug_output )  THEN
             WRITE( debug_string, * )  'open joint restart file "' // TRIM( io_file_name ) //      &
                                       '" for write with MPI-IO'
             CALL debug_message( debug_string, 'end' )
          ENDIF

       ENDIF

    ENDIF
#else
    IF ( rd_flag )  THEN

       IF ( debug_output )  THEN
          WRITE( debug_string, * )  'open restart file "' // TRIM( io_file_name ) //               &
                                    '" for read in serial mode (posix)'
          CALL debug_message( debug_string, 'start' )
       ENDIF

       fh = posix_open( TRIM( io_file_name ), .TRUE. )

       IF ( debug_output )  THEN
          WRITE( debug_string, * )  'open restart file "' // TRIM( io_file_name ) //               &
                                    '" for read in serial mode (posix)'
          CALL debug_message( debug_string, 'end' )
       ENDIF

    ELSEIF ( wr_flag )  THEN

       IF ( debug_output )  THEN
          WRITE( debug_string, * )  'open restart file "' // TRIM( io_file_name ) //               &
                                    '" for write in serial mode (posix)'
          CALL debug_message( debug_string, 'start' )
       ENDIF

       fh = posix_open( TRIM( io_file_name ), .FALSE. )

       IF ( debug_output )  THEN
          WRITE( debug_string, * )  'open restart file "' // TRIM( io_file_name ) //               &
                                    '" for write in serial mode (posix)'
          CALL debug_message( debug_string, 'end' )
       ENDIF

    ENDIF

    IF ( fh < 0 )  THEN
       message_string = 'error opening restart file for posix I/O'
       CALL message( 'restart_data_mpi_io_mod', 'PAC0297', 1, 2, 0, 6, 0 )
    ENDIF
#endif

    array_position  = 65536          !> Start offset for writing 2-D and 3.D arrays at 64 k
    header_position = 0

    header_int_index   = 1
    header_char_index  = 1
    header_real_index  = 1
    header_array_index = 1

    int_names    = ' '
    int_values   = 0
    text_lines   = ' '
    real_names   = ' '
    real_values  = 0.0
    array_names  = ' '
    array_offset = 0

    int_names(1)     = 'nx'
    int_values(1)    = nx
    int_names(2)     = 'ny'
    int_values(2)    = ny
    int_names(3)     = 'nz'
    int_values(3)    = nz
    header_int_index = header_int_index+3

    DO  i = 1, max_nr_arrays
       array_offset(i) = 0
       array_names(i)  = ' '
    ENDDO

    gh_size = STORAGE_SIZE( tgh ) / 8

    IF ( rd_flag )  THEN
       IF ( sm_io%iam_io_pe )  THEN
!
!--       File is open for reading
#if defined( __parallel )
!--       Set the default view
          CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
!
!--       Read the file header size
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, tgh, gh_size, MPI_BYTE, status, ierr )
#else
          CALL posix_lseek( fh, header_position )
          buf_ptr = C_LOC( tgh )
          CALL posix_read( fh, buf_ptr, gh_size )
#endif
       ENDIF
#if defined( __parallel )
       IF ( sm_io%is_sm_active() )  THEN
          CALL MPI_BCAST( tgh, gh_size, MPI_BYTE, 0, sm_io%comm_shared, ierr )
       ENDIF
#endif
       header_position = header_position + gh_size

       total_domain_boundaries_on_file = ( tgh%outer_ghost_layers_included == 1 )

!
!--    File types depend on boundaries of the total domain being included in data. This has been
!--    checked with the previous statement.
       IF ( set_filetype )  THEN
          CALL rd_mpi_io_create_filetypes
          filetypes_created = .TRUE.
       ENDIF

       IF ( sm_io%iam_io_pe )  THEN
#if defined( __parallel )
!
!--       Read INTEGER values
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, int_names, SIZE( int_names ) * 32, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( int_names ) * 32

          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ (fh, int_values, SIZE( int_values ), MPI_INT, status, ierr )
          header_position = header_position + SIZE( int_values ) * iwp
!
!--       Character entries
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, text_lines, SIZE( text_lines ) * 128, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE ( text_lines ) * 128
!
!--       REAL values
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, real_names, SIZE( real_names ) * 32, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( real_names ) * 32

          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, real_values, SIZE( real_values ), MPI_REAL, status, ierr )
          header_position = header_position + SIZE( real_values ) * wp
!
!--       2d- and 3d-array headers
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, array_names, SIZE( array_names ) * 32, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( array_names ) * 32

          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_READ( fh, array_offset, SIZE( array_offset ) * MPI_OFFSET_KIND, MPI_BYTE,  &
                              status,ierr )
          header_position = header_position + SIZE( array_offset ) * rd_offset_kind
#else
          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, int_names )
          header_position = header_position + SIZE( int_names ) * 32

          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, int_values, SIZE( int_values ) )
          header_position = header_position + SIZE( int_values ) * iwp
!
!--       Character entries
          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, text_lines )
          header_position = header_position + SIZE( text_lines ) * 128
!
!--       REAL values
          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, real_names )
          header_position = header_position + SIZE( real_names ) * 32

          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, real_values, SIZE( real_values ) )
          header_position = header_position + SIZE( real_values ) * wp
!
!--       2d- and 3d-array headers
          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, array_names )
          header_position = header_position + SIZE( array_names ) * 32

          CALL posix_lseek( fh, header_position )
          CALL posix_read( fh, array_offset, SIZE( array_offset ) )
          header_position = header_position + SIZE( array_offset ) * rd_offset_kind
#endif
          IF ( debug_output )  CALL rd_mpi_io_print_header

       ENDIF

#if defined( __parallel )
!
!--    Broadcast header to all remaining cores that are not involved in I/O
       IF ( sm_io%is_sm_active() )  THEN
!
!--        Not sure, that it is possible to broadcast CHARACTER array in one MPI_Bcast call
           DO  i = 1, SIZE( int_names )
              CALL MPI_BCAST( int_names(i), 32, MPI_CHARACTER, 0, sm_io%comm_shared, ierr )
           ENDDO
           CALL MPI_BCAST( int_values, SIZE( int_values ), MPI_INTEGER, 0, sm_io%comm_shared, ierr )

           DO  i = 1, SIZE( text_lines )
              CALL MPI_BCAST( text_lines(i), 128, MPI_CHARACTER, 0, sm_io%comm_shared, ierr )
           ENDDO

           DO  i = 1, SIZE( real_names )
              CALL MPI_BCAST( real_names(i), 32, MPI_CHARACTER, 0, sm_io%comm_shared, ierr )
           ENDDO
           CALL MPI_BCAST( real_values, SIZE( real_values ), MPI_REAL, 0, sm_io%comm_shared, ierr )

           DO  i = 1, SIZE( array_names )
              CALL MPI_BCAST( array_names(i), 32, MPI_CHARACTER, 0, sm_io%comm_shared, ierr )
           ENDDO
           CALL MPI_BCAST( array_offset, SIZE( array_offset )*8, MPI_BYTE, 0, sm_io%comm_shared,   &
                           ierr )

           CALL MPI_BCAST( header_position, rd_offset_kind, MPI_BYTE, 0, sm_io%comm_shared, ierr )

       ENDIF
#endif


    ENDIF

!
!-- TODO: describe in more detail what is happening here
!-- activate model grid
    IF ( cyclic_fill_mode )  CALL mainrun_grid%activate_grid_from_this_class()

 CONTAINS

    SUBROUTINE setup_cyclic_fill

       IMPLICIT NONE

       INTEGER      ::  color  !< used to set the IO PEs for MPI_COMM_SPLIT
       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
#if defined( __parallel )
       INTEGER      ::  ierr   !<
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window
#else
       INTEGER(idp) ::  winsize
#endif

!
!--    Save grid information of the mainrun, i.e. grid variables like nxl, nxr, nys, nyn and other
!--    values are stored within the mainrun_grid structure
       CALL mainrun_grid%save_grid_into_this_class()

       ALLOCATE( remote_pe(0:nx_on_file,0:ny_on_file) )
       ALLOCATE( remote_pe_s(0:nx_on_file,0:ny_on_file) )
       ALLOCATE( rma_offset(0:nx_on_file,0:ny_on_file) )
       ALLOCATE( rma_offset_s(0:nx_on_file,0:ny_on_file) )

       remote_pe_s  = 0
       rma_offset_s = 0
!
!--    Determine, if gridpoints of the prerun are located on this PE.
!--    Set the (cyclic) prerun grid.
       nxr = MIN( nxr, nx_on_file )
       IF ( nxl > nx_on_file )  THEN
          nxl = -99
          nxr = -99
          nnx = 0
       ELSE
          nnx =nxr-nxl+1
       ENDIF

       nyn = MIN( nyn, ny_on_file )
       IF ( nys > ny_on_file )  THEN
          nys = -99
          nyn = -99
          nny = 0
       ELSE
          nny = nyn-nys+1
       ENDIF

       nx = nx_on_file
       ny = ny_on_file
!
!--    Determine, if this PE is doing IO
       IF ( nnx > 0  .AND.  nny > 0 )  THEN
          color = 1
          pe_active_for_read = .TRUE.
          remote_pe_s(nxl:nxr,nys:nyn) = myid   ! myid from comm2d
          DO  j = nys, nyn
             DO  i = nxl, nxr
                rma_offset_s(i,j) = ( j-nys ) + ( i-nxl ) * nny
             ENDDO
          ENDDO
       ELSE
#if defined( __parallel )
          color = MPI_UNDEFINED
#endif
          pe_active_for_read = .FALSE.
       ENDIF

#if defined( __parallel )
       CALL MPI_ALLREDUCE( remote_pe_s,  remote_pe,  SIZE(remote_pe_s),  MPI_INTEGER, MPI_SUM,     &
                           comm2d, ierr )
       CALL MPI_ALLREDUCE( rma_offset_s, rma_offset, SIZE(rma_offset_s), MPI_INTEGER, MPI_SUM,     &
                           comm2d, ierr )
       CALL MPI_COMM_SPLIT( comm2d, color, 0, comm_cyclic_fill, ierr )

       IF ( pe_active_for_read )  THEN
          CALL MPI_COMM_SIZE( comm_cyclic_fill, numprocs, ierr )
          CALL MPI_COMM_RANK( comm_cyclic_fill, myid, ierr )
       ENDIF
#else
       remote_pe  = remote_pe_s
       rma_offset = rma_offset_s
       myid       = 0
       numprocs   = 1
#endif
!
!--    Allocate 2d buffers as RMA window, accessible on all PEs
       IF ( pe_active_for_read )  THEN
          ALLOCATE( rmabuf_2di(nys:nyn,nxl:nxr) )
          ALLOCATE( rmabuf_2di8(nys:nyn,nxl:nxr) )
       ELSE
          ALLOCATE( rmabuf_2di(1,1) )
          ALLOCATE( rmabuf_2di8(1,1) )
       ENDIF
       winsize = SIZE( rmabuf_2di ) * iwp

#if defined( __parallel )
       CALL MPI_WIN_CREATE( rmabuf_2di, winsize, iwp, MPI_INFO_NULL, comm2d, rmawin_2di, ierr )
       CALL MPI_WIN_FENCE( 0, rmawin_2di, ierr )
#endif
       winsize = SIZE( rmabuf_2di8 ) * idp

#if defined( __parallel )
       CALL MPI_WIN_CREATE( rmabuf_2di8, winsize, iwp, MPI_INFO_NULL, comm2d, rmawin_2di8, ierr )
       CALL MPI_WIN_FENCE( 0, rmawin_2di8, ierr )
#endif

       IF ( pe_active_for_read )  THEN
          ALLOCATE( rmabuf_2d(nys:nyn,nxl:nxr) )
       ELSE
          ALLOCATE( rmabuf_2d(1,1) )
       ENDIF
       winsize = SIZE( rmabuf_2d ) * wp

#if defined( __parallel )
       CALL MPI_WIN_CREATE( rmabuf_2d, winsize, wp, MPI_INFO_NULL, comm2d, rmawin_2d, ierr )
       CALL MPI_WIN_FENCE( 0, rmawin_2d, ierr )
#endif

!
!--    Allocate 3d buffer as RMA window, accessable on all PEs
       IF ( pe_active_for_read )  THEN
          ALLOCATE( rmabuf_3d(nzb:nzt+1,nys:nyn,nxl:nxr) )
       ELSE
          ALLOCATE( rmabuf_3d(1,1,1) )
       ENDIF
       winsize = SIZE( rmabuf_3d ) * wp

#if defined( __parallel )
       CALL MPI_WIN_CREATE( rmabuf_3d, winsize, wp, MPI_INFO_NULL, comm2d, rmawin_3d, ierr )
       CALL MPI_WIN_FENCE( 0, rmawin_3d, ierr )
#endif

!
!--    Save grid of the prerun, i.e. grid variables like nxl, nxr, nys, nyn and other values
!--    are stored within the prerun_grid structure.
!--    The prerun grid can later be activated by calling prerun_grid%activate_grid_from_this_class()
       CALL prerun_grid%save_grid_into_this_class()
       prerun_grid%comm2d = comm_cyclic_fill

       DEALLOCATE( remote_pe_s, rma_offset_s )

    END SUBROUTINE setup_cyclic_fill

 END SUBROUTINE rd_mpi_io_open



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check, if array exists in restart file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_check_array( name, found )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name  !<

    INTEGER(iwp) ::  i  !<

    LOGICAl ::  found  !<


    DO  i = 1, tgh%nr_arrays
       IF ( TRIM( array_names(i) ) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          RETURN
       ENDIF
    ENDDO

    found = .FALSE.

 END SUBROUTINE rd_mpi_io_check_array



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read INTEGER with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: name  !<

    INTEGER(iwp)                   ::  i      !<
    INTEGER(KIND=iwp), INTENT(OUT) ::  value  !<

    LOGICAL ::  found  !<


    found = .FALSE.
    value = 0

    DO  i = 1, tgh%nr_int
       IF ( TRIM(int_names(i)) == TRIM( name ) )  THEN
          value = int_values(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( .NOT. found )  THEN
       message_string = 'INTEGER variable "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int', 'PAC0298', 3, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE rrd_mpi_io_int



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read INTEGER*8 with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int8( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   ::  name   !<

    INTEGER(KIND=idp), INTENT(OUT) ::  value  !<

    REAL(KIND=wp)                  ::  tmp  !<

!
!-- Preliminary Version, read INTEGER*8 as REAL.
!-- Later create an INTEGER*8 header section.
    CALL rrd_mpi_io ( name, tmp )

    value = NINT( tmp, idp )

 END SUBROUTINE rrd_mpi_io_int8



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read REAL with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_real( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name   !<

    INTEGER(iwp)                 ::  i      !<

    LOGICAL                      ::  found  !<

    REAL(KIND=wp), INTENT(OUT)   ::  value  !<


    found = .FALSE.
    value = 0.0

    DO  i = 1, tgh%nr_real
       IF ( TRIM(real_names(i)) == TRIM( name ) )  THEN
          value = real_values(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( .NOT. found )  THEN
       message_string = 'REAL variable "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int', 'PAC0298', 3, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE rrd_mpi_io_real



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-real array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_real_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name       !< name of array to be read

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status     !< returned status of MPI call
#endif
    INTEGER(iwp)                       ::  i          !< loop index
    INTEGER(iwp)                       ::  shift_i    !< shift of i index
    INTEGER(iwp)                       ::  shift_j    !< shift of j index

    LOGICAL                            ::  found      !< true if array found in restart file

    REAL(wp), INTENT(INOUT), DIMENSION(:,:) ::  data  !< array to be read


    found = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

       IF ( SIZE( data, 1 ) == ( nyn - nys + 1 + 2 * nbgp )  .AND.                                 &
            SIZE( data, 2 ) == ( nxr - nxl + 1 + 2 * nbgp )  )  THEN
!
!--       REAL arrays with ghost boundaries.
!--       This kind of array would be dimensioned in the caller subroutine like this:
!--       REAL, DIMENSION(nys:nyn,nxl:nxr) ::  data
!
!--       Because the dimension bounds of 'data' are not explicitly given, they start at 1. To match
!--       the indices of data to the correct indices of array_2d, they must be shifted.
!--       Index nysg in array_2d corresponds to index 1 in data. Hence, to get the correct position
!--       inside the data array, the index "nysg" must be shifted by "-nysg+1". The same for x.
          shift_i = - nxlg + 1
          shift_j = - nysg + 1

          IF ( cyclic_fill_mode )  THEN

             CALL rrd_mpi_io_real_2d_cyclic_fill

          ELSE

#if defined( __parallel )
             CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited # of cores is inactive
             IF ( sm_io%iam_io_pe )  THEN
                CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_2d, 'native',             &
                                        MPI_INFO_NULL, ierr )
                CALL MPI_FILE_READ_ALL( fh, array_2d, SIZE( array_2d ), MPI_REAL, status, ierr )
             ENDIF
             CALL sm_io%sm_node_barrier()
#else
             CALL posix_lseek( fh, array_position )
             CALL posix_read( fh, array_2d, SIZE( array_2d ) )
#endif

             IF ( total_domain_boundaries_on_file )  THEN
                DO  i = iog%nxl, iog%nxr
                   data(iog%nys-nbgp+shift_j:iog%nyn-nbgp+shift_j,i-nbgp+shift_i) =                &
                     array_2d(i,iog%nys:iog%nyn)
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   data(nys+shift_j:nyn+shift_j,i+shift_i) = array_2d(i,nys:nyn)
                ENDDO
             ENDIF

          ENDIF

          CALL exchange_horiz_2d( data )
          CALL set_lateral_neumann_bc( data )

       ELSE
          message_string = '2d-REAL array "' // TRIM( name ) // '" to be read ' //                 &
                           'from restart file is defined with illegal dimensions'
          CALL message( 'rrd_mpi_io_real_2d', 'PAC0299', 3, 2, 0, 6, 0 )
       ENDIF

    ELSE
       message_string = '2d-REAL array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int', 'PAC0298', 3, 2, 0, 6, 0 )
    ENDIF


 CONTAINS

    SUBROUTINE rrd_mpi_io_real_2d_cyclic_fill

       IMPLICIT NONE

       INTEGER(iwp)    :: i         !<
       INTEGER(iwp)    :: ie        !<
       INTEGER(iwp)    :: is        !<
       INTEGER(iwp)    :: i_remote  !<
       INTEGER(iwp)    :: j         !<
       INTEGER(iwp)    :: je        !<
       INTEGER(iwp)    :: js        !<
       INTEGER(iwp)    :: j_remote  !<
       INTEGER(iwp)    :: nval      !<
       INTEGER(iwp)    :: rem_pe    !<

#if defined( __parallel )
       INTEGER(iwp)    :: ierr      !<
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  rem_offs  !<
#else
       INTEGER(idp) ::  rem_offs
#endif


!
!--    Reading 2d real array on prerun grid
       CALL prerun_grid%activate_grid_from_this_class()

       IF ( pe_active_for_read )  THEN
#if defined( __parallel )
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_2d, 'native', MPI_INFO_NULL,    &
                                  ierr )
          CALL MPI_FILE_READ_ALL( fh, array_2d, SIZE( array_2d ), MPI_REAL, status, ierr )
#endif
          DO  i = nxl, nxr
             rmabuf_2d(nys:nyn,i) = array_2d(i,nys:nyn)
          ENDDO
!
!--       Copy prerund data directly into output array data
          data(nys+shift_j:nyn+shift_j,nxl+shift_i:nxr+shift_i) = rmabuf_2d
       ENDIF

       CALL mainrun_grid%activate_grid_from_this_class()

#if defined( __parallel )
!
!--    Close RMA window to allow remote access
       CALL MPI_WIN_FENCE( 0, rmawin_2d, ierr )
#endif

!
!--    TODO: describe in more detail what is happening in this IF/ELSE clause
       IF ( .NOT. pe_active_for_read )  THEN

          is = nxl
          ie = nxr
          js = nys
          je = nyn

       ELSE
!
!--       Extra get for cyclic data north of prerun data
          is = nxl
          ie = nxr
          js = prerun_grid%nys+1
          je = nyn
          DO  i = is, ie
             DO  j = js, je
                i_remote = MOD(i,nx_on_file+1)
                j_remote = MOD(j,ny_on_file+1)
                rem_pe   = remote_pe(i_remote,j_remote)
                rem_offs = rma_offset(i_remote,j_remote)
                nval     = 1

#if defined( __parallel )
                IF ( rem_pe /= myid )  THEN
                   CALL MPI_GET( data(j+shift_j,i+shift_i), nval, MPI_REAL, rem_pe, rem_offs,      &
                                 nval, MPI_REAL, rmawin_2d, ierr )
                ELSE
                   data(j+shift_j,i+shift_i) = rmabuf_2d(j_remote,i_remote)
                ENDIF
#else
                data(j+shift_j,i+shift_i) = array_2d(i_remote,j_remote)
#endif
             ENDDO
          ENDDO
!
!--       Prepare setup for stripe right of prerun data
          is = prerun_grid%nxr+1
          ie = nxr
          js = nys
          je = nyn

       ENDIF

       DO  i = is, ie
          DO j = js, je
             i_remote = MOD(i,nx_on_file+1)
             j_remote = MOD(j,ny_on_file+1)
             rem_pe   = remote_pe(i_remote,j_remote)
             rem_offs = rma_offset(i_remote,j_remote)
             nval     = 1

#if defined( __parallel )
             IF ( rem_pe /= myid )  THEN
                CALL MPI_GET( data(j+shift_j,i+shift_i), nval, MPI_REAL, rem_pe, rem_offs, nval,   &
                              MPI_REAL, rmawin_2d, ierr )
             ELSE
                data(j+shift_j,i+shift_i) = rmabuf_2d(j_remote,i_remote)
             ENDIF
#else
             data(j+shift_j,i+shift_i) = array_2d(i_remote,j_remote)
#endif
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Reopen RMA window to allow filling
       CALL MPI_WIN_FENCE( 0, rmawin_2d, ierr )
#endif

    END SUBROUTINE rrd_mpi_io_real_2d_cyclic_fill

 END SUBROUTINE rrd_mpi_io_real_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-INTEGER array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<
    INTEGER(iwp)                       ::  j       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    INTEGER(KIND=iwp), INTENT(INOUT), DIMENSION(:,:) ::  data  !<

    LOGICAL ::  found  !<


    found = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

       IF ( ( nxr - nxl + 1 + 2 * nbgp ) == SIZE( data, 2 ) )  THEN
!
!--       Output array with Halos.
!--       ATTENTION: INTEGER arrays with ghost boundaries are not implemented yet. This kind of
!--                  array would be dimensioned in the caller subroutine like this:
!--                  INTEGER, DIMENSION(nysg:nyng,nxlg:nxrg)::  data
          message_string = '2d-INTEGER array "' // TRIM( name ) // '" to be read ' //              &
                           'from restart file is defined with illegal dimensions'
          CALL message( 'rrd_mpi_io_int_2d', 'PAC0301', 3, 2, 0, 6, 0 )

       ELSEIF ( (nxr - nxl + 1) == SIZE( data, 2 ) )  THEN
!
!--       INTEGER input array without Halos.
!--       This kind of array is dimensioned in the caller subroutine
!--       INTEGER, DIMENSION(nys:nyn,nxl:nxr) ::  data
          IF ( cyclic_fill_mode )  THEN

             CALL rrd_mpi_io_int_2d_cyclic_fill

          ELSE

#if defined( __parallel )
             CALL sm_io%sm_node_barrier() ! Has no effect if I/O on limited # of cores is inactive
             IF ( sm_io%iam_io_pe )  THEN
                CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER, ft_2di_nb, 'native',      &
                                        MPI_INFO_NULL, ierr )
                CALL MPI_FILE_READ_ALL( fh, array_2di, SIZE( array_2di ), MPI_INTEGER, status,     &
                                        ierr )
             ENDIF
             CALL sm_io%sm_node_barrier()
#else
             CALL posix_lseek( fh, array_position )
             CALL posix_read( fh, array_2di, SIZE( array_2di ) )
#endif
             DO  j = nys, nyn
                DO  i = nxl, nxr
                   data(j-nys+1,i-nxl+1) = array_2di(i,j)
                ENDDO
             ENDDO

          ENDIF

       ELSE

          message_string = '2d-INTEGER array "' // TRIM( name ) // '" to be read from restart ' // &
                           'file is defined with illegal dimensions'
          CALL message( 'rrd_mpi_io_int_2d', 'PAC0301', 3, 2, 0, 6, 0 )

       ENDIF

    ELSE

       message_string = '2d-INTEGER array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int_2d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF


 CONTAINS

    SUBROUTINE rrd_mpi_io_int_2d_cyclic_fill

       IMPLICIT NONE

       INTEGER(iwp)    :: i         !<
       INTEGER(iwp)    :: ie        !<
       INTEGER(iwp)    :: is        !<
       INTEGER(iwp)    :: i_remote  !<
       INTEGER(iwp)    :: j         !<
       INTEGER(iwp)    :: je        !<
       INTEGER(iwp)    :: js        !<
       INTEGER(iwp)    :: j_remote  !<
       INTEGER(iwp)    :: nval      !<
       INTEGER(iwp)    :: rem_pe    !<

#if defined( __parallel )
       INTEGER(iwp)    :: ierr      !<
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  rem_offs  !<
#else
       INTEGER(idp) ::  rem_offs
#endif


       IF ( pe_active_for_read )  THEN
          CALL prerun_grid%activate_grid_from_this_class()

#if defined( __parallel )
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER, ft_2di_nb, 'native',            &
                                  MPI_INFO_NULL, ierr )
          CALL MPI_FILE_READ_ALL( fh, array_2di, SIZE( array_2di ), MPI_INTEGER, status, ierr )
#else
          CALL posix_lseek( fh, array_position )
          CALL posix_read( fh, array_2di, SIZE( array_2di ) )
#endif
          DO  i = nxl, nxr
             rmabuf_2di(nys:nyn,i) = array_2di(i,nys:nyn)
          ENDDO
          data(1:nny,1:nnx) = rmabuf_2di

          CALL mainrun_grid%activate_grid_from_this_class()
       ENDIF

#if defined( __parallel )
!
!--    Close RMA window to allow remote access
       CALL MPI_WIN_FENCE( 0, rmawin_2di, ierr )
#endif

       is = nxl
       ie = nxr
       js = nys
       je = nyn

       DO  i = is, ie
          DO  j = js, je
             i_remote = MOD(i,nx_on_file+1)
             j_remote = MOD(j,ny_on_file+1)
             rem_pe   = remote_pe(i_remote,j_remote)
             rem_offs = rma_offset(i_remote,j_remote)
             nval     = 1
#if defined( __parallel )
             IF ( rem_pe /= myid )  THEN
                CALL MPI_GET( data(j-nys+1,i-nxl+1), nval, MPI_INTEGER, rem_pe, rem_offs, nval,    &
                              MPI_INTEGER, rmawin_2di, ierr)
             ELSE
                data(j-nys+1,i-nxl+1) = rmabuf_2di(j_remote,i_remote)
             ENDIF
#else
             data(j-nys+1,i-nxl+1) = array_2di(i_remote,j_remote)
#endif
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Reopen RMA window to allow filling
       CALL MPI_WIN_FENCE( 0, rmawin_2di, ierr )
#endif

    END SUBROUTINE rrd_mpi_io_int_2d_cyclic_fill

 END SUBROUTINE rrd_mpi_io_int_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-INTEGER array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int8_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<
    INTEGER(iwp)                       ::  j       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:) ::  data  !<

    LOGICAL ::  found  !<


    found = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM( array_names(i) ) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

       IF ( ( nxr - nxl + 1 + 2 * nbgp ) == SIZE( data, 2 ) )  THEN
!
!--       Output array with ghost layers.
!--       ATTENTION: INTEGER arrays with ghost boundaries are not implemented yet. This kind of
!--                  array would be dimensioned in the caller subroutine like this:
!--                  INTEGER, DIMENSION(nysg:nyng,nxlg:nxrg)::  data
          message_string = '2d-INTEGER*8 array "' // TRIM( name ) // '" to be read ' //            &
                           'from restart file is defined with illegal dimensions'
          CALL message( 'rrd_mpi_io_int8_2d', 'PAC0301', 3, 2, 0, 6, 0 )

       ELSEIF ( (nxr - nxl + 1) == SIZE( data, 2 ) )  THEN
!
!--       INTEGER input array without ghost layers.
!--       This kind of array is dimensioned in the caller subroutine as
!--       INTEGER, DIMENSION(nys:nyn,nxl:nxr) ::  data
          IF ( cyclic_fill_mode )  THEN

             CALL rrd_mpi_io_int8_2d_cyclic_fill

          ELSE
#if defined( __parallel )
             CALL sm_io%sm_node_barrier() ! Has no effect if I/O on limited # of cores is inactive
             IF ( sm_io%iam_io_pe )  THEN
                CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER8, ft_2di8_nb, 'native',    &
                                        MPI_INFO_NULL, ierr )
                CALL MPI_FILE_READ_ALL( fh, array_2di8, SIZE( array_2di8 ), MPI_INTEGER8, status,  &
                                        ierr )
             ENDIF
             CALL sm_io%sm_node_barrier()
#else
             CALL posix_lseek( fh, array_position )
             CALL posix_read( fh, array_2di8, SIZE( array_2di8 ) )
#endif
             DO  j = nys, nyn
                DO  i = nxl, nxr
                   data(j-nys+1,i-nxl+1) = array_2di8(i,j)
                 ENDDO
             ENDDO

          ENDIF

       ELSE

          message_string = '2d-INTEGER*8 array "' // TRIM( name ) // '" to be read from ' //       &
                           'restart file is defined with illegal dimensions'
          CALL message( 'rrd_mpi_io_int8_2d', 'PAC0301', 3, 2, 0, 6, 0 )

       ENDIF

    ELSE

       message_string = '2d-INTEGER*8 array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int8_2d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF

 CONTAINS

    SUBROUTINE rrd_mpi_io_int8_2d_cyclic_fill

       IMPLICIT NONE

       INTEGER(iwp)    :: i         !<
       INTEGER(iwp)    :: ie        !<
       INTEGER(iwp)    :: is        !<
       INTEGER(iwp)    :: i_remote  !<
       INTEGER(iwp)    :: j         !<
       INTEGER(iwp)    :: je        !<
       INTEGER(iwp)    :: js        !<
       INTEGER(iwp)    :: j_remote  !<
       INTEGER(iwp)    :: nval      !<
       INTEGER(iwp)    :: rem_pe    !<

#if defined( __parallel )
       INTEGER(iwp)    :: ierr      !<
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  rem_offs  !<
#else
       INTEGER(idp) ::  rem_offs
#endif


       IF ( pe_active_for_read )  THEN
          CALL prerun_grid%activate_grid_from_this_class()

#if defined( __parallel )
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER8, ft_2di8_nb, 'native',            &
                                  MPI_INFO_NULL, ierr )
          CALL MPI_FILE_READ_ALL( fh, array_2di8, SIZE( array_2di8 ), MPI_INTEGER8, status, ierr )
#else
          CALL posix_lseek( fh, array_position )
          CALL posix_read( fh, array_2di8, SIZE( array_2di8 ) )
#endif
          DO  i = nxl, nxr
             rmabuf_2di8(nys:nyn,i) = array_2di8(i,nys:nyn)
          ENDDO
          data(1:nny,1:nnx) = rmabuf_2di8

          CALL mainrun_grid%activate_grid_from_this_class()
       ENDIF

#if defined( __parallel )
!
!--    Close RMA window to allow remote access
       CALL MPI_WIN_FENCE( 0, rmawin_2di8, ierr )
#endif

       is = nxl
       ie = nxr
       js = nys
       je = nyn

       DO  i = is, ie
          DO  j = js, je
             i_remote = MOD(i,nx_on_file+1)
             j_remote = MOD(j,ny_on_file+1)
             rem_pe   = remote_pe(i_remote,j_remote)
             rem_offs = rma_offset(i_remote,j_remote)
             nval     = 1
#if defined( __parallel )
             IF ( rem_pe /= myid )  THEN
                CALL MPI_GET( data(j-nys+1,i-nxl+1), nval, MPI_INTEGER8, rem_pe, rem_offs, nval,    &
                              MPI_INTEGER8, rmawin_2di8, ierr)
             ELSE
                data(j-nys+1,i-nxl+1) = rmabuf_2di8(j_remote,i_remote)
             ENDIF
#else
             data(j-nys+1,i-nxl+1) = array_2di8(i_remote,j_remote)
#endif
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Reopen RMA window to allow filling
       CALL MPI_WIN_FENCE( 0, rmawin_2di8, ierr )
#endif

    END SUBROUTINE rrd_mpi_io_int8_2d_cyclic_fill

 END SUBROUTINE rrd_mpi_io_int8_2d




!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 3d-INTEGER*4 array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int4_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    LOGICAL                            ::  found   !<

    INTEGER(isp), INTENT(INOUT), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data  !<


    found = .FALSE.
    data  = -1.0

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

#if defined( __parallel )
       CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited # of cores is inactive
       IF( sm_io%iam_io_pe )  THEN
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER, ft_3di4, 'native',              &
                                  MPI_INFO_NULL, ierr )
          CALL MPI_FILE_READ_ALL( fh, array_3di4, SIZE( array_3di4 ), MPI_INTEGER, status, ierr )
       ENDIF
       CALL sm_io%sm_node_barrier()
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_read(fh, array_3di4, SIZE( array_3di4 ) )
#endif
       IF ( total_domain_boundaries_on_file )  THEN
          DO  i = iog%nxl, iog%nxr
             data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp) = array_3di4(:,i,iog%nys:iog%nyn)
          ENDDO
       ELSE
          DO  i = nxl, nxr
             data(:,nys:nyn,i) = array_3di4(:,i,nys:nyn)
          ENDDO
       ENDIF

    ELSE

       message_string = '3d-INTEGER*4 array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int4_3d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF


 END SUBROUTINE rrd_mpi_io_int4_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 3d-INTEGER*4 array with MPI-IO;
!> Dimension order: y, x, special
!> Dimension x and y do not contain ghost layers. The 3rd dimension is of size 'size_3rd_dimension'
!> and must be smaller than (nzt+1) - nzb + 1.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int4_3d_special( name, data, size_3rd_dimension )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name            !< name of variable to be read

    INTEGER(iwp)             ::  n                   !< loop index
    INTEGER(iwp), INTENT(IN) ::  size_3rd_dimension  !< size of 3rd dimension of array 'data'

    INTEGER(isp), INTENT(OUT), DIMENSION(nys:nyn,nxl:nxr,1:size_3rd_dimension) ::  data  !< array to be read

    INTEGER(isp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  out_buf  !< reordered data array

    IF ( size_3rd_dimension > nzt - nzb + 2 )  THEN
       WRITE(message_string, '(A,I6,A,I6,A)')  '3rd dimension of array "' // TRIM( name ) //       &
                                             '" is larger than number of domain height levels:&',  &
                                             size_3rd_dimension, ' > ', nzt - nzb + 1,             &
                                             ', writing impossible'
       CALL message( 'rrd_mpi_io_int4_3d_special', 'PAC0302', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Read data into array of default size (entire subdomain)
    CALL rrd_mpi_io_int4_3d( name, out_buf )

!
!-- Reorder data from file to output array
    DO  n = 1, size_3rd_dimension
       data(nys:nyn,nxl:nxr,n) = out_buf(nzb+n-1,nys:nyn,nxl:nxr)
    ENDDO

END SUBROUTINE rrd_mpi_io_int4_3d_special



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 3d-INTEGER*8 array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_int8_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    LOGICAL                            ::  found   !<

    INTEGER(idp), INTENT(INOUT), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data  !<


    found = .FALSE.
    data  = -1.0

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

#if defined( __parallel )
          CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited # of cores is inactive
          IF( sm_io%iam_io_pe )  THEN
             CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER8, ft_3di8, 'native', MPI_INFO_NULL, &
                                     ierr )
             CALL MPI_FILE_READ_ALL( fh, array_3di8, SIZE( array_3di8 ), MPI_INTEGER8, status, ierr )
          ENDIF
          CALL sm_io%sm_node_barrier()
#else
          CALL posix_lseek( fh, array_position )
          CALL posix_read(fh, array_3di8, SIZE( array_3di8 ) )
#endif
          IF ( total_domain_boundaries_on_file )  THEN
             DO  i = iog%nxl, iog%nxr
                data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp) = array_3di8(:,i,iog%nys:iog%nyn)
             ENDDO
          ELSE
             DO  i = nxl, nxr
                data(:,nys:nyn,i) = array_3di8(:,i,nys:nyn)
             ENDDO
          ENDIF

    ELSE

       message_string = '3d-INTEGER*8 array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_int8_3d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF


 END SUBROUTINE rrd_mpi_io_int8_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-REAL array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_real_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    LOGICAL                            ::  found   !<

    REAL(wp), INTENT(INOUT), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data  !<


    found = .FALSE.
    data  = -1.0

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

       IF ( cyclic_fill_mode )  THEN

          CALL rrd_mpi_io_real_3d_cyclic_fill
!
!--       Cyclic fill mode requires to use the "cyclic" communicator, in order to initialize
!--       grid points at the outer boundaries (ghost layers) of the total domain. These points
!--       are not contained in the prerun data, because the prerun used cyclic boundary conditions.
          CALL exchange_horiz( data, nbgp, alternative_communicator = 1 )

       ELSE
#if defined( __parallel )
          CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited # of cores is inactive
          IF( sm_io%iam_io_pe )  THEN
             CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_3d, 'native', MPI_INFO_NULL, &
                                     ierr )
             CALL MPI_FILE_READ_ALL( fh, array_3d, SIZE( array_3d ), MPI_REAL, status, ierr )
          ENDIF
          CALL sm_io%sm_node_barrier()
#else
          CALL posix_lseek( fh, array_position )
          CALL posix_read(fh, array_3d, SIZE( array_3d ) )
#endif
          IF ( total_domain_boundaries_on_file )  THEN
             DO  i = iog%nxl, iog%nxr
                data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp) = array_3d(:,i,iog%nys:iog%nyn)
             ENDDO
          ELSE
             DO  i = nxl, nxr
                data(:,nys:nyn,i) = array_3d(:,i,nys:nyn)
             ENDDO
          ENDIF

          CALL exchange_horiz( data, nbgp )

       ENDIF

    ELSE

       message_string = '3d-REAL array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_real_3d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF


 CONTAINS

    SUBROUTINE rrd_mpi_io_real_3d_cyclic_fill

       IMPLICIT NONE

       INTEGER(iwp)    :: i         !<
       INTEGER(iwp)    :: ie        !<
       INTEGER(iwp)    :: is        !<
       INTEGER(iwp)    :: i_remote  !<
       INTEGER(iwp)    :: j         !<
       INTEGER(iwp)    :: je        !<
       INTEGER(iwp)    :: js        !<
       INTEGER(iwp)    :: j_remote  !<
       INTEGER(iwp)    :: nval      !<
       INTEGER(iwp)    :: rem_pe    !<

#if defined( __parallel )
       INTEGER(iwp)    :: ierr      !<
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  rem_offs  !<
#else
       INTEGER(idp) ::  rem_offs
#endif


       CALL prerun_grid%activate_grid_from_this_class()

       IF ( pe_active_for_read )  THEN
#if defined( __parallel )
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_3d, 'native', MPI_INFO_NULL,    &
                                  ierr )
          CALL MPI_FILE_READ_ALL( fh, array_3d, SIZE( array_3d ), MPI_REAL, status, ierr )
#else
          CALL posix_lseek( fh, array_position )
          CALL posix_read( fh, array_3d, SIZE( array_3d ) )
#endif
          DO  i = nxl, nxr
             rmabuf_3d(:,nys:nyn,i) = array_3d(:,i,nys:nyn)
          ENDDO
          data(:,nys:nyn,nxl:nxr) = rmabuf_3d
       ENDIF
       CALL mainrun_grid%activate_grid_from_this_class ()

#if defined( __parallel )
!
!--    Close RMA window to allow remote access
       CALL MPI_WIN_FENCE( 0, rmawin_3d, ierr )
#endif

       is = nxl
       ie = nxr
       js = nys
       je = nyn

       DO  i = is, ie
          DO  j = js, je
             i_remote = MOD( i, nx_on_file+1 )
             j_remote = MOD( j, ny_on_file+1 )
             rem_pe   = remote_pe(i_remote,j_remote)
             rem_offs = rma_offset(i_remote,j_remote) * ( nzt-nzb+2 )
             nval     = nzt-nzb+2

#if defined( __parallel )
             IF ( rem_pe /= myid )  THEN
                CALL MPI_GET( data(nzb,j,i), nval, MPI_REAL, rem_pe, rem_offs, nval, MPI_REAL,     &
                              rmawin_3d, ierr)
             ELSE
                data(:,j,i) = rmabuf_3d(:,j_remote,i_remote)
             ENDIF
#else
             data(:,j,i) = array_3d(:,i_remote,j_remote)
#endif
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Reopen RMA window to allow filling
       CALL MPI_WIN_FENCE( 0, rmawin_3d, ierr )
#endif

    END SUBROUTINE rrd_mpi_io_real_3d_cyclic_fill

 END SUBROUTINE rrd_mpi_io_real_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 3d-REAL soil array with MPI-IO
!> nzb_soil, nzt_soil are located in the module land_surface_model_mod. Since Fortran does not allow
!> cross referencing of module variables, it is required to pass these variables as arguments.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_real_3d_soil( name, data, nzb_soil, nzt_soil )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name      !<

    INTEGER(iwp)                       ::  i         !<
    INTEGER, INTENT(IN)                ::  nzb_soil  !<
    INTEGER, INTENT(IN)                ::  nzt_soil  !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status    !<
    INTEGER(iwp)                       ::  ierr      !<
#endif

    LOGICAL                            ::  found     !<

    REAL(wp), INTENT(INOUT), DIMENSION(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) ::  data  !<


!
!-- Prerun data is not allowed to contain soil information so far
    IF ( cyclic_fill_mode )  THEN
       message_string = 'prerun data is not allowed to contain soil information'
       CALL message( 'rrd_mpi_io_real_3d_soil', 'PAC0303', 3, 2, -1, 6, 0 )
    ENDIF

    found = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN
       CALL rd_mpi_io_create_filetypes_3dsoil( nzb_soil, nzt_soil )
#if defined( __parallel )
       CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
       IF ( sm_io%iam_io_pe )  THEN
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_3dsoil, 'native',               &
                                  MPI_INFO_NULL, ierr )
          CALL MPI_FILE_READ_ALL( fh, array_3d_soil, SIZE( array_3d_soil ), MPI_REAL, status, ierr )
          CALL MPI_TYPE_FREE( ft_3dsoil, ierr )
       ENDIF
       CALL sm_io%sm_node_barrier()
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_read( fh, array_3d_soil, SIZE( array_3d_soil ) )
#endif
       IF ( total_domain_boundaries_on_file )  THEN
          DO  i = iog%nxl, iog%nxr
             data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp) = array_3d_soil(:,i,iog%nys:iog%nyn)
          ENDDO
       ELSE
          DO  i = nxl, nxr
             data(:,nys:nyn,i) = array_3d_soil(:,i,nys:nyn)
          ENDDO
       ENDIF

#if defined( __parallel )
      IF ( sm_io%is_sm_active() )  THEN
         CALL sm_io%sm_free_shared( win_3ds )
      ELSE
         DEALLOCATE( array_3d_soil )
      ENDIF
#else
      DEALLOCATE( array_3d_soil )
#endif

    ELSE

       message_string = '3d-REAL soil array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_real_3d_soil', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE rrd_mpi_io_real_3d_soil


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read CHARACTER with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_char( name, text )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  ::  name   !<
    CHARACTER(LEN=*), INTENT(OUT) ::  text   !<
    CHARACTER(LEN=128)            ::  line   !<

    INTEGER(iwp)                  ::  i      !<

    LOGICAL                       ::  found  !<


    found = .FALSE.
    text = ' '

    DO  i = 1, tgh%nr_char
       line = text_lines(i)
       IF ( TRIM( line(1:32) ) == TRIM( name ) )  THEN
          text = line(33:)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( .NOT. found )  THEN
       message_string = 'CHARACTER variable "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_char', 'PAC0298', 3, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE rrd_mpi_io_char



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read LOGICAL with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_logical( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name                !<

    INTEGER(iwp)                 ::  logical_as_integer  !<

    LOGICAL, INTENT(OUT)         ::  value               !<


    CALL rrd_mpi_io_int( name, logical_as_integer )
    value = ( logical_as_integer == 1 )

 END SUBROUTINE rrd_mpi_io_logical



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 3D-LOGICAL with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_logical_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name                !<

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  logical_as_integer  !<

    LOGICAL, INTENT(OUT), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data  !<


    CALL rrd_mpi_io_int4_3d( name, logical_as_integer )
    data(:,:,:) = ( logical_as_integer(:,:,:) == 1 )

 END SUBROUTINE rrd_mpi_io_logical_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write INTEGER with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  ::  name   !<

    INTEGER(KIND=iwp), INTENT(IN) ::  value  !<


    IF ( header_int_index == max_nr_int )  THEN
       STOP '+++ maximum number of INTEGER entries in restart file header exceeded'
    ENDIF

    int_names(header_int_index)  = name
    int_values(header_int_index) = value
    header_int_index = header_int_index + 1

 END SUBROUTINE wrd_mpi_io_int



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write INTEGER*8 with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int8( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  ::  name   !<

    INTEGER(KIND=idp), INTENT(IN) ::  value  !<

    REAL(KIND=wp)                 ::  tmp  !<

!
!-- Preliminary Version, write INTEGER*8 as REAL.
!-- Later create an INTEGER*8 header section.
    tmp = value

    CALL wrd_mpi_io( name, tmp )

 END SUBROUTINE wrd_mpi_io_int8



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> To do: Description missing!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_real( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name   !<

    REAL(wp), INTENT(IN)         ::  value  !<


    IF ( header_real_index == max_nr_real )  THEN
       STOP '+++ maximum number of REAL entries in restart file header exceeded'
    ENDIF

    real_names(header_real_index)  = name
    real_values(header_real_index) = value
    header_real_index = header_real_index + 1

 END SUBROUTINE wrd_mpi_io_real



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 2d-REAL array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_real_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name     !< name of array to be written

    INTEGER(iwp)                       ::  i        !< loop index
    INTEGER(iwp)                       ::  shift_i  !< shift of i index
    INTEGER(iwp)                       ::  shift_j  !< shift of j index

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status   !< returned status of MPI call
#endif

    REAL(wp), INTENT(IN), DIMENSION(:,:) ::  data   !< array to be written


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_real_2d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    IF ( SIZE( data, 1 ) == ( nyn - nys + 1 + 2 * nbgp )  .AND.                                    &
         SIZE( data, 2 ) == ( nxr - nxl + 1 + 2 * nbgp ) )  THEN
!
!--    Array with ghost layers.
!--    This kind of array is dimensioned in the caller subroutine as
!--    REAL, DIMENSION(nysg:nyng,nxlg:nxrg) ::  data
!
!--    Because the dimension bounds of 'data' are not explicitly given, they start at 1. To match
!--    the indices of data to the correct indices of array_2d, they must be shifted.
!--    Index nysg in array_2d corresponds to index 1 in data. Hence, to get the correct position
!--    inside the data array, the index "nysg" must be shifted by "-nysg+1". The same applies for x.
       shift_i = - nxlg + 1
       shift_j = - nysg + 1

       IF ( total_domain_boundaries_on_file )  THEN
!
!--       Prepare output with outer boundaries
          DO  i = iog%nxl, iog%nxr
             array_2d(i,iog%nys:iog%nyn) =                                                         &
                data(iog%nys-nbgp+shift_j:iog%nyn-nbgp+shift_j,i-nbgp+shift_i)
          ENDDO

       ELSE
!
!--       Prepare output without outer boundaries
          DO  i = nxl, nxr
             array_2d(i,iog%nys:iog%nyn) = data(nys+shift_j:nyn+shift_j,i+shift_i)
          ENDDO

       ENDIF

#if defined( __parallel )
       CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
       IF ( sm_io%iam_io_pe )  THEN
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_2d, 'native', MPI_INFO_NULL,    &
                                  ierr )
          CALL MPI_FILE_WRITE_ALL( fh, array_2d, SIZE( array_2d), MPI_REAL, status, ierr )
       ENDIF
       CALL sm_io%sm_node_barrier()
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_write( fh, array_2d, SIZE( array_2d ) )
#endif
!
!--    Type conversion required, otherwise right hand side brackets are calculated assuming 4 byte
!--    INT. Maybe a compiler problem.
       array_position = array_position + ( INT( iog%ny, KIND=rd_offset_kind ) + 1 ) *              &
                                         ( INT( iog%nx, KIND=rd_offset_kind ) + 1 ) * wp

!     ELSEIF ( SIZE( data, 1 ) == ( nyn - nys + 1 )  .AND.  &
!              SIZE( data, 2 ) == ( nxr - nxl + 1 ) )  THEN
!
!--    Array without ghost layers are not implemented yet.
!--    This kind of array is dimensioned in the caller subroutine as
!--    REAL, DIMENSION(nys:nyn,nxl:nxr) ::  data

    ELSE

       message_string = '2d-REAL array "' // TRIM( name ) // '" to be written to restart ' //   &
                        'file is defined with illegal dimensions'
       CALL message( 'wrd_mpi_io_real_2d', 'PAC0299', 3, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE wrd_mpi_io_real_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 2d-INTEGER array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                  ::  name    !<

    INTEGER(iwp)                                  ::  i       !<
    INTEGER(iwp)                                  ::  j       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size)            ::  status  !<
#endif
    INTEGER(KIND=iwp), INTENT(IN), DIMENSION(:,:) ::  data    !<


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_int_2d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    IF ( ( nxr - nxl + 1 + 2 * nbgp ) == SIZE( data, 2 ) )  THEN
!
!--    Integer arrays with ghost layers are not implemented yet. These kind of arrays would be
!--    dimensioned in the caller subroutine as
!--    INTEGER, DIMENSION(nysg:nyng,nxlg:nxrg) ::  data
       message_string = '2d-INTEGER array "' // TRIM( name ) // '" to be written to restart ' //   &
                        'file is defined with illegal dimensions'
       CALL message( 'wrd_mpi_io_int_2d', 'PAC0301', 3, 2, 0, 6, 0 )

    ELSEIF ( ( nxr-nxl+1 ) == SIZE( data, 2 ) )  THEN
!
!--    INTEGER input array without ghost layers.
!--    This kind of array is dimensioned in the caller subroutine as
!--    INTEGER, DIMENSION(nys:nyn,nxl:nxr) ::  data
       DO  j = nys, nyn
          DO  i = nxl, nxr
             array_2di(i,j) = data(j-nys+1,i-nxl+1)
          ENDDO
       ENDDO
#if defined( __parallel )
       CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
       IF ( sm_io%iam_io_pe )  THEN
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER, ft_2di_nb, 'native',            &
                                  MPI_INFO_NULL, ierr )
          CALL MPI_FILE_WRITE_ALL( fh, array_2di, SIZE( array_2di ), MPI_INTEGER, status, ierr )
       ENDIF
       CALL sm_io%sm_node_barrier()
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_write( fh, array_2di, SIZE( array_2di ) )
#endif
!
!--    Type conversion required, otherwise rigth hand side brackets are calculated assuming 4 byte
!--    INT. Maybe a compiler problem.
       array_position = array_position + INT( (ny+1), KIND=rd_offset_kind ) *                      &
                                         INT( (nx+1), KIND=rd_offset_kind ) * 4

    ELSE

       message_string = '2d-INTEGER array "' // TRIM( name ) // '" to be written to restart ' //   &
                        'file is defined with illegal dimensions'
       CALL message( 'wrd_mpi_io_int_2d', 'PAC0301', 3, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE wrd_mpi_io_int_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 2d-INTEGER array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int8_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                  ::  name    !<

    INTEGER(iwp)                                  ::  i       !<
    INTEGER(iwp)                                  ::  j       !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size)            ::  status  !<
#endif
    INTEGER(KIND=idp), INTENT(IN), DIMENSION(:,:) ::  data    !<


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_int8_2d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    IF ( ( nxr - nxl + 1 + 2 * nbgp ) == SIZE( data, 2 ) )  THEN
!
!--    Integer arrays with ghost layers are not implemented yet. These kind of arrays would be
!--    dimensioned in the caller subroutine as
!--    INTEGER, DIMENSION(nysg:nyng,nxlg:nxrg) ::  data
       message_string = '2d-INTEGER*8 array "' // TRIM( name ) // '" to be written to restart ' // &
                        'file is defined with illegal dimensions'
       CALL message( 'wrd_mpi_io_int8_2d', 'PAC0301', 3, 2, 0, 6, 0 )

    ELSEIF ( ( nxr-nxl+1 ) == SIZE( data, 2 ) )  THEN
!
!--    INTEGER input array without ghost layers.
!--    This kind of array is dimensioned in the caller subroutine as
!--    INTEGER, DIMENSION(nys:nyn,nxl:nxr) ::  data
       DO  j = nys, nyn
          DO  i = nxl, nxr
             array_2di8(i,j) = data(j-nys+1,i-nxl+1)
          ENDDO
       ENDDO
#if defined( __parallel )
       CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
       IF ( sm_io%iam_io_pe )  THEN
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER8, ft_2di8_nb, 'native',          &
                                  MPI_INFO_NULL, ierr )
          CALL MPI_FILE_WRITE_ALL( fh, array_2di8, SIZE( array_2di8 ), MPI_INTEGER8, status, ierr )
       ENDIF
       CALL sm_io%sm_node_barrier()
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_write( fh, array_2di8, SIZE( array_2di8 ) )
#endif
!
!--    Type conversion required, otherwise rigth hand side brackets are calculated assuming 4 byte
!--    INT. Maybe a compiler problem.
       array_position = array_position + INT( (ny+1), KIND=rd_offset_kind ) *                      &
                                         INT( (nx+1), KIND=rd_offset_kind ) * idp

    ELSE

       message_string = '2d-INTEGER*8 array "' // TRIM( name ) // '" to be written to restart ' // &
                        'file is defined with illegal dimensions'
       CALL message( 'wrd_mpi_io_int8_2d', 'PAC0301', 3, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE wrd_mpi_io_int8_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3d-INTEGER*4 array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int4_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<
#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif
    INTEGER(isp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data !<


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_int4_3d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    IF ( total_domain_boundaries_on_file )  THEN
!
!--    Prepare output of 3d-REAL-array with ghost layers. In the virtual PE grid, the first
!--    dimension is PEs along x, and the second along y. For MPI-IO it is recommended to have the
!--    index order of the array in the same way, i.e. the first dimension should be along x and the
!--    second along y. For this reason, the original PALM data need to be swaped.
       DO  i = iog%nxl, iog%nxr
          array_3di4(:,i,iog%nys:iog%nyn) = data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp)
       ENDDO

    ELSE
!
!--    Prepare output of 3d-REAL-array without ghost layers
       DO  i = nxl, nxr
           array_3di4(:,i,iog%nys:iog%nyn) = data(:,nys:nyn,i)
       ENDDO

    ENDIF
#if defined( __parallel )
    CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER, ft_3di4, 'native', MPI_INFO_NULL, ierr )
       CALL MPI_FILE_WRITE_ALL( fh, array_3di4, SIZE( array_3di4 ), MPI_INTEGER, status, ierr )
    ENDIF
    CALL sm_io%sm_node_barrier()
#else
    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, array_3di4, SIZE( array_3di4 ) )
#endif
!
!-- Type conversion required, otherwise right hand side brackets are calculated assuming 4 byte INT.
!-- Maybe a compiler problem.
    array_position = array_position + INT(     (nz+2), KIND = rd_offset_kind ) *                   &
                                      INT( (iog%ny+1), KIND = rd_offset_kind ) *                   &
                                      INT( (iog%nx+1), KIND = rd_offset_kind ) * isp

 END SUBROUTINE wrd_mpi_io_int4_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3d-INTEGER*4 array with MPI-IO
!> Dimension order: y, x, special
!> Dimension x and y do not contain ghost layers. The 3rd dimension is of size 'size_3rd_dimension'
!> and must be smaller than (nzt+1) - nzb + 1.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int4_3d_special( name, data, size_3rd_dimension )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name            !< name of variable to be written

    INTEGER(iwp), INTENT(IN) ::  size_3rd_dimension  !< size of 3rd dimension of array 'data'
    INTEGER(iwp)             ::  n                   !< loop index

    INTEGER(isp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr,1:size_3rd_dimension) ::  data  !< array to be written

    INTEGER(isp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  out_buf  !< resorted data array


    IF ( size_3rd_dimension > nzt - nzb + 2 )  THEN
       WRITE(message_string, '(A,I6,A,I6,A)')  '3rd dimension of array "' // TRIM( name ) //       &
                                             '" is larger than number of domain height levels:&',  &
                                             size_3rd_dimension, ' > ', nzt - nzb + 1,             &
                                             ', reading impossible'
       CALL message( 'wrd_mpi_io_int4_3d_special', 'PAC0304', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Reorder data from input array to default output buffer
    out_buf = 0_isp

    DO  n = 1, size_3rd_dimension
       out_buf(nzb+n-1,nys:nyn,nxl:nxr) = data(nys:nyn,nxl:nxr,n)
    ENDDO

!
!-- Write data to file
    CALL wrd_mpi_io_int4_3d( name, out_buf )

 END SUBROUTINE wrd_mpi_io_int4_3d_special



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3d-INTEGER*8 array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_int8_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<
#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif
    INTEGER(idp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data !<


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_int8_3d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    IF ( total_domain_boundaries_on_file )  THEN
!
!--    Prepare output of 3d-REAL-array with ghost layers. In the virtual PE grid, the first
!--    dimension is PEs along x, and the second along y. For MPI-IO it is recommended to have the
!--    index order of the array in the same way, i.e. the first dimension should be along x and the
!--    second along y. For this reason, the original PALM data need to be swaped.
       DO  i = iog%nxl, iog%nxr
          array_3di8(:,i,iog%nys:iog%nyn) = data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp)
       ENDDO

    ELSE
!
!--    Prepare output of 3d-REAL-array without ghost layers
       DO  i = nxl, nxr
           array_3di8(:,i,iog%nys:iog%nyn) = data(:,nys:nyn,i)
       ENDDO

    ENDIF
#if defined( __parallel )
    CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_INTEGER8, ft_3di8, 'native', MPI_INFO_NULL, ierr )
       CALL MPI_FILE_WRITE_ALL( fh, array_3di8, SIZE( array_3di8 ), MPI_INTEGER8, status, ierr )
    ENDIF
    CALL sm_io%sm_node_barrier()
#else
    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, array_3di8, SIZE( array_3di8 ) )
#endif
!
!-- Type conversion required, otherwise right hand side brackets are calculated assuming 4 byte INT.
!-- Maybe a compiler problem.
    array_position = array_position + INT(     (nz+2), KIND = rd_offset_kind ) *                   &
                                      INT( (iog%ny+1), KIND = rd_offset_kind ) *                   &
                                      INT( (iog%nx+1), KIND = rd_offset_kind ) * dp

 END SUBROUTINE wrd_mpi_io_int8_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3d-REAL array with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_real_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<
#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif
    REAL(wp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data !<


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_real_3d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    IF ( total_domain_boundaries_on_file )  THEN
!
!--    Prepare output of 3d-REAL-array with ghost layers. In the virtual PE grid, the first
!--    dimension is PEs along x, and the second along y. For MPI-IO it is recommended to have the
!--    index order of the array in the same way, i.e. the first dimension should be along x and the
!--    second along y. For this reason, the original PALM data need to be swaped.
       DO  i = iog%nxl, iog%nxr
          array_3d(:,i,iog%nys:iog%nyn) = data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp)
       ENDDO

    ELSE
!
!--    Prepare output of 3d-REAL-array without ghost layers
       DO  i = nxl, nxr
           array_3d(:,i,iog%nys:iog%nyn) = data(:,nys:nyn,i)
       ENDDO

    ENDIF
#if defined( __parallel )
    CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_3d, 'native', MPI_INFO_NULL, ierr )
       CALL MPI_FILE_WRITE_ALL( fh, array_3d, SIZE( array_3d ), MPI_REAL, status, ierr )
    ENDIF
    CALL sm_io%sm_node_barrier()
#else
    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, array_3d, SIZE( array_3d ) )
#endif
!
!-- Type conversion required, otherwise right hand side brackets are calculated assuming 4 byte INT.
!-- Maybe a compiler problem.
    array_position = array_position + INT(     (nz+2), KIND = rd_offset_kind ) *                   &
                                      INT( (iog%ny+1), KIND = rd_offset_kind ) *                   &
                                      INT( (iog%nx+1), KIND = rd_offset_kind ) * wp

    IF ( debug_output )  THEN
       WRITE( debug_string, '(3A,I10)' ) 'array_position real3d ', TRIM( name ), ' ', array_position
       CALL debug_message( debug_string, 'info' )
    ENDIF

 END SUBROUTINE wrd_mpi_io_real_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3d-REAL soil array with MPI-IO.
!> nzb_soil, nzt_soil are located in the module land_surface_model_mod. Since Fortran does not allow
!> cross referencing of module variables, it is required to pass these variables as arguments.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_real_3d_soil( name, data, nzb_soil, nzt_soil )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name      !<

    INTEGER(iwp)                       ::  i         !<
    INTEGER, INTENT(IN)                ::  nzb_soil  !<
    INTEGER, INTENT(IN)                ::  nzt_soil  !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    REAL(wp), INTENT(IN), DIMENSION(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) ::  data  !<


    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_real_3d_soil', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    CALL rd_mpi_io_create_filetypes_3dsoil( nzb_soil, nzt_soil )

    IF ( total_domain_boundaries_on_file)  THEN
!
!--    Prepare output of 3d-REAL-array with ghost layers. In the virtual PE grid, the first
!--    dimension is PEs along x, and the second along y. For MPI-IO it is recommended to have the
!--    index order of the array in the same way, i.e. the first dimension should be along x and the
!--    second along y. For this reason, the original PALM data need to be swaped.
       DO  i = iog%nxl, iog%nxr
          array_3d_soil(:,i,iog%nys:iog%nyn) = data(:,iog%nys-nbgp:iog%nyn-nbgp,i-nbgp)
       ENDDO

    ELSE
!
!--    Prepare output of 3d-REAL-array without ghost layers
       DO  i = nxl, nxr
          array_3d_soil(:,i,iog%nys:iog%nyn) = data(:,nys:nyn,i)
       ENDDO

    ENDIF
#if defined( __parallel )
    CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_3dsoil, 'native', MPI_INFO_NULL,   &
                               ierr )
       CALL MPI_FILE_WRITE_ALL( fh, array_3d_soil, SIZE( array_3d_soil ), MPI_REAL, status, ierr )
    ENDIF
    CALL sm_io%sm_node_barrier()

    IF ( sm_io%is_sm_active() )  THEN
       CALL sm_io%sm_free_shared( win_3ds )
    ELSE
       DEALLOCATE( array_3d_soil )
    ENDIF
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_FREE( ft_3dsoil, ierr )
    ENDIF
#else
    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, array_3d_soil, SIZE( array_3d_soil ) )
    DEALLOCATE( array_3d_soil )
#endif
!
!-- Type conversion required, otherwise right hand side brackets are calculated assuming 4 byte INT.
!-- Maybe a compiler problem.
    array_position = array_position + INT( (nzt_soil-nzb_soil+1), KIND = rd_offset_kind ) *        &
                                      INT( (iog%ny+1),            KIND = rd_offset_kind ) *        &
                                      INT( (iog%nx+1),            KIND = rd_offset_kind ) * wp

 END SUBROUTINE wrd_mpi_io_real_3d_soil


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write CHARATCTER with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_char( name, text )

    IMPLICIT NONE

    CHARACTER(LEN=128)           ::  lo_line  !<
    CHARACTER(LEN=*), INTENT(IN) ::  name     !<
    CHARACTER(LEN=*), INTENT(IN) ::  text     !<


    IF ( header_char_index == max_nr_char )  THEN
       STOP '+++ maximum number of CHARACTER entries in restart file header exceeded'
    ENDIF

    lo_line      = name
    lo_line(33:) = text
    text_lines(header_char_index) = lo_line
    header_char_index = header_char_index + 1

 END SUBROUTINE wrd_mpi_io_char



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write LOGICAL with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_logical( name, value )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name                !<

    INTEGER(iwp)                 ::  logical_as_integer  !<

    LOGICAL, INTENT(IN)          ::  value               !<


    IF ( value )  THEN
       logical_as_integer = 1
    ELSE
       logical_as_integer = 0
    ENDIF

    CALL wrd_mpi_io_int( name, logical_as_integer )

 END SUBROUTINE wrd_mpi_io_logical



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3D-LOGICAL with MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_logical_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name                !<

    INTEGER(iwp) ::  i  !< loop index
    INTEGER(iwp) ::  j  !< loop index
    INTEGER(iwp) ::  k  !< loop index

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  logical_as_integer  !<

    LOGICAL, INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  data   !<


    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          DO  k = nzb, nzt+1
             IF ( data(k,j,i) )  THEN
                logical_as_integer(k,j,i) = 1
             ELSE
                logical_as_integer(k,j,i) = 0
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    CALL wrd_mpi_io_int4_3d( name, logical_as_integer )

 END SUBROUTINE wrd_mpi_io_logical_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 1d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_global_array_real_1d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)       ::  name    !<

    INTEGER(iwp)                       ::  i       !<
    INTEGER(KIND=rd_offset_kind)       ::  offset  !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status  !<
#endif

    LOGICAL                            ::  found   !<

    REAL(KIND=wp), INTENT(INOUT), DIMENSION(:) ::  data  !<


    offset = 0
    found  = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO


    IF ( found )  THEN

!
!--    Set default view
#if defined( __parallel )
       IF ( cyclic_fill_mode )  THEN
          IF ( pe_active_for_read )  THEN
             CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
             CALL MPI_FILE_SEEK( fh, array_position, MPI_SEEK_SET, ierr )
             CALL MPI_FILE_READ_ALL( fh, data, SIZE( data ), MPI_REAL, status, ierr )
         ENDIF
         CALL MPI_BCAST( data, SIZE( data ), MPI_REAL, 0, comm2d, ierr )
       ELSE
          IF( sm_io%iam_io_pe )  THEN
             CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
          ENDIF
          IF ( myid == 0 )  THEN
             CALL MPI_FILE_SEEK( fh, array_position, MPI_SEEK_SET, ierr )
             CALL MPI_FILE_READ( fh, data, SIZE( data ), MPI_REAL, status, ierr )
          ENDIF
          CALL MPI_BCAST( data, SIZE( data ), MPI_REAL, 0, comm2d, ierr )
       ENDIF
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_read( fh, data, SIZE( data ) )
#endif

    ELSE

       message_string = '1d/2d/3d/4d-REAL global array "' // TRIM( name ) //                       &
                        '" not found in restart file'
       CALL message( 'rrd_mpi_io_global_array_real_1d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE rrd_mpi_io_global_array_real_1d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_global_array_real_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                      ::  name      !<

    INTEGER, DIMENSION(1)                             ::  bufshape  !<

    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), TARGET ::  data      !<
    REAL(KIND=wp), POINTER, DIMENSION(:)              ::  buf       !<

    TYPE(C_PTR)                                       ::  c_data    !<


    c_data = C_LOC( data )
    bufshape(1) = SIZE( data )
    CALL C_F_POINTER( c_data, buf, bufshape )

    CALL rrd_mpi_io_global_array_real_1d( name, buf )

 END SUBROUTINE rrd_mpi_io_global_array_real_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 3d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_global_array_real_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                        ::  name      !<

    INTEGER, DIMENSION(1)                               ::  bufshape  !<

    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), TARGET ::  data      !<
    REAL(KIND=wp), POINTER, DIMENSION(:)                ::  buf       !<

    TYPE(C_PTR)                                         ::  c_data    !<


    c_data = C_LOC( data )
    bufshape(1) = SIZE( data )
    CALL C_F_POINTER( c_data, buf, bufshape )

    CALL rrd_mpi_io_global_array_real_1d( name, buf )

 END SUBROUTINE rrd_mpi_io_global_array_real_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 4d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_global_array_real_4d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                          ::  name      !<

    INTEGER, DIMENSION(1)                                 ::  bufshape  !<

    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:,:), TARGET ::  data      !<
    REAL(KIND=wp), POINTER, DIMENSION(:)                  ::  buf       !<

    TYPE(C_PTR)                                           ::  c_data    !<


    c_data = C_LOC( data )
    bufshape(1) = SIZE( data)
    CALL C_F_POINTER( c_data, buf, bufshape )

    CALL rrd_mpi_io_global_array_real_1d( name, buf )

 END SUBROUTINE rrd_mpi_io_global_array_real_4d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 1d-INTEGER global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_global_array_int_1d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                   ::  name    !<

    INTEGER(iwp)                                   ::  i       !<
    INTEGER(KIND=rd_offset_kind)                   ::  offset  !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size)             ::  status  !<
#endif
    INTEGER(KIND=iwp), INTENT(INOUT), DIMENSION(:) ::  data    !<

    LOGICAL                                        ::  found   !<


    offset = 0
    found  = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM( array_names(i) ) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN
!
!--    Set default view
#if defined( __parallel )
       IF ( cyclic_fill_mode )  THEN
          IF ( pe_active_for_read )  THEN
             CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
             CALL MPI_FILE_SEEK( fh, array_position, MPI_SEEK_SET, ierr )
             CALL MPI_FILE_READ_ALL( fh, data, SIZE( data), MPI_INTEGER, status, ierr )
          ENDIF
          CALL MPI_BCAST( data, SIZE( data ), MPI_INTEGER, 0, comm2d, ierr )
       ELSE
          IF( sm_io%iam_io_pe )  THEN
             CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
          ENDIF
          IF ( myid == 0 )  THEN
             CALL MPI_FILE_SEEK( fh, array_position, MPI_SEEK_SET, ierr )
             CALL MPI_FILE_READ( fh, data, SIZE( data), MPI_INTEGER, status, ierr )
          ENDIF
          CALL MPI_BCAST( data, SIZE( data ), MPI_INTEGER, 0, comm2d, ierr )
        ENDIF
#else
       CALL posix_lseek( fh, array_position )
       CALL posix_read( fh, data, SIZE( data ) )
#endif
    ELSE

       message_string = '1d-INTEGER global array "' // TRIM( name ) // '" not found in restart file'
       CALL message( 'rrd_mpi_io_global_array_int_1d', 'PAC0298', 3, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE rrd_mpi_io_global_array_int_1d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 1d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_global_array_real_1d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)            ::  name    !<

    INTEGER(KIND=rd_offset_kind)            ::  offset  !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size)      ::  status  !<
#endif

    REAL(KIND=wp), INTENT(IN), DIMENSION(:) ::  data    !<


    offset = 0

    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_global_array_real_1d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

!
!-- Set default view
#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
    ENDIF
!
!-- Only PE 0 writes replicated data
    IF ( myid == 0 )  THEN
       CALL MPI_FILE_SEEK( fh, array_position, MPI_SEEK_SET, ierr )
       CALL MPI_FILE_WRITE( fh, data, SIZE( data ), MPI_REAL, status, ierr )
    ENDIF
#else
    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, data, SIZE( data ) )
#endif
    array_position = array_position + SIZE( data ) * wp

 END SUBROUTINE wrd_mpi_io_global_array_real_1d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 2d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_global_array_real_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                      ::  name      !<

    INTEGER, DIMENSION(1)                             ::  bufshape  !<

    REAL(KIND=wp), POINTER, DIMENSION(:)              ::  buf       !<
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), TARGET ::  data      !<

    TYPE(C_PTR)                                       ::  c_data    !<


    c_data = C_LOC( data )
    bufshape(1) = SIZE( data)
    CALL C_F_POINTER( c_data, buf, bufshape )

    CALL wrd_mpi_io_global_array_real_1d( name, buf )

 END SUBROUTINE wrd_mpi_io_global_array_real_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 3d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_global_array_real_3d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                        ::  name      !<

    INTEGER, DIMENSION(1)                               ::  bufshape  !<

    REAL(KIND=wp), POINTER, DIMENSION(:)                ::  buf       !<
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), TARGET ::  data      !<

    TYPE(C_PTR)                                         ::  c_data    !<


    c_data = C_LOC( data )
    bufshape(1) = SIZE( data )
    CALL C_F_POINTER( c_data, buf, bufshape )

    CALL wrd_mpi_io_global_array_real_1d( name, buf )

 END SUBROUTINE wrd_mpi_io_global_array_real_3d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 4d-REAL global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_global_array_real_4d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                          ::  name      !<

    INTEGER, DIMENSION(1)                                 ::  bufshape  !<

    REAL(KIND=wp), POINTER, DIMENSION(:)                  ::  buf       !<
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:,:), TARGET ::  data      !<

    TYPE(C_PTR)                                           ::  c_data    !<


    c_data = C_LOC( data )
    bufshape(1) = SIZE( data)
    CALL C_F_POINTER( c_data, buf, bufshape )

    CALL wrd_mpi_io_global_array_real_1d( name, buf )

 END SUBROUTINE wrd_mpi_io_global_array_real_4d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 1d-INTEGER global array with MPI-IO.
!> Array contains identical data on all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_global_array_int_1d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                ::  name    !<

    INTEGER(KIND=rd_offset_kind)                ::  offset  !<

    INTEGER(KIND=iwp), INTENT(IN), DIMENSION(:) ::  data    !<
#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size)          ::  status  !<
#endif

    IF ( header_array_index == max_nr_arrays )  THEN
       message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
       CALL message( 'wrd_mpi_io_global_array_int_1d', 'PAC0300', 1, 2, 0, 6, 0 )
    ENDIF

    offset = 0
    array_names(header_array_index)  = name
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

!
!-- Set default view
#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
    ENDIF
!
!-- Only PE 0 writes replicated data
    IF ( myid == 0 )  THEN
       CALL MPI_FILE_SEEK( fh, array_position, MPI_SEEK_SET, ierr )
       CALL MPI_FILE_WRITE( fh, data, SIZE( data), MPI_INTEGER, status, ierr )
    ENDIF
#else
    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, data, SIZE( data ) )
#endif
    array_position = array_position + SIZE( data ) * 4

 END SUBROUTINE wrd_mpi_io_global_array_int_1d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read particle data array with MPI-IO.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_particles( name, prt_global_index )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                       ::  name            !<
    INTEGER(idp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  prt_global_index      !<

    INTEGER(iwp)                       ::  array_size      !<
    INTEGER(iwp)                       ::  byte_column     !<
    INTEGER(iwp)                       ::  i               !<
    INTEGER(iwp)                       ::  ind             !<
    INTEGER(iwp)                       ::  j               !<
    INTEGER(iwp)                       ::  k               !<
    INTEGER(iwp)                       ::  n               !<
    INTEGER(iwp)                       ::  particle_size   !<

    INTEGER(KIND=rd_offset_kind)       ::  disp            !<
    INTEGER(KIND=rd_offset_kind)       ::  offset          !<
    INTEGER(KIND=rd_offset_kind)       ::  prt_nr_bytes    !<

    LOGICAL                            ::  found           !<

    REAL(dp)                           :: rr               !< there is no data type INTEGER*8 in MPI
    REAL(dp)                           :: rs               !< use REAL*8 to compute max offset

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE, TARGET :: prt_data   !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status          !<
#else
    TYPE(C_PTR)                        ::  buf
#endif

    found = .FALSE.

    DO  i = 1, tgh%nr_arrays
       IF ( TRIM(array_names(i)) == TRIM( name ) )  THEN
          array_position = array_offset(i)
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF ( found )  THEN

       offset = 0

       particle_size = STORAGE_SIZE(zero_particle) / 8  ! 8 here means number of bits/byte, NOT wp

       array_size = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             array_size = MAX( array_size, SUM(prt_count(:,j,i)) )
          ENDDO
       ENDDO

       ALLOCATE( prt_data(MAX(array_size,1)) )

!
!--    Write columns of particle
#if defined( __parallel )
       CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
#endif
       prt_nr_bytes = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             disp         = array_position + prt_global_index(nzb,j,i) * particle_size
             byte_column  = SUM( prt_count(:,j,i) ) * particle_size
             prt_nr_bytes = MAX( disp+byte_column, prt_nr_bytes )

#if defined( __parallel )
             CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
             IF ( byte_column > 0 )  THEN
                CALL MPI_FILE_SEEK( fh, disp, MPI_SEEK_SET, ierr )
                CALL MPI_FILE_READ( fh, prt_data, byte_column, MPI_BYTE, status, ierr )
             ENDIF
             CALL sm_io%sm_node_barrier()
#else
             buf = C_LOC(prt_data)     ! use C_PTR to avaid another overlay in posix interface
             CALL posix_lseek( fh, disp )
             CALL posix_read( fh, buf, byte_column )
#endif
             ind = 1
             DO  k = nzb, nzt+1
                DO  n = 1, prt_count(k,j,i)
                   grid_particles(k,j,i)%particles(n) = prt_data(ind)
                   ind = ind+1
                ENDDO
             ENDDO

          ENDDO
       ENDDO

#if defined( __parallel )
       rs = prt_nr_bytes
       CALL MPI_ALLREDUCE( rs, rr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm2d, ierr )
       prt_nr_bytes = rr
#else
       rr = rs
#endif
       array_position = prt_nr_bytes

       DEALLOCATE( prt_data )

    ENDIF

 END SUBROUTINE rrd_mpi_io_particles



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 1d-REAL surface data array with MPI-IO.
!> This is a recursive subroutine. In case of cyclic fill mode it may call itself for reading parts
!> of the prerun grid.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE rrd_mpi_io_surface( name, data, first_index )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name            !<

    INTEGER(iwp), OPTIONAL       ::  first_index     !<
    INTEGER(idp)                 ::  i               !<
    INTEGER(idp)                 ::  j               !<
    INTEGER(iwp)                 ::  lo_first_index  !<

#if defined( __parallel )
    INTEGER(idp)                 ::  buf_start       !<
    INTEGER(KIND=rd_offset_kind) ::  disp            !< displacement of actual indices
    INTEGER(idp)                 ::  ie              !<
    INTEGER(idp)                 ::  ind_gb          !<
    INTEGER(idp)                 ::  ind_out         !<
    INTEGER(idp)                 ::  is              !<
    INTEGER(iwp)                 ::  n               !<
    INTEGER(iwp)                 ::  n_trans         !<
    INTEGER(KIND=rd_offset_kind) ::  offset          !<

    INTEGER(iwp),DIMENSION(0:numprocs-1) ::  lo_index  !<
    INTEGER, DIMENSION(rd_status_size)   ::  status    !<
#endif
    LOGICAL                      ::  found  !<

    REAL(wp), INTENT(OUT), DIMENSION(:), TARGET ::  data  !<
#if defined( __parallel )
    REAL(wp),DIMENSION(:),ALLOCATABLE    ::  put_buffer  !<
#endif


    found = .FALSE.
    lo_first_index = 1

    IF ( PRESENT( first_index ) )   THEN
       lo_first_index = first_index
    ENDIF

    DO  i = 1, tgh%nr_arrays
        IF ( TRIM( array_names(i) ) == TRIM( name ) )  THEN
!
!--        ATTENTION: The total_number_of_surface_elements and wp MUST be INTERGER*8.
!--        The compiler (at least Intel) first computes total_number_of_surface_elements*wp
!--        and then does the conversion to INTEGER(8).
!--        This may lead to wrong results when total_number_of_surface_elements*wp is > 2*10**6
           array_position = array_offset(i) + ( lo_first_index - 1 ) *                             &
                            INT( total_number_of_surface_elements, idp ) * INT( wp, idp )
           found = .TRUE.
           EXIT
        ENDIF
    ENDDO

!
!-- In case of 2d-data, name is written only once
    IF ( lo_first_index == 1 )  THEN

       array_names(header_array_index)  = name
       array_offset(header_array_index) = array_position
       header_array_index = header_array_index + 1

    ENDIF

    IF ( found )  THEN
       IF ( cyclic_fill_mode )  THEN

          CALL rrd_mpi_io_surface_cyclic_fill
          RETURN

       ELSE
#if defined( __parallel )
!
!--       Read data from restart file
          CALL sm_io%sm_node_barrier() ! has no effect if I/O on limited number of cores is inactive
          IF ( sm_io%iam_io_pe )  THEN
             IF ( less_than_2g_elements )  THEN
                CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_surf, 'native',           &
                                        MPI_INFO_NULL, ierr )
             ELSE
                offset =  array_position+local_offset
                CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native',                  &
                                        MPI_INFO_NULL, ierr )
             ENDIF
             CALL MPI_FILE_READ_ALL( fh, array_out, SIZE( array_out ), MPI_REAL, status, ierr )
          ENDIF
          CALL sm_io%sm_node_barrier()

!
!--       Copy data into transfer buffer. Data is organized in a way that only one MPI_PUT to the
!--       respective PE ist required.
          ALLOCATE( put_buffer(SUM( transfer_index(4,:) )) )
!
!--       Copy into to transfer buffer. Here, array_out matches the indices global_start and
!--       global_end, reduced by array_out_reduction.
!--       The use of array_out_reduction allows the use 32 bit indices in sm_allocate_shared.
          ind_gb = 1
          DO  i = 1, SIZE( local_indices, 2 )
             ind_out = local_indices(1,i) - array_out_reduction
             DO  j = 1, local_indices(2,i)
                put_buffer(ind_gb) = array_out(ind_out)
                ind_out = ind_out + 1
                ind_gb  = ind_gb  + 1
             ENDDO
          ENDDO
!
!--       Transfer data from I/O PEs to the respective PEs to which they belong.
          CALL MPI_WIN_FENCE( 0, win_surf, ierr )

          buf_start = 1
          DO  n = 0, numprocs-1
             n_trans = transfer_index(4,n)
             IF ( n_trans > 0 )  THEN
                disp = transfer_index(3,n) - 1
                CALL MPI_PUT( put_buffer(buf_start), n_trans, MPI_REAL, n, disp, n_trans,          &
                              MPI_REAL, win_surf, ierr)
                buf_start = buf_start + n_trans
             ENDIF
          ENDDO

          CALL MPI_WIN_FENCE( 0, win_surf, ierr )
          DEALLOCATE( put_buffer )
!
!--       Copy from RMA window into output array (data) to allow transfering data to target PEs.
!--       Check, if the number of surface elements per horizontal gridbox match the index setup.
          lo_index = nr_surfaces_in_tb
          DO  i = nxl, nxr
             DO  j = nys, nyn
                is = lo_index(s_index_in_window(j,i)) + 1
                ie = is + m_end_index(j,i) - m_start_index(j,i)
                data(m_start_index(j,i):m_end_index(j,i)) = array_1d(is:ie)
                lo_index(s_index_in_window(j,i)) = lo_index(s_index_in_window(j,i)) +              &
                                                   e_end_index(j,i) - e_start_index(j,i) + 1
!
!--             Internal check. Should not happen and would mean programming error.
                IF ( e_end_index(j,i) - e_start_index(j,i) + 1 /= NINT( array_1d(is-1) ) )  THEN
                   WRITE( message_string, '(A,6I8)' ) 'number of surface elements in array_1d ' // &
                                    'does not match ', j, i, e_start_index(j,i), e_end_index(j,i), &
                                    e_end_index(j,i)-e_start_index(j,i)+1 , NINT( array_1d(is-1) )
                   CALL message( 'rrd_mpi_io_surface', 'PAC0305', 2, 2, 0, 6, 1 )
                ENDIF
             ENDDO
          ENDDO


#else
          CALL posix_lseek( fh, array_position )
          CALL posix_read( fh, array_out, SIZE(array_out) )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                data(m_start_index(j,i):m_end_index(j,i)) =                                        &
                                                   array_out(e_start_index(j,i)+1:e_end_index(j,i))
!
!--             Internal check. Should not happen and would mean programming error.
                IF ( e_end_index(j,i) - e_start_index(j,i) + 1 /=                                  &
                     NINT( array_out(e_start_index(j,i)) ) )                                       &
                THEN
                   WRITE( message_string, '(A,6I8)' ) 'number of surface elements in array_out ' //&
                      'does not match ', j, i, e_start_index(j,i), e_end_index(j,i),               &
                      e_end_index(j,i)-e_start_index(j,i)+1 , NINT( array_out(e_start_index(j,i)) )
                   CALL message( 'rrd_mpi_io_surface', 'PAC0305', 2, 2, 0, 6, 1 )
                ENDIF
             ENDDO
          ENDDO
#endif
       ENDIF

    ENDIF

 CONTAINS

    SUBROUTINE rrd_mpi_io_surface_cyclic_fill

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !<
       INTEGER(iwp) ::  ie        !<
       INTEGER(iwp) ::  is        !<
       INTEGER(iwp) ::  i_remote  !<
       INTEGER(iwp) ::  j         !<
       INTEGER(iwp) ::  je        !<
       INTEGER(iwp) ::  js        !<
       INTEGER(iwp) ::  j_remote  !<
       INTEGER(iwp) ::  nval      !<
       INTEGER(iwp) ::  rem_pe    !<

#if defined( __parallel )
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  rem_offs  !<
#else
       INTEGER(idp) ::  rem_offs                    !<
#endif

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_data  !<


!
!--    ATTENTION: This version allows only 1 surface element per grid cell.
!
!--    Activate grid of the smaller prerun, i.e. grid variables like nxl, nxr, nys, nyn and other
!--    values are set according to the prerun settings.
       CALL prerun_grid%activate_grid_from_this_class()

       IF ( pe_active_for_read )  THEN

          IF ( MAXVAL( m_end_index ) <= 0 )  THEN
             CALL mainrun_grid%activate_grid_from_this_class()
             IF ( debug_output )  THEN
                CALL debug_message( 'PE inactive for reading restart or prerun data', 'start' )
             ENDIF
             RETURN
          ENDIF

          ALLOCATE( c_data(MAXVAL( m_end_index )) )

!
!--       Recursive CALL of rrd_mpi_io_surface.
!--       rrd_mpi_io_surface is called with cyclic_fill_mode = .FALSE. on the smaller
!--       prerun grid.
          cyclic_fill_mode = .FALSE.
          CALL rrd_mpi_io_surface( name, c_data )
          cyclic_fill_mode = .TRUE.

          DO  i = nxl, nxr
             DO  j = nys, nyn
                rmabuf_2d(j,i) = c_data(c_start_index(j,i))
             ENDDO
          ENDDO

       ENDIF
!
!--    Activate grid of the mainrun, i.e. grid variables like nxl, nxr, nys, nyn and other values
!--    are set according to the mainrun settings.
       CALL mainrun_grid%activate_grid_from_this_class()

#if defined( __parallel )
!
!--    Close RMA window to allow remote access
       CALL MPI_WIN_FENCE( 0, rmawin_2d, ierr )
#endif

!
!--   After reading surface data on the small grid, map these data in a cyclic way to all respective
!--   grid points of the main run.
      is = nxl
      ie = nxr
      js = nys
      je = nyn

       DO  i = is, ie
          DO  j = js, je
             i_remote = MOD( i, nx_on_file+1 )
             j_remote = MOD( j, ny_on_file+1 )
             rem_pe   = remote_pe(i_remote,j_remote)
             rem_offs = rma_offset(i_remote,j_remote)
             nval     = 1

#if defined( __parallel )
             IF ( rem_pe /= myid )  THEN
                CALL MPI_GET( data(o_start_index(j,i)), nval, MPI_REAL, rem_pe, rem_offs, nval,    &
                              MPI_REAL, rmawin_2d, ierr)
             ELSE
                data(o_start_index(j,i)) = rmabuf_2d(j_remote,i_remote)
             ENDIF
#else
             data(o_start_index(j,i)) = array_2d(i_remote,j_remote)
#endif
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Reopen RMA window to allow filling
       CALL MPI_WIN_FENCE( 0, rmawin_2d, ierr )
#endif

       IF ( ALLOCATED( c_data ) )  DEALLOCATE( c_data )

    END SUBROUTINE rrd_mpi_io_surface_cyclic_fill

 END SUBROUTINE rrd_mpi_io_surface



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-REAL surface data array with MPI-IO.
!> These consist of multiple 1d-REAL surface data arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rrd_mpi_io_surface_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)          ::  name  !<

    INTEGER(iwp)                          ::  i     !<

    REAL(wp), INTENT(OUT), DIMENSION(:,:) ::  data  !<
    REAL(wp), DIMENSION(SIZE( data,2))    ::  tmp   !<


    DO  i = 1, SIZE( data, 1 )
       CALL rrd_mpi_io_surface( name, tmp, first_index = i )
       data(i,:) = tmp
    ENDDO

!
!-- Comment from Klaus Ketelsen (September 2018)
!-- The intention of the following loop was to let the compiler do the copying on return.
!-- In my understanding it is standard conform to pass the second dimension to a 1d-array inside a
!-- subroutine and the compiler is responsible to generate code for copying. Acually this works fine
!-- for INTENT(IN) variables (wrd_mpi_io_surface_2d). For INTENT(OUT) like in this case the code
!-- works on pgi compiler. But both, the Intel 16 and the Cray compiler show wrong answers using
!-- this loop. That is the reason why the above auxiliary array tmp was introduced.
!    DO  i = 1, SIZE(  data,1)
!       CALL rrd_mpi_io_surface( name, data(i,:), first_index = i )
!    ENDDO

 END SUBROUTINE rrd_mpi_io_surface_2d



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write particle data with MPI-IO.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_particles( name, prt_global_index )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                       ::  name            !<
    INTEGER(idp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  prt_global_index      !<

    INTEGER(iwp)                       ::  array_size      !<
    INTEGER(iwp)                       ::  byte_column     !<
    INTEGER(iwp)                       ::  i               !<
    INTEGER(iwp)                       ::  ind             !<
    INTEGER(iwp)                       ::  j               !<
    INTEGER(iwp)                       ::  k               !<
    INTEGER(iwp)                       ::  n               !<
    INTEGER(iwp)                       ::  particle_size   !<

    INTEGER(KIND=rd_offset_kind)       ::  disp            !<
    INTEGER(KIND=rd_offset_kind)       ::  offset          !<
    INTEGER(KIND=rd_offset_kind)       ::  prt_nr_bytes    !<

    REAL(dp)                           :: rs               !< use REAL*8 to compute max offset
    REAL(dp)                           :: rr               !< there is no data type INTEGER*8 in MPI


    TYPE(particle_type), DIMENSION(:), ALLOCATABLE, TARGET ::  prt_data   !<

#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status          !<
#else
    TYPE(C_PTR)                        ::  buf
#endif


    offset = 0

    array_names(header_array_index)  = TRIM(name)
    array_offset(header_array_index) = array_position
    header_array_index = header_array_index + 1

    particle_size = STORAGE_SIZE( zero_particle ) / 8

    array_size = 0
    DO  i = nxl, nxr
      DO  j = nys, nyn
         array_size = MAX( array_size, SUM(prt_count(:,j,i)) )
       ENDDO
    ENDDO

    ALLOCATE( prt_data(MAX(array_size,1)) )

!
!-- Write columns of particles.
!-- The particles of a column are stored sequentially in the first dimension of the particles array.
!-- Store only the particles of one cell would be possible with this setup, but the I/O portions
!-- for a maximum of 100 particles are not big enough.
#if defined( __parallel )
    CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
#endif
    prt_nr_bytes = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          disp         = array_position + prt_global_index(nzb,j,i) * particle_size
          byte_column  = SUM( prt_count(:,j,i) ) * particle_size
          prt_nr_bytes = MAX( disp+byte_column, prt_nr_bytes )

          ind = 1
          DO  k = nzb, nzt+1
             DO  n = 1, prt_count(k,j,i)
                prt_data(ind) = grid_particles(k,j,i)%particles(n)
                ind = ind+1
             ENDDO
          ENDDO

#if defined( __parallel )
          CALL sm_io%sm_node_barrier()  ! Has no effect if I/O on limited number of cores is inactive
          IF ( byte_column > 0 )  THEN
             CALL MPI_FILE_SEEK( fh, disp, MPI_SEEK_SET, ierr )
             CALL MPI_FILE_WRITE( fh, prt_data, byte_column, MPI_BYTE, status, ierr )
          ENDIF
          CALL sm_io%sm_node_barrier()
#else
          buf = C_LOC(prt_data)  ! use C_PTR to avoid another overlay in posix interface
          CALL posix_lseek( fh, disp )
          CALL posix_write( fh, buf, byte_column )
#endif
       ENDDO
    ENDDO

#if defined( __parallel )
    rs = prt_nr_bytes
    CALL MPI_ALLREDUCE( rs, rr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm2d, ierr )
    prt_nr_bytes = rr
#else
    rr = rs
#endif
    array_position = prt_nr_bytes

    DEALLOCATE( prt_data )

 END SUBROUTINE wrd_mpi_io_particles



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write 1d-REAL surface data array with MPI-IO.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_surface( name, data, first_index )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name  !<

    INTEGER(iwp), OPTIONAL ::  first_index     !<
    INTEGER(idp) ::  i               !<
    INTEGER(idp) ::  j               !<
    INTEGER(iwp) ::  lo_first_index  !<
#if defined( __parallel )
    INTEGER(idp) ::  buf_start       !<
    INTEGER(iwp) ::  ie              !<
    INTEGER(idp) ::  is              !<
    INTEGER(idp) ::  ind_gb          !<
    INTEGER(idp) ::  ind_out         !<
    INTEGER(iwp) ::  n               !<
    INTEGER(iwp) ::  n_trans         !<
#endif

#if defined( __parallel )
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  disp    !< displacement in RMA window
    INTEGER(KIND=rd_offset_kind)   ::  offset  !<

    INTEGER(iwp), DIMENSION(0:numprocs-1)   ::  lo_index  !<
    INTEGER(iwp), DIMENSION(rd_status_size) ::  status    !<
#endif

    REAL(wp), INTENT(IN), DIMENSION(:), TARGET ::  data  !<
#if defined( __parallel )
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  get_buffer  !<
#endif


    lo_first_index = 1

    IF ( PRESENT( first_index ) )  THEN
       lo_first_index = first_index
    ENDIF

!
!-- In case of 2d-data, name is written only once
    IF ( lo_first_index == 1 )  THEN

       IF ( header_array_index == max_nr_arrays )  THEN
          message_string = 'maximum number of 2d/3d-array entries in restart file header exceeded'
          CALL message( 'wrd_mpi_io_surface', 'PAC0300', 1, 2, 0, 6, 0 )
       ENDIF

       array_names(header_array_index)  = name
       array_offset(header_array_index) = array_position
       header_array_index = header_array_index + 1

    ENDIF

#if defined( __parallel )
    offset = 0

    ALLOCATE( get_buffer(SUM( transfer_index(4,:) )) )
!
!-- Copy from input array (data) to RMA window to allow the target PEs to get the appropiate data.
!-- At this point, a dummy surface element is added. This makes sure that every x-y gridbox owns
!-- at least one surface element. This way, bookkeeping becomes much easier.
    lo_index = nr_surfaces_in_tb
    DO  i = nxl, nxr
       DO  j = nys, nyn
          is = lo_index(s_index_in_window(j,i)) + 1
          ie = is + m_end_index(j,i) - m_start_index(j,i)
!
!--       Store number of surface elements in dummy additional surface element
          array_1d(is-1)  = e_end_index(j,i) - e_start_index(j,i) + 1
          array_1d(is:ie) = data(m_start_index(j,i):m_end_index(j,i))
          lo_index(s_index_in_window(j,i)) = lo_index(s_index_in_window(j,i)) +                    &
                                             e_end_index(j,i) - e_start_index(j,i) + 1
       ENDDO
    ENDDO
!
!-- On target PE, get data from source PEs which are assigned for output on this PE.
    CALL MPI_WIN_FENCE( 0, win_surf, ierr )

    buf_start = 1
    DO  n = 0, numprocs-1
       n_trans = transfer_index(4,n)
       IF ( n_trans > 0 )  THEN
          disp = transfer_index(3,n) - 1
          CALL MPI_GET( get_buffer(buf_start), n_trans, MPI_REAL, n, disp, n_trans, MPI_REAL,      &
                        win_surf, ierr )
          buf_start = buf_start + n_trans
       ENDIF
    ENDDO

    CALL MPI_WIN_FENCE( 0, win_surf, ierr )
!
!-- Copy data to output buffer. Here, the outpuf buffer matches the indices global_start and
!-- global_end, reduced by array_out_reduction.
!-- The use of array_out_reduction allows to use 32 bit indices in sm_allocate_shared.
    ind_gb  = 1
    DO  i = 1, SIZE( local_indices, 2 )
       ind_out = local_indices(1,i) - array_out_reduction
       DO  j = 1, local_indices(2,i)
          array_out(ind_out) = get_buffer(ind_gb)
          ind_out = ind_out+1
          ind_gb  = ind_gb+1
       ENDDO
    ENDDO

    DEALLOCATE( get_buffer )

!
!-- Write data to disk.
    CALL sm_io%sm_node_barrier()  ! has no effect if I/O on limited number of cores is inactive
    IF ( sm_io%iam_io_pe )  THEN
       IF ( less_than_2g_elements )  THEN
          CALL MPI_FILE_SET_VIEW( fh, array_position, MPI_REAL, ft_surf, 'native', MPI_INFO_NULL,  &
                                  ierr )
       ELSE
          offset =  array_position + local_offset
          CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )

       ENDIF

       CALL MPI_FILE_WRITE_ALL( fh, array_out, SIZE( array_out ), MPI_REAL, status, ierr )
    ENDIF
    CALL sm_io%sm_node_barrier()
#else
    DO  i = nxl, nxr
       DO  j = nys, nyn
          array_out(e_start_index(j,i)) = e_end_index(j,i) - e_start_index(j,i) + 1
          array_out(e_start_index(j,i)+1:e_end_index(j,i)) =                                       &
                                                          data(m_start_index(j,i):m_end_index(j,i))
       ENDDO
    ENDDO

    CALL posix_lseek( fh, array_position )
    CALL posix_write( fh, array_out, SIZE(array_out) )
#endif

    array_position = array_position + total_number_of_surface_elements * wp

 END SUBROUTINE wrd_mpi_io_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read 2d-REAL surface data array with MPI-IO.
!> This consist of multiple 1d-REAL surface data arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wrd_mpi_io_surface_2d( name, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)         ::  name  !<

    INTEGER(iwp)                         ::  i     !<

    REAL(wp), INTENT(IN), DIMENSION(:,:) ::  data  !<


    DO  i = 1, SIZE( data, 1 )
       CALL wrd_mpi_io_surface( name, data(i,:), first_index = i )
    ENDDO

 END SUBROUTINE wrd_mpi_io_surface_2d




!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Close restart file for MPI-IO
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_close

    IMPLICIT NONE

    INTEGER(iwp)                       ::  gh_size  !<
    INTEGER(KIND=rd_offset_kind)       ::  offset   !<
#if defined( __parallel )
    INTEGER, DIMENSION(rd_status_size) ::  status   !<
#endif

#if ! defined( __parallel )
    TYPE(C_PTR)                        ::  buf_ptr  !<
#endif


    offset = 0

    IF ( wr_flag  .AND.  sm_io%iam_io_pe )  THEN

       tgh%nr_int      = header_int_index - 1
       tgh%nr_char     = header_char_index - 1
       tgh%nr_real     = header_real_index - 1
       tgh%nr_arrays   = header_array_index - 1
       tgh%total_nx    = iog%nx + 1
       tgh%total_ny    = iog%ny + 1
       tgh%pes_along_x = npex
       tgh%pes_along_y = npey
       IF ( total_domain_boundaries_on_file )  THEN   ! Not sure, if LOGICAL interpretation is the same for all compilers,
          tgh%outer_ghost_layers_included = 1         ! therefore store as INTEGER in general header
       ELSE
          tgh%outer_ghost_layers_included = 0
       ENDIF
!
!--    Check for big/little endian format. This check is currently not used, and could be removed
!--    if we can assume little endian as the default on all machines.
       CALL rd_mpi_io_check_endian( tgh%endian )

!
!--    Set default view
#if defined( __parallel )
       CALL MPI_FILE_SET_VIEW( fh, offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr )
#endif
!
!--    Write header file
       gh_size = storage_size(tgh) / 8
       IF ( myid == 0 )  THEN   ! myid = 0 always performs I/O, even if I/O is limited to some cores
#if defined( __parallel )
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, tgh, gh_size, MPI_INT, status, ierr )
          header_position = header_position + gh_size
!
!--       INTEGER values
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, int_names, SIZE( int_names ) * 32, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( int_names ) * 32

          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, int_values, SIZE( int_values ), MPI_INT, status, ierr )
          header_position = header_position + SIZE( int_values ) * iwp
!
!--       Character entries
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, text_lines, SIZE( text_lines ) * 128, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( text_lines ) * 128
!
!---      REAL values
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, real_names, SIZE( real_names ) * 32, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( real_names ) * 32

          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, real_values, SIZE( real_values ), MPI_REAL, status, ierr )
          header_position = header_position + SIZE( real_values ) * wp
!
!--       2d- and 3d- distributed array headers, all replicated array headers
          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, array_names, SIZE( array_names ) * 32, MPI_CHAR, status, ierr )
          header_position = header_position + SIZE( array_names ) * 32

          CALL MPI_FILE_SEEK( fh, header_position, MPI_SEEK_SET, ierr )
          CALL MPI_FILE_WRITE( fh, array_offset, SIZE( array_offset ) * MPI_OFFSET_KIND, MPI_BYTE, &
                               status, ierr )
          header_position = header_position + SIZE( array_offset ) * rd_offset_kind
#else
          CALL posix_lseek( fh, header_position )
          buf_ptr = C_LOC( tgh )
          CALL posix_write( fh, buf_ptr, gh_size )
          header_position = header_position + gh_size
!
!--       INTEGER values
          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, int_names )
          header_position = header_position + SIZE( int_names ) * 32

          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, int_values, SIZE( int_values ) )
          header_position = header_position + SIZE( int_values ) * iwp
!
!--       Character entries
          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, text_lines )
          header_position = header_position + SIZE( text_lines ) * 128
!
!--       REAL values
          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, real_names )
          header_position = header_position + SIZE( real_names ) * 32

          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, real_values, SIZE( real_values ) )
          header_position = header_position + SIZE( real_values ) * wp
!
!--       2d- and 3d-distributed array headers, all replicated array headers
          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, array_names )
          header_position = header_position + SIZE( array_names ) * 32

          CALL posix_lseek( fh, header_position )
          CALL posix_write( fh, array_offset, SIZE( array_offset ) )
          header_position = header_position + SIZE( array_offset ) * rd_offset_kind
#endif
          IF ( debug_output )  CALL rd_mpi_io_print_header
       ENDIF

    ENDIF

!
!-- Free file types
    CALL rd_mpi_io_free_filetypes

!
!-- Close MPI-IO files
#if defined( __parallel )
!
!-- Free RMA windows
    IF ( cyclic_fill_mode )  THEN
       CALL MPI_WIN_FREE( rmawin_2di, ierr )
       CALL MPI_WIN_FREE( rmawin_2di8, ierr )
       CALL MPI_WIN_FREE( rmawin_2d, ierr )
       CALL MPI_WIN_FREE( rmawin_3d, ierr )
    ENDIF
#endif

    IF ( ALLOCATED( e_start_index )     )   DEALLOCATE( e_start_index )
    IF ( ALLOCATED( e_end_index )       )   DEALLOCATE( e_end_index )
    IF ( ALLOCATED( m_start_index )     )   DEALLOCATE( m_start_index )
    IF ( ALLOCATED( m_end_index )       )   DEALLOCATE( m_end_index )
    IF ( ALLOCATED( nr_surfaces_in_tb ) )   DEALLOCATE( nr_surfaces_in_tb )
    IF ( ALLOCATED( s_index_in_tb )     )   DEALLOCATE( s_index_in_tb )
    IF ( ALLOCATED( s_index_in_window ) )   DEALLOCATE( s_index_in_window )
    IF ( ALLOCATED( transfer_index )    )   DEALLOCATE( transfer_index )

    IF ( .NOT. pe_active_for_read )  RETURN
!
!-- TODO: better explain the following message
!-- In case on non cyclic read, pe_active_for_read is set .TRUE.
    IF ( sm_io%iam_io_pe )  THEN

#if defined( __parallel )
       CALL MPI_FILE_CLOSE( fh, ierr )
#else
       CALL posix_close( fh )
#endif

    ENDIF

    fh = -1

    restart_file_size = array_position / ( 1024.0_dp * 1024.0_dp )

 END SUBROUTINE rd_mpi_io_close



 FUNCTION rd_mpi_io_check_open()  RESULT ( isopen )

    LOGICAL ::  isopen

    isopen = ( fh /= -1 )

 END FUNCTION rd_mpi_io_check_open



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine prepares a filetype and some variables for the I/O of surface data.
!> Whenever a new set of start_index and end_index is used, rd_mpi_io_surface_filetypes has to be
!> called. A main feature of this subroutine is computing the global start indices of the 1d- and
!> 2d- surface arrays.
!> A reorganization of data is done to get nearly equally distributed number of cells on all PEs.
!> The total number of surface cells can be greater than 2**31-1 so an INTEGER(KIND=8) index space
!> is required for global indices.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write,          &
                                                   global_start, global_end )

    IMPLICIT NONE

    INTEGER(idp) ::  e_lo_start        !<
    INTEGER(iwp) ::  i                 !<  loop index
    INTEGER(iwp) ::  j                 !<  loop index
    INTEGER(iwp) ::  index_offset      !<
    INTEGER(iwp) ::  last_end_index    !<
    INTEGER(idp) ::  lo_start          !<
    INTEGER(idp) ::  nr_surfcells_pe   !<
    INTEGER(idp) ::  rest_cells_pe     !<
    INTEGER(idp) ::  rest_bound        !<
#if defined( __parallel )
    INTEGER(iwp) ::  io_end_index      !<
    INTEGER(iwp) ::  io_start_index    !<
    INTEGER(iwp) ::  n                 !<  loop index
    INTEGER(idp) ::  nr_previous       !<
    INTEGER(iwp) ::  sum_local_ind     !<  sum of local indices (I4)
#endif

    INTEGER(idp), DIMENSION(0:numprocs-1,2) ::  nr_surfcells_all_r  !<
    INTEGER(idp), DIMENSION(0:numprocs-1,2) ::  nr_surfcells_all_s  !<
#if defined( __parallel )
    INTEGER(iwp), DIMENSION(1)              ::  dims1               !< global dimension for MPI_TYPE_CREATE_SUBARRAY
    INTEGER(idp), DIMENSION(1)              ::  dims8               !< global dimension for MPI_TYPE_CREATE_SUBARRAY   (I8)
    INTEGER(iwp), DIMENSION(1)              ::  lsize1              !< local size for MPI_TYPE_CREATE_SUBARRAY
    INTEGER(iwp), DIMENSION(0:numprocs-1)   ::  nr_gp_with_surfaces_for_pe  !< number of horizontal gridpoints containing surface elements for the respective pe
    INTEGER(iwp), DIMENSION(0:numprocs-1)   ::  nr_surface_elements_for_pe !< total number of surface elements for the respective pe
    INTEGER(idp), DIMENSION(0:npex)         ::  nr_surf_cells_x     !<
    INTEGER(idp), DIMENSION(0:npex)         ::  nr_surf_cells_x_s   !<
    INTEGER(iwp), DIMENSION(1)              ::  start1              !< start index for MPI_TYPE_CREATE_SUBARRAY
    INTEGER(idp), DIMENSION(1)              ::  start8              !< start index for MPI_TYPE_CREATE_SUBARRAY   (I8)

    INTEGER(idp), DIMENSION(nxl:nxr)        ::  sum_y               !<
#endif

    INTEGER(iwp), INTENT(INOUT), DIMENSION(nys:nyn,nxl:nxr) ::  end_index     !< local end indx
    INTEGER(idp), INTENT(INOUT), DIMENSION(nys:nyn,nxl:nxr) ::  global_start  !< global start index
    INTEGER(idp), INTENT(INOUT), DIMENSION(nys:nyn,nxl:nxr) ::  global_end    !< global end index
    INTEGER(iwp), INTENT(INOUT), DIMENSION(nys:nyn,nxl:nxr) ::  start_index   !< local start index
#if defined( __parallel )
    INTEGER(idp), DIMENSION(0:myidy,nxl:nxr) ::  nr_previous_y      !<
    INTEGER(idp), DIMENSION(0:npey,nxl:nxr)  ::  nr_surf_cells_y    !<
    INTEGER(idp), DIMENSION(0:npey,nxl:nxr)  ::  nr_surf_cells_y_s  !<
    INTEGER(idp), DIMENSION(4,0:numprocs-1)  ::  transfer_index_s   !<
#endif

    LOGICAL, INTENT(OUT) ::  data_to_write        !< returns .TRUE., if surface data have been written
    LOGICAL              ::  only_dummy_elements  !< only dummy elements, i.e. no data to write

    data_to_write = .FALSE.
!
!-- Clean up previous calls.
#if defined( __parallel )
    IF ( win_surf /= -1 )  THEN
       CALL MPI_WIN_FREE( win_surf, ierr )
       DEALLOCATE( array_1d )
       win_surf = -1
    ENDIF
    IF ( ft_surf /= -1  .AND.  sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_FREE( ft_surf, ierr )
    ENDIF
    ft_surf = -1
    IF ( sm_io%is_sm_active() )  THEN
       IF ( win_out /= -1 )  THEN
          CALL sm_io%sm_free_shared( win_out )
          win_out = -1
       ENDIF
    ELSE
       IF ( ASSOCIATED( array_out ) )  DEALLOCATE( array_out )
    ENDIF
#else
    IF ( ASSOCIATED( array_out ) )  DEALLOCATE( array_out )
#endif

    IF ( cyclic_fill_mode )  THEN
       CALL cyclic_fill_surface_filetype
       RETURN
    ELSE
       IF ( .NOT. ALLOCATED( e_end_index )        )  ALLOCATE( e_end_index(nys:nyn,nxl:nxr) )
       IF ( .NOT. ALLOCATED( e_start_index )      )  ALLOCATE( e_start_index(nys:nyn,nxl:nxr) )
       IF ( .NOT. ALLOCATED( m_end_index )        )  ALLOCATE( m_end_index(nys:nyn,nxl:nxr) )
       IF ( .NOT. ALLOCATED( m_start_index )      )  ALLOCATE( m_start_index(nys:nyn,nxl:nxr) )
       IF ( .NOT. ALLOCATED( nr_surfaces_in_tb )  )  ALLOCATE( nr_surfaces_in_tb(0:numprocs-1) )
       IF ( .NOT. ALLOCATED( s_index_in_tb )      )  ALLOCATE( s_index_in_tb(0:numprocs-1) )
       IF ( .NOT. ALLOCATED( s_index_in_window )  )  ALLOCATE( s_index_in_window(nys:nyn,nxl:nxr) )
       IF ( .NOT. ALLOCATED( transfer_index )     )  ALLOCATE( transfer_index(4,0:numprocs-1) )
    ENDIF

    IF ( wr_flag)  THEN
!
!--    Add one dummy value at every grid box.
!--    This allows to use MPI_FILE_WRITE_ALL and MPI_FILE_READ_ALL with subarray file type.
       index_offset   = 0
       last_end_index = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             e_start_index(j,i) = start_index (j,i) + index_offset
             IF ( end_index (j,i) - start_index(j,i) < 0 )  THEN
                e_end_index (j,i) = last_end_index+1
                last_end_index    = last_end_index+1
             ELSE
                e_end_index (j,i) = end_index(j,i) + index_offset + 1
                last_end_index    = e_end_index (j,i)
             ENDIF
             index_offset = index_offset + 1
           ENDDO
       ENDDO
#if defined( __parallel )
!
!--    Compute indices for global, PE independent 1-d surface element array.
       nr_surf_cells_y_s = 0
!
!--    Count number of surface elements in y-direction for every x.
       DO  i = nxl, nxr
          nr_surf_cells_y_s(myidy,i) = SUM( e_end_index (:,i) - e_start_index (:,i) + 1 )
       ENDDO
!
!--    Distribute these elements to all PEs along y.
       CALL MPI_ALLREDUCE( nr_surf_cells_y_s, nr_surf_cells_y, SIZE( nr_surf_cells_y ),            &
                           MPI_INTEGER8, MPI_SUM, comm1dy, ierr )

!
!--    Sum all surface elements along y for individual x PEs
       nr_surf_cells_x_s        = 0
       nr_surf_cells_x_s(myidx) = SUM( nr_surf_cells_y )
!
!--    Distribute to all PEs along x.
       CALL MPI_ALLREDUCE( nr_surf_cells_x_s, nr_surf_cells_x, SIZE( nr_surf_cells_x ),            &
                           MPI_INTEGER8, MPI_SUM, comm1dx, ierr )
       DO  i = nxl, nxr
          nr_previous_y(:,i) = 0
          DO  n = 1, myidy
             nr_previous_y(n,i) = nr_previous_y(n-1,i) + nr_surf_cells_y(n-1,i)
          ENDDO
       ENDDO

       sum_y(nxl) = SUM( nr_surf_cells_y(:,nxl) )
       DO  i = nxl, nxr
          IF ( i > nxl )  THEN
             sum_y(i) = sum_y(i-1) + SUM( nr_surf_cells_y(:,i) )
          ENDIF
       ENDDO

       nr_previous = 0
       IF ( myidx >= 1 )  THEN
          nr_previous = SUM(nr_surf_cells_x(0:myidx-1))
       ENDIF

       global_start(nys,nxl) = 1 + nr_previous + nr_previous_y(myidy,nxl)
       DO  j = nys+1, nyn
          global_start(j,nxl) = global_start(j-1,nxl) +                                            &
                                e_end_index(j-1,nxl) - e_start_index(j-1,nxl) + 1
       ENDDO

       DO  i = nxl+1, nxr
          global_start(nys,i) = 1 + nr_previous + nr_previous_y(myidy,i) + sum_y(i-1)
          DO  j = nys+1, nyn
             global_start(j,i) = global_start(j-1,i) +                                             &
                                 e_end_index(j-1,i) - e_start_index(j-1,i) + 1
          ENDDO
       ENDDO
#else
       global_start = e_start_index
#endif
       DO  i = nxl, nxr
          DO  j = nys, nyn
             global_end(j,i) = global_start(j,i) + e_end_index (j,i) - e_start_index (j,i)
          ENDDO
       ENDDO

    ELSE
!
!--    In case of read, compute e_start_index and e_end_index for current processor grid.
!--    This data contains one extra value for every i and j.
       e_lo_start = 1
       lo_start   = 1
       DO  i = nxl, nxr
          DO  j = nys, nyn
             e_start_index(j,i) = e_lo_start
             e_end_index(j,i)   = e_lo_start + global_end(j,i) - global_start(j,i)
             e_lo_start         = e_lo_start + global_end(j,i) - global_start(j,i) + 1
             start_index(j,i)   = lo_start
             end_index(j,i)     = lo_start + global_end(j,i) - global_start(j,i) - 1
             lo_start           = lo_start + global_end(j,i) - global_start(j,i)
          ENDDO
       ENDDO

    ENDIF

    nr_surfcells_all_s = 0
    nr_surfcells_all_s(myid,1) = MAXVAL( e_end_index ) ! don't split surface elements of one gridbox
    nr_surfcells_all_s(myid,2) = MAXVAL( e_end_index - e_start_index )

#if defined( __parallel )
    CALL MPI_ALLREDUCE( nr_surfcells_all_s, nr_surfcells_all_r, SIZE( nr_surfcells_all_s ),        &
                        MPI_INTEGER8, MPI_SUM, comm2d, ierr )
#else
    nr_surfcells_all_r = nr_surfcells_all_s
#endif

    total_number_of_surface_elements = 0
    DO  i = 0, numprocs-1
       IF ( i == myid )  THEN
          glo_start = total_number_of_surface_elements + 1
       ENDIF
       total_number_of_surface_elements = total_number_of_surface_elements + nr_surfcells_all_r(i,1)
    ENDDO
    only_dummy_elements   = ( MAXVAL( nr_surfcells_all_r(:,2) ) <= 0 )
    less_than_2g_elements = ( total_number_of_surface_elements < INT( HUGE( i4_rep ), 8 ) )

!
!-- Compute indices of equally distributed surface elements.
!-- Number of surface elements scheduled for ouput on this PE:
    nr_surfcells_pe  = total_number_of_surface_elements  / INT( numprocs, idp )
    rest_cells_pe    = MOD( total_number_of_surface_elements, INT( numprocs, idp ) )
    rest_bound       = rest_cells_pe * ( nr_surfcells_pe + 1 )
    m_start_index    = start_index
    m_end_index      = end_index

!
!-- Compute number of surface elements on source PE, which have to be send to the corresponding
!-- target PE.
#if defined( __parallel )
    nr_gp_with_surfaces_for_pe = 0
    nr_surface_elements_for_pe = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( rest_cells_pe == 0 )  THEN
             s_index_in_window(j,i) = ( global_start(j,i) - 1 ) / nr_surfcells_pe
          ELSE
             IF ( global_start(j,i) <= rest_bound )  THEN
                s_index_in_window(j,i) = ( global_start(j,i) - 1 ) / ( nr_surfcells_pe + 1 )
             ELSE
                s_index_in_window(j,i) = ( global_start(j,i) - rest_bound - 1 ) / nr_surfcells_pe
                s_index_in_window(j,i) = s_index_in_window(j,i) + rest_cells_pe
             ENDIF
          ENDIF
          nr_gp_with_surfaces_for_pe(s_index_in_window(j,i)) =                                     &
                                             nr_gp_with_surfaces_for_pe(s_index_in_window(j,i)) + 1
          nr_surface_elements_for_pe(s_index_in_window(j,i)) =                                     &
                                             nr_surface_elements_for_pe(s_index_in_window(j,i)) +  &
                                             e_end_index(j,i) - e_start_index(j,i) + 1
       ENDDO
    ENDDO

!
!-- Compute start index in the transfer buffer on the source side for the corresponding target PE.
    s_index_in_tb(0)     = 1
    nr_surfaces_in_tb(0) = 1
    DO  n = 1, numprocs-1
       s_index_in_tb(n)     = s_index_in_tb(n-1)     + nr_gp_with_surfaces_for_pe(n-1)
       nr_surfaces_in_tb(n) = nr_surfaces_in_tb(n-1) + nr_surface_elements_for_pe(n-1)
    ENDDO

!
!-- Buffer distribution on the source side.
    DO  n = 0, numprocs-1
       transfer_index_s(1,n) = s_index_in_tb(n)
       transfer_index_s(2,n) = nr_gp_with_surfaces_for_pe(n)
       transfer_index_s(3,n) = nr_surfaces_in_tb(n)
       transfer_index_s(4,n) = nr_surface_elements_for_pe(n)
    ENDDO

    CALL MPI_ALLTOALL( transfer_index_s, 4, MPI_INTEGER8, transfer_index, 4, MPI_INTEGER8, comm2d, &
                       ierr)
!
!-- Buffer distribution on the target side side.
    CALL get_remote_indices()
!
!-- Create surface element file type.
    IF ( total_number_of_surface_elements > 0 .AND. .NOT. only_dummy_elements)  THEN
        data_to_write = .TRUE.
    ELSE
        data_to_write = .FALSE.
    ENDIF

    start8(1) = MINVAL( local_indices(1,:) ) - 1
    IF ( sm_io%is_sm_active() )  THEN
       sum_local_ind = SUM( local_indices(2,:) )
       CALL MPI_ALLREDUCE( sum_local_ind, lsize1(1), 1, MPI_INTEGER, MPI_SUM, sm_io%comm_shared,   &
                           ierr )
    ELSE
       lsize1(1) = SUM( local_indices(2,:) )
    ENDIF

    IF ( less_than_2g_elements )  THEN
       CALL MPI_ALLREDUCE( global_end(nyn,nxr), dims8(1), 1, MPI_INTEGER8, MPI_MAX, comm2d, ierr )
       dims1(1)  = dims8(1)
       start1(1) = start8(1)

       IF ( sm_io%iam_io_pe )  THEN
          IF ( total_number_of_surface_elements > 0 )  THEN
             CALL MPI_TYPE_CREATE_SUBARRAY( 1, dims1, lsize1, start1, MPI_ORDER_FORTRAN, MPI_REAL, &
                                            ft_surf, ierr )
             CALL MPI_TYPE_COMMIT( ft_surf, ierr )
          ENDIF
       ENDIF
    ELSE
       local_offset =  start8(1) * wp
    ENDIF

!
!-- Allocate rma window to supply surface data to other PEs.
    CALL rd_alloc_rma_mem( array_1d, SUM( nr_surface_elements_for_pe ), win_surf )

!
!-- Allocate shared array on IO-PE to supply data for MPI-IO (write or read).
    array_out_reduction = start8(1)
    IF ( sm_io%is_sm_active() )  THEN
       IF ( sm_io%iam_io_pe )  THEN
          io_start_index = start8(1) + 1_dp - array_out_reduction
          io_end_index   = start8(1) + INT( lsize1(1), 8 ) - array_out_reduction
       ENDIF
       CALL MPI_BCAST( io_start_index,        1, MPI_INTEGER,  0, sm_io%comm_shared, ierr )
       CALL MPI_BCAST( io_end_index,          1, MPI_INTEGER,  0, sm_io%comm_shared, ierr )
       CALL MPI_BCAST( array_out_reduction,   1, MPI_INTEGER8, 0, sm_io%comm_shared, ierr )
       CALL sm_io%sm_allocate_shared( array_out, io_start_index, io_end_index, win_out )
    ELSE
       io_start_index = start8(1) + 1_dp - array_out_reduction
       io_end_index   = start8(1) + INT( lsize1(1) , 8 ) - array_out_reduction
       ALLOCATE( array_out(io_start_index:io_end_index) )
    ENDIF
#else
    IF ( total_number_of_surface_elements > 0  .AND.  .NOT. only_dummy_elements )  THEN
        data_to_write = .TRUE.
    ELSE
        data_to_write = .FALSE.
    ENDIF
    ALLOCATE( array_out(1:total_number_of_surface_elements) )
#endif

 CONTAINS

    SUBROUTINE cyclic_fill_surface_filetype

       INTEGER(iwp) ::  i  !<  loop index
       INTEGER(iwp) ::  j  !<  loop index


       IF ( .NOT. ALLOCATED( o_start_index ) )  ALLOCATE( o_start_index(nys:nyn,nxl:nxr) )
       IF ( .NOT. ALLOCATED( o_end_index )   )  ALLOCATE( o_end_index(nys:nyn,nxl:nxr)   )

       lo_start   = 1
       DO  i = nxl, nxr
           DO  j = nys, nyn
               o_start_index(j,i) = lo_start
               o_end_index(j,i)   = lo_start
               lo_start           = lo_start + 1
           ENDDO
       ENDDO
       start_index = o_start_index
       end_index   = o_end_index

!
!--    Activate grid of the smaller prerun, i.e. grid variables like nxl, nxr, nys, nyn and others
!--    are set according to the prerun layout.
       CALL prerun_grid%activate_grid_from_this_class()

       IF ( pe_active_for_read )  THEN

          IF ( .NOT. ALLOCATED( c_global_start ) )  ALLOCATE( c_global_start(nys:nyn,nxl:nxr) )
          IF ( .NOT. ALLOCATED( c_global_end )   )  ALLOCATE( c_global_end(nys:nyn,nxl:nxr)   )
          IF ( .NOT. ALLOCATED( c_start_index )  )  ALLOCATE( c_start_index(nys:nyn,nxl:nxr)  )
          IF ( .NOT. ALLOCATED( c_end_index )    )  ALLOCATE( c_end_index(nys:nyn,nxl:nxr)    )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                c_global_start(j,i) = global_start(j,i)
                c_global_end(j,i)   = global_end(j,i)
             ENDDO
          ENDDO

!
!--       Recursive call of rd_mpi_io_surface_filetypes.
!--       Prerun data are read, but they are treated as if they are mainrun data, just on a smaller
!--       grid.
          cyclic_fill_mode = .FALSE.
          CALL rd_mpi_io_surface_filetypes( c_start_index, c_end_index, data_to_write,             &
                                            c_global_start, c_global_end )
          cyclic_fill_mode = .TRUE.

       ENDIF
!
!--    Activate grid of the mainrun, i.e. grid variables like nxl, nxr, nys, nyn and others
!--    are set according to the mainrun layout.
       CALL mainrun_grid%activate_grid_from_this_class()

#if defined( __parallel )
       CALL MPI_BCAST( data_to_write, 1, MPI_LOGICAL, 0, comm2d, ierr )
#endif

    END SUBROUTINE cyclic_fill_surface_filetype

#if defined( __parallel )
!
!-- Get the indices of the surface elements inside the RMA window on the remote PE.
!-- This information is required to fetch the surface element data on remote PEs
!-- in rrd_mpi_io_surface and wrd_mpi_io_surface.
    SUBROUTINE get_remote_indices

       INTEGER(idp) ::  buf_start  !<
       INTEGER(iwp) ::  i          !<
       INTEGER(iwp) ::  j          !<
       INTEGER(iwp) ::  n          !<
       INTEGER(iwp) ::  n_trans    !<
       INTEGER(iwp) ::  win_ind    !<

       INTEGER(KIND=MPI_ADDRESS_KIND) ::  disp     !< displacement in RMA window
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window

       INTEGER(idp), DIMENSION(0:numprocs-1) ::  lo_index         !<

       INTEGER(idp), POINTER, DIMENSION(:,:) ::  surf_val_index  !<


       IF ( ALLOCATED( local_indices ) )  DEALLOCATE( local_indices )
       ALLOCATE( local_indices(2,MAX( SUM( transfer_index(2,:) ), 2_idp )))

       local_indices(1,:) = HUGE(local_indices)
       local_indices(2,:) = 0

       winsize = MAX( 2 * SUM( nr_gp_with_surfaces_for_pe ), 2 )

       ALLOCATE( surf_val_index(2,winsize) )
       winsize = winsize * idp
       CALL MPI_WIN_CREATE( surf_val_index, winsize, idp, MPI_INFO_NULL, comm2d, win_ind, ierr )
       CALL MPI_WIN_FENCE( 0, win_ind, ierr )

       lo_index = s_index_in_tb
       DO  i = nxl, nxr
          DO  j = nys, nyn
             surf_val_index(1,lo_index(s_index_in_window(j,i))) = global_start(j,i)
             surf_val_index(2,lo_index(s_index_in_window(j,i))) = global_end(j,i) -                &
                                                                  global_start(j,i) + 1
             lo_index(s_index_in_window(j,i)) = lo_index(s_index_in_window(j,i)) + 1
          ENDDO
       ENDDO

       CALL MPI_WIN_FENCE( 0, win_ind, ierr )

       buf_start = 1
       DO  n = 0, numprocs-1
          n_trans = transfer_index(2,n)
          IF ( n_trans > 0 )  THEN
             disp = 2 * ( transfer_index(1,n) - 1 )
             CALL MPI_GET( local_indices(1,buf_start), 2*n_trans, MPI_INTEGER8, n, disp, 2*n_trans, &
                           MPI_INTEGER8, win_ind, ierr )
             buf_start = buf_start + n_trans
          ENDIF
       ENDDO

       CALL MPI_WIN_FENCE( 0, win_ind, ierr )

       CALL MPI_WIN_FREE( win_ind, ierr )

       DEALLOCATE( surf_val_index )

    END SUBROUTINE get_remote_indices

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate memory and create window for one-sided communication (1-d INTEGER array)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE rd_alloc_rma_mem( array, idim, win )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)       ::  idim     !< Dimension of this 1-D array
       INTEGER                        ::  ierr     !< MPI error code
       INTEGER(iwp), INTENT(OUT)      ::  win      !< MPI window
       INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window

       REAL(wp), DIMENSION(:), POINTER, INTENT(INOUT) ::  array  !< array to access RMA window locally


       winsize = MAX( idim, 2 )
       ALLOCATE( array(winsize) )
       winsize = winsize * wp
       CALL MPI_WIN_CREATE( array, winsize, wp, MPI_INFO_NULL, comm2d, win, ierr )
       array = -1
       CALL MPI_WIN_FENCE( 0, win, ierr )

    END SUBROUTINE rd_alloc_rma_mem
#endif

 END SUBROUTINE rd_mpi_io_surface_filetypes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates file types to access 2d-/3d-REAL arrays and 2d-INTEGER arrays
!> distributed in blocks among processes to a single file that contains the global arrays.
!--------------------------------------------------------------------------------------------------!
  SUBROUTINE rd_mpi_io_create_filetypes

    IMPLICIT NONE

    INTEGER, DIMENSION(2) ::  dims2   !<
    INTEGER, DIMENSION(2) ::  lize2   !<
    INTEGER, DIMENSION(2) ::  start2  !<

    INTEGER, DIMENSION(3) ::  dims3   !<
    INTEGER, DIMENSION(3) ::  lize3   !<
    INTEGER, DIMENSION(3) ::  start3  !<

    TYPE(domain_decomposition_grid_features) ::  save_io_grid  !< temporary variable to store grid settings


    IF ( sm_io%is_sm_active() )  THEN
       save_io_grid = sm_io%io_grid
    ENDIF

    IF( .NOT. pe_active_for_read )  RETURN

    IF ( cyclic_fill_mode )  CALL prerun_grid%activate_grid_from_this_class()

!
!-- Decision, if storage with or without ghost layers.
!-- Please note that the indexing of the global array always starts at 0, even in Fortran.
!-- Therefore the PE local indices have to be shifted by nbgp in the case with ghost layers.
    IF ( total_domain_boundaries_on_file )  THEN

       iog%nxl = nxl + nbgp
       iog%nxr = nxr + nbgp
       iog%nys = nys + nbgp
       iog%nyn = nyn + nbgp
       iog%nnx = nnx
       iog%nny = nny
       iog%nx  = nx + 2 * nbgp
       iog%ny  = ny + 2 * nbgp
       IF ( myidx == 0 )  THEN
          iog%nxl = iog%nxl - nbgp
          iog%nnx = iog%nnx + nbgp
       ENDIF
       IF ( myidx == npex-1 )  THEN
          iog%nxr = iog%nxr + nbgp
          iog%nnx = iog%nnx + nbgp
       ENDIF
       IF ( myidy == 0 )  THEN
          iog%nys = iog%nys - nbgp
          iog%nny = iog%nny + nbgp
       ENDIF
       IF ( myidy == npey-1 )  THEN
          iog%nyn = iog%nyn + nbgp
          iog%nny = iog%nny + nbgp
       ENDIF

       CALL sm_io%sm_adjust_outer_boundary()

    ELSE

       iog%nxl = nxl
       iog%nxr = nxr
       iog%nys = nys
       iog%nyn = nyn
       iog%nnx = nnx
       iog%nny = nny
       iog%nx  = nx
       iog%ny  = ny

    ENDIF

    IF ( sm_io%is_sm_active() )  THEN
#if defined( __parallel )
       CALL sm_io%sm_allocate_shared( array_2d,  sm_io%io_grid%nxl, sm_io%io_grid%nxr,             &
                                      sm_io%io_grid%nys, sm_io%io_grid%nyn, win_2dr )
       CALL sm_io%sm_allocate_shared( array_2di, save_io_grid%nxl, save_io_grid%nxr,               &
                                      save_io_grid%nys, save_io_grid%nyn, win_2di )
       CALL sm_io%sm_allocate_shared( array_2di8, save_io_grid%nxl, save_io_grid%nxr,              &
                                      save_io_grid%nys, save_io_grid%nyn, win_2di8 )
       CALL sm_io%sm_allocate_shared( array_3d, nzb, nzt+1, sm_io%io_grid%nxl, sm_io%io_grid%nxr,  &
                                      sm_io%io_grid%nys, sm_io%io_grid%nyn, win_3dr )
       CALL sm_io%sm_allocate_shared( array_3di4, nzb, nzt+1, sm_io%io_grid%nxl, sm_io%io_grid%nxr,&
                                      sm_io%io_grid%nys, sm_io%io_grid%nyn, win_3di4 )
#endif
    ELSE
       ALLOCATE( array_2d(iog%nxl:iog%nxr,iog%nys:iog%nyn) )
       ALLOCATE( array_2di(nxl:nxr,nys:nyn) )
       ALLOCATE( array_2di8(nxl:nxr,nys:nyn) )
       ALLOCATE( array_3d(nzb:nzt+1,iog%nxl:iog%nxr,iog%nys:iog%nyn) )
       ALLOCATE( array_3di4(nzb:nzt+1,iog%nxl:iog%nxr,iog%nys:iog%nyn) )
       sm_io%io_grid = iog
    ENDIF

!
!-- Create filetype for 2d-REAL array with ghost layers around the total domain
    dims2(1)  = iog%nx + 1
    dims2(2)  = iog%ny + 1

    lize2(1)  = sm_io%io_grid%nnx
    lize2(2)  = sm_io%io_grid%nny

    start2(1) = sm_io%io_grid%nxl
    start2(2) = sm_io%io_grid%nys

#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_CREATE_SUBARRAY( 2, dims2, lize2, start2, MPI_ORDER_FORTRAN, MPI_REAL,        &
                                      ft_2d, ierr )
       CALL MPI_TYPE_COMMIT( ft_2d, ierr )
    ENDIF
#endif
!
!-- Create filetype for 2d-INTEGER array without ghost layers around the total domain
    dims2(1)  = nx + 1
    dims2(2)  = ny + 1

    IF ( sm_io%is_sm_active() )  THEN

       lize2(1)  = save_io_grid%nnx
       lize2(2)  = save_io_grid%nny

       start2(1) = save_io_grid%nxl
       start2(2) = save_io_grid%nys

    ELSE

       lize2(1)  = nnx
       lize2(2)  = nny

       start2(1) = nxl
       start2(2) = nys

    ENDIF

#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_CREATE_SUBARRAY( 2, dims2, lize2, start2, MPI_ORDER_FORTRAN, MPI_INTEGER,     &
                                      ft_2di_nb, ierr )
       CALL MPI_TYPE_COMMIT( ft_2di_nb, ierr )

       CALL MPI_TYPE_CREATE_SUBARRAY( 2, dims2, lize2, start2, MPI_ORDER_FORTRAN, MPI_INTEGER8,    &
                                      ft_2di8_nb, ierr )
       CALL MPI_TYPE_COMMIT( ft_2di8_nb, ierr )
    ENDIF
#endif
!
!-- Create filetype for 3d-REAL array
    dims3(1)  = nz + 2
    dims3(2)  = iog%nx + 1
    dims3(3)  = iog%ny + 1

    lize3(1)  = dims3(1)
    lize3(2)  = sm_io%io_grid%nnx
    lize3(3)  = sm_io%io_grid%nny

    start3(1) = nzb
    start3(2) = sm_io%io_grid%nxl
    start3(3) = sm_io%io_grid%nys

#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_CREATE_SUBARRAY( 3, dims3, lize3, start3, MPI_ORDER_FORTRAN, MPI_INTEGER,     &
                                      ft_3di4, ierr )
       CALL MPI_TYPE_COMMIT( ft_3di4, ierr )

       CALL MPI_TYPE_CREATE_SUBARRAY( 3, dims3, lize3, start3, MPI_ORDER_FORTRAN, MPI_REAL, ft_3d, &
                                      ierr )
       CALL MPI_TYPE_COMMIT( ft_3d, ierr )
    ENDIF
#endif

    IF ( cyclic_fill_mode )  CALL mainrun_grid%activate_grid_from_this_class()

 END SUBROUTINE rd_mpi_io_create_filetypes



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates file types to access 3d-INTEGER*4 arrays and 3d-INTEGER*8 arrays
!> distributed in blocks among processes to a single file that contains the global arrays.
!> These types are only used for particle data.
!--------------------------------------------------------------------------------------------------!
  SUBROUTINE rd_mpi_io_particle_filetypes

    IMPLICIT NONE

    INTEGER, DIMENSION(3) ::  dims3   !<
    INTEGER, DIMENSION(3) ::  lize3   !<
    INTEGER, DIMENSION(3) ::  start3  !<

    TYPE(domain_decomposition_grid_features) ::  save_io_grid  !< temporary variable to store grid settings

!
!-- MPI_INTEGER8 is not standard MPI, but is supported on most MPI distibutions
!-- If not suppported, a workaround could be enabled with the following preprocessor directive
!#if defined( __NO_INTEGER8)
!    CALL MPI_TYPE_CONTIGUOUS( 2, MPI_INTEGER, MPI_INTEGER8, ierr )
!    CALL MPI_TYPE_COMMIT( MPI_INTEGER8, ierr )
!#endif

    IF ( sm_io%is_sm_active() )  THEN
       save_io_grid = sm_io%io_grid
    ENDIF

    IF( .NOT. pe_active_for_read )  RETURN

!
!-- Decision, if storage with or without ghost layers.
!-- Please note that the indexing of the global array always starts at 0, even in Fortran.
!-- Therefore the PE local indices have to be shifted by nbgp in the case with ghost layers.
    IF ( total_domain_boundaries_on_file )  THEN

       iog%nxl = nxl + nbgp
       iog%nxr = nxr + nbgp
       iog%nys = nys + nbgp
       iog%nyn = nyn + nbgp
       iog%nnx = nnx
       iog%nny = nny
       iog%nx  = nx + 2 * nbgp
       iog%ny  = ny + 2 * nbgp
       IF ( myidx == 0 )  THEN
          iog%nxl = iog%nxl - nbgp
          iog%nnx = iog%nnx + nbgp
       ENDIF
       IF ( myidx == npex-1 )  THEN
          iog%nxr = iog%nxr + nbgp
          iog%nnx = iog%nnx + nbgp
       ENDIF
       IF ( myidy == 0 )  THEN
          iog%nys = iog%nys - nbgp
          iog%nny = iog%nny + nbgp
       ENDIF
       IF ( myidy == npey-1 )  THEN
          iog%nyn = iog%nyn + nbgp
          iog%nny = iog%nny + nbgp
       ENDIF

       CALL sm_io%sm_adjust_outer_boundary()

    ELSE

       iog%nxl = nxl
       iog%nxr = nxr
       iog%nys = nys
       iog%nyn = nyn
       iog%nnx = nnx
       iog%nny = nny
       iog%nx  = nx
       iog%ny  = ny

    ENDIF

    IF ( sm_io%is_sm_active() )  THEN
#if defined( __parallel )
       CALL sm_io%sm_allocate_shared( array_3di8, nzb, nzt+1, sm_io%io_grid%nxl, sm_io%io_grid%nxr,&
                                      sm_io%io_grid%nys, sm_io%io_grid%nyn, win_3di8 )
#endif
    ELSE
       ALLOCATE( array_3di8(nzb:nzt+1,iog%nxl:iog%nxr,iog%nys:iog%nyn) )

       sm_io%io_grid = iog
    ENDIF

!
!-- Create filetype for 3d INTEGER array
    dims3(1)  = nz + 2
    dims3(2)  = iog%nx + 1
    dims3(3)  = iog%ny + 1

    lize3(1)  = dims3(1)
    lize3(2)  = sm_io%io_grid%nnx
    lize3(3)  = sm_io%io_grid%nny

    start3(1) = nzb
    start3(2) = sm_io%io_grid%nxl
    start3(3) = sm_io%io_grid%nys

#if defined( __parallel )
    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_CREATE_SUBARRAY( 3, dims3, lize3, start3, MPI_ORDER_FORTRAN, MPI_INTEGER8,    &
                                      ft_3di8, ierr )
       CALL MPI_TYPE_COMMIT( ft_3di8, ierr )
    ENDIF
#endif

 END SUBROUTINE rd_mpi_io_particle_filetypes



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates file types to access 3d-soil arrays distributed in blocks among processes
!> to a single file that contains the global arrays. It is not required for the serial mode.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_create_filetypes_3dsoil( nzb_soil, nzt_soil )

    IMPLICIT NONE

    INTEGER, INTENT(IN)   ::  nzb_soil  !<
    INTEGER, INTENT(IN)   ::  nzt_soil  !<

#if defined( __parallel )
    INTEGER, DIMENSION(3) ::  dims3     !<
    INTEGER, DIMENSION(3) ::  lize3     !<
    INTEGER, DIMENSION(3) ::  start3    !<


    IF ( sm_io%is_sm_active() )  THEN
       CALL sm_io%sm_allocate_shared( array_3d_soil, nzb_soil, nzt_soil, sm_io%io_grid%nxl,        &
                                      sm_io%io_grid%nxr, sm_io%io_grid%nys, sm_io%io_grid%nyn,     &
                                      win_3ds )
    ELSE
       ALLOCATE( array_3d_soil(nzb_soil:nzt_soil,iog%nxl:iog%nxr,iog%nys:iog%nyn) )
       sm_io%io_grid = iog
    ENDIF

!
!-- Create filetype for 3d-soil array
    dims3(1)  = nzt_soil - nzb_soil + 1
    dims3(2)  = iog%nx + 1
    dims3(3)  = iog%ny + 1

    lize3(1)  = dims3(1)
    lize3(2)  = sm_io%io_grid%nnx
    lize3(3)  = sm_io%io_grid%nny

    start3(1) = nzb_soil
    start3(2) = sm_io%io_grid%nxl
    start3(3) = sm_io%io_grid%nys

    IF ( sm_io%iam_io_pe )  THEN
       CALL MPI_TYPE_CREATE_SUBARRAY( 3, dims3, lize3, start3, MPI_ORDER_FORTRAN, MPI_REAL,        &
                                      ft_3dsoil, ierr )
       CALL MPI_TYPE_COMMIT( ft_3dsoil, ierr )
    ENDIF
#else
    ALLOCATE( array_3d_soil(nzb_soil:nzt_soil,iog%nxl:iog%nxr,iog%nys:iog%nyn) )
    sm_io%io_grid = iog
#endif

 END SUBROUTINE rd_mpi_io_create_filetypes_3dsoil

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Free all file types that have been created for MPI-IO.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_free_filetypes

    IMPLICIT NONE

#if defined( __parallel )
    IF ( filetypes_created )  THEN

       IF ( sm_io%iam_io_pe )  THEN
          CALL MPI_TYPE_FREE( ft_2d, ierr )
          CALL MPI_TYPE_FREE( ft_2di_nb, ierr )
          CALL MPI_TYPE_FREE( ft_2di8_nb, ierr )
          CALL MPI_TYPE_FREE( ft_3d, ierr )
       ENDIF

       IF ( sm_io%is_sm_active() )  THEN
          CALL sm_io%sm_free_shared( win_2dr )
          CALL sm_io%sm_free_shared( win_2di )
          CALL sm_io%sm_free_shared( win_2di8 )
          CALL sm_io%sm_free_shared( win_3dr )
       ELSE
          DEALLOCATE( array_2d, array_2di, array_3d )
       ENDIF

    ENDIF

!
!-- Free last surface filetype
    IF ( sm_io%iam_io_pe .AND. ft_surf /= -1 )  THEN
       CALL MPI_TYPE_FREE( ft_surf, ierr )
    ENDIF

    IF ( win_surf /= -1 )  THEN
       CALL MPI_WIN_FREE( win_surf, ierr )
       DEALLOCATE( array_1d )
       win_surf = -1
    ENDIF

    IF ( sm_io%is_sm_active() )  THEN
       IF ( win_out /= -1 )  THEN
          CALL sm_io%sm_free_shared( win_out )  ! frees shared memory window and releases memory
          win_out = -1
       ENDIF
       NULLIFY( array_out )  ! disassociate Fortran pointer
    ELSE
       IF ( ASSOCIATED( array_out ) )  DEALLOCATE( array_out )
    ENDIF
!
!-- Free last particle filetypes
    IF ( sm_io%iam_io_pe .AND. ft_3di4 /= -1 )  THEN
       CALL MPI_TYPE_FREE( ft_3di4, ierr )
       ft_3di4 = -1
    ENDIF
    IF ( sm_io%iam_io_pe .AND. ft_3di8 /= -1 )  THEN
       CALL MPI_TYPE_FREE( ft_3di8, ierr )
       ft_3di8 = -1
    ENDIF

    IF ( sm_io%is_sm_active() .AND.  win_3di4 /= -1 )  THEN
       CALL sm_io%sm_free_shared( win_3di4 )
       win_3di4 = -1
    ENDIF
    IF ( sm_io%is_sm_active() .AND.  win_3di8 /= -1 )  THEN
       CALL sm_io%sm_free_shared( win_3di8 )
       win_3di8 = -1
    ENDIF

    IF ( win_start /= -1 )  THEN
       CALL sm_io%sm_free_shared( win_start)
       CALL sm_io%sm_free_shared( win_end)
       CALL sm_io%sm_free_shared( win_glost)
       win_start = -1
       win_end   = -1
       win_glost = -1
    ENDIF

    ft_surf  = -1
    win_surf = -1
#else
    IF ( ASSOCIATED( array_2d )   )  DEALLOCATE( array_2d )
    IF ( ASSOCIATED( array_2di )  )  DEALLOCATE( array_2di )
    IF ( ASSOCIATED( array_3d )   )  DEALLOCATE( array_3d )
    IF ( ASSOCIATED( array_3di4 ) )  DEALLOCATE( array_3di4 )
    IF ( ASSOCIATED( array_3di8 ) )  DEALLOCATE( array_3di8 )
#endif

 END SUBROUTINE rd_mpi_io_free_filetypes



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print the restart data file header (MPI-IO format) for debugging.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_print_header

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<


    WRITE (9,*)  'header position after reading the restart file header: ', header_position
    WRITE (9,*)  ' '
    WRITE (9,*)  'restart file header content:'
    WRITE (9,*)  '----------------------------'
    WRITE (9,*)  ' '

    WRITE (9,*)  ' CHARACTER header values   Total number: ', tgh%nr_char
    WRITE (9,*)  ' '
    DO  i = 1, tgh%nr_char
       WRITE ( 9, '(I3,A,1X,A)' )  i, ': ', text_lines(i)(1:80)
    ENDDO
    WRITE (9,*)  ' '

    WRITE (9,*) ' INTEGER header variables and values   Total number: ', tgh%nr_int
    WRITE (9,*)  ' '
    DO  i = 1, tgh%nr_int
       WRITE(9,*)  ' variable: ', int_names(i), '  value: ', int_values(i)
    ENDDO
    WRITE (9,*)  ' '

    WRITE (9,*)  ' REAL header variables and values   Total number: ', tgh%nr_real
    WRITE (9,*)  ' '
    DO  i = 1, tgh%nr_real
       WRITE(9,*)  ' variable: ', real_names(i), '  value: ', real_values(i)
    ENDDO
    WRITE (9,*)  ' '

    WRITE (9,*)  ' Header entries with offset (2d/3d arrays)   Total number: ', tgh%nr_arrays
    WRITE (9,*)  ' '
    DO  i = 1, tgh%nr_arrays
       WRITE(9,*)  ' variable: ', array_names(i), '  offset: ', array_offset(i)
    ENDDO
    WRITE (9,*)  ' '

 END SUBROUTINE rd_mpi_io_print_header



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check if big/little endian data format is used.
!> An int*4 pointer is set to a int*8 variable, the int*8 is set to 1, and then it is checked, if
!> the first 4 bytes of the pointer are equal 1 (little endian) or not.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rd_mpi_io_check_endian( i_endian )

    IMPLICIT NONE

    INTEGER, INTENT(out)                   ::  i_endian  !<
    INTEGER(KIND=8), TARGET                ::  int8      !<

    INTEGER, DIMENSION(1)                  ::  bufshape  !<
    INTEGER(KIND=4), POINTER, DIMENSION(:) ::  int4      !<

    TYPE(C_PTR)                            ::  ptr       !<


    ptr = C_LOC( int8 )
    bufshape(1) = 2
    CALL C_F_POINTER( ptr, int4, bufshape )

    int8 = 1

    IF ( int4(1) == 1 )  THEN
       i_endian = 1    ! Little endian
    ELSE
       i_endian = 2    ! Big endian
    ENDIF

 END SUBROUTINE rd_mpi_io_check_endian

 END MODULE restart_data_mpi_io_mod
