!> @file shared_memory_io_mod.f90
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
!> Handle MPI-IO or NetCDF-IO shared memory arrays.
!> This module performs the organization of new communicators, adapted PE-grids and allocation of
!> shared memory arrays. The IO itself is not done here.
!> TODO: !kk this module will also be used by poisfft_sm. Maybe it should be renamed in
!>       shared_memory_util_mod
!--------------------------------------------------------------------------------------------------!
 MODULE shared_memory_io_mod

#if defined( __parallel )
    USE MPI
#endif

    USE, INTRINSIC ::  ISO_C_BINDING

    USE control_parameters,                                                                        &
        ONLY:  maximum_grid_level,                                                                 &
               message_string,                                                                     &
               mg_switch_to_pe0_level

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nnx,                                                                                &
               nny,                                                                                &
               nnz,                                                                                &
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
               nzb,                                                                                &
               nzt

    USE kinds,                                                                                     &
        ONLY:  dp,                                                                                 &
               idp,                                                                                &
               isp,                                                                                &
               iwp,                                                                                &
               sp,                                                                                 &
               wp

    USE pegrid,                                                                                    &
        ONLY:  comm1dx,                                                                            &
               comm1dy,                                                                            &
               comm2d,                                                                             &
               ierr,                                                                               &
               myid,                                                                               &
               myidx,                                                                              &
               myidy,                                                                              &
               npex,                                                                               &
               npey,                                                                               &
               numprocs,                                                                           &
               pleft,                                                                              &
               pnorth,                                                                             &
               pright,                                                                             &
               psouth,                                                                             &
               sendrecvcount_xy

    IMPLICIT NONE

    PRIVATE

    SAVE

!
!-- Type to store information about the domain decomposition grid
    TYPE, PUBLIC ::  domain_decomposition_grid_features  !<

       INTEGER(iwp) ::  comm2d    !<
       INTEGER(iwp) ::  myid      !<
       INTEGER(iwp) ::  nnx       !<
       INTEGER(iwp) ::  nny       !<
       INTEGER(iwp) ::  nx        !<
       INTEGER(iwp) ::  nxl       !<
       INTEGER(iwp) ::  nxr       !<
       INTEGER(iwp) ::  ny        !<
       INTEGER(iwp) ::  nyn       !<
       INTEGER(iwp) ::  nys       !<
       INTEGER(iwp) ::  numprocs  !<

       CONTAINS

          PROCEDURE, PASS(this), PUBLIC :: activate_grid_from_this_class
          PROCEDURE, PASS(this), PUBLIC :: save_grid_into_this_class

    END TYPE domain_decomposition_grid_features

    TYPE, PUBLIC ::  sm_remote_array

       TYPE(C_PTR)  ::  rem_ptr  !<
       INTEGER(iwp) ::  d1e      !<
       INTEGER(iwp) ::  d1s      !<
       INTEGER(iwp) ::  d2e      !<
       INTEGER(iwp) ::  d2s      !<
       INTEGER(iwp) ::  d3e      !<
       INTEGER(iwp) ::  d3s      !<
       INTEGER(iwp) ::  d4e      !<
       INTEGER(iwp) ::  d4s      !<

    END TYPE sm_remote_array

!
!-- Class definition for shared memory instances.
!-- For every use of shared memory IO, one instance of this class is created.
    TYPE, PUBLIC ::  sm_class  !<

       INTEGER(iwp) ::  sm_blocks_per_node            !< typical configuration, 2 sockets per node
       LOGICAL      ::  no_shared_Memory_in_this_run  !<

       INTEGER(iwp) ::  comm_model            !< communicator of this model run
!
!--    Variables for the shared memory communicator
       INTEGER(iwp), PUBLIC ::  comm_shared   !< communicator for processes with shared array
       INTEGER(iwp), PUBLIC ::  sh_npes       !<
       INTEGER(iwp), PUBLIC ::  sh_rank       !<

!
!--    Variables for the I/O virtual grid
       INTEGER(iwp), PUBLIC ::  comm_io  !< communicator for all IO processes
       INTEGER(iwp), PUBLIC ::  io_npes  !<
       INTEGER(iwp), PUBLIC ::  io_rank  !<
!
!--    Variables for the node local communicator
       INTEGER(iwp) ::  comm_node          !< communicator for all processes of current node
       INTEGER(iwp) ::  io_pe_global_rank  !<
       INTEGER(iwp) ::  n_npes             !<
       INTEGER(iwp) ::  n_rank             !<

       LOGICAL, PUBLIC ::  is_root_pe          !<
       LOGICAL, PUBLIC ::  iam_io_pe = .TRUE.  !< this PE is an IO-PE

       TYPE(domain_decomposition_grid_features), PUBLIC ::  io_grid  !< io grid features, depending on reading from prerun or main run

       CONTAINS

          PRIVATE

          PROCEDURE, PASS(this), PUBLIC ::  is_sm_active
          PROCEDURE, PASS(this), PUBLIC ::  sm_adjust_outer_boundary
          PROCEDURE, PASS(this), PUBLIC ::  sm_free_shared
          PROCEDURE, PASS(this), PUBLIC ::  sm_init_comm
          PROCEDURE, PASS(this), PUBLIC ::  sm_init_data_output_particles
          PROCEDURE, PASS(this), PUBLIC ::  sm_node_barrier
#if defined( __parallel )
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_1d_32
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_1d_64
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_1di
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_2d_32
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_2d_64
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_2di
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_2di_64
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_3d_32
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_3d_64
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_4d_32
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_4d_64
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_3di_32
          PROCEDURE, PASS(this), PUBLIC ::  sm_allocate_shared_3di_64
          PROCEDURE, PASS(this), PUBLIC ::  sm_all_allocate_shared_3d_64

          GENERIC, PUBLIC ::  sm_allocate_shared =>                                                &
                                              sm_allocate_shared_1d_64, sm_allocate_shared_1d_32,  &
                                              sm_allocate_shared_2d_64, sm_allocate_shared_2d_32,  &
                                              sm_allocate_shared_2di,   sm_allocate_shared_2di_64, &
                                              sm_allocate_shared_3d_64, sm_allocate_shared_4d_64,  &
                                              sm_allocate_shared_4d_32, sm_allocate_shared_3d_32,  &
                                              sm_allocate_shared_1di,   sm_allocate_shared_3di_32, &
                                              sm_allocate_shared_3di_64

          GENERIC, PUBLIC ::  sm_all_allocate_shared => sm_all_allocate_shared_3d_64
#endif
    END TYPE sm_class


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create the communicator for shared memory groups and IO-PEs.
!> Setup the grid for shared memory IO.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_init_comm( this, sm_active, comm_input )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(INOUT) ::  this  !< pointer to access internal variables of this call
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  comm_input  !< main model communicator (comm2d)

#if defined( __parallel )
    INTEGER ::  color              !<
    INTEGER ::  max_npes_per_node  !< maximum number of PEs/node
#endif

    LOGICAL, INTENT(IN) ::  sm_active  !< flag to activate shared-memory IO


    this%sm_blocks_per_node = 2

!
!-- Next line is just to avoid compile errors in serial mode because of unused arguments
    IF ( PRESENT( comm_input )  .AND.  sm_active )  CONTINUE

#if defined( __parallel )
    IF ( PRESENT( comm_input ) )  THEN
       this%comm_model = comm_input
    ELSE
       this%comm_model = comm2d
    ENDIF

    this%no_shared_memory_in_this_run = .NOT. sm_active
    this%comm_io = this%comm_model      ! preset in case of non shared-memory-IO

    IF ( this%no_shared_memory_in_this_run )  THEN
       this%iam_io_pe = .TRUE.
       this%sh_rank   = 0
       this%sh_npes   = 1
       RETURN
    ENDIF

!
!-- Determine, how many PEs are running on a node.
    this%iam_io_pe = .FALSE.
    CALL MPI_COMM_SPLIT_TYPE( this%comm_model, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,             &
                              this%comm_node, ierr )
    CALL MPI_COMM_SIZE( this%comm_node, this%n_npes, ierr )
    CALL MPI_COMM_RANK( this%comm_node, this%n_rank, ierr )

    CALL MPI_ALLREDUCE( this%n_npes, max_npes_per_node, 1, MPI_INTEGER, MPI_MAX, this%comm_model,  &
                        ierr )
!
!-- Special configuration on the HLRN-IV system with 4 shared memory blocks/node.
    IF ( max_npes_per_node > 64 )  this%sm_blocks_per_node = 4

!
!-- Divide a node into shared memory groups, depending on the virtual x-y grid
    CALL compute_color( color )
!
!-- If no shared memory IO possible, nothing is left to be done here.
    IF ( this%no_shared_memory_in_this_run )  RETURN

!
!-- Setup the shared memory area
    CALL MPI_COMM_SPLIT( this%comm_node, color, 0, this%comm_shared, ierr )
    CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
    CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

!
!-- Setup the communicator across the nodes depending on the shared memory rank.
!-- All PEs with shared memory rank 0 will be I/O PEs.
    color = this%sh_rank
    CALL MPI_COMM_SPLIT( this%comm_model, color, 0, this%comm_io, ierr )

    IF ( this%comm_io /= MPI_COMM_NULL )  THEN
       CALL MPI_COMM_SIZE( this%comm_io, this%io_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_io, this%io_rank, ierr )
    ELSE
       this%io_npes = -1
       this%io_rank = -1
    ENDIF

    IF ( this%sh_rank == 0 )  THEN
       this%iam_io_pe = .TRUE.
       this%io_pe_global_rank = myid
    ENDIF
    CALL MPI_BCAST( this%io_pe_global_rank, 1, MPI_INTEGER, 0, this%comm_shared, ierr )
#else
    this%iam_io_pe  = .TRUE.
    this%comm_model = comm2d
    this%sh_rank    = 0
    this%sh_npes    = 1
    this%no_shared_memory_in_this_run = .TRUE.
#endif

#if defined( __parallel )
 CONTAINS

 SUBROUTINE compute_color( color )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(OUT) ::  color  !<

    INTEGER(iwp) ::  group_start    !<
    INTEGER(iwp) ::  my_color       !<
    INTEGER(iwp) ::  n              !<
    INTEGER(iwp) ::  pe             !<
    INTEGER(iwp) ::  sh_group_size  !<

    INTEGER(iwp), DIMENSION(4,0:this%n_npes-1) ::  local_dim_s   !<
    INTEGER(iwp), DIMENSION(4,0:this%n_npes-1) ::  local_dim_r   !<

    TYPE(domain_decomposition_grid_features), DIMENSION(32) ::  node_grid  !<

!
!-- No shared memory I/O on one node jobs
    IF ( numprocs < this%n_npes )  THEN
       this%no_shared_memory_in_this_run = .TRUE.
       RETURN
    ENDIF

    local_dim_s = 0
    local_dim_s(1,this%n_rank) = nxl
    local_dim_s(2,this%n_rank) = nxr
    local_dim_s(3,this%n_rank) = nys
    local_dim_s(4,this%n_rank) = nyn

    node_grid%nyn = -1
!
!-- Distribute the x-y layout of all cores of a node to all node processes
    CALL MPI_ALLREDUCE( local_dim_s, local_dim_r, SIZE( local_dim_s ), MPI_INTEGER, MPI_SUM,       &
                        this%comm_node, ierr )
    sh_group_size = ( max_npes_per_node + this%sm_blocks_per_node - 1 ) / this%sm_blocks_per_node

    pe       = 0
!
!-- color is used to split the shared memory communicator into a communicator for io groups.
    my_color = 1

    group_start = pe
    node_grid(my_color)%nxl = local_dim_r(1,group_start)
    node_grid(my_color)%nxr = local_dim_r(2,group_start)
    node_grid(my_color)%nys = local_dim_r(3,group_start)

    DO  n = 1, this%n_npes-1

       pe =  n
       IF ( n > 0  .AND.  MOD( n,sh_group_size ) == 0 )  THEN
!
!--       If group boundary, start new IO group
          node_grid(my_color)%nyn = local_dim_r(4,pe-1)
          my_color = my_color + 1
          group_start = pe
          node_grid(my_color)%nxl = local_dim_r(1,group_start)
          node_grid(my_color)%nxr = local_dim_r(2,group_start)
          node_grid(my_color)%nys = local_dim_r(3,group_start)

       ELSEIF ( local_dim_r(1,pe) /= node_grid(my_color)%nxl )  THEN
!
!--       If nxl changes, start new IO group
          node_grid(my_color)%nyn = local_dim_r(4,pe-1)
          my_color = my_color + 1
          group_start = pe
          node_grid(my_color)%nxl = local_dim_r(1,group_start)
          node_grid(my_color)%nxr = local_dim_r(2,group_start)
          node_grid(my_color)%nys = local_dim_r(3,group_start)
       ENDIF
!
!--    Save values for local PE
       IF ( this%n_rank == pe )  THEN                                 !
          color = my_color
       ENDIF
       IF ( n == this%n_npes-1 )  node_grid(my_color)%nyn = local_dim_r(4,pe)

    ENDDO

    IF ( this%n_rank == 0 )  THEN
       color = 1
    ENDIF

    this%io_grid = node_grid(color)
    this%io_grid%nnx = this%io_grid%nxr - this%io_grid%nxl + 1
    this%io_grid%nny = this%io_grid%nyn - this%io_grid%nys + 1
!
!-- Consider special case where only a single core is placed on a node.
    IF ( this%n_npes == 1 )  THEN
       this%io_grid%nys = nys
       this%io_grid%nyn = nyn
       this%io_grid%nxl = nxl
       this%io_grid%nxr = nxr
       this%io_grid%nnx = this%io_grid%nxr - this%io_grid%nxl + 1
       this%io_grid%nny = this%io_grid%nyn - this%io_grid%nys + 1
    ENDIF

 END SUBROUTINE compute_color
#endif

 END SUBROUTINE sm_init_comm


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializing setup for output of particle time series.
!> This output always uses a shared memory to reduce the number of particle transfers.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_init_data_output_particles( this )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(INOUT) ::  this  !< pointer to access internal variables of this call

#if defined( __parallel )
    INTEGER(iwp) ::  color              !<
    INTEGER(iwp) ::  ierr               !<
    INTEGER(iwp) ::  max_npes_per_node  !< maximum number of PEs/node
#endif

    LOGICAL :: sm_active  !<


    this%sm_blocks_per_node = 2

    sm_active       = .TRUE.   ! particle IO always uses shared memory
    this%comm_model = comm2d

    this%no_shared_memory_in_this_run = .NOT. sm_active
    this%comm_io = this%comm_model  ! preset in case of non shared-memory-IO

    IF ( this%no_shared_memory_in_this_run )  THEN
       this%iam_io_pe = .TRUE.
       RETURN
    ENDIF

#if defined( __parallel )
!
!-- Determine, how many PEs are running on a node.
    this%iam_io_pe = .FALSE.
    IF ( numprocs > 4 )  THEN
       CALL MPI_COMM_SPLIT_TYPE( this%comm_model, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,          &
                                 this%comm_node, ierr )
    ELSE
!
!--    This ELSE branch is a workaraound for non reproducible MPI aborts in MPI_COMM_SPLIT_TYPE
!--    running on 2 PEs. The problem was seen in the combination OpenMPI/gfortran.
!--    For small configurations <= 4 PEs, every MPI process is treated as its own node.
       color = myid
       CALL MPI_COMM_SPLIT( this%comm_model, color, 0, this%comm_node, ierr )
       this%sm_blocks_per_node = 1
    ENDIF

    CALL MPI_COMM_SIZE( this%comm_node, this%n_npes, ierr )
    CALL MPI_COMM_RANK( this%comm_node, this%n_rank, ierr )

    CALL MPI_ALLREDUCE( this%n_npes, max_npes_per_node, 1, MPI_INTEGER, MPI_MAX, this%comm_model,  &
                        ierr )
!
!-- TODO: better explanation
!-- It has to be testet, if using memory blocks for an IO process (MPI shared Memory), or if it is
!-- even better to use the complete node for MPI shared memory (this%sm_blocks_per_node = 1).
!-  In the latter case, the access to the MPI shared memory buffer is slower, the number of
!-- particles to move between PEs will be much smaller.
    IF ( max_npes_per_node > 64 )  THEN
!
!--    Special configuration on the HLRN-IV system with 4 shared memory blocks/node
       this%sm_blocks_per_node = 4
    ENDIF

    IF ( this%sm_blocks_per_node == 1 )  THEN
!
!--    This branch is not realized so far
       this%iam_io_pe   = ( this%n_rank == 0 )
       this%comm_shared = this%comm_node
       CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

    ELSEIF( this%sm_blocks_per_node == 2 )  THEN

       this%iam_io_pe = ( this%n_rank == 0  .OR.  this%n_rank == this%n_npes/2 )
       IF ( this%n_rank < this%n_npes/2 )  THEN
          color = 1
       ELSE
          color = 2
       ENDIF
       CALL MPI_COMM_SPLIT( this%comm_node, color, 0, this%comm_shared, ierr )
       CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

    ELSEIF( this%sm_blocks_per_node == 4 )  THEN

       this%iam_io_pe = ( this%n_rank == 0  .OR.  this%n_rank == this%n_npes/4  .OR.               &
                          this%n_rank == this%n_npes/2  .OR.  this%n_rank == (3*this%n_npes)/4 )
       IF ( this%n_rank < this%n_npes/4 )  THEN
          color = 1
       ELSEIF( this%n_rank < this%n_npes/2 )  THEN
          color = 2
       ELSEIF( this%n_rank < (3*this%n_npes)/4 )  THEN
          color = 3
       ELSE
          color = 4
       ENDIF
       CALL MPI_COMM_SPLIT( this%comm_node, color, 0, this%comm_shared, ierr )
       CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

    ELSE

       WRITE( *, * ) 'shared_memory_io_mod: internal error'
       WRITE( *, * ) 'only 1, 2 or 4 shared memory groups per node are allowed '
       WRITE( *, * ) 'here, ', this%sm_blocks_per_node, ' groups have been set'
       STOP

    ENDIF

!
!-- Setup the shared memory area
    CALL MPI_COMM_SPLIT( this%comm_node, color, 0, this%comm_shared, ierr )
    CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
    CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

!
!-- Setup the communicator across the nodes depending on the shared memory rank.
!-- All PEs with shared memory rank 0 will be I/O PEs.
    color = this%sh_rank
    CALL MPI_COMM_SPLIT( this%comm_model, color, 0, this%comm_io, ierr )

    IF ( this%comm_io /= MPI_COMM_NULL )  THEN
       CALL MPI_COMM_SIZE( this%comm_io, this%io_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_io, this%io_rank, ierr )
    ELSE
       this%io_npes = -1
       this%io_rank = -1
    ENDIF

    IF ( this%sh_rank == 0 )  THEN
       this%iam_io_pe = .TRUE.
       this%io_pe_global_rank = myid
    ENDIF
    CALL MPI_BCAST( this%io_pe_global_rank, 1, MPI_INTEGER, 0, this%comm_shared, ierr )

#else
    this%iam_io_pe = .FALSE.
#endif

 END SUBROUTINE sm_init_data_output_particles

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Function to return if shared Memory IO is active.
!--------------------------------------------------------------------------------------------------!
 FUNCTION is_sm_active( this ) RESULT( ac )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout) ::  this  !<

    LOGICAL ::  ac  !<

    ac = .NOT. this%no_shared_memory_in_this_run

 END FUNCTION is_sm_active


#if defined( __parallel )

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 1d-REAL (64 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE sm_allocate_shared_1d_64( this, p1, d1, d2, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)  ::  this

    INTEGER(iwp), INTENT(IN)        ::  d1
    INTEGER(iwp), INTENT(IN)        ::  d2
    INTEGER(iwp), INTENT(OUT)       ::  win

    INTEGER(iwp)                    ::  disp_unit
    INTEGER(iwp), SAVE              ::  pe_from = 0

    INTEGER(KIND=MPI_ADDRESS_KIND)  ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)  ::  wsize

    INTEGER(idp), DIMENSION(1)      ::  buf_shape

    REAL(dp), DIMENSION(:), POINTER ::  buf
    REAL(dp), DIMENSION(:), POINTER ::  p1

    TYPE(C_PTR), SAVE               ::  base_ptr
    TYPE(C_PTR), SAVE               ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = d2 - d1 + 1
    ELSE
       wsize = 1
    ENDIF
    wsize = wsize * dp  ! please note, size is always in bytes, independently of the displacement
                        ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, dp, MPI_INFO_NULL, this%comm_shared,base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(1) = d2 - d1 + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p1(d1:) => buf

!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_1d_64


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 1d-REAL (32 bit) array on PE 0 and pass address to all PEs
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_1d_32( this, p1, d1, d2, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)  ::  this

    INTEGER(iwp)                    ::  disp_unit
    INTEGER(iwp), INTENT(IN)        ::  d1
    INTEGER(iwp), INTENT(IN)        ::  d2
    INTEGER(iwp), SAVE              ::  pe_from = 0
    INTEGER(iwp), INTENT(OUT)       ::  win

    INTEGER(KIND=MPI_ADDRESS_KIND)  ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)  ::  wsize

    INTEGER(iwp), DIMENSION(1)      ::  buf_shape

    REAL(sp), DIMENSION(:), POINTER ::  buf
    REAL(sp), DIMENSION(:), POINTER ::  p1

    TYPE(C_PTR), SAVE               ::  base_ptr
    TYPE(C_PTR), SAVE               ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = d2 - d1 + 1
    ELSE
       wsize = 1
    ENDIF
    wsize = wsize * sp  ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, sp, MPI_INFO_NULL, this%comm_shared,base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(1) = d2 - d1 + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p1(d1:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_1d_32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 1d-INTEGER array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_1di( this, p1, d1, d2, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)  ::  this

    INTEGER(iwp)                    ::  disp_unit
    INTEGER(iwp), INTENT(IN)        ::  d1
    INTEGER(iwp), INTENT(IN)        ::  d2
    INTEGER(iwp), SAVE              ::  pe_from = 0
    INTEGER(iwp), INTENT(OUT)       ::  win

    INTEGER(KIND=MPI_ADDRESS_KIND)  ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)  ::  wsize

    INTEGER(iwp), DIMENSION(1)          ::  buf_shape

    INTEGER(iwp), DIMENSION(:), POINTER ::  buf
    INTEGER(iwp), DIMENSION(:), POINTER ::  p1

    TYPE(C_PTR), SAVE                   ::  base_ptr
    TYPE(C_PTR), SAVE                   ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = d2 - d1 + 1
    ELSE
       wsize = 1
    ENDIF
    wsize = wsize * iwp  ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, iwp, MPI_INFO_NULL, this%comm_shared,base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(1) = d2 - d1 + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p1(d1:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_1di


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 2d-REAL array (64 bit) on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_2d_64( this, p2, n_nxlg, n_nxrg, n_nysg, n_nyng, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(INOUT)    ::  this

    INTEGER(iwp)                      ::  disp_unit
    INTEGER(iwp), INTENT(IN)          ::  n_nxlg
    INTEGER(iwp), INTENT(IN)          ::  n_nxrg
    INTEGER(iwp), INTENT(IN)          ::  n_nyng
    INTEGER(iwp), INTENT(IN)          ::  n_nysg
    INTEGER(iwp), SAVE                ::  pe_from = 0
    INTEGER(iwp), INTENT(OUT)         ::  win

    INTEGER(KIND=MPI_ADDRESS_KIND)    ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)    ::  wsize

    INTEGER(iwp), DIMENSION(2)        ::  buf_shape

    REAL(dp), DIMENSION(:,:), POINTER ::  buf
    REAL(dp), DIMENSION(:,:), POINTER ::  p2

    TYPE(C_PTR), SAVE                 ::  base_ptr
    TYPE(C_PTR), SAVE                 ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( n_nyng - n_nysg + 1 ) * ( n_nxrg - n_nxlg + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * dp  ! Please note, size is always in bytes, independently of the displacement
                        ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, 8, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(2) = n_nyng - n_nysg + 1
    buf_shape(1) = n_nxrg - n_nxlg + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p2(n_nxlg:, n_nysg:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_2d_64


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 2d-REAL (32 Bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_2d_32( this, p2, n_nxlg, n_nxrg, n_nysg, n_nyng, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(INOUT)    ::  this

    INTEGER(iwp)                      ::  disp_unit
    INTEGER(iwp), INTENT(IN)          ::  n_nxlg
    INTEGER(iwp), INTENT(IN)          ::  n_nxrg
    INTEGER(iwp), INTENT(IN)          ::  n_nyng
    INTEGER(iwp), INTENT(IN)          ::  n_nysg
    INTEGER(iwp), SAVE                ::  pe_from = 0
    INTEGER(iwp), INTENT(OUT)         ::  win

    INTEGER(KIND=MPI_ADDRESS_KIND)    ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)    ::  wsize

    INTEGER(iwp), DIMENSION(2)        ::  buf_shape

    REAL(sp), DIMENSION(:,:), POINTER ::  buf
    REAL(sp), DIMENSION(:,:), POINTER ::  p2

    TYPE(C_PTR), SAVE                 ::  base_ptr
    TYPE(C_PTR), SAVE                 ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( n_nyng - n_nysg + 1 ) * ( n_nxrg - n_nxlg + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * sp  ! Please note, size is always in bytes, independently of the displacement
                        ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, dp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(2) = n_nyng - n_nysg + 1
    buf_shape(1) = n_nxrg - n_nxlg + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p2(n_nxlg:, n_nysg:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_2d_32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 2d-INTEGER array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_2di( this, p2i, n_nxlg, n_nxrg, n_nysg, n_nyng, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)        ::  this         !<

    INTEGER(iwp)                          ::  disp_unit    !<
    INTEGER(iwp), INTENT(IN)              ::  n_nxlg       !<
    INTEGER(iwp), INTENT(IN)              ::  n_nxrg       !<
    INTEGER(iwp), INTENT(IN)              ::  n_nyng       !<
    INTEGER(iwp), INTENT(IN)              ::  n_nysg       !<
    INTEGER(iwp), SAVE                    ::  pe_from = 0  !<
    INTEGER(iwp), INTENT(OUT)             ::  win          !<

    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  wsize        !<

    INTEGER(iwp), DIMENSION(2)            ::  buf_shape    !<

    INTEGER(iwp), DIMENSION(:,:), POINTER ::  buf          !<
    INTEGER(iwp), DIMENSION(:,:), POINTER ::  p2i          !<

    TYPE(C_PTR), SAVE                     ::  base_ptr     !<
    TYPE(C_PTR), SAVE                     ::  rem_ptr      !<


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( n_nyng - n_nysg + 1 ) * ( n_nxrg - n_nxlg + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * 4  ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, 4, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(2) = n_nyng - n_nysg + 1
    buf_shape(1) = n_nxrg - n_nxlg + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p2i(n_nxlg:, n_nysg:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_2di



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 2d-INTEGER*8 array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_2di_64( this, p2i, n_nxlg, n_nxrg, n_nysg, n_nyng, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)        ::  this         !<

    INTEGER(iwp)                          ::  disp_unit    !<
    INTEGER(iwp), INTENT(IN)              ::  n_nxlg       !<
    INTEGER(iwp), INTENT(IN)              ::  n_nxrg       !<
    INTEGER(iwp), INTENT(IN)              ::  n_nyng       !<
    INTEGER(iwp), INTENT(IN)              ::  n_nysg       !<
    INTEGER(iwp), SAVE                    ::  pe_from = 0  !<
    INTEGER(iwp), INTENT(OUT)             ::  win          !<

    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  wsize        !<

    INTEGER(iwp), DIMENSION(2)            ::  buf_shape    !<

    INTEGER(idp), DIMENSION(:,:), POINTER ::  buf          !<
    INTEGER(idp), DIMENSION(:,:), POINTER ::  p2i          !<

    TYPE(C_PTR), SAVE                     ::  base_ptr     !<
    TYPE(C_PTR), SAVE                     ::  rem_ptr      !<


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( n_nyng - n_nysg + 1 ) * ( n_nxrg - n_nxlg + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * idp  ! Please note, size is always in bytes, independently of the displacement
                         ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, idp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(2) = n_nyng - n_nysg + 1
    buf_shape(1) = n_nxrg - n_nxlg + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p2i(n_nxlg:, n_nysg:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_2di_64



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 3d-REAL (64 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_3d_64( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)      ::  this         !<

    INTEGER(iwp)                        ::  disp_unit    !<
    INTEGER(iwp), INTENT(IN)            ::  d1e          !<
    INTEGER(iwp), INTENT(IN)            ::  d1s          !<
    INTEGER(iwp), INTENT(IN)            ::  d2e          !<
    INTEGER(iwp), INTENT(IN)            ::  d2s          !<
    INTEGER(iwp), INTENT(IN)            ::  d3e          !<
    INTEGER(iwp), INTENT(IN)            ::  d3s          !<
    INTEGER(iwp), SAVE                  ::  pe_from = 0  !<
    INTEGER(iwp), INTENT(OUT)           ::  win          !<

    INTEGER(KIND=MPI_ADDRESS_KIND)      ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)      ::  wsize        !<

    INTEGER(iwp), DIMENSION(3)          ::  buf_shape    !<

    REAL(dp), DIMENSION(:,:,:), POINTER ::  buf          !<
    REAL(dp), DIMENSION(:,:,:), POINTER ::  p3           !<

    TYPE(C_PTR), SAVE                   ::  base_ptr     !<
    TYPE(C_PTR), SAVE                   ::  rem_ptr      !<


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * dp ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, dp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_3d_64


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 3d-REAL (32 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_3d_32( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)      ::  this

    INTEGER(iwp)                        ::  disp_unit
    INTEGER(iwp), INTENT(IN)            ::  d1e
    INTEGER(iwp), INTENT(IN)            ::  d1s
    INTEGER(iwp), INTENT(IN)            ::  d2e
    INTEGER(iwp), INTENT(IN)            ::  d2s
    INTEGER(iwp), INTENT(IN)            ::  d3e
    INTEGER(iwp), INTENT(IN)            ::  d3s
    INTEGER(iwp), SAVE                  ::  pe_from = 0
    INTEGER(iwp), INTENT(OUT)           ::  win

    INTEGER(KIND=MPI_ADDRESS_KIND)      ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)      ::  wsize

    INTEGER(iwp), DIMENSION(3)          ::  buf_shape

    REAL(sp), DIMENSION(:,:,:), POINTER ::  buf
    REAL(sp), DIMENSION(:,:,:), POINTER ::  p3

    TYPE(C_PTR), SAVE                   ::  base_ptr
    TYPE(C_PTR), SAVE                   ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * sp ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, sp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_3d_32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 4d-REAL (64 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_4d_64( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, d4s, d4e, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)        ::  this         !<

    INTEGER                               ::  disp_unit    !<
    INTEGER(iwp), INTENT(IN)              ::  d1e          !<
    INTEGER(iwp), INTENT(IN)              ::  d1s          !<
    INTEGER(iwp), INTENT(IN)              ::  d2e          !<
    INTEGER(iwp), INTENT(IN)              ::  d2s          !<
    INTEGER(iwp), INTENT(IN)              ::  d3e          !<
    INTEGER(iwp), INTENT(IN)              ::  d3s          !<
    INTEGER(iwp), INTENT(IN)              ::  d4e          !<
    INTEGER(iwp), INTENT(IN)              ::  d4s          !<
    INTEGER(iwp), SAVE                    ::  pe_from = 0  !<
    INTEGER(iwp), INTENT(OUT)             ::  win          !<

    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  wsize        !<

    INTEGER(iwp), DIMENSION(4)            ::  buf_shape    !<

    REAL(dp), DIMENSION(:,:,:,:), POINTER ::  buf          !<
    REAL(dp), DIMENSION(:,:,:,:), POINTER ::  p3           !<

    TYPE(C_PTR), SAVE                     ::  base_ptr     !<
    TYPE(C_PTR), SAVE                     ::  rem_ptr      !<


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = (d4e - d4s +1) * ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * dp ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, dp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(4) = d4e - d4s + 1
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:,d4s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_4d_64


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 4d-REAL (32 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_4d_32( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, d4s, d4e, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)        ::  this         !<

    INTEGER                               ::  disp_unit    !<
    INTEGER(iwp), INTENT(IN)              ::  d1e          !<
    INTEGER(iwp), INTENT(IN)              ::  d1s          !<
    INTEGER(iwp), INTENT(IN)              ::  d2e          !<
    INTEGER(iwp), INTENT(IN)              ::  d2s          !<
    INTEGER(iwp), INTENT(IN)              ::  d3e          !<
    INTEGER(iwp), INTENT(IN)              ::  d3s          !<
    INTEGER(iwp), INTENT(IN)              ::  d4e          !<
    INTEGER(iwp), INTENT(IN)              ::  d4s          !<
    INTEGER(iwp), SAVE                    ::  pe_from = 0  !<
    INTEGER(iwp), INTENT(OUT)             ::  win          !<

    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)        ::  wsize        !<

    INTEGER(iwp), DIMENSION(4)            ::  buf_shape    !<

    REAL(sp), DIMENSION(:,:,:,:), POINTER ::  buf          !<
    REAL(sp), DIMENSION(:,:,:,:), POINTER ::  p3           !<

    TYPE(C_PTR), SAVE                     ::  base_ptr     !<
    TYPE(C_PTR), SAVE                     ::  rem_ptr      !<


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = (d4e - d4s +1) * ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * sp ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, dp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(4) = d4e - d4s + 1
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:,d4s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_4d_32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 3d-INTEGER (32 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_3di_32( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)          ::  this

    INTEGER                                 ::  disp_unit
    INTEGER(iwp), INTENT(IN)                ::  d1e
    INTEGER(iwp), INTENT(IN)                ::  d1s
    INTEGER(iwp), INTENT(IN)                ::  d2e
    INTEGER(iwp), INTENT(IN)                ::  d2s
    INTEGER(iwp), INTENT(IN)                ::  d3e
    INTEGER(iwp), INTENT(IN)                ::  d3s
    INTEGER(iwp), SAVE                      ::  pe_from = 0
    INTEGER(iwp), INTENT(OUT)               ::  win

    INTEGER(KIND=MPI_ADDRESS_KIND)          ::  rem_size
    INTEGER(KIND=MPI_ADDRESS_KIND)          ::  wsize

    INTEGER(iwp), DIMENSION(3)              ::  buf_shape

    INTEGER(isp), DIMENSION(:,:,:), POINTER ::  buf
    INTEGER(isp), DIMENSION(:,:,:), POINTER ::  p3

    TYPE(C_PTR), SAVE                       ::  base_ptr
    TYPE(C_PTR), SAVE                       ::  rem_ptr


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * isp ! Please note, size is always in bytes, independently of the displacement
                       ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, isp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_3di_32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 3d-INTEGER (64 bit) array on PE 0 and pass address to all PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_allocate_shared_3di_64( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)          ::  this         !<

    INTEGER                                 ::  disp_unit    !<
    INTEGER(iwp), INTENT(IN)                ::  d1e          !<
    INTEGER(iwp), INTENT(IN)                ::  d1s          !<
    INTEGER(iwp), INTENT(IN)                ::  d2e          !<
    INTEGER(iwp), INTENT(IN)                ::  d2s          !<
    INTEGER(iwp), INTENT(IN)                ::  d3e          !<
    INTEGER(iwp), INTENT(IN)                ::  d3s          !<
    INTEGER(iwp), SAVE                      ::  pe_from = 0  !<
    INTEGER(iwp), INTENT(OUT)               ::  win          !<

    INTEGER(KIND=MPI_ADDRESS_KIND)          ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)          ::  wsize        !<

    INTEGER(iwp), DIMENSION(3)              ::  buf_shape    !<

    INTEGER(idp), DIMENSION(:,:,:), POINTER ::  buf          !<
    INTEGER(idp), DIMENSION(:,:,:), POINTER ::  p3           !<

    TYPE(C_PTR), SAVE                       ::  base_ptr     !<
    TYPE(C_PTR), SAVE                       ::  rem_ptr      !<


    IF ( this%no_shared_memory_in_this_run )  RETURN
!
!-- Allocate shared memory on node rank 0 PEs.
    IF ( this%sh_rank == pe_from )  THEN
       wsize = ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )
    ELSE
       wsize = 1
    ENDIF

    wsize = wsize * idp ! Please note, size is always in bytes, independently of the displacement
                        ! unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, idp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)
    CALL MPI_WIN_SHARED_QUERY( win, pe_from, rem_size, disp_unit, rem_ptr, ierr )
!
!-- Convert C- to Fortran-pointer
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_allocate_shared_3di_64


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate shared 3d-REAL (64 Bit) array on ALL PEs.
!>
!> Every PE allocates the local part of a node-shared array.
!> The C-Pointer of this array and the local limits are broadcasted to all PEs of the node
!> The information is store in an array of type sm_remote_array and can be retrieved
!> by sm_remote_array to access remote data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_all_allocate_shared_3d_64( this, p3, d1s, d1e, d2s, d2e, d3s, d3e, remote_arrays, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout)      ::  this         !< class pointer
    REAL(dp), DIMENSION(:,:,:), POINTER ::  p3           !< return local array pointer

    INTEGER(iwp), INTENT(IN)            ::  d1e          !< end index dimension 1
    INTEGER(iwp), INTENT(IN)            ::  d1s          !< start index dimension 1
    INTEGER(iwp), INTENT(IN)            ::  d2e          !< end index dimension 2
    INTEGER(iwp), INTENT(IN)            ::  d2s          !< start index dimension 2
    INTEGER(iwp), INTENT(IN)            ::  d3e          !< end index dimension 3
    INTEGER(iwp), INTENT(IN)            ::  d3s          !< start index dimension 3
    INTEGER(iwp), INTENT(OUT)           ::  win          !< MPI Window

    INTEGER(iwp), DIMENSION(3)          ::  buf_shape    !<
    INTEGER(iwp)                        ::  disp_unit    !<
    INTEGER(iwp)                        ::  i            !<
    INTEGER(iwp), SAVE                  ::  pe_from = 0  !<

    INTEGER(KIND=MPI_ADDRESS_KIND)      ::  rem_size     !<
    INTEGER(KIND=MPI_ADDRESS_KIND)      ::  wsize        !<

    REAL(dp), DIMENSION(:,:,:), POINTER ::  buf          !<

    TYPE(sm_remote_array),INTENT(INOUT), DIMENSION(0:this%sh_npes-1) :: remote_arrays !< info about all remote arrays

    TYPE(C_PTR), SAVE                   ::  base_ptr     !<

    INTEGER(iwp),DIMENSION(6,0:this%sh_npes-1)              ::  all_indices_s
    INTEGER(iwp),DIMENSION(6,0:this%sh_npes-1)              ::  all_indices


    IF ( this%no_shared_memory_in_this_run )  RETURN

    all_indices_s = 0


    wsize = ( d3e - d3s + 1 ) * ( d2e - d2s + 1 ) * ( d1e - d1s + 1 )

    wsize = wsize * dp   ! Please note, size is always in bytes, independently of the displacement unit

    CALL MPI_WIN_ALLOCATE_SHARED( wsize, dp, MPI_INFO_NULL, this%comm_shared, base_ptr, win, ierr )
!
!-- Get C-pointer of the memory located on node-rank pe_from (sh_rank == pe_from)

    all_indices_s(1,this%sh_rank) = d1s
    all_indices_s(2,this%sh_rank) = d1e
    all_indices_s(3,this%sh_rank) = d2s
    all_indices_s(4,this%sh_rank) = d2e
    all_indices_s(5,this%sh_rank) = d3s
    all_indices_s(6,this%sh_rank) = d3e

    CALL MPI_ALLREDUCE (all_indices_s ,all_indices, SIZE(all_indices_s), MPI_INTEGER, MPI_SUM, this%comm_shared, ierr)

    DO i=0,this%sh_npes-1
       CALL MPI_WIN_SHARED_QUERY( win, i, rem_size, disp_unit, remote_arrays(i)%rem_ptr, ierr )
       remote_arrays(i)%d1s = all_indices(1,i)
       remote_arrays(i)%d1e = all_indices(2,i)
       remote_arrays(i)%d2s = all_indices(3,i)
       remote_arrays(i)%d2e = all_indices(4,i)
       remote_arrays(i)%d3s = all_indices(5,i)
       remote_arrays(i)%d3e = all_indices(6,i)
    END DO

!
!-- Convert C- to Fortran-pointer
    buf_shape(3) = d3e - d3s + 1
    buf_shape(2) = d2e - d2s + 1
    buf_shape(1) = d1e - d1s + 1
    CALL C_F_POINTER( remote_arrays(this%sh_rank)%rem_ptr, buf, buf_shape )
    p3(d1s:,d2s:,d3s:) => buf
!
!-- Allocate shared memory in round robin on all PEs of a node.
    pe_from = MOD( pe_from, this%sh_npes )

 END SUBROUTINE sm_all_allocate_shared_3d_64
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> ???
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_adjust_outer_boundary( this )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout) ::  this  !<


    IF ( this%no_shared_memory_in_this_run )  RETURN

    IF ( this%io_grid%nxl == 0 )  THEN
       this%io_grid%nxl = this%io_grid%nxl - nbgp
       this%io_grid%nnx = this%io_grid%nnx + nbgp
    ENDIF

    IF ( this%io_grid%nxr == nx  .OR.  npex == -1 )  THEN   ! npex == -1 if -D__parallel not set
       this%io_grid%nxr = this%io_grid%nxr + nbgp
       this%io_grid%nnx = this%io_grid%nnx + nbgp
    ENDIF

    IF ( this%io_grid%nys == 0 )  THEN
       this%io_grid%nys = this%io_grid%nys - nbgp
       this%io_grid%nny = this%io_grid%nny + nbgp
    ENDIF

    IF ( this%io_grid%nyn == ny .OR.  npey == -1 )  THEN   ! npey == -1 if -D__parallel not set
       this%io_grid%nyn = this%io_grid%nyn + nbgp
       this%io_grid%nny = this%io_grid%nny + nbgp
    ENDIF

    this%io_grid%nxl = this%io_grid%nxl + nbgp
    this%io_grid%nxr = this%io_grid%nxr + nbgp
    this%io_grid%nys = this%io_grid%nys + nbgp
    this%io_grid%nyn = this%io_grid%nyn + nbgp
    this%io_grid%nnx = this%io_grid%nnx
    this%io_grid%nny = this%io_grid%nny

 END SUBROUTINE sm_adjust_outer_boundary


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocate shared aray and free related window.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_free_shared( this, win )

    IMPLICIT NONE

    CLASS(sm_class), INTENT(inout) ::  this  !<

    INTEGER(iwp), INTENT(INOUT)    ::  win   !<

    IF ( this%no_shared_memory_in_this_run )  RETURN
#if defined( __parallel )
    CALL MPI_WIN_FREE( win, ierr )
#endif
    win = -1

 END SUBROUTINE sm_free_shared


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> ...
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_node_barrier( this, win )

    IMPLICIT NONE

    INTEGER(iwp), OPTIONAL         ::  win   !<

    CLASS(sm_class), INTENT(inout) ::  this  !<


    IF ( this%no_shared_memory_in_this_run )  RETURN

#if defined( __parallel )
    CALL MPI_BARRIER( this%comm_shared, ierr )
#endif
    IF ( PRESENT(win) )  THEN
#if defined( __parallel )
       CALL MPI_WIN_FENCE(0, win, ierr )
#else
       CONTINUE
#endif
    ENDIF

 END SUBROUTINE sm_node_barrier


 SUBROUTINE save_grid_into_this_class( this )

    IMPLICIT NONE

    CLASS(domain_decomposition_grid_features), INTENT(inout) ::  this  !<

       this%myid     = myid      !<
       this%nnx      = nnx       !<
       this%nny      = nny       !<
       this%nx       = nx        !<
       this%nxl      = nxl       !<
       this%nxr      = nxr       !<
       this%ny       = ny        !<
       this%nyn      = nyn       !<
       this%nys      = nys       !<
       this%numprocs = numprocs  !<
       this%comm2d   = comm2d    !<

 END SUBROUTINE save_grid_into_this_class


 SUBROUTINE activate_grid_from_this_class( this )

    IMPLICIT NONE

    CLASS(domain_decomposition_grid_features), INTENT(inout) ::  this  !<

       myid     = this%myid      !<
       nnx      = this%nnx       !<
       nny      = this%nny       !<
       nx       = this%nx        !<
       nxl      = this%nxl       !<
       nxr      = this%nxr       !<
       ny       = this%ny        !<
       nyn      = this%nyn       !<
       nys      = this%nys       !<
       numprocs = this%numprocs  !<
       comm2d   = this%comm2d    !<

 END SUBROUTINE activate_grid_from_this_class

 END MODULE shared_memory_io_mod
