!> @file sm_poisfft_mod.f90
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
! Copyright 1997-2022 Leibniz Universitaet Hannover, Klaus Ketelsen
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Routines to setup and activate settings for the shared memory Poisson-FFT-solver.
!--------------------------------------------------------------------------------------------------!
 MODULE sm_poisfft_mod

#if defined( __parallel )
    USE MPI

    USE control_parameters,                                                                        &
        ONLY:  debug_output,                                                                       &
               message_string

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

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  comm1dx,                                                                            &
               comm1dy,                                                                            &
               comm2d,                                                                             &
               comm_palm,                                                                          &
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

   USE transpose_mod,                                                                              &
       ONLY:  nxl_y,                                                                               &
              nxl_z,                                                                               &
              nxr_y,                                                                               &
              nxr_z,                                                                               &
              nys_x,                                                                               &
              nys_z,                                                                               &
              nyn_x,                                                                               &
              nyn_z,                                                                               &
              nzb_x,                                                                               &
              nzb_y,                                                                               &
              nzt_x,                                                                               &
              nzt_y

   USE shared_memory_io_mod,                                                                       &
       ONLY:  domain_decomposition_grid_features,                                                  &
              sm_class

   IMPLICIT NONE

   TYPE, PUBLIC, EXTENDS(domain_decomposition_grid_features) ::  sm_setup  !<

      INTEGER(iwp) ::  npex              !<
      INTEGER(iwp) ::  npey              !<
      INTEGER(iwp) ::  nys_x             !<
      INTEGER(iwp) ::  nyn_x             !<
      INTEGER(iwp) ::  nys_z             !<
      INTEGER(iwp) ::  nyn_z             !<
      INTEGER(iwp) ::  nxl_z             !<
      INTEGER(iwp) ::  nxr_z             !<
      INTEGER(iwp) ::  sendrecvcount_xy  !<

      CONTAINS

          PROCEDURE, PASS(this), PUBLIC ::  activate_grid_from_this_class_poisfft
          PROCEDURE, PASS(this), PUBLIC ::  save_grid_into_this_class_poisfft

   END TYPE sm_setup


   TYPE, PUBLIC, EXTENDS(sm_class) ::  sm_poisfft_class  !<

      INTEGER ::  node_pe0  !<

      CONTAINS

         PRIVATE

         PROCEDURE, PASS(this), PUBLIC ::  sm_check_layout
         PROCEDURE, PASS(this), PUBLIC ::  sm_init_pegrid_poisfft

   END TYPE sm_poisfft_class

   TYPE(sm_poisfft_class), PUBLIC ::  sm_poisfft  !<


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines the virtual PE grid (npex,npey) in case that the shared memory Poisson-FFT-solver
!> is used.
!> In such a case the number of PEs along y (npey) must exactly match the number of cores belonging
!> to a shared memory block. The largest possible shared memory block is given by the number of
!> cores used on a node, but 2 or 4 times smaller blocks will be automatically selected for
!> performance reasons, depending on the total number of cores that are used by a run.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_init_pegrid_poisfft( this, sm_active )

    IMPLICIT NONE

    CLASS(sm_poisfft_class), INTENT(INOUT) ::  this  !< pointer to access internal variables of this call

    LOGICAL,INTENT(IN) ::  sm_active

    INTEGER(iwp) ::  color        !<
    INTEGER(iwp) ::  comm_tmp     !< temporary MPI communicator
    INTEGER(iwp) ::  ierr         !< MPI error tag
    INTEGER(iwp) ::  key          !<
    INTEGER(iwp) ::  max_sh_npes  !<
    INTEGER(iwp) ::  new_npes     !<
    INTEGER(iwp) ::  new_rank     !<


!
!-- Set flag and return, if the shared memory Poisson solver has not been selected.
    this%no_shared_memory_in_this_run = .NOT. sm_active
    IF ( this%no_shared_memory_in_this_run )  RETURN

!
!-- Create subcommunicator that spans the shared-memoy domain of the largest possible size.
    CALL  MPI_COMM_SPLIT_TYPE( comm2d, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, this%comm_shared,   &
                               ierr )
!
!-- Calculate of number of cores belonging to this subcommunicator and their respective rank.
!-- This values will be re-defined further below, if nodes are split into several (2 or 4)
!-- smaller shared-memory domains.
    CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
    CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

    CALL MPI_ALLREDUCE( this%sh_npes, max_sh_npes, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )

!
!-- In order to split the shared memory of a node into a maximum of 4 smaller blocks, the number of
!-- PEs per node must be a multiple of 4.
    IF ( MOD( max_sh_npes, 4 ) /= 0 )  THEN
       WRITE( message_string, * ) 'psolver = "poisfft_sm" requires the number of PEs per node',    &
                                  ' to be a multiple of 4 but this run uses ', max_sh_npes,         &
                                  ' PEs per node'
       CALL message( 'sm_init_pegrid_poisfft', 'PAC0306', 3, 2, 0, 6, 1 )
    ENDIF

    IF ( debug_output )  THEN
       WRITE( 9, * )  'sm_init_pegrid_poisfft after MPI_COMM_SPLIT_TYPE: this%sh_npes = ',         &
                      this%sh_npes, ' this%sh_rank = ', this%sh_rank
    ENDIF

!
!-- Decide, for performance reasons, if the shared memory of a node is split into 2 or 4 smaller
!-- blocks. Run will be aborted if no virtual PE grid with good performance can be found.
    CALL determine_sm_blocks_per_node

!
!-- Re-define comm_shared if it has been decided above to split the shared memory of a node
!-- into smaller blocks.
    IF ( this%sm_blocks_per_node == 2 )  THEN

       comm_tmp = this%comm_shared
       IF ( this%sh_npes == max_sh_npes/2 )  THEN
!
!--       TODO: unclear, if the extra check for this%sh_npes == max_sh_npes/2 is really
!--             required. Can perhaps be removed because there is no respective check in the
!--             below case for 4 shared memory blocks.
          color = 1
       ELSE
          IF ( this%sh_rank < this%sh_npes/2 )  THEN
             color = 1
          ELSE
             color = 2
          ENDIF
       ENDIF

       CALL MPI_COMM_SPLIT( comm_tmp, color, 0, this%comm_shared, ierr )
       CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

    ELSEIF ( this%sm_blocks_per_node == 4 )  THEN

       comm_tmp = this%comm_shared
       IF ( this%sh_rank < this%sh_npes/4 )  THEN
          color = 1
       ELSEIF ( this%sh_rank < this%sh_npes/2 )  THEN
          color = 2
       ELSEIF ( this%sh_rank < (3*this%sh_npes)/4 )  THEN
          color = 3
       ELSE
          color = 4
       ENDIF

       CALL MPI_COMM_SPLIT( comm_tmp, color, 0, this%comm_shared, ierr )
       CALL MPI_COMM_SIZE( this%comm_shared, this%sh_npes, ierr )
       CALL MPI_COMM_RANK( this%comm_shared, this%sh_rank, ierr )

    ENDIF

!
!-- Set root PE flag (PEs with rank 0) for the above generated shared memory blocks..
    IF ( this%sh_rank == 0 )  this%is_root_pe = .TRUE.

!
!-- Create a communicator, which mimics the communication pattern between the shared memory groups.
!-- For the real communication in the shared memory Poisson-FFT-solver comm1dx is used lateron.
!-- sh_rank is the rank within one of the shared memory groups, which means that the size of
!-- comm_node is equal numprocs / sh_npes = npex.
    color = this%sh_rank
    CALL MPI_COMM_SPLIT( comm2d, color, 0, this%comm_node, ierr )

    IF ( this%comm_node == MPI_COMM_NULL )  THEN
       message_string = 'cannot create communicator for shared memory Poisson-FFT-solver'
       CALL message( 'sm_init_pegrid_poisfft', 'PAC0307', 3, 2, 0, 6, 1 )
    ENDIF

!
!-- Size n_npes equals npex, and n_rank is the rank of the shared memory block, i.e. we have
!-- npex shared memory blocks.
    CALL MPI_COMM_SIZE( this%comm_node, this%n_npes, ierr )
    CALL MPI_COMM_RANK( this%comm_node, this%n_rank, ierr )

!
!-- Check, if the calculated topology ratio npex/npey is too extreme, and then stop for
!-- performance reasons. The factor 8 below is estimated empirically and may need to be adjusted.
    IF ( REAL( this%n_npes ) / REAL( this%sh_npes ) > 8.0_wp )  THEN
       WRITE( message_string, * ) 'calculated PE grid topology ratio npex/npey > 8',               &
                                  'may result in poor performance & npex =', this%n_npes,          &
                                  ' npey = ', this%sh_npes
       CALL message( 'sm_init_pegrid_poisfft', 'PAC0308', 3, 2, 0, 6, 1 )
    ENDIF

    IF ( debug_output )  THEN
       IF ( this%sh_rank == 0 )  this%node_pe0 = myid
       CALL MPI_BCAST(this%node_pe0, 1, MPI_INTEGER, 0, this%comm_shared, ierr )
       WRITE( 9, '(A,8I7)' )  'sm_init_pegrid_poisfft ', this%sh_rank, this%sh_npes, this%n_rank,  &
                              this%n_npes, this%node_pe0
    ENDIF

!
!-- Finally set the virtual PE-grid to be used for the 2d domain decomposition. Be aware that the
!-- shared memory Poisson solver itself is using a 1d decomposition where npey = 1.
    npex = this%n_npes
    npey = this%sh_npes

!
!-- Make sure that the PE order of the PALM model communicator comm_palm (later comm2d) fits the
!-- order given by the shared memory layout. This might e.g. be violated if MPI_CART_CREATE is
!-- allowed a re-ordering of processes.
!-- The order must be first along y, and second along x (see calculation of key).
    comm_tmp = comm_palm
    color = 1
    key = this%n_rank*this%sh_npes + this%sh_rank
    CALL MPI_COMM_SPLIT( comm_tmp, color, key, comm_palm, ierr )

    IF ( debug_output )  THEN
       CALL MPI_COMM_SIZE( comm_palm, new_npes, ierr )
       CALL MPI_COMM_RANK( comm_palm, new_rank, ierr )
       WRITE( 9, '(A,8I6)' )  'PE layout ', this%sh_rank, this%sh_npes, this%n_rank, this%n_npes,  &
                              new_rank, new_npes
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines the number of shared memory blocks to be used on one node. This number depends on the
!> number of cores that physically exist on one node, as well as on the total number of cores that
!> are used in the run.
!> The calculation itself is a general one, but it is optimized here for a system with
!> 96 cores/node. Other systems may profit from fine tuning.
!> The calculation below results in
!>
!>    this%sm_blocks_per_node = 4  for    1  -   6 full nodes
!>    this%sm_blocks_per_node = 2  for    7  -  48 full nodes
!>    this%sm_blocks_per_node = 1  for   49  - 384 full nodes
!>
!> If no full nodes are used, the number of shared memory blocks per node is adjusted (see code).
!> If more than 384 full nodes are used, a switch back to the default poisfft (without shared
!> memory usage) is done.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE determine_sm_blocks_per_node

    INTEGER(iwp) ::  n1   !<
    INTEGER(iwp) ::  n2   !<
    INTEGER(iwp) ::  n4   !<
    INTEGER(iwp) ::  nx1  !<
    INTEGER(iwp) ::  nx2  !<
    INTEGER(iwp) ::  nx4  !<


!
!-- Calculate the number of cores per shared memory block, for the case that a node would be
!-- split into 1, 2 or 4 shared memory blocks.
    n1 = max_sh_npes
    n2 = max_sh_npes / 2
    n4 = max_sh_npes / 4
!
!-- Calculate the total number of shared memory blocks in this run, for the case that 1, 2,
!-- or 4 shared memory blocks per node would be used.
    nx1 = numprocs / n1
    nx2 = numprocs / n2
    nx4 = numprocs / n4

!
!-- Setup for small systems. Thresholds of 32 and 16 are estimated but not that important, since
!-- the Poisson solver does not take much time at all for such small setups at all.
!-- Two shared memory block per node are chosen to allow for a testing on the testserver.
    IF ( numprocs <= 32  .OR.  n1 <= 16 )  THEN
       this%sm_blocks_per_node = 2
       RETURN
    ENDIF

!
!-- The below conditions try to make the virtual PE grid as quadratic as possible.
!-- The number of cores per shared memory block is equivalent to npey, the total number of shared
!-- memory block is equivalnet to npex.
    IF ( nx4 < n4 )  THEN
       this%sm_blocks_per_node = 4
    ELSEIF ( nx2 < n2*2 )  THEN
       this%sm_blocks_per_node = 2
    ELSEIF ( nx1 < n1*4 )  THEN
       this%sm_blocks_per_node = 1
    ELSE
!
!--    The total number of cores is too large, or the size of the largest shared memory block
!--    (which is the number of cores per node) is too small.
       message_string = 'no virtual PE grid with sufficient performance can be generated'
       CALL message( 'determine_sm_blocks_per_node', 'PAC0309', 3, 2, 0, 6, 1 )
    ENDIF

    IF ( debug_output )  THEN
       WRITE( 9,'(A,8I6,I8)' ) 'determine_sm_blocks_per_node ', n1, n2, n4, nx1, nx2, nx4,         &
                               max_sh_npes, this%sm_blocks_per_node, numprocs
       FLUSH( 9 )
    ENDIF

 END SUBROUTINE determine_sm_blocks_per_node

 END SUBROUTINE sm_init_pegrid_poisfft


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compare, if the virtual shared memory solver layout fits the virtual PE grid.
!> The PE rank of a shared memory block must fit myidy, the node rank of the shared memory block
!> must fit myidx.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sm_check_layout( this )

    IMPLICIT NONE

    CLASS(sm_poisfft_class), INTENT(inout) ::  this  !<


    IF ( this%sh_rank /= myidy  .OR.  this%n_rank /= myidx )  THEN
       WRITE( message_string,'(A,4(A4,I6))' ) 'shared memory virtual layout error &',              &
                                              'shared rank / myidy: ', this%sh_rank, ' / ', myidy, &
                                              'node rank / myidx: ', this%n_rank, ' / ', myidx
       CALL message( 'init_pegrid', 'PAC0310', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE sm_check_layout


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE save_grid_into_this_class_poisfft( this )

    IMPLICIT NONE

    CLASS(sm_setup), INTENT(inout) ::  this  !<


    CALL this%save_grid_into_this_class()

    this%npex = npex
    this%npey = npey
    this%nys_x = nys_x
    this%nyn_x = nyn_x
    this%nys_z = nys_z
    this%nyn_z = nyn_z
    this%nxl_z = nxl_z
    this%nxr_z = nxr_z
    this%sendrecvcount_xy = sendrecvcount_xy

 END SUBROUTINE save_grid_into_this_class_poisfft


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE activate_grid_from_this_class_poisfft( this )

    IMPLICIT NONE

    CLASS(sm_setup), INTENT(inout) ::  this  !<


    CALL this%activate_grid_from_this_class()

    npex = this%npex
    npey = this%npey
    nys_x = this%nys_x
    nyn_x = this%nyn_x
    nys_z = this%nys_z
    nyn_z = this%nyn_z
    nxl_z = this%nxl_z
    nxr_z = this%nxr_z
    sendrecvcount_xy = this%sendrecvcount_xy

 END SUBROUTINE activate_grid_from_this_class_poisfft
#endif

END MODULE sm_poisfft_mod
