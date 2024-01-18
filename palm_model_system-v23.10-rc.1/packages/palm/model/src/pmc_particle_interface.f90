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
! Introduce particle transfer in nested models.
! Interface to palm lpm model to handel particle transfer between parent and child model.
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_particle_interface
  
#if defined( __parallel )

    USE, INTRINSIC ::  ISO_C_BINDING

    USE MPI

    USE arrays_3d,                                                                                 &
        ONLY:  zw

    USE control_parameters,                                                                        &
        ONLY:  debug_output,                                                                       &
               debug_string,                                                                       &
               message_string

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               nxlg,                                                                               &
               nxrg,                                                                               &
               ny,                                                                                 &
               nys,                                                                                &
               nyn,                                                                                &
               nysg,                                                                               &
               nyng,                                                                               &
               nzb,                                                                                &
               nzt

    USE kinds

    USE particle_attributes,                                                                       &
        ONLY:  alloc_factor,                                                                       &
               grid_particles,                                                                     &
               ibc_par_lr,                                                                         &
               ibc_par_ns,                                                                         &
               ibc_par_t,                                                                          &
               particles,                                                                          &
               particle_type,                                                                      &
               prt_count,                                                                          &
               number_of_particles,                                                                &
               zero_particle

    USE pegrid,                                                                                    &
        ONLY:  myid,                                                                               &
               myidx,                                                                              &
               myidy,                                                                              &
               numprocs

#if defined( __parallel )
    USE pmc_child,                                                                                 &
        ONLY:  me,                                                                                 &
               pmc_c_getbuffer,                                                                    &
               pmc_c_putbuffer

    USE pmc_general,                                                                               &
        ONLY:  pedef

    USE pmc_handle_communicator,                                                                   &
        ONLY:  pmc_parent_for_child,                                                               &
               pmc_get_model_info

    USE pmc_interface,                                                                             &
        ONLY:  coord_x,                                                                            &
               coord_y,                                                                            &
               cpl_id,                                                                             &
               get_childid,                                                                        &
               get_child_edges,                                                                    &
               get_child_gridspacing,                                                              &
               get_number_of_children,                                                             &
               lower_left_coord_x,                                                                 &
               lower_left_coord_y,                                                                 &
               nr_part,                                                                            &
               nr_partc,                                                                           &
               parent_bound,                                                                       &
               particle_coupling,                                                                  &
               part_adr,                                                                           &
               part_adrc,                                                                          &
               pg,                                                                                 &
               comm_world_nesting,                                                                 &
               root_model

    USE pmc_mpi_wrapper,                                                                           &
        ONLY:  pmc_recv_from_child,                                                                &
               pmc_send_to_parent

    USE pmc_parent,                                                                                &
        ONLY:  children,                                                                           &
               pmc_s_fillbuffer,                                                                   &
               pmc_s_getdata_from_buffer,                                                          &
               pmc_s_get_child_npes
#endif

    IMPLICIT NONE

    PRIVATE
    SAVE

    TYPE  coarse_particle_def
       INTEGER(iwp) ::  nr_particle  !<

       TYPE(particle_type),ALLOCATABLE,DIMENSION(:) ::  parent_particles  !<
    END TYPE  coarse_particle_def

    INTEGER(iwp),PARAMETER ::  max_nr_prts_in_rma_win = 100000  !< maximum size of RMA window to store particles for transfer
    INTEGER(iwp),PARAMETER ::  min_particles_per_column   = 100     !<


    INTEGER(iwp) ::  nr_fine_in_coarse   !< Number of fine grid cells in coarse grid (one direction)
    INTEGER(iwp) ::  particle_win_child  !<

    INTEGER(iwp),ALLOCATABLE,DIMENSION(:) ::  particle_win_parent  !<

    TYPE(C_PTR), ALLOCATABLE,DIMENSION(:) ::  buf_ptr  !<

    TYPE(particle_type), DIMENSION(:),POINTER ::  particle_in_win  !<

    TYPE(coarse_particle_def),ALLOCATABLE,DIMENSION(:,:) ::  coarse_particles  !<

    TYPE communicators_for_nested_pts_output
       INTEGER(iwp) ::  my_npes              !<
       INTEGER(iwp) ::  my_rank              !<
       INTEGER(iwp) ::  to_child_comm        !<
       INTEGER(iwp) ::  to_parent_comm       !<
       INTEGER(iwp) ::  rem_size             !<
       INTEGER(iwp) ::  rem_out_win = -1     !<
       INTEGER(iwp) ::  buf_win = -1         !<
       INTEGER(iwp) ::  nr_rem_out_particle  !< total number of remote output particle on this PE
       INTEGER(iwp), POINTER, DIMENSION(:) ::  index_and_size_buffer  !<
       TYPE(particle_type), ALLOCATABLE, DIMENSION(:) ::  transfer_buffer  !< buffer for output particles to be transferred to the root model
    END TYPE communicators_for_nested_pts_output

    TYPE(communicators_for_nested_pts_output), PUBLIC ::  cnpo  !< communicator infos for output of particle time series in nested runs


    INTERFACE pmcp_c_get_particle_from_parent
       MODULE PROCEDURE pmcp_c_get_particle_from_parent
    END INTERFACE pmcp_c_get_particle_from_parent

    INTERFACE pmcp_c_send_particle_to_parent
       MODULE PROCEDURE pmcp_c_send_particle_to_parent
    END INTERFACE pmcp_c_send_particle_to_parent

    INTERFACE pmcp_g_alloc_win
       MODULE PROCEDURE pmcp_g_alloc_win
    END INTERFACE pmcp_g_alloc_win

    INTERFACE pmcp_p_delete_particles_in_fine_grid_area
       MODULE PROCEDURE pmcp_p_delete_particles_in_fine_grid_area
    END INTERFACE pmcp_p_delete_particles_in_fine_grid_area

    INTERFACE pmcp_g_free_win
       MODULE PROCEDURE pmcp_g_free_win
    END INTERFACE pmcp_g_free_win

    INTERFACE pmcp_g_init
       MODULE PROCEDURE pmcp_g_init
    END INTERFACE pmcp_g_init

    INTERFACE pmcp_g_nested_init
       MODULE PROCEDURE pmcp_g_nested_init
    END INTERFACE pmcp_g_nested_init

    INTERFACE pmcp_g_print_number_of_particles
       MODULE PROCEDURE pmcp_g_print_number_of_particles
    END INTERFACE pmcp_g_print_number_of_particles

    INTERFACE pmcp_p_empty_particle_win
       MODULE PROCEDURE pmcp_p_empty_particle_win
    END INTERFACE pmcp_p_empty_particle_win

    INTERFACE pmcp_p_fill_particle_win
       MODULE PROCEDURE pmcp_p_fill_particle_win
    END INTERFACE pmcp_p_fill_particle_win

    PUBLIC pmcp_c_get_particle_from_parent
    PUBLIC pmcp_c_send_particle_to_parent
    PUBLIC pmcp_g_alloc_win
    PUBLIC pmcp_p_delete_particles_in_fine_grid_area
    PUBLIC pmcp_g_free_win
    PUBLIC pmcp_g_init
    PUBLIC pmcp_g_nested_init
    PUBLIC pmcp_g_print_number_of_particles
    PUBLIC pmcp_p_empty_particle_win
    PUBLIC pmcp_p_fill_particle_win

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> General routine:
!> Initializing actions of the particle interface check particle boundary conditions for the child
!> models.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_g_init

    IMPLICIT NONE

    INTEGER(iwp) ::  nr_childs  !< Number of child models of the current model

#if defined( __parallel )

    nr_childs = get_number_of_children()
!
!-- Check if the current model has child models
    IF ( nr_childs > 0 )  THEN
       ALLOCATE( nr_part(nysg:nyng, nxlg:nxrg) )
       ALLOCATE( part_adr(nysg:nyng, nxlg:nxrg) )
       nr_part  = 0
       part_adr = 0
    ENDIF
!
!-- Set the boundary conditions to nested for all non root (i.e child) models
    IF ( cpl_id > 1 )  THEN
       IF ( particle_coupling )  THEN
          ibc_par_t  = 3
          ibc_par_lr = 3
          ibc_par_ns = 3
       ENDIF

    ENDIF

#endif
 END SUBROUTINE pmcp_g_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize MPI nesting setup for particle output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_g_nested_init

    INTEGER(iwp) ::  color       !<
    INTEGER(iwp) ::  comm        !<
    INTEGER(iwp) ::  ierr        !<
    INTEGER(iwp) ::  nest_npes   !<
    INTEGER(iwp) ::  nest_rank   !<
    INTEGER(iwp) ::  peer_comm   !<
    INTEGER(iwp) ::  tag         !<
    INTEGER(iwp) ::  world_rank  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window


    CALL MPI_COMM_RANK( comm_world_nesting, world_rank, ierr )

    IF ( root_model )  THEN
       color = 1
    ELSE
       color = 2
    ENDIF

    CALL MPI_COMM_SPLIT( comm_world_nesting, color, world_rank, comm, ierr )
!
!-- Get size and rank of the model running on this PE.
    CALL  MPI_COMM_RANK( comm, nest_rank, ierr )
    CALL  MPI_COMM_SIZE( comm, nest_npes, ierr )

!
!-- Create intercommunicator between the root model and all children.
!-- MPI_INTERCOMM_CREATE creates an intercommunicator between 2 groups of different colors. The
!-- grouping was done above with MPI_COMM_SPLIT. A duplicate of comm_world_nesting is created and
!-- used as peer communicator (peer_comm) for MPI_INTERCOMM_CREATE.
    CALL MPI_COMM_DUP( comm_world_nesting, peer_comm, ierr )

    IF ( root_model )  THEN
!
!--    Identify all childs of the current model and create inter-communicator to connect between the
!--    root model and all childs.
       tag = 600
       CALL MPI_INTERCOMM_CREATE( comm, 0, peer_comm, numprocs, tag, cnpo%to_child_comm, ierr )
       CALL MPI_COMM_RANK( cnpo%to_child_comm, cnpo%my_rank, ierr )
       CALL MPI_COMM_SIZE( cnpo%to_child_comm, cnpo%my_npes, ierr )
       CALL MPI_COMM_REMOTE_SIZE( cnpo%to_child_comm, cnpo%rem_size, ierr )

       ALLOCATE( cnpo%index_and_size_buffer(cnpo%rem_size) )

       winsize = cnpo%rem_size * iwp
       CALL MPI_WIN_CREATE( cnpo%index_and_size_buffer, winsize, isp, MPI_INFO_NULL,               &
                            comm_world_nesting, cnpo%rem_out_win, ierr )

    ELSE
!
!--    Create an inter-communicator to connect between all child models and the root model.
       tag = 600
       CALL MPI_INTERCOMM_CREATE( comm, 0, peer_comm, 0, tag, cnpo%to_parent_comm, ierr )
       CALL MPI_COMM_RANK( cnpo%to_parent_comm, cnpo%my_rank, ierr )
       CALL MPI_COMM_SIZE( cnpo%to_parent_comm, cnpo%my_npes, ierr )
       CALL MPI_COMM_REMOTE_SIZE( cnpo%to_parent_comm, cnpo%rem_size, ierr )

       ALLOCATE( cnpo%index_and_size_buffer(1) )
!
!--    Dummy window on child processes.
       winsize = iwp
       CALL MPI_WIN_CREATE( cnpo%index_and_size_buffer, winsize, iwp, MPI_INFO_NULL,               &
                            comm_world_nesting, cnpo%rem_out_win, ierr )

    ENDIF

    CALL MPI_WIN_FENCE( 0, cnpo%rem_out_win, ierr )

 END SUBROUTINE pmcp_g_nested_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate the MPI windows used for particle output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_g_alloc_win

    IMPLICIT NONE

    INTEGER(iwp) ::  child_id   !< Id of a child model
    INTEGER(iwp) ::  ierr       !< error code
    INTEGER(iwp) ::  ipl        !< left boundary in coarse(parent) index space
    INTEGER(iwp) ::  ipr        !< right boundary in coarse(parent) index space
    INTEGER(iwp) ::  jps        !< south boundary in coarse(parent) index space
    INTEGER(iwp) ::  jpn        !< north boundary in coarse(parent) index space
    INTEGER(iwp) ::  m          !< loop index
    INTEGER(iwp) ::  nr_childs  !< Number of child models of the current model

    INTEGER ::  parsize  !<

    INTEGER(iwp),DIMENSION(1) ::  buf_shape  !<

    TYPE(C_PTR), SAVE ::  ptr  !<

    TYPE(particle_type),DIMENSION(:),POINTER ::  win_buffer  !<

#if defined( __parallel )
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  parsize_mpi_address_kind  !<
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize                   !<


!
!-- If the model has a parent model prepare the structures for transfer.
    IF ( cpl_id > 1 )  THEN

       parsize_mpi_address_kind = STORAGE_SIZE(zero_particle) / 8

       CALL MPI_ALLOC_MEM( parsize_mpi_address_kind , MPI_INFO_NULL, ptr, ierr )
       parsize = parsize_mpi_address_kind
       buf_shape(1) = 1
       CALL C_F_POINTER( ptr, win_buffer, buf_shape )
       CALL MPI_WIN_CREATE( win_buffer, parsize_mpi_address_kind, parsize, MPI_INFO_NULL,          &
                            me%intra_comm, particle_win_child, ierr )

!
!--    Child domain boundaries in the parent index space.
       ipl = parent_bound(1)
       ipr = parent_bound(2)
       jps = parent_bound(3)
       jpn = parent_bound(4)

       ALLOCATE( coarse_particles(jps:jpn,ipl:ipr) )

       coarse_particles(:,:)%nr_particle = 0

    ENDIF

!
!-- If the model has child models prepare the structures for transfer.
    nr_childs = get_number_of_children()
    IF ( nr_childs > 0 )   THEN
       ALLOCATE( particle_win_parent(nr_childs) )
       ALLOCATE( buf_ptr(nr_childs) )
       DO  m = 1, nr_childs
          child_id = get_childid(m)
          parsize_mpi_address_kind = STORAGE_SIZE(zero_particle) / 8
          parsize = parsize_mpi_address_kind

          winsize = max_nr_prts_in_rma_win * parsize_mpi_address_kind
          CALL MPI_ALLOC_MEM( winsize , MPI_INFO_NULL, buf_ptr(m), ierr )
          buf_shape(1) = max_nr_prts_in_rma_win
          CALL C_F_POINTER( buf_ptr(m), win_buffer, buf_shape )
          CALL MPI_WIN_CREATE( win_buffer, winsize, parsize, MPI_INFO_NULL,                        &
                               children(child_id)%intra_comm, particle_win_parent(m), ierr )
          ENDDO
    ENDIF

#endif
 END SUBROUTINE pmcp_g_alloc_win


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Free open MPI windows used for particle exchange.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_g_free_win

    IMPLICIT NONE

    INTEGER(iwp) ::  ierr       !< error code
    INTEGER(iwp) ::  m          !< loop index
    INTEGER(iwp) ::  nr_childs  !< Number of child models of the current model


    IF ( cpl_id > 1 )  THEN
       CALL MPI_WIN_FREE( particle_win_child, ierr )
    ENDIF

    nr_childs = get_number_of_children()
    IF ( nr_childs > 0 )   THEN
       DO  m = 1, nr_childs
          CALL MPI_WIN_FREE( particle_win_parent(m), ierr )
       ENDDO
    ENDIF

 END SUBROUTINE pmcp_g_free_win


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read/get particles out of the parent MPI window. Routine is used on childs only.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_c_get_particle_from_parent

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< x grid index
    INTEGER(iwp) ::  ierr  !< error code
    INTEGER(iwp) ::  ij    !< combined xy index for the buffer array
    INTEGER(iwp) ::  ip    !< loop index (child PEs)
    INTEGER(iwp) ::  ipl   !< left boundary in coarse(parent) index space
    INTEGER(iwp) ::  j     !< y grid index
    INTEGER(iwp) ::  jps   !< south boundary in coarse(parent) index space
    INTEGER(iwp) ::  nr    !< number of particles to receive from a parent box

    INTEGER ::  parsize !<

#if defined( __parallel )
    TYPE(pedef), POINTER ::  ape  !< TO_DO Klaus: give a description and better name of the variable

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  parsize_mpi_address_kind  !<
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_disp               !<


    IF ( cpl_id > 1 )  THEN
!
!--    Get number of particle/column and offset in RMA window.
!--    The parent side (in pmc_s_fillbuffer) is filling the buffer in the MPI RMA window. When the
!--    filling is complete, a MPI_BARRIER is called. The child is not allowd to access the
!--    parent-buffer before it is completely filled. Synchronization is done implicitely in
!--    pmc_c_getbuffer and pmc_s_fillbuffer on the parent side.
       CALL pmc_c_getbuffer( particle_transfer = .TRUE. )

       ipl = parent_bound(1)
       jps = parent_bound(3)

       DO  ip = 1, me%inter_npes

          ape => me%pes(ip)

          DO  ij = 1, ape%nrele
              j  = ape%locind(ij)%j + jps - 1
              i  = ape%locind(ij)%i + ipl - 1
              nr = nr_partc(j,i)
              IF ( nr > 0 )  THEN

                 CALL check_and_alloc_coarse_particle( i, j, nr )
                 parsize_mpi_address_kind = STORAGE_SIZE(zero_particle) / 8
                 parsize = parsize_mpi_address_kind
                 target_disp = part_adrc(j,i) - 1
                 CALL MPI_WIN_LOCK( MPI_LOCK_SHARED , ip - 1, 0, particle_win_child, ierr )
                 CALL MPI_GET( coarse_particles(j,i)%parent_particles, nr * parsize, MPI_BYTE,     &
                               ip - 1, target_disp, nr * parsize, MPI_BYTE, particle_win_child,    &
                               ierr )
                 CALL MPI_WIN_UNLOCK( ip - 1, particle_win_child, ierr )
             ENDIF
             coarse_particles(j,i)%nr_particle = nr
          ENDDO
       ENDDO

       CALL c_copy_particle_to_child_grid
    ENDIF

#endif
 END SUBROUTINE pmcp_c_get_particle_from_parent


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write/put particles into the parent MPI window. Routine is used on childs only.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_c_send_particle_to_parent

    IMPLICIT NONE

    INTEGER(iwp) ::  disp_offset             !<
    INTEGER(iwp) ::  i                       !< x loop index
    INTEGER(iwp) ::  ierr                    !< error code
    INTEGER(iwp) ::  ij                      !< combined xy index for the buffer array
    INTEGER(iwp) ::  ip                      !< loop index (child PEs)
    INTEGER(iwp) ::  ipl                     !< left boundary in coarse(parent) index space
    INTEGER(iwp) ::  ipr                     !< right boundary in coarse(parent) index space
    INTEGER(iwp) ::  j                       !< y loop index
    INTEGER(iwp) ::  jpn                     !< north boundary in coarse(parent) index space
    INTEGER(iwp) ::  jps                     !< south boundary in coarse(parent) index space
    INTEGER(iwp) ::  max_nr_particle_per_pe  !< maximum number of particles per PE (depending on grid apect ratio)
    INTEGER(iwp) ::  pe_offset               !< offset index of the current PE

    INTEGER ::  parsize  !<

    REAL(wp) ::  eps=0.00001  !< used in calculations to avoid rounding errors
    REAL(wp) ::  xx           !< number of fine grid cells inside a coarse grid cell in x-direction
    REAL(wp) ::  yy           !< number of fine grid cells inside a coarse grid cell in y-direction

 !   TYPE(particle_type) ::  dummy_part !< dummy particle (needed for size calculations)

#if defined( __parallel )
    TYPE(pedef), POINTER ::  ape  !< TO_DO Klaus: give a description and better name of the variable

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  parsize_mpi_address_kind  !<
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_disp               !<


    IF ( cpl_id > 1 )  THEN
       CALL c_copy_particle_to_coarse_grid

!
!--    Child domain boundaries in the parent index space.
       ipl = parent_bound(1)
       ipr = parent_bound(2)
       jps = parent_bound(3)
       jpn = parent_bound(4)

       nr_partc = 0

       DO i = ipl, ipr
          DO j = jps, jpn
             nr_partc(j,i) = coarse_particles(j,i)%nr_particle
          ENDDO
       ENDDO
       part_adrc = 0

!
!--    Compute number of fine grid cells in coarse grid (one direction). eps is added to avoid
!--    rounding errors.
       xx = ( pg%dx + eps ) / dx
       yy = ( pg%dy + eps ) / dy
       nr_fine_in_coarse = MAX( INT( xx ), INT( yy ) )

       IF ( MOD( coord_x(0), pg%dx ) /= 0.0  .OR.  MOD( coord_y(0), pg%dy ) /= 0.0 )  THEN
          nr_fine_in_coarse = nr_fine_in_coarse + 1
       ENDIF

!
!--    Assign a number to my child PE to select different areas in the RMA window on server side
!--    With this number a square of child PEs is defined which share the same coarse grid cells.
       pe_offset   = MOD( myidx, nr_fine_in_coarse ) * nr_fine_in_coarse +                         &
                     MOD( myidy, nr_fine_in_coarse )
       max_nr_particle_per_pe = max_nr_prts_in_rma_win / ( nr_fine_in_coarse * nr_fine_in_coarse )
       disp_offset            = pe_offset * max_nr_particle_per_pe
       parsize_mpi_address_kind = STORAGE_SIZE( zero_particle ) / 8
       parsize = parsize_mpi_address_kind
       DO  ip = 1, me%inter_npes

          ape => me%pes(ip)
          target_disp = disp_offset
          DO  ij = 1, ape%nrele

             j  = ape%locind(ij)%j + jps - 1
             i  = ape%locind(ij)%i + ipl - 1

             IF ( nr_partc(j,i) > 0 )  THEN

                IF ( target_disp + nr_partc(j,i) - disp_offset >= max_nr_particle_per_pe )  THEN
                   IF ( debug_output )  THEN
                      WRITE( debug_string, '(A,3(1X,I10))' )  'RMA window too small on child ',    &
                             target_disp + nr_partc(j,i) - disp_offset, max_nr_particle_per_pe,    &
                             max_nr_prts_in_rma_win
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   message_string = 'RMA window too small on child'
                   CALL message( 'pmcp_c_send_particle_to_parent', 'PMC0028', 3, 2, 0, 6, 0 )
                ENDIF

                CALL MPI_WIN_LOCK( MPI_LOCK_SHARED , ip - 1, 0, particle_win_child, ierr )
                CALL MPI_PUT( coarse_particles(j,i)%parent_particles,                              &
                              INT( nr_partc(j,i), iwp ) * parsize, MPI_BYTE, ip-1, target_disp,    &
                              INT( nr_partc(j,i), iwp ) * parsize, MPI_BYTE,                       &
                              particle_win_child, ierr )
                CALL MPI_WIN_UNLOCK( ip - 1, particle_win_child, ierr )
                part_adrc(j,i) = target_disp + 1
                target_disp    = target_disp + nr_partc(j,i)

             ENDIF

          ENDDO

       ENDDO
!
!--    Send new number of particles/column and offset to parent.
       CALL pmc_c_putbuffer( particle_transfer = .TRUE. )

    ENDIF

#endif
 END SUBROUTINE pmcp_c_send_particle_to_parent


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write particles into the MPI window. Routine is used on parents only.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_p_fill_particle_win

    IMPLICIT NONE

    INTEGER(iwp) ::  child_id             !< id of the child model
    INTEGER(iwp) ::  i                    !< x grid box index
    INTEGER(iwp) ::  ij                   !< combined xy index for the buffer array
    INTEGER(iwp) ::  ip                   !< loop index (child PEs)
    INTEGER(iwp) ::  j                    !< y grid box index
    INTEGER(iwp) ::  k                    !< z grid box index
    INTEGER(iwp) ::  m                    !< loop index (number of childs)
    INTEGER(iwp) ::  n                    !< loop index (number of particles)
    INTEGER(iwp) ::  nr_part_col          !< Number of particles to transfer per column
    INTEGER(iwp) ::  number_of_particles  !<
    INTEGER(iwp) ::  pindex               !<
    INTEGER(iwp) ::  tot_particle_count   !< Total number of particles per child

    INTEGER(iwp),DIMENSION(1) ::  buf_shape  !<

    LOGICAL       ::  do_transfer      !< flag to indicate that particle is located in the fine/child grid and requires transfer to parent
    LOGICAL, SAVE ::  lfirst = .TRUE.  !<

    REAL(wp) ::  dx_child    !< child grid spacing
    REAL(wp) ::  dy_child    !< child grid spacing
    REAL(wp) ::  dz_child    !< child grid spacing
    REAL(wp) ::  ny_coord    !< north coordinate of child grid
    REAL(wp) ::  ny_coord_b  !< north coordinate of child grid boundary
    REAL(wp) ::  lx_coord    !< left coordinate of child grid
    REAL(wp) ::  lx_coord_b  !< left coordinate of child grid boundary
    REAL(wp) ::  rx_coord    !< right coordinate of child grid
    REAL(wp) ::  rx_coord_b  !< right coordinate of child grid boundary
    REAL(wp) ::  sy_coord    !< south coordinate of child grid
    REAL(wp) ::  sy_coord_b  !< south coordinate of child grid boundary
    REAL(wp) ::  uz_coord    !< upper coordinate of child grid
    REAL(wp) ::  uz_coord_b  !< upper coordinate of child grid boundary
    REAL(wp) ::  x           !< particle position
    REAL(wp) ::  y           !< particle position
    REAL(wp) ::  z           !< particle position

#if defined( __parallel )
    TYPE(pedef), POINTER ::  ape  !< TO_DO Klaus: give a description and better name of the variable


    DO  m = 1, get_number_of_children()

       child_id = get_childid(m)

       CALL get_child_edges( m, lx_coord, lx_coord_b, rx_coord, rx_coord_b, sy_coord, sy_coord_b,  &
                             ny_coord, ny_coord_b, uz_coord, uz_coord_b)

       CALL get_child_gridspacing( m, dx_child, dy_child, dz_child )

       IF ( lfirst  .AND.  debug_output )  THEN
          WRITE( debug_string, '(3(A,5(1X,F10.2)/))' )                                             &
             'edges          ', lx_coord,rx_coord, sy_coord, ny_coord, uz_coord,                   &
             'edges boundary ', lx_coord_b, rx_coord_b, sy_coord_b, ny_coord_b, uz_coord_b,        &
             'child spacing  ', dx_child, dy_child, dz_child, lower_left_coord_x, lower_left_coord_y
          CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Reset values for every child.
       tot_particle_count = 0
       nr_part  = 0
       part_adr = 0
       pindex   = 1

       buf_shape(1) = max_nr_prts_in_rma_win
       CALL C_F_POINTER( buf_ptr(m), particle_in_win , buf_shape )

       DO  ip = 1, children(child_id)%inter_npes

          ape => children(child_id)%pes(ip)

          nr_part_col   = 0

          DO  ij = 1, ape%nrele

!
!--          Inside the PMC adressing of 3d arrays starts with 1.
             i = ape%locind(ij)%i + nxl - nbgp - 1
             j = ape%locind(ij)%j + nys - nbgp - 1
             nr_part_col = 0   ! Number of particles to transfer per column
             part_adr(j,i) = pindex

             DO  k = nzb + 1, nzt
                number_of_particles = prt_count(k,j,i)

                IF ( number_of_particles <= 0 )  CYCLE

                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
!
!--             Select particles within boundary area.
                DO  n = 1, number_of_particles
                   x = particles(n)%x + lower_left_coord_x
                   y = particles(n)%y + lower_left_coord_y
                   z = particles(n)%z
!
!--                Check if the particle is located in the fine grid area and therefore requires
!--                to be transferred to the parent.
                   do_transfer = ( ( x > lx_coord  .AND.  x < rx_coord )  .AND.                    &
                                   ( y > sy_coord  .AND.  y < ny_coord )  .AND.                    &
                                   ( z > 0.000000001  .AND.  z < uz_coord )                        &
                                 )

                   IF ( do_transfer .AND. particles(n)%particle_mask )  THEN

                      particle_in_win(pindex) = particles(n)
!
!--                   Change particle positions and origin relative to global origin.
                      particle_in_win(pindex)%x = particle_in_win(pindex)%x + lower_left_coord_x
                      particle_in_win(pindex)%y = particle_in_win(pindex)%y + lower_left_coord_y
                      particle_in_win(pindex)%origin_x = particle_in_win(pindex)%origin_x          &
                                                         + lower_left_coord_x
                      particle_in_win(pindex)%origin_y = particle_in_win(pindex)%origin_y          &
                                                         + lower_left_coord_y

                      tot_particle_count = tot_particle_count + 1
                      nr_part_col        = nr_part_col + 1
                      pindex             = pindex + 1
                      IF ( pindex > max_nr_prts_in_rma_win )  THEN
                         IF ( debug_output )  THEN
                            WRITE( debug_string, '(A,2(1X,I10))' )                                 &
                                   'RMA window too small on parent ', pindex, max_nr_prts_in_rma_win
                            CALL debug_message( debug_string, 'info' )
                         ENDIF
                         message_string = 'RMA window too small on parent'
                         CALL message( 'pmcp_p_fill_particle_win', 'PMC0028', 3, 2, 0, 6, 0 )
                     ENDIF

                   ENDIF
                ENDDO
             ENDDO
             nr_part(j,i) = nr_part_col
          ENDDO
       ENDDO

       CALL pmc_s_fillbuffer( child_id, particle_transfer = .TRUE. )

    ENDDO

    lfirst = .FALSE.

#endif
 END SUBROUTINE pmcp_p_fill_particle_win


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Delete particles from the MPI window. Routine is used on parents only.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_p_empty_particle_win

    IMPLICIT NONE

    INTEGER(iwp) ::  child_id  !< model id of the child
    INTEGER(iwp) ::  ip        !< loop index (child PEs)
    INTEGER(iwp) ::  m         !< loop index (number of childs)

    INTEGER(iwp),DIMENSION(1) ::  buf_shape  !<

#if defined( __parallel )
    DO  m = 1, get_number_of_children()

       child_id = get_childid(m)

       buf_shape(1) = max_nr_prts_in_rma_win
       CALL C_F_POINTER( buf_ptr(m), particle_in_win , buf_shape )

!
!--    In some cells of the coarse grid, there are contributions from more than one child process.
!--    Therefore p_copy_particle_to_org_grid is done for one child process per call.
       DO  ip = 1, pmc_s_get_child_npes( child_id )

          nr_part  = 0
          part_adr = 0

          CALL pmc_s_getdata_from_buffer( child_id, particle_transfer = .TRUE.,                    &
                                          child_process_nr = ip )
          CALL p_copy_particle_to_org_grid( m )
       ENDDO

    ENDDO

#endif
 END SUBROUTINE pmcp_p_empty_particle_win


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> After the transfer mark all parent particles that are still inside on of the child areas for
!> deletion. Routine is used on parents only.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_p_delete_particles_in_fine_grid_area

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< loop index (x grid)
    INTEGER(iwp) ::  j  !< loop index (y grid)
    INTEGER(iwp) ::  k  !< loop index (z grid)
    INTEGER(iwp) ::  m  !< loop index (number of particles)
    INTEGER(iwp) ::  n  !< loop index (number of childs)

    LOGICAL ::  to_delete  !< particles outside of model domain are marked as to_delete

    REAL(wp) ::  dx_child    !< child grid spacing
    REAL(wp) ::  dy_child    !< child grid spacing
    REAL(wp) ::  dz_child    !< child grid spacing
    REAL(wp) ::  ny_coord    !< north coordinate of child grid
    REAL(wp) ::  ny_coord_b  !< north coordinate of child grid boundary
    REAL(wp) ::  lx_coord    !< left coordinate of child grid
    REAL(wp) ::  lx_coord_b  !< left coordinate of child grid boundary
    REAL(wp) ::  rx_coord    !< right coordinate of child grid
    REAL(wp) ::  rx_coord_b  !< right coordinate of child grid boundary
    REAL(wp) ::  sy_coord    !< south coordinate of child grid
    REAL(wp) ::  sy_coord_b  !< south coordinate of child grid boundary
    REAL(wp) ::  uz_coord    !< upper coordinate of child grid
    REAL(wp) ::  uz_coord_b  !< upper coordinate of child grid boundary
    REAL(wp) ::  x           !< particle position
    REAL(wp) ::  y           !< particle position
    REAL(wp) ::  z           !< particle position

#if defined( __parallel )
    DO  m = 1, get_number_of_children()
       CALL get_child_edges( m, lx_coord, lx_coord_b, rx_coord, rx_coord_b, sy_coord, sy_coord_b,  &
                             ny_coord, ny_coord_b, uz_coord, uz_coord_b )


       CALL get_child_gridspacing( m, dx_child, dy_child, dz_child )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
                number_of_particles = prt_count(k,j,i)

                IF ( number_of_particles == 0 )  CYCLE

                particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                DO  n = 1, number_of_particles
                   x = particles(n)%x + lower_left_coord_x
                   y = particles(n)%y + lower_left_coord_y
                   z = particles(n)%z

                   to_delete = ( (x > lx_coord .AND. x < rx_coord) .AND.                           &
                                 (y > sy_coord .AND. y < ny_coord) .AND.                           &
                                 (z > 0.000000001 .AND. z < uz_coord) )

                   IF ( to_delete )  THEN
                      particles(n)%particle_mask = .FALSE.
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDDO

#endif
 END SUBROUTINE pmcp_p_delete_particles_in_fine_grid_area


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print the total number of particles of the current model and its childs into a file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmcp_g_print_number_of_particles( my_time, local_nr_particles )

    USE pegrid,                                                                                    &
        ONLY:  myid

    IMPLICIT NONE

    CHARACTER(LEN=16) ::  fname  !< filename

    INTEGER(iwp) ::  child_id            !<
    INTEGER(iwp) ::  child_nr_particles  !< total number of particles in all child models
    INTEGER(iwp) ::  ierr                !< error code
    INTEGER(iwp) ::  m                   !< loop index (number of childs

    INTEGER(iwp),INTENT(IN) ::  local_nr_particles  !<

    INTEGER(iwp),DIMENSION(2) ::  ivalr  !< integer value to be received
    INTEGER(iwp),DIMENSION(2) ::  ivals  !< integer value to be send

    LOGICAL, SAVE ::  is_file_open=.FALSE.  !<

    REAL(wp),INTENT(IN) ::  my_time  !<


#if defined( __parallel )
    child_nr_particles = 0
    IF ( myid == 0 )  THEN
       IF ( cpl_id > 1 )  THEN
          ivals(1) = local_nr_particles
          CALL pmc_send_to_parent( ivals, 1, 0, 400, ierr )
       ENDIF
       DO  m = 1, SIZE( pmc_parent_for_child ) - 1

          child_id = pmc_parent_for_child(m)
          CALL pmc_recv_from_child( child_id, ivalr, 1, 0, 400, ierr )
          child_nr_particles = child_nr_particles + ivalr(1)
       ENDDO

       IF ( SIZE( pmc_parent_for_child ) > 1 )  THEN
          IF ( .NOT. is_file_open )  THEN  !kk muss noch auf file_open umgestellt werden
             WRITE( fname, '(A,I2.2)' ) 'nr_particles_', cpl_id
             OPEN ( 333, FILE = fname )
             is_file_open = .TRUE.
          ENDIF
          WRITE( 333, '(F12.2,3I10)' ) my_time, local_nr_particles + child_nr_particles,           &
                                       local_nr_particles, child_nr_particles
       ENDIF
    ENDIF

#endif
 END SUBROUTINE pmcp_g_print_number_of_particles


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Update the size of the local buffer (coarse_particles).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE check_and_alloc_coarse_particle( ic, jc, nr, with_copy )

    IMPLICIT NONE

    INTEGER(iwp),INTENT(IN) ::  ic  !< coarse x grid index
    INTEGER(iwp),INTENT(IN) ::  jc  !< coarse y grid index
    INTEGER(iwp),INTENT(IN) ::  nr  !< number of particles to be transferred/stored in local buffer

    INTEGER(iwp) ::  new_size  !< new size of the local buffer
    INTEGER(iwp) ::  old_size  !< old size of the local buffer

    LOGICAL :: with_copy_lo  !< local variable of with copy

    LOGICAL,INTENT(IN),OPTIONAL ::  with_copy  !< copy particles in buffer? or reallocate empty buffer

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  tmp_particles_d  !<

#if defined( __parallel )
    with_copy_lo = .FALSE.
    IF ( PRESENT( with_copy ) )  with_copy_lo = with_copy

    IF ( .NOT. ALLOCATED( coarse_particles(jc,ic)%parent_particles ) )  THEN
       new_size = MAX( nr, min_particles_per_column )
       ALLOCATE( coarse_particles(jc,ic)%parent_particles(new_size) )
    ELSEIF ( nr > SIZE( coarse_particles(jc,ic)%parent_particles ) )  THEN

       old_size = SIZE( coarse_particles(jc,ic)%parent_particles )
       new_size = old_size * ( 1.0_wp + alloc_factor / 100.0_wp )
       new_size = MAX( nr, new_size, old_size + min_particles_per_column )

!
!--    Copy existing particles to new particle buffer
       IF ( with_copy_lo )  THEN
          ALLOCATE( tmp_particles_d(old_size) )
          tmp_particles_d(1:old_size) = coarse_particles(jc,ic)%parent_particles

          DEALLOCATE( coarse_particles(jc,ic)%parent_particles )
          ALLOCATE( coarse_particles(jc,ic)%parent_particles(new_size) )

          coarse_particles(jc,ic)%parent_particles(1:old_size)          = tmp_particles_d(1:old_size)
          coarse_particles(jc,ic)%parent_particles(old_size+1:new_size) = zero_particle

          DEALLOCATE( tmp_particles_d )
!
!--    allocate or reallocate an empty buffer
       ELSE
          DEALLOCATE( coarse_particles(jc,ic)%parent_particles )
          ALLOCATE( coarse_particles(jc,ic)%parent_particles(new_size) )
       ENDIF
    ENDIF

#endif
 END SUBROUTINE check_and_alloc_coarse_particle


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Copy/sort particles out of the local buffer into the respective grid boxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE c_copy_particle_to_child_grid

    IMPLICIT NONE

    INTEGER(iwp) ::  ic   !< coarse x grid index
    INTEGER(iwp) ::  ipl  !< left boundary in coarse(parent) index space
    INTEGER(iwp) ::  ipr  !< right boundary in coarse(parent) index space
    INTEGER(iwp) ::  ip   !< x grid index
    INTEGER(iwp) ::  jc   !< coarse y grid index
    INTEGER(iwp) ::  jpn  !< north boundary in coarse(parent) index space
    INTEGER(iwp) ::  jps  !< south boundary in coarse(parent) index space
    INTEGER(iwp) ::  jp   !< y grid index
    INTEGER(iwp) ::  kp   !< z grid index
    INTEGER(iwp) ::  n    !< loop index (number of particles)
    INTEGER(iwp) ::  nr   !< short variable for name number or particles

    REAL(wp) ::  xc   !< child x coordinate
    REAL(wp) ::  xoc  !< child x origin
    REAL(wp) ::  yc   !< child y coordinate
    REAL(wp) ::  yoc  !< child y origin
    REAL(wp) ::  zc   !< child z coordinate

#if defined( __parallel )
!
!-- Child domain boundaries in the parent index space
    ipl = parent_bound(1)
    ipr = parent_bound(2)
    jps = parent_bound(3)
    jpn = parent_bound(4)

    DO  ic = ipl, ipr
       DO  jc = jps, jpn
          nr = coarse_particles(jc,ic)%nr_particle

          IF ( nr > 0 )  THEN
             DO  n = 1, nr
                xc =  coarse_particles(jc,ic)%parent_particles(n)%x-lower_left_coord_x
                yc =  coarse_particles(jc,ic)%parent_particles(n)%y-lower_left_coord_y
                zc =  coarse_particles(jc,ic)%parent_particles(n)%z
                xoc = coarse_particles(jc,ic)%parent_particles(n)%origin_x-lower_left_coord_x
                yoc = coarse_particles(jc,ic)%parent_particles(n)%origin_y-lower_left_coord_y
                ip = xc / dx
                jp = yc / dy
                kp = nzt
                DO WHILE ( zw(kp-1) > zc .AND. kp > nzb + 1 )  ! kk search loop has to be optimzed !!!
                   kp = kp - 1
                ENDDO

                IF ( ip >= nxl  .AND.  ip <= nxr  .AND.  jp >= nys  .AND.  jp <= nyn )  THEN
                   prt_count(kp,jp,ip) = prt_count(kp,jp,ip) + 1
                   IF ( prt_count(kp,jp,ip) > SIZE( grid_particles(kp,jp,ip)%particles ) )  THEN
                      CALL pmc_realloc_particles_array( ip, jp, kp )
                   ENDIF
                   coarse_particles(jc,ic)%parent_particles(n)%x = xc         ! Adjust coordinates to child grid
                   coarse_particles(jc,ic)%parent_particles(n)%y = yc
                   coarse_particles(jc,ic)%parent_particles(n)%origin_x = xoc ! Adjust origins to child grid
                   coarse_particles(jc,ic)%parent_particles(n)%origin_y = yoc
                   grid_particles(kp,jp,ip)%particles(prt_count(kp,jp,ip))                         &
                                                       = coarse_particles(jc,ic)%parent_particles(n)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO

#endif
 END SUBROUTINE c_copy_particle_to_child_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Copy particles which left the model area into the local buffer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE c_copy_particle_to_coarse_grid

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< loop index (x grid)
    INTEGER(iwp) ::  ic    !< loop index (coarse x grid)
    INTEGER(iwp) ::  ipl   !< left boundary in coarse(parent) index space
    INTEGER(iwp) ::  ipr   !< left boundary in coarse(parent) index space
    INTEGER(iwp) ::  j     !< loop index (y grid)
    INTEGER(iwp) ::  jc    !< loop index (coarse y grid)
    INTEGER(iwp) ::  jps   !< south boundary in coarse(parent) index space
    INTEGER(iwp) ::  jpn   !< north boundary in coarse(parent) index space
    INTEGER(iwp) ::  k     !< loop index (z grid)
    INTEGER(iwp) ::  n     !< loop index (number of particles)

    LOGICAL ::  boundary_particle  !<

    REAL(wp) ::  parent_lower_left_x  !< lower left x-coordinate of respective parent
    REAL(wp) ::  parent_lower_left_y  !< lower left y-coordinate of respective parent
    REAL(wp) ::  x                    !< x coordinate
    REAL(wp) ::  xo                   !< x origin
    REAL(wp) ::  x_left               !< absolute left boundary
    REAL(wp) ::  x_right              !< absolute right boundary
    REAL(wp) ::  y                    !< left boundary in coarse(parent) index space
    REAL(wp) ::  yo                   !< y origin
    REAL(wp) ::  y_north              !< absolute north boundary
    REAL(wp) ::  y_south              !< absoulte south boundary
    REAL(wp) ::  z                    !< z coordinate

#if defined( __parallel )
!
!-- Child domain boundaries in the parent index space
    ipl = parent_bound(1)
    ipr = parent_bound(2)
    jps = parent_bound(3)
    jpn = parent_bound(4)

!
!-- Absolute coordinates
    x_left  = coord_x(0)
    x_right = coord_x(nx) + dx
    y_south = coord_y(0)
    y_north = coord_y(ny) + dy
!
!-- Clear Particle Buffer
    DO  ic = ipl, ipr
       DO  jc = jps, jpn
          coarse_particles(jc,ic)%nr_particle = 0
       ENDDO
    ENDDO

    CALL pmc_get_model_info( parent_lower_left_x = parent_lower_left_x,                            &
                             parent_lower_left_y = parent_lower_left_y)

!
!-- Determine particles which leave the inner area in east or west dirextion.
!-- Compute only first (nxl) and last (nxr) loop iterration.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb + 1, nzt
             number_of_particles = prt_count(k,j,i)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
                IF ( particles(n)%particle_mask )  THEN
                   x  = particles(n)%x+ lower_left_coord_x
                   y  = particles(n)%y+ lower_left_coord_y
                   xo = particles(n)%origin_x + lower_left_coord_x
                   yo = particles(n)%origin_y + lower_left_coord_y
                   z  = particles(n)%z

                   boundary_particle = .FALSE.
                   boundary_particle = ( boundary_particle .OR. ( x < x_left  .OR. x > x_right ) )
                   boundary_particle = ( boundary_particle .OR. ( y < y_south .OR. y > y_north ) )
                   boundary_particle = ( boundary_particle .OR. ( z > zw(nzt) ) )

                   IF ( boundary_particle )  THEN
                      ic = ( x - parent_lower_left_x ) / pg%dx
                      jc = ( y - parent_lower_left_y ) / pg%dy
                      
                      IF ( ic >= ipl  .AND.  ic <= ipr  .AND.  jc >= jps  .AND.  jc <= jpn )  THEN
                         coarse_particles(jc,ic)%nr_particle = coarse_particles(jc,ic)%nr_particle &
                                                               + 1
                         CALL check_and_alloc_coarse_particle( ic, jc,                             &
                                                               coarse_particles(jc,ic)%nr_particle,&
                                                               with_copy = .TRUE. )

                         coarse_particles(jc,ic)%parent_particles(                                 &
                         coarse_particles(jc,ic)%nr_particle) = particles(n)
                         coarse_particles(jc,ic)%parent_particles(                                 &
                         coarse_particles(jc,ic)%nr_particle)%x = x   !Adapt to global coordinates
                         coarse_particles(jc,ic)%parent_particles(                                 &
                         coarse_particles(jc,ic)%nr_particle)%y = y
                         coarse_particles(jc,ic)%parent_particles(                                 &
                         coarse_particles(jc,ic)%nr_particle)%origin_x = xo
                         coarse_particles(jc,ic)%parent_particles(                                 &
                         coarse_particles(jc,ic)%nr_particle)%origin_y = yo
!
!--                      Mark particle as deleted after copying it to the transfer buffer
                         grid_particles(k,j,i)%particles(n)%particle_mask = .FALSE.
                      ELSE
!
!--                      For security reasons. Error shouldn't appear.
                         WRITE( message_string, '(A,10(1X,I6))' )                                  &
                              'particle outside parent bounds:', i, j, k, ic, jc, ipl, ipr, jps, jpn
                         CALL message( 'c_copy_particle_to_coarse_grid', 'PMC0029', 3, 2, 0, 6, 0 )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- Pack particles (eliminate those marked for deletion), determine new number of particles
!    CALL lpm_sort_in_subboxes

#endif
 END SUBROUTINE c_copy_particle_to_coarse_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Copy/sort particles from the MPI window into the respective grid boxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE p_copy_particle_to_org_grid( m )

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index (x grid)
    INTEGER(iwp) ::  j       !< loop index (y grid)
    INTEGER(iwp) ::  k       !< loop index (z grid)
    INTEGER(iwp) ::  n       !< loop index (nr part)
    INTEGER(iwp) ::  nr      !< short variable name for nr_part
    INTEGER(iwp) ::  pindex  !< short variable name part_adr

    INTEGER(iwp),INTENT(IN) ::  m  !<

    INTEGER(iwp),DIMENSION(1) ::  buf_shape  !<

    REAL(wp) ::  x   !< x coordinate
    REAL(wp) ::  xo  !< x origin
    REAL(wp) ::  y   !< y coordinate
    REAL(wp) ::  yo  !< y origin
    REAL(wp) ::  z   !< z coordinate


#if defined( __parallel )
    buf_shape(1) = max_nr_prts_in_rma_win
    CALL C_F_POINTER( buf_ptr(m), particle_in_win , buf_shape )

    DO  i = nxl, nxr
       DO  j = nys, nyn
          nr = nr_part(j,i)
          IF ( nr > 0 )  THEN
             pindex = part_adr(j,i)
             DO  n = 1, nr
                x = particle_in_win(pindex)%x-lower_left_coord_x
                y = particle_in_win(pindex)%y-lower_left_coord_y
                z = particle_in_win(pindex)%z
                xo = particle_in_win(pindex)%origin_x-lower_left_coord_x
                yo = particle_in_win(pindex)%origin_y-lower_left_coord_y
                k = nzt + 1
                DO WHILE ( zw(k-1) > z .AND. k > nzb + 1 )  ! kk search loop has to be optimzed !!!
                   k = k - 1
                END DO

                prt_count(k,j,i) = prt_count(k,j,i) + 1
                IF ( prt_count(k,j,i) > SIZE( grid_particles(k,j,i)%particles ) )  THEN
                   CALL pmc_realloc_particles_array( i, j, k )
                ENDIF
                grid_particles(k,j,i)%particles(prt_count(k,j,i)) = particle_in_win(pindex)

!
!--             Update particle positions and origins relative to parent domain
                grid_particles(k,j,i)%particles(prt_count(k,j,i))%x = x
                grid_particles(k,j,i)%particles(prt_count(k,j,i))%y = y
                grid_particles(k,j,i)%particles(prt_count(k,j,i))%origin_x = xo
                grid_particles(k,j,i)%particles(prt_count(k,j,i))%origin_y = yo
                pindex = pindex + 1
             ENDDO
          ENDIF
       ENDDO
    ENDDO

#endif
 END SUBROUTINE p_copy_particle_to_org_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> If the allocated memory for the particle array do not suffice to add arriving particles from
!> neighbour grid cells, this subroutine reallocates the particle array to assure enough memory is
!> available.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_realloc_particles_array( i, j, k, size_in )

    INTEGER(iwp) ::  new_size  !<
    INTEGER(iwp) ::  old_size  !<

    INTEGER(iwp), INTENT(IN) ::  i  !<
    INTEGER(iwp), INTENT(IN) ::  j  !<
    INTEGER(iwp), INTENT(IN) ::  k  !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  size_in  !<


    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  tmp_particles_d  !<

    TYPE(particle_type), DIMENSION(500)            ::  tmp_particles_s  !<

    old_size = SIZE( grid_particles(k,j,i)%particles )

    IF ( PRESENT( size_in ) )  THEN
       new_size = size_in
    ELSE
       new_size = old_size * ( 1.0_wp + alloc_factor / 100.0_wp )
    ENDIF

    new_size = MAX( new_size, 1, old_size + 1 )

    IF ( old_size <= 500 )  THEN

       tmp_particles_s(1:old_size) = grid_particles(k,j,i)%particles(1:old_size)

       DEALLOCATE( grid_particles(k,j,i)%particles )
       ALLOCATE( grid_particles(k,j,i)%particles(new_size) )

       grid_particles(k,j,i)%particles(1:old_size)          = tmp_particles_s(1:old_size)
       grid_particles(k,j,i)%particles(old_size+1:new_size) = zero_particle

    ELSE

       ALLOCATE( tmp_particles_d(new_size) )
       tmp_particles_d(1:old_size) = grid_particles(k,j,i)%particles

       DEALLOCATE( grid_particles(k,j,i)%particles )
       ALLOCATE( grid_particles(k,j,i)%particles(new_size) )

       grid_particles(k,j,i)%particles(1:old_size)          = tmp_particles_d(1:old_size)
       grid_particles(k,j,i)%particles(old_size+1:new_size) = zero_particle

       DEALLOCATE( tmp_particles_d )

    ENDIF
    particles => grid_particles(k,j,i)%particles(1:new_size)

    RETURN

 END SUBROUTINE pmc_realloc_particles_array

#endif
END MODULE pmc_particle_interface
