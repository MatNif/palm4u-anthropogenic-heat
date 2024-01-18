!> @file pmc_parent_mod.f90
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
! Authors:
! --------
!> @author Klaus Ketelsen (no affiliation)
!
! Description:
! ------------
!> Parent part of Palm Model Coupler
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_parent

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

    USE MPI

    USE kinds

    USE pmc_general,                                                                               &
        ONLY:  arraydef,                                                                           &
               childdef,                                                                           &
               da_namedef,                                                                         &
               da_namelen,                                                                         &
               pedef,                                                                              &
               pmc_g_setname,                                                                      &
               pmc_max_array,                                                                      &
               pmc_max_models

    USE pmc_handle_communicator,                                                                   &
        ONLY:  m_model_comm,                                                                       &
               m_model_rank,                                                                       &
               m_model_npes,                                                                       &
               m_to_child_comm,                                                                    &
               m_world_rank,                                                                       &
               pmc_parent_for_child

    USE pmc_mpi_wrapper,                                                                           &
        ONLY:  pmc_alloc_mem,                                                                      &
               pmc_bcast,                                                                          &
               pmc_time

    IMPLICIT NONE

    INTEGER ::  next_array_in_list = 0  !<

    TYPE childindexdef
       INTEGER ::  nrpoints  !<

       INTEGER, DIMENSION(:,:), ALLOCATABLE ::  index_list_2d  !<
    END TYPE childindexdef

    TYPE(childdef), DIMENSION(pmc_max_models) ::  children  !<

    TYPE(childindexdef), DIMENSION(pmc_max_models) ::  indchildren  !<

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC pmc_parent_for_child

!
!-- Public variables, constants and types
    PUBLIC children,                                                                               &
           pmc_parentinit,                                                                         &
           pmc_s_clear_next_array_list,                                                            &
           pmc_s_fillbuffer,                                                                       &
           pmc_s_finalize,                                                                         &
           pmc_s_getdata_from_buffer,                                                              &
           pmc_s_getnextarray,                                                                     &
           pmc_s_setind_and_allocmem,                                                              &
           pmc_s_set_active_data_array,                                                            &
           pmc_s_set_dataarray,                                                                    &
           pmc_s_set_2d_index_list,                                                                &
           pmc_s_get_child_npes

    INTERFACE pmc_parentinit
       MODULE PROCEDURE  pmc_parentinit
    END INTERFACE pmc_parentinit

    INTERFACE pmc_s_clear_next_array_list
       MODULE PROCEDURE pmc_s_clear_next_array_list
    END INTERFACE pmc_s_clear_next_array_list

    INTERFACE pmc_s_fillbuffer
       MODULE PROCEDURE pmc_s_fillbuffer
    END INTERFACE pmc_s_fillbuffer

    INTERFACE pmc_s_finalize
       MODULE PROCEDURE  pmc_s_finalize
    END INTERFACE pmc_s_finalize

    INTERFACE pmc_s_getdata_from_buffer
       MODULE PROCEDURE pmc_s_getdata_from_buffer
    END INTERFACE pmc_s_getdata_from_buffer

    INTERFACE pmc_s_getnextarray
       MODULE PROCEDURE pmc_s_getnextarray
    END INTERFACE pmc_s_getnextarray

    INTERFACE pmc_s_get_child_npes
       MODULE PROCEDURE pmc_s_get_child_npes
    END INTERFACE pmc_s_get_child_npes

    INTERFACE pmc_s_setind_and_allocmem
       MODULE PROCEDURE pmc_s_setind_and_allocmem
    END INTERFACE pmc_s_setind_and_allocmem

    INTERFACE pmc_s_set_active_data_array
       MODULE PROCEDURE pmc_s_set_active_data_array
    END INTERFACE pmc_s_set_active_data_array

    INTERFACE pmc_s_set_dataarray
       MODULE PROCEDURE pmc_s_set_dataarray_2d
       MODULE PROCEDURE pmc_s_set_dataarray_3d
       MODULE PROCEDURE pmc_s_set_dataarray_ip2d
    END INTERFACE pmc_s_set_dataarray

    INTERFACE pmc_s_set_2d_index_list
       MODULE PROCEDURE pmc_s_set_2d_index_list
    END INTERFACE pmc_s_set_2d_index_list

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> If this thread is intended as parent, initialize parent part of parent-client data transfer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_parentinit

    INTEGER(iwp) ::  childid  !<
    INTEGER(iwp) ::  i        !<
    INTEGER(iwp) ::  istat    !<
    INTEGER(iwp) ::  j        !<


    DO  i = 1, SIZE( pmc_parent_for_child ) - 1
       childid = pmc_parent_for_child( i )

       children(childid)%model_comm = m_model_comm
       children(childid)%inter_comm = m_to_child_comm(childid)

!
!--    Get rank and size
       CALL MPI_COMM_RANK( children(childid)%model_comm, children(childid)%model_rank, istat )
       CALL MPI_COMM_SIZE( children(childid)%model_comm, children(childid)%model_npes, istat )
       CALL MPI_COMM_REMOTE_SIZE( children(childid)%inter_comm, children(childid)%inter_npes,      &
                                  istat )

!
!--    Intra communicator is used for MPI_GET
       CALL MPI_INTERCOMM_MERGE( children(childid)%inter_comm, .FALSE.,                            &
                                 children(childid)%intra_comm, istat )
       CALL MPI_COMM_RANK( children(childid)%intra_comm, children(childid)%intra_rank, istat )

       ALLOCATE( children(childid)%pes(children(childid)%inter_npes) )

!
!--    Allocate array of TYPE arraydef for all child PEs to store information of the transfer array
       DO  j = 1, children(childid)%inter_npes
         ALLOCATE( children(childid)%pes(j)%array_list(pmc_max_array) )
       ENDDO

       CALL get_da_names_from_child( childid )
    ENDDO

 END SUBROUTINE pmc_parentinit


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> thread 0 transfers the index list, which contains all parent grid cells involved in
!> parent client data transfer to the thread, on which this grid cell is located
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_set_2d_index_list( childid, index_list )

     INTEGER(iwp) ::  ian        !<
     INTEGER(iwp) ::  i          !<
     INTEGER(iwp) ::  ip         !<
     INTEGER(iwp) ::  istat      !<
     INTEGER(iwp) ::  max_cells  !<

     INTEGER(iwp), INTENT(IN) ::  childid  !<

     INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  index_list  !<

     INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  cells_on_pe  !<

     INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  lo_ind_list  !<


     IF ( m_model_rank == 0 )  THEN
!
!--     Compute maximum number of grid cells located on one parent thread
        ALLOCATE( cells_on_pe(0:m_model_npes-1) )
        cells_on_pe = 0

        DO  i = 1, SIZE( index_list, 2 )
           cells_on_pe(index_list(6,i )) = cells_on_pe(index_list(6,i ))+1
        ENDDO

        max_cells = MAXVAL( cells_on_pe )

!
!--     Allocate temp array for thread dependent transfer of index_list
        ALLOCATE( lo_ind_list(SIZE( index_list,1 ),max_cells) )

        DO  ip = 0, m_model_npes-1
!
!--        Split into parent processes
           ian = 0

           DO  i = 1, SIZE( index_list, 2 )
              IF ( index_list(6,i) == ip )  THEN
                 ian = ian+1
                 lo_ind_list(:,ian) = index_list(:,i)
              ENDIF
           ENDDO

!
!--        Send data to other parent processes
           IF ( ip == 0 )  THEN
              indchildren(childid)%nrpoints = ian
!
!--           Allocate array for index_list_2d. Note, the array will also be allocated in case
!--           ian = 0, in order to avoid errors when array bounds are checked.
              ALLOCATE( indchildren(childid)%index_list_2d(6,1:ian) )
              IF ( ian > 0 )  THEN
                  indchildren(childid)%index_list_2d(:,1:ian) = lo_ind_list(:,1:ian)
              ENDIF
           ELSE
              CALL MPI_SEND( ian, 1, MPI_INTEGER, ip, 1000, m_model_comm, istat )
              IF ( ian > 0 )  THEN
                  CALL MPI_SEND( lo_ind_list, 6*ian, MPI_INTEGER, ip, 1001, m_model_comm, istat )
              ENDIF
           ENDIF
        ENDDO

        DEALLOCATE( lo_ind_list )
        DEALLOCATE( cells_on_pe )
     ELSE
        CALL MPI_RECV( indchildren(childid)%nrpoints, 1, MPI_INTEGER, 0, 1000, m_model_comm,       &
                       MPI_STATUS_IGNORE, istat )
        ian = indchildren(childid)%nrpoints
!
!--     Allocate array for index_list_2d. Note, the array will also be allocated in case ian=0, in
!--     order to avoid errors when array bounds are checked.
        ALLOCATE( indchildren(childid)%index_list_2d(6,1:ian) )
        IF ( ian > 0 )  THEN
           CALL MPI_RECV( indchildren(childid)%index_list_2d, 6*ian, MPI_INTEGER, 0, 1001,         &
                          m_model_comm, MPI_STATUS_IGNORE, istat)
        ENDIF
     ENDIF

     CALL set_pe_index_list( children(childid), indchildren(childid)%index_list_2d,                &
                             indchildren(childid)%nrpoints )

 END SUBROUTINE pmc_s_set_2d_index_list


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Before creating an array list with arrays schedule for parent client transfer
!> make sure that the list is empty
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_clear_next_array_list

    next_array_in_list = 0

 END SUBROUTINE pmc_s_clear_next_array_list


 LOGICAL FUNCTION pmc_s_getnextarray( childid, myname )

!
!-- Althoug there are no linked lists any more in PMC, this call still looks like working with a list
    CHARACTER(LEN=*), INTENT(OUT) ::  myname  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    TYPE(pedef),    POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar  !<


    next_array_in_list = next_array_in_list + 1

!
!-- Array names are the same on all children processes, so take first process to get the name
    ape => children(childid)%pes(1)

    IF ( next_array_in_list > ape%nr_arrays )  THEN
!
!--    All arrays are done
       pmc_s_getnextarray = .FALSE.
       RETURN
    ENDIF

    ar => ape%array_list(next_array_in_list)
    myname = ar%name

!
!-- Return true if there is still an array in the list
    pmc_s_getnextarray = .TRUE.

 END FUNCTION pmc_s_getnextarray


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> add 2D real array to list of arrays scheduled for parent-client transfer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_set_dataarray_2d( childid, array, array_2 )

    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER ::  array  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER, OPTIONAL ::  array_2  !<

    TYPE(C_PTR) ::  array_adr   !<
    TYPE(C_PTR) ::  second_adr  !<


    dims      = 1
    nrdims    = 2
    dims(1)   = SIZE( array, 1 )
    dims(2)   = SIZE( array, 2 )
    array_adr = C_LOC( array )

    IF ( PRESENT( array_2 ) )  THEN
       second_adr = C_LOC( array_2 )
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr, second_adr = second_adr )
    ELSE
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr )
    ENDIF

 END SUBROUTINE pmc_s_set_dataarray_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> add 2D integer array to list of arrays scheduled for parent-client transfer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_set_dataarray_ip2d( childid, array )

    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(idp), INTENT(IN), DIMENSION(:,:), POINTER ::  array  !<

    TYPE(C_PTR) ::  array_adr  !<


    dims      = 1
    nrdims    = 2
    dims(1)   = SIZE( array, 1 )
    dims(2)   = SIZE( array, 2 )
    array_adr = C_LOC( array )

    CALL pmc_s_setarray( childid, nrdims, dims, array_adr , dimkey = 22 )

 END SUBROUTINE pmc_s_set_dataarray_ip2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> add 3D real array to list of arrays scheduled for parent-client transfer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_set_dataarray_3d( childid, array, nz_cl, nz, array_2, ks_cl, ke_cl )

    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp), INTENT(IN) ::  nz       !<
    INTEGER(iwp), INTENT(IN) ::  nz_cl    !<
    INTEGER(iwp), INTENT(IN) ::  ks_cl    !<
    INTEGER(iwp), INTENT(IN) ::  ke_cl    !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER ::  array  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER, OPTIONAL ::  array_2  !<

    TYPE(C_PTR) ::  array_adr   !<
    TYPE(C_PTR) ::  second_adr  !<


    nrdims  = 3
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )
    dims(3) = SIZE( array, 3 )
    dims(4) = nz_cl+dims(1)-nz  ! Works for first dimension 1:nz and 0:nz+1

    array_adr = C_LOC( array )
!
!-- In PALM's pointer version, two indices have to be stored internally.
!-- The active address of the data array is set in swap_timelevel.
    IF ( PRESENT( array_2 ) )  THEN
       second_adr = C_LOC( array_2 )
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr, second_adr = second_adr, ks=ks_cl, ke=ke_cl)
    ELSE
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr, ks=ks_cl, ke=ke_cl )
    ENDIF

 END SUBROUTINE pmc_s_set_dataarray_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Naming convention for appendices:   _pc  -> parent to child transfer
!>                                     _cp  -> child to parent transfer
!>                                     send -> parent to child transfer
!>                                     recv -> child to parent transfer
!>
!> @todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_setind_and_allocmem( childid )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    INTEGER(iwp) ::  arlen         !<
    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  ierr          !<
    INTEGER(iwp) ::  j             !<
    INTEGER(iwp) ::  lo_nr_arrays  !< store number of arrays in  local variiab le
    INTEGER(iwp) ::  myindex       !<
    INTEGER(iwp) ::  total_npes    !< Total Number of PEs Parent and Child
    INTEGER(idp) ::  bufsize       !< size of MPI data window

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_s  !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_r  !<

    REAL(wp),DIMENSION(:), POINTER, SAVE ::  base_array_cp  !< base array for child to parent transfer
    REAL(wp),DIMENSION(:), POINTER, SAVE ::  base_array_pc  !< base array for parent to child transfer

    TYPE(C_PTR) ::  base_ptr  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<


    call MPI_COMM_SIZE( children(childid)%intra_comm, total_npes, ierr )
!
!-- Parent to child direction
    myindex = 1
    bufsize = 8

!
!-- All Child processes get the same number of arrays.
!-- Therfore the number of arrays form the first Child process can be used for Dimension.
    lo_nr_arrays = children(childid)%pes(1)%nr_arrays

    ALLOCATE( myindex_s(lo_nr_arrays,0:total_npes-1) )
    ALLOCATE( myindex_r(lo_nr_arrays,0:total_npes-1) )

    myindex_s = 0

!
!-- First stride: compute size and set index
    DO  i = 1, children(childid)%inter_npes
       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar  => ape%array_list(j)
          IF ( ar%nrdims == 2 )  THEN
             arlen = ape%nrele
          ELSEIF ( ar%nrdims == 3 )  THEN
             arlen = ape%nrele * ar%a_dim(4)
          ELSE
             arlen = -1
          ENDIF
          ar%sendindex = myindex

!
!--       Using intra communicator for MPU_Alltoall, the numbers of the child processes are after
!--       the parent ones.
          myindex_s(j,i-1+children(childid)%model_npes) = myindex

          myindex = myindex + arlen
          bufsize = bufsize + arlen
          ar%sendsize = arlen
       ENDDO
    ENDDO

!
!-- Using MPI_Alltoall to send indices from  Parent to Child
!-- The data comming back from the child processes are ignored.
    CALL MPI_ALLTOALL( myindex_s, lo_nr_arrays, MPI_INTEGER, myindex_r, lo_nr_arrays, MPI_INTEGER, &
                       children(childid)%intra_comm, ierr )

!
!-- Using MPI_Alltoall to receive indices from Child
    myindex_s = 0
    myindex_r = 0

    CALL MPI_ALLTOALL( myindex_s, lo_nr_arrays, MPI_INTEGER, myindex_r, lo_nr_arrays, MPI_INTEGER, &
                       children(childid)%intra_comm, ierr )

!
!-- Create RMA (One Sided Communication) window for data buffer parent to child transfer.
!-- The buffer of MPI_GET (counterpart of transfer) can be PE-local, i.e. it can but must not be
!-- part of the MPI RMA window. Only one RMA window is required to prepare the data for:
!--          parent -> child transfer on the parent side
!-- and for:
!--          child -> parent transfer on the child side
    CALL pmc_alloc_mem( base_array_pc, bufsize )
    children(childid)%totalbuffersize = bufsize * wp

    winsize = bufsize * wp
    CALL MPI_WIN_CREATE( base_array_pc, winsize, wp, MPI_INFO_NULL, children(childid)%intra_comm,  &
                         children(childid)%win_parent_child, ierr )

!
!-- Open window to set data
    CALL MPI_WIN_FENCE( 0, children(childid)%win_parent_child, ierr )

!
!-- Second stride: set buffer pointer
    DO  i = 1, children(childid)%inter_npes
       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%sendbuf = C_LOC( base_array_pc(ar%sendindex) )
          IF ( ar%sendindex + ar%sendsize > bufsize )  THEN
             WRITE( message_string, '(a,i4,4i7,1x,a)' ) 'parent buffer too small ',i ,             &
                    ar%sendindex, ar%sendsize, ar%sendindex + ar%sendsize, bufsize, TRIM( ar%name )
             CALL message( 'pmc_s_setind_and_allocmem', 'PMC0004', 3, 2, 0, 6, 0 )
          ENDIF
       ENDDO
    ENDDO

!
!-- Child to parent direction
    bufsize = 8

!
!-- First stride: compute size and set index
    DO  i = 1, children(childid)%inter_npes
       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
!
!--       Receive index from child
          IF ( ar%nrdims == 3 )  THEN
             bufsize = MAX( bufsize, INT( ape%nrele * ar%a_dim(4), MPI_ADDRESS_KIND ) )
          ELSE
             bufsize = MAX( bufsize, INT( ape%nrele, MPI_ADDRESS_KIND ) )
          ENDIF
          ar%recvindex = myindex_r(j,i-1+children(childid)%model_npes)
        ENDDO
    ENDDO

    DEALLOCATE( myindex_s )
    DEALLOCATE( myindex_r )

!
!-- Create RMA (one sided communication, RMA = Remote Memory Access) data buffer.
!-- The buffer for MPI_GET can be PE local, i.e. it can but must not be part of the MPI RMA window.
    CALL pmc_alloc_mem( base_array_cp, bufsize, base_ptr )
    children(childid)%totalbuffersize = bufsize * wp

    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

!
!-- Second stride: set buffer pointer
    DO  i = 1, children(childid)%inter_npes
       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvbuf = base_ptr
       ENDDO
    ENDDO

 END SUBROUTINE pmc_s_setind_and_allocmem


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fill buffer in RMA window to enable the client to fetch the dat with MPI_Get
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_fillbuffer( childid, waittime, particle_transfer )

    INTEGER(iwp) ::  ierr     !<
    INTEGER(iwp) ::  ij       !<
    INTEGER(iwp) ::  ip       !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  myindex  !<

    INTEGER(iwp), INTENT(IN) ::  childid   !<

    INTEGER(iwp), DIMENSION(1) ::  buf_shape  !<

    INTEGER(idp), POINTER, DIMENSION(:) ::  ibuf  !<

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  idata_2d  !<

    LOGICAL ::  lo_ptrans  !<

    LOGICAL, INTENT(IN), OPTIONAL   ::  particle_transfer  !<

    REAL(wp) ::  t1  !<
    REAL(wp) ::  t2  !<

    REAL(wp), INTENT(OUT), OPTIONAL ::  waittime  !<

    REAL(wp), POINTER, DIMENSION(:) ::  buf  !<

    REAL(wp), POINTER, DIMENSION(:,:) ::  data_2d  !<

    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<

!
!-- Synchronization of the model is done in pmci_synchronize. Therefor the RMA window can be filled
!-- without sychronization at this point and a barrier is not necessary.
!-- Please note that waittime has to be set in pmc_s_fillbuffer AND pmc_c_getbuffer.
    IF ( PRESENT( waittime ) )  THEN
      t1 = pmc_time()
      CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
      t2 = pmc_time()
      waittime = t2 - t1
    ENDIF

    lo_ptrans = .FALSE.
    IF ( PRESENT( particle_transfer ) )    lo_ptrans = particle_transfer

    DO  ip = 1, children(childid)%inter_npes
       ape => children(childid)%pes(ip)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          myindex = 1
          IF ( ar%dimkey == 2 .AND. .NOT.lo_ptrans  )  THEN         ! PALM 2D REAL*8 Array
             buf_shape(1) = ape%nrele
             CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
             CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                buf(myindex) = data_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + 1
             ENDDO
          ELSEIF ( ar%dimkey == 3  .AND. .NOT.lo_ptrans  )  THEN      ! PALM 3D REAL*8 Array
             buf_shape(1) = ape%nrele*ar%a_dim(4)
             CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
             CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
             DO  ij = 1, ape%nrele
                buf(myindex:myindex+ar%a_dim(4)-1) =                                               &
                data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + ar%a_dim(4)
             ENDDO
          ELSEIF ( ar%dimkey == 22 .AND. lo_ptrans  )  THEN  ! 2D INTEGER*8 Array for particle Transfer
             buf_shape(1) = ape%nrele
             CALL C_F_POINTER( ar%sendbuf, ibuf, buf_shape )
             CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                ibuf(myindex) = idata_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + 1
             ENDDO
          ENDIF
        ENDDO
    ENDDO

!
!-- Buffer is filled
    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

 END SUBROUTINE pmc_s_fillbuffer


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Get client data from RMM window
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_getdata_from_buffer( childid, waittime , particle_transfer, child_process_nr )

    INTEGER(iwp) ::  ierr       !<
    INTEGER(iwp) ::  ij         !<
    INTEGER(iwp) ::  ip         !<
    INTEGER(iwp) ::  ip_start   !<
    INTEGER(iwp) ::  ip_end     !<
    INTEGER(iwp) ::  j          !<
    INTEGER(iwp) ::  myindex    !<
    INTEGER(iwp) ::  nr         !<
    INTEGER(iwp) ::  target_pe  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  child_process_nr  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_disp  !<

    INTEGER(iwp), DIMENSION(1) ::  buf_shape  !<

    INTEGER(idp), POINTER, DIMENSION(:) ::  ibuf  !<

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  idata_2d  !<

    LOGICAL ::  lo_ptrans  !<

    LOGICAL, INTENT(IN), OPTIONAL ::  particle_transfer  !<

    REAL(wp) ::  t1  !<
    REAL(wp) ::  t2  !<

    REAL(wp), INTENT(OUT), OPTIONAL ::  waittime  !<

    REAL(wp), POINTER, DIMENSION(:) ::  buf  !<

    REAL(wp), POINTER, DIMENSION(:,:) ::  data_2d  !<

    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(pedef), POINTER  ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<


    t1 = pmc_time()

    IF( PRESENT( child_process_nr ) )  THEN
       ip_start = child_process_nr
       ip_end   = child_process_nr
    ELSE
       ip_start = 1
       ip_end   = children(childid)%inter_npes
    END IF

    lo_ptrans = .FALSE.
    IF ( PRESENT( particle_transfer ) )  lo_ptrans = particle_transfer

    IF(ip_start == 1)  THEN
!
!--    Wait for child to fill buffer
       CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
       t2 = pmc_time() - t1
       IF ( PRESENT( waittime ) )  waittime = t2

       CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
    ENDIF

    DO  ip = ip_start, ip_end
       ape => children(childid)%pes(ip)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)

          IF ( ar%recvindex < 0 )  CYCLE

          IF ( ar%dimkey == 2  .AND.  .NOT. lo_ptrans  )  THEN
             nr = ape%nrele
          ELSEIF ( ar%dimkey == 3  .AND.  .NOT. lo_ptrans  )  THEN
             nr = ape%nrele * ar%a_dim(4)
          ELSE IF ( ar%dimkey == 22 .AND. lo_ptrans)  THEN
             nr = ape%nrele
          ELSE
             CYCLE            ! Particle arrays are not transfered here
          ENDIF

          buf_shape(1) = nr

          IF(lo_ptrans)   THEN
             CALL C_F_POINTER( ar%recvbuf, ibuf, buf_shape )
          ELSE
             CALL C_F_POINTER( ar%recvbuf, buf, buf_shape )
          ENDIF

!
!--       MPI passive target RMA
          IF ( nr > 0 )  THEN
             target_disp = ar%recvindex - 1
!
!--          Child processes are located behind parent process
             target_pe = ip - 1 + m_model_npes
             CALL MPI_WIN_LOCK( MPI_LOCK_SHARED, target_pe, 0, children(childid)%win_parent_child, &
                                ierr )
             IF ( lo_ptrans )  THEN
                CALL MPI_GET( ibuf, nr * 8, MPI_BYTE, target_pe, target_disp, nr * 8, MPI_BYTE,    &
                              children(childid)%win_parent_child, ierr )
             ELSE
                CALL MPI_GET( buf, nr, MPI_REAL, target_pe, target_disp, nr, MPI_REAL,             &
                              children(childid)%win_parent_child, ierr )
             ENDIF
             CALL MPI_WIN_UNLOCK( target_pe, children(childid)%win_parent_child, ierr )
          ENDIF

          myindex = 1

          IF ( ar%dimkey == 2  .AND.  .NOT. lo_ptrans )  THEN
             CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                data_2d(ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex)
                myindex = myindex + 1
             ENDDO
          ELSE IF ( ar%dimkey == 3  .AND.  .NOT. lo_ptrans )  THEN
             CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
             DO  ij = 1, ape%nrele
                data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i) =                         &
                buf(myindex:myindex+ar%a_dim(4)-1)
                myindex = myindex + ar%a_dim(4)
             ENDDO
          ELSE IF ( ar%dimkey == 22 .AND. lo_ptrans)  THEN
             CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                idata_2d(ape%locind(ij)%j,ape%locind(ij)%i) = ibuf(myindex)
                myindex = myindex + 1
             ENDDO
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE pmc_s_getdata_from_buffer


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> broadcast name of transfer arrays from child thread 0 to parent threads
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_da_names_from_child( childid )

!
!-- Get data array description and name from child
    INTEGER(iwp), INTENT(IN) ::  childid  !<

    TYPE(da_namedef) ::  myname  !<


    DO
       CALL pmc_bcast( myname%couple_index, 0, comm=m_to_child_comm(childid) )

       IF ( myname%couple_index == -1 )  EXIT

       CALL pmc_bcast( myname%parentdesc,   0, comm=m_to_child_comm(childid) )
       CALL pmc_bcast( myname%nameonparent, 0, comm=m_to_child_comm(childid) )
       CALL pmc_bcast( myname%childdesc,    0, comm=m_to_child_comm(childid) )
       CALL pmc_bcast( myname%nameonchild,  0, comm=m_to_child_comm(childid) )

       CALL pmc_g_setname( children(childid), myname%couple_index, myname%nameonparent )
   ENDDO

 END SUBROUTINE get_da_names_from_child


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_setarray( childid, nrdims, dims, array_adr, second_adr, dimkey, ks, ke )

!
!-- Set array for child inter process 0
    INTEGER(iwp) ::  i  !< local counter

    INTEGER(iwp), INTENT(IN) :: childid  !<
    INTEGER(iwp), INTENT(IN) :: nrdims   !<
    INTEGER(iwp), INTENT(IN), OPTIONAL :: ks       !<
    INTEGER(iwp), INTENT(IN), OPTIONAL :: ke       !<

    INTEGER(iwp), INTENT(IN), OPTIONAL :: dimkey  !<

    INTEGER(iwp), INTENT(IN), DIMENSION(:) :: dims  !<

    TYPE(C_PTR), INTENT(IN) :: array_adr  !<

    TYPE(C_PTR), INTENT(IN), OPTIONAL ::  second_adr  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<


    DO  i = 1, children(childid)%inter_npes
       ape => children(childid)%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%dimkey = nrdims

       IF( PRESENT( dimkey ) )  ar%dimkey = dimkey
       ar%a_dim  = dims
       ar%data   = array_adr
       IF ( PRESENT( second_adr ) )  THEN
          ar%po_data(1) = array_adr
          ar%po_data(2) = second_adr
       ELSE
          ar%po_data(1) = C_NULL_PTR
          ar%po_data(2) = C_NULL_PTR
       ENDIF

       IF( PRESENT (ks) ) THEN
          ar%ks = ks
       END IF
       IF( PRESENT (ke) ) THEN
          ar%ke = ke
       ENDIF
    ENDDO

 END SUBROUTINE pmc_s_setarray


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_set_active_data_array( childid, iactive )

    INTEGER(iwp) :: ip  !<
    INTEGER(iwp) :: j   !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp), INTENT(IN) ::  iactive  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<

    DO  ip = 1, children(childid)%inter_npes
       ape => children(childid)%pes(ip)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          IF ( MOD( ar%dimkey, 10 ) == 2 )  CYCLE  !Not for 2D array
          IF ( iactive == 1  .OR.  iactive == 2 )  THEN
             ar%data = ar%po_data(iactive)
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE pmc_s_set_active_data_array


 !--------------------------------------------------------------------------------------------------!
 ! Description:
 ! ------------
 !> @todo: Missing function description.
 !--------------------------------------------------------------------------------------------------!
 INTEGER FUNCTION pmc_s_get_child_npes( child_id )

   INTEGER(iwp),INTENT(IN) ::  child_id  !<

   pmc_s_get_child_npes = children(child_id)%inter_npes

   RETURN

 END FUNCTION pmc_s_get_child_npes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE set_pe_index_list( mychild, index_list, nrp )

    INTEGER(iwp) :: i        !<
    INTEGER(iwp) :: ierr     !<
    INTEGER(iwp) :: ind      !<
    INTEGER(iwp) :: indwin   !<
    INTEGER(iwp) :: indwin2  !<
    INTEGER(iwp) :: i2       !<
    INTEGER(iwp) :: j        !<
    INTEGER(iwp) :: rempe    !<

    TYPE(childdef), INTENT(INOUT) ::  mychild  !<

    INTEGER(iwp), DIMENSION(mychild%inter_npes) ::  remind  !<

    INTEGER(iwp), INTENT(IN) ::  nrp  !<

    INTEGER(iwp), INTENT(IN), DIMENSION(:,:) ::  index_list  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<

    INTEGER(iwp), DIMENSION(:), POINTER ::  remindw  !<
    INTEGER(iwp), DIMENSION(:), POINTER ::  rldef    !<

    TYPE(pedef), POINTER ::  ape  !<


!
!-- First, count entries for every remote child process
    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ape%nrele = 0
    ENDDO

!
!-- Loop over number of parent grid cells.
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1   ! Process number on remote process
       ape => mychild%pes(rempe)
       ape%nrele = ape%nrele + 1     ! Increment number of elements for this child process
    ENDDO

    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ALLOCATE( ape%locind(ape%nrele) )
    ENDDO

    remind = 0

!
!-- Second, create lists
!-- Loop over number of parent grid cells.
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1
       ape => mychild%pes(rempe)
       remind(rempe)     = remind(rempe) + 1
       ind               = remind(rempe)
       ape%locind(ind)%i = index_list(1,j)
       ape%locind(ind)%j = index_list(2,j)
    ENDDO

!
!-- Prepare number of elements for children processes
    CALL pmc_alloc_mem( rldef, mychild%inter_npes * 2 )

!
!-- Number of child processes * size of INTEGER (i just arbitrary INTEGER)
    winsize = mychild%inter_npes * STORAGE_SIZE( i ) / 8 * 2

    CALL MPI_WIN_CREATE( rldef, winsize, iwp, MPI_INFO_NULL, mychild%intra_comm, indwin, ierr )

!
!-- Open window to set data
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

    rldef(1) = 0            ! Index on remote process 0
    rldef(2) = remind(1)    ! Number of elements on remote process 0

!
!-- Reserve buffer for index array
    DO  i = 2, mychild%inter_npes
       i2          = ( i - 1 ) * 2 + 1
       rldef(i2)   = rldef(i2-2) + rldef(i2-1) * 2  ! Index on remote process
       rldef(i2+1) = remind(i)                      ! Number of elements on remote process
    ENDDO

!
!-- Close window to allow child to access data
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

!
!-- Child has retrieved data
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

    i2 = 2 * mychild%inter_npes - 1
    winsize = ( rldef(i2) + rldef(i2+1) ) * 2

!
!-- Make sure, MPI_ALLOC_MEM works
    winsize = MAX( winsize, INT( 1, MPI_ADDRESS_KIND ) )

    CALL pmc_alloc_mem( remindw, INT( winsize ) )

    CALL MPI_BARRIER( m_model_comm, ierr )
    CALL MPI_WIN_CREATE( remindw, winsize * STORAGE_SIZE( i ) / 8, iwp, MPI_INFO_NULL,             &
                         mychild%intra_comm, indwin2, ierr )

!
!-- Open window to set data
    CALL MPI_WIN_FENCE( 0, indwin2, ierr )

!
!-- Create the 2D index list
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1    ! Process number on remote process
       ape => mychild%pes(rempe)
       i2    = rempe * 2 - 1
       ind   = rldef(i2) + 1
       remindw(ind)   = index_list(3,j)
       remindw(ind+1) = index_list(4,j)
       rldef(i2)      = rldef(i2)+2
    ENDDO

!
!-- All data are set
    CALL MPI_WIN_FENCE( 0, indwin2, ierr )

!
!-- Don't know why, but this barrier is necessary before windows can be freed
!-- TODO: find out why this is required
    CALL MPI_BARRIER( mychild%intra_comm, ierr )

    CALL MPI_WIN_FREE( indwin, ierr )
    CALL MPI_WIN_FREE( indwin2, ierr )

!
!-- TODO: check if the following idea needs to be done
!-- Should work, Problem with MPI implementation
!-- https://www.lrz.de/services/software/parallel/mpi/onesided
!-- CALL MPI_Free_mem (remindw, ierr)

 END SUBROUTINE set_pe_index_list


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> free window of pmc data buffer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_s_finalize( childid )

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp) :: ierr  !<


    CALL MPI_WIN_FREE( children(childid)%win_parent_child, ierr )

 END SUBROUTINE pmc_s_finalize
#endif

 END MODULE pmc_parent
