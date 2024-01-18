!> @file poisfft_sm_mod.f90
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
! Copyright 1997-2021 Leibniz Universitaet Hannover, Klaus Ketelsen
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Solves the Poisson equation with a 2D spectral method
!>        d^2 p / dx^2 + d^2 p / dy^2 + d^2 p / dz^2 = s
!>
!> Input:
!> real   ar   contains (nnz,nny,nnx) elements of the velocity divergence, starting from (1,nys,nxl)
!>
!> Output:
!> real   ar   contains the solution for perturbation pressure p
!>
!> This version of the Poisson-FFT-solver is based on a 1d-decomposition using shared memory,
!> while the remaining part of the model (may) use a 2d-decomposition.
!--------------------------------------------------------------------------------------------------!
 MODULE poisfft_sm_mod

#if defined( __parallel )
    USE MPI

    USE control_parameters,                                                                        &
        ONLY:  debug_output

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE fft_xy,                                                                                    &
        ONLY:  fft_init,                                                                           &
               fft_y,                                                                              &
               fft_y_1d,                                                                           &
               fft_x,                                                                              &
               fft_x_1d

    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nnx_pe,                                                                             &
               nny,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxl_pe,                                                                             &
               nxr,                                                                                &
               nxr_pe,                                                                             &
               ny,                                                                                 &
               nys,                                                                                &
               nyn,                                                                                &
               nz

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  collective_wait,                                                                    &
               comm1dx,                                                                            &
               comm2d,                                                                             &
               myid,                                                                               &
               myidx,                                                                              &
               non_uniform_data_for_transpose,                                                     &
               npex,                                                                               &
               npey,                                                                               &
               sendrecvcount_xy,                                                                   &
               threads_per_task

    USE sm_poisfft_mod,                                                                            &
        ONLY:  sm_poisfft,                                                                         &
               sm_setup

    USE transpose_mod,                                                                             &
        ONLY:  nxl_y,                                                                              &
               nxl_z,                                                                              &
               nxr_y,                                                                              &
               nxr_z,                                                                              &
               nyn_x,                                                                              &
               nyn_z,                                                                              &
               nyn_z_max,                                                                          &
               nys_x,                                                                              &
               nys_z,                                                                              &
               nzb_x,                                                                              &
               nzb_y,                                                                              &
               nzt_x,                                                                              &
               nzt_y

    USE tridia_solver,                                                                             &
        ONLY:  tridia_init,                                                                        &
               tridia_substi

    IMPLICIT NONE

    PRIVATE
    SAVE

    TYPE(sm_setup) ::  pe_grid    !<
    TYPE(sm_setup) ::  node_grid  !<

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  nys_x_pe  !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  nyn_x_pe  !<

    LOGICAL ::  poisfft_sm_initialized = .FALSE.  !<

!
!-- Shared arrays.
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  work_trix  !<
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  yz_buf     !<

    REAL(wp), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::  xz_buf  !<

!
!-- Window variables for shared arrays.
    INTEGER(iwp) ::  nnx_max      !<
    INTEGER(iwp) ::  nny_x        !<
    INTEGER(iwp) ::  nny_x_max    !<
    INTEGER(iwp) ::  nyn_x_work   !<
    INTEGER(iwp) ::  nys_x_work   !<
    INTEGER(iwp) ::  ny_dim       !<
    INTEGER(iwp) ::  nz_max       !<
    INTEGER(iwp) ::  trix_win     !<
    INTEGER(iwp) ::  xend         !<
    INTEGER(iwp) ::  xstart       !<
    INTEGER(iwp) ::  xzbuftype    !<
    INTEGER(iwp) ::  xz_win       !<
    INTEGER(iwp) ::  yz_win       !<
    INTEGER(iwp) ::  zend         !<
    INTEGER(iwp) ::  zend_max     !<
    INTEGER(iwp) ::  zend_work    !<
    INTEGER(iwp) ::  zstart       !<
    INTEGER(iwp) ::  zstart_work  !<

    INTERFACE ffty_tr_yx
       MODULE PROCEDURE ffty_tr_yx
    END INTERFACE ffty_tr_yx

    INTERFACE poisfft_sm
       MODULE PROCEDURE poisfft_sm
    END INTERFACE poisfft_sm

    INTERFACE poisfft_sm_init
       MODULE PROCEDURE poisfft_sm_init
    END INTERFACE poisfft_sm_init

    INTERFACE tr_xy_ffty
       MODULE PROCEDURE tr_xy_ffty
    END INTERFACE tr_xy_ffty

    PUBLIC poisfft_sm_init, poisfft_sm

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Setup coefficients for FFT and the tridiagonal solver.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poisfft_sm_init

    IMPLICIT NONE

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  ierr
    INTEGER(iwp) ::  irest_x  !<
    INTEGER(iwp) ::  irest_y  !<
    INTEGER(iwp) ::  irest_z  !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  k        !<
    INTEGER(iwp) ::  xinc     !<
    INTEGER(iwp) ::  zinc     !<

    INTEGER(iwp) ::  shp(4)  !<

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  xend_pe    !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  xstart_pe  !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  zend_pe    !<
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  zstart_pe  !<


    CALL fft_init

    poisfft_sm_initialized = .TRUE.

!
!-- Set the shared memory communication features / variables.
    CALL pe_grid%save_grid_into_this_class_poisfft()
    CALL node_grid%save_grid_into_this_class_poisfft()

    node_grid%npex = sm_poisfft%n_npes
    node_grid%npey = 1
    node_grid%myid = sm_poisfft%n_rank
    node_grid%nny  = ny + 1

    ierr = 0

    IF ( non_uniform_data_for_transpose )  THEN

       ALLOCATE( nyn_x_pe(0:npex-1) )
       ALLOCATE( nys_x_pe(0:npex-1) )

       nny_x   = node_grid%nny / node_grid%npex
       irest_y = MOD( node_grid%nny, node_grid%npex )

        DO  j = 0, npex-1
           IF ( j < irest_y )  THEN
              nys_x_pe(j) = j * ( nny_x + 1 )
              nyn_x_pe(j) = nys_x_pe(j) + ( nny_x + 1 ) - 1
           ELSE
              nys_x_pe(j) = irest_y * ( nny_x + 1 ) + ( j - irest_y ) * nny_x
              nyn_x_pe(j) = nys_x_pe(j) + nny_x - 1
           ENDIF
        ENDDO

        node_grid%nys_x = nys_x_pe(myidx)
        node_grid%nyn_x = nyn_x_pe(myidx)
        nny_x = node_grid%nyn_x - node_grid%nys_x + 1
        nny_x_max = nyn_x_pe(0) - nys_x_pe(0) + 1
        ny_dim = ny + npex

        nnx_max = nxr_pe(0) - nxl_pe(0) + 1

        node_grid%sendrecvcount_xy = nnx_max * nny_x_max

     ELSE

        nny_x = node_grid%nny / node_grid%npex
        node_grid%nys_x = node_grid%myid * nny_x
        node_grid%nyn_x = ( node_grid%myid + 1 ) * nny_x - 1
        node_grid%sendrecvcount_xy = nnx * nny_x
        nny_x_max = nny_x
        nnx_max   = nnx
        ny_dim    = ny

     ENDIF

     nys_x_work = myidx * nny_x_max
     nyn_x_work = nys_x_work + nny_x     - 1
     nyn_z_max  = nys_x_work + nny_x_max - 1

!
!--  Switch indices to the set that is used for the 1d-decomposition (shared memory).
     CALL node_grid%activate_grid_from_this_class_poisfft()

     IF ( non_uniform_data_for_transpose )  THEN
!
!--     Compute x-distribution inside shared memory blocks.
        ALLOCATE( xstart_pe(0:sm_poisfft%sh_npes-1) )
        ALLOCATE( xend_pe(0:sm_poisfft%sh_npes-1) )
        xinc    = ( nx + 1 ) / sm_poisfft%sh_npes
        irest_x = MOD( ( nx + 1 ), sm_poisfft%sh_npes )

!
!--     Calculate left and right index bounds along x for each core of the shared memory block.
!--     Required for domain decomposition along x.
        DO  i = 0, sm_poisfft%sh_npes-1
           IF ( i < irest_x )  THEN
!
!--           Remaining grid points (irest_x) are distributed among the lower rank cores.
!--           Therefore,  xinc+1 grid points are assigned.
              xstart_pe(i) = i * ( xinc + 1 )
              xend_pe(i)   = xstart_pe(i) + ( xinc + 1 ) - 1
           ELSE
              xstart_pe(i) = irest_x * ( xinc + 1 ) + ( i - irest_x ) * xinc
              xend_pe(i)   = xstart_pe(i) + xinc - 1
           ENDIF
        ENDDO
        xstart = xstart_pe(sm_poisfft%sh_rank)
        xend   = xend_pe(sm_poisfft%sh_rank)

!
!--     Compute z-distribution inside shared memory blocks.
        ALLOCATE( zstart_pe(0:sm_poisfft%sh_npes-1) )
        ALLOCATE( zend_pe(0:sm_poisfft%sh_npes-1) )
        zinc    = nz / sm_poisfft%sh_npes
        irest_z = MOD( nz, sm_poisfft%sh_npes )

!
!--     Calculate lower and upper index bounds along z for each core of the shared memory block.
!--     Required for domain decomposition along z.
        DO  k = 0, sm_poisfft%sh_npes-1
           IF ( k < irest_z )  THEN
!
!--           Remaining grid points (irest_z) are distributed among the lower rank cores.
!--           Therefore,  zinc+1 grid points are assigned.
              zstart_pe(k) = k * ( zinc + 1 ) + 1
              zend_pe(k)   = zstart_pe(k) + ( zinc + 1 ) - 1
           ELSE
              zstart_pe(k) = irest_z * ( zinc + 1 ) + ( k - irest_z ) * zinc + 1
              zend_pe(k)   = zstart_pe(k) + zinc - 1
           ENDIF
        ENDDO
!
!--     Set the lower and upper index bound along z for this specific PE.
        zstart = zstart_pe(sm_poisfft%sh_rank)
        zend   = zend_pe(sm_poisfft%sh_rank)

!
!--     The z areas in equal sized arrays must not overlap.
!--     Following indices are e.g. used to allocate arrays. Although the actual number of
!--     grid points may vary by one from core to core, the number of array elements must be the same
!--     on each core (required by MPI_ALLTOALL), i.e. the arrays need to be dimensioned with respect
!--     to the number of grid points on core 0 (which has one point more than the larger rank PEs).
!--     zend_max is either zend or zend+1.
        zend_max    = zstart + zend_pe(0) - 1
        zstart_work = sm_poisfft%sh_rank * ( zend_pe(0) - zstart_pe(0) + 1 ) + 1
        zend_work   = zstart_work + zend_pe(0) - zstart_pe(0)
!
!--     Cores with rank > irest_z have one more (dummy) grid point.
        nz_max = nz + sm_poisfft%sh_npes - irest_z

     ELSE

        xinc   = ( nx + 1 ) / sm_poisfft%sh_npes
        xstart = sm_poisfft%sh_rank * xinc
        xend   = xstart + xinc - 1

        zinc   = nz / sm_poisfft%sh_npes
        zstart = sm_poisfft%sh_rank * zinc + 1
        zend   = zstart + zinc - 1
!
!--     Since domains are uniform, all arrays / indices are the same.
        zend_max    = zend
        zstart_work = zstart
        zend_work   = zend
        nz_max      = nz

     ENDIF

     nys_z = nys_x
     nyn_z = nyn_x
     nxl_z = xstart
     nxr_z = xend

     node_grid%nys_z = nys_z
     node_grid%nyn_z = nyn_z
     node_grid%nxl_z = nxl_z
     node_grid%nxr_z = nxr_z

!
!--  Allocate the shared memory arrays.
!--  yz_buf holds data before transposition y -> x (contains all grid points along y, 0:ny).
!--  xz_buf holds data after transposition y -> x (has 4 dimensions, where the 4th dimension
!--  displays the data along x on the respective PE). Will be collected to a 0:nx array in
!--  fftx_tri_fftx.
     CALL sm_poisfft%sm_allocate_shared( yz_buf, 1, nz, 0, ny, nxl, nxr, yz_win)
     CALL sm_poisfft%sm_allocate_shared( xz_buf, 1, nnx_max, 1, nz_max, nys_x_work,                &
                                         nys_x_work+nny_x_max-1, 1, node_grid%npex, xz_win )

     shp = SHAPE( xz_buf )
     shp(4) = 1
     CALL MPI_TYPE_CREATE_SUBARRAY( 4, shp, (/nnx_max,1,nny_x_max,1/), (/0,0,0,0/),                &
                                    MPI_ORDER_FORTRAN, MPI_REAL, xzbuftype, ierr )
     CALL MPI_TYPE_COMMIT( xzbuftype, ierr )

     CALL sm_poisfft%sm_allocate_shared( work_trix, 0, nx, nys_z, nyn_z, 1, nz, trix_win )

!
!--  Initialize/calculate coefficients required for solving the tridiagonal set of equations.
     CALL tridia_init

!
!--  Following output can be removed after shared memory Poisson-FFT-solver has been sufficiently
!--  used in operational mode.
     IF ( debug_output )  THEN
!
!--     Node-grid printout. npex and npey here are node grid values!
        WRITE( 9, '(A,10I5)' )  'init_poisfft_no1 ', npex, npey, nys_x, nyn_x, xstart, xend, xinc, &
                                zstart, zend, zinc
        WRITE( 9, '(A,10I6)' )  'init_poisfft_no2 ', nny_x, node_grid%nys_x, node_grid%nyn_x,      &
                                node_grid%sendrecvcount_xy
        WRITE( 9, '(A,10I6)' )  'init_poisfft_no3 ', nys_x, nyn_x, sendrecvcount_xy,               &
                                sm_poisfft%sh_npes, nys_z, nyn_z, nxl_z, nxr_z, nx
        FLUSH( 9 )
     ENDIF

!
!--  Finally, switch indices back to the set that is used for the 2d-decomposition.
     CALL pe_grid%activate_grid_from_this_class_poisfft()

!
!--  Following output can be removed after shared memory Poisson-FFT-solver has been sufficiently
!--  used in operational mode.
     IF ( debug_output )  THEN
!
!--     PE-grid printout.
        WRITE( 9, '(A,10I6)' )  'init_poisfft_pe  ', nys_x, nyn_x, sendrecvcount_xy, npex, npey,   &
                                nys_z, nyn_z, nxl_z, nxr_z, nx
        FLUSH( 9 )
     ENDIF

 END SUBROUTINE poisfft_sm_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Two-dimensional Fourier Transformation in x- and y-direction.
!> This is the version using a 1d-decomposition based on shared memory. comm1dx is used as
!> communicator.
!> Node barriers are always used after a shared array has been modified.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poisfft_sm( ar )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT), DIMENSION(1:nz,nys:nyn,nxl:nxr) ::  ar  !<


    CALL cpu_log( log_point_s(102), 'poisfft_sm', 'start' )

    IF ( .NOT. poisfft_sm_initialized )  CALL poisfft_sm_init

!
!-- Activate settings for the 1d-decomposition along x. Indices will be redefined, e.g. nxl = 0,
!-- nxr = nx.
    CALL node_grid%activate_grid_from_this_class_poisfft()

!
!-- Fill the shared array. 2nd dimension of yz_buf is 0:ny.
    yz_buf(1:nz,nys:nyn,nxl:nxr) = ar
    CALL sm_poisfft%sm_node_barrier

!
!-- FFT along y and transposition y --> x.
    CALL ffty_tr_yx
    CALL sm_poisfft%sm_node_barrier

!
!-- FFT along x, solving the tridiagonal system and backward FFT.
    CALL fftx_tri_fftx( xz_buf )
    CALL sm_poisfft%sm_node_barrier

!
!-- Transposition x --> y and backward FFT along y.
    CALL tr_xy_ffty
    CALL sm_poisfft%sm_node_barrier

!
!-- Restore the solution from the shared array.
    ar = yz_buf(1:nz,nys:nyn,nxl:nxr)
    CALL sm_poisfft%sm_node_barrier

!
!-- Re-activate the settings for the 2d-domain decomposition.
    CALL pe_grid%activate_grid_from_this_class_poisfft()

    CALL cpu_log( log_point_s(102), 'poisfft_sm', 'stop' )

 END SUBROUTINE poisfft_sm


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y with subsequent transposition y --> x for a 1d-decomposition
!> along x (based on shared memory).
!>
!> The transposition is done for each k-layer separately, to allow for some overlap of communication
!> and calculation.
!>
!> After the transposition, data is not resorted to a 3d-array, i.e. the transposed array xz_buf
!> has four dimensions with the last dimension (1:npex).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ffty_tr_yx

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  ierr  !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  jj    !<
    INTEGER(iwp) ::  je    !<
    INTEGER(iwp) ::  js    !<
    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  kw    !<
    INTEGER(iwp) ::  l     !<

    INTEGER(iwp), DIMENSION(zstart:zend_max) ::  req  !<

    REAL(wp), DIMENSION(nxl:nxl+nnx_max-1,0:ny_dim,zstart_work:zend_work) ::  work       !<
    REAL(wp), DIMENSION(0:ny,zstart:zend,nxl:nxr)                         ::  work_ffty  !<


!
!-- Carry out the FFT along y, where all data are present due to the 1d-decomposition along x.
!-- Resort the data in a way that y becomes the first index. This does not require a transposition
!-- because all data along y are present in the shared memory block (yz_buf).
    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'start' )
    DO  i = nxl, nxr
       DO  k = zstart, zend
          DO  j = 0, ny
             work_ffty(j,k,i) = yz_buf(k,j,i)
          ENDDO
       ENDDO
    ENDDO
    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'pause' )

    kw = zstart_work
    DO  k = zstart, zend_max
!
!--    FFTs are only required for k levels that contain data. Levels beyond zend have been added
!--    only to deal with non-uniform subdomains.
       IF ( k <= zend )  THEN

          CALL cpu_log( log_point_s(7), 'fft_y_1d', 'continue' )
          DO  i = nxl, nxr
!
!--          FFT along y.
             CALL fft_y_1d( work_ffty(:,k,i), 'forward' )
!
!--          Resort data for x as first index.
             IF ( non_uniform_data_for_transpose )  THEN
                jj = 0
                DO  l = 0, npex - 1
                   js = l * nny_x_max
                   je = js + nyn_x_pe(l) - nys_x_pe(l)
                   DO  j = js, je
                      work(i,j,kw) = work_ffty(jj,k,i)
                      jj = jj + 1
                   ENDDO
                ENDDO
             ELSE
                DO  j = 0, ny
                   work(i,j,kw) = work_ffty(j,k,i)
                ENDDO
             ENDIF
          ENDDO
          CALL cpu_log( log_point_s(7), 'fft_y_1d', 'pause' )

       ENDIF
!
!--    A wait for the MPI_IALLTOALL is carried out from the second vertical level on to allow for an
!--    overlap of FFT calculation and data exchange. While the FFT for level k is calculated, the
!--    data exchange via MPI_IALLTOALL for level k-1 is still on its way.
       IF ( k >= zstart + 1 )  THEN
          CALL cpu_log( log_point_s(101), 'mpi_ialltoall_2', 'start' )
          CALL MPI_WAIT( req(k-1), MPI_STATUS_IGNORE, ierr )
          CALL cpu_log( log_point_s(101), 'mpi_ialltoall_2', 'stop' )
       ENDIF
!
!--    Transpose data y --> x for level k.
!--    MPI_IALLTOALL needs the same number of calls on every PE. Therefore the z dimension in xz_buf
!--    must be dimensioned with nnz taken from PE0.
       CALL cpu_log( log_point_s(100), 'mpi_ialltoall_1', 'start' )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_IALLTOALL( work(nxl,0,kw), sendrecvcount_xy, MPI_REAL,                             &
                           xz_buf(1,kw,nys_x_work,1), 1, xzbuftype, comm1dx, req(k), ierr )
       CALL cpu_log( log_point_s(100), 'mpi_ialltoall_1', 'stop' )

       kw = kw + 1

    ENDDO

!
!-- Wait for the final remaining k level zend_max.
    CALL cpu_log( log_point_s(101), 'mpi_ialltoall_2', 'start' )
    CALL MPI_WAIT( req(zend_max), MPI_STATUS_IGNORE, ierr )
    CALL cpu_log( log_point_s(101), 'mpi_ialltoall_2', 'stop' )

 END SUBROUTINE ffty_tr_yx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> FFT along x, solution of the tridiagonal system and backward FFT for a 1d-decomposition along x,
!> based on shared memory.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fftx_tri_fftx( ar )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT),                                                                       &
              DIMENSION(nnx_max,1:nz_max,nys_x_work:nys_x_work+nny_x_max-1,npex) ::  ar  !<

    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  j   !<
    INTEGER(iwp) ::  jj  !<
    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  kw  !<
    INTEGER(iwp) ::  m   !<
    INTEGER(iwp) ::  n   !<

    REAL(wp), DIMENSION(0:nx) ::  work_fftx  !<


    CALL cpu_log( log_point_s(33), 'fft_x_1d ', 'start' )

    jj = nys_x_work

!--   OpenMP Version not implemented yet, jj needs action
!--   !$OMP  PARALLEL DO PRIVATE ( i, j, k, m, n, work_fftx )
!
!-- Carry out the fast Fourier transformation along x.
    DO  j = nys_z, nyn_z
       kw = zstart_work
       DO  k = zstart, zend
!
!--       Collect the data for the FFT along x in a 1d-array. The source array has four dimensions
!--       because it hasn't been resorted after the transposition y --> x.
          m = 0
          IF ( non_uniform_data_for_transpose )  THEN
             DO  n = 1, npex
                DO  i = 1, nnx_pe(n-1)
                   work_fftx(m) = ar(i,kw,jj,n)
                   m = m + 1
                ENDDO
             ENDDO
          ELSE
             DO  n = 1, npex
                DO  i = 1, nnx
                   work_fftx(m) = ar(i,k,j,n)
                   m = m + 1
                ENDDO
             ENDDO
          ENDIF
!
!--       FFT along x.
          CALL fft_x_1d( work_fftx, 'forward' )
!
!--       Store result in a work array to be used by the tri-diagonal solver below.
          DO  i = 0, nx
             work_trix(i,j,k) = work_fftx(i)
          ENDDO
          kw = kw + 1
       ENDDO
       jj = jj + 1
    ENDDO
    CALL cpu_log( log_point_s(33), 'fft_x_1d ', 'pause' )

    CALL sm_poisfft%sm_node_barrier
!
!-- Solve the linear equation system.
    CALL cpu_log( log_point_s(6), 'tridia', 'start' )
    CALL tridia_substi (work_trix, 0, nx)
    CALL cpu_log( log_point_s(6), 'tridia', 'stop' )

    CALL sm_poisfft%sm_node_barrier

    CALL cpu_log( log_point_s(33), 'fft_x_1d ', 'continue' )

    jj = nys_x_work

!--  !$OMP  PARALLEL DO PRIVATE ( i, j, k, m, n, work_fftx )
!
!-- Carry out the backward fast Fourier transformation along x.
    DO  j = nys_z, nyn_z
       kw = zstart_work
       DO  k = zstart, zend
!
!--       Collect the data for the backward FFT along x in a 1d-array.
          DO  i = 0, nx
             work_fftx(i) = work_trix(i,j,k)
          ENDDO
!
!--       Backward FFT along x.
          CALL fft_x_1d( work_fftx, 'backward' )
!
!--       Store result in a 4-d array as input for the transposition x --> y as next step.
          m = 0
          IF ( non_uniform_data_for_transpose )  THEN
             DO  n = 1, npex
                DO  i = 1, nnx_pe(n-1)
                   ar(i,kw,jj,n) = work_fftx(m)
                   m = m + 1
                ENDDO
             ENDDO
          ELSE
             DO  n = 1, npex
                DO  i = 1, nnx
                   ar(i,k,j,n) = work_fftx(m)
                   m = m + 1
                ENDDO
             ENDDO
          ENDIF
          kw = kw + 1
       ENDDO
       jj = jj + 1
    ENDDO

    CALL cpu_log( log_point_s(33), 'fft_x_1d ', 'stop' )

 END SUBROUTINE fftx_tri_fftx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition x --> y with subsequent backward FFT along y for a 1d-decomposition along x,
!> based on shared memory.
!> The transposition is done for each k-layer separately, to allow for some overlap of communication
!> and calculation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE tr_xy_ffty

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  ierr  !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  jj    !<
    INTEGER(iwp) ::  js    !<
    INTEGER(iwp) ::  je    !<
    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  kw    !<
    INTEGER(iwp) ::  l     !<

    INTEGER(iwp), DIMENSION(zstart:zend_max) ::  req  !<

    REAL(wp), DIMENSION(0:ny,zstart:zend,nxl:nxr)                         ::  work_ffty  !<
    REAL(wp), DIMENSION(nxl:nxl+nnx_max-1,0:ny_dim,zstart_work:zend_work) ::  work       !<


    req = 0
    kw = zstart_work - 2
!
!-- Carry out transposition and FFT for each vertical layer separately.
    DO  k = zstart-1, zend_max

       kw = kw + 1
!
!--    For the last vertical level transposition is not required because it has been initiated
!--    in the last loop iteration.
       IF ( k <= zend_max-1 )  THEN
!
!--       Initiate the transposition for vertical level k.
          CALL cpu_log( log_point_s(100), 'mpi_ialltoall_1', 'start' )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_IALLTOALL( xz_buf(1,kw+1,nys_x_work,1), 1, xzbuftype, work(nxl,0,kw+1),         &
                              sendrecvcount_xy, MPI_REAL, comm1dx, req(k+1), ierr )
          CALL cpu_log( log_point_s(100), 'mpi_ialltoall_1', 'stop' )
       ENDIF
!
!--    Skip FFT calculation for the first vertical level to allow for an overlap of data
!--    communication and transposition.
       IF ( k == zstart-1 )  CYCLE
!
!--    Wait for completion of data communication before FFT calculation can be started.
       CALL cpu_log( log_point_s(101), 'mpi_ialltoall_2', 'start' )
       CALL MPI_wait (req(k), MPI_STATUS_IGNORE, ierr)
       CALL cpu_log( log_point_s(101), 'mpi_ialltoall_2', 'stop' )
!
!--    Resort the data in a way that y becomes the first index and carry out the
!--    backward fft along y.
       CALL cpu_log( log_point_s(7), 'fft_y_1d', 'continue' )
       IF ( k <= zend )  THEN
!--        !$OMP  PARALLEL DO PRIVATE ( i, j )
          DO  i = nxl, nxr
!
!--          Resort the data.
             IF ( non_uniform_data_for_transpose )  THEN
                jj = 0
                DO  l = 0, npex - 1
                   js = l * nny_x_max
                   je = js + nyn_x_pe(l) - nys_x_pe(l)
                   DO  j = js, je
                      work_ffty(jj,k,i) = work(i,j,kw)
                      jj = jj + 1
                   ENDDO
                ENDDO
             ELSE
                DO  j = 0, ny
                    work_ffty(j,k,i) = work(i,j,kw)
                ENDDO
             ENDIF

!
!--          Backward FFT along y.
             CALL fft_y_1d( work_ffty(:,k,i), 'backward' )

          ENDDO

       ENDIF
       CALL cpu_log( log_point_s(7), 'fft_y_1d', 'pause' )

    ENDDO

!--  !$OMP  PARALLEL DO PRIVATE ( i, j ,k)
!
!-- Resort data to get k as first index. All data along y are present in the shared memory block
!-- (yz_buf).
    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'continue' )
    DO  i = nxl, nxr
       DO  j = 0, ny
          DO  k = zstart,zend
             yz_buf(k,j,i) = work_ffty(j,k,i)
          ENDDO
       ENDDO
    ENDDO
    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'stop' )

 END SUBROUTINE tr_xy_ffty
#endif

 END MODULE poisfft_sm_mod
