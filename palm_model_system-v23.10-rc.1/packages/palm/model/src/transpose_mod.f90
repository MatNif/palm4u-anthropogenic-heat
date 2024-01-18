!> @file transpose.f90
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
!> Resorting data for the transposition from x to y. The transposition itself is carried out in
!> transpose_xy.
!--------------------------------------------------------------------------------------------------!

#define __acc_fft_device ( defined( _OPENACC ) && ( defined ( __cuda_fft ) ) )

 MODULE transpose_mod

#if defined( __parallel )
    USE MPI
#endif

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  collective_wait,                                                                    &
               comm2d,                                                                             &
               comm1dx,                                                                            &
               comm1dy,                                                                            &
               ierr,                                                                               &
               non_uniform_data_for_transpose,                                                     &
               npex,                                                                               &
               npey,                                                                               &
               numprocs,                                                                           &
               sendrecvcount_xy,                                                                   &
               sendrecvcount_yz,                                                                   &
               sendrecvcount_zx,                                                                   &
               sendrecvcount_zyd

    IMPLICIT NONE

    INTEGER(iwp) ::  nxl_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nxl_yd  !< internal index bound for transpositions
    INTEGER(iwp) ::  nxl_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nxr_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nxr_yd  !< internal index bound for transpositions
    INTEGER(iwp) ::  nxr_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nyn_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nyn_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nys_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nys_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzb_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzb_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzb_yd  !< internal index bound for transpositions
    INTEGER(iwp) ::  nzt_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzt_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzt_yd  !< internal index bound for transpositions
!
!-- Variables to handle non-uniform subdomains
    INTEGER(iwp) ::  nnx_x_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nnx_y_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nny_yd_max  !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nny_z_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nnz_x_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nnz_yd_max  !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nnz_z_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nxr_x_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nxr_y_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nx_y_max    !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nyn_x_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nyn_yd_max  !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nyn_z_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  ny_z_max    !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nzt_x_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nzt_y_max   !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nzt_yd_max  !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nz_x_max    !< internal index bound for allocating transpose arrays
    INTEGER(iwp) ::  nz_yd_max   !< internal index bound for allocating transpose arrays

    INTEGER, DIMENSION(:), ALLOCATABLE ::  nxl_y_pe   !< index bounds for all PEs for transposition x --> y
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nxr_y_pe   !< index bounds for all PEs for transposition x --> y
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nyn_z_pe   !< index bounds for all PEs for transposition y --> z
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nys_z_pe   !< index bounds for all PEs for transposition y --> z
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nzb_x_pe   !< index bounds for all PEs for transposition z --> x
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nzb_yd_pe  !< index bounds for all PEs for direct transposition z --> y
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nzt_x_pe   !< index bounds for all PEs for transposition z --> x
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nzt_yd_pe  !< index bounds for all PEs for direct transposition z --> y

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f_vec_x  !< array used for vectorized Temperton FFT, allocated in fft_xy

!
!-- Variables for transpose_zx with z dimensioned different than (1:nz)
    INTEGER(iwp) ::  nnx_x_max_diffnz
    INTEGER(iwp) ::  nnz_x_diffnz
    INTEGER(iwp) ::  nnz_x_max_diffnz
    INTEGER(iwp) ::  nzb_x_diffnz
    INTEGER(iwp) ::  nzb_x_max_diffnz
    INTEGER(iwp) ::  nzt_x_diffnz
    INTEGER(iwp) ::  nzt_x_max_diffnz
    INTEGER(iwp) ::  nz_x_max_diffnz
    INTEGER(iwp) ::  sendrecvcount_zx_diffnz

    INTEGER, DIMENSION(:), ALLOCATABLE ::  nzb_x_pe_diffnz  !< index bounds for all PEs for transposition z --> x
    INTEGER, DIMENSION(:), ALLOCATABLE ::  nzt_x_pe_diffnz  !< index bounds for all PEs for transposition z --> x

    SAVE

    INTERFACE resort_for_xy
       MODULE PROCEDURE resort_for_xy
    END INTERFACE resort_for_xy

    INTERFACE transpose_xy
       MODULE PROCEDURE transpose_xy
    END INTERFACE transpose_xy

    INTERFACE resort_for_xz
       MODULE PROCEDURE resort_for_xz
    END INTERFACE resort_for_xz

    INTERFACE transpose_xz
       MODULE PROCEDURE transpose_xz
    END INTERFACE transpose_xz

    INTERFACE resort_for_yx
       MODULE PROCEDURE resort_for_yx
    END INTERFACE resort_for_yx

    INTERFACE transpose_yx
       MODULE PROCEDURE transpose_yx
    END INTERFACE transpose_yx

    INTERFACE resort_for_yz
       MODULE PROCEDURE resort_for_yz
    END INTERFACE resort_for_yz

    INTERFACE transpose_yz
       MODULE PROCEDURE transpose_yz
    END INTERFACE transpose_yz

    INTERFACE resort_for_zx
       MODULE PROCEDURE resort_for_zx_diffnz
       MODULE PROCEDURE resort_for_zx
    END INTERFACE resort_for_zx

    INTERFACE transpose_zx
       MODULE PROCEDURE transpose_zx
#if defined( __parallel )
       MODULE PROCEDURE transpose_zx_diffnz
#endif
    END INTERFACE transpose_zx

#if defined( __parallel )
    INTERFACE setup_transpose_indices_zx_diffnz
       MODULE PROCEDURE setup_transpose_indices_zx_diffnz
    END INTERFACE setup_transpose_indices_zx_diffnz
#endif

    INTERFACE resort_for_zy
       MODULE PROCEDURE resort_for_zy
    END INTERFACE resort_for_zy

    INTERFACE transpose_zy
       MODULE PROCEDURE transpose_zy
    END INTERFACE transpose_zy

#if defined( __parallel )
    INTERFACE transpose_zyd
       MODULE PROCEDURE transpose_zyd
    END INTERFACE transpose_zyd
#endif

    PUBLIC

 CONTAINS


 SUBROUTINE resort_for_xy( f_in, f_inv )

    USE indices,                                                                                   &
        ONLY:  nx

    IMPLICIT NONE

    INTEGER(iwp) ::  i        !<
    INTEGER(iwp) ::  ii       !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  k        !<
    INTEGER(iwp) ::  l        !<

    REAL(wp) ::  f_in(0:nx,nys_x:nyn_x,nzb_x:nzt_x)   !<
    REAL(wp) ::  f_inv(nys_x:nyn_x_max,nzb_x:nzt_x,0:nx_y_max)  !<

!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    IF ( non_uniform_data_for_transpose )  THEN

       !$OMP  PARALLEL PRIVATE ( i, j, k, ii, l )
       !$OMP  DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,ii,l) &
       !$ACC PRESENT(f_inv, f_in)
#endif
       DO  k = nzb_x, nzt_x
          DO  j = nys_x, nyn_x
             ii = 0
             DO  l = 0, npey - 1
                DO  i = l*nnx_y_max, l*nnx_y_max+nxr_y_pe(l)-nxl_y_pe(l)
                   f_inv(j,k,i) = f_in(ii,j,k)
                   ii = ii+1
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !$OMP  END PARALLEL

    ELSE

       !$OMP  PARALLEL PRIVATE ( i, j, k )
       !$OMP  DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_inv, f_in)
#endif
       DO  k = nzb_x, nzt_x
          DO  j = nys_x, nyn_x
             DO  i = 0, nx
                f_inv(j,k,i) = f_in(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       !$OMP  END PARALLEL
    ENDIF

 END SUBROUTINE resort_for_xy


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from x to y. For the input array, all elements along x reside
!> on the same PE, while after transposition, all elements along y reside on the same PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_xy( f_inv, f_out )

#if defined( __parallel )
    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s
#endif
    USE indices,                                                                                   &
        ONLY:  ny
#if defined( __parallel )
    USE indices,                                                                                   &
        ONLY:  nny_pe
#endif

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

#if defined( __parallel )
    INTEGER(iwp) ::  jj !<
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  ys  !<
#endif

    REAL(wp) ::  f_inv(nys_x:nyn_x_max,nzb_x:nzt_x,0:nx_y_max)  !<
    REAL(wp) ::  f_out(0:ny,nxl_y:nxr_y,nzb_y:nzt_y)  !<

#if defined( __parallel )
    REAL(wp), DIMENSION(nyn_x_max-nys_x+1,nzb_y:nzt_y,nxl_y:nxr_y_max,0:npey-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif
#endif


    IF ( numprocs /= 1 )  THEN

#if defined( __parallel )
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE HOST(f_inv)
#else
       !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_inv(nys_x,nzb_x,0),  sendrecvcount_xy, MPI_REAL,                       &
                          work(1,nzb_y,nxl_y,0), sendrecvcount_xy, MPI_REAL, comm1dy, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE DEVICE(work)
#else
       !$ACC END HOST_DATA
#endif
#endif

       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!--    Reorder transposed array.
       IF ( non_uniform_data_for_transpose )  THEN

          !$OMP  PARALLEL PRIVATE ( i, j, k, l, jj )
#if __acc_fft_device
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,jj) &
          !$ACC PRESENT(f_out, work)
#endif
          !$OMP DO
          DO  i = nxl_y, nxr_y
             DO  k = nzb_y, nzt_y
                jj = 0
                DO  l = 0, npey - 1
                   DO  j = 1, nny_pe(l)
                      f_out(jj,i,k) = work(j,k,i,l)
                      jj = jj+1
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

       ELSE

          !$OMP  PARALLEL PRIVATE ( i, j, k, l, ys )
          DO  l = 0, npey - 1
             ys = 0 + l * ( nyn_x - nys_x + 1 )
#if __acc_fft_device
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
             !$ACC PRESENT(f_out, work)
#endif
             !$OMP DO
             DO  i = nxl_y, nxr_y
                DO  k = nzb_y, nzt_y
                   DO  j = ys, ys + nyn_x - nys_x
                      f_out(j,i,k) = work(j-ys+1,k,i,l)
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END DO NOWAIT
          ENDDO
          !$OMP END PARALLEL

       ENDIF
#endif

    ELSE

!
!--    Reorder transposed array
       !$OMP  PARALLEL PRIVATE ( i, j, k )
       !$OMP  DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_out, f_inv)
#endif
       DO  k = nzb_y, nzt_y
          DO  i = nxl_y, nxr_y
             DO  j = 0, ny
                f_out(j,i,k) = f_inv(j,k,i)
             ENDDO
          ENDDO
       ENDDO
       !$OMP  END PARALLEL

    ENDIF

 END SUBROUTINE transpose_xy


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data after the transposition from x to z. The transposition itself is carried out in
!> transpose_xz.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resort_for_xz( f_inv, f_out )

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  kk !<
    INTEGER(iwp) ::  l  !<

    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr_x_max,1:nz_x_max)  !<
    REAL(wp) ::  f_out(1:nz,nys:nyn,nxl:nxr)  !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
!-- In case of parallel fft/transposition, scattered store is faster in backward direction.
    IF ( non_uniform_data_for_transpose )  THEN

       !$OMP PARALLEL PRIVATE ( i, j, k, kk, l )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_in, f_inv)
#endif
       DO  i = nxl, nxr
          DO  j = nys, nyn
             kk = 1
             DO  l = 0, npex - 1
                DO  k = l*nnz_x_max+1,l*nnz_x_max+nzt_x_pe(l)-nzb_x_pe(l)+1
                   f_out(kk,j,i) = f_inv(j,i,k)
                   kk = kk+1
                END DO
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ELSE

       !$OMP  PARALLEL PRIVATE ( i, j, k )
       !$OMP  DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_out, f_inv)
#endif
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 1, nz
                f_out(k,j,i) = f_inv(j,i,k)
             ENDDO
          ENDDO
       ENDDO
       !$OMP  END PARALLEL

    ENDIF

 END SUBROUTINE resort_for_xz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from x to z. For the input array, all elements along x reside
!> on the same PE, while after transposition, all elements along z reside on the same PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_xz( f_in, f_inv )

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  temperton_fft_vec

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s
#endif

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz

#if defined( __parallel )
    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nxl_pe,                                                                             &
               nxr_pe
#endif

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

#if defined( __parallel )
    INTEGER(iwp) ::  ii  !<
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  mm  !<‚
    INTEGER(iwp) ::  xs  !<
#endif

    REAL(wp) ::  f_in(0:nx,nys_x:nyn_x,nzb_x:nzt_x)  !<
    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr_x_max,1:nz_x_max)          !<


#if defined( __parallel )
    REAL(wp), DIMENSION(nys_x:nyn_x,nnx_x_max,nzb_x:nzt_x_max,0:npex-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif
#endif


!
!-- If the PE grid is one-dimensional along y, the array has only to be reordered locally and
!-- therefore no transposition has to be done.
    IF ( npex /= 1 )  THEN

#if defined( __parallel )
!
!--    Reorder input array for transposition. Data from the vectorized Temperton-fft is stored in
!--    different array format (f_vec_x).
       IF ( temperton_fft_vec )  THEN

          DO  l = 0, npex - 1
             xs = 0 + l * nnx
             DO  k = nzb_x, nzt_x
                IF ( non_uniform_data_for_transpose )  THEN
                   ii = 1
                   DO  i = nxl_pe(l), nxr_pe(l)
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         work(j,ii,k,l) = f_vec_x(mm,i)
                      ENDDO
                      ii = ii+1
                   ENDDO
                ELSE
                   DO  i = xs, xs + nnx - 1
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         work(j,i-xs+1,k,l) = f_vec_x(mm,i)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

       ELSE

          !$OMP  PARALLEL PRIVATE ( i, j, k, l, xs, ii )
          DO  l = 0, npex - 1
             xs = 0 + l * nnx
#if __acc_fft_device
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,ii) &
             !$ACC PRESENT(work, f_in)
#endif
             !$OMP DO
             DO  k = nzb_x, nzt_x
                IF ( non_uniform_data_for_transpose )  THEN
                   ii = 1
                   DO  i = nxl_pe(l), nxr_pe(l)
                      DO  j = nys_x, nyn_x
                         work(j,ii,k,l) = f_in(i,j,k)
                      ENDDO
                      ii = ii+1
                   ENDDO
                ELSE
                   DO  i = xs, xs + nnx - 1
                      DO  j = nys_x, nyn_x
                         work(j,i-xs+1,k,l) = f_in(i,j,k)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             !$OMP END DO NOWAIT
          ENDDO
          !$OMP  END PARALLEL

       ENDIF
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE HOST(work)
#else
       !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(nys_x,1,nzb_x,0), sendrecvcount_zx, MPI_REAL,                       &
                          f_inv(nys,nxl,1),      sendrecvcount_zx, MPI_REAL, comm1dx, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE DEVICE(f_inv)
#else
       !$ACC END HOST_DATA
#endif
#endif
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    ELSE

!
!--    Reorder the array in a way that the z index is in first position
       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_inv, f_in)
#endif
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 1, nz
                f_inv(j,i,k) = f_in(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ENDIF

 END SUBROUTINE transpose_xz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data after the transposition from y to x. The transposition itself is carried out in
!> transpose_yx.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resort_for_yx( f_inv, f_out )

    USE indices,                                                                                   &
        ONLY:  nx

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  ii !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<

    REAL(wp) ::  f_inv(nys_x:nyn_x_max,nzb_x:nzt_x,0:nx_y_max)  !<
    REAL(wp) ::  f_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x)                    !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    IF ( non_uniform_data_for_transpose )  THEN

       !$OMP  PARALLEL PRIVATE ( i, j, k, ii, l )
       !$OMP  DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,ii,l) &
       !$ACC PRESENT(f_inv, f_in)
#endif
       DO  k = nzb_x, nzt_x
          DO  j = nys_x, nyn_x
             ii = 0
             DO  l = 0, npey - 1
                DO  i = l*nnx_y_max, l*nnx_y_max+nxr_y_pe(l)-nxl_y_pe(l)
                   f_out(ii,j,k) = f_inv(j,k,i)
                   ii = ii+1
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !$OMP  END PARALLEL

    ELSE

       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_out, f_inv)
#endif
       DO  k = nzb_x, nzt_x
          DO  j = nys_x, nyn_x
             DO  i = 0, nx
                f_out(i,j,k) = f_inv(j,k,i)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ENDIF

 END SUBROUTINE resort_for_yx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from y to x. For the input array, all  elements along y
!> reside on the same PE, while after transposition, all elements along x reside on the same PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_yx( f_in, f_inv )

#if defined( __parallel )
    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s
#endif
    USE indices,                                                                                   &
        ONLY:  ny
#if defined( __parallel )
    USE indices,                                                                                   &
        ONLY:  nny_pe
#endif

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

#if defined( __parallel )
    INTEGER(iwp) ::  jj !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  ys  !<
#endif

    REAL(wp) ::  f_in(0:ny,nxl_y:nxr_y,nzb_y:nzt_y)             !<
    REAL(wp) ::  f_inv(nys_x:nyn_x_max,nzb_x:nzt_x,0:nx_y_max)  !<


#if defined( __parallel )
    REAL(wp), DIMENSION(nyn_x_max-nys_x+1,nzb_y:nzt_y,nxl_y:nxr_y_max,0:npey-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif
#endif


    IF ( numprocs /= 1 )  THEN

#if defined( __parallel )
!
!--    Reorder input array for transposition
       IF ( non_uniform_data_for_transpose )  THEN

          !$OMP PARALLEL PRIVATE ( i, j, k, l, jj)
#if __acc_fft_device
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,l,jj) &
          !$ACC PRESENT(work, f_in)
#endif
          !$OMP DO
          DO  i = nxl_y, nxr_y
             DO  k = nzb_y, nzt_y
                jj = 0
                DO  l = 0, npey - 1
                   ys = 0 + l * ( nyn_x - nys_x + 1 )
                   DO  j = 1, nny_pe(l)
                      work(j,k,i,l )= f_in(jj,i,k)
                      jj = jj+1
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

       ELSE

          !$OMP PARALLEL PRIVATE ( i, j, k, l, ys )
          DO  l = 0, npey - 1
             ys = 0 + l * ( nyn_x - nys_x + 1 )
#if __acc_fft_device
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
             !$ACC PRESENT(work, f_in)
#endif
             !$OMP DO
             DO  i = nxl_y, nxr_y
                DO  k = nzb_y, nzt_y
                   DO  j = ys, ys + nyn_x - nys_x
                      work(j-ys+1,k,i,l) = f_in(j,i,k)
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END DO NOWAIT
          ENDDO
          !$OMP END PARALLEL

       ENDIF

!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE HOST(work)
#else
       !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(1,nzb_y,nxl_y,0), sendrecvcount_xy, MPI_REAL,                       &
                          f_inv(nys_x,nzb_x,0),  sendrecvcount_xy, MPI_REAL, comm1dy, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE DEVICE(f_inv)
#else
       !$ACC END HOST_DATA
#endif
#endif
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    ELSE

!
!--    Reorder array f_in the same way as ALLTOALL did it.
       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_inv, f_in)
#endif
       DO  i = nxl_y, nxr_y
          DO  k = nzb_y, nzt_y
             DO  j = 0, ny
                f_inv(j,k,i) = f_in(j,i,k)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ENDIF

 END SUBROUTINE transpose_yx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data for the transposition from y to z. The transposition itself is carried out in
!> transpose_yz.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resort_for_yz( f_in, f_inv )

    USE indices,                                                                                   &
        ONLY:  ny

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  jj !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<

    REAL(wp) ::  f_in(0:ny,nxl_y:nxr_y,nzb_y:nzt_y)                     !<
    REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y_max,0:ny_z_max)  !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    !$OMP PARALLEL PRIVATE ( i, j, k, l, jj )
    !$OMP DO
#if __acc_fft_device
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,l,jj) &
    !$ACC PRESENT(f_inv, f_in)
#endif
    DO  k = nzb_y, nzt_y
       DO  i = nxl_y, nxr_y
          IF ( non_uniform_data_for_transpose )  THEN
             jj = 0
             DO  l = 0, npex - 1
                DO  j = l*nny_z_max, l*nny_z_max+nyn_z_pe(l)-nys_z_pe(l)
                   f_inv(i,k,j) = f_in(jj,i,k)
                   jj = jj+1
                ENDDO
             ENDDO
          ELSE
             DO  j = 0, ny
                f_inv(i,k,j) = f_in(j,i,k)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !$OMP END PARALLEL

 END SUBROUTINE resort_for_yz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from y to z. For the input array, all elements along y reside
!> on the same PE, while after transposition, all elements along z reside on the same PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_yz( f_inv, f_out )

#if defined( __parallel )
    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s
#endif

    USE indices,                                                                                   &
        ONLY:  ny,                                                                                 &
               nz

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

#if defined( __parallel )
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  kk !<
    INTEGER(iwp) ::  zs  !<
#endif

    REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y_max,0:ny_z_max)  !<
    REAL(wp) ::  f_out(nxl_z:nxr_z,nys_z:nyn_z,1:nz)                    !<

#if defined( __parallel )
    REAL(wp), DIMENSION(nxl_z:nxr_z,nnz_z_max,nys_z:nyn_z_max,0:npex-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif
#endif


!
!-- If the PE grid is one-dimensional along y, only local reordering of the data is necessary and no
!-- transposition has to be done.
    IF ( npex == 1 )  THEN

       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_out, f_inv)
#endif
       DO  j = 0, ny
          DO  k = nzb_y, nzt_y
             DO  i = nxl_y, nxr_y
                f_out(i,j,k) = f_inv(i,k,j)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ELSE

#if defined( __parallel )
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE HOST(f_inv)
#else
       !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_inv(nxl_y,nzb_y,0),  sendrecvcount_yz, MPI_REAL,                       &
                          work(nxl_z,1,nys_z,0), sendrecvcount_yz, MPI_REAL, comm1dx, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE DEVICE(work)
#else
       !$ACC END HOST_DATA
#endif
#endif

       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!--    Reorder transposed array
       IF ( non_uniform_data_for_transpose )  THEN

          !$OMP PARALLEL PRIVATE ( i, j, k, l, kk )
#if __acc_fft_device
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,l,kk) &
          !$ACC PRESENT(f_out, work)
#endif
          !$OMP DO
          DO  j = nys_z, nyn_z
             kk = 1
             DO  l = 0, npex - 1
                DO  k = 1, nzt_x_pe(l)-nzb_x_pe(l)+1
                   DO  i = nxl_z, nxr_z
                      f_out(i,j,kk) = work(i,k,j,l)
                   ENDDO
                   kk = kk+1
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

       ELSE

          !$OMP PARALLEL PRIVATE ( i, j, k, l, zs )
          DO  l = 0, npex - 1
             zs = 1 + l * ( nzt_y - nzb_y + 1 )
#if __acc_fft_device
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
             !$ACC PRESENT(f_out, work)
#endif
             !$OMP DO
             DO  j = nys_z, nyn_z
                DO  k = zs, zs + nzt_y - nzb_y
                   DO  i = nxl_z, nxr_z
                      f_out(i,j,k) = work(i,k-zs+1,j,l)
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END DO NOWAIT
          ENDDO
          !$OMP END PARALLEL

       ENDIF
#endif

   ENDIF

 END SUBROUTINE transpose_yz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data for the transposition from z to x. The transposition itself is carried out in
!> transpose_zx.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resort_for_zx( f_in, f_inv )

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  kk !<
    INTEGER(iwp) ::  l  !<

    REAL(wp) ::  f_in(1:nz,nys:nyn,nxl:nxr)   !<
    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr_x_max,1:nz_x_max)  !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    IF ( non_uniform_data_for_transpose )  THEN

       !$OMP PARALLEL PRIVATE ( i, j, k, kk, l )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,kk,l)   &
       !$ACC PRESENT(f_in, f_inv)
#endif
       DO  i = nxl, nxr
          DO  j = nys, nyn
             kk = 1
             DO  l = 0, npex - 1
                DO  k = l*nnz_x_max+1, l*nnz_x_max+nzt_x_pe(l)-nzb_x_pe(l)+1
                   f_inv(j,i,k) = f_in(kk,j,i)
                   kk = kk+1
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ELSE

       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_in, f_inv)
#endif
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 1,nz
                f_inv(j,i,k) = f_in(k,j,i)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ENDIF

 END SUBROUTINE resort_for_zx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data for the transposition from z to x. The transposition itself is carried out in
!> transpose_zx.
!> This routine assumes that subdomains are non-uniform, which includes the case of uniform
!> subdomains.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resort_for_zx_diffnz( f_in, f_inv, nzb_do, nzt_do )

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN)  ::  nzb_do
    INTEGER(iwp), INTENT(IN)  ::  nzt_do
    REAL(wp), INTENT(IN), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do)        ::  f_in   !<
    REAL(sp), INTENT(OUT),DIMENSION(nys:nyn,nxl:nxr_x_max,1:nz_x_max_diffnz) ::  f_inv  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  kk !<
    INTEGER(iwp) ::  l  !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    !$OMP PARALLEL PRIVATE ( i, j, k, kk, l )
    !$OMP DO
#if __acc_fft_device
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,kk,l)   &
    !$ACC PRESENT(f_in, f_inv)
#endif
    DO  i = nxl, nxr
       DO  j = nys, nyn
          kk = nzb_do
          DO  l = 0, npex - 1
             DO  k = l * nnz_x_max_diffnz + 1,                                                     &
                     l * nnz_x_max_diffnz + nzt_x_pe_diffnz(l) - nzb_x_pe_diffnz(l) + 1
                f_inv(j,i,k) = f_in(i,j,kk)
                kk = kk+1
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL

 END SUBROUTINE resort_for_zx_diffnz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from z to x. For the input array, all elements along z reside
!> on the same PE, while after transposition, all elements along x reside on the same PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_zx( f_inv, f_out )

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               temperton_fft_vec

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s

#endif

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz

#if defined( __parallel )
    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nxl_pe,                                                                             &
               nxr_pe
#endif

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

#if defined( __parallel )
    INTEGER(iwp) ::  ii  !<
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  mm  !<
    INTEGER(iwp) ::  xs  !<
#endif

    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr_x_max,1:nz_x_max)          !<
    REAL(wp) ::  f_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x)  !<

#if defined( __parallel )
    REAL(wp), DIMENSION(nys_x:nyn_x,nnx_x_max,nzb_x:nzt_x_max,0:npex-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif
#endif


!
!-- If the PE grid is one-dimensional along y, only local reordering of the data is necessary and no
!-- transposition has to be done.
    IF ( npex == 1 )  THEN

       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_out, f_inv)
#endif
       DO  k = 1, nz
          DO  i = nxl, nxr
             DO  j = nys, nyn
                f_out(i,j,k) = f_inv(j,i,k)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ELSE

#if defined( __parallel )
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE HOST(f_inv)
#else
       !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_inv(nys,nxl,1),      sendrecvcount_zx, MPI_REAL,                       &
                          work(nys_x,1,nzb_x,0), sendrecvcount_zx, MPI_REAL, comm1dx, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE DEVICE(work)
#else
       !$ACC END HOST_DATA
#endif
#endif

       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!--    Reorder transposed array.
!--    Data for the vectorized Temperton-fft is stored in different array format (f_vec_x) which
!--    saves additional data copy in fft_x.
!--    INFO: Using f_vec_x for data transfer does not work for temperton_fft_vec,
!--          non_uniform_data_for_transpose, and lateral Neumann boundary conditions.
!--          With these conditions, data transfer via f_out is activated as workaround.
!--          The problem only was only observed with gfortran, with Intel compiler transfer
!--          via f_vec_x works!
       IF ( temperton_fft_vec  .AND.  bc_lr_cyc )  THEN

          DO  l = 0, npex - 1
             xs = 0 + l * nnx
             DO  k = nzb_x, nzt_x
                IF ( non_uniform_data_for_transpose )  THEN
                   ii = 1
                   DO  i = nxl_pe(l), nxr_pe(l)
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         f_vec_x(mm,i) = work(j,ii,k,l)
                      ENDDO
                      ii = ii+1
                   ENDDO
                ELSE
                   DO  i = xs, xs + nnx - 1
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         f_vec_x(mm,i) = work(j,i-xs+1,k,l)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

       ELSE

          !$OMP PARALLEL PRIVATE ( i, j, k, l, xs, ii )
          DO  l = 0, npex - 1
             xs = 0 + l * nnx
#if __acc_fft_device
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,ii) &
             !$ACC PRESENT(f_out, work)
#endif
             !$OMP DO
             DO  k = nzb_x, nzt_x
                IF ( non_uniform_data_for_transpose )  THEN
                   ii = 1
                   DO  i = nxl_pe(l), nxr_pe(l)
                      DO  j = nys_x, nyn_x
                         f_out(i,j,k) = work(j,ii,k,l)
                      ENDDO
                      ii = ii+1
                   ENDDO
                ELSE
                   DO  i = xs, xs + nnx - 1
                      DO  j = nys_x, nyn_x
                         f_out(i,j,k) = work(j,i-xs+1,k,l)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             !$OMP END DO NOWAIT
          ENDDO
          !$OMP END PARALLEL

       ENDIF
#endif

    ENDIF

 END SUBROUTINE transpose_zx


#if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from z to x. For the input array, all elements along z reside
!> on the same PE, while after transposition, all elements along x reside on the same PE.
!> This version is for an arbitrary z-dimension (nzb_do:nzt_do) and assumes arrays to be of single
!> precision.
!> This routine assumes that subdomains are non-uniform, which includes the case of uniform
!> subdomains.
!> Attention: Routine does not work for a 1d-decomposition along y (npex=1).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_zx_diffnz( f_inv, f_out, diffnz )

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nyn,                                                                                &
               nys

    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nxl_pe,                                                                             &
               nxr_pe

    IMPLICIT NONE

    LOGICAL, INTENT(IN) ::  diffnz  !< dummy parameter to decide for this routine (overloading)

    REAL(sp), INTENT(IN)  ::  f_inv(nys:nyn,nxl:nxr_x_max,1:nz_x_max_diffnz)     !<
    REAL(sp), INTENT(OUT) ::  f_out(0:nx,nys_x:nyn_x,nzb_x_diffnz:nzt_x_diffnz)  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    INTEGER(iwp) ::  ii  !<
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  xs  !<

    REAL(sp), DIMENSION(nys_x:nyn_x,nnx_x_max,nzb_x_diffnz:nzt_x_max_diffnz,0:npex-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif


!
!-- Dummy statement to avoid compiler warning.
    IF ( diffnz )  CONTINUE

!
!-- Transpose array
    CALL cpu_log( log_point_s(65), 'mpi_alltoall (nzt_do)', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
    !$ACC UPDATE HOST(f_inv)
#else
    !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif

    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
!
!-- 32 Bit version. mpi_real MUST be in small letters here, because MPI_REAL in big letter would
!-- be substituted by MPI_DOUBLE_PRECISION by cpp preprocessor.
    CALL MPI_ALLTOALL( f_inv(nys,nxl,1),             sendrecvcount_zx_diffnz, mpi_real,         &
                       work(nys_x,1,nzb_x_diffnz,0), sendrecvcount_zx_diffnz, mpi_real,         &
                       comm1dx, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
    !$ACC UPDATE DEVICE(work)
#else
    !$ACC END HOST_DATA
#endif
#endif

    CALL cpu_log( log_point_s(65), 'mpi_alltoall (nzt_do)', 'stop' )

    !$OMP PARALLEL PRIVATE ( i, j, k, l, xs, ii )
    DO  l = 0, npex - 1
       xs = 0 + l * nnx
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,ii) &
       !$ACC PRESENT(f_out, work)
#endif
       !$OMP DO
       DO  k = nzb_x_diffnz, nzt_x_diffnz
          ii = 1
          DO  i = nxl_pe(l), nxr_pe(l)
             DO  j = nys_x, nyn_x
                f_out(i,j,k) = work(j,ii,k,l)
             ENDDO
             ii = ii+1
          ENDDO
       ENDDO
       !$OMP END DO NOWAIT
    ENDDO
    !$OMP END PARALLEL

 END SUBROUTINE transpose_zx_diffnz
#endif


#if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Setup the transposition indices for the transposition from z to x (transpose_zx_diffnz) with
!> arbitrary z-dimension (nzb_do:nzt_do).
!> This routine assumes that subdomains are non-uniform, which includes the case of uniform
!> subdomains.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE setup_transpose_indices_zx_diffnz( nzb_do, nzt_do )

    USE indices,                                                                                   &
        ONLY:  nnx_pe,                                                                             &
               nny

    USE pegrid,                                                                                    &
        ONLY:  myidx

    INTEGER(iwp), INTENT(IN) ::  nzb_do  !< lower index along z
    INTEGER(iwp), INTENT(IN) ::  nzt_do  !< upper index along z

    INTEGER(iwp) ::  irest_z
    INTEGER(iwp) ::  k
    INTEGER(iwp) ::  nz_do


    IF ( .NOT. ALLOCATED( nzb_x_pe_diffnz ) )  ALLOCATE( nzb_x_pe_diffnz(0:npex-1) )
    IF ( .NOT. ALLOCATED( nzt_x_pe_diffnz ) )  ALLOCATE( nzt_x_pe_diffnz(0:npex-1) )

!
!-- Always assume that subdomains may be non-uniform. nnx then must have the value of the largest
!-- subdomain, which is the one on PE0 anyway.
    nnx_x_max_diffnz = nnx_pe(0)

!
!-- nz_do is normally not equal to nz (i.e. the z-axis is not dimensioned 1:nzt).
    nz_do        = nzt_do - nzb_do + 1
    nnz_x_diffnz = nz_do / npex
    irest_z      = MOD( nz_do, npex )

!
!-- Always calculate for non-uniform subdomains (k dimension is almost never 1:nzt).
    DO  k = 0, npex-1
       IF ( k < irest_z )  THEN
          nzb_x_pe_diffnz(k) = k * ( nnz_x_diffnz + 1 ) + 1
          nzt_x_pe_diffnz(k) = nzb_x_pe_diffnz(k) + ( nnz_x_diffnz + 1 ) - 1
       ELSE
          nzb_x_pe_diffnz(k) = irest_z * ( nnz_x_diffnz + 1 ) + ( k - irest_z ) * nnz_x_diffnz + 1
          nzt_x_pe_diffnz(k) = nzb_x_pe_diffnz(k) + nnz_x_diffnz - 1
       ENDIF
    ENDDO

    nzb_x_diffnz = nzb_x_pe_diffnz(myidx)
    nzt_x_diffnz = nzt_x_pe_diffnz(myidx)
    nnz_x_diffnz = nzt_x_diffnz - nzb_x_diffnz + 1

!
!-- Calculate index bounds and number of grid points that are required for allocating the
!-- 3d-arrays used by the transpose routine. Arrays on all PEs need to have identical size.
!-- Therefore, the dimensions defined on PE0 are used, because in case of non-uniform
!-- subdomains, at least on PE0 the respective arrays contain one more grid point than on the
!-- other PEs.
    nnz_x_max_diffnz = nzt_x_pe_diffnz(0) - nzb_x_pe_diffnz(0) + 1
    nz_x_max_diffnz  = npex * nnz_x_max_diffnz
    nzt_x_max_diffnz = nzb_x_diffnz + nnz_x_max_diffnz - 1

    sendrecvcount_zx_diffnz = nnx_x_max_diffnz * nny * nnz_x_max_diffnz

 END SUBROUTINE setup_transpose_indices_zx_diffnz
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data after the transposition from z to y. The transposition itself is carried out in
!> transpose_zy.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resort_for_zy( f_inv, f_out )


    USE indices,                                                                                   &
        ONLY:  ny

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  jj !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<

    REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y_max,0:ny_z_max)  !<
    REAL(wp) ::  f_out(0:ny,nxl_y:nxr_y,nzb_y:nzt_y)  !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    !$OMP PARALLEL PRIVATE ( i, j, k, l, jj )
    !$OMP DO
#if __acc_fft_device
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,l,jj) &
    !$ACC PRESENT(f_out, f_inv)
#endif
    DO  k = nzb_y, nzt_y
       DO  i = nxl_y, nxr_y
          IF ( non_uniform_data_for_transpose )  THEN
             jj = 0
             DO  l = 0, npex - 1
                DO  j = l*nny_z_max, l*nny_z_max+nyn_z_pe(l)-nys_z_pe(l)
                   f_out(jj,i,k) = f_inv(i,k,j)
                   jj = jj+1
                ENDDO
             ENDDO
          ELSE
             DO  j = 0, ny
                f_out(j,i,k) = f_inv(i,k,j)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !$OMP END PARALLEL

 END SUBROUTINE resort_for_zy


!--------------------------------------------------------------------------------------------------!
! Description:cpu_log_nowait
! ------------
!> Transposition of input array (f_in) from z to y. For the input array, all elements along z reside
!> on the same PE, while after transposition, all elements along y reside on the same PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_zy( f_in, f_inv )

#if defined( __parallel )
    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               cpu_log_nowait,                                                                     &
               log_point_s
#endif

    USE indices,                                                                                   &
        ONLY:  ny,                                                                                 &
               nz

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

#if defined( __parallel )
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  kk
    INTEGER(iwp) ::  zs  !<
#endif

    REAL(wp) ::  f_in(nxl_z:nxr_z,nys_z:nyn_z,1:nz)   !<
    REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y_max,0:ny_z_max)  !<

#if defined( __parallel )
    REAL(wp), DIMENSION(nxl_z:nxr_z,nnz_z_max,nys_z:nyn_z_max,0:npex-1) ::  work  !<
#if __acc_fft_device
    !$ACC DECLARE CREATE(work)
#endif
#endif


!
!-- If the PE grid is one-dimensional along y, the array has only to be reordered locally and
!-- therefore no transposition has to be done.
    IF ( npex /= 1 )  THEN

#if defined( __parallel )
!
!--    Reorder input array for transposition.
       IF ( non_uniform_data_for_transpose )  THEN

          !$OMP PARALLEL PRIVATE ( i, j, k, l, kk )
#if __acc_fft_device
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k,l,kk) &
          !$ACC PRESENT(f_out, work)
#endif
          !$OMP DO
          DO  j = nys_z, nyn_z
             kk = 1
             DO  l = 0, npex - 1
                DO  k = 1, nzt_x_pe(l) - nzb_x_pe(l) + 1
                   DO  i = nxl_z, nxr_z
                      work(i,k,j,l) = f_in(i,j,kk)
                   ENDDO
                   kk = kk+1
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

       ELSE

          !$OMP PARALLEL PRIVATE ( i, j, k, l, zs )
          DO  l = 0, npex - 1
             zs = 1 + l * ( nzt_y - nzb_y + 1 )
#if __acc_fft_device
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
             !$ACC PRESENT(work, f_in)
#endif
             !$OMP DO
             DO  j = nys_z, nyn_z
                DO  k = zs, zs + nzt_y - nzb_y
                   DO  i = nxl_z, nxr_z
                      work(i,k-zs+1,j,l) = f_in(i,j,k)
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END DO NOWAIT
          ENDDO
          !$OMP END PARALLEL

       ENDIF
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE HOST(work)
#else
       !$ACC HOST_DATA USE_DEVICE(work, f_inv)
#endif
#endif

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(nxl_z,1,nys_z,0), sendrecvcount_yz, MPI_REAL,                       &
                          f_inv(nxl_y,nzb_y,0),  sendrecvcount_yz, MPI_REAL, comm1dx, ierr )

#if __acc_fft_device
#ifndef __cuda_aware_mpi
       !$ACC UPDATE DEVICE(f_inv)
#else
       !$ACC END HOST_DATA
#endif
#endif

       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    ELSE
!
!--    Reorder the array in the same way like ALLTOALL did it
       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
#if __acc_fft_device
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC PRESENT(f_inv, f_in)
#endif
       DO  k = nzb_y, nzt_y
          DO  j = 0, ny
             DO  i = nxl_y, nxr_y
                f_inv(i,k,j) = f_in(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ENDIF

 END SUBROUTINE transpose_zy


#if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from z to y. For the input array, all elements along z reside
!> on the same PE, while after transposition, all elements along y reside on the same PE. This is a
!> direct transposition for arrays with indices in regular order (k,j,i) (cf. transpose_zy).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transpose_zyd( f_in, f_out )

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nny,                                                                                &
               nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nyn_pe,                                                                             &
               nys,                                                                                &
               nys_pe,                                                                             &
               ny,                                                                                 &
               nz

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  j   !<
    INTEGER(iwp) ::  jj  !<
    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  kk  !<
    INTEGER(iwp) ::  l   !<
    INTEGER(iwp) ::  ys  !<

    REAL(wp) ::  f_in(1:nz,nys:nyn,nxl:nxr)                            !<
    REAL(wp) ::  f_out(0:ny,nxl_yd:nxr_yd,nzb_yd:nzt_yd)               !<

    REAL(wp), DIMENSION(nys:nyn_yd_max,nxl_yd:nxr_yd,1:nz_yd_max)      ::  f_inv  !<
    REAL(wp), DIMENSION(nny_yd_max,nxl:nxr,nzb_yd:nzt_yd_max,0:npey-1) ::  work   !<


!
!-- Rearrange indices of input array in order to make data to be send by MPI contiguous.
    IF ( non_uniform_data_for_transpose )  THEN

       DO  i = nxl, nxr
          DO  j = nys, nyn
             kk = 1
             DO  l = 0, npey - 1
                DO  k = l*nnz_yd_max+1, l*nnz_yd_max+nzt_yd_pe(l)-nzb_yd_pe(l)+1
                   f_inv(j,i,k) = f_in(kk,j,i)
                   kk = kk+1
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSE

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 1, nz
                f_inv(j,i,k) = f_in(k,j,i)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- Move data to different array, because memory location of work1 is needed further below
!-- (work1 = work2). If the PE grid is one-dimensional along x, only local reordering of the data is
!-- necessary and no transposition has to be done.
    IF ( npey == 1 )  THEN
       DO  k = 1, nz
          DO  i = nxl, nxr
             DO  j = nys, nyn
                f_out(j,i,k) = f_inv(j,i,k)
             ENDDO
          ENDDO
       ENDDO
       RETURN
    ENDIF

!
!-- Transpose array
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( f_inv(nys,nxl,1), sendrecvcount_zyd, MPI_REAL,                              &
       work(1,nxl,nzb_yd,0), sendrecvcount_zyd, MPI_REAL, comm1dy, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!-- Reorder transposed array
    DO  l = 0, npey - 1
       ys = 0 + l * nny
       DO  k = nzb_yd, nzt_yd
          DO  i = nxl,nxr
             IF ( non_uniform_data_for_transpose )  THEN
                jj = 1
                DO  j = nys_pe(l), nyn_pe(l)
                   f_out(j,i,k) = work(jj,i,k,l)
                   jj = jj+1
                ENDDO
             ELSE
                DO  j = ys, ys + nny - 1
                   f_out(j,i,k) = work(j-ys+1,i,k,l)
                ENDDO
             end if
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE transpose_zyd
#endif

 END MODULE transpose_mod
