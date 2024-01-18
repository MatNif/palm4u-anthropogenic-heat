!> @file poisfft_mod.f90
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
!> Solves the Poisson equation with a 2D spectral method
!>        d^2 p / dx^2 + d^2 p / dy^2 + d^2 p / dz^2 = s
!>
!> Input:
!> real   ar   contains (nnz,nny,nnx) elements of the velocity divergence, starting from (1,nys,nxl)
!>
!> Output:
!> real   ar   contains the solution for perturbation pressure p
!--------------------------------------------------------------------------------------------------!
 MODULE poisfft_mod

#if defined( __parallel )
    USE MPI
#endif

    USE control_parameters,                                                                        &
        ONLY:  message_string,                                                                     &
               temperton_fft_vec

    USE fft_xy,                                                                                    &
        ONLY:  fft_init,                                                                           &
               fft_y,                                                                              &
               fft_x

    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nny,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               ny,                                                                                 &
               nys,                                                                                &
               nyn,                                                                                &
               nz

    USE pegrid,                                                                                    &
        ONLY:  non_uniform_data_for_transpose

    USE transpose_mod,                                                                             &
        ONLY:  nnx_x_max,                                                                          &
               nnx_y_max,                                                                          &
               nnz_x_max,                                                                          &
               nnz_z_max,                                                                          &
               nxl_y,                                                                              &
               nxl_z,                                                                              &
               nxr_x_max,                                                                          &
               nxr_y_max,                                                                          &
               nxr_y,                                                                              &
               nxr_z,                                                                              &
               nx_y_max,                                                                           &
               nys_x,                                                                              &
               nys_z,                                                                              &
               nyn_x,                                                                              &
               nyn_x_max,                                                                          &
               nyn_z,                                                                              &
               nyn_z_max,                                                                          &
               ny_z_max,                                                                           &
               nzb_x,                                                                              &
               nzb_y,                                                                              &
               nzt_x,                                                                              &
               nzt_x_max,                                                                          &
               nzt_y_max,                                                                          &
               nzt_y,                                                                              &
               nz_x_max,                                                                           &
               resort_for_xy, transpose_xy,                                                        &
               resort_for_xz, transpose_xz,                                                        &
               resort_for_yx, transpose_yx,                                                        &
               resort_for_yz, transpose_yz,                                                        &
               resort_for_zx, transpose_zx,                                                        &
               resort_for_zy, transpose_zy

    USE tridia_solver,                                                                             &
        ONLY:  tridia_1dd,                                                                         &
               tridia_init,                                                                        &
               tridia_substi


    IMPLICIT NONE

    LOGICAL, SAVE ::  poisfft_initialized = .FALSE.  !<

    PRIVATE

    PUBLIC  poisfft, poisfft_init

    INTERFACE poisfft
       MODULE PROCEDURE poisfft
    END INTERFACE poisfft

    INTERFACE poisfft_init
       MODULE PROCEDURE poisfft_init
    END INTERFACE poisfft_init


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Setup coefficients for FFT and the tridiagonal solver
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poisfft_init

    IMPLICIT NONE


    CALL fft_init
    CALL tridia_init

    poisfft_initialized = .TRUE.

 END SUBROUTINE poisfft_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Two-dimensional Fourier Transformation in x- and y-direction.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poisfft( ar )

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz

    USE kinds

    USE pegrid

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nz,nys:nyn,nxl:nxr) ::  ar      !<

#define __acc_fft_device ( defined( _OPENACC ) && ( defined ( __cuda_fft ) ) )
#if __acc_fft_device
    !$ACC DECLARE CREATE(ar_inv)
#endif

    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  xy_in    !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  xy_out   !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  yz_in    !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  yz_out   !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  zx_in    !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  zx_out   !<


    IF ( .NOT. poisfft_initialized )  CALL poisfft_init

    CALL cpu_log( log_point_s(3), 'poisfft', 'start' )

#if !__acc_fft_device
    !$ACC UPDATE HOST(ar)
#endif

!
!-- Two-dimensional Fourier Transformation in x- and y-direction.
    ALLOCATE( xy_in(nys_x:nyn_x_max,nzb_x:nzt_x,0:nx_y_max) )
    ALLOCATE( zx_in(nys:nyn,nxl:nxr_x_max,1:nz_x_max) )
    ALLOCATE( zx_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x) )
!
!-- 2d-domain-decomposition or no decomposition (1 PE run).
!-- Transposition z --> x
    CALL cpu_log( log_point_s(5), 'transpo forward', 'start' )
    CALL resort_for_zx( ar, zx_in )
!
!-- In case of temperton_fft_vec, zx_out is bypassed by f_vec_x.
    CALL transpose_zx( zx_in, zx_out )
    CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

    CALL cpu_log( log_point_s(4), 'fft_x', 'start' )
    IF ( temperton_fft_vec  .AND. .NOT.  non_uniform_data_for_transpose)  THEN
!
!--    Vector version outputs a transformed array ar_inv that does not require resorting
!--    (which is done for ar further below).
       CALL fft_x( zx_out, 'forward',  ar_inv=xy_in)
    ELSE
       CALL fft_x( zx_out, 'forward')
    ENDIF
    CALL cpu_log( log_point_s(4), 'fft_x', 'pause' )

!
!-- Transposition x --> y.
    ALLOCATE( xy_out(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) )
    ALLOCATE( yz_in(nxl_y:nxr_y,nzb_y:nzt_y_max,0:ny_z_max) )

    CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
    IF ( .NOT. temperton_fft_vec  .OR.  non_uniform_data_for_transpose )  THEN
       CALL resort_for_xy( zx_out, xy_in )
    ENDIF
    CALL transpose_xy( xy_in, xy_out )
    CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

    CALL cpu_log( log_point_s(7), 'fft_y', 'start' )
    IF ( temperton_fft_vec .AND. .NOT. non_uniform_data_for_transpose)  THEN
!
!--    Input array ar_inv from fft_x can be directly used here.
!--    The output (also in array ar_inv) does not require resorting below.
!--    This is currently only programmed for uniform sundomains.
!--    TODO: Please check performance on NEC to decide if this branch should also be
!--          implemented for nonuniform subdomains.
!--    This branch saves one resort call, the vector version of Temperton-fft is active in
!--    both cases of this IF condition.
       CALL fft_y( xy_out, 'forward', ar_inv = yz_in, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y,    &
                   nxl_y_l = nxl_y, nxr_y_l = nxr_y )
    ELSE
       CALL fft_y( xy_out, 'forward', ar_tr = xy_out, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y,    &
                   nxl_y_l = nxl_y, nxr_y_l = nxr_y )
    ENDIF
    CALL cpu_log( log_point_s(7), 'fft_y', 'pause' )

!
!-- Transposition y --> z.
    ALLOCATE( yz_out(nxl_z:nxr_z,nys_z:nyn_z,1:nz) )

    CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
    IF ( .NOT. temperton_fft_vec  .OR.  non_uniform_data_for_transpose )  THEN
       CALL resort_for_yz( xy_out, yz_in )
    ENDIF
    CALL transpose_yz( yz_in, yz_out )
    CALL cpu_log( log_point_s(5), 'transpo forward', 'stop' )

!
!-- Solve the tridiagonal equation system along z.
    CALL cpu_log( log_point_s(6), 'tridia', 'start' )
    CALL tridia_substi( yz_out, nxl_z, nxr_z )
    CALL cpu_log( log_point_s(6), 'tridia', 'stop' )
!
!
!-- Inverse Fourier Transformation.
!-- Transposition z --> y.
    CALL cpu_log( log_point_s(8), 'transpo invers', 'start' )
    CALL transpose_zy( yz_out, yz_in )
!
!-- The fft_y below (vector branch) can directly process ar_inv (i.e. does not require a resorting).
    IF ( .NOT. temperton_fft_vec  .OR. non_uniform_data_for_transpose )  THEN
       CALL resort_for_zy( yz_in, xy_out )
    ENDIF
    CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

    DEALLOCATE( yz_out )

    CALL cpu_log( log_point_s(7), 'fft_y', 'continue' )
    IF ( temperton_fft_vec  .AND.  .NOT. non_uniform_data_for_transpose )  THEN
!
!--    Output array ar_inv can be used as input to the below fft_x routine without resorting.
       CALL fft_y( xy_out, 'backward', ar_inv = yz_in, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y,   &
                   nxl_y_l = nxl_y, nxr_y_l = nxr_y )
    ELSE
       CALL fft_y( xy_out, 'backward', ar_tr = xy_out, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y,   &
                   nxl_y_l = nxl_y, nxr_y_l = nxr_y )
    ENDIF

    DEALLOCATE( yz_in )

    CALL cpu_log( log_point_s(7), 'fft_y', 'stop' )

!
!-- Transposition y --> x.
    CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
    CALL transpose_yx( xy_out, xy_in )
    IF ( .NOT. temperton_fft_vec  .OR.  non_uniform_data_for_transpose )  THEN
       CALL resort_for_yx( xy_in, zx_out )
    ENDIF
    CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

    DEALLOCATE( xy_out )

    CALL cpu_log( log_point_s(4), 'fft_x', 'continue' )
    IF ( temperton_fft_vec  .AND.  .NOT. non_uniform_data_for_transpose )  THEN
       CALL fft_x( zx_out, 'backward',  ar_inv=xy_in )
    ELSE
       CALL fft_x( zx_out, 'backward' )
    ENDIF

    DEALLOCATE( xy_in )

    CALL cpu_log( log_point_s(4), 'fft_x', 'stop' )

!
!-- Transposition x --> z.
    CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
    CALL transpose_xz( zx_out, zx_in )
    CALL resort_for_xz( zx_in, ar )
    CALL cpu_log( log_point_s(8), 'transpo invers', 'stop' )

    DEALLOCATE( zx_in )
    DEALLOCATE( zx_out )

#if !__acc_fft_device
    !$ACC UPDATE DEVICE(ar)
#endif

    CALL cpu_log( log_point_s(3), 'poisfft', 'stop' )

 END SUBROUTINE poisfft

 END MODULE poisfft_mod
