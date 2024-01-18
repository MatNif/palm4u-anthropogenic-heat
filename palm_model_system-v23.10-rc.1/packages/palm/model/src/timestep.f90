!> @file timestep.f90
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
!> Compute the time step under consideration of the FCL and diffusion criterion.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE timestep

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzu,                                                                                &
               dzw,                                                                                &
               kh,                                                                                 &
               km,                                                                                 &
               u,                                                                                  &
               u_stokes_zu,                                                                        &
               v,                                                                                  &
               v_stokes_zu,                                                                        &
               w

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  dt_precipitation,                                                                   &
               bulk_cloud_model

    USE control_parameters,                                                                        &
        ONLY:  cfl_factor,                                                                         &
               dt_3d,                                                                              &
               dt_fixed,                                                                           &
               dt_max,                                                                             &
               galilei_transformation,                                                             &
               message_string,                                                                     &
               nested_run,                                                                         &
               rans_mode,                                                                          &
               stop_dt,                                                                            &
               timestep_reason,                                                                    &
               u_gtrans,                                                                           &
               use_ug_for_galilei_tr,                                                              &
               v_gtrans


    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dx2,                                                                                &
               dy,                                                                                 &
               dy2

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzt

    USE interfaces

    USE kinds

    USE land_surface_model_mod,                                                                    &
        ONLY:  dt_lsm

    USE module_interface,                                                                          &
        ONLY:  module_interface_timestep

    USE pegrid

    USE statistics,                                                                                &
        ONLY:  flow_statistics_called,                                                             &
               hom,                                                                                &
               u_max,                                                                              &
               u_max_ijk,                                                                          &
               v_max,                                                                              &
               v_max_ijk,                                                                          &
               w_max,                                                                              &
               w_max_ijk

    USE urban_surface_mod,                                                                         &
        ONLY:  dt_usm

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  km_max_ijk(3) = -1  !< index values (i,j,k) of location where km_max occurs
    INTEGER(iwp) ::  kh_max_ijk(3) = -1  !< index values (i,j,k) of location where kh_max occurs

    LOGICAL ::  stop_dt_local  !< local switch for controlling the time stepping

    REAL(wp) ::  div         !<
    REAL(wp) ::  dt_diff     !<
    REAL(wp) ::  dt_diff_l   !<
    REAL(wp) ::  dt_u        !<
    REAL(wp) ::  dt_u_l      !<
    REAL(wp) ::  dt_v        !<
    REAL(wp) ::  dt_v_l      !<
    REAL(wp) ::  dt_w        !<
    REAL(wp) ::  dt_w_l      !<
    REAL(wp) ::  km_max      !< maximum of Km in entire domain
    REAL(wp) ::  kh_max      !< maximum of Kh in entire domain
    REAL(wp) ::  u_gtrans_l  !<
    REAL(wp) ::  v_gtrans_l  !<

    REAL(wp), DIMENSION(2)         ::  uv_gtrans_l  !<
#if defined( __parallel )
    REAL(wp), DIMENSION(2)         ::  uv_gtrans    !<
    REAL(wp), DIMENSION(3)         ::  reduce       !<
    REAL(wp), DIMENSION(3)         ::  reduce_l     !<
#endif
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dxyz2_min    !<
    !$ACC DECLARE CREATE(dxyz2_min)


    CALL cpu_log( log_point(12), 'calculate_timestep', 'start' )

!
!-- In case of Galilei-transform not using the geostrophic wind as translation velocity, compute the
!-- volume-averaged horizontal velocity components, which will then be subtracted from the
!-- horizontal wind for the time step and horizontal advection routines.
    IF ( galilei_transformation  .AND. .NOT.  use_ug_for_galilei_tr )  THEN
       IF ( flow_statistics_called )  THEN
!
!--       Horizontal averages already existent, just need to average them vertically.
          u_gtrans = 0.0_wp
          v_gtrans = 0.0_wp
          DO  k = nzb+1, nzt
             u_gtrans = u_gtrans + hom(k,1,1,0)
             v_gtrans = v_gtrans + hom(k,1,2,0)
          ENDDO
          u_gtrans = u_gtrans / REAL( nzt - nzb, KIND = wp )
          v_gtrans = v_gtrans / REAL( nzt - nzb, KIND = wp )
       ELSE
!
!--       Averaging over the entire model domain.
          u_gtrans_l = 0.0_wp
          v_gtrans_l = 0.0_wp
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   u_gtrans_l = u_gtrans_l + u(k,j,i)
                   v_gtrans_l = v_gtrans_l + v(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          uv_gtrans_l(1) = u_gtrans_l / REAL( (nxr-nxl+1) * (nyn-nys+1) * (nzt-nzb), KIND = wp )
          uv_gtrans_l(2) = v_gtrans_l / REAL( (nxr-nxl+1) * (nyn-nys+1) * (nzt-nzb), KIND = wp )
#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( uv_gtrans_l, uv_gtrans, 2, MPI_REAL, MPI_SUM, comm2d, ierr )
          u_gtrans = uv_gtrans(1) / REAL( numprocs, KIND = wp )
          v_gtrans = uv_gtrans(2) / REAL( numprocs, KIND = wp )
#else
          u_gtrans = uv_gtrans_l(1)
          v_gtrans = uv_gtrans_l(2)
#endif
       ENDIF
    ENDIF

!
!-- Determine the maxima of the velocity components, including their grid index positions.
    CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, u, 'abs', 0.0_wp, u_max, u_max_ijk )
    CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, v, 'abs', 0.0_wp, v_max, v_max_ijk )
    CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, w, 'abs', 0.0_wp, w_max, w_max_ijk )

    IF ( .NOT. dt_fixed )  THEN
!
!--    Variable time step:
!--    Calculate the maximum time step according to the CFL-criterion
       dt_u_l = 999999.9_wp
       dt_v_l = 999999.9_wp
       dt_w_l = 999999.9_wp

       IF ( .NOT. rans_mode )  THEN
!
!--       Consider each velocity component individually

          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
          !$ACC COPY(dt_u_l, dt_v_l, dt_w_l, u_stokes_zu, v_stokes_zu) &
          !$ACC REDUCTION(MIN: dt_u_l, dt_v_l, dt_w_l) &
          !$ACC PRESENT(u, v, w, dzu)
          !$OMP PARALLEL DO PRIVATE(i,j,k) &
          !$OMP REDUCTION(MIN: dt_u_l, dt_v_l, dt_w_l)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   dt_u_l = MIN( dt_u_l, ( dx / ( ABS( u(k,j,i) - u_gtrans + u_stokes_zu(k) )      &
                                                  + 1.0E-10_wp ) ) )
                   dt_v_l = MIN( dt_v_l, ( dy / ( ABS( v(k,j,i) - v_gtrans + v_stokes_zu(k) )      &
                                                  + 1.0E-10_wp ) ) )
                   dt_w_l = MIN( dt_w_l, ( dzu(k) / ( ABS( w(k,j,i) ) + 1.0E-10_wp ) ) )
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Consider the wind speed at the scalar-grid point
!--       !> @note Considering the wind speed instead of each individual wind component is only a
!--       !>       workaround so far. This has to be changed in the future.

          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
          !$ACC COPY(dt_u_l, u_stokes_zu, v_stokes_zu) &
          !$ACC REDUCTION(MIN: dt_u_l) &
          !$ACC PRESENT(u, v, w, dzu)
          !$OMP PARALLEL DO PRIVATE(i,j,k) &
          !$OMP REDUCTION(MIN: dt_u_l)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   dt_u_l = MIN( dt_u_l, ( MIN( dx, dy, dzu(k) ) / ( SQRT(                         &
                            ( 0.5 * ( u(k,j,i) + u(k,j,i+1) ) - u_gtrans + u_stokes_zu(k) )**2     &
                          + ( 0.5 * ( v(k,j,i) + v(k,j+1,i) ) - v_gtrans + v_stokes_zu(k) )**2     &
                          + ( 0.5 * ( w(k,j,i) + w(k-1,j,i) ) )**2 ) + 1.0E-10_wp ) ) )
                ENDDO
             ENDDO
          ENDDO

          dt_v_l = dt_u_l
          dt_w_l = dt_u_l

       ENDIF

#if defined( __parallel )
       reduce_l(1) = dt_u_l
       reduce_l(2) = dt_v_l
       reduce_l(3) = dt_w_l
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( reduce_l, reduce, 3, MPI_REAL, MPI_MIN, comm2d, ierr )
       dt_u = reduce(1)
       dt_v = reduce(2)
       dt_w = reduce(3)
#else
       dt_u = dt_u_l
       dt_v = dt_v_l
       dt_w = dt_w_l
#endif

!
!--    Compute time step according to the diffusion criterion.
!--    First calculate minimum grid spacing which only depends on index k. When using the dynamic
!--    subgrid model, negative km are possible.
       dt_diff_l = 999999.0_wp

       !$ACC PARALLEL LOOP PRESENT(dxyz2_min, dzw)
       DO  k = nzb+1, nzt
           dxyz2_min(k) = MIN( dx2, dy2, dzw(k) * dzw(k) ) * 0.125_wp
       ENDDO

       !$OMP PARALLEL DO private(i,j,k) reduction(MIN: dt_diff_l)
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC COPY(dt_diff_l) REDUCTION(MIN: dt_diff_l) &
       !$ACC PRESENT(dxyz2_min, kh, km)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                dt_diff_l = MIN( dt_diff_l, dxyz2_min(k) / ( MAX( kh(k,j,i), 2.0_wp *              &
                                 ABS( km(k,j,i) ) ) + 1E-20_wp ) )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( dt_diff_l, dt_diff, 1, MPI_REAL, MPI_MIN, comm2d, ierr )
#else
       dt_diff = dt_diff_l
#endif
!
!--    Determine further timestep-limitations. This, e.g., includes limitations by the energy-
!--    balance models or precipitation.
       CALL module_interface_timestep
!
!--    The time step is the minimum of the 3 velocity components according to the CFL
!--    criterion, the diffusion time step, the precipitation time step, as well as
!--    the maximum allowed timestep of the land- and urban-surface model where a
!--    diffusion criterion is applied.
!--    The resulting time step is further reducted by scaling with cfl_factor, just to be on
!--    the safe side. Further on, the time step must not exceed the maximum allowed value.
!--    Note, only in very rare cases the timestep is limited by the land- or urban-surface
!--    models, e.g. when metal surfaces appear.
       dt_3d = cfl_factor * MIN( dt_diff, dt_u, dt_v, dt_w, dt_precipitation, dt_lsm, dt_usm )
!
!--    Remember the restricting time step criterion for later output.
!--    Limitation by user-defined (or default) setting of dt_max.
       IF ( dt_max <= dt_3d )  THEN
          timestep_reason = 'X'
!
!--    Advection.
       ELSEIF ( MIN( dt_u, dt_v, dt_w ) < MIN( dt_diff, dt_precipitation, dt_lsm, dt_usm ) )  THEN
          timestep_reason = 'A'
!
!--    Diffusion.
       ELSEIF ( dt_diff < MIN( dt_u, dt_v, dt_w, dt_precipitation, dt_lsm, dt_usm ) )  THEN
          timestep_reason = 'D'
!
!--    Surface energy balance models.
       ELSEIF ( MIN( dt_lsm, dt_usm ) < MIN( dt_u, dt_v, dt_w, dt_precipitation, dt_diff ) )  THEN
          timestep_reason = 'S'
!
!--    Precipitation.
       ELSEIF ( bulk_cloud_model  .AND.                                                            &
                dt_precipitation < MIN( dt_u, dt_v, dt_w, dt_diff, dt_lsm, dt_usm ) )  THEN
          timestep_reason = 'P'
       ENDIF
!
!--    Finally, limit the timestep by the user-defined (or default) setting of dt_max.
       dt_3d = MIN( dt_3d, dt_max )
!
!--    Set flag if the time step becomes too small.
       IF ( dt_3d < ( 0.00001_wp * dt_max ) )  THEN
          stop_dt = .TRUE.

!
!--       Determine the maxima of the diffusion coefficients, including their grid index positions.
          CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, km, 'abs', 0.0_wp, km_max,      &
                               km_max_ijk )
          CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, kh, 'abs', 0.0_wp, kh_max,      &
                               kh_max_ijk )

          WRITE( message_string, * ) 'Time step has reached minimum limit.',                       &
               '&dt              = ', dt_3d, ' s  Simulation is terminated.',                      &
               '&dt_u            = ', dt_u, ' s',                                                  &
               '&dt_v            = ', dt_v, ' s',                                                  &
               '&dt_w            = ', dt_w, ' s',                                                  &
               '&dt_diff         = ', dt_diff, ' s',                                               &
               '&dt_lsm          = ', dt_lsm, ' s',                                                &
               '&dt_usm          = ', dt_usm, ' s',                                                &
               '&u_max           = ', u_max, ' m/s    k=', u_max_ijk(1),                           &
               '  j=', u_max_ijk(2), '  i=', u_max_ijk(3),                                         &
               '&v_max           = ', v_max, ' m/s    k=', v_max_ijk(1),                           &
               '  j=', v_max_ijk(2), '  i=', v_max_ijk(3),                                         &
               '&w_max           = ', w_max, ' m/s    k=', w_max_ijk(1),                           &
               '  j=', w_max_ijk(2), '  i=', w_max_ijk(3),                                         &
               '&km_max          = ', km_max, ' m2/s2  k=', km_max_ijk(1),                         &
               '  j=', km_max_ijk(2), '  i=', km_max_ijk(3),                                       &
               '&kh_max          = ', kh_max, ' m2/s2  k=', kh_max_ijk(1),                         &
                '  j=', kh_max_ijk(2), '  i=', kh_max_ijk(3)
          CALL message( 'timestep', 'PAC0318', 0, 1, 0, 6, 0 )
       ENDIF

!
!--    In case of nested runs all parent/child processes have to terminate if one process has set
!--    the stop flag, i.e. they need to set the stop flag too.
       IF ( nested_run )  THEN
          stop_dt_local = stop_dt
#if defined( __parallel )
          CALL MPI_ALLREDUCE( stop_dt_local, stop_dt, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
#endif
       ENDIF

!
!--    Ensure a smooth value (two significant digits) of the timestep.
       div = 1000.0_wp
       DO  WHILE ( dt_3d < div )
          div = div / 10.0_wp
       ENDDO
       dt_3d = NINT( dt_3d * 100.0_wp / div ) * div / 100.0_wp

    ELSE
!
!--    This is the branch for a fixed time step given by the user.
       timestep_reason = 'F'
    ENDIF

    CALL cpu_log( log_point(12), 'calculate_timestep', 'stop' )

 END SUBROUTINE timestep
