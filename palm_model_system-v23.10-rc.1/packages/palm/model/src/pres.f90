!> @file pres.f90
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
!> Compute the divergence of the provisional velocity field. Solve the Poisson equation for the
!> perturbation pressure. Compute the final velocities using this perturbation pressure. Compute the
!> remaining divergence.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pres

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  d,                                                                                  &
               ddzu,                                                                               &
               ddzu_pres,                                                                          &
               ddzw,                                                                               &
               drho_air,                                                                           &
               dzw,                                                                                &
               p,                                                                                  &
               p_loc,                                                                              &
               rho_air,                                                                            &
               rho_air_zw,                                                                         &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               w

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               child_domain,                                                                       &
               conserve_volume_flow,                                                               &
               dt_3d,                                                                              &
               gathered_size,                                                                      &
               ibc_p_b,                                                                            &
               ibc_p_t,                                                                            &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               mg_switch_to_pe0_level,                                                             &
               nesting_offline,                                                                    &
               psolver,                                                                            &
               subdomain_size,                                                                     &
               topography,                                                                         &
               use_sm_for_poisfft,                                                                 &
               volume_flow,                                                                        &
               volume_flow_area,                                                                   &
               volume_flow_initial

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE diagnostic_output_quantities_mod,                                                          &
        ONLY:  div_new,                                                                            &
               div_old

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               ngp_2dh_wgrid,                                                                      &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxl_mg,                                                                             &
               nxr,                                                                                &
               nxrg,                                                                               &
               nxr_mg,                                                                             &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nys_mg,                                                                             &
               nyn,                                                                                &
               nyng,                                                                               &
               nyn_mg,                                                                             &
               nzb,                                                                                &
               nzt,                                                                                &
               nzt_mg,                                                                             &
               topo_flags

    USE kinds

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  nesting_bounds

    USE poisfft_mod,                                                                               &
        ONLY:  poisfft

#if defined( __parallel )
    USE poisfft_sm_mod,                                                                            &
        ONLY:  poisfft_sm
#endif

    USE poismg_mod

    USE poismg_noopt_mod

    USE statistics,                                                                                &
        ONLY:  statistic_regions,                                                                  &
               sums_divnew_l,                                                                      &
               sums_divold_l,                                                                      &
               weight_pres,                                                                        &
               weight_substep

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    REAL(wp) ::  ddt_3d            !<
    REAL(wp) ::  d_weight_pres     !<
    REAL(wp) ::  localsum          !<
    REAL(wp) ::  threadsum         !<
    REAL(wp) ::  weight_pres_l     !<
    REAL(wp) ::  weight_substep_l  !<

    REAL(wp), DIMENSION(1:3)   ::  volume_flow_l       !<
    REAL(wp), DIMENSION(1:3)   ::  volume_flow_offset  !<
    REAL(wp), DIMENSION(nzb:nzt) ::  w_l               !<
    REAL(wp), DIMENSION(nzb:nzt) ::  w_l_l             !<


    CALL cpu_log( log_point(8), 'pres', 'start' )

!
!-- Calculate quantities to be used locally
    ddt_3d = 1.0_wp / dt_3d
    IF ( intermediate_timestep_count == 0 )  THEN
!
!--    If pres is called before initial time step
       weight_pres_l    = 1.0_wp
       d_weight_pres    = 1.0_wp
       weight_substep_l = 1.0_wp
    ELSE
       weight_pres_l    = weight_pres(intermediate_timestep_count)
       d_weight_pres    = 1.0_wp / weight_pres(intermediate_timestep_count)
       weight_substep_l = weight_substep(intermediate_timestep_count)
    ENDIF

!
!-- Multigrid method expects array d to have one ghost layer.
!--
    IF ( psolver(1:9) == 'multigrid' )  THEN

       DEALLOCATE( d )
       ALLOCATE( d(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) )

!
!--    Since p is later used to hold the weighted average of the substeps, it cannot be used in the
!--    iterative solver. Therefore, its initial value is stored on p_loc, which is then iteratively
!--    advanced in every substep. Attention: Because PALM solves the anelastic system of equations,
!--    p_loc is defined as the perturbation pressure divided by the density.
       IF ( intermediate_timestep_count <= 1 )  THEN
          DO  i = nxl-1, nxr+1
             DO  j = nys-1, nyn+1
                DO  k = nzb, nzt+1
                   p_loc(k,j,i) = p(k,j,i) * drho_air(k)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

    ELSEIF ( psolver == 'sor'  .AND.  intermediate_timestep_count <= 1 )  THEN

!
!--    Since p is later used to hold the weighted average of the substeps, it cannot be used in the
!--    iterative solver. Therefore, its initial value is stored on p_loc, which is then iteratively
!--    advanced in every substep. Attention: Because PALM solves the anelastic system of equations,
!--    p_loc is defined as the perturbation pressure divided by the density.
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                p_loc(k,j,i) = p(k,j,i) * drho_air(k)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- Conserve the volume flow at the outflow in case of non-cyclic lateral boundary conditions
!-- WARNING: so far, this conservation does not work at the left/south boundary if the topography at
!--          the inflow differs from that at the outflow! For this case, volume_flow_area needs
!--          adjustment!
!
!-- Left/right
    IF ( conserve_volume_flow  .AND.  ( bc_radiation_l .OR. bc_radiation_r ) )  THEN

       volume_flow(1)   = 0.0_wp
       volume_flow_l(1) = 0.0_wp

       IF ( bc_radiation_l )  THEN
          i = 0
       ELSEIF ( bc_radiation_r )  THEN
          i = nx+1
       ENDIF

       DO  j = nys, nyn
!
!--       Sum up the volume flow through the south/north boundary
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1) + u(k,j,i) * dzw(k)                               &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l(1), volume_flow(1), 1, MPI_REAL, MPI_SUM, comm1dy, ierr )
#else
       volume_flow = volume_flow_l
#endif
       volume_flow_offset(1) = ( volume_flow_initial(1) - volume_flow(1) ) / volume_flow_area(1)

       DO  j = nysg, nyng
          DO  k = nzb+1, nzt
             u(k,j,i) = u(k,j,i) + volume_flow_offset(1)                                           &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO

    ENDIF

!
!-- South/north
    IF ( conserve_volume_flow  .AND.  ( bc_radiation_n .OR. bc_radiation_s ) )  THEN

       volume_flow(2)   = 0.0_wp
       volume_flow_l(2) = 0.0_wp

       IF ( bc_radiation_s )  THEN
          j = 0
       ELSEIF ( bc_radiation_n )  THEN
          j = ny+1
       ENDIF

       DO  i = nxl, nxr
!
!--       Sum up the volume flow through the south/north boundary
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2) + v(k,j,i) * dzw(k)                               &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l(2), volume_flow(2), 1, MPI_REAL, MPI_SUM, comm1dx, ierr )
#else
       volume_flow = volume_flow_l
#endif
       volume_flow_offset(2) = ( volume_flow_initial(2) - volume_flow(2) ) / volume_flow_area(2)

       DO  i = nxlg, nxrg
          DO  k = nzb+1, nzt
             v(k,j,i) = v(k,j,i) + volume_flow_offset(2)                                           &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
       ENDDO

    ENDIF

!
!-- Remove mean vertical velocity in case that Neumann conditions are used both at bottom and top
!-- boundary. With Neumann conditions at both vertical boundaries, the solver cannot remove mean
!-- vertical velocities. They should be removed, because incompressibility requires that the
!-- vertical gradient of vertical velocity is zero. Since w=0 at the solid surface, it must be zero
!-- everywhere.
!-- This must not be done in case of a 3d-nesting child domain or in case of pure vertical nesting
!-- with non-cyclic conditions, because a mean vertical velocity can physically exist in such a
!-- domain.
!-- Also in case of offline nesting, mean vertical velocities may exist (and must not be removed),
!-- caused by horizontal divergence/convergence of the large scale flow that is prescribed at the
!-- side boundaries.
!-- The removal cannot be done before the first initial time step because ngp_2dh_wgrid is not yet
!-- known then.
    IF ( ibc_p_b == 1  .AND.  ibc_p_t == 1  .AND.  .NOT. nesting_offline                           &
         .AND. .NOT. ( child_domain .AND. ( nesting_bounds /= 'vertical_only'  .OR.  &
              ( nesting_bounds == 'vertical_only' .AND. .NOT. ( bc_lr_cyc .AND. bc_ns_cyc ) ) ) )  &
         .AND. intermediate_timestep_count /= 0 )  THEN
       w_l = 0.0_wp;  w_l_l = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
                w_l_l(k) = w_l_l(k) + w(k,j,i) *                                                   &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( w_l_l(1), w_l(1), nzt, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       w_l = w_l_l
#endif
       DO  k = 1, nzt
          w_l(k) = w_l(k) / ngp_2dh_wgrid(k)
       ENDDO
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
                w(k,j,i) = w(k,j,i) - w_l(k) *                                                     &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
             ENDDO
          ENDDO
       ENDDO
!
!--    Instead of running the above loop over ghost points, they are set via exchange_horiz,
!--    in order to correctly consider non-cyclic boundary conditions, where ghost boundaries
!--    of the total domain must not be set. Otherwise w may continuously increase/decrease
!--    at these points.
       CALL exchange_horiz( w, nbgp )
    ENDIF

!
!-- Compute the divergence of the provisional velocity field.
    CALL cpu_log( log_point_s(1), 'divergence', 'start' )

    IF ( psolver(1:9) == 'multigrid' )  THEN
       !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                d(k,j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDDO
    ELSE
       !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(d)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                d(k,j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    localsum  = 0.0_wp
    threadsum = 0.0_wp

    !$OMP PARALLEL PRIVATE (i,j,k)
    !$OMP DO SCHEDULE( STATIC )
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
    !$ACC PRESENT(u, v, w, rho_air, rho_air_zw, ddzw, topo_flags) &
    !$ACC PRESENT(d)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = 1, nzt
             d(k,j,i) = ( ( u(k,j,i+1) - u(k,j,i) ) * rho_air(k) * ddx +                           &
                          ( v(k,j+1,i) - v(k,j,i) ) * rho_air(k) * ddy +                           &
                          ( w(k,j,i)   * rho_air_zw(k) - w(k-1,j,i) * rho_air_zw(k-1) )            &
                          * ddzw(k) ) * ddt_3d * d_weight_pres                                     &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL

!
!-- Compute possible PE-sum of divergences for flow_statistics. Carry out computation only at last
!-- Runge-Kutta substep.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.                      &
         intermediate_timestep_count == 0 )  THEN
       !$OMP PARALLEL PRIVATE (i,j,k) FIRSTPRIVATE(threadsum) REDUCTION(+:localsum)
       !$OMP DO SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i,j,k) &
       !$ACC REDUCTION(+:threadsum) COPY(threadsum) &
       !$ACC PRESENT(d)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                threadsum = threadsum + ABS( d(k,j,i) )
             ENDDO
          ENDDO
       ENDDO
       localsum = localsum + threadsum * dt_3d * weight_pres_l
       !$OMP END PARALLEL
    ENDIF

!
!-- For completeness, set the divergence sum of all statistic regions to those of the total domain
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.                      &
         intermediate_timestep_count == 0 )                                                        &
    THEN
       sums_divold_l(0:statistic_regions) = localsum
    ENDIF

!
!-- Store the "old" divergence for diagnostic output, only for the last Runge-Kutta timestep.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.                      &
         intermediate_timestep_count == 0 )                                                        &
    THEN
       IF ( ALLOCATED( div_old ) )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, nzt
                   div_old(k,j,i) = d(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point_s(1), 'divergence', 'stop' )

!
!-- Compute the pressure perturbation solving the Poisson equation.
    IF ( psolver(1:7) == 'poisfft' )  THEN

!
!--    Solve Poisson equation via FFT and solution of tridiagonal matrices.
       IF ( use_sm_for_poisfft )  THEN
#if defined( __parallel )
          CALL poisfft_sm( d )
#endif
       ELSE
          CALL poisfft( d )
       ENDIF

!
!--    Store computed perturbation pressure and set boundary condition in z-direction.
       !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(d, tend)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = d(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Set Neumann boundary conditions for pressure in non-cyclic case
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  tend(:,:,nxl-1) = tend(:,:,nxl)
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  tend(:,:,nxr+1) = tend(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  tend(:,nyn+1,:) = tend(:,nyn,:)
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  tend(:,nys-1,:) = tend(:,nys,:)
       ENDIF
!
!--    Bottom boundary:
!--    This condition is only required for internal output. The pressure gradient
!--    (dp(nzb+1)-dp(nzb))/dz is not used anywhere else.
!--    Neumann
       IF ( ibc_p_b == 1 )  THEN
          !$OMP PARALLEL DO PRIVATE (i,j) SCHEDULE( STATIC )
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                tend(nzb,j,i) = tend(nzb+1,j,i)
             ENDDO
          ENDDO
!
!--    Dirichlet
       ELSE
          !$OMP PARALLEL DO PRIVATE (i,j) SCHEDULE( STATIC )
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                tend(nzb,j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary
       IF ( ibc_p_t == 1 )  THEN
!
!--       Neumann
          !$OMP PARALLEL DO PRIVATE (i,j) SCHEDULE( STATIC )
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                tend(nzt+1,j,i) = tend(nzt,j,i)
             ENDDO
          ENDDO

       ELSE
!
!--       Dirichlet
          !$OMP PARALLEL DO PRIVATE (i,j) SCHEDULE( STATIC )
          !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
          !$ACC PRESENT(tend)
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                tend(nzt+1,j,i) = 0.0_wp
             ENDDO
          ENDDO

       ENDIF

!
!--    Exchange boundaries for p
       CALL exchange_horiz( tend, nbgp )

    ELSEIF ( psolver == 'sor' )  THEN

!
!--    Solve Poisson equation for perturbation pressure using SOR-Red/Black scheme
       CALL sor( d, ddzu_pres, ddzw, p_loc )
       tend = p_loc

    ELSEIF ( psolver(1:9) == 'multigrid' )  THEN

!
!--    Solve Poisson equation for perturbation pressure using Multigrid scheme, array tend is used
!--    to store the residuals.

!--    If the number of grid points of the gathered grid, which is collected on PE0, is larger than
!--    the number of grid points of an PE, than array tend will be enlarged.
       IF ( gathered_size > subdomain_size )  THEN
          DEALLOCATE( tend )
          ALLOCATE( tend(nzb:nzt_mg(mg_switch_to_pe0_level)+1,nys_mg(                              &
                    mg_switch_to_pe0_level)-1:nyn_mg(mg_switch_to_pe0_level)+1,                    &
                    nxl_mg(mg_switch_to_pe0_level)-1:nxr_mg(mg_switch_to_pe0_level)+1) )
       ENDIF

       IF ( psolver == 'multigrid' )  THEN
          CALL poismg( tend )
       ELSE
          CALL poismg_noopt( tend )
       ENDIF

       IF ( gathered_size > subdomain_size )  THEN
          DEALLOCATE( tend )
          ALLOCATE( tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

!
!--    Restore perturbation pressure on tend because this array is used further below to correct the
!--    velocity fields
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                tend(k,j,i) = p_loc(k,j,i)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- Store perturbation pressure on array p, used for pressure data output.
!-- Attention: Because PALM solves the anelastic system of equations, tend
!-- contains the perturbation pressure divided by the density.
!-- Ghost layers are added in the output routines (except sor-method: see below)
    IF ( intermediate_timestep_count <= 1 )  THEN
       !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(p, tend)
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                p(k,j,i) = tend(k,j,i)  * rho_air(k) * weight_substep_l
             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( intermediate_timestep_count > 1 )  THEN
       !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(p, tend)
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                p(k,j,i) = p(k,j,i) + tend(k,j,i) * rho_air(k) * weight_substep_l
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- SOR-method needs ghost layers for the next timestep
    IF ( psolver == 'sor' )  CALL exchange_horiz( p, nbgp )

!
!-- Correction of the provisional velocities with the current perturbation pressure just computed
    IF ( conserve_volume_flow  .AND.  ( bc_lr_cyc .OR. bc_ns_cyc ) )  THEN
       volume_flow_l(1) = 0.0_wp
       volume_flow_l(2) = 0.0_wp
    ENDIF
!
!-- Add pressure gradients to the velocity components. Note, the loops are running over the entire
!-- model domain, even though, in case of non-cyclic boundaries u- and v-component are not
!-- prognostic at i=0 and j=0, respectiveley. However, in case of Dirichlet boundary conditions for
!-- the velocities, zero-gradient conditions for the pressure are set, so that no modification is
!-- imposed at the boundaries.
    !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k) &
    !$ACC PRESENT(u, v, w, tend, ddzu, topo_flags)
    DO  i = nxl, nxr
       DO  j = nys, nyn

          DO  k = nzb+1, nzt
             w(k,j,i) = w(k,j,i) - dt_3d * ( tend(k+1,j,i) - tend(k,j,i) ) * ddzu(k+1)             &
                        * weight_pres_l                                                            &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
          ENDDO

          DO  k = nzb+1, nzt
             u(k,j,i) = u(k,j,i) - dt_3d * ( tend(k,j,i) - tend(k,j,i-1) ) * ddx * weight_pres_l   &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO

          DO  k = nzb+1, nzt
             v(k,j,i) = v(k,j,i) - dt_3d * ( tend(k,j,i) - tend(k,j-1,i) ) * ddy * weight_pres_l   &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO

       ENDDO
    ENDDO

!
!-- The vertical velocity may be non zero at the top of nested child domains. w(nzt+1) is set to
!-- w(nzt) in routine pmci_ensure_nest_mass_conservation BEFORE calling the pressure solver.
!-- Same is done below, because w(nzt) has been changed above, to avoid jumps in the profile output
!-- of w. Hint: w level nzt+1 does not impact results.
    IF ( child_domain )  THEN
       w(nzt+1,:,:) = w(nzt,:,:)
    ENDIF

!
!-- Sum up the volume flow through the right and north boundary
    IF ( conserve_volume_flow  .AND.  bc_lr_cyc  .AND.  bc_ns_cyc  .AND.  nxr == nx )  THEN
       threadsum = 0.0_wp
!
!--    Summation of the volume flow is done on threadsum rather than on volumen_flow itself.
!--    This is because intel compiler, when compiled with openmp, do not allow reduction
!--    operation on array elements.
       !$OMP PARALLEL DO PRIVATE (j,k) REDUCTION (+:threadsum)
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             threadsum = threadsum + u(k,j,nxr) * dzw(k)                                           &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,nxr), 1 ) )
          ENDDO
       ENDDO
       volume_flow_l(1) = threadsum

    ENDIF

    IF ( conserve_volume_flow  .AND.  bc_ns_cyc  .AND.  bc_lr_cyc  .AND. nyn == ny )  THEN
       threadsum = 0.0_wp
       !$OMP PARALLEL DO PRIVATE (j,k) REDUCTION (+:threadsum)
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             threadsum = threadsum + v(k,nyn,i) * dzw(k)                                           &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,nyn,i), 2 ) )
           ENDDO
       ENDDO
       volume_flow_l(2) = threadsum

    ENDIF

!
!-- Conserve the volume flow
    IF ( conserve_volume_flow  .AND.  ( bc_lr_cyc  .AND.  bc_ns_cyc ) )  THEN

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l(1), volume_flow(1), 2, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       volume_flow = volume_flow_l
#endif

       volume_flow_offset(1:2) = ( volume_flow_initial(1:2) - volume_flow(1:2) )                   &
                                 / volume_flow_area(1:2)

       !$OMP PARALLEL DO PRIVATE (i,j,k) SCHEDULE( STATIC )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                u(k,j,i) = u(k,j,i) + volume_flow_offset(1)                                        &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
             ENDDO
             DO  k = nzb+1, nzt
                v(k,j,i) = v(k,j,i) + volume_flow_offset(2)                                        &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- Exchange of boundaries for the velocities
    CALL exchange_horiz( u, nbgp )
    CALL exchange_horiz( v, nbgp )
    CALL exchange_horiz( w, nbgp )

!
!-- Compute the divergence of the corrected velocity field.
!-- A possible PE-sum is computed in flow_statistics. Carry out computation only at last
!-- Runge-Kutta step.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.                      &
         intermediate_timestep_count == 0 )  THEN
       CALL cpu_log( log_point_s(1), 'divergence', 'start' )
       sums_divnew_l = 0.0_wp

!
!--    d must be reset to zero because it can contain nonzero values below the topography
       IF ( topography /= 'flat' )  d = 0.0_wp

       localsum  = 0.0_wp
       threadsum = 0.0_wp

       !$OMP PARALLEL PRIVATE (i,j,k) FIRSTPRIVATE(threadsum) REDUCTION(+:localsum)
       !$OMP DO SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(u, v, w, rho_air, rho_air_zw, ddzw, topo_flags) &
       !$ACC PRESENT(d)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                d(k,j,i) = ( ( u(k,j,i+1) - u(k,j,i) ) * rho_air(k) * ddx +    &
                             ( v(k,j+1,i) - v(k,j,i) ) * rho_air(k) * ddy +    &
                             ( w(k,j,i)   * rho_air_zw(k) - w(k-1,j,i) * rho_air_zw(k-1) )         &
                             * ddzw(k) )                                                           &
                             * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             ENDDO
          ENDDO
       ENDDO
!
!--    Compute possible PE-sum of divergences for flow_statistics
       !$OMP DO SCHEDULE( STATIC )
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC REDUCTION(+:threadsum) COPY(threadsum) &
       !$ACC PRESENT(d)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                threadsum = threadsum + ABS( d(k,j,i) )
             ENDDO
          ENDDO
       ENDDO

       localsum = localsum + threadsum
       !$OMP END PARALLEL

!
!--    For completeness, set the divergence sum of all statistic regions to those of the total
!--    domain
       sums_divnew_l(0:statistic_regions) = localsum

!
!--    Store the "new" divergence for diagnostic output.
       IF ( ALLOCATED( div_new ) )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, nzt
                   div_new(k,j,i) = d(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       CALL cpu_log( log_point_s(1), 'divergence', 'stop' )

    ENDIF

    CALL cpu_log( log_point(8), 'pres', 'stop' )

 END SUBROUTINE pres
